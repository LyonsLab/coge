require({cache:{
'JBrowse/Store/TabixIndexedFile':function(){
define("JBrowse/Store/TabixIndexedFile", [
           'dojo/_base/declare',
           'dojo/_base/array',
           'JBrowse/Util',
           'JBrowse/Store/LRUCache',
           'JBrowse/Errors',
           'JBrowse/Model/XHRBlob',
           'JBrowse/Model/BGZip/BGZBlob',
           'JBrowse/Model/TabixIndex'
       ],
       function(
           declare,
           array,
           Util,
           LRUCache,
           Errors,
           XHRBlob,
           BGZBlob,
           TabixIndex
       ) {

return declare( null, {

    constructor: function( args ) {
        this.browser = args.browser;
        this.index = new TabixIndex({ blob: new BGZBlob( args.tbi ), browser: args.browser } );
        this.data  = new BGZBlob( args.file );
        this.indexLoaded = this.index.load();

        this.chunkSizeLimit = args.chunkSizeLimit || 15000000;
    },

    getLines: function( ref, min, max, itemCallback, finishCallback, errorCallback ) {
        var thisB = this;
        var args = Array.prototype.slice.call(arguments);
        this.indexLoaded.then(function() {
            thisB._fetch.apply( thisB, args );
        }, errorCallback);
    },

    _fetch: function( ref, min, max, itemCallback, finishCallback, errorCallback ) {
        errorCallback = errorCallback || function(e) { console.error(e, e.stack); };

        var chunks = this.index.blocksForRange( ref, min, max);
        if ( ! chunks ) {
            errorCallback('Error in index fetch ('+[ref,min,max].join(',')+')');
            return;
        }

        // toString function is used by the cache for making cache keys
        chunks.toString = chunks.toUniqueString = function() {
            return this.join(', ');
        };

        // check the chunks for any that are over the size limit.  if
        // any are, don't fetch any of them
        for( var i = 0; i<chunks.length; i++ ) {
            var size = chunks[i].fetchedSize();
            if( size > this.chunkSizeLimit ) {
                errorCallback( new Errors.DataOverflow('Too much data. Chunk size '+Util.commifyNumber(size)+' bytes exceeds chunkSizeLimit of '+Util.commifyNumber(this.chunkSizeLimit)+'.' ) );
                return;
            }
        }

        var fetchError;
        try {
            this._fetchChunkData(
                chunks,
                ref,
                min,
                max,
                itemCallback,
                finishCallback,
                errorCallback
            );
        } catch( e ) {
            errorCallback( e );
        }
    },

    _fetchChunkData: function( chunks, ref, min, max, itemCallback, endCallback, errorCallback ) {
        var thisB = this;

        if( ! chunks.length ) {
            endCallback();
            return;
        }

        var allItems = [];
        var chunksProcessed = 0;

        var cache = this.chunkCache = this.chunkCache || new LRUCache({
            name: 'TabixIndexedFileChunkedCache',
            fillCallback: dojo.hitch( this, '_readChunkItems' ),
            sizeFunction: function( chunkItems ) {
                return chunkItems.length;
            },
            maxSize: 100000 // cache up to 100,000 items
        });

        var regRef = this.browser.regularizeReferenceName( ref );

        var haveError;
        array.forEach( chunks, function( c ) {
            cache.get( c, function( chunkItems, e ) {
                if( e && !haveError )
                    errorCallback( e );
                if(( haveError = haveError || e )) {
                    return;
                }

                for( var i = 0; i< chunkItems.length; i++ ) {
                    var item = chunkItems[i];
                    if( item._regularizedRef == regRef ) {
                        // on the right ref seq
                        if( item.start > max ) // past end of range, can stop iterating
                            break;
                        else if( item.end >= min ) // must be in range
                            itemCallback( item );
                    }
                }
                if( ++chunksProcessed == chunks.length ) {
                    endCallback();
                }
            });
        });
    },

    _readChunkItems: function( chunk, callback ) {
        var thisB = this;
        var items = [];

        thisB.data.read(chunk.minv.block, chunk.maxv.block - chunk.minv.block + 1, function( data ) {
            data = new Uint8Array(data);

            // throw away the first (probably incomplete) line
            var parseStart = array.indexOf( data, thisB._newlineCode, 0 ) + 1;

            try {
                thisB.parseItems(
                    data,
                    parseStart,
                    function(i) { items.push(i); },
                    function() { callback(items); }
                );
            } catch( e ) {
                callback( null, e );
            }
        },
        function(e) {
            callback( null, e );
        });
    },

    parseItems: function( data, blockStart, itemCallback, finishCallback ) {
        var that = this;
        var itemCount = 0;

        var maxItemsWithoutYielding = 300;
        var parseState = { data: data, offset: blockStart };

        while ( true ) {
            // if we've read no more than a certain number of items this cycle, read another one
            if( itemCount <= maxItemsWithoutYielding ) {
                var item = this.parseItem( parseState ); //< increments parseState.offset
                if( item ) {
                    itemCallback( item );
                    itemCount++;
                }
                else {
                    finishCallback();
                    return;
                }
            }
            // if we're not done but we've read a good chunk of
            // items, schedule the rest of our work in a timeout to continue
            // later, avoiding blocking any UI stuff that needs to be done
            else {
                window.setTimeout( function() {
                    that.parseItems( data, parseState.offset, itemCallback, finishCallback );
                }, 1);
                return;
            }
        }
    },

    // stub method, override in subclasses or instances
    parseItem: function( parseState ) {
        var metaChar = this.index.metaChar;

        var line;
        do {
            line = this._getline( parseState );
        } while( line && line[0] == metaChar );

        if( !line )
            return null;

        // function extractColumn( colNum ) {
        //     var skips = '';
        //     while( colNum-- > 1 )
        //         skips += '^[^\t]*\t';
        //     var match = (new Regexp( skips+'([^\t]*)' )).exec( line );
        //     if( ! match )
        //         return null;
        //     return match[1];
        // }
        var fields = line.split( "\t" );
        var item = { // note: index column numbers are 1-based
            ref:   fields[this.index.columnNumbers.ref-1],
            _regularizedRef: this.browser.regularizeReferenceName( fields[this.index.columnNumbers.ref-1] ),
            start: parseInt(fields[this.index.columnNumbers.start-1]),
            end:   parseInt(fields[this.index.columnNumbers.end-1]),
            fields: fields
        };
        return item;
    },

    _newlineCode: "\n".charCodeAt(0),

    _getline: function( parseState ) {
        var data = parseState.data;
        var newlineIndex = array.indexOf( data, this._newlineCode, parseState.offset );

        if( newlineIndex == -1 ) // no more lines
            return null;

        var line = '';
        for( var i = parseState.offset; i < newlineIndex; i++ )
            line += String.fromCharCode( data[i] );
        parseState.offset = newlineIndex+1;
        return line;
    }
});
});

},
'JBrowse/Model/BGZip/BGZBlob':function(){
/**
 * File blob in Heng Li's `bgzip` format.
 */
define( [
            'dojo/_base/declare',
            'jszlib/inflate',
            'jszlib/arrayCopy'
        ],
        function(
            declare,
            inflate,
            arrayCopy
        ) {

var BGZBlob = declare( null,
{
    constructor: function( blob ) {
        this.blob = blob;
    },

    blockSize: 1<<16,

    slice: function(s, l) {
        return new BGZBlob( this.blob.slice( s, l ) );
    },

    fetch: function( callback, failCallback ) {
        this.blob.fetch(
            this._wrap( callback ),
            failCallback
        );
    },

    read: function( offset, length, callback, failCallback ) {
        this.blob.read( offset,
                        length + this.blockSize, //< need to over-fetch by a whole block size
                        this._wrap( callback, length ),
                        failCallback
                      );
    },

    _wrap: function( callback, maxLen ) {
        var thisB = this;
        return function( bgzData ) {
            callback( thisB.unbgzf( bgzData, maxLen ) );
        };
    },

    readInt: function(ba, offset) {
        return (ba[offset + 3] << 24) | (ba[offset + 2] << 16) | (ba[offset + 1] << 8) | (ba[offset]);
    },

    readShort: function(ba, offset) {
        return (ba[offset + 1] << 8) | (ba[offset]);
    },

    readFloat: function(ba, offset) {
        var temp = new Uint8Array( 4 );
        for( var i = 0; i<4; i++ ) {
            temp[i] = ba[offset+i];
        }
        var fa = new Float32Array( temp.buffer );
        return fa[0];
    },

    unbgzf: function(data, lim) {
        lim = Math.min( lim || Infinity, data.byteLength - 27);
        var oBlockList = [];
        var totalSize = 0;

        for( var ptr = [0]; ptr[0] < lim; ptr[0] += 8) {

            var ba = new Uint8Array( data, ptr[0], 18 );

            // check the bgzf block magic
            if( !( ba[0] == 31 && ba[1] == 139 ) ) {
                console.error( 'invalid BGZF block header, skipping', ba );
                break;
            }

            var xlen = this.readShort( ba, 10 );
            var compressedDataOffset = ptr[0] + 12 + xlen;

            // var inPtr = ptr[0];
            // var bSize = Utils.readShort( ba, 16 );
            // var logLength = Math.min(data.byteLength-ptr[0], 40);
            // console.log( xlen, bSize, bSize - xlen - 19, new Uint8Array( data, ptr[0], logLength ), logLength );

            var unc;
            try {
                unc = inflate(
                    data,
                    compressedDataOffset,
                    data.byteLength - compressedDataOffset,
                    ptr
                );
            } catch( inflateError ) {
                // if we have a buffer error and we have already
                // inflated some data, there is probably just an
                // incomplete BGZF block at the end of the data, so
                // ignore it and stop inflating
                if( /^Z_BUF_ERROR/.test(inflateError.statusString) && oBlockList.length ) {
                    break;
                }
                // otherwise it's some other kind of real error
                else {
                    throw inflateError;
                }
            }
            if( unc.byteLength ) {
                totalSize += unc.byteLength;
                oBlockList.push( unc );
            }
            // else {
            //     console.error( 'BGZF decompression failed for block ', compressedDataOffset, data.byteLength-compressedDataOffset, [inPtr] );
            // }
        }

        if (oBlockList.length == 1) {
            return oBlockList[0];
        } else {
            var out = new Uint8Array(totalSize);
            var cursor = 0;
            for (var i = 0; i < oBlockList.length; ++i) {
                var b = new Uint8Array(oBlockList[i]);
                arrayCopy(b, 0, out, cursor, b.length);
                cursor += b.length;
            }
            return out.buffer;
        }
    }



});

return BGZBlob;
});
},
'JBrowse/Model/TabixIndex':function(){
define([
           'dojo/_base/declare',
           'dojo/_base/array',
           'dojo/_base/Deferred',
           'JBrowse/has',
           'jDataView',
           'JBrowse/Util',
           'JBrowse/Model/BGZip/VirtualOffset'
       ],
       function(
           declare,
           array,
           Deferred,
           has,
           jDataView,
           Util,
           VirtualOffset
       ) {

// inner class representing a chunk
var Chunk = Util.fastDeclare({
    constructor: function(minv,maxv,bin) {
        this.minv = minv;
        this.maxv = maxv;
        this.bin = bin;
    },
    toUniqueString: function() {
        return this.minv+'..'+this.maxv+' (bin '+this.bin+')';
    },
    toString: function() {
        return this.toUniqueString();
    },
    compareTo: function( b ) {
        return this.minv - b.minv || this.maxv - b.maxv || this.bin - b.bin;
    },
    compare: function( b ) {
        return this.compareTo( b );
    },
    fetchedSize: function() {
        return this.maxv.block + (1<<16) - this.minv.block + 1;
    }
});

return declare( null, {

   constructor: function( args ) {
       this.browser = args.browser;
       this.blob = args.blob;
       this.load();
   },

   load: function() {
       var thisB = this;
       return this._loaded = this._loaded || function() {
           var d = new Deferred();
           if( ! has('typed-arrays') )
               d.reject( 'This web browser does not support JavaScript typed arrays.' );
           else
               this.blob.fetch( function( data) {
                                    thisB._parseIndex( data, d );
                                }, dojo.hitch( d, 'reject' ) );
           return d;
       }.call(this);
   },

   // fetch and parse the index
   _parseIndex: function( bytes, deferred ) {

       this._littleEndian = true;
       var data = new jDataView( bytes, 0, undefined, this._littleEndian );

       // check TBI magic numbers
       if( data.getInt32() != 21578324 /* "TBI\1" */) {
           // try the other endianness if no magic
           this._littleEndian = false;
           data = new jDataView( bytes, 0, undefined, this._littleEndian );
           if( data.getInt32() != 21578324 /* "TBI\1" */) {
               console.error('Not a TBI file');
               deferred.reject('Not a TBI file');
               return;
           }
       }

       // number of reference sequences in the index
       var refCount = data.getInt32();
       this.presetType = data.getInt32();
       this.columnNumbers = {
           ref:   data.getInt32(),
           start: data.getInt32(),
           end:   data.getInt32()
       };
       this.metaValue = data.getInt32();
       this.metaChar = this.metaValue ? String.fromCharCode( this.metaValue ) : null;
       this.skipLines = data.getInt32();

       // read sequence dictionary
       this._refIDToName = new Array( refCount );
       this._refNameToID = {};
       var nameSectionLength = data.getInt32();
       this._parseNameBytes( data.getBytes( nameSectionLength, undefined, false ) );

       // read the per-reference-sequence indexes
       this._indices = new Array( refCount );
       for (var i = 0; i < refCount; ++i) {
           // the binning index
           var binCount = data.getInt32();
           var idx = this._indices[i] = { binIndex: {} };
           for (var j = 0; j < binCount; ++j) {
               var bin        = data.getInt32();
               var chunkCount = data.getInt32();
               var chunks = new Array( chunkCount );
               for (var k = 0; k < chunkCount; ++k) {
                   var u = new VirtualOffset( data.getBytes(8) );
                   var v = new VirtualOffset( data.getBytes(8) );
                   this._findFirstData( u );
                   chunks[k] = new Chunk( u, v, bin );
               }
               idx.binIndex[bin] = chunks;
           }
           // the linear index
           var linearCount = data.getInt32();
           var linear = idx.linearIndex = new Array( linearCount );
           for (var k = 0; k < linearCount; ++k) {
               linear[k] = new VirtualOffset( data.getBytes(8) );
               this._findFirstData( linear[k] );
           }
       }
       deferred.resolve({ success: true });
   },

   _findFirstData: function( virtualOffset ) {
       var fdl = this.firstDataLine;
       this.firstDataLine = fdl ? fdl.compareTo( virtualOffset ) > 0 ? virtualOffset
                                                                     : fdl
                                : virtualOffset;
   },

   _parseNameBytes: function( namesBytes ) {
       var offset = 0;

       function getChar() {
           var b = namesBytes[ offset++ ];
           return b ? String.fromCharCode( b ) : null;
       }

       function getString() {
           var c, s = '';
           while(( c = getChar() ))
               s += c;
           return s.length ? s : null;
       }

       var refName, refID = 0;
       for( ; refName = getString(); refID++ ) {
           this._refIDToName[refID] = refName;
           this._refNameToID[ this.browser.regularizeReferenceName( refName ) ] = refID;
       }
   },

    /**
     * Interrogate whether a store has data for a given reference
     * sequence.  Calls the given callback with either true or false.
     *
     * Implemented as a binary interrogation because some stores are
     * smart enough to regularize reference sequence names, while
     * others are not.
     */
    hasRefSeq: function( seqName, callback, errorCallback ) {
       var thisB = this;
       seqName = thisB.browser.regularizeReferenceName( seqName );
       thisB.load().then( function() {
           if( seqName in thisB._refNameToID ) {
               callback(true);
               return;
           }
           callback( false );
       });
   },

   getRefId: function( refName ) {
       refName = this.browser.regularizeReferenceName( refName );
       return this._refNameToID[refName];
   },

   TAD_LIDX_SHIFT: 14,

   blocksForRange: function( refName, beg, end ) {
       if( beg < 0 )
           beg = 0;

       var tid = this.getRefId( refName );
       var indexes = this._indices[tid];
       if( ! indexes )
           return [];

       var linearIndex = indexes.linearIndex,
            binIndex   = indexes.binIndex;

       var bins = this._reg2bins(beg, end);

       var min_off = linearIndex.length
           ? linearIndex[
                 ( beg >> this.TAD_LIDX_SHIFT >= linearIndex.length )
                     ?  linearIndex.length - 1
                     :  beg >> this.TAD_LIDX_SHIFT
               ]
           : new VirtualOffset( 0, 0 );

       var i, l, n_off = 0;
       for( i = 0; i < bins.length; ++i ) {
           n_off += ( binIndex[ bins[i] ] || [] ).length;
       }

       if( n_off == 0 )
           return [];

       var off = [];

       var chunks;
       for (i = n_off = 0; i < bins.length; ++i)
           if (( chunks = binIndex[ bins[i] ] ))
               for (var j = 0; j < chunks.length; ++j)
                   if( min_off.compareTo( chunks[j].maxv ) < 0 )
                       off[n_off++] = new Chunk( chunks[j].minv, chunks[j].maxv, chunks[j].bin );

       if( ! off.length )
           return [];

       off = off.sort( function(a,b) {
                           return a.compareTo(b);
                       });

       // resolve completely contained adjacent blocks
       for (i = 1, l = 0; i < n_off; ++i) {
           if( off[l].maxv.compareTo( off[i].maxv ) < 0 ) {
               ++l;
               off[l].minv = off[i].minv;
               off[l].maxv = off[i].maxv;
           }
       }
       n_off = l + 1;

       // resolve overlaps between adjacent blocks; this may happen due to the merge in indexing
       for (i = 1; i < n_off; ++i)
           if ( off[i-1].maxv >= off[i].minv )
               off[i-1].maxv = off[i].minv;
       // merge adjacent blocks
       for (i = 1, l = 0; i < n_off; ++i) {
           if( off[l].maxv.block == off[i].minv.block )
               off[l].maxv = off[i].maxv;
           else {
               ++l;
               off[l].minv = off[i].minv;
               off[l].maxv = off[i].maxv;
           }
       }
       n_off = l + 1;

       return off.slice( 0, n_off );
   },

    /* calculate bin given an alignment covering [beg,end) (zero-based, half-close-half-open) */
    _reg2bin: function(beg, end) {
        --end;
        if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
        if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
        if (beg>>20 == end>>20) return ((1<<9)-1)/7 + (beg>>20);
        if (beg>>23 == end>>23) return ((1<<6)-1)/7 + (beg>>23);
        if (beg>>26 == end>>26) return ((1<<3)-1)/7 + (beg>>26);
        return 0;
    },

    /* calculate the list of bins that may overlap with region [beg,end) (zero-based) */
    _reg2bins: function(beg, end) {
        var k, list = [];
        --end;
        list.push(0);
        for (k = 1 + (beg>>26); k <= 1 + (end>>26); ++k) list.push(k);
        for (k = 9 + (beg>>23); k <= 9 + (end>>23); ++k) list.push(k);
        for (k = 73 + (beg>>20); k <= 73 + (end>>20); ++k) list.push(k);
        for (k = 585 + (beg>>17); k <= 585 + (end>>17); ++k) list.push(k);
        for (k = 4681 + (beg>>14); k <= 4681 + (end>>14); ++k) list.push(k);
        return list;
    }

});
});

},
'jDataView/jdataview':function(){
define([], function() {
var scope = {};

//
// jDataView by Vjeux - Jan 2010
//
// A unique way to read a binary file in the browser
// http://github.com/vjeux/jDataView
// http://blog.vjeux.com/ <vjeuxx@gmail.com>
//

(function (global) {

var compatibility = {
	ArrayBuffer: typeof ArrayBuffer !== 'undefined',
	DataView: typeof DataView !== 'undefined' &&
		('getFloat64' in DataView.prototype ||				// Chrome
		 'getFloat64' in new DataView(new ArrayBuffer(1))), // Node
	// NodeJS Buffer in v0.5.5 and newer
	NodeBuffer: typeof Buffer !== 'undefined' && 'readInt16LE' in Buffer.prototype
};

var dataTypes = {
	'Int8': 1,
	'Int16': 2,
	'Int32': 4,
	'Uint8': 1,
	'Uint16': 2,
	'Uint32': 4,
	'Float32': 4,
	'Float64': 8
};

var nodeNaming = {
	'Int8': 'Int8',
	'Int16': 'Int16',
	'Int32': 'Int32',
	'Uint8': 'UInt8',
	'Uint16': 'UInt16',
	'Uint32': 'UInt32',
	'Float32': 'Float',
	'Float64': 'Double'
};

var jDataView = function (buffer, byteOffset, byteLength, littleEndian) {
	if (!(this instanceof jDataView)) {
		throw new Error("jDataView constructor may not be called as a function");
	}

	this.buffer = buffer;

	// Handle Type Errors
	if (!(compatibility.NodeBuffer && buffer instanceof Buffer) &&
		!(compatibility.ArrayBuffer && buffer instanceof ArrayBuffer) &&
		typeof buffer !== 'string') {
		throw new TypeError('jDataView buffer has an incompatible type');
	}

	// Check parameters and existing functionnalities
	this._isArrayBuffer = compatibility.ArrayBuffer && buffer instanceof ArrayBuffer;
	this._isDataView = compatibility.DataView && this._isArrayBuffer;
	this._isNodeBuffer = compatibility.NodeBuffer && buffer instanceof Buffer;

	// Default Values
	this._littleEndian = Boolean(littleEndian);

	var bufferLength = this._isArrayBuffer ? buffer.byteLength : buffer.length;
	if (byteOffset === undefined) {
		byteOffset = 0;
	}
	this.byteOffset = byteOffset;

	if (byteLength === undefined) {
		byteLength = bufferLength - byteOffset;
	}
	this.byteLength = byteLength;

	if (!this._isDataView) {
		// Do additional checks to simulate DataView
		if (typeof byteOffset !== 'number') {
			throw new TypeError('jDataView byteOffset is not a number');
		}
		if (typeof byteLength !== 'number') {
			throw new TypeError('jDataView byteLength is not a number');
		}
		if (byteOffset < 0) {
			throw new Error('jDataView byteOffset is negative');
		}
		if (byteLength < 0) {
			throw new Error('jDataView byteLength is negative');
		}
	}

	// Instanciate
	if (this._isDataView) {
		this._view = new DataView(buffer, byteOffset, byteLength);
		this._start = 0;
	}
	this._start = byteOffset;
	if (byteOffset + byteLength > bufferLength) {
		throw new Error("jDataView (byteOffset + byteLength) value is out of bounds");
	}

	this._offset = 0;

	// Create uniform reading methods (wrappers) for the following data types

	if (this._isDataView) { // DataView: we use the direct method
		for (var type in dataTypes) {
			if (!dataTypes.hasOwnProperty(type)) {
				continue;
			}
			(function(type, view){
				var size = dataTypes[type];
				view['get' + type] = function (byteOffset, littleEndian) {
					// Handle the lack of endianness
					if (littleEndian === undefined) {
						littleEndian = view._littleEndian;
					}

					// Handle the lack of byteOffset
					if (byteOffset === undefined) {
						byteOffset = view._offset;
					}

					// Move the internal offset forward
					view._offset = byteOffset + size;

					return view._view['get' + type](byteOffset, littleEndian);
				};
			})(type, this);
		}
	} else if (this._isNodeBuffer && compatibility.NodeBuffer) {
		for (var type in dataTypes) {
			if (!dataTypes.hasOwnProperty(type)) {
				continue;
			}

			var name;
			if (type === 'Int8' || type === 'Uint8') {
				name = 'read' + nodeNaming[type];
			} else if (littleEndian) {
				name = 'read' + nodeNaming[type] + 'LE';
			} else {
				name = 'read' + nodeNaming[type] + 'BE';
			}

			(function(type, view, name){
				var size = dataTypes[type];
				view['get' + type] = function (byteOffset, littleEndian) {
					// Handle the lack of endianness
					if (littleEndian === undefined) {
						littleEndian = view._littleEndian;
					}

					// Handle the lack of byteOffset
					if (byteOffset === undefined) {
						byteOffset = view._offset;
					}

					// Move the internal offset forward
					view._offset = byteOffset + size;

					return view.buffer[name](view._start + byteOffset);
				};
			})(type, this, name);
		}
	} else {
		for (var type in dataTypes) {
			if (!dataTypes.hasOwnProperty(type)) {
				continue;
			}
			(function(type, view){
				var size = dataTypes[type];
				view['get' + type] = function (byteOffset, littleEndian) {
					// Handle the lack of endianness
					if (littleEndian === undefined) {
						littleEndian = view._littleEndian;
					}

					// Handle the lack of byteOffset
					if (byteOffset === undefined) {
						byteOffset = view._offset;
					}

					// Move the internal offset forward
					view._offset = byteOffset + size;

					if (view._isArrayBuffer && (view._start + byteOffset) % size === 0 && (size === 1 || littleEndian)) {
						// ArrayBuffer: we use a typed array of size 1 if the alignment is good
						// ArrayBuffer does not support endianess flag (for size > 1)
						return new global[type + 'Array'](view.buffer, view._start + byteOffset, 1)[0];
					} else {
						// Error checking:
						if (typeof byteOffset !== 'number') {
							throw new TypeError('jDataView byteOffset is not a number');
						}
						if (byteOffset + size > view.byteLength) {
							throw new Error('jDataView (byteOffset + size) value is out of bounds');
						}

						return view['_get' + type](view._start + byteOffset, littleEndian);
					}
				};
			})(type, this);
		}
	}
};

if (compatibility.NodeBuffer) {
	jDataView.createBuffer = function () {
		return new Buffer(arguments);
	};
} else if (compatibility.ArrayBuffer) {
	jDataView.createBuffer = function () {
		return new Uint8Array(arguments).buffer;
	};
} else {
	jDataView.createBuffer = function () {
		return String.fromCharCode.apply(null, arguments);
	};
}

jDataView.prototype = {
	compatibility: compatibility,

	// Helpers

	_getBytes: function (length, byteOffset, littleEndian) {
		var result;

		// Handle the lack of endianness
		if (littleEndian === undefined) {
			littleEndian = this._littleEndian;
		}

		// Handle the lack of byteOffset
		if (byteOffset === undefined) {
			byteOffset = this._offset;
		}

		// Error Checking
		if (typeof byteOffset !== 'number') {
			throw new TypeError('jDataView byteOffset is not a number');
		}
		if (length < 0 || byteOffset + length > this.byteLength) {
			throw new Error('jDataView length or (byteOffset+length) value is out of bounds');
		}

		byteOffset += this._start;

		if (this._isArrayBuffer) {
			result = new Uint8Array(this.buffer, byteOffset, length);
		}
		else {
			result = this.buffer.slice(byteOffset, byteOffset + length);

			if (!this._isNodeBuffer) {
				result = Array.prototype.map.call(result, function (ch) {
					return ch.charCodeAt(0) & 0xff;
				});
			}
		}

		if (littleEndian && length > 1) {
			if (!(result instanceof Array)) {
				result = Array.prototype.slice.call(result);
			}

			result.reverse();
		}

		this._offset = byteOffset - this._start + length;

		return result;
	},

	// wrapper for external calls (do not return inner buffer directly to prevent it's modifying)
	getBytes: function (length, byteOffset, littleEndian) {
		var result = this._getBytes.apply(this, arguments);

		if (!(result instanceof Array)) {
			result = Array.prototype.slice.call(result);
		}

		return result;
	},

	getString: function (length, byteOffset) {
		var value;

		if (this._isNodeBuffer) {
			// Handle the lack of byteOffset
			if (byteOffset === undefined) {
				byteOffset = this._offset;
			}

			// Error Checking
			if (typeof byteOffset !== 'number') {
				throw new TypeError('jDataView byteOffset is not a number');
			}
			if (length < 0 || byteOffset + length > this.byteLength) {
				throw new Error('jDataView length or (byteOffset+length) value is out of bounds');
			}

			value = this.buffer.toString('ascii', this._start + byteOffset, this._start + byteOffset + length);
			this._offset = byteOffset + length;
		}
		else {
			value = String.fromCharCode.apply(null, this._getBytes(length, byteOffset, false));
		}

		return value;
	},

	getChar: function (byteOffset) {
		return this.getString(1, byteOffset);
	},

	tell: function () {
		return this._offset;
	},

	seek: function (byteOffset) {
		if (typeof byteOffset !== 'number') {
			throw new TypeError('jDataView byteOffset is not a number');
		}
		if (byteOffset < 0 || byteOffset > this.byteLength) {
			throw new Error('jDataView byteOffset value is out of bounds');
		}

		return this._offset = byteOffset;
	},

	// Compatibility functions on a String Buffer

	_getFloat64: function (byteOffset, littleEndian) {
		var b = this._getBytes(8, byteOffset, littleEndian),

			sign = 1 - (2 * (b[0] >> 7)),
			exponent = ((((b[0] << 1) & 0xff) << 3) | (b[1] >> 4)) - (Math.pow(2, 10) - 1),

		// Binary operators such as | and << operate on 32 bit values, using + and Math.pow(2) instead
			mantissa = ((b[1] & 0x0f) * Math.pow(2, 48)) + (b[2] * Math.pow(2, 40)) + (b[3] * Math.pow(2, 32)) +
						(b[4] * Math.pow(2, 24)) + (b[5] * Math.pow(2, 16)) + (b[6] * Math.pow(2, 8)) + b[7];

		if (exponent === 1024) {
			if (mantissa !== 0) {
				return NaN;
			} else {
				return sign * Infinity;
			}
		}

		if (exponent === -1023) { // Denormalized
			return sign * mantissa * Math.pow(2, -1022 - 52);
		}

		return sign * (1 + mantissa * Math.pow(2, -52)) * Math.pow(2, exponent);
	},

	_getFloat32: function (byteOffset, littleEndian) {
		var b = this._getBytes(4, byteOffset, littleEndian),

			sign = 1 - (2 * (b[0] >> 7)),
			exponent = (((b[0] << 1) & 0xff) | (b[1] >> 7)) - 127,
			mantissa = ((b[1] & 0x7f) << 16) | (b[2] << 8) | b[3];

		if (exponent === 128) {
			if (mantissa !== 0) {
				return NaN;
			} else {
				return sign * Infinity;
			}
		}

		if (exponent === -127) { // Denormalized
			return sign * mantissa * Math.pow(2, -126 - 23);
		}

		return sign * (1 + mantissa * Math.pow(2, -23)) * Math.pow(2, exponent);
	},

	_getInt32: function (byteOffset, littleEndian) {
		var b = this._getUint32(byteOffset, littleEndian);
		return b > Math.pow(2, 31) - 1 ? b - Math.pow(2, 32) : b;
	},

	_getUint32: function (byteOffset, littleEndian) {
		var b = this._getBytes(4, byteOffset, littleEndian);
		return (b[0] * Math.pow(2, 24)) + (b[1] << 16) + (b[2] << 8) + b[3];
	},

	_getInt16: function (byteOffset, littleEndian) {
		var b = this._getUint16(byteOffset, littleEndian);
		return b > Math.pow(2, 15) - 1 ? b - Math.pow(2, 16) : b;
	},

	_getUint16: function (byteOffset, littleEndian) {
		var b = this._getBytes(2, byteOffset, littleEndian);
		return (b[0] << 8) + b[1];
	},

	_getInt8: function (byteOffset) {
		var b = this._getUint8(byteOffset);
		return b > Math.pow(2, 7) - 1 ? b - Math.pow(2, 8) : b;
	},

	_getUint8: function (byteOffset) {
		return this._getBytes(1, byteOffset)[0];
	}
};

if (typeof jQuery !== 'undefined' && jQuery.fn.jquery >= "1.6.2") {
	var convertResponseBodyToText = function (byteArray) {
		// http://jsperf.com/vbscript-binary-download/6
		var scrambledStr;
		try {
			scrambledStr = IEBinaryToArray_ByteStr(byteArray);
		} catch (e) {
			// http://stackoverflow.com/questions/1919972/how-do-i-access-xhr-responsebody-for-binary-data-from-javascript-in-ie
			// http://miskun.com/javascript/internet-explorer-and-binary-files-data-access/
			var IEBinaryToArray_ByteStr_Script =
				"Function IEBinaryToArray_ByteStr(Binary)\r\n"+
				"	IEBinaryToArray_ByteStr = CStr(Binary)\r\n"+
				"End Function\r\n"+
				"Function IEBinaryToArray_ByteStr_Last(Binary)\r\n"+
				"	Dim lastIndex\r\n"+
				"	lastIndex = LenB(Binary)\r\n"+
				"	if lastIndex mod 2 Then\r\n"+
				"		IEBinaryToArray_ByteStr_Last = AscB( MidB( Binary, lastIndex, 1 ) )\r\n"+
				"	Else\r\n"+
				"		IEBinaryToArray_ByteStr_Last = -1\r\n"+
				"	End If\r\n"+
				"End Function\r\n";

			// http://msdn.microsoft.com/en-us/library/ms536420(v=vs.85).aspx
			// proprietary IE function
			window.execScript(IEBinaryToArray_ByteStr_Script, 'vbscript');

			scrambledStr = IEBinaryToArray_ByteStr(byteArray);
		}

		var lastChr = IEBinaryToArray_ByteStr_Last(byteArray),
		result = "",
		i = 0,
		l = scrambledStr.length % 8,
		thischar;
		while (i < l) {
			thischar = scrambledStr.charCodeAt(i++);
			result += String.fromCharCode(thischar & 0xff, thischar >> 8);
		}
		l = scrambledStr.length;
		while (i < l) {
			result += String.fromCharCode(
				(thischar = scrambledStr.charCodeAt(i++), thischar & 0xff), thischar >> 8,
				(thischar = scrambledStr.charCodeAt(i++), thischar & 0xff), thischar >> 8,
				(thischar = scrambledStr.charCodeAt(i++), thischar & 0xff), thischar >> 8,
				(thischar = scrambledStr.charCodeAt(i++), thischar & 0xff), thischar >> 8,
				(thischar = scrambledStr.charCodeAt(i++), thischar & 0xff), thischar >> 8,
				(thischar = scrambledStr.charCodeAt(i++), thischar & 0xff), thischar >> 8,
				(thischar = scrambledStr.charCodeAt(i++), thischar & 0xff), thischar >> 8,
				(thischar = scrambledStr.charCodeAt(i++), thischar & 0xff), thischar >> 8);
		}
		if (lastChr > -1) {
			result += String.fromCharCode(lastChr);
		}
		return result;
	};

	jQuery.ajaxSetup({
		converters: {
			'* dataview': function(data) {
				return new jDataView(data);
			}
		},
		accepts: {
			dataview: "text/plain; charset=x-user-defined"
		},
		responseHandler: {
			dataview: function (responses, options, xhr) {
				// Array Buffer Firefox
				if ('mozResponseArrayBuffer' in xhr) {
					responses.text = xhr.mozResponseArrayBuffer;
				}
				// Array Buffer Chrome
				else if ('responseType' in xhr && xhr.responseType === 'arraybuffer' && xhr.response) {
					responses.text = xhr.response;
				}
				// Internet Explorer (Byte array accessible through VBScript -- convert to text)
				else if ('responseBody' in xhr) {
					responses.text = convertResponseBodyToText(xhr.responseBody);
				}
				// Older Browsers
				else {
					responses.text = xhr.responseText;
				}
			}
		}
	});

	jQuery.ajaxPrefilter('dataview', function(options, originalOptions, jqXHR) {
		// trying to set the responseType on IE 6 causes an error
		if (jQuery.support.ajaxResponseType) {
			if (!options.hasOwnProperty('xhrFields')) {
				options.xhrFields = {};
			}
			options.xhrFields.responseType = 'arraybuffer';
		}
		options.mimeType = 'text/plain; charset=x-user-defined';
	});
}

global.jDataView = (global.module || {}).exports = jDataView;
if (typeof module !== 'undefined') {
	module.exports = jDataView;
}

})(scope);

return scope.jDataView;
});
},
'JBrowse/Model/BGZip/VirtualOffset':function(){
/**
 * a virtual offset into a bgzipped file
 */
define([
         'JBrowse/Util'
       ],
       function( Util ) {

var VirtualOffset = Util.fastDeclare({
    constructor: function(b, o) {
        if( arguments.length >= 2 ) {
            this.block  = b;
            this.offset = o;
        }
        else {
            this._fromBytes( b );
        }
    },

    _fromBytes: function( ba, offset ) {
        offset = offset || 0;

        //console.log( 'readVob', offset );
        var block =
              ba[offset  ] * 0x10000000000
            + ba[offset+1] * 0x100000000
            + ba[offset+2] * 0x1000000
            + ba[offset+3] * 0x10000
            + ba[offset+4] * 0x100
            + ba[offset+5];
        var bint = (ba[offset+6] << 8) | ba[offset+7];
        if (block == 0 && bint == 0) {
            this.block = this.offset = null;
        } else {
            this.block = block;
            this.offset = bint;
        }
    },
    toString: function() {
        return '' + this.block + ':' + this.offset;
    },
    compareTo: function(b) {
        return this.block - b.block || this.offset - b.offset;
    },
    cmp: function(b) {
        return this.compareTo( b );
    }
});

return VirtualOffset;

});
},
'JBrowse/Store/SeqFeature/GlobalStatsEstimationMixin':function(){
/**
 * Mixin that adds _estimateGlobalStats method to a store, which
 * samples a section of the features in the store and uses those to
 * esimate the statistics of the whole data set.
 */

define([
           'dojo/_base/declare',
           'dojo/_base/array'
       ],
       function( declare, array ) {

return declare( null, {

    /**
     * Fetch a region of the current reference sequence and use it to
     * estimate the feature density of the store.
     * @private
     */
    _estimateGlobalStats: function( finishCallback, errorCallback ) {

        var statsFromInterval = function( refSeq, length, callback ) {
            var thisB = this;
            var sampleCenter = refSeq.start*0.75 + refSeq.end*0.25;
            var start = Math.max( 0, Math.round( sampleCenter - length/2 ) );
            var end = Math.min( Math.round( sampleCenter + length/2 ), refSeq.end );
            var features = [];
            this._getFeatures({ ref: refSeq.name, start: start, end: end},
                              function( f ) { features.push(f); },
                              function( error ) {
                                  features = array.filter( features, function(f) { return f.get('start') >= start && f.get('end') <= end; } );
                                  callback.call( thisB, length,
                                                 {
                                                     featureDensity: features.length / length,
                                                     _statsSampleFeatures: features.length,
                                                     _statsSampleInterval: { ref: refSeq.name, start: start, end: end, length: length }
                                                 });
                              },
                              function( error ) {
                                      console.error( error );
                                      callback.call( thisB, length,  null, error );
                              });
        };

        var maybeRecordStats = function( interval, stats, error ) {
            if( error ) {
                finishCallback( null, error );
            } else {
                var refLen = this.refSeq.end - this.refSeq.start;
                 if( stats._statsSampleFeatures >= 300 || interval * 2 > refLen || error ) {
                     console.log( 'Store statistics: '+(this.source||this.name), stats );
                     finishCallback( stats );
                 } else {
                     statsFromInterval.call( this, this.refSeq, interval * 2, maybeRecordStats );
                 }
            }
        };

        statsFromInterval.call( this, this.refSeq, 100, maybeRecordStats );
    }

});
});
},
'JBrowse/Store/SeqFeature/VCFTabix/Parser':function(){
define([
           'dojo/_base/declare',
           'dojo/_base/array',
           'dojo/json',
           'JBrowse/Digest/Crc32',
           './LazyFeature'
       ],
       function(
           declare,
           array,
           JSON,
           Digest,
           LazyFeature
       ) {

return declare( null, {

    /**
     * Parse the bytes that contain the VCF header, storing the parsed
     * data in this.header.
     */
    parseHeader: function( headerBytes ) {

        // parse the header lines
        var headData = {};
        var parseState = { data: headerBytes, offset: 0 };
        var line;
        while(( line = this._getlineFromBytes( parseState ))) {
            // only interested in meta and header lines
            if( line[0] != '#' )
                continue;

            var match = /^##([^\s#=]+)=(.+)/.exec( line);
            // parse meta line
            if( match && match[1] ) {
                var metaField = match[1].toLowerCase();
                var metaData = (match[2]||'');

                // TODO: do further parsing for some fields
                if( metaData.match(/^<.+>$/) ) {
                    metaData = this._parseGenericHeaderLine( metaData );
                }

                if( ! headData[metaField] )
                    headData[metaField] = [];

                headData[metaField].push( metaData );
            }
            else if( /^#CHROM\t/.test( line ) ) {
                var f = line.split("\t");
                if( f[8] == 'FORMAT' && f.length > 9 )
                    headData.samples = f.slice(9);
            }
        }
        //console.log(headData);

        // index some of the headers by ID
        for( var headerType in headData ) {
            if( dojo.isArray( headData[headerType] ) && typeof headData[headerType][0] == 'object' && 'id' in headData[headerType][0] )
                headData[headerType] = this._indexUniqObjects( headData[headerType], 'id' );
        }

        this.header = headData;
        return headData;
    },

    /**
     * Given a line from a TabixIndexedFile, convert it into a feature
     * and return it.  Assumes that the header has already been parsed
     * and stored (i.e. _parseHeader has already been called.)
     */
    lineToFeature: function( line ) {
        var fields = line.fields;
        var ids = [];
        for( var i=0; i<fields.length; i++ )
            if( fields[i] == '.' )
                fields[i] = null;

        var ref = fields[3];
        var alt = fields[4];

        var SO_term = this._find_SO_term( ref, alt );
        var featureData = {
            start:  line.start,
            end:    line.start+ref.length,
            seq_id: line.ref,
            description: this._makeDescriptionString( SO_term, ref, alt ),
            type:   SO_term,
            reference_allele:    ref
        };

        if( fields[2] !== null ) {
            ids = (fields[2]||'').split(';');
            featureData.name = ids[0];
            if( ids.length > 1 )
                featureData.aliases = ids.slice(1).join(',');
        }

        if( fields[5] !== null )
            featureData.score = parseFloat( fields[5] );
        if( fields[6] !== null )
            featureData.filter = fields[6];

        if( alt && alt[0] != '<' )
            featureData.alternative_alleles = {
                meta: {
                    description: 'VCF ALT field, list of alternate non-reference alleles called on at least one of the samples'
                },
                values: alt
            };

        // parse the info field and store its contents as attributes in featureData
        this._parseInfoField( featureData, fields );

        var f = new LazyFeature({
            id: ids[0] || fields.slice( 0, 9 ).join('/'),
            data: featureData,
            fields: fields,
            parser: this
        });

        return f;
    },

    _newlineCode: "\n".charCodeAt(0),

    /**
     *  helper method that parses the next line from a Uint8Array or similar.
     *  @param parseState.data the byte array
     *  @param parseState.offset the offset to start parsing.  <MODIFIED AS A SIDE EFFECT OF THI SMETHOD
     */
    _getlineFromBytes: function( parseState ) {
        if( ! parseState.offset )
            parseState.offset = 0;

        var newlineIndex = array.indexOf( parseState.data, this._newlineCode, parseState.offset );

        if( newlineIndex == -1 ) // no more lines
            return null;

        var line = String.fromCharCode.apply( String, Array.prototype.slice.call( parseState.data, parseState.offset, newlineIndex ));
        parseState.offset = newlineIndex+1;
        return line;
    },


    _parseGenericHeaderLine: function( metaData ) {
        metaData = metaData.replace(/^<|>$/g,'');
        return this._parseKeyValue( metaData, ',', ';', 'lowercase' );
    },

    _vcfReservedInfoFields: {
        // from the VCF4.1 spec, http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
        AA:    { description: "ancestral allele" },
        AC:    { description: "allele count in genotypes, for each ALT allele, in the same order as listed" },
        AF:    { description: "allele frequency for each ALT allele in the same order as listed: use this when estimated from primary data, not called genotypes" },
        AN:    { description: "total number of alleles in called genotypes" },
        BQ:    { description: "RMS base quality at this position" },
        CIGAR: { description: "cigar string describing how to align an alternate allele to the reference allele" },
        DB:    { description: "dbSNP membership" },
        DP:    { description: "combined depth across samples, e.g. DP=154" },
        END:   { description: "end position of the variant described in this record (esp. for CNVs)" },
        H2:    { description: "membership in hapmap2" },
        MQ:    { description: "RMS mapping quality, e.g. MQ=52" },
        MQ0:   { description: "Number of MAPQ == 0 reads covering this record" },
        NS:    { description: "Number of samples with data" },
        SB:    { description: "strand bias at this position" },
        SOMATIC: { description: "indicates that the record is a somatic mutation, for cancer genomics" },
        VALIDATED: { description: "validated by follow-up experiment" },

        //specifically for structural variants
        "IMPRECISE": { number: 0, type: 'Flag', description: "Imprecise structural variation" },
        "NOVEL":     { number: 0, type: 'Flag',description: "Indicates a novel structural variation" },
        "END":       { number: 1, type: 'Integer', description: "End position of the variant described in this record" },

        // For precise variants, END is POS + length of REF allele -
        // 1, and the for imprecise variants the corresponding best
        // estimate.

        "SVTYPE": { number: 1, type: 'String',description: "Type of structural variant" },

        // Value should be one of DEL, INS, DUP, INV, CNV, BND. This
        // key can be derived from the REF/ALT fields but is useful
        // for filtering.

        "SVLEN": { number:'.',type: 'Integer', description: 'Difference in length between REF and ALT alleles' },

        // One value for each ALT allele. Longer ALT alleles
        // (e.g. insertions) have positive values, shorter ALT alleles
        // (e.g. deletions) have negative values.

        "CIPOS":  { number: 2, "type": 'Integer', "description": 'Confidence interval around POS for imprecise variants' },
        "CIEND":  { number: 2, "type": 'Integer', "description": "Confidence interval around END for imprecise variants" },
        "HOMLEN": {            "type": "Integer", "description": "Length of base pair identical micro-homology at event breakpoints" },
        "HOMSEQ": {            "type": "String",  "description": "Sequence of base pair identical micro-homology at event breakpoints" },
        "BKPTID": {            "type": "String",  "description": "ID of the assembled alternate allele in the assembly file" },

        // For precise variants, the consensus sequence the alternate
        // allele assembly is derivable from the REF and ALT
        // fields. However, the alternate allele assembly file may
        // contain additional information about the characteristics of
        // the alt allele contigs.

        "MEINFO":  { number:4, "type": "String", "description": "Mobile element info of the form NAME,START,END,POLARITY" },
        "METRANS": { number:4, "type": "String", "description": "Mobile element transduction info of the form CHR,START,END,POLARITY" },
        "DGVID":   { number:1, "type": "String", "description": "ID of this element in Database of Genomic Variation"},
        "DBVARID": { number:1, "type": "String", "description": "ID of this element in DBVAR"},
        "DBRIPID": { number:1, "type": "String", "description": "ID of this element in DBRIP"},
        "MATEID":  {           "type": "String", "description": "ID of mate breakends"},
        "PARID":   { number:1, "type": "String", "description": "ID of partner breakend"},
        "EVENT":   { number:1, "type": "String", "description": "ID of event associated to breakend"},
        "CILEN":   { number:2, "type": "Integer","description": "Confidence interval around the length of the inserted material between breakends"},
        "DP":      { number:1, "type": "Integer","description": "Read Depth of segment containing breakend"},
        "DPADJ":   {           "type": "Integer","description": "Read Depth of adjacency"},
        "CN":      { number:1, "type": "Integer","description": "Copy number of segment containing breakend"},
        "CNADJ":   {           "type": "Integer","description": "Copy number of adjacency"},
        "CICN":    { number:2, "type": "Integer","description": "Confidence interval around copy number for the segment"},
        "CICNADJ": {           "type": "Integer","description": "Confidence interval around copy number for the adjacency"}
    },

    _vcfStandardGenotypeFields: {
        // from the VCF4.1 spec, http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
        GT : { description: "genotype, encoded as allele values separated by either of '/' or '|'. The allele values are 0 for the reference allele (what is in the REF field), 1 for the first allele listed in ALT, 2 for the second allele list in ALT and so on. For diploid calls examples could be 0/1, 1|0, or 1/2, etc. For haploid calls, e.g. on Y, male non-pseudoautosomal X, or mitochondrion, only one allele value should be given; a triploid call might look like 0/0/1. If a call cannot be made for a sample at a given locus, '.' should be specified for each missing allele in the GT field (for example './.' for a diploid genotype and '.' for haploid genotype). The meanings of the separators are as follows (see the PS field below for more details on incorporating phasing information into the genotypes): '/' meaning genotype unphased, '|' meaning genotype phased" },
        DP : { description: "read depth at this position for this sample (Integer)" },
        FT : { description: "sample genotype filter indicating if this genotype was \"called\" (similar in concept to the FILTER field). Again, use PASS to indicate that all filters have been passed, a semi-colon separated list of codes for filters that fail, or \".\" to indicate that filters have not been applied. These values should be described in the meta-information in the same way as FILTERs (String, no white-space or semi-colons permitted)" },
        GL : { description: "genotype likelihoods comprised of comma separated floating point log10-scaled likelihoods for all possible genotypes given the set of alleles defined in the REF and ALT fields. In presence of the GT field the same ploidy is expected and the canonical order is used; without GT field, diploidy is assumed. If A is the allele in REF and B,C,... are the alleles as ordered in ALT, the ordering of genotypes for the likelihoods is given by: F(j/k) = (k*(k+1)/2)+j.  In other words, for biallelic sites the ordering is: AA,AB,BB; for triallelic sites the ordering is: AA,AB,BB,AC,BC,CC, etc.  For example: GT:GL 0/1:-323.03,-99.29,-802.53 (Floats)" },
        GLE : { description: "genotype likelihoods of heterogeneous ploidy, used in presence of uncertain copy number. For example: GLE=0:-75.22,1:-223.42,0/0:-323.03,1/0:-99.29,1/1:-802.53 (String)" },
        PL : { description: "the phred-scaled genotype likelihoods rounded to the closest integer (and otherwise defined precisely as the GL field) (Integers)" },
        GP : { description: "the phred-scaled genotype posterior probabilities (and otherwise defined precisely as the GL field); intended to store imputed genotype probabilities (Floats)" },
        GQ : { description: "conditional genotype quality, encoded as a phred quality -10log_10p(genotype call is wrong, conditioned on the site's being variant) (Integer)" },
        HQ : { description: "haplotype qualities, two comma separated phred qualities (Integers)" },
        PS : { description: "phase set.  A phase set is defined as a set of phased genotypes to which this genotype belongs.  Phased genotypes for an individual that are on the same chromosome and have the same PS value are in the same phased set.  A phase set specifies multi-marker haplotypes for the phased genotypes in the set.  All phased genotypes that do not contain a PS subfield are assumed to belong to the same phased set.  If the genotype in the GT field is unphased, the corresponding PS field is ignored.  The recommended convention is to use the position of the first variant in the set as the PS identifier (although this is not required). (Non-negative 32-bit Integer)" },
        PQ : { description: "phasing quality, the phred-scaled probability that alleles are ordered incorrectly in a heterozygote (against all other members in the phase set).  We note that we have not yet included the specific measure for precisely defining \"phasing quality\"; our intention for now is simply to reserve the PQ tag for future use as a measure of phasing quality. (Integer)" },
        EC : { description: "comma separated list of expected alternate allele counts for each alternate allele in the same order as listed in the ALT field (typically used in association analyses) (Integers)" },
        MQ : { description: "RMS mapping quality, similar to the version in the INFO field. (Integer)" }
    },

    _vcfReservedAltTypes: {
        "DEL": { description: "Deletion relative to the reference", so_term: 'deletion' },
        "INS": { description: "Insertion of novel sequence relative to the reference", so_term: 'insertion' },
        "DUP": { description: "Region of elevated copy number relative to the reference", so_term: 'copy_number_gain' },
        "INV": { description: "Inversion of reference sequence", so_term: 'inversion' },
        "CNV": { description: "Copy number variable region (may be both deletion and duplication)", so_term: 'copy_number_variation' },
        "DUP:TANDEM": { description: "Tandem duplication", so_term: 'copy_number_gain' },
        "DEL:ME": { description: "Deletion of mobile element relative to the reference" },
        "INS:ME": { description: "Insertion of a mobile element relative to the reference" }
    },

    /**
     * parse a VCF line's INFO field, storing the contents as
     * attributes in featureData
     */
    _parseInfoField: function( featureData, fields ) {
        if( !fields[7] || fields[7] == '.' )
            return;
        var info = this._parseKeyValue( fields[7] );

        // decorate the info records with references to their descriptions
        for( var field in info ) {
            if( info.hasOwnProperty( field ) ) {
                    var i = info[field] = {
                        values: info[field],
                        toString: function() { return (this.values || []).join(','); }
                    };
                    var meta = this._getInfoMeta( field );
                    if( meta )
                        i.meta = meta;
            }
        }

        dojo.mixin( featureData, info );
    },

    _getAltMeta: function( alt ) {
        return (this.header.alt||{})[alt] || this._vcfReservedAltTypes[alt];
    },
    _getInfoMeta: function( id ) {
        return (this.header.info||{})[id] || this._vcfReservedInfoFields[id];
    },
    _getFormatMeta: function( fieldname ) {
        return (this.header.format||{})[fieldname] || this._vcfStandardGenotypeFields[fieldname];
    },

    /**
     * Take an array of objects and make another object that indexes
     * them into another object for easy lookup by the given field.
     * WARNING: Values of the field must be unique.
     */
    _indexUniqObjects: function( entries, indexField, lowerCase ) {
        // index the info fields by field ID
        var items = {};
        array.forEach( entries, function( rec ) {
            var k = rec[indexField];
            if( dojo.isArray(k) )
                k = k[0];
            if( lowerCase )
                k = k.toLowerCase();
            items[ rec[indexField] ]= rec;
        });
        return items;
    },

    /**
     * Parse a VCF key-value string like DP=154;Foo="Bar; baz";MQ=52;H2 into an object like
     *  { DP: [154], Foo:['Bar',' baz'], ... }
     *
     * Done in a low-level style to properly support quoted values.  >:-{
     */
    _parseKeyValue: function( str, pairSeparator, valueSeparator, lowercaseKeys ) {
        pairSeparator  = pairSeparator  || ';';
        valueSeparator = valueSeparator || ',';

        var data = {};
        var currKey = '';
        var currValue = '';
        var state = 1;  // states: 1: read key to =, 2: read value to comma or sep, 3: read value to quote
        for( var i = 0; i < str.length; i++ ) {
            if( state == 1 ) { // read key
                if( str[i] == '=' ) {
                    if( lowercaseKeys )
                        currKey = currKey.toLowerCase();
                    data[currKey] = [];
                    state = 2;
                }
                else if( str[i] == pairSeparator ) {
                    if( lowercaseKeys )
                        currKey = currKey.toLowerCase();
                    data[currKey] = [true];
                    currKey = '';
                    state = 1;
                }
                else {
                    currKey += str[i];
                }
            }
            else if( state == 2 ) { // read value to value sep or pair sep
                if( str[i] == valueSeparator ) {
                    data[currKey].push( currValue );
                    currValue = '';
                }
                else if( str[i] == pairSeparator ) {
                    data[currKey].push( currValue );
                    currKey = '';
                    state = 1;
                    currValue = '';
                } else if( str[i] == '"' ) {
                    state = 3;
                    currValue = '';
                }
                else
                    currValue += str[i];
            }
            else if( state == 3 ) { // read value to quote
                if( str[i] != '"' )
                    currValue += str[i];
                else
                    state = 2;
            }
        }

        if( state == 2 || state == 3) {
            data[currKey].push( currValue );
        }

        return data;
    },

    _find_SO_term: function( ref, alt ) {
        // it's just a remark if there are no alternate alleles
        if( ! alt || alt == '.' )
            return 'remark';

        var types = array.filter( array.map( alt.split(','), function( alt ) {
                                                 return this._find_SO_term_from_alt_definitions( alt );
                                             }, this ),
                                  function( t ) { return t; } );

        if( types[0] )
            return types.join(',');


        return this._find_SO_term_by_examination( ref, alt );
    },

    /**
     * Given an ALT string, return a string suitable for appending to
     * the feature description, if available.
     */
    _makeDescriptionString: function( SO_term, ref, alt ) {
        if( ! alt )
            return 'no alternative alleles';

        alt = alt.replace(/^<|>$/g,'');

        var def = this._getAltMeta( alt );
        return def && def.description ? alt+' - '+def.description : SO_term+" "+ref+" -> "+ alt;
    },

    _find_SO_term_from_alt_definitions: function( alt ) {
        // not a symbolic ALT if doesn't begin with '<', so we'll have no definition
        if( alt[0] != '<' )
            return null;

        alt = alt.replace(/^<|>$/g,''); // trim off < and >

        // look for a definition with an SO type for this
        var def = (this.header.alt||{})[alt] || this._vcfReservedAltTypes[alt];
        if( def && def.so_term )
            return def.so_term;

        // try to look for a definition for a parent term if we can
        alt = alt.split(':');
        if( alt.length > 1 )
            return this._find_SO_term_from_alt_definitions( '<'+alt.slice( 0, alt.length-1 ).join(':')+'>' );
        else // no parent
            return null;
    },

    _find_SO_term_by_examination: function( ref, alt ) {
        alt = alt.split(',');

        var minAltLen = Infinity;
        var maxAltLen = -Infinity;
        var altLen = array.map( alt, function(a) {
            var l = a.length;
            if( l < minAltLen )
                minAltLen = l;
            if( l > maxAltLen )
                maxAltLen = l;
            return a.length;
        });

        if( ref.length == 1 && minAltLen == 1 && maxAltLen == 1 )
            return 'SNV'; // use SNV because SO definition of SNP says
                          // abundance must be at least 1% in
                          // population, and can't be sure we meet
                          // that

        if( ref.length == minAltLen && ref.length == maxAltLen )
            if( alt.length == 1 && ref.split('').reverse().join('') == alt[0] )
                return 'inversion';
            else
                return 'substitution';

        if( ref.length <= minAltLen && ref.length < maxAltLen )
            return 'insertion';

        if( ref.length > minAltLen && ref.length >= maxAltLen )
            return 'deletion';

        return 'indel';
    }

});
});

},
'JBrowse/Store/SeqFeature/VCFTabix/LazyFeature':function(){
/**
 * Lazy-parsing feature implementation for VCF stores.
 */

define( ['dojo/_base/array',
         'JBrowse/Util'
        ],
        function( array, Util ) {

var Feature = Util.fastDeclare(
{
    constructor: function( args ) {
        this.parser  = args.parser;
        this.data    = args.data;
        this._id = args.id;
        this.fields  = args.fields;
    },

    get: function( field) {
        return this._get( field ) || this._get( field.toLowerCase() );
    },

    // same as get(), except requires lower-case arguments.  used
    // internally to save lots of calls to field.toLowerCase()
    _get: function( field ) {
        return field in this.data ? this.data[field] : // have we already parsed it out?
            function(field) {
                var v = this.data[field] =
                    this['_parse_'+field] ? this['_parse_'+field]()            : // maybe we have a special parser for it
                                            undefined;
                return v;
            }.call(this,field);
    },

    parent: function() {
        return null;
    },

    children: function() {
        return null;
    },

    tags: function() {
        var t = [];
        var d = this.data;
        for( var k in d ) {
            if( d.hasOwnProperty( k ) )
                t.push( k );
        }
        if( ! d.genotypes )
            t.push('genotypes');
        return t;
    },

    id: function() {
        return this._id;
    },

    _parse_genotypes: function() {
        var fields = this.fields;
        var parser = this.parser;
        delete this.fields; // TODO: remove these deletes if we add other laziness
        delete this.parser;

        if( fields.length < 10 )
            return null;

        // parse the genotype data fields
        var genotypes = [];
        var format = array.map( fields[8].split(':'), function( fieldID ) {
                         return { id: fieldID, meta: parser._getFormatMeta( fieldID ) };
                     }, this );
        for( var i = 9; i < fields.length; ++i ) {
            var g = (fields[i]||'').split(':');
            var gdata = {};
            for( var j = 0; j<format.length; ++j ) {
                var gData = g[j] || '';
                gdata[ format[j].id ] = {
                    // don't split on commas if it looks like a string
                    values: gData.charAt(0) == '"' ? [ gData ] : gData.split(','),
                    meta: format[j].meta
                };
            }
            genotypes.push( gdata );
        }

        // index the genotypes by sample ID
        var bySample = {};
        for( var i = 0; i<genotypes.length; i++ ) {
            var sname = (parser.header.samples||{})[i];
            if( sname ) {
                bySample[sname] = genotypes[i];
            }
        }

        // add a toString to it that serializes it to JSON without its metadata
        bySample.toString = this._stringifySample;

        return bySample;
    },

    _stringifySample: function() {
        var ex = {};
        for( var sample in this ) {
            var srec = ex[sample] = {};
            for( var field in this[sample] ) {
                    srec[field] = this[sample][field].values;
            }
        }
        return JSON.stringify( ex );
    }

});

return Feature;
});
}}});
define("JBrowse/Store/SeqFeature/VCFTabix", [
           'dojo/_base/declare',
           'dojo/_base/Deferred',
           'JBrowse/Store/SeqFeature',
           'JBrowse/Store/DeferredStatsMixin',
           'JBrowse/Store/DeferredFeaturesMixin',
           'JBrowse/Store/TabixIndexedFile',
           'JBrowse/Store/SeqFeature/GlobalStatsEstimationMixin',
           'JBrowse/Model/XHRBlob',
           './VCFTabix/Parser'
       ],
       function(
           declare,
           Deferred,
           SeqFeatureStore,
           DeferredStatsMixin,
           DeferredFeaturesMixin,
           TabixIndexedFile,
           GlobalStatsEstimationMixin,
           XHRBlob,
           VCFParser
       ) {


// subclass the TabixIndexedFile to modify the parsed items a little
// bit so that the range filtering in TabixIndexedFile will work.  VCF
// files don't actually have an end coordinate, so we have to make it
// here.  also convert coordinates to interbase.
var VCFIndexedFile = declare( TabixIndexedFile, {
    parseItem: function() {
        var i = this.inherited( arguments );
        if( i ) {
            i.start--;
            i.end = i.start + i.fields[3].length;
        }
        return i;
    }
});

return declare( [ SeqFeatureStore, DeferredStatsMixin, DeferredFeaturesMixin, GlobalStatsEstimationMixin, VCFParser ],
{

    constructor: function( args ) {
        var thisB = this;

        var tbiBlob = args.tbi ||
            new XHRBlob(
                this.resolveUrl(
                    this.getConf('tbiUrlTemplate',[]) || this.getConf('urlTemplate',[])+'.tbi'
                )
            );

        var fileBlob = args.file ||
            new XHRBlob(
                this.resolveUrl( this.getConf('urlTemplate',[]) )
            );

        this.indexedData = new VCFIndexedFile(
            {
                tbi: tbiBlob,
                file: fileBlob,
                browser: this.browser,
                chunkSizeLimit: args.chunkSizeLimit
            });

        this._loadHeader().then( function() {
            thisB._estimateGlobalStats( function( stats, error ) {
                if( error )
                    thisB._failAllDeferred( error );
                else {
                    thisB.globalStats = stats;
                    thisB._deferred.stats.resolve({success:true});
                    thisB._deferred.features.resolve({success:true});
                }
            });
        },
        function(error) {
            thisB._failAllDeferred( error );
        });
    },

    /** fetch and parse the VCF header lines */
    _loadHeader: function() {
        var thisB = this;
        return this._parsedHeader = this._parsedHeader || function() {
            var d = new Deferred();

            thisB.indexedData.indexLoaded.then( function() {
                var maxFetch = thisB.indexedData.index.firstDataLine
                    ? thisB.indexedData.index.firstDataLine.block + thisB.indexedData.data.blockSize - 1
                    : null;

                thisB.indexedData.data.read(
                    0,
                    maxFetch,
                    function( bytes ) {

                        thisB.parseHeader( new Uint8Array( bytes ) );

                        d.resolve({ success:true});
                    },
                    dojo.hitch( d, 'reject' )
                );
             },
             dojo.hitch( d, 'reject' )
            );

            return d;
        }.call();
    },

    _getFeatures: function( query, featureCallback, finishedCallback, errorCallback ) {
        var thisB = this;
        thisB._loadHeader().then( function() {
            thisB.indexedData.getLines(
                query.ref || thisB.refSeq.name,
                query.start,
                query.end,
                function( line ) {
                    var f = thisB.lineToFeature( line );
                    //console.log(f);
                    featureCallback( f );
                    //return f;
                },
                finishedCallback,
                errorCallback
            );
        }, errorCallback );
    },

    /**
     * Interrogate whether a store has data for a given reference
     * sequence.  Calls the given callback with either true or false.
     *
     * Implemented as a binary interrogation because some stores are
     * smart enough to regularize reference sequence names, while
     * others are not.
     */
    hasRefSeq: function( seqName, callback, errorCallback ) {
        return this.indexedData.index.hasRefSeq( seqName, callback, errorCallback );
    }

});
});