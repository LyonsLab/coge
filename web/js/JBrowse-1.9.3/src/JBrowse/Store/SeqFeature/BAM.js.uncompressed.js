require({cache:{
'JBrowse/Store/SeqFeature/GlobalStatsEstimationMixin':function(){
/**
 * Mixin that adds _estimateGlobalStats method to a store, which
 * samples a section of the features in the store and uses those to
 * esimate the statistics of the whole data set.
 */

define("JBrowse/Store/SeqFeature/GlobalStatsEstimationMixin", [
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
'JBrowse/Store/SeqFeature/BAM/File':function(){
define( [
            'dojo/_base/declare',
            'dojo/_base/array',
            'JBrowse/has',
            'JBrowse/Util',
            'JBrowse/Errors',
            'JBrowse/Store/LRUCache',
            './Util',
            './LazyFeature'
        ],
        function( declare, array, has, Util, Errors, LRUCache, BAMUtil, BAMFeature ) {

var BAM_MAGIC = 21840194;
var BAI_MAGIC = 21578050;

var dlog = function(){ console.error.apply(console, arguments); };

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
    fetchedSize: function() {
        return this.maxv.block + (1<<16) - this.minv.block + 1;
    }
});

var readInt   = BAMUtil.readInt;
var readVirtualOffset = BAMUtil.readVirtualOffset;

var BamFile = declare( null,


/**
 * @lends JBrowse.Store.SeqFeature.BAM.File
 */
{

    /**
     * Low-level BAM file reading code.
     *
     * Adapted by Robert Buels from bam.js in the Dalliance Genome
     * Explorer which is copyright Thomas Down 2006-2010
     * @constructs
     */
    constructor: function( args ) {
        this.store = args.store;
        this.data  = args.data;
        this.bai   = args.bai;

        this.chunkSizeLimit = args.chunkSizeLimit || 5000000;
    },

    init: function( args ) {
        var bam = this;
        var successCallback = args.success || function() {};
        var failCallback = args.failure || function(e) { console.error(e, e.stack); };

        this._readBAI( dojo.hitch( this, function() {
            this._readBAMheader( function() {
                successCallback();
            }, failCallback );
        }), failCallback );
    },

    _readBAI: function( successCallback, failCallback ) {
        // Do we really need to fetch the whole thing? :-(
        this.bai.fetch( dojo.hitch( this, function(header) {
            if (!header) {
                dlog("No data read from BAM index (BAI) file");
                failCallback("No data read from BAM index (BAI) file");
                return;
            }

            if( ! has('typed-arrays') ) {
                dlog('Web browser does not support typed arrays');
                failCallback('Web browser does not support typed arrays');
                return;
            }

            var uncba = new Uint8Array(header);
            if( readInt(uncba, 0) != BAI_MAGIC) {
                dlog('Not a BAI file');
                failCallback('Not a BAI file');
                return;
            }

            var nref = readInt(uncba, 4);

            this.indices = [];

            var p = 8;
            for (var ref = 0; ref < nref; ++ref) {
                var blockStart = p;
                var nbin = readInt(uncba, p); p += 4;
                for (var b = 0; b < nbin; ++b) {
                    var bin = readInt(uncba, p);
                    var nchnk = readInt(uncba, p+4);
                    p += 8;
                    for( var chunkNum = 0; chunkNum < nchnk; chunkNum++ ) {
                        var vo = readVirtualOffset( uncba, p );
                        this._findMinAlignment( vo );
                        p += 16;
                    }
                }
                var nintv = readInt(uncba, p); p += 4;
                // as we're going through the linear index, figure out
                // the smallest virtual offset in the indexes, which
                // tells us where the BAM header ends
                this._findMinAlignment( nintv ? readVirtualOffset(uncba,p) : null );

                p += nintv * 8;
                if( nbin > 0 || nintv > 0 ) {
                    this.indices[ref] = new Uint8Array(header, blockStart, p - blockStart);
                }
            }

            this.empty = ! this.indices.length;

            successCallback( this.indices, this.minAlignmentVO );
        }), failCallback );
    },

    _findMinAlignment: function( candidate ) {
        if( candidate && ( ! this.minAlignmentVO || this.minAlignmentVO.cmp( candidate ) < 0 ) )
            this.minAlignmentVO = candidate;
    },

    _readBAMheader: function( successCallback, failCallback ) {
        var thisB = this;
        // We have the virtual offset of the first alignment
        // in the file.  Cannot completely determine how
        // much of the first part of the file to fetch to get just
        // up to that, since the file is compressed.  Thus, fetch
        // up to the start of the BGZF block that the first
        // alignment is in, plus 64KB, which should get us that whole
        // BGZF block, assuming BGZF blocks are no bigger than 64KB.
        thisB.data.read(
            0,
            thisB.minAlignmentVO ? thisB.minAlignmentVO.block + 65535 : null,
            function(r) {
                var unc = BAMUtil.unbgzf(r);
                var uncba = new Uint8Array(unc);

                if( readInt(uncba, 0) != BAM_MAGIC) {
                    dlog('Not a BAM file');
                    failCallback( 'Not a BAM file' );
                    return;
                }

                var headLen = readInt(uncba, 4);

                thisB._readRefSeqs( headLen+8, 65536*4, successCallback, failCallback );
            },
            failCallback
        );
    },

    _readRefSeqs: function( start, refSeqBytes, successCallback, failCallback ) {
        var thisB = this;
        // have to do another request, because sometimes
        // minAlignment VO is just flat wrong.
        // if headLen is not too big, this will just be in the
        // RemoteBinaryFile cache
        thisB.data.read( 0, start+refSeqBytes,
                         function(r) {
            var unc = BAMUtil.unbgzf(r);
            var uncba = new Uint8Array(unc);

            var nRef = readInt(uncba, start );
            var p = start + 4;

            thisB.chrToIndex = {};
            thisB.indexToChr = [];
            for (var i = 0; i < nRef; ++i) {
                var lName = readInt(uncba, p);
                var name = '';
                for (var j = 0; j < lName-1; ++j) {
                    name += String.fromCharCode(uncba[p + 4 + j]);
                }

                var lRef = readInt(uncba, p + lName + 4);
                //console.log(name + ': ' + lRef);
                thisB.chrToIndex[ thisB.store.browser.regularizeReferenceName( name ) ] = i;
                thisB.indexToChr.push({ name: name, length: lRef });

                p = p + 8 + lName;
                if( p > uncba.length ) {
                    // we've gotten to the end of the data without
                    // finishing reading the ref seqs, need to fetch a
                    // bigger chunk and try again.  :-(
                    refSeqBytes *= 2;
                    console.warn( 'BAM header is very big.  Re-fetching '+refSeqBytes+' bytes.' );
                    thisB._readRefSeqs( start, refSeqBytes, successCallback, failCallback );
                    return;
                }
            }

            successCallback();

        }, failCallback );
    },

    /**
     * Get an array of Chunk objects for the given ref seq id and range.
     */
    blocksForRange: function(refId, min, max) {
        var index = this.indices[refId];
        if (!index) {
            return [];
        }

        // object as { <binNum>: true, ... } containing the bin numbers
        // that overlap this range
        var overlappingBins = function() {
            var intBins = {};
            var intBinsL = this._reg2bins(min, max);
            for (var i = 0; i < intBinsL.length; ++i) {
                intBins[intBinsL[i]] = true;
            }
            return intBins;
        }.call(this);

        // parse the chunks for the overlapping bins out of the index
        // for this ref seq, keeping a distinction between chunks from
        // leaf (lowest-level, smallest) bins, and chunks from other,
        // larger bins
        var leafChunks  = [];
        var otherChunks = [];
        var nbin = readInt(index, 0);
        var p = 4;
        for (var b = 0; b < nbin; ++b) {
            var bin   = readInt(index, p  );
            var nchnk = readInt(index, p+4);
            p += 8;
            if( overlappingBins[bin] ) {
                for (var c = 0; c < nchnk; ++c) {
                    var cs = readVirtualOffset( index, p     );
                    var ce = readVirtualOffset( index, p + 8 );
                    ( bin < 4681 ? otherChunks : leafChunks ).push( new Chunk(cs, ce, bin) );
                    p += 16;
                }
            } else {
                p += nchnk * 16;
            }
        }

        // parse the linear index to find the lowest virtual offset
        var lowest = function() {
            var lowest = null;
            var nintv  = readInt(index, p);
            var minLin = Math.min(min>>14, nintv - 1);
            var maxLin = Math.min(max>>14, nintv - 1);
            for (var i = minLin; i <= maxLin; ++i) {
                var lb =  readVirtualOffset(index, p + 4 + (i * 8));
                if( !lb )
                    continue;

                if ( ! lowest || lb.cmp( lowest ) < 0  )
                    lowest = lb;
            }
            return lowest;
        }();

        // discard any chunks that come before the lowest
        // virtualOffset that we got from the linear index
        otherChunks = function( otherChunks ) {
            var relevantOtherChunks = [];
            if (lowest != null) {
                for (var i = 0; i < otherChunks.length; ++i) {
                    var chnk = otherChunks[i];
                    if (chnk.maxv.block >= lowest.block && chnk.maxv.offset >= lowest.offset) {
                        relevantOtherChunks.push(chnk);
                    }
                }
            }
            return relevantOtherChunks;
        }(otherChunks);

        // add the leaf chunks in, and sort the chunks ascending by virtual offset
        var allChunks = otherChunks
            .concat( leafChunks )
            .sort( function(c0, c1) {
                      return c0.minv.block - c1.minv.block || c0.minv.offset - c1.minv.offset;
                   });

        // merge chunks from the same block together
        var mergedChunks = [];
        if( allChunks.length ) {
            var cur = allChunks[0];
            for (var i = 1; i < allChunks.length; ++i) {
                var nc = allChunks[i];
                if (nc.minv.block == cur.maxv.block /* && nc.minv.offset == cur.maxv.offset */) { // no point splitting mid-block
                    cur = new Chunk(cur.minv, nc.maxv, 'merged');
                } else {
                    mergedChunks.push(cur);
                    cur = nc;
                }
            }
            mergedChunks.push(cur);
        }

        return mergedChunks;
    },

    fetch: function(chr, min, max, featCallback, endCallback, errorCallback ) {

        chr = this.store.browser.regularizeReferenceName( chr );

        var chrId = this.chrToIndex && this.chrToIndex[chr];
        var chunks;
        if( !( chrId >= 0 ) ) {
            chunks = [];
        } else {
            chunks = this.blocksForRange(chrId, min, max);
            if (!chunks) {
                callback(null, new Errors.Fatal('Error in index fetch') );
            }
        }

        // toString function is used by the cache for making cache keys
        chunks.toString = function() {
            return this.join(', ');
        };

        //console.log( chr, min, max, chunks.toString() );

        try {
            this._fetchChunkFeatures(
                chunks,
                chrId,
                min,
                max,
                featCallback,
                endCallback,
                errorCallback
            );
        } catch( e ) {
            errorCallback( e );
        }
    },

    _fetchChunkFeatures: function( chunks, chrId, min, max, featCallback, endCallback, errorCallback ) {
        var thisB = this;

        if( ! chunks.length ) {
            endCallback();
            return;
        }

        var chunksProcessed = 0;

        var cache = this.featureCache = this.featureCache || new LRUCache({
            name: 'bamFeatureCache',
            fillCallback: dojo.hitch( this, '_readChunk' ),
            sizeFunction: function( features ) {
                return features.length;
            },
            maxSize: 100000 // cache up to 100,000 BAM features
        });

        // check the chunks for any that are over the size limit.  if
        // any are, don't fetch any of them
        for( var i = 0; i<chunks.length; i++ ) {
            var size = chunks[i].fetchedSize();
            if( size > this.chunkSizeLimit ) {
                errorCallback( new Errors.DataOverflow('Too many BAM features. BAM chunk size '+Util.commifyNumber(size)+' bytes exceeds chunkSizeLimit of '+Util.commifyNumber(this.chunkSizeLimit)+'.' ) );
                return;
            }
        }

        var haveError;
        var pastStart;
        array.forEach( chunks, function( c ) {
            cache.get( c, function( f, e ) {
                if( e && !haveError )
                    errorCallback(e);
                if(( haveError = haveError || e )) {
                    return;
                }

                for( var i = 0; i<f.length; i++ ) {
                    var feature = f[i];
                    if( feature._refID == chrId ) {
                        // on the right ref seq
                        if( feature.get('start') > max ) // past end of range, can stop iterating
                            break;
                        else if( feature.get('end') >= min ) // must be in range
                            featCallback( feature );
                    }
                }
                if( ++chunksProcessed == chunks.length ) {
                    endCallback();
                }
            });
        });

    },

    _readChunk: function( chunk, callback ) {
        var thisB = this;
        var features = [];
        // console.log('chunk '+chunk+' size ',Util.humanReadableNumber(size));

        thisB.data.read( chunk.minv.block, chunk.fetchedSize(), function(r) {
            try {
                var data = BAMUtil.unbgzf(r, chunk.maxv.block - chunk.minv.block + 1);
                thisB.readBamFeatures( new Uint8Array(data), chunk.minv.offset, features, callback );
            } catch( e ) {
                callback( null, new Errors.Fatal(e) );
            }
        }, function( e ) {
            callback( null, new Errors.Fatal(e) );
        });
    },

    readBamFeatures: function(ba, blockStart, sink, callback ) {
        var that = this;
        var featureCount = 0;

        var maxFeaturesWithoutYielding = 300;

        while ( true ) {
            if( blockStart >= ba.length ) {
                // if we're done, call the callback and return
                callback( sink );
                return;
            }
            else if( featureCount <= maxFeaturesWithoutYielding ) {
                // if we've read no more than 200 features this cycle, read another one
                var blockSize = readInt(ba, blockStart);
                var blockEnd = blockStart + 4 + blockSize - 1;

                // only try to read the feature if we have all the bytes for it
                if( blockEnd < ba.length ) {
                    var feature = new BAMFeature({
                        store: this.store,
                        file: this,
                        bytes: { byteArray: ba, start: blockStart, end: blockEnd }
                     });
                    sink.push(feature);
                    featureCount++;
                }

                blockStart = blockEnd+1;
            }
            else {
                // if we're not done but we've read a good chunk of
                // features, put the rest of our work into a timeout to continue
                // later, avoiding blocking any UI stuff that's going on
                window.setTimeout( function() {
                    that.readBamFeatures( ba, blockStart, sink, callback );
                }, 1);
                return;
            }
        }
    },
    //
    // Binning (transliterated from SAM1.3 spec)
    //

    /* calculate bin given an alignment covering [beg,end) (zero-based, half-close-half-open) */
    _reg2bin: function( beg, end ) {
        --end;
        if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
        if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
        if (beg>>20 == end>>20) return ((1<<9)-1)/7 + (beg>>20);
        if (beg>>23 == end>>23) return ((1<<6)-1)/7 + (beg>>23);
        if (beg>>26 == end>>26) return ((1<<3)-1)/7 + (beg>>26);
        return 0;
    },

    /* calculate the list of bins that may overlap with region [beg,end) (zero-based) */
    MAX_BIN: (((1<<18)-1)/7),
    _reg2bins: function( beg, end ) {
        var k, list = [ 0 ];
        --end;
        for (k = 1    + (beg>>26); k <= 1    + (end>>26); ++k) list.push(k);
        for (k = 9    + (beg>>23); k <= 9    + (end>>23); ++k) list.push(k);
        for (k = 73   + (beg>>20); k <= 73   + (end>>20); ++k) list.push(k);
        for (k = 585  + (beg>>17); k <= 585  + (end>>17); ++k) list.push(k);
        for (k = 4681 + (beg>>14); k <= 4681 + (end>>14); ++k) list.push(k);
        return list;
    }

});

return BamFile;

});

},
'JBrowse/Store/SeqFeature/BAM/Util':function(){
define( [ 'jszlib/inflate',
          'jszlib/arrayCopy',
          'JBrowse/Util'
        ],
        function( inflate, arrayCopy, Util ) {

var VirtualOffset = Util.fastDeclare({
    constructor: function(b, o) {
        this.block = b;
        this.offset = o;
    },
    toString: function() {
        return '' + this.block + ':' + this.offset;
    },
    cmp: function(b) {
        var a = this;
        return b.block - a.block || b.offset - a.offset;
    }
});

/**
 * @lends JBrowse.Store.SeqFeature.BAM.Util
 * Package of utility functions used in various places in the BAM code.
 */
var Utils = {

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

    readVirtualOffset: function(ba, offset) {
        //console.log( 'readVob', offset );
        var block = (ba[offset+6] & 0xff) * 0x100000000
            + (ba[offset+5] & 0xff) * 0x1000000
            + (ba[offset+4] & 0xff) * 0x10000
            + (ba[offset+3] & 0xff) * 0x100
            + (ba[offset+2] & 0xff);
        var bint = (ba[offset+1] << 8) | ba[offset];
        if (block == 0 && bint == 0) {
            return null;  // Should only happen in the linear index?
        } else {
            return new VirtualOffset(block, bint);
        }
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

            var xlen = Utils.readShort( ba, 10 );
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
};

return Utils;

});
},
'JBrowse/Store/SeqFeature/BAM/LazyFeature':function(){
define( ['dojo/_base/array',
         'JBrowse/Util',
         './Util',
         'JBrowse/Model/SimpleFeature'
        ],
        function( array, Util, BAMUtil, SimpleFeature ) {

var SEQRET_DECODER = ['=', 'A', 'C', 'x', 'G', 'x', 'x', 'x', 'T', 'x', 'x', 'x', 'x', 'x', 'x', 'N'];
var CIGAR_DECODER  = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', '?', '?', '?', '?', '?', '?', '?'];

var readInt   = BAMUtil.readInt;
var readShort = BAMUtil.readShort;
var readFloat = BAMUtil.readFloat;

var Feature = Util.fastDeclare(
{
    constructor: function( args ) {
        this.store = args.store;
        this.file  = args.file;
        this.data  = {
            type: 'match',
            source: args.store.source
        };
        this.bytes = {
            start: args.bytes.start,
            end: args.bytes.end,
            byteArray: args.bytes.byteArray
        };

        this._coreParse();
    },

    get: function( field) {
        return this._get( field.toLowerCase() );
    },

    // same as get(), except requires lower-case arguments.  used
    // internally to save lots of calls to field.toLowerCase()
    _get: function( field ) {
        return field in this.data ? this.data[field] : // have we already parsed it out?
            function(field) {
                var v = this.data[field] =
                    this[field]            ? this[field]()            : // maybe we have a special parser for it
                    this._flagMasks[field] ? this._parseFlag( field ) : // or is it a flag?
                                             this._parseTag( field );   // otherwise, look for it in the tags
                return v;
            }.call(this,field);
    },

    tags: function() {
        return this._get('_tags');
    },

    _tags: function() {
        this._parseAllTags();

        var tags = [ 'seq', 'seq_reverse_complemented', 'unmapped' ];
        if( ! this._get('unmapped') )
            tags.push( 'start', 'end', 'strand', 'score', 'qual', 'MQ', 'CIGAR', 'length_on_ref' );
        if( this._get('multi_segment_template') ) {
            tags.push( 'multi_segment_all_aligned',
                       'multi_segment_next_segment_unmapped',
                       'multi_segment_next_segment_reversed',
                       'multi_segment_first',
                       'multi_segment_last',
                       'secondary_alignment',
                       'qc_failed',
                       'duplicate',
                       'next_segment_position'
                     );
        }
        tags = tags.concat( this._tagList || [] );

        var d = this.data;
        for( var k in d ) {
            if( d.hasOwnProperty( k ) && k[0] != '_' )
                tags.push( k );
        }

        var seen = {};
        tags = array.filter( tags, function(t) {
            if( t in this.data && this.data[t] === undefined )
                return false;

            var lt = t.toLowerCase();
            var s = seen[lt];
            seen[lt] = true;
            return ! s;
        },this);

        return tags;
    },

    parent: function() {
        return undefined;
    },

    children: function() {
        return this._get('subfeatures');
    },

    id: function() {
        return this._get('name')+'/'+this._get('md')+'/'+this._get('cigar')+'/'+this._get('start');
    },

    // special parsers
    /**
     * Mapping quality score.
     */
    mq: function() {
        var mq = (this._get('_bin_mq_nl') & 0xff00) >> 8;
        return mq == 255 ? undefined : mq;
    },
    score: function() {
        return this._get('mq');
    },
    qual: function() {
        if( this._get('unmapped') )
            return undefined;

        var qseq = [];
        var byteArray = this.bytes.byteArray;
        var p = this.bytes.start + 36 + this._get('_l_read_name') + this._get('_n_cigar_op')*4 + this._get('_seq_bytes');
        var lseq = this._get('seq_length');
        for (var j = 0; j < lseq; ++j) {
            qseq.push( byteArray[p + j] );
        }
        return qseq.join(' ');
    },
    strand: function() {
        var xs = this._get('xs');
        return xs ? ( xs == '-' ? -1 : 1 ) :
               this._get('seq_reverse_complemented') ? -1 :  1;
    },
    /**
     * Length in characters of the read name.
     */
    _l_read_name: function() {
        return this._get('_bin_mq_nl') & 0xff;
    },
    /**
     * number of bytes in the sequence field
     */
    _seq_bytes: function() {
        return (this._get('seq_length') + 1) >> 1;
    },
    seq: function() {
        var seq = '';
        var byteArray = this.bytes.byteArray;
        var p = this.bytes.start + 36 + this._get('_l_read_name') + this._get('_n_cigar_op')*4;
        var seqBytes = this._get('_seq_bytes');
        for (var j = 0; j < seqBytes; ++j) {
            var sb = byteArray[p + j];
            seq += SEQRET_DECODER[(sb & 0xf0) >> 4];
            seq += SEQRET_DECODER[(sb & 0x0f)];
        }
        return seq;
    },
    name: function() {
        return this._get('_read_name');
    },
    _read_name: function() {
        var byteArray = this.bytes.byteArray;
        var readName = '';
        var nl = this._get('_l_read_name');
        var p = this.bytes.start + 36;
        for (var j = 0; j < nl-1; ++j) {
            readName += String.fromCharCode(byteArray[p+j]);
        }
        return readName;
    },
    _n_cigar_op: function() {
        return this._get('_flag_nc') & 0xffff;
    },
    cigar: function() {
        if( this._get('unmapped') )
            return undefined;

        var byteArray   = this.bytes.byteArray;
        var numCigarOps = this._get('_n_cigar_op');
        var p = this.bytes.start + 36 + this._get('_l_read_name');
        var cigar = '';
        var lref = 0;
        for (var c = 0; c < numCigarOps; ++c) {
            var cigop = readInt(byteArray, p);
            var lop = cigop >> 4;
            var op = CIGAR_DECODER[cigop & 0xf];
            cigar += lop + op;

            // soft clip, hard clip, and insertion don't count toward
            // the length on the reference
            if( op != 'H' && op != 'S' && op != 'I' )
                lref += lop;

            p += 4;
        }

        this.data.length_on_ref = lref;
        return cigar;
    },
    next_segment_position: function() {
        var nextRefID = this._get('_next_refid');
        var nextSegment = this.file.indexToChr[nextRefID];
        if( nextSegment )
            return nextSegment.name+':'+this._get('_next_pos');
        else
            return undefined;
    },
    subfeatures: function() {
        if( ! this.store.createSubfeatures )
            return undefined;

        var cigar = this._get('cigar');
        if( cigar )
            return this._cigarToSubfeats( cigar );

        return undefined;
    },
    length_on_ref: function() {
        var c = this._get('cigar'); // the length_on_ref is set as a
                                   // side effect of the CIGAR parsing
        return this.data.length_on_ref;
    },
    _flags: function() {
        return (this.get('_flag_nc') & 0xffff0000) >> 16;
    },
    end: function() {
        return this._get('start') + ( this._get('length_on_ref') || this._get('seq_length') || undefined );
    },

    seq_id: function() {
        if( this._get('unmapped') )
            return undefined;

        return ( this.file.indexToChr[ this._refID ] || {} ).name;
    },

    _bin_mq_nl: function() {
        with( this.bytes )
            return readInt( byteArray, start + 12  );
    },
    _flag_nc: function() {
        with( this.bytes )
            return readInt( byteArray, start + 16 );
    },
    seq_length: function() {
        with( this.bytes )
            return readInt( byteArray, start + 20 );
    },
    _next_refid: function() {
        with( this.bytes )
            return readInt( byteArray, start + 24 );
    },
    _next_pos: function() {
        with( this.bytes )
            return readInt( byteArray, start + 28 );
    },
    template_length: function() {
        with( this.bytes )
            return readInt( byteArray, start + 32 );
    },

    /**
     * parse the core data: ref ID and start
     */
    _coreParse: function() {
        with( this.bytes ) {
            this._refID      = readInt( byteArray, start + 4 );
            this.data.start  = readInt( byteArray, start + 8 );
        }
    },

    /**
     * Get the value of a tag, parsing the tags as far as necessary.
     * Only called if we have not already parsed that field.
     */
    _parseTag: function( tagName ) {
        // if all of the tags have been parsed and we're still being
        // called, we already know that we have no such tag, because
        // it would already have been cached.
        if( this._allTagsParsed )
            return undefined;

        this._tagList = this._tagList || [];
        var byteArray = this.bytes.byteArray;
        var p = this._tagOffset || this.bytes.start + 36 + this._get('_l_read_name') + this._get('_n_cigar_op')*4 + this._get('_seq_bytes') + this._get('seq_length');

        var blockEnd = this.bytes.end;
        while( p < blockEnd && lcTag != tagName ) {
            var tag      = String.fromCharCode( byteArray[p], byteArray[ p+1 ] );
            var lcTag    = tag.toLowerCase();
            var type = String.fromCharCode( byteArray[ p+2 ] );
            p += 3;

            var value;
            switch( type.toLowerCase() ) {
            case 'a':
                value = String.fromCharCode( byteArray[p] );
                p += 1;
                break;
            case 'i':
                value = readInt(byteArray, p );
                p += 4;
                break;
            case 'c':
                value = byteArray[p];
                p += 1;
                break;
            case 's':
                value = readShort(byteArray, p);
                p += 2;
                break;
            case 'f':
                value = readFloat( byteArray, p );
                p += 4;
                break;
            case 'z':
            case 'h':
                value = '';
                while( p <= blockEnd ) {
                    var cc = byteArray[p++];
                    if( cc == 0 ) {
                        break;
                    }
                    else {
                        value += String.fromCharCode(cc);
                    }
                }
                break;
            default:
                console.warn( "Unknown BAM tag type '"+type
                              +"', tags may be incomplete"
                            );
                value = undefined;
                p = blockEnd; // stop parsing tags
            }

            this._tagOffset = p;

            this._tagList.push( tag );
            if( lcTag == tagName )
                return value;
            else {
                this.data[ lcTag ] = value;
            }
        }
        this._allTagsParsed = true;
        return undefined;
    },
    _parseAllTags: function() {
        this._parseTag(); // calling _parseTag with no arg just parses
        // all the tags and returns the last one
    },

    _flagMasks: {
        multi_segment_template:              0x1,
        multi_segment_all_aligned:           0x2,
        unmapped:                            0x4,
        multi_segment_next_segment_unmapped: 0x8,
        seq_reverse_complemented:            0x10,
        multi_segment_next_segment_reversed: 0x20,
        multi_segment_first:                 0x40,
        multi_segment_last:                  0x80,
        secondary_alignment:                 0x100,
        qc_failed:                           0x200,
        duplicate:                           0x400
    },

    _parseFlag: function( flagName ) {
        return !!( this._get('_flags') & this._flagMasks[flagName] );
    },

    _parseCigar: function( cigar ) {
        return array.map( cigar.match(/\d+\D/g), function( op ) {
           return [ op.match(/\D/)[0].toUpperCase(), parseInt( op ) ];
        });
    },

    /**
     *  take a cigar string, and initial position, return an array of subfeatures
     */
    _cigarToSubfeats: function(cigar)    {
        var subfeats = [];
        var min = this._get('start');
        var max;
        var ops = this._parseCigar( cigar );
        for (var i = 0; i < ops.length; i++)  {
            var lop = ops[i][1];
            var op = ops[i][0];  // operation type
            // converting "=" to "E" to avoid possible problems later with non-alphanumeric type name
            if (op === "=")  { op = "E"; }

            switch (op) {
            case 'M':
            case 'D':
            case 'N':
            case 'E':
            case 'X':
                max = min + lop;
                break;
            case 'I':
                max = min;
                break;
            case 'P':  // not showing padding deletions (possibly change this later -- could treat same as 'I' ?? )
            case 'H':  // not showing hard clipping (since it's unaligned, and offset arg meant to be beginning of aligned part)
            case 'S':  // not showing soft clipping (since it's unaligned, and offset arg meant to be beginning of aligned part)
                break;
                // other possible cases
            }
            if( op !== 'N' ) {
                var subfeat = new SimpleFeature({
                    data: {
                    type: op,
                        start: min,
                        end: max,
                        strand: this._get('strand'),
                        cigar_op: lop+op
                    },
                    parent: this
                });
                subfeats.push(subfeat);
            }
            min = max;
        }
        return subfeats;
    }

});

return Feature;
});
},
'JBrowse/Model/SimpleFeature':function(){
/**
 * Simple implementation of a feature object.
 */
define([
        'JBrowse/Util'
       ],
       function( Util ) {

var counter = 0;

var SimpleFeature = Util.fastDeclare({

    /**
     * @param args.data {Object} key-value data, must include 'start' and 'end'
     * @param args.parent {Feature} optional parent feature
     * @param args.id {String} optional unique identifier.  can also be in data.uniqueID.
     *
     * Note: args.data.subfeatures can be an array of these same args,
     * which will be inflated to more instances of this class.
     */
    constructor: function( args ) {
        args = args || {};
        this.data = args.data || {};
        this._parent = args.parent;
        this._uniqueID = args.id || this.data.uniqueID || (
            this._parent ? this._parent.id()+'_'+(counter++) : 'SimpleFeature_'+(counter++)
        );

        // inflate any subfeatures that are not already feature objects
        var subfeatures;
        if(( subfeatures = this.data.subfeatures )) {
            for( var i = 0; i < subfeatures.length; i++ ) {
                if( typeof subfeatures[i].get != 'function' ) {
                    subfeatures[i] = new SimpleFeature(
                        { data: subfeatures[i],
                          parent: this
                        });
                }
            }
        }
    },

    /**
     * Get a piece of data about the feature.  All features must have
     * 'start' and 'end', but everything else is optional.
     */
    get: function(name) {
        return this.data[ name ];
    },

    /**
     * Set an item of data.
     */
    set: function( name, val ) {
        this.data[ name ] = val;
    },

    /**
     * Get an array listing which data keys are present in this feature.
     */
    tags: function() {
        var t = [];
        var d = this.data;
        for( var k in d ) {
            if( d.hasOwnProperty( k ) )
                t.push( k );
        }
        return t;
    },

    /**
     * Get the unique ID of this feature.
     */
    id: function( newid ) {
        if( newid )
            this._uniqueID = newid;
        return this._uniqueID;
    },

    /**
     * Get this feature's parent feature, or undefined if none.
     */
    parent: function() {
        return this._parent;
    },

    /**
     * Get an array of child features, or undefined if none.
     */
    children: function() {
        return this.get('subfeatures');
    }

});

return SimpleFeature;
});
}}});
define( "JBrowse/Store/SeqFeature/BAM", [
            'dojo/_base/declare',
            'dojo/_base/array',
            'dojo/_base/Deferred',
            'dojo/_base/lang',
            'JBrowse/has',
            'JBrowse/Util',
            'JBrowse/Store/SeqFeature',
            'JBrowse/Store/DeferredStatsMixin',
            'JBrowse/Store/DeferredFeaturesMixin',
            'JBrowse/Model/XHRBlob',
            'JBrowse/Store/SeqFeature/GlobalStatsEstimationMixin',
            './BAM/File'
        ],
        function(
            declare,
            array,
            Deferred,
            lang,
            has,
            Util,
            SeqFeatureStore,
            DeferredStatsMixin,
            DeferredFeaturesMixin,
            XHRBlob,
            GlobalStatsEstimationMixin,
            BAMFile
        ) {

var BAMStore = declare( [ SeqFeatureStore, DeferredStatsMixin, DeferredFeaturesMixin, GlobalStatsEstimationMixin ],

/**
 * @lends JBrowse.Store.SeqFeature.BAM
 */
{
    /**
     * Data backend for reading feature data directly from a
     * web-accessible BAM file.
     *
     * @constructs
     */
    constructor: function( args ) {

        this.createSubfeatures = args.subfeatures;

        var bamBlob = args.bam ||
            new XHRBlob( this.resolveUrl(
                             args.urlTemplate || 'data.bam'
                         )
                       );

        var baiBlob = args.bai ||
            new XHRBlob( this.resolveUrl(
                             args.baiUrlTemplate || ( args.urlTemplate ? args.urlTemplate+'.bai' : 'data.bam.bai' )
                         )
                       );

        this.bam = new BAMFile({
                store: this,
                data: bamBlob,
                bai: baiBlob,
                chunkSizeLimit: args.chunkSizeLimit
        });

        this.source = ( bamBlob.url  ? bamBlob.url.match( /\/([^/\#\?]+)($|[\#\?])/ )[1] :
                        bamBlob.blob ? bamBlob.blob.name : undefined ) || undefined;

        if( ! has( 'typed-arrays' ) ) {
            this._failAllDeferred( 'Web browser does not support typed arrays');
            return;
        }

        this.bam.init({
            success: dojo.hitch( this, '_estimateGlobalStats',
                                 dojo.hitch( this, function( stats, error ) {
                                     if( error )
                                         this._failAllDeferred( error );
                                     else {
                                         this.globalStats = stats;
                                         this._deferred.stats.resolve({success:true});
                                         this._deferred.features.resolve({success:true});
                                     }

                                 }),
                                 dojo.hitch( this, '_failAllDeferred' )
                               ),
            failure: dojo.hitch( this, '_failAllDeferred' )
        });
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
        this._deferred.stats.then( function() {
            callback( seqName in thisB.bam.chrToIndex );
        }, errorCallback );
    },

    // called by getFeatures from the DeferredFeaturesMixin
    _getFeatures: function( query, featCallback, endCallback, errorCallback ) {
        this.bam.fetch( this.refSeq.name, query.start, query.end, featCallback, endCallback, errorCallback );
    }

});

return BAMStore;
});