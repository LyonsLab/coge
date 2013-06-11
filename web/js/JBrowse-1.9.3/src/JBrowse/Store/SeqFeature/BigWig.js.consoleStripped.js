require({cache:{
'JBrowse/Store/SeqFeature/BigWig/Window':function(){
define( [
            'dojo/_base/declare',
            'dojo/_base/lang',
            './RequestWorker'
        ],
        function( declare, lang, RequestWorker ) {

var dlog = function(){ console.log.apply(console, arguments); };

return declare( null,
 /**
  * @lends JBrowse.Store.BigWig.Window.prototype
  */
{

    /**
     * View into a subset of the data in a BigWig file.
     *
     * Adapted by Robert Buels from bigwig.js in the Dalliance Genome
     * Explorer by Thomas Down.
     * @constructs
     */
    constructor: function(bwg, cirTreeOffset, cirTreeLength, isSummary) {
        this.bwg = bwg;
        this.cirTreeOffset = cirTreeOffset;
        this.cirTreeLength = cirTreeLength;
        this.isSummary = isSummary;
    },

    BED_COLOR_REGEXP: /^[0-9]+,[0-9]+,[0-9]+/,

    readWigData: function(chrName, min, max, callback, errorCallback ) {
        // 0 && console.log( 'reading wig data from '+chrName+':'+min+'..'+max);
        var chr = this.bwg.refsByName[chrName];
        if ( ! chr ) {
            // Not an error because some .bwgs won't have data for all chromosomes.

            // dlog("Couldn't find chr " + chrName);
            // dlog('Chroms=' + miniJSONify(this.bwg.refsByName));
            callback([]);
        } else {
            this.readWigDataById( chr.id, min, max, callback, errorCallback );
        }
    },

    readWigDataById: function(chr, min, max, callback, errorCallback ) {
        if( !this.cirHeader ) {
            var readCallback = lang.hitch( this, 'readWigDataById', chr, min, max, callback );
            if( this.cirHeaderLoading ) {
                this.cirHeaderLoading.push( readCallback );
            }
            else {
                this.cirHeaderLoading = [ readCallback ];
                // dlog('No CIR yet, fetching');
                this.bwg.data
                    .slice(this.cirTreeOffset, 48)
                    .fetch( lang.hitch( this, function(result) {
                                this.cirHeader = result;
                                var la = new Int32Array( this.cirHeader, 0, 2 );
                                this.cirBlockSize = la[1];
                                dojo.forEach( this.cirHeaderLoading, function(c) { c(); });
                                delete this.cirHeaderLoading;
                            }), errorCallback );
            }
            return;
        }

        //dlog('_readWigDataById', chr, min, max, callback);

        var worker = new RequestWorker( this, chr, min, max, callback, errorCallback );
        worker.cirFobRecur([this.cirTreeOffset + 48], 1);
    }
});

});

},
'JBrowse/Store/SeqFeature/BigWig/RequestWorker':function(){
define( [
            'dojo/_base/declare',
            'dojo/_base/lang',
            'JBrowse/Model/Range',
            'jszlib/inflate'
        ],
        function( declare, dlang, Range, inflate ) {

var dlog = function(){ console.log.apply(console, arguments); };

var gettable = declare( null, {
    get: function(name) {
        return this[ { start: 'min', end: 'max', seq_id: 'segment' }[name] || name ];
    },
    tags: function() {
        return ['start','end','seq_id','score','type','source'];
    }
});
var Feature = declare( gettable, {} );
var Group = declare( gettable, {} );

var RequestWorker = declare( null,
 /**
  * @lends JBrowse.Store.BigWig.Window.RequestWorker.prototype
  */
 {

    BIG_WIG_TYPE_GRAPH: 1,
    BIG_WIG_TYPE_VSTEP: 2,
    BIG_WIG_TYPE_FSTEP: 3,

    /**
     * Worker object for reading data from a bigwig or bigbed file.
     * Manages the state necessary for traversing the index trees and
     * so forth.
     *
     * Adapted by Robert Buels from bigwig.js in the Dalliance Genome
     * Explorer by Thomas Down.
     * @constructs
     */
    constructor: function( window, chr, min, max, callback, errorCallback ) {
        this.window = window;
        this.source = window.bwg.name || undefined;

        this.blocksToFetch = [];
        this.outstanding = 0;

        this.chr = chr;
        this.min = min;
        this.max = max;
        this.callback = callback;
        this.errorCallback = errorCallback || function(e) { console.error( e, e.stack, arguments.caller ); };
    },

    cirFobRecur: function(offset, level) {
        this.outstanding += offset.length;

        var maxCirBlockSpan = 4 +  (this.window.cirBlockSize * 32);   // Upper bound on size, based on a completely full leaf node.
        var spans;
        for (var i = 0; i < offset.length; ++i) {
            var blockSpan = new Range(offset[i], Math.min(offset[i] + maxCirBlockSpan, this.window.cirTreeOffset + this.window.cirTreeLength));
            spans = spans ? spans.union( blockSpan ) : blockSpan;
        }

        var fetchRanges = spans.ranges();
        //dlog('fetchRanges: ' + fetchRanges);
        for (var r = 0; r < fetchRanges.length; ++r) {
            var fr = fetchRanges[r];
            this.cirFobStartFetch(offset, fr, level);
        }
    },

    cirFobStartFetch: function(offset, fr, level, attempts) {
        var length = fr.max() - fr.min();
        //dlog('fetching ' + fr.min() + '-' + fr.max() + ' (' + (fr.max() - fr.min()) + ')');
        //0 && console.log('cirfobstartfetch');
        this.window.bwg.data
            .slice(fr.min(), fr.max() - fr.min())
            .fetch( dlang.hitch( this,function(resultBuffer) {
                                     for (var i = 0; i < offset.length; ++i) {
                                         if (fr.contains(offset[i])) {
                                             this.cirFobRecur2(resultBuffer, offset[i] - fr.min(), level);
                                             --this.outstanding;
                                             if (this.outstanding == 0) {
                                                 this.cirCompleted();
                                             }
                                         }
                                     }
                                 }), this.errorCallback );
    },

    cirFobRecur2: function(cirBlockData, offset, level) {
        var ba = new Int8Array(cirBlockData);
        var sa = new Int16Array(cirBlockData);
        var la = new Int32Array(cirBlockData);

        var isLeaf = ba[offset];
        var cnt = sa[offset/2 + 1];
        // dlog('cir level=' + level + '; cnt=' + cnt);
        offset += 4;

        if (isLeaf != 0) {
            for (var i = 0; i < cnt; ++i) {
                var lo = offset/4;
                var startChrom = la[lo];
                var startBase = la[lo + 1];
                var endChrom = la[lo + 2];
                var endBase = la[lo + 3];
                var blockOffset = (la[lo + 4]<<32) | (la[lo + 5]);
                var blockSize = (la[lo + 6]<<32) | (la[lo + 7]);
                if ((startChrom < this.chr || (startChrom == this.chr && startBase <= this.max)) &&
                    (endChrom   > this.chr || (endChrom == this.chr && endBase >= this.min)))
                {
                    // dlog('Got an interesting block: startBase=' + startBase + '; endBase=' + endBase + '; offset=' + blockOffset + '; size=' + blockSize);
                    this.blocksToFetch.push({offset: blockOffset, size: blockSize});
                }
                offset += 32;
            }
        } else {
            var recurOffsets = [];
            for (var i = 0; i < cnt; ++i) {
                var lo = offset/4;
                var startChrom = la[lo];
                var startBase = la[lo + 1];
                var endChrom = la[lo + 2];
                var endBase = la[lo + 3];
                var blockOffset = (la[lo + 4]<<32) | (la[lo + 5]);
                if ((startChrom < this.chr || (startChrom == this.chr && startBase <= this.max)) &&
                    (endChrom   > this.chr || (endChrom == this.chr && endBase >= this.min)))
                {
                    recurOffsets.push(blockOffset);
                }
                offset += 24;
            }
            if (recurOffsets.length > 0) {
                this.cirFobRecur(recurOffsets, level + 1);
            }
        }
    },

    cirCompleted: function() {
        this.blocksToFetch.sort(function(b0, b1) {
                               return (b0.offset|0) - (b1.offset|0);
                           });

        if (this.blocksToFetch.length == 0) {
            this.callback([]);
        } else {
            this.features = [];
            this.tramp();
        }
    },

    createFeature: function(fmin, fmax, opts) {
        // dlog('createFeature(' + fmin +', ' + fmax + ')');

        if (!opts) {
            opts = {};
        }

        var f = new Feature();
        f.segment = (this.window.bwg.refsByNumber[this.chr]||{}).name;
        f.min = fmin;
        f.max = fmax;
        f.type = 'remark';
        f.source = this.source;

        for (k in opts) {
            f[k] = opts[k];
        }

        this.features.push(f);
    },
    maybeCreateFeature: function(fmin, fmax, opts) {
        if (fmin <= this.max && fmax >= this.min) {
            this.createFeature( fmin, fmax, opts );
        }
    },
    tramp: function() {
        if (this.blocksToFetch.length == 0) {
            //var afterBWG = new Date();
            // dlog('BWG fetch took ' + (afterBWG - beforeBWG) + 'ms');
            this.callback( this.features );
            return;  // just in case...
        } else {
            var block = this.blocksToFetch[0];
            if (block.data) {
                var ba = new Uint8Array(block.data);

                if (this.window.isSummary) {
                    var sa = new Int16Array(block.data);
                    var la = new Int32Array(block.data);
                    var fa = new Float32Array(block.data);

                    var itemCount = block.data.byteLength/32;
                    for (var i = 0; i < itemCount; ++i) {
                        var chromId =   la[(i*8)];
                        var start =     la[(i*8)+1];
                        var end =       la[(i*8)+2];
                        var validCnt =  la[(i*8)+3];
                        var minVal    = fa[(i*8)+4];
                        var maxVal    = fa[(i*8)+5];
                        var sumData   = fa[(i*8)+6];
                        var sumSqData = fa[(i*8)+7];

                        if (chromId == this.chr) {
                            var summaryOpts = {score: sumData/validCnt};
                            if (this.window.bwg.type == 'bigbed') {
                                summaryOpts.type = 'density';
                            }
                            this.maybeCreateFeature( start, end, summaryOpts);
                        }
                    }
                } else if (this.window.bwg.type == 'bigwig') {
                    var sa = new Int16Array(block.data);
                    var la = new Int32Array(block.data);
                    var fa = new Float32Array(block.data);

                    var chromId = la[0];
                    var blockStart = la[1];
                    var blockEnd = la[2];
                    var itemStep = la[3];
                    var itemSpan = la[4];
                    var blockType = ba[20];
                    var itemCount = sa[11];

                    // dlog('processing bigwig block, type=' + blockType + '; count=' + itemCount);

                    if (blockType == this.BIG_WIG_TYPE_FSTEP) {
                        for (var i = 0; i < itemCount; ++i) {
                            var score = fa[i + 6];
                            this.maybeCreateFeature( blockStart + (i*itemStep), blockStart + (i*itemStep) + itemSpan, {score: score});
                        }
                    } else if (blockType == this.BIG_WIG_TYPE_VSTEP) {
                        for (var i = 0; i < itemCount; ++i) {
                            var start = la[(i*2) + 6];
                            var score = fa[(i*2) + 7];
                            this.maybeCreateFeature( start, start + itemSpan, {score: score});
                        }
                    } else if (blockType == this.BIG_WIG_TYPE_GRAPH) {
                        for (var i = 0; i < itemCount; ++i) {
                            var start = la[(i*3) + 6];
                            var end   = la[(i*3) + 7];
                            var score = fa[(i*3) + 8];
                            if (start > end) {
                                start = end;
                            }
                            this.maybeCreateFeature( start, end, {score: score});
                        }
                    } else {
                        dlog('Currently not handling bwgType=' + blockType);
                    }
                } else if (this.window.bwg.type == 'bigbed') {
                    var offset = 0;
                    while (offset < ba.length) {
                        var chromId = (ba[offset+3]<<24) | (ba[offset+2]<<16) | (ba[offset+1]<<8) | (ba[offset+0]);
                        var start = (ba[offset+7]<<24) | (ba[offset+6]<<16) | (ba[offset+5]<<8) | (ba[offset+4]);
                        var end = (ba[offset+11]<<24) | (ba[offset+10]<<16) | (ba[offset+9]<<8) | (ba[offset+8]);
                        offset += 12;
                        var rest = '';
                        while (true) {
                            var ch = ba[offset++];
                            if (ch != 0) {
                                rest += String.fromCharCode(ch);
                            } else {
                                break;
                            }
                        }

                        var featureOpts = {};

                        var bedColumns = rest.split('\t');
                        if (bedColumns.length > 0) {
                            featureOpts.label = bedColumns[0];
                        }
                        if (bedColumns.length > 1) {
                            featureOpts.score = stringToInt(bedColumns[1]);
                        }
                        if (bedColumns.length > 2) {
                            featureOpts.orientation = bedColumns[2];
                        }
                        if (bedColumns.length > 5) {
                            var color = bedColumns[5];
                            if (this.window.BED_COLOR_REGEXP.test(color)) {
                                featureOpts.override_color = 'rgb(' + color + ')';
                            }
                        }

                        if (bedColumns.length < 9) {
                            if (chromId == this.chr) {
                                this.maybeCreateFeature( start, end, featureOpts);
                            }
                        } else if (chromId == this.chr && start <= this.max && end >= this.min) {
                            // Complex-BED?
                            // FIXME this is currently a bit of a hack to do Clever Things with ensGene.bb

                            var thickStart = bedColumns[3]|0;
                            var thickEnd   = bedColumns[4]|0;
                            var blockCount = bedColumns[6]|0;
                            var blockSizes = bedColumns[7].split(',');
                            var blockStarts = bedColumns[8].split(',');

                            featureOpts.type = 'bb-transcript';
                            var grp = new Group();
                            grp.id = bedColumns[0];
                            grp.type = 'bb-transcript';
                            grp.notes = [];
                            featureOpts.groups = [grp];

                            if (bedColumns.length > 10) {
                                var geneId = bedColumns[9];
                                var geneName = bedColumns[10];
                                var gg = new Group();
                                gg.id = geneId;
                                gg.label = geneName;
                                gg.type = 'gene';
                                featureOpts.groups.push(gg);
                            }

                            var spans = null;
                            for (var b = 0; b < blockCount; ++b) {
                                var bmin = (blockStarts[b]|0) + start;
                                var bmax = bmin + (blockSizes[b]|0);
                                var span = new Range(bmin, bmax);
                                if (spans) {
                                    spans = spans.union( span );
                                } else {
                                    spans = span;
                                }
                            }

                            var tsList = spans.ranges();
                            for (var s = 0; s < tsList.length; ++s) {
                                var ts = tsList[s];
                                this.createFeature( ts.min(), ts.max(), featureOpts);
                            }

                            if (thickEnd > thickStart) {
                                var tl = spans.intersection( new Range(thickStart, thickEnd) );
                                if (tl) {
                                    featureOpts.type = 'bb-translation';
                                    var tlList = tl.ranges();
                                    for (var s = 0; s < tlList.length; ++s) {
                                        var ts = tlList[s];
                                        this.createFeature( ts.min(), ts.max(), featureOpts);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    dlog("Don't know what to do with " + this.window.bwg.type);
                }
                this.blocksToFetch.splice(0, 1);
                this.tramp();
            } else {
                var fetchStart = block.offset;
                var fetchSize = block.size;
                var bi = 1;
                while (bi < this.blocksToFetch.length && this.blocksToFetch[bi].offset == (fetchStart + fetchSize)) {
                    fetchSize += this.blocksToFetch[bi].size;
                    ++bi;
                }

                //dlog('other thing');
                this.window.bwg.data
                    .slice(fetchStart, fetchSize)
                    .fetch(dlang.hitch( this, function(result) {
                                            var offset = 0;
                                            var bi = 0;
                                            while (offset < fetchSize) {
                                                var fb = this.blocksToFetch[bi];

                                                var data;
                                                if (this.window.bwg.uncompressBufSize > 0) {
                                                    // var beforeInf = new Date();
                                                    data = inflate(result, offset + 2, fb.size - 2);
                                                    // var afterInf = new Date();
                                                    // dlog('inflate: ' + (afterInf - beforeInf) + 'ms');
                                                } else {
                                                    var tmp = new Uint8Array(fb.size);    // FIXME is this really the best we can do?
                                                    arrayCopy(new Uint8Array(result, offset, fb.size), 0, tmp, 0, fb.size);
                                                    data = tmp.buffer;
                                                }
                                                fb.data = data;

                                                offset += fb.size;
                                                ++bi;
                                            }
                                            this.tramp();
                                        }), this.errorCallback );
            }
        }
    }
});

return RequestWorker;

});
},
'JBrowse/Model/Range':function(){
define( [
            'dojo/_base/declare'
        ],
        function( declare ) {

var Range = declare( null,
/**
 * @lends JBrowse.Model.Range.prototype
 */
{

    /**
     * Adapted from a combination of Range and _Compound in the
     * Dalliance Genome Explorer, (c) Thomas Down 2006-2010.
     */
    constructor: function() {
        this._ranges =
            arguments.length == 2 ? [ { min: arguments[0], max: arguments[1] } ] :
            0 in arguments[0]     ? dojo.clone( arguments[0] )                   :
                                    [ arguments[0] ];
    },

    min: function() {
        return this._ranges[0].min;
    },

    max: function() {
        return this._ranges[this._ranges.length - 1].max;
    },

    contains: function(pos) {
        for (var s = 0; s < this._ranges.length; ++s) {
            var r = this._ranges[s];
            if ( r.min <= pos && r.max >= pos ) {
                return true;
            }
        }
        return false;
    },

    isContiguous: function() {
        return this._ranges.length > 1;
    },

    ranges: function() {
        return this._ranges.map( function(r) {
            return new Range( r.min, r.max );
        });
    },

    toString: function() {
        return this._ranges
            .map(function(r) { return '['+r.min+'-'+r.max+']'; })
            .join(',');
    },

    union: function(s1) {
        var s0 = this;
        var ranges = s0.ranges().concat(s1.ranges()).sort( this.rangeOrder );
        var oranges = [];
        var current = ranges[0];

        for (var i = 1; i < ranges.length; ++i) {
            var nxt = ranges[i];
            if (nxt.min() > (current.max() + 1)) {
                oranges.push(current);
                current = nxt;
            } else {
                if (nxt.max() > current.max()) {
                    current = new Range(current.min(), nxt.max());
                }
            }
        }
        oranges.push(current);

        if (oranges.length == 1) {
            return oranges[0];
        } else {
            return new _Compound(oranges);
        }
    },

    intersection: function( s1 ) {
        var s0 = this;
        var r0 = s0.ranges();
        var r1 = s1.ranges();
        var l0 = r0.length, l1 = r1.length;
        var i0 = 0, i1 = 0;
        var or = [];

        while (i0 < l0 && i1 < l1) {
            var s0 = r0[i0], s1 = r1[i1];
            var lapMin = Math.max(s0.min(), s1.min());
            var lapMax = Math.min(s0.max(), s1.max());
            if (lapMax >= lapMin) {
                or.push(new Range(lapMin, lapMax));
            }
            if (s0.max() > s1.max()) {
                ++i1;
            } else {
                ++i0;
            }
        }

        if (or.length == 0) {
            return null; // FIXME
        } else if (or.length == 1) {
            return or[0];
        } else {
            return new _Compound(or);
        }
    },

    coverage: function() {
        var tot = 0;
        var rl = this.ranges();
        for (var ri = 0; ri < rl.length; ++ri) {
            var r = rl[ri];
            tot += (r.max() - r.min() + 1);
        }
        return tot;
    },

    rangeOrder: function( a, b ) {
        if( arguments.length < 2 ) {
            b = a;
            a = this;
        }

        if (a.min() < b.min()) {
            return -1;
        } else if (a.min() > b.min()) {
            return 1;
        } else if (a.max() < b.max()) {
            return -1;
        } else if (b.max() > a.max()) {
            return 1;
        } else {
            return 0;
        }
    }
});

return Range;
});


}}});
define( "JBrowse/Store/SeqFeature/BigWig", [
            'dojo/_base/declare',
            'dojo/_base/lang',
            'dojo/_base/array',
            'dojo/_base/url',
            'JBrowse/has',
            'JBrowse/Store/SeqFeature',
            'JBrowse/Store/DeferredStatsMixin',
            'JBrowse/Store/DeferredFeaturesMixin',
            './BigWig/Window',
            'JBrowse/Util',
            'JBrowse/Model/XHRBlob'
        ],
        function( declare, lang, array, urlObj, has, SeqFeatureStore, DeferredFeaturesMixin, DeferredStatsMixin, Window, Util, XHRBlob ) {
return declare([ SeqFeatureStore, DeferredFeaturesMixin, DeferredStatsMixin ],

 /**
  * @lends JBrowse.Store.BigWig
  */
{

    BIG_WIG_MAGIC: -2003829722,
    BIG_BED_MAGIC: -2021002517,

    BIG_WIG_TYPE_GRAPH: 1,
    BIG_WIG_TYPE_VSTEP: 2,
    BIG_WIG_TYPE_FSTEP: 3,

    /**
     * Data backend for reading wiggle data from BigWig or BigBed files.
     *
     * Adapted by Robert Buels from bigwig.js in the Dalliance Genome
     * Explorer which is copyright Thomas Down 2006-2010
     * @constructs
     */
    constructor: function( args ) {

        this.data = args.blob ||
            new XHRBlob( this.resolveUrl(
                             args.urlTemplate || 'data.bigwig'
                         )
                       );

        this.name = args.name || ( this.data.url && new urlObj( this.data.url ).path.replace(/^.+\//,'') ) || 'anonymous';

        this._load();
    },

    _getGlobalStats: function( successCallback, errorCallback ) {
        var s = this._globalStats || {};

        // calc mean and standard deviation if necessary
        if( !( 'scoreMean' in s ))
            s.scoreMean = s.basesCovered ? s.scoreSum / s.basesCovered : 0;
        if( !( 'scoreStdDev' in s ))
            s.scoreStdDev = this._calcStdFromSums( s.scoreSum, s.scoreSumSquares, s.basesCovered );

        successCallback( s );
    },

    _load: function() {
        var bwg = this;
        var headerSlice = bwg.data.slice(0, 512);
        headerSlice.fetch( function( result ) {
            if( ! result ) {
                this._failAllDeferred( 'BBI header not readable' );
                return;
            }

            bwg.fileSize = result.fileSize;;
            var header = result;
            var sa = new Int16Array(header);
            var la = new Int32Array(header);
            if (la[0] == bwg.BIG_WIG_MAGIC) {
                bwg.type = 'bigwig';
            } else if (la[0] == bwg.BIG_BED_MAGIC) {
                bwg.type = 'bigbed';
            } else {
                console.error( 'Format '+la[0]+' not supported' );
                bwg._failAllDeferred( 'Format '+la[0]+' not supported' );
                return;
            }
            //        dlog('magic okay');

            bwg.version = sa[2];             // 4
            bwg.numZoomLevels = sa[3];       // 6
            bwg.chromTreeOffset = (la[2] << 32) | (la[3]);     // 8
            bwg.unzoomedDataOffset = (la[4] << 32) | (la[5]);  // 16
            bwg.unzoomedIndexOffset = (la[6] << 32) | (la[7]); // 24
            bwg.fieldCount = sa[16];         // 32
            bwg.definedFieldCount = sa[17];  // 34
            bwg.asOffset = (la[9] << 32) | (la[10]);    // 36 (unaligned longlong)
            bwg.totalSummaryOffset = (la[11] << 32) | (la[12]);    // 44 (unaligned longlong)
            bwg.uncompressBufSize = la[13];  // 52

            // dlog('bigType: ' + bwg.type);
            // dlog('chromTree at: ' + bwg.chromTreeOffset);
            // dlog('uncompress: ' + bwg.uncompressBufSize);
            // dlog('data at: ' + bwg.unzoomedDataOffset);
            // dlog('index at: ' + bwg.unzoomedIndexOffset);
            // dlog('field count: ' + bwg.fieldCount);
            // dlog('defined count: ' + bwg.definedFieldCount);

            bwg.zoomLevels = [];
            for (var zl = 0; zl < bwg.numZoomLevels; ++zl) {
                var zlReduction = la[zl*6 + 16];
                var zlData = (la[zl*6 + 18]<<32)|(la[zl*6 + 19]);
                var zlIndex = (la[zl*6 + 20]<<32)|(la[zl*6 + 21]);
                //          dlog('zoom(' + zl + '): reduction=' + zlReduction + '; data=' + zlData + '; index=' + zlIndex);
                bwg.zoomLevels.push({reductionLevel: zlReduction, dataOffset: zlData, indexOffset: zlIndex});
            }

            // parse the totalSummary if present (summary of all data in the file)
            if( bwg.totalSummaryOffset ) {
                if( Float64Array ) {
                    (function() {
                        var ua = new Uint32Array( header, bwg.totalSummaryOffset, 2 );
                        var da = new Float64Array( header, bwg.totalSummaryOffset+8, 4 );
                        var s = {
                            basesCovered: ua[0]<<32 | ua[1],
                            scoreMin: da[0],
                            scoreMax: da[1],
                            scoreSum: da[2],
                            scoreSumSquares: da[3]
                        };
                        bwg._globalStats = s;
                        // rest of these will be calculated on demand in getGlobalStats
                    }).call();
                } else {
                    console.warn("BigWig "+bwg.data.url+ " total summary not available, this web browser is not capable of handling this data type.");
                }
            } else {
                    console.warn("BigWig "+bwg.data.url+ " has no total summary data.");
            }

            bwg._readChromTree(
                function() {
                    bwg._deferred.features.resolve({success: true});
                    bwg._deferred.stats.resolve({success: true});
                },
                dojo.hitch( bwg, '_failAllDeferred' )
            );
        },
        dojo.hitch( this, '_failAllDeferred' )
       );
    },


    _readInt: function(ba, offset) {
        return (ba[offset + 3] << 24) | (ba[offset + 2] << 16) | (ba[offset + 1] << 8) | (ba[offset]);
    },

    _readShort: function(ba, offset) {
        return (ba[offset + 1] << 8) | (ba[offset]);
    },

    /**
     * @private
     */
    _readChromTree: function( callback, errorCallback ) {
        var thisB = this;
        this.refsByNumber = {};
        this.refsByName = {};

        var udo = this.unzoomedDataOffset;
        while ((udo % 4) != 0) {
            ++udo;
        }

        var readInt   = this._readInt;
        var readShort = this._readShort;

        this.data.slice( this.chromTreeOffset, udo - this.chromTreeOffset )
            .fetch(function(bpt) {
                       if( ! has('typed-arrays') ) {
                           thisB._failAllDeferred( 'Web browser does not support typed arrays' );
                           return;
                       }
                       var ba = new Uint8Array(bpt);
                       var la = new Int32Array(bpt, 0, 6);
                       var bptMagic = la[0];
                       if( bptMagic !== 2026540177 )
                           throw "parse error: not a Kent bPlusTree";
                       var blockSize = la[1];
                       var keySize = la[2];
                       var valSize = la[3];
                       var itemCount = (la[4] << 32) | (la[5]);
                       var rootNodeOffset = 32;

                       //dlog('blockSize=' + blockSize + '    keySize=' + keySize + '   valSize=' + valSize + '    itemCount=' + itemCount);

                       var bptReadNode = function(offset) {
                           if( offset >= ba.length )
                               throw "reading beyond end of buffer";
                           var isLeafNode = ba[offset];
                           var cnt = readShort( ba, offset+2 );
                           //dlog('ReadNode: ' + offset + '     type=' + isLeafNode + '   count=' + cnt);
                           offset += 4;
                           for (var n = 0; n < cnt; ++n) {
                               if( isLeafNode ) {
                                   // parse leaf node
                                   var key = '';
                                   for (var ki = 0; ki < keySize; ++ki) {
                                       var charCode = ba[offset++];
                                       if (charCode != 0) {
                                           key += String.fromCharCode(charCode);
                                       }
                                   }
                                   var refId = readInt( ba, offset );
                                   var refSize = readInt( ba, offset+4 );
                                   offset += 8;

                                   var refRec = { name: key, id: refId, length: refSize };

                                   //dlog(key + ':' + refId + ',' + refSize);
                                   thisB.refsByName[ thisB.browser.regularizeReferenceName(key) ] = refRec;
                                   thisB.refsByNumber[refId] = refRec;
                               } else {
                                   // parse index node
                                   offset += keySize;
                                   var childOffset = (readInt( ba, offset+4 ) << 32) | readInt( ba, offset );
                                   offset += 8;
                                   childOffset -= thisB.chromTreeOffset;
                                   bptReadNode(childOffset);
                               }
                           }
                       };
                       bptReadNode(rootNodeOffset);

                       callback.call( thisB, thisB );
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
        var thisB = this;
        seqName = thisB.browser.regularizeReferenceName( seqName );
        this._deferred.features.then(function() {
            callback( seqName in thisB.refsByName );
        }, errorCallback );
    },

    _getFeatures: function( query, featureCallback, endCallback, errorCallback ) {

        var chrName = this.browser.regularizeReferenceName( query.ref );
        var min = query.start;
        var max = query.end;

        var v = query.basesPerSpan ? this.getView( 1/query.basesPerSpan ) :
                       query.scale ? this.getView( scale )                :
                                     this.getView( 1 );

        if( !v ) {
            endCallback();
            return;
        }

        v.readWigData( chrName, min, max, dojo.hitch( this, function( features ) {
            array.forEach( features || [], featureCallback );
            endCallback();
        }), errorCallback );
    },

    getUnzoomedView: function() {
        if (!this.unzoomedView) {
            var cirLen = 4000;
            var nzl = this.zoomLevels[0];
            if (nzl) {
                cirLen = this.zoomLevels[0].dataOffset - this.unzoomedIndexOffset;
            }
            this.unzoomedView = new Window( this, this.unzoomedIndexOffset, cirLen, false );
        }
        return this.unzoomedView;
    },

    getView: function( scale ) {
        if( ! this.zoomLevels || ! this.zoomLevels.length )
            return null;

        if( !this._viewCache || this._viewCache.scale != scale ) {
            this._viewCache = {
                scale: scale,
                view: this._getView( scale )
            };
        }
        return this._viewCache.view;
    },

    _getView: function( scale ) {
        var basesPerSpan = 1/scale;
        //0 && console.log('getting view for '+basesPerSpan+' bases per span');
        for( var i = this.zoomLevels.length - 1; i > 0; i-- ) {
            var zh = this.zoomLevels[i];
            if( zh && zh.reductionLevel <= basesPerSpan ) {
                var indexLength = i < this.zoomLevels.length - 1
                    ? this.zoomLevels[i + 1].dataOffset - zh.indexOffset
                    : this.fileSize - 4 - zh.indexOffset;
                //0 && console.log( 'using zoom level '+i);
                return new Window( this, zh.indexOffset, indexLength, true );
            }
        }
        //0 && console.log( 'using unzoomed level');
        return this.getUnzoomedView();
    }
});

});