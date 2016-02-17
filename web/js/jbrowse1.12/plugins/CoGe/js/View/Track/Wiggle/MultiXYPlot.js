var coge_xyplot;
define( [
            'dojo/_base/declare',
            'dojo/_base/array',
            'dojo/_base/Color',
            'dojo/dom-construct',
            'dijit/Dialog',
            'JBrowse/View/InfoDialog',
            'JBrowse/View/Track/WiggleBase',
            'JBrowse/View/Track/_YScaleMixin',
            'JBrowse/Util',
            './_Scale',
            'CoGe/View/ColorDialog'
        ],
        function( declare, array, Color, domConstruct, Dialog, InfoDialog, WiggleBase, YScaleMixin, Util, Scale, ColorDialog ) {

var XYPlot = declare( [WiggleBase, YScaleMixin], // mdb: this file is a copy of XYPlot, extend that class instead?

/**
 * Wiggle track that shows data with an X-Y plot along the reference.
 *
 * @lends JBrowse.View.Track.Wiggle.XYPlot
 * @extends JBrowse.View.Track.WiggleBase
 */
{
    // Load cookie params - mdb added 1/13/14, issue 279
    constructor: function() {
        this.inherited(arguments); // call superclass constructor
        coge_xyplot = this;

        if (!this.config.style.featureColor)
            this.config.style.featureColor = {};
        if (typeof(this.config.showHoverScores) == "undefined")
        	this.config.showHoverScores = 1;
        if (typeof(this.config.showLabels) == "undefined")
        	this.config.showLabels = 1;
        if (typeof(this.config.showBackground) == "undefined")
        	this.config.showBackground = 0;
        if (typeof(this.config.disableZoomLimit) == "undefined")
        	this.config.disableZoomLimit = 0;
    },

    _defaultConfig: function() {
        return Util.deepUpdate(
            dojo.clone( this.inherited(arguments) ),
            {
                style: {
                    //pos_color: 'blue',
                    //neg_color: 'red',
                    origin_color: '#888',
                    variance_band_color: 'rgba(0,0,0,0.3)',
                },
            }
        );
    },

    _getScaling: function( viewArgs, successCallback, errorCallback ) {

        this._getScalingStats( viewArgs, dojo.hitch(this, function( stats ) {

            //calculate the scaling if necessary
            if( ! this.lastScaling || ! this.lastScaling.sameStats( stats ) ) {

                var scaling = new Scale( this.config, stats );

                // bump minDisplayed to 0 if it is within 0.5% of it
                if( Math.abs( scaling.min / scaling.max ) < 0.005 )
                    scaling.min = 0;

                // update our track y-scale to reflect it
                this.makeYScale({
                    fixBounds: true,
                    min: scaling.min,
                    max: scaling.max
                });

                // and finally adjust the scaling to match the ruler's scale rounding
                scaling.min = this.ruler.scaler.bounds.lower;
                scaling.max = this.ruler.scaler.bounds.upper;
                scaling.range = scaling.max - scaling.min;

                this.lastScaling = scaling;
            }

            successCallback( this.lastScaling );
        }), errorCallback );
    },

    updateStaticElements: function( coords ) {
        //console.log('updateStaticElements');
        this.inherited( arguments );
        this.updateYScaleFromViewDimensions( coords );
        _adjust_nav(this.config.coge.id)
    },

    fillTooManyFeaturesMessage: function( blockIndex, block, scale ) {
        this.fillMessage(
            blockIndex,
            block,
            'Too much data to show'
                + (scale >= this.browser.view.maxPxPerBp ? '': '; zoom in to see detail')
                + '.'
        );
    },

    fillMessage: function( blockIndex, block, message, class_ ) {
        domConstruct.empty( block.domNode );
        var msgDiv = dojo.create(
            'div', {
                className: class_ || 'message',
                innerHTML: message
            }, block.domNode );
        this.heightUpdate( dojo.position(msgDiv).h, blockIndex );
    },

    renderBlock: function( args ) {
        var featureScale = this.config.style.featureScale;
        var scale = args.block.scale;
        if (scale <= featureScale) { // don't draw, too zoomed-out, modeled after HTMLFeatures
            this.fillTooManyFeaturesMessage(args.blockIndex, args.block, scale);
        }
        else { // render features
            this.inherited( arguments );
        }
    },

    _draw: function(scale, leftBase, rightBase, block, canvas, features, featureRects, dataScale, pixels, spans) {
        this._preDraw(      scale, leftBase, rightBase, block, canvas, features, featureRects, dataScale );

        this._drawFeatures( scale, leftBase, rightBase, block, canvas, features, featureRects, dataScale );

        if ( spans ) {
            this._maskBySpans( scale, leftBase, rightBase, block, canvas, pixels, dataScale, spans );
        }
        this._postDraw(     scale, leftBase, rightBase, block, canvas, features, featureRects, dataScale );
    },
    /**
     * Draw a set of features on the canvas.
     * @private
     */
    _drawFeatures: function( scale, leftBase, rightBase, block, canvas, features, featureRects, dataScale ) {
        //console.log('_drawFeatures');
        var config = this.config;
        var context = canvas.getContext('2d');
        var canvasHeight = canvas.height;
        var toY = dojo.hitch( this, function( val ) {
           return canvasHeight * ( 1-dataScale.normalize.call(this, val) );
        });
        var originY = toY( dataScale.origin );
        var disableClipMarkers = this.config.disable_clip_markers;

        // mdb added 5/7/14 issue 374 - add background bar so that zero values are visible
        if (config.showBackground) {
            dojo.forEach( features, function(f,i) {
                var fRect = featureRects[i];
                var score = ( f.get('score') >= 0 ? 1 : -1 );
                fRect.t = toY(score);
                context.fillStyle = 'lightgray';
                if( fRect.t <= canvasHeight ) { // if the rectangle is visible at all
                    if (fRect.t <= originY) // bar goes upward
                        context.fillRect( fRect.l, fRect.t, fRect.w, originY-fRect.t+1);
                    else // downward
                        context.fillRect( fRect.l, originY, fRect.w, fRect.t-originY+1 );
                }
            }, this );
        }

        // Note: transform cases below can be consolidated/optimized
        if (config.transformAverage) {
            var sum_f = [];
            var sum_r = [];
            var count_f = [];
            var count_r = [];
            var width = [];
            dojo.forEach( features, function(f,i) {
                var score = f.get('score');
                var l = featureRects[i].l;
                var w = featureRects[i].w;
                if (score >= 0) {
                    sum_f[l] = l in sum_f ? sum_f[l] + score : score;
                    count_f[l] = l in count_f ? count_f[l] + 1 : 1;
                }
                else {
                    sum_r[l] = l in sum_r ? sum_r[l] + score : score;
                    count_r[l] = l in count_r ? count_r[l] + 1 : 1;
                }
                width[l] = l in width ? Math.max(width[l], w) : w;
            });
            sum_f.forEach( function(x,l) { // bar goes upward
                var avg = sum_f[l]/count_f[l];
                var height = toY( avg );
                context.fillStyle = 'gray';
                if( height <= canvasHeight ) { // if the rectangle is visible at all
                    context.fillRect( l, height, width[l], originY-height+1);
                }
            });
            sum_r.forEach( function(x,l) { // bar goes downward
                var avg = sum_r[l]/count_r[l];
                var height = toY( avg );
                context.fillStyle = 'gray';
                if( height <= canvasHeight ) { // if the rectangle is visible at all
                    context.fillRect( l, originY, width[l], height-originY+1 );
                }
            });
        }
        else if (config.transformDifference) {
            var width   = [];
            var max_f   = []; // forward strand
            var max_r   = []; // reverse strand
            var min_f   = [];
            var min_r   = [];
            var color_f = [];
            var color_r = [];
            var count_f = [];
            var count_r = [];
            dojo.forEach( features, dojo.hitch(this, function(f,i) {
                var score = f.get('score');
                var l = featureRects[i].l;
                var w = featureRects[i].w;
                width[l] = l in width ? Math.max(width[l], w) : w;

                var color = this._getFeatureColor( f.get('id') );

                if (score >= 0) {
                    count_f[l] = l in count_f ? count_f[l] + 1 : 1;
                    max_f[l] = l in max_f ? Math.max(max_f[l], score) : score;
                    min_f[l] = l in min_f ? Math.min(min_f[l], score) : score;
                    if (score >= max_f[l])
                        color_f[l] = color;
                }
                else {
                    count_r[l] = l in count_r ? count_r[l] + 1 : 1;
                    score = Math.abs(score);
                    max_r[l] = l in max_r ? Math.max(max_r[l], score) : score;
                    min_r[l] = l in min_r ? Math.min(min_r[l], score) : score;
                    if (score >= max_r[l])
                        color_r[l] = color;
                }
            }));
            max_f.forEach( function(x,l) { // bar goes upward
                var diff = max_f[l] - min_f[l];
                if (count_f[l] == 1)
                    diff = max_f[l];
                var height = toY(diff);
                if( height <= canvasHeight ) { // if the rectangle is visible at all
                    context.fillStyle = color_f[l];
                    context.fillRect( l, height, width[l], originY-height+1);
                }
            });
            max_r.forEach( function(x,l) { // bar goes downward
                var diff = max_r[l] - min_r[l];
                if (count_r[l] == 1)
                    diff = max_r[l];
                var height = toY(-1*diff);
                if( height <= canvasHeight ) { // if the rectangle is visible at all
                    context.fillStyle = color_r[l];
                    context.fillRect( l, originY, width[l], height-originY+1 );
                }
            });
        }
        else if (config.transformInflate) {
            // sort features by score
            var sorted = [];
            dojo.forEach( features, function(f,i) {
                sorted.push({ feature: f, featureRect: featureRects[i] });
            });
            sorted.sort( sortByScore );

            dojo.forEach( sorted, function(pair,i) {
                var f = pair.feature;
                var fRect = pair.featureRect;
                var score = ( f.get('score') > 0 ? 1 : -1 );

                fRect.t = toY( score );
                if( fRect.t <= canvasHeight ) { // if the rectangle is visible at all
                    var id = f.get('id');
                    context.fillStyle = this._getFeatureColor(id);
                    if (fRect.t <= originY) // bar goes upward
                        context.fillRect( fRect.l, fRect.t, fRect.w, originY-fRect.t+1);
                    else // downward
                        context.fillRect( fRect.l, originY, fRect.w, fRect.t-originY+1 );
                }
            }, this );
        }
        else {
            // sort features by score
            var sorted = [];
            dojo.forEach( features, function(f,i) {
                sorted.push({ feature: f, featureRect: featureRects[i] });
            });
            sorted.sort( sortByScore );

            dojo.forEach( sorted, function(pair,i) {
                var f = pair.feature;
                var fRect = pair.featureRect;
                var score = f.get('score');

                if (config.transformLog10) {
                    if (score >= 0)
                        score = log10(Math.abs(score)+1);
                    else
                        score = -1*log10(Math.abs(score)+1);
                }
                else if (config.transformLog2) {
                    if (score >= 0)
                        score = log2(Math.abs(score)+1);
                    else
                        score = -1*log2(Math.abs(score)+1);
                }

                fRect.t = toY( score );
                if( fRect.t <= canvasHeight ) { // if the rectangle is visible at all
                    var id = f.get('id');
                    context.fillStyle = this._getFeatureColor(id);
                    if (fRect.t <= originY) // bar goes upward
                        context.fillRect( fRect.l, fRect.t, fRect.w, originY-fRect.t+1);
                    else // downward
                        context.fillRect( fRect.l, originY, fRect.w, fRect.t-originY+1 );
                }
            }, this );
        }

        // mdb added 4/2/14 - draw labels on top of bars, issue 346
        var prevStart, prevEnd;
        if (config.showLabels && scale > config.style.labelScale) {
            dojo.forEach( sorted, function(pair,i) {
                var f = pair.feature;
                var fRect = pair.featureRect;
                var isUpward = (fRect.t <= originY); // bar goes upward
                var start = f.get('start');
                if (start >= block.startBase && start <= block.endBase) { // print label only for first spanning block
	                var label = f.get('label');
	                if (label && label != '.') {
	                	if (!(start >= prevStart && start <= prevEnd)) { // mdb added 4/15/14 - don't allow overlapping labels, only print the first one
		                	var topOffset = ( isUpward ? fRect.t-12 : fRect.t );
		                    var rulerdiv =
		                        dojo.create('div',
		                    		{   style: {
		                                	width: '100px',
		                                    position: 'absolute',
		                                    left: fRect.l,
		                                    top: topOffset,
		                                    //zIndex: 10,
		                                },
		                                innerHTML: label
		                            }, canvas.parentNode );
		                    prevStart = start;
		                    prevEnd = f.get('end');
	                	}
	                }
                }
            }, this );
        }
    },

    /**
     * Draw anything needed after the features are drawn.
     */
    _postDraw: function( scale, leftBase, rightBase, block, canvas, features, featureRects, dataScale ) {
//      console.log('_postDraw');
        var context = canvas.getContext('2d');
        var canvasHeight = canvas.height;
        var toY = dojo.hitch( this, function( val ) {
           return canvasHeight * (1-dataScale.normalize.call(this, val));
        });

        // draw the variance_band if requested
        if( this.config.variance_band ) {
            var bandPositions =
                typeof this.config.variance_band == 'object'
                    ? array.map( this.config.variance_band, function(v) { return parseFloat(v); } ).sort().reverse()
                    : [ 2, 1 ];
            this.getGlobalStats( dojo.hitch( this, function( stats ) {
                if( ('scoreMean' in stats) && ('scoreStdDev' in stats) ) {
                    var drawVarianceBand = function( plusminus, fill, label ) {
                        context.fillStyle = fill;
                        var varTop = toY( stats.scoreMean + plusminus );
                        var varHeight = toY( stats.scoreMean - plusminus ) - varTop;
                        varHeight = Math.max( 1, varHeight );
                        context.fillRect( 0, varTop, canvas.width, varHeight );
                        context.font = '12px sans-serif';
                        if( plusminus > 0 ) {
                            context.fillText( '+'+label, 2, varTop );
                            context.fillText( '-'+label, 2, varTop+varHeight );
                        }
                        else {
                            context.fillText( label, 2, varTop );
                        }
                    };

                    var maxColor = new Color( this.config.style.variance_band_color );
                    var minColor = new Color( this.config.style.variance_band_color );
                    minColor.a /= bandPositions.length;

                    var bandOpacityStep = 1/bandPositions.length;
                    var minOpacity = bandOpacityStep;

                    array.forEach( bandPositions, function( pos,i ) {
                        drawVarianceBand( pos*stats.scoreStdDev,
                                          Color.blendColors( minColor, maxColor, (i+1)/bandPositions.length).toCss(true),
                                          pos+'σ');
                    });
                    drawVarianceBand( 0, 'rgba(255,255,0,0.7)', 'mean' );
                }
            }));
        }

        // draw the origin line if it is not disabled
        var originColor = this.config.style.origin_color;
        if( typeof originColor == 'string' && !{'none':1,'off':1,'no':1,'zero':1}[originColor] ) {
            var originY = toY( dataScale.origin );
            context.fillStyle = originColor;
            context.fillRect( 0, originY, canvas.width-1, 1 );
        }
    },

    _calculatePixelScores: function( canvasWidth, features, featureRects ) {
        var pixelValues = new Array( canvasWidth );

        if (this.config.showHoverScores) {
            // sort features by score
            var sorted = [];
            dojo.forEach( features, function(f,i) {
                sorted.push({ feature: f, featureRect: featureRects[i] });
            });
            sorted.sort( sortByScore );

            // make an array of the max score at each pixel on the canvas
            dojo.forEach( sorted, function( f, i ) {
                var fRect = f.featureRect;
                var jEnd = fRect.r;
                var score = f.feature.get('score');
                var score2 = f.feature.get('score2');
                var id = f.feature.get('id');
                var fLabel = f.feature.get('label');
                if (!fLabel || fLabel == '.') fLabel = '';
                var name = this._getFeatureName(f.feature);
                var color = this._getFeatureColor(id);
                for( var j = Math.round(fRect.l); j < jEnd; j++ ) {
                    var label = '<div style="background-color:'+color+';">' +
                        nbspPad(score.toPrecision(6).toString(), 11) +
                        (score2 ? nbspPad(score2.toPrecision(6).toString(), 11) : '') +
                        fLabel+ ' ' + name + '</div>';
                    pixelValues[j] = j in pixelValues ? pixelValues[j] + label : label;
                }
            },this);

            // compute transform scores - FIXME dup'ed in _drawFeatures
            if (this.config.transformAverage) {
                var sum_f = [];
                var sum_r = [];
                var count_f = [];
                var count_r = [];
                var width = [];
                dojo.forEach( features, function(f,i) {
                    var score = f.get('score');
                    var l = featureRects[i].l;
                    var w = featureRects[i].w;
                    if (score >= 0) {
                        sum_f[l] = l in sum_f ? sum_f[l] + score : score;
                        count_f[l] = l in count_f ? count_f[l] + 1 : 1;
                    }
                    else {
                        sum_r[l] = l in sum_r ? sum_r[l] + score : score;
                        count_r[l] = l in count_r ? count_r[l] + 1 : 1;
                    }
                    width[l] = l in width ? Math.max(width[l], w) : w;
                });

                sum_f.forEach( function(x,l) {
                    var avg = sum_f[l]/count_f[l];
                    for( var j = Math.round(l); j < l+width[l]; j++ ) {
                        var label = '<div style="background-color:gray;">' +
                            nbspPad(avg.toPrecision(6).toString(), 11)
                            + 'Average (+)' + '</div>';
                        pixelValues[j] = j in pixelValues ? pixelValues[j] + label : label;
                    }
                });
                sum_r.forEach( function(x,l) {
                    var avg = sum_r[l]/count_r[l];
                    for( var j = Math.round(l); j < l+width[l]; j++ ) {
                        var label = '<div style="background-color:gray;">' +
                            nbspPad(avg.toPrecision(6).toString(), 11)
                            + 'Average (-)' + '</div>';
                        pixelValues[j] = j in pixelValues ? pixelValues[j] + label : label;
                    }
                });
            }
            else if (this.config.transformDifference) {
                var max_f = [];
                var max_r = [];
                var min_f = [];
                var min_r = [];
                var count_f = [];
                var count_r = [];
                var width = [];
                dojo.forEach( features, function(f,i) {
                    var score = f.get('score');
                    var l = featureRects[i].l;
                    var w = featureRects[i].w;
                    width[l] = l in width ? Math.max(width[l], w) : w;
                    if (score >= 0) {
                        count_f[l] = l in count_f ? count_f[l] + 1 : 1;
                        max_f[l] = l in max_f ? Math.max(max_f[l], score) : score;
                        min_f[l] = l in min_f ? Math.min(min_f[l], score) : score;
                    }
                    else {
                        score = Math.abs(score);
                        count_r[l] = l in count_r ? count_r[l] + 1 : 1;
                        max_r[l] = l in max_r ? Math.max(max_r[l], score) : score;
                        min_r[l] = l in min_r ? Math.min(min_r[l], score) : score;
                    }
                });

                max_f.forEach( function(x,l) {
                    var diff = max_f[l] - min_f[l];
                    if (count_f[l] == 1)
                        diff = max_f[l];
                    for( var j = Math.round(l); j < l+width[l]; j++ ) {
                        var label = '<div style="background-color:gray;">' +
                            nbspPad(diff.toPrecision(6).toString(), 11)
                            + 'Difference (+)' + '</div>';
                        pixelValues[j] = j in pixelValues ? pixelValues[j] + label : label;
                    }
                });
                max_r.forEach( function(x,l) {
                    var diff = max_r[l] - min_r[l];
                    if (count_r[l] == 1)
                        diff = max_r[l];
                    for( var j = Math.round(l); j < l+width[l]; j++ ) {
                        var label = '<div style="background-color:gray;">' +
                            nbspPad(diff.toPrecision(6).toString(), 11)
                            + 'Difference (-)' + '</div>';
                        pixelValues[j] = j in pixelValues ? pixelValues[j] + label : label;
                    }
                });
            }
        }

        return pixelValues;
    },

    _error: function(title, content) {
    	if (content.responseText)
    		content = content.responseText;
    	else if (content.error)
    		content = JSON.stringify(content.error);
    	new InfoDialog({
    		title: title,
    		content: content,
            onHide: function(){this.destroyRecursive()}
    	}).show();	
    },

    _search_track: function() {
		var eid = dojo.byId('eid').value;
		var div = dojo.byId('coge-track-search');
		var params;

    	if (dojo.byId('highest').checked)
    		params = 'type=max';
    	else if (dojo.byId('lowest').checked)
    		params = 'type=min';
    	else {
    		if (this._brush.empty()) {
    			this._error('Unspecified Range', 'Please drag on the histogram to select the range of values you wish to search for');
    			return;
    		}
    		var params = 'type=range';
    		var extent = this._brush.extent();
    		var domain = this._brush.x().domain();
    		if (extent[0] != domain[0])
    			params += '&gte=' + extent[0];
    		if (extent[1] != domain[1])
    			params += '&lte=' + extent[1];
    	}
    	var ref_seq = dojo.byId('coge_ref_seq');
    	if (ref_seq.selectedIndex > 0)
    		params += '&chr=' + ref_seq.options[ref_seq.selectedIndex].innerHTML;

		dojo.empty(div);
		div.innerHTML = '<img src="picts/ajax-loader.gif">';
     	dojo.xhrGet({
    		url: api_base_url + '/experiment/' + eid + '/query?' + params,
    		handleAs: 'json',
	  		load: dojo.hitch(this, function(data) {
	  			if (data.length == 0) {
	  				this._error('Search', 'Search returned zero hits');
	  				this._track_search_dialog.hide();
	  				return;
	  			}
	  			this._new_nav(eid, data);
 				this._track_search_dialog.hide();
    		}),
    		error: dojo.hitch(this, function(data) {
    			this._error('Search', data);
 				this._track_search_dialog.hide();
    		})
    	});
    },
    
    _new_nav: function(eid, data) {
		var first = JSON.parse('[' + data[0] + ']');
		this._track.browser.navigateToLocation({
			ref: first[0],
			start: first[1],
			end: first[2]
		});
		var nav = dojo.byId('nav_' + eid);
		if (nav)
			dojo.destroy(nav);
		nav = dojo.create('div', { id: 'nav_' + eid, style: { background: 'white', opacity: 0.7, position: 'absolute' } }, dojo.byId('container'));
		_adjust_nav(eid);
		nav.hits = data;
		nav.hit = 0;
		nav.browser = this._track.browser;
		dojo.create('span', { className: 'glyphicon glyphicon-step-backward', onclick: function() { _go_to_hit(nav, 0) }, style: { cursor: 'pointer' } }, nav);
		dojo.create('span', { className: 'glyphicon glyphicon-chevron-left', onclick: function() { if (nav.hit > 0) _go_to_hit(nav, nav.hit - 1) }, style: { cursor: 'pointer' } }, nav);
		nav.num_span = dojo.create('span', { innerHTML: '1', style: { cursor: 'default' } }, nav);
		dojo.create('span', { innerHTML: ' of ' + data.length + ' hit' + (data.length != 1 ? 's ' : ' '), style: { cursor: 'default', marginRight: 5 } }, nav);
		dojo.create('span', { className: 'glyphicon glyphicon-chevron-right', onclick: function() { if (nav.hit < nav.hits.length - 1) _go_to_hit(nav, nav.hit + 1) }, style: { cursor: 'pointer' } }, nav);
		dojo.create('span', { className: 'glyphicon glyphicon-step-forward', onclick: function() { _go_to_hit(nav, nav.hits.length - 1) }, style: { cursor: 'pointer' } }, nav);
        this.browser.subscribe('/jbrowse/v1/v/tracks/hide', function(configs) {
        	for (var i=0; i<configs.length; i++)
        		if (configs[i].coge.id == eid)
        			dojo.destroy(dojo.byId('nav_' + eid));
        });
    },

    _showPixelValue: function( scoreDisplay, score ) {
        var scoreType = typeof score;
        if( scoreType == 'number' ) {
            // display the score with only 6
            // significant digits, avoiding
            // most confusion about the
            // approximative properties of
            // IEEE floating point numbers
            // parsed out of BigWig files
            scoreDisplay.innerHTML = parseFloat( score.toPrecision(6) );
            return true;
        }
        else if( scoreType == 'string' ) {
            scoreDisplay.innerHTML = score;
            return true;
        }
        else {
            return false;
        }
    },

    _trackMenuOptions: function() {
        var options = this.inherited(arguments);
        var track = this;
        var config = this.config;
        if (config.coge.menuOptions) {
            config.coge.menuOptions.forEach( function(e) {
                options.push(e);
            });
        }

        options.push.apply(
            options,
            [
                { type: 'dijit/MenuSeparator' },
                {
                    label: 'Show scores on hover',
                    type: 'dijit/CheckedMenuItem',
                    checked: this.config.showHoverScores,
                    onClick: function(event) {
                        track.config.showHoverScores = this.checked;
                        track.changed();
                    }
                },
                {
                    label: 'Show labels',
                    type: 'dijit/CheckedMenuItem',
                    checked: this.config.showLabels,
                    onClick: function(event) {
                        track.config.showLabels = this.checked;
                        track.changed();
                    }
                },
                {
                    label: 'Show background',
                    type: 'dijit/CheckedMenuItem',
                    checked: this.config.showBackground,
                    onClick: function(event) {
                        track.config.showBackground = this.checked;
                        track.changed();
                    }
                },
                {
                    label: 'Change colors',
                    onClick: function(event) {
                        if (!track.colorDialog) {
                            track.colorDialog = new ColorDialog({
                                title: "Change colors",
                                style: {
                                    width: '230px',
                            },
                            items: track.config.coge.experiments || [track.config.coge],
                            featureColor: track.config.style.featureColor,
                            callback: function(id, color) {
                                var curColor = track.config.style.featureColor[id];
                                if (!curColor || curColor != color) {
                                    // Save color choice
                                    track.config.style.featureColor[id] = color;

                                    track.updateUserStyles({ featureColor : track.config.style.featureColor });

                                    // Repaint track
                                    track.changed();

                                    //FIXME TrackList should update itself
                                    track.browser.publish('/jbrowse/v1/c/tracks/show', [track.config]);
                                    }
                                }
                            });
                        }
                        track.colorDialog.show();
                    }
                },
                { // mdb added 11/6/15 COGE-678
                    label: 'Disable zoom limit',
                    type: 'dijit/CheckedMenuItem',
                    checked: this.config.disableZoomLimit,
                    onClick: function(event) {
                    	var config = track.config;
                        config.disableZoomLimit = this.checked;
                        if (config.disableZoomLimit) {
                        	config.savedFeatureScale = config.style.featureScale;
                        	config.style.featureScale = 0;
                        }
                        else {
                        	config.style.featureScale = config.savedFeatureScale;
                        }
                        track.changed();
                    }
                }
            ]);

        options.push.apply(options, [{
	        label: 'Search',
	        onClick: function(event) {
	        	coge_xyplot._track = track;
	        	var content = '<div id="coge-track-search"><input type="hidden" id="eid" value="' + track.config.coge.id + '"><table align="center"><tr><td>RefSeq:</td><td><select id="coge_ref_seq"><option>Any</option>';
	        	coge_xyplot.browser.refSeqOrder.forEach(function(rs){
	        		content += '<option';
	        		if (rs == coge_xyplot.browser.refSeq.name)
	        			content += ' selected';
	        		content += '>' + rs + '</option>';
	        	})
	        	content += '</select></td></tr>' +
	        		'<tr><td>Values:</td><td style="white-space:nowrap"><input id="highest" type="radio" name="type" checked="checked"> highest</td></tr>' +
	        		'<tr><td></td><td style="white-space:nowrap"><input id="lowest" type="radio" name="type"> lowest</td></tr>' +
	        		'<tr><td></td><td style="white-space:nowrap" valign="top"><input id="range" type="radio" name="type"> range: <span id="selected_range">&nbsp;</span></td></tr>' +
	        		'<tr><td></td><td><div id="coge-hist"><img src="picts/ajax-loader.gif"></div></td></tr></table><div class="dijitDialogPaneActionBar"><button data-dojo-type="dijit/form/Button" type="button" onClick="coge_xyplot._search_track()">OK</button><button data-dojo-type="dijit/form/Button" type="button" onClick="coge_xyplot._track_search_dialog.hide()">Cancel</button></div></div>';
	        	coge_xyplot._track_search_dialog = new Dialog({
                    title: 'Search Track',
                    content: content,
                    onHide: function(){this.destroyRecursive()},
                    style: "width: 300px"
                });
	        	dojo.xhrGet({
	        		url: api_base_url + '/experiment/' + track.config.coge.id + '/histogram',
	        		handleAs: 'json',
	    	  		load: function(data) {
	    	  			if (data.error) {
	    	  				coge_xyplot._error('Search', data);
	    	  				coge_xyplot._track_search_dialog.hide();
	    	  				return;
	    	  			}
	    	  			dojo.destroy(dojo.byId('coge-hist').firstChild);
	    	  			coge_xyplot._brush = chart(d3.select('#coge-hist'), data.first, data.gap, data.counts);
	        		},
	        		error: function(data) {
	        			coge_xyplot._error('Search', data);
	        			coge_xyplot._track_search_dialog.hide();
	        		}
	        	});
	        	coge_xyplot._track_search_dialog.show();
	        }
	    }]);

        if (config.coge.type == 'notebook') {
            options.push.apply(
                    options,
                    [
                        // Note: would prefer a radio submenu but this is
                        // dojo 1.8 and RadioMenuItem doesn't exist until
                        // dojo 1.9.
                        {   label: 'Transform',
                            type: 'dijit/DropDownMenu',
                            children: [
                                {   label: 'None',
                                    onClick: function(event) {
                                        clearTransforms(config);
                                        track.changed();
                                    }
                                },
                                {   label: 'Average',
                                    onClick: function(event) {
                                        clearTransforms(config);
                                        track.config.transformAverage = true;
                                        track.changed();
                                    }
                                },
                                {   label: 'Difference',
                                    onClick: function(event) {
                                        clearTransforms(config);
                                        track.config.transformDifference = true;
                                        track.changed();
                                    }
                                },
                                {   label: 'Log10',
                                    onClick: function(event) {
                                        clearTransforms(config);
                                        track.config.transformLog10 = true;
                                        track.changed();
                                    }
                                },
                                {   label: 'Log2',
                                    onClick: function(event) {
                                        clearTransforms(config);
                                        track.config.transformLog2 = true;
                                        track.changed();
                                    }
                                }
                            ]
                        }
                    ]
                );
        }
        else { // individual experiment track
            options.push.apply(
                    options,
                    [
                        // Note: would prefer a radio submenu but this is
                        // dojo 1.8 and RadioMenuItem doesn't exist until
                        // dojo 1.9.
                        {   label: 'Transform',
                            type: 'dijit/DropDownMenu',
                            children: [
                                {   label: 'None',
                                    onClick: function(event) {
                                        clearTransforms(config);
                                        track.changed();
                                    }
                                },
                                {   label: 'Log10',
                                    onClick: function(event) {
                                        clearTransforms(config);
                                        track.config.transformLog10 = true;
                                        track.changed();
                                    }
                                },
                                {   label: 'Log2',
                                    onClick: function(event) {
                                        clearTransforms(config);
                                        track.config.transformLog2 = true;
                                        track.changed();
                                    }
                                },
                                {   label: 'Inflate',
                                    onClick: function(event) {
                                        clearTransforms(config);
                                        track.config.transformInflate = true;
                                        track.changed();
                                    }
                                }
                            ]
                        }
                    ]
                );
        }

        return options;
    },

    _getFeatureName: function(f) {
        var id = f.get('id');
        var coge = this.config.coge;

        if (coge.type == 'experiment') {
            return coge.name;
        }
        else if (coge.type == 'notebook') {
            var experiments = coge.experiments || [];
            var name = '';
            experiments.every( function(e) {
                if (e.id == id) {
                    name = e.name;
                    if (e.type == 'snp')
                        name += ' ' + f.get('name');
                    return false;
                }
                return true;
            });
            return name;
        }
    },

    _getFeatureColor: function(id) {
        if (this.config.style.featureColor && this.config.style.featureColor[id]) {
            return this.config.style.featureColor[id];
        }
        return '#' + ((((id * 1234321) % 0x1000000) | 0x444444) & 0xe7e7e7 ).toString(16); //FIXME: dup'ed in CoGe.js
    }

});

return XYPlot;
});

function log10(x) {
    return Math.log(x) / Math.log(10);
}

function log2(x) {
    return Math.log(x) / Math.log(2);
}

function clearTransforms(config) {
	// FIXME change to config.transforms.average, etc...
    config.transformAverage = false;
    config.transformDifference = false;
    config.transformLog10 = false;
    config.transformLog2 = false;
    config.transformInflate = false;
}

function nbspPad(s, padLength) {
    for (var i = s.length;  i < padLength;  i++) {
        s += '&nbsp;';
    }
    return s;
}

function sortByScore(a,b) { // sort features by score
    return Math.abs( b.feature.get('score') ) - Math.abs( a.feature.get('score') );
}

function _adjust_nav(eid) {
 	var l = dojo.byId('label_experiment' + eid);
 	if (l) {
 		var nav = dojo.byId('nav_' + eid);
 		if (nav) {
	 		var track = dojo.byId('track_experiment' + eid);
	     	dojo.style(nav, 'left', dojo.style(l, 'left') + 10);
	     	dojo.style(nav, 'top', dojo.style(track, 'top') + 32);
 		}
 	}
}

function _go_to_hit(nav, hit) {
	nav.hit = hit % nav.hits.length;
	nav.num_span.innerHTML = nav.hit + 1;
	var loc = JSON.parse('[' + nav.hits[nav.hit] + ']');
	nav.browser.navigateToLocation({
		ref: loc[0],
		start: loc[1],
		end: loc[2]
	});	
}

function chart(div, first, gap, counts) {
	var margin = {top: 10, right: 10, bottom: 20, left: 10},
		x = d3.scale.linear().range([0, 200]).domain([first, first + gap * counts.length]),
		y = d3.scale.linear().range([100, 0]).domain([0, d3.max(counts)]),
		brush = d3.svg.brush(),
	    axis = d3.svg.axis().orient("bottom").scale(x),
	    width = x.range()[1],
      	height = y.range()[0];
  
	var g = div.append("svg")
      	.attr("width", width + margin.left + margin.right)
      	.attr("height", height + margin.top + margin.bottom)
      .append("g")
      	.attr("transform", "translate(" + margin.left + "," + margin.top + ")");

	g.append("clipPath")
      	.attr("id", "clip")
      .append("rect")
      	.attr("width", width)
      	.attr("height", height);

	g.selectAll(".bar")
      	.data(["background", "foreground"])
      .enter().append("path")
      	.attr("class", function(d) { return d + " bar"; })
      	.datum(function(d,i){return counts[i]});//group.all());

	g.selectAll(".foreground.bar")
      	.attr("clip-path", "url(#clip)");

	g.append("g")
      	.attr("class", "axis")
      	.attr("transform", "translate(0," + height + ")")
      	.call(axis);

	// Initialize the brush component with pretty resize handles.
	brush.x(x);
	var gBrush = g.append("g").attr("class", "brush").call(brush);
	gBrush.selectAll("rect").attr("height", height);
	gBrush.selectAll(".resize").append("path").attr("d", resizePath);

    g.selectAll(".bar").attr("d", barPath);

	function barPath() {
    	var path = [],
            i = -1,
            n = counts.length,
            X = first;
        while (++i < n) {
        	path.push("M", x(X), ",", height, "V", y(counts[i]), "h9V", height);
        	X += gap;
        }
        return path.join("");
	}

	function resizePath(d) {
        var e = +(d == "e"),
            x = e ? 1 : -1,
            y = height / 3;
        return "M" + (.5 * x) + "," + y
            + "A6,6 0 0 " + e + " " + (6.5 * x) + "," + (y + 6)
            + "V" + (2 * y - 6)
            + "A6,6 0 0 " + e + " " + (.5 * x) + "," + (2 * y)
            + "Z"
            + "M" + (2.5 * x) + "," + (y + 8)
            + "V" + (2 * y - 8)
            + "M" + (4.5 * x) + "," + (y + 8)
            + "V" + (2 * y - 8);
	}

    brush.on("brushstart.chart", function() {
    	var div = d3.select(this.parentNode.parentNode.parentNode);
    	div.select(".title a").style("display", null);
    	dojo.byId('range').checked = true;
    });

    brush.on("brush.chart", function() {
    	var g = d3.select(this.parentNode),
          	extent = brush.extent();
//      if (round) g.select(".brush")
//          .call(brush.extent(extent = extent.map(round)))
//        .selectAll(".resize")
//          .style("display", null);
    	g.select("#clip rect")
          	.attr("x", x(extent[0]))
          	.attr("width", x(extent[1]) - x(extent[0]));
    	dojo.byId('selected_range').innerHTML = parseFloat(Math.round(extent[0] * 100) / 100).toFixed(2) + ' .. ' + parseFloat(Math.round(extent[1] * 100) / 100).toFixed(2);
    });

    brush.on("brushend.chart", function() {
    	if (brush.empty()) {
    		var div = d3.select(this.parentNode.parentNode.parentNode);
    		div.select(".title a").style("display", "none");
    		div.select("#clip rect").attr("x", null).attr("width", "100%");
    		dojo.byId('selected_range').innerHTML = '&nbsp;';
    	}
    });
    d3.rebind(chart, brush, "on");
    return brush;
}
