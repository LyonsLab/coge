var coge_xyplot;
define( [
			'dojo/_base/declare',
			'dojo/_base/array',
			'dojo/_base/Color',
			'dojo/dom-construct',
			'dijit/Dialog',
			'JBrowse/View/Track/Wiggle/XYPlot',
			'JBrowse/Util',
			'./_Scale',
			'CoGe/View/ColorDialog',
            'JBrowse/Store/LRUCache'
		],
		function( declare, array, Color, domConstruct, Dialog, XYPlotBase, Util, Scale, ColorDialog, LRUCache ) {

var XYPlot = declare( [XYPlotBase], {
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

	// ----------------------------------------------------------------

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
				if (!fLabel || fLabel == '.')
					fLabel = '';
				var name = this._getFeatureName(f.feature);
				var color = this._getFeatureColor(id);
				for( var j = Math.round(fRect.l); j < jEnd; j++ ) {
					var label = '<div style="background-color:' + color + ';">' +
						nbspPad(score.toPrecision(6).toString(), 11) +
						(score2 ? nbspPad(score2.toPrecision(6).toString(), 11) : '') +
						fLabel+ ' ' + name + '&nbsp;&nbsp;' + f.feature.get('start') + '..' + f.feature.get('end') + '</div>';
					pixelValues[j] = j in pixelValues ? pixelValues[j] + label : label;
				}
			}, this);

			// compute transform scores - FIXME dup'ed in _drawFeatures
			if (this.config.coge.transform == 'Average') {
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
						var label = '<div style="background-color:gray;">' + nbspPad(avg.toPrecision(6).toString(), 11) + 'Average (+)' + '</div>';
						pixelValues[j] = j in pixelValues ? pixelValues[j] + label : label;
					}
				});
				sum_r.forEach( function(x,l) {
					var avg = sum_r[l]/count_r[l];
					for( var j = Math.round(l); j < l+width[l]; j++ ) {
						var label = '<div style="background-color:gray;">' + nbspPad(avg.toPrecision(6).toString(), 11) + 'Average (-)' + '</div>';
						pixelValues[j] = j in pixelValues ? pixelValues[j] + label : label;
					}
				});
			}
			else if (this.config.coge.transform == 'Difference') {
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
						var label = '<div style="background-color:gray;">' + nbspPad(diff.toPrecision(6).toString(), 11) + 'Difference (+)' + '</div>';
						pixelValues[j] = j in pixelValues ? pixelValues[j] + label : label;
					}
				});
				max_r.forEach( function(x,l) {
					var diff = max_r[l] - min_r[l];
					if (count_r[l] == 1)
						diff = max_r[l];
					for( var j = Math.round(l); j < l+width[l]; j++ ) {
						var label = '<div style="background-color:gray;">' + nbspPad(diff.toPrecision(6).toString(), 11) + 'Difference (-)' + '</div>';
						pixelValues[j] = j in pixelValues ? pixelValues[j] + label : label;
					}
				});
			}
		}

		return pixelValues;
	},

	// ----------------------------------------------------------------
    // utility method that calculates standard deviation from sum and sum of squares

    _calcStdFromSums: function( sum, sumSquares, n ) {
        if( n == 0 )
            return 0;

        var variance = sumSquares - sum*sum/n;
        if (n > 1) {
            variance /= n-1;
        }
        return variance < 0 ? 0 : Math.sqrt(variance);
    },

	// ----------------------------------------------------------------

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

	// ----------------------------------------------------------------

	_draw: function(scale, leftBase, rightBase, block, canvas, features, featureRects, dataScale, pixels, spans) {
		this._preDraw(scale, leftBase, rightBase, block, canvas, features, featureRects, dataScale);
		this._drawFeatures(scale, leftBase, rightBase, block, canvas, features, featureRects, dataScale);
		if (spans)
			this._maskBySpans(scale, leftBase, rightBase, block, canvas, pixels, dataScale, spans);
		this._postDraw(scale, leftBase, rightBase, block, canvas, features, featureRects, dataScale);
	},

	// ----------------------------------------------------------------
	/**
	 * Draw a set of features on the canvas.
	 * @private
	 */
	_drawFeatures: function( scale, leftBase, rightBase, block, canvas, features, featureRects, dataScale ) {
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

		if (config.coge.transform == 'Average') {
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
		else if (config.coge.transform == 'Difference') {
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

				if (config.coge.transform == 'Inflate')
					score = score > 0 ? dataScale.max : dataScale.min;
				else if (config.coge.transform == 'Log10') {
					if (score >= 0)
						score = log10(Math.abs(score)+1);
					else
						score = -1*log10(Math.abs(score)+1);
				}
				else if (config.coge.transform == 'Log2') {
					if (score >= 0)
						score = log2(Math.abs(score)+1);
					else
						score = -1*log2(Math.abs(score)+1);
				}
				else if (config.coge.transform == 'Normalize')
					score /= config.coge.max;

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
							dojo.create('div',
								{   style: {
										width: '100px',
										position: 'absolute',
										left: fRect.l,
										top: topOffset
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

	// ----------------------------------------------------------------

	fillMessage: function( blockIndex, block, message, class_ ) {
		domConstruct.empty( block.domNode );
		var msgDiv = dojo.create(
			'div', {
				className: class_ || 'message',
				innerHTML: message,
				style: { 'margin-top': '30px' }
			}, block.domNode );
		this.heightUpdate( 100, //dojo.position(msgDiv).h,
			blockIndex );
	},

	// ----------------------------------------------------------------

	fillTooManyFeaturesMessage: function( blockIndex, block, scale ) {
		this.fillMessage(
			blockIndex,
			block,
			'Too much data to show' + (scale >= this.browser.view.maxPxPerBp ? '.': '; zoom in to see detail.')
		);
	},

	// ----------------------------------------------------------------

	_getFeatureColor: function(id) {
		if (this.config.coge.type == 'search')
			id = this.config.coge.eid;
		if (this.config.style.featureColor && this.config.style.featureColor[id])
			return this.config.style.featureColor[id];
		 return coge_plugin.calc_color(id);
	},

	// ----------------------------------------------------------------

	_getFeatureName: function(f) {
		var id = f.get('id');
		var coge = this.config.coge;

		if (coge.type == 'experiment')
			return coge.name;
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
 
	// ----------------------------------------------------------------

   getGlobalStats: function( successCallback, errorCallback ) {
		if (this.config.coge.transform == 'Normalize') {
			if (!this.config.coge.max) {
				var coge = this.config.coge;
				this.store.getGlobalStats( function(stats) {
					coge.max = Math.max(stats.scoreMax, Math.abs(stats.scoreMin));
					successCallback({ scoreMin: -1, scoreMax: 1 });
				}, errorCallback);
			} else
				successCallback({ scoreMin: -1, scoreMax: 1 });
		} else
	        this.store.getGlobalStats( successCallback, errorCallback );
    },

	// ----------------------------------------------------------------

	_getHistogram: function(chr) {
		dojo.byId('coge-hist').innerHTML='loading histogram <img src="picts/ajax-loader.gif">';
		dojo.xhrGet({
			url: api_base_url + '/experiment/' + this._track.config.coge.id + '/histogram/' + chr,
			handleAs: 'json',
			load: function(data) {
				if (data.error) {
					coge_plugin.error('Search', data);
					if (coge_xyplot._search_dialog)
						coge_xyplot._search_dialog.hide();
					return;
				}
				if (chr != coge_xyplot._hist_chr) {
					var hist = dojo.byId('coge-hist');
					if (hist) {
						coge_xyplot._hist_chr = chr;
						dojo.empty(hist);
						coge_xyplot._brush = chart(d3.select('#coge-hist'), data.first, data.gap, data.counts);
					}
				}
			},
			error: function(data) {
				coge_plugin.error('Search', data);
				if (coge_xyplot._search_dialog)
					coge_xyplot._search_dialog.hide();
			}
		});
	},

	// ----------------------------------------------------------------

    getRegionStats: function( query, successCallback, errorCallback ) {
        return this._getRegionStats.apply( this, arguments );
    },

    _getRegionStats: function( query, successCallback, errorCallback ) {
        var thisB = this;
        var cache = thisB._regionStatsCache = thisB._regionStatsCache || new LRUCache({
            name: 'regionStatsCache',
            maxSize: 1000, // cache stats for up to 1000 different regions
            sizeFunction: function( stats ) { return 1; },
            fillCallback: function( query, callback ) {
                //console.log( '_getRegionStats', query );
                var s = {
                    scoreMax: -Infinity,
                    scoreMin: Infinity,
                    scoreSum: 0,
                    scoreSumSquares: 0,
                    basesCovered: query.end - query.start,
                    featureCount: 0
                };
                thisB.getFeatures( query,
                                  function( feature ) {
                                      var score = feature.get('score') || 0;
										if (thisB.config.coge.transform == 'Log10') {
											if (score >= 0)
												score = log10(Math.abs(score)+1);
											else
												score = -1*log10(Math.abs(score)+1);
										}
										else if (thisB.config.coge.transform == 'Log2') {
											if (score >= 0)
												score = log2(Math.abs(score)+1);
											else
												score = -1*log2(Math.abs(score)+1);
										}
										else if (thisB.config.coge.transform == 'Inflate')
											score = score > 0 ? 1 : -1;
                                      s.scoreMax = Math.max( score, s.scoreMax );
                                      s.scoreMin = Math.min( score, s.scoreMin );
                                      s.scoreSum += score;
                                      s.scoreSumSquares += score*score;
                                      s.featureCount++;
                                  },
                                  function() {
                                      s.scoreMean = s.featureCount ? s.scoreSum / s.featureCount : 0;
                                      s.scoreStdDev = thisB._calcStdFromSums( s.scoreSum, s.scoreSumSquares, s.featureCount );
                                      s.featureDensity = s.featureCount / s.basesCovered;
                                      //console.log( '_getRegionStats done', s );
                                      callback( s );
                                  },
                                  function(error) {
                                      callback( null, error );
                                  }
                                );
            }
         });

         cache.get( query,
                    function( stats, error ) {
                        if( error )
                            errorCallback( error );
                        else
                            successCallback( stats );
                    });

    },

	// ----------------------------------------------------------------

	 _getScaling: function( viewArgs, successCallback, errorCallback ) {

		this._getScalingStats( viewArgs, dojo.hitch(this, function( stats ) {

			//calculate the scaling if necessary
			if( ! this.lastScaling || ! this.lastScaling.sameStats( stats ) || this.trackHeightChanged ) {

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
				this.trackHeightChanged=false; //reset flag
			}

			successCallback( this.lastScaling );
		}), errorCallback );
	},

	// ----------------------------------------------------------------

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

	// ----------------------------------------------------------------

	_search_track: function() {
		var div = dojo.byId('coge-track-search');
		var search = {};

		if (dojo.byId('max').checked)
			search.type = 'max';
		else if (dojo.byId('min').checked)
			search.type = 'min';
		else {
			search.type = 'range';
			search.gte = dojo.byId('hist_from').value;
			search.lte = dojo.byId('hist_to').value;
			if (!search.gte && !search.lte) {
				coge_plugin.error('Unspecified Range', 'Please drag on the histogram or enter the range of values you wish to search for');
				return;
			}
		}
		var ref_seq = dojo.byId('coge_ref_seq');
		if (ref_seq.selectedIndex > 0)
			search.chr = ref_seq.options[ref_seq.selectedIndex].innerHTML;

		dojo.empty(div);
		div.innerHTML = '<img src="picts/ajax-loader.gif">';
		this._track.config.coge.search = search;
		dojo.xhrGet({
			url: api_base_url + '/search/query/' + this._track.config.coge.id + '?' + coge_plugin.search_to_params(search),
			handleAs: 'json',
			load: dojo.hitch(this, function(data) {
				if (this._search_dialog)
					this._search_dialog.hide();
				if (data.length == 0)
					coge_plugin.error('Search', 'Search returned zero hits');
				else
					coge_plugin.new_search_track(this._track, data);
			}),
			error: dojo.hitch(this, function(data) {
				if (this._search_dialog)
					this._search_dialog.hide();
				coge_plugin.error('Search', data);
			})
		});
	},

	// ----------------------------------------------------------------

	_search_track_dialog: function(track) {
		this._track = track;
		var content = '<div id="coge-track-search"><table align="center"><tr><td>Chromosome:</td><td>';
		content += coge_plugin.build_chromosome_select('Any', 'coge_xyplot._getHistogram(this.options[this.selectedIndex].text)');
		content += '</td></tr>';
		content += '<tr><td>Values:</td><td style="white-space:nowrap"><input id="max" type="radio" name="type" checked="checked"> max</td></tr>';
		content += '<tr><td></td><td style="white-space:nowrap"><input id="min" type="radio" name="type"> min</td></tr>';
		content += '<tr><td></td><td style="white-space:nowrap" valign="top"><input id="range" type="radio" name="type"> range: from <input id="hist_from" size="4" onfocus="dojo.byId(\'range\').checked=true;" /> to <input id="hist_to" size="4" onfocus="dojo.byId(\'range\').checked=true;" /></td></tr>';
		content += '<tr><td></td><td><div id="coge-hist"></div></td></tr></table>';
		content += coge_plugin.build_buttons('coge_xyplot._search_track()', 'coge_xyplot._search_dialog.hide()');
		content += '</div>';
		this._search_dialog = new Dialog({
			title: 'Search Track',
			content: content,
			onHide: function() {
				this.destroyRecursive();
				coge_xyplot._search_dialog = null;
				coge_xyplot._hist_chr = null;
			}
		});
		this._getHistogram(coge_plugin.browser.refSeq.name);
		this._search_dialog.show();
	},

	// ----------------------------------------------------------------

	_setTransform: function(track, transform) {
		this.config.coge.transform = transform;
		track.changed();
	},

	// ----------------------------------------------------------------

	_showPixelValue: function( scoreDisplay, score ) {
		var scoreType = typeof score;
		if( scoreType == 'number' ) {
			// display the score with only 6 significant digits, avoiding
			// most confusion about the approximative properties of
			// IEEE floating point numbers parsed out of BigWig files
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

	// ----------------------------------------------------------------

	_trackMenuOptions: function() {
		var options = this.inherited(arguments);
		var track = this;
		var config = this.config;

		options.push({ type: 'dijit/MenuSeparator' });

		if (config.coge.type == 'experiment')
			options.push({
				label: 'ExperimentView',
				onClick: function(){window.open( 'ExperimentView.pl?eid=' + config.coge.id );}
			});
		else if (config.coge.type == 'notebook')
			options.push({
				label: 'NotebookView',
				onClick: function(){window.open( 'NotebookView.pl?nid=' + config.coge.id );}
			});

		options.push.apply(
			options,
			[
				{
					label: 'Autoscale',
					type: 'dijit/CheckedMenuItem',
					checked: this.config.autoscale == 'local',
					onClick: function(event) {
						track.config.autoscale = this.checked ? 'local' : 'global';
						track.changed();
					}
				},
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
								style: { width: '230px', },
								items: track.config.coge.experiments || [{id: track.config.coge.id, name: track.config.coge.name}],
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

		if (config.coge.type != 'search' && config.coge.type != 'notebook')
            options.push({
                label: 'Find Data that Overlaps Features',
                onClick: function(){coge_plugin.features_overlap_search_dialog(track, 'Data');}
            });
			options.push({
				label: 'Search Experiment Data',
				onClick: function(){coge_xyplot._search_track_dialog(track);}
			});

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
										config.coge.transform = null;
										track.changed();
									}
								},
								{   label: 'Average',
									onClick: function(event) {
										config.coge.transform = 'Average';
										track.changed();
									}
								},
								{   label: 'Difference',
									onClick: function(event) {
										config.coge.transform = 'Difference';
										track.changed();
									}
								},
								{   label: 'Log10',
									onClick: function(event) {
										config.coge.transform = 'Log10';
										track.changed();
									}
								},
								{   label: 'Log2',
									onClick: function(event) {
										config.coge.transform = 'Log2';
										track.changed();
									}
								},
								{   label: 'Normalize',
									onClick: function(event) {
										config.coge.transform = 'Normalize';

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
										config.coge.transform = null;
										track.changed();
									}
								},
								{   label: 'Inflate',
									onClick: function(event) {
										config.coge.transform = 'Inflate';
										track.changed();
									}
								},
								{   label: 'Log10',
									onClick: function(event) {
										config.coge.transform = 'Log10';
										track.changed();
									}
								},
								{   label: 'Log2',
									onClick: function(event) {
										config.coge.transform = 'Log2';
										track.changed();
									}
								},
								{   label: 'Normalize',
									onClick: function(event) {
										config.coge.transform = 'Normalize';
										track.changed();
									}
								}
							]
						}
					]
				);
		}

		if (config.coge.type != 'notebook')
			options.push({
				label: 'Export Track Data',
				onClick: function(){coge_plugin.export_dialog(track);}
			});
		if (config.coge.search && config.coge.data_type == 1) {
			options.push({
				label: 'Convert to Marker Track',
				onClick: function(){coge_plugin.convert_to_marker_dialog(track)}
			});
			options.push({
				label: 'Save Results as New Experiment',
				onClick: function(){coge_plugin.save_as_experiment_dialog(track)}
			});
		}
		if (config.coge.type == 'merge') {
			options.push({
				label: 'Create New Notebook with Merged Tracks',
				onClick: function(){coge_plugin.create_notebook_dialog(track)}
			});
			if (config.coge.keys.length > 1) { // can't remove last track
				var tracks = [];
				config.coge.keys.forEach(function(key){
					tracks.push({
						label: key,
						onClick: function(event) {
							var index = config.coge.keys.indexOf(key);
							config.coge.keys.splice(index, 1);
							config.coge.eids.splice(index, 1);
							track.browser.getStore(config.store, function(store){
								store.baseUrl = store.config.baseUrl = api_base_url + '/experiment/' + config.coge.eids.join(',') + '/';
							});
							track.changed();
							track.makeTrackMenu();
						}
					});
				});
				options.push({
					label: 'Remove Track',
					type: 'dijit/DropDownMenu',
					children: tracks
				});
			}
		}

		return options;
	},

	// ----------------------------------------------------------------

	updateStaticElements: function( coords ) {
		this.inherited( arguments );
		coge_plugin.adjust_nav(this.config.coge.id)
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

function nbspPad(s, padLength) {
	for (var i = s.length;  i < padLength;  i++) {
		s += '&nbsp;';
	}
	return s;
}

function sortByScore(a,b) { // sort features by score
	return Math.abs( b.feature.get('score') ) - Math.abs( a.feature.get('score') );
}

function chart(div, first, gap, counts) {
	var margin = {top: 10, right: 10, bottom: 20, left: 10},
		x = d3.scale.linear().range([0, 200]).domain([0, first + gap * counts.length]),
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
		g.select("#clip rect")
			.attr("x", x(extent[0]))
			.attr("width", x(extent[1]) - x(extent[0]));
		dojo.byId('hist_from').value = parseFloat(Math.round(extent[0] * 100) / 100).toFixed(2);
		dojo.byId('hist_to').value = parseFloat(Math.round(extent[1] * 100) / 100).toFixed(2);
	});

	brush.on("brushend.chart", function() {
		if (brush.empty()) {
			var div = d3.select(this.parentNode.parentNode.parentNode);
			div.select(".title a").style("display", "none");
			div.select("#clip rect").attr("x", null).attr("width", "100%");
			dojo.byId('hist_from').value = '';
			dojo.byId('hist_to').value = '';
		}
	});
	d3.rebind(chart, brush, "on");
	return brush;
}
