define( [
            'dojo/_base/declare',
            'dojo/_base/array',
            'dojo/_base/Color',
            'dojo/on',
            'JBrowse/View/Track/WiggleBase',
            'JBrowse/View/Track/YScaleMixin',
            'JBrowse/Util',
            './_Scale'
        ],
        function( declare, array, Color, on, WiggleBase, YScaleMixin, Util, Scale ) {

var XYPlot = declare( [WiggleBase, YScaleMixin], // mdb: this file is a copy of XYPlot, extend that class instead?

/**
 * Wiggle track that shows data with an X-Y plot along the reference.
 *
 * @lends JBrowse.View.Track.Wiggle.XYPlot
 * @extends JBrowse.View.Track.WiggleBase
 */
{
    _defaultConfig: function() {
//    	console.log('_defaultConfig');
        return Util.deepUpdate(
            dojo.clone( this.inherited(arguments) ),
            {
                style: {
                    pos_color: 'blue',
                    neg_color: 'red',
                    origin_color: '#888',
                    variance_band_color: 'rgba(0,0,0,0.3)'
                }
            }
        );
    },

    _getScaling: function( successCallback, errorCallback ) {
//    	console.log('_getScaling');
    	
        this._getScalingStats( dojo.hitch(this, function( stats ) {

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
    },

    /**
     * Draw a set of features on the canvas.
     * @private
     */
    _drawFeatures: function( scale, leftBase, rightBase, block, canvas, features, featureRects, dataScale ) {
//    	console.log('_drawFeatures');
        var context = canvas.getContext('2d');
        var canvasHeight = canvas.height;
        var toY = dojo.hitch( this, function( val ) {
           return canvasHeight * ( 1-dataScale.normalize.call(this, val) );
        });
        var originY = toY( dataScale.origin );

        var disableClipMarkers = this.config.disable_clip_markers;
        
        // sort features by score
        var sorted = [];
        dojo.forEach( features, function(f,i) {
        	sorted.push({ feature: f, featureRect: featureRects[i] });
        });
        sorted.sort( sortByScore );
        
        // compute average scores
        var sum = [];
        var count = [];
        var width = [];
        dojo.forEach( features, function(f,i) {
        	var score = f.get('score');
        	var l = featureRects[i].l;
        	var w = featureRects[i].w;
        	sum[l] = l in sum ? sum[l] + score : score;
        	count[l] = l in count ? count[l] + 1 : 1;
        	width[l] = l in width ? Math.max(width[l], w) : w;
        });
        
        if (this.config.showAverage) {
        	sum.forEach( function(x,l) {
        		var avg = sum[l]/count[l];
        		var height = toY( avg );
        		context.fillStyle = 'gray';
        		
        		if( height <= canvasHeight ) { // if the rectangle is visible at all
        			if( height <= originY ) { // bar goes upward
        				context.fillRect( l, height, width[l], originY-height+1);
        			}
        			else { // bar goes downward
        				context.fillRect( l, originY, width[l], height-originY+1 );
        			}
        		}
        	});
        }
        else {
            dojo.forEach( sorted, function(pair,i) {
            	var f = pair.feature;
                var fRect = pair.featureRect;//featureRects[i];
                var score = f.get('score');
                fRect.t = toY( score );

                // draw the background color if we are configured to do so
//                if( fRect.t >= 0 ) {
//                    var bgColor = this.getConfForFeature('style.bg_color', f );
//                    if( bgColor ) {
//                        context.fillStyle = bgColor;
//                        context.fillRect( fRect.l, 0, fRect.w, canvasHeight );
//                    }
//                }
                
                if( fRect.t <= canvasHeight ) { // if the rectangle is visible at all
                	var id = f.get('id');
                	context.fillStyle = this._getFeatureColor(id);

                	if( fRect.t <= originY ) {
                        // bar goes upward
//                        context.fillStyle = this.getConfForFeature('style.pos_color',f);
                        context.fillRect( fRect.l, fRect.t, fRect.w, originY-fRect.t+1);
//                        if( !disableClipMarkers && fRect.t < 0 ) { // draw clip marker if necessary
//                            context.fillStyle = this.getConfForFeature('style.clip_marker_color',f) || this.getConfForFeature('style.neg_color',f);
//                            context.fillRect( fRect.l, 0, fRect.w, 2 );
//                        }
                    }
                    else {
                        // bar goes downward
//                        context.fillStyle = this.getConfForFeature('style.neg_color',f);
                        context.fillRect( fRect.l, originY, fRect.w, fRect.t-originY+1 );
//                        if( !disableClipMarkers && fRect.t >= canvasHeight ) { // draw clip marker if necessary
//                            context.fillStyle = this.getConfForFeature('style.clip_marker_color',f) || this.getConfForFeature('style.pos_color',f);
//                            context.fillRect( fRect.l, canvasHeight-3, fRect.w, 2 );
//                        }
                    }
                }
            }, this );        	
        }

    },

    /**
     * Draw anything needed after the features are drawn.
     */
    _postDraw: function( scale, leftBase, rightBase, block, canvas, features, featureRects, dataScale ) {
//    	console.log('_postDraw');
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
                                          pos+'Ïƒ');
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
        
        //Êdraw average line - mdb added 5/15/13 ... not working
//        if (this.config.coge.showAverage) {
//        	var originY = toY( dataScale.origin );
//    		var sum = [];
//    		var count = [];
//    		featureRects.forEach( function (rect) {
//				var height = originY-rect.t+1;
//				sum[rect.l] = rect.l in sum ? sum[rect.l] + height : height;
//				count[rect.l] = rect.l in count ? count[rect.l] + 1 : 1;
//    		});
//    		featureRects.forEach( function (rect) {
//    			var avg = sum[rect.l]/count[rect.l];
//	    		context.fillStyle = 'red';
//	        	context.fillRect( rect.l-5, avg, rect.w+10, 1 );
//    		});
//        }

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
	            var id = f.feature.get('id');
	            var name = this._getFeatureName(id);
	            var color = this._getFeatureColor(id);
	            for( var j = Math.round(fRect.l); j < jEnd; j++ ) {
	            	var label = '<div style="background-color:'+color+';">' + nbspPad(score.toString(), 11) + name + '</div>';
	                pixelValues[j] = j in pixelValues ? pixelValues[j] + label : label;
	            }
	        },this);
	        
	        
	        // compute average scores - FIXME dup'ed in _drawFeatures
	        if (this.config.showAverage) {
		        var sum = [];
		        var count = [];
		        var width = [];
		        dojo.forEach( features, function(f,i) {
		        	var score = f.get('score');
		        	var l = featureRects[i].l;
		        	var w = featureRects[i].w;
		        	sum[l] = l in sum ? sum[l] + score : score;
		        	count[l] = l in count ? count[l] + 1 : 1;
		        	width[l] = l in width ? Math.max(width[l], w) : w;
		        });
	        
		        sum.forEach( function(x,l) {
	        		var avg = sum[l]/count[l];
	        		for( var j = Math.round(l); j < l+width[l]; j++ ) {
	                	var label = '<div style="background-color:gray;">' + nbspPad(avg.toPrecision(6).toString(), 11) + 'Average' + '</div>';
	                    pixelValues[j] = j in pixelValues ? pixelValues[j] + label : label;
	                }
	        	});
	        }
        }
        
        return pixelValues;
    }, 
    
    _trackMenuOptions: function() {
    	var options = this.inherited(arguments); 
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
                    { label: 'Show scores on hover',
                      type: 'dijit/CheckedMenuItem',
                      checked: this.config.showHoverScores,
                      onClick: function(event) {
                          track.config.showHoverScores = this.checked;
                          track.changed();
                      }
                    }
                ]
            );
    	
    	if (config.coge.type == 'notebook') {
	    	var track = this;
	        options.push.apply(
	                options,
	                [
	                    { label: 'Show average',
	                      type: 'dijit/CheckedMenuItem',
	                      checked: this.config.showAverage,
	                      onClick: function(event) {
	                          track.config.showAverage = this.checked;
	                          track.changed();
	                      }
	                    }
	                ]
	            );
    	}
    	
    	return options;
    },

    _getFeatureName: function(id) {
    	var cogeConfig = this.config.coge;
    	
    	if (cogeConfig.type == 'experiment') {
    		return cogeConfig.name;
    	}
    	else if (cogeConfig.type == 'notebook') {
	    	var experiments = cogeConfig.experiments || [];
	    	var name;
	    	experiments.forEach( function(e) {
	    		if (e.id == id) {
	    			name = e.name;
	    		}
	    	});
	    	return name;
    	}
    },
    
    _getFeatureColor: function(id) {
    	return '#' + ((((id * 1234321) % 0x1000000) | 0x444444) & 0xe7e7e7 ).toString(16); //FIXME: dup'ed in CoGe.js
    }

});

return XYPlot;
});

function nbspPad(s, padLength) {
	for (var i = s.length;  i < padLength;  i++) {
		s += '&nbsp;';
	}
	return s;
}

function sortByScore(a,b) {
	return Math.abs( b.feature.get('score') ) - Math.abs( a.feature.get('score') );
}
