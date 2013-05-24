define(['dojo/_base/declare',
        'dojo/_base/array',
        'dojo/dom-construct',
        'dijit/layout/ContentPane',
        'dojo/dnd/Source',
        'dojo/fx/easing',
        'dijit/form/TextBox'
       ],
       function( declare, array, dom, ContentPane, dndSource, animationEasing, dijitTextBox ) {
return declare( 'JBrowse.View.TrackList.CoGe', null,

    /** @lends JBrowse.View.TrackList.CoGe.prototype */
    {
    /**
     * CoGe drag-and-drop track selector.
     * @constructs
     */
    constructor: function( args ) {
        this.browser = args.browser;

        // make the track list DOM nodes and widgets
        this.createTrackList( args.browser.container, args.trackConfigs );

        // maintain a list of the HTML nodes of filtered tracks
        this.filteredNodes = {};        
        
        // subscribe to drop events for tracks being DND'ed
        this.browser.subscribe(
            "/dnd/drop",
            dojo.hitch( this,
                        function( source, nodes, copy, target ) {
			            	if ( source == target) { // both source and target
			            		console.log('source = target');
			            		return;
			            	}
			            	
			            	var isSource = this.trackListWidgets.indexOf(source) != -1;
			            	var isTarget = this.trackListWidgets.indexOf(target) != -1;
			            	
			            	if( isSource && !isTarget ) { // source
                            	console.log('/dnd/drop/source');
                            	// get the configs from the tracks being dragged in
	                            var confs = dojo.filter(
	                                dojo.map( nodes, function(n) {
	                                              return source.map[n.id].data;
	                                        }),
	                                function(c) {return c;}
	                            );
                            	
	                            // highlight track to show it is enabled
	                            this.dndDrop = true;
	                            this.browser.publish( '/jbrowse/v1/v/tracks/show', confs ); // mdb: why not just call setTrackActive directly?
	                            this.dndDrop = false;
                            }
			            	
                            if( this.trackListWidgets.indexOf(target) != -1 ) { // target
                            	console.log('/dnd/drop/target');
                            	
	                            // get the configs from the tracks being dragged in
	                            var confs = dojo.filter(
	                                dojo.map( nodes, function(n) {
	                                              return n.track && n.track.config;
	                                        }),
	                                function(c) {return c;}
	                            );
	                            
	                            // return if no confs; whatever was dragged here probably wasn't a track
	                            if( ! confs.length )
	                                return;

	                            
	                            // un-highlight track to show it is disabled
	                            this.dndDrop = true;
	                            this.browser.publish( '/jbrowse/v1/v/tracks/hide', confs ); // mdb: why not just call setTrackInactive directly?
	                            this.dndDrop = false;
            				}
                        }
                      ));

        // subscribe to commands coming from the the controller
        this.browser.subscribe( '/jbrowse/v1/c/tracks/show',
                                dojo.hitch( this, 'setTracksActive' ));
        this.browser.subscribe( '/jbrowse/v1/c/tracks/hide',
                                dojo.hitch( this, 'setTracksInactive' ));
        this.browser.subscribe( '/jbrowse/v1/c/tracks/new',
                                dojo.hitch( this, 'addTracks' ));
        this.browser.subscribe( '/jbrowse/v1/c/tracks/replace',
                                dojo.hitch( this, 'replaceTracks' ));
        this.browser.subscribe( '/jbrowse/v1/c/tracks/delete',
                                dojo.hitch( this, 'deleteTracks' ));
    },

    addTracks: function( trackConfigs ) {
    	console.log('addTracks');
//        // note that new tracks are, by default, hidden, so we just put them in the list
//        this.trackListWidget.insertNodes(
//            false,
//            trackConfigs
//        );
//
//        this._blinkTracks( trackConfigs );
    },

    // mdb: unused now
    replaceTracks: function( trackConfigs ) {
    	console.log('replaceTracks');
//        // for each one
//        array.forEach( trackConfigs, function( conf ) {
//            // figure out its position in the genome view and delete it
//            var oldNode = this.inactiveTrackNodes[ conf.label ];
//            if( ! oldNode )
//                return;
//            delete this.inactiveTrackNodes[ conf.label ];
//
//            this.trackListWidget.delItem( oldNode.id );
//            if( oldNode.parentNode )
//                oldNode.parentNode.removeChild( oldNode );
//
//           // insert the new track config into the trackListWidget after the 'before'
//           this.trackListWidget.insertNodes( false, [conf], false, oldNode.previousSibling );
//       },this);
    },

    /** @private */
    createTrackList: function( renderTo, trackConfigs ) {
        var trackPane = this.pane = dojo.create(
            'div',
            { id: 'trackPane',
              style: { width: '12em' }
            },
            renderTo
        );

        // splitter on right side
        var trackWidget = new ContentPane({region: "right", splitter: true}, trackPane);

        this.div = dojo.create(
            'div',
            { id: 'tracksAvail',
              className: 'container handles',
              style: { width: '100%', height: '100%', overflowX: 'hidden', overflowY: 'auto' },
              //innerHTML: '<h2>Available Tracks</h2>',
              //onclick: dojo.hitch( this, function() { this.trackListWidget.selectNone(); } )
            },
            trackPane
        );

        // create text filter input
        this._createTextFilter();
        this._updateTextFilterControl();

        // create a DnD source for sequence, annotation, etc.
        this.trackListWidgets = [];
        
        var temp = this._createDnDSource();
        temp.insertNodes(
            	false,
            	trackConfigs.filter( function(e) {
            		return ( !e.coge.type ||
            				 (e.coge.type == 'sequence' || e.coge.type == 'annotation') );
            	})
            );
        this.trackListWidgets.push( temp );
        
        // create a DnD source for each notebook
        var notebooks = trackConfigs.filter( function(e) {
    		return (e.coge.type && e.coge.type == 'notebook');
    	});
    	var experiments = trackConfigs.filter( function(e) {
    		return (e.coge.type && e.coge.type == 'experiment');
    	});
    	var that = this;
        notebooks.forEach( function(n) {
    		temp = that._createDnDSource();
    		temp.insertNodes(
                	false,
                	[n].concat(experiments.filter( function(e) {
                		return e.coge.notebooks && dojo.indexOf(e.coge.notebooks, n.coge.id) != -1;
                	}))
                );
    		that.trackListWidgets.push( temp );
    	});
        
        // create a DnD source for remaining experiments
        temp = this._createDnDSource();
        temp.insertNodes(
            	false,
            	experiments.filter( function(e) {
            		return !e.coge.notebooks;
            	})
            );
        this.trackListWidgets.push( temp );
        
        return this.div;
    },
    
    _createDnDSource: function() {
        var div = dojo.create( 'div', {}, this.div );
    	
    	return new dndSource( // modifies div to be DnD-able
            div,
            {
                accept: ["track"], // accepts only tracks
                withHandles: false,
                copyOnly: true,
                checkAcceptance: function( source, nodes ) {
                	console.log('checkAcceptance');
                	var accept = true;
                	var target = this;
                	nodes.forEach( function (n) {
                		var type = source.map[n.id].data.coge.type;
	                	if (!type || type != 'experiment' || n.id in target.map) {
	                		accept = false;
	                	}
                	});
                	return accept;
                },
//	            onDrop: function( source, nodes, copy ) {
//	                console.log('onDrop');
//	                return true;
//	            },
//                onDropExternal: dojo.hitch( this, function( source, nodes, copy ) {
//                	console.log('onDropExternal');
//                	return this.inherited(arguments);
//                }),
//                onDropInternal: dojo.hitch( this, function( source, nodes, copy ) {
//	            	console.log('onDropInternal');
//	            	return this.inherited(arguments);
//                }),
                creator: dojo.hitch( this, function( trackConfig, hint ) {
                	console.log('creator');
                	var id = trackConfig.coge.type + '_' + trackConfig.coge.id;
                	var name = ('name' in trackConfig.coge ? trackConfig.coge.name : trackConfig.key);
                	var node = this._createLabelNode( trackConfig );
                	if (trackConfig.coge.type == 'notebook') {
                		node.innerHTML = '<img style="padding-right:3px" src="picts/notebook-icon-small.png">' + 
                			'<span class="tracklist-text">' + name + '</span>';
                	}
                	else if (trackConfig.coge.type == 'experiment') {
                		node.innerHTML = (trackConfig.coge.notebooks ? '<span style="padding-left:25px;"></span>' : '') 
                              	+ '<img src="picts/testtube-icon-small.png" />' + ' ' + 
                              	'<span class="tracklist-text">' + name + '</span>';
                	}
                	else {
                		node.innerHTML = '<span class="tracklist-text">' + name + '</span>';
                	}
                	node.id = id;
                	
                    // in the list, wrap the list item in a container for border drag-insertion-point monkeying
                    dojo.connect( node, "click", dojo.hitch(this, function() {
                    	//console.log('click ' + node.id);
                    	if (dojo.hasClass(node, 'selected')) {
                    		this.browser.publish( '/jbrowse/v1/v/tracks/hide', [trackConfig] );
                    	}
                    	else {
                    		this.browser.publish( '/jbrowse/v1/v/tracks/show', [trackConfig] );
                    	}
                    }));
                    
                    var container = dojo.create( 'div', { className: 'tracklist-container' });
                    
                    if (trackConfig.coge.type == 'experiment') { 
                    	if (trackConfig.coge.notebooks) {
                    		dojo.addClass(container, 'collapsed');
                    	}
                    }
                    else if (trackConfig.coge.type == 'notebook') {
                		var button = dom.create(
                			'button',
                    	    {	innerHTML: '+',
                				style: { float: 'right' }
                    	    }, 
                    	    container
                    	);
                    	dojo.connect( button, "click", dojo.hitch(this, function() {
                            var expanded = (button.innerHTML == '-');
                            if (expanded) {
                            	button.innerHTML = '+';
                            	var children = div.children;
                            	for (var i = 1;  i < children.length;  i++) {
		                            dojo.addClass(children[i], 'collapsed');
                            	}
                            }
                            else {
                            	button.innerHTML = '-';
                            	var children = div.children;
                            	for (var i = 1;  i < children.length;  i++) {
		                            dojo.removeClass(children[i], 'collapsed');
                            	}
                            }
                        }));
                    }
                    
                    container.appendChild(node);
                    container.id = dojo.dnd.getUniqueId();
                    return {node: container, data: trackConfig, type: ["track", trackConfig.coge.type]};
                })
            }
        ); 
    },
    
    _createTextFilter: function( ) {
        this.textFilterDiv = dom.create( 'div', {
            className: 'textfilter',
            style: {
                width: '100%',
                position: 'relative',
                overflow: 'hidden'
            }
        }, this.div );
        
		this.textFilterInput = dom.create(
			'input',
			{	type: 'text',
				style: {
					paddingLeft: '18px',
					height: '16px',
					width: '80%'
				},
				placeholder: 'filter by text',
				onkeypress: dojo.hitch( this, function( evt ) {
					if( this.textFilterTimeout )
						window.clearTimeout( this.textFilterTimeout );
					this.textFilterTimeout = window.setTimeout(
						dojo.hitch( this, function() {
						      this._updateTextFilterControl();
						      this._textFilter( this.textFilterInput.value, this.filteredNodes );
						  }),
						500
					);
					this._updateTextFilterControl();
					
					evt.stopPropagation();
				})
			},
			dom.create('div',{ style: 'overflow: show;' }, this.textFilterDiv )
		);
		
		this.textFilterClearButton = dom.create('div', {
			className: 'jbrowseIconCancel',
			onclick: dojo.hitch( this, function() {
				this._clearTextFilterControl();
				this._textFilter( this.textFilterInput.value, this.filteredNodes );
			}),
			style: {
				position: 'absolute',
				left: '4px',
				top: '2px'
			}
		}, this.textFilterDiv );
    },
    
    _createLabelNode: function( trackConfig ) {
    	return dojo.create(
    				'div',
	                { className: 'coge-tracklist-label coge-' + trackConfig.coge.type,
	                  title: capitalize(trackConfig.coge.type) + " id" + trackConfig.coge.id + "\n" +
	                  		 "Name: " + trackConfig.coge.name + "\n" +
	                  		 "Description: " + trackConfig.coge.description + "\n" +
	                  		 (trackConfig.coge.annotations ?
	                  				trackConfig.coge.annotations.
	                  				map(
                  						function(a) {
                  							return a.type + ': ' + a.text
                  						}
	                  				)
	                  				.join("\n")
	                  				: '') +
	                  		 "\n\nDrag or click to activate"
	                }
	        	);
    },

    _textFilter: function( text, filteredNodes ) {
        if( text && /\S/.test(text) ) {
            text = text.toLowerCase();
            dojo.query( '.tracklist-text', this.div )
                .forEach( function( labelNode, i ) {
                	var container = labelNode.parentNode.parentNode;
                    if( labelNode.innerHTML.toLowerCase().indexOf( text ) != -1 ) {
                        dojo.removeClass( container, 'collapsed');
                        delete filteredNodes[container.id];
                    } 
                    else if (!dojo.hasClass(container, 'collapsed')) { // check if already hidden in collapsed notebook
                        dojo.addClass( container, 'collapsed');
                        filteredNodes[container.id] = container;
                    }
                 });
        } 
        else { // empty string, show all
//            dojo.query( '.tracklist-container.collapsed', this.div )
//                .removeClass('collapsed');
        	for (var id in filteredNodes) {
        		dojo.removeClass(filteredNodes[id], 'collapsed');
        	}
        	filteredNodes = {};
        }
    },

   /**
    * Clear the text filter control input.
    * @private
    */
    _clearTextFilterControl: function() {
        this.textFilterInput.value = '';
        this._updateTextFilterControl();
    },
    /**
     * Update the display of the text filter control based on whether
     * it has any text in it.
     * @private
     */
    _updateTextFilterControl: function() {
        if( this.textFilterInput.value.length )
            dojo.removeClass( this.textFilterDiv, 'dijitDisabled' );
        else
            dojo.addClass( this.textFilterDiv, 'dijitDisabled' );
    },

    /**
     * Given an array of track configs, update the track list to show
     * that they are turned on.
     */
    setTracksActive: function( /**Array[Object]*/ trackConfigs ) {
    	console.log('setTracksActive ');
        dojo.query( '.coge-tracklist-label', this.div )
	        .forEach( function( labelNode, i ) {
	        	trackConfigs.forEach( function (trackConfig) {
	        		var trackId = trackConfig.coge.type + '_' + trackConfig.coge.id;
	        		if (labelNode.id == trackId) {
	    				dojo.addClass(labelNode, 'selected');
	        			if (dojo.hasClass(labelNode, 'coge-experiment')) {
	        				var id = trackConfig.coge.id;
	        				var color = getFeatureColor(id);
	        				dojo.style(labelNode, 'background', color);
	        			}
	        			else {
	        				dojo.style(labelNode, 'background', 'lightgray');
	        			}
	        		}
	        	});
	        });
    },

    // mdb: unused now
    deleteTracks: function( /**Array[Object]*/ trackConfigs ) { // mdb: unused now ...?
//    	console.log('deleteTracks');
//        // remove any tracks in our track list that are being set as visible
//        array.forEach( trackConfigs || [], function( conf ) {
//            var oldNode = this.inactiveTrackNodes[ conf.label ];
//            if( ! oldNode )
//                return;
//            delete this.inactiveTrackNodes[ conf.label ];
//
//            if( oldNode.parentNode )
//                oldNode.parentNode.removeChild( oldNode );
//
//            this.trackListWidget.delItem( oldNode.id );
//        },this);
    },

    /**
     * Given an array of track configs, update the track list to show
     * that they are turned off.
     */
    setTracksInactive: function( /**Array[Object]*/ trackConfigs ) {
    	console.log('setTracksInactive');
        dojo.query( '.coge-tracklist-label', this.div )
	        .forEach( function( labelNode, i ) {
	        	trackConfigs.forEach(function (trackConfig) {
	        		var trackId = trackConfig.coge.type + '_' + trackConfig.coge.id;
	        		if (labelNode.id == trackId) {
	    				dojo.style(labelNode, 'background', '');
	        			dojo.removeClass(labelNode, 'selected');
	        		}
	        	});
	        });  
    	
        // remove any tracks in our track list that are being set as visible
//        if( ! this.dndDrop ) {
//            var n = this.trackListWidget.insertNodes( false, trackConfigs );
//
//            // blink the track(s) that we just turned off to make it
//            // easier for users to tell where they went.
//            // note that insertNodes will have put its html element in
//            // inactivetracknodes
//            this._blinkTracks( trackConfigs );
//        }
    },
    
    _blinkTracks: function( trackConfigs ) {
    	console.log('_blinkTracks');
    	
        // scroll the tracklist all the way to the bottom so we can see the blinking nodes
//        this.trackListWidget.node.scrollTop = this.trackListWidget.node.scrollHeight;
//
//        array.forEach( trackConfigs, function(c) {
//            var label = this.inactiveTrackNodes[c.label].firstChild;
//            if( label ) {
//                dojo.animateProperty({
//                                         node: label,
//                                         duration: 400,
//                                         properties: {
//                                             backgroundColor: { start: '#DEDEDE', end:  '#FFDE2B' }
//                                         },
//                                         easing: animationEasing.sine,
//                                         repeat: 2,
//                                         onEnd: function() {
//                                             label.style.backgroundColor = null;
//                                         }
//                                     }).play();
//            }
//        },this);
    },
    
    /**
     * Make the track selector visible.
     * This does nothing for the Simple track selector, since it is always visible.
     */
    show: function() {
    },

    /**
     * Make the track selector invisible.
     * This does nothing for the Simple track selector, since it is always visible.
     */
    hide: function() {
    },

    /**
     * Toggle visibility of this track selector.
     * This does nothing for the Simple track selector, since it is always visible.
     */
    toggle: function() {
    }

});
});

function getFeatureColor(id) { //FIXME: dup'ed in MultiXYPlot.js
	return '#' + ((((id * 1234321) % 0x1000000) | 0x444444) & 0xe7e7e7 ).toString(16); 
}

function capitalize( string ) { //FIXME: doesn't go here, extend String class instead
    return string.charAt(0).toUpperCase() + string.slice(1);
}
