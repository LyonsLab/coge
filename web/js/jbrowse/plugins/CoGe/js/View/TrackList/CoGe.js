var coge_track_list;
var create_notebook_dialog;

define(['dojo/_base/declare',
        'dojo/_base/array',
        'dojo/dom-construct',
        'dojo/dom-geometry',
        'dojo/aspect',
        'dijit/layout/ContentPane',
        'dojo/dnd/Source',
//        'dojo/fx/easing',
        'dijit/form/TextBox',
        'dojo/mouse',
        'dijit/Dialog',
        "dijit/form/DropDownButton",
//        "dijit/DropDownMenu",
        "dijit/Menu",
        "dijit/MenuItem",
        'JBrowse/View/ConfirmDialog',
        'JBrowse/View/InfoDialog'
       ],
       function( declare, array, dom, domGeom, aspect, ContentPane, dndSource, /*animationEasing,*/ dijitTextBox, mouse, Dialog, DropDownButton, /*DropDownMenu,*/ Menu, MenuItem, ConfirmDialog, InfoDialog ) {
	return declare( 'JBrowse.View.TrackList.CoGe', null,

    /** @lends JBrowse.View.TrackList.CoGe.prototype */
    {
    /**
     * CoGe drag-and-drop track selector.
     * @constructs
     */
    constructor: function( args ) {
    	coge_track_list = this;
        this.browser = args.browser;

        // make the track list DOM nodes and widgets
        this.createTrackList( args.browser.container, args.trackConfigs );

        // maximum tracks that can be added via "+" button
        this.maxTracksToAdd = 20;

        // subscribe to drop events for tracks being DND'ed
        this.browser.subscribe( '/dnd/drop',
        						dojo.hitch( this, 'moveTracks' ));

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

    //----------------------------------------------------------------

    _add_expander: function (container, collapsed) {
    	var button = dom.create('img', {
    			className: (collapsed ? 'coge-tracklist-expandIcon' : 'coge-tracklist-collapseIcon'),
    			src: (collapsed ? 'js/jbrowse/plugins/CoGe/img/arrow-right-icon.png' : 'js/jbrowse/plugins/CoGe/img/arrow-down-icon.png')
        	},
        	container
        );
    	dojo.connect(button, 'click', dojo.hitch(this, function() {
    		dojo.toggleClass(button, 'coge-tracklist-expandIcon coge-tracklist-collapseIcon');
            if (dojo.hasClass(button, 'coge-tracklist-expandIcon')) {
            	button.src = 'js/jbrowse/plugins/CoGe/img/arrow-right-icon.png';
            	// Hide child nodes
            	var n = container.nextSibling;
            	while (n) {
            		dojo.addClass(n, 'collapsed');
            		n = n.nextSibling;
            	}
            }
            else {
            	button.src = 'js/jbrowse/plugins/CoGe/img/arrow-down-icon.png';
            	// Show child nodes
            	var n = container.nextSibling;
            	while (n) {
            		dojo.removeClass(n, 'collapsed');
            		n = n.nextSibling;
            	}
            }
        }));
    },

    //----------------------------------------------------------------

    addTracks: function( trackConfigs ) {
        // note that new tracks are, by default, hidden, so we just put them in the list
        this.trackListWidget.insertNodes(
            false,
            trackConfigs
        );

        this._blinkTracks( trackConfigs );
    },

    //----------------------------------------------------------------

    _capitalize: function (string) {
        return string.charAt(0).toUpperCase() + string.substring(1);
    },

    //----------------------------------------------------------------

   /**
    * Clear the text filter control input.
    * @private
    */
    _clearTextFilterControl: function() {
        this.textFilterInput.value = '';
        this._updateTextFilterControl();
    },

    //----------------------------------------------------------------

    _createDnDSource: function(isStatic) {
        var div = dojo.create( 'div', {}, this.div );

    	return new dndSource( // modifies div to be DnD-able
            div,
            {
                accept: ["track"], // accepts only tracks
                withHandles: false,
                copyOnly: true,
                checkAcceptance: function( source, nodes ) {
                	console.log('checkAcceptance');
                	if (isStatic)
                		return false;
                	var target = this;
                	var nbNode = target.getAllNodes().shift();
                	var nbConfig = target.map[nbNode.id].data.coge;
                	return ( nbConfig.editable && !this._hasLabelNode(div, nodes) );
                },
                creator: dojo.hitch( this, function( trackConfig, hint ) {
                	var node = this._createLabelNode( trackConfig );
                	var coge = trackConfig.coge;

                	node.id = coge.type + coge.id;

                    if (coge.classes) {
                    	dojo.addClass( node, coge.classes.join(' ') );
                    }

                    this._setTrackTitle(trackConfig, node);

                    dojo.connect( node, "click", dojo.hitch(this, function(e) {
                    	//console.log('click ' + node.id + ' ' + e.shiftKey);
                    	if (dojo.hasClass(node, 'selected')) {
                    		this.browser.publish( '/jbrowse/v1/v/tracks/hide', [trackConfig] );
                    	}
                    	else {
                    		this.browser.publish( '/jbrowse/v1/v/tracks/show', [trackConfig] );
                    	}
                    }));

                    // in the list, wrap the list item in a container for border drag-insertion-point monkeying
                    var container = dojo.create( 'div', { className: 'coge-tracklist-container' });
                    dojo.connect(container, "onmouseenter", dojo.hitch(this, function(){
                    	if (this._menu_popped_up())
                    		return;
                    	var b = dijit.byId('coge_track_menu_button');
                    	if (b)
                    		b.destroy();
                    	var deletable = dojo.hasClass(node, 'coge-tracklist-deletable');
                    	var info = dojo.hasClass(node, 'coge-tracklist-info');
                    	if (deletable || info) {
	                		var menu_button = dom.create('div', {id: 'coge_track_menu_button'}, container);
	                        var menu = new Menu();
	                        var type = this._capitalize(coge.type);
	                        if (info)
		                        menu.addChild(new MenuItem({
		                            label: type + ' Info',
		                            onClick: dojo.hitch( this, function() {
		                            	this._info(trackConfig, coge.type);
		                			})
		                        }));
	                        if (deletable)
		                        menu.addChild(new MenuItem({
		                            label: 'Delete ' + type,
		                            onClick: dojo.hitch( this, function() {
		                            	this._delete(trackConfig, coge.type, coge.id, div);
		                			})
		                        }));
	                        var btn = new DropDownButton({
	                            dropDown: menu
	                        }, menu_button);
	                        menu.startup();
	                        btn.startup();
                    	}
                    }));
                    dojo.connect(container, "onmouseleave", dojo.hitch(this, function(){
                    	if (!this._menu_popped_up()) {
                    		var b = dijit.byId('coge_track_menu_button');
                    		if (b)
                    			b.destroy();
                    	}
                    }));

                    if (dojo.hasClass(node, 'coge-tracklist-collapsible'))
                    	this._add_expander(container, coge.collapsed);
                    else if (coge.collapsed)
                    	dojo.addClass(container, 'collapsed');

                    container.appendChild(node);
                    container.id = node.id;//dojo.dnd.getUniqueId();
                    return {node: container, data: trackConfig, type: ["track", coge.type]};
                })
            }
        );
    },

    //----------------------------------------------------------------

    _createLabelNode: function( trackConfig ) {
    	var coge = trackConfig.coge;
    	var title = this._capitalize(coge.type) + " id" + coge.id;
    	if (coge.name)
    		title += "\nName: " + coge.name;
    	if (coge.description)
    		title += "\nDescription: " + coge.description;
    	if (coge.annotations)
    		title += coge.annotations;
    	return dojo.create('div',
            {
				className: 'coge-tracklist-label coge-' + coge.type,
				title:  title
            }
    	);
    },

    //----------------------------------------------------------------

    _createTextFilter: function( trackConfigs, parent ) {
        this.textFilterDiv = dom.create( 'div', {
            className: 'coge-textfilter' //className: 'textfilter', // replace jbrowse styling
        }, parent); //this.div );

        var d = dom.create('div',{ style: 'display:inline;overflow:show;position:relative' }, this.textFilterDiv );
		this.textFilterInput = dom.create(
			'input',
			{	type: 'text',
				id: 'textFilterInput',
				placeholder: 'filter by text',
				onkeypress: dojo.hitch( this, function( evt ) {
					if( this.textFilterTimeout )
						window.clearTimeout( this.textFilterTimeout );
					this.textFilterTimeout = window.setTimeout(
						dojo.hitch( this, function() {
						      this._updateTextFilterControl();
						      this._textFilter( this.textFilterInput.value );
						  }),
						500
					);
					this._updateTextFilterControl();

					evt.stopPropagation();
				})
			},
			d
		);

		this.textFilterCancel = dom.create('div', {
			className: 'jbrowseIconCancel',
			id: 'textFilterCancel',
			onclick: dojo.hitch( this, function() {
				this._clearTextFilterControl();
				this._textFilter( this.textFilterInput.value );
			})
		}, d );

		var menu_button = dom.create('div', {id: 'coge_menu_button'}, this.textFilterDiv);
        var menu = new Menu();
        menu.addChild(new MenuItem({
            label: "Add All Tracks",
            onClick: dojo.hitch( this, function() {
				var configs = this._getVisibleConfigs(this.div, trackConfigs);
				if (configs.length) {
					if (configs.length > this.maxTracksToAdd) {
					    var myDialog = new dijit.Dialog({
					        title: "Warning",
					        content: "There are too many tracks to add (>" + this.maxTracksToAdd + "), please filter them further.",
					        style: "width: 300px"
					    });
					    myDialog.show();
					}
					else {
						this.browser.publish( '/jbrowse/v1/v/tracks/show', configs );
					}
				}
			})
        }));
        menu.addChild(new MenuItem({
            label: "Clear All Tracks",
            onClick: dojo.hitch( this, function() {
				var configs = this._getVisibleConfigs(this.div, trackConfigs);
				if (configs.length) {
					this.browser.publish( '/jbrowse/v1/v/tracks/hide', configs );
				}
			})
        }));
        menu.addChild(new MenuItem({
            label: "Create New Notebook",
            onClick: dojo.hitch( this, function() {
            	create_notebook_dialog = new Dialog({
                    title: "Create New Notebook",
                    content: '<table><tr><td><label>Name:</label></td><td><input id="notebook_name"></td></tr><tr><td><label>Description:</label></td><td><input id="notebook_description"></td></tr><tr><td><label>Restricted:</label></td><td><input type="checkbox" id="notebook_restricted"></td></tr></table><div class="dijitDialogPaneActionBar"><button data-dojo-type="dijit/form/Button" type="button" onClick="create_notebook()">OK</button><button data-dojo-type="dijit/form/Button" type="button" onClick="create_notebook_dialog.hide()">Cancel</button></div>',
                    style: "width: 300px"
                });
            	create_notebook_dialog.show();
            })
        }));
        var btn = new DropDownButton({
            dropDown: menu
        }, menu_button);
        menu.startup();
        btn.startup();

		dom.create('div', {
			style: { clear: 'both' }
		}, this.textFilterDiv );

		this.textFilterLabel = dom.create('div', {
			id: 'textFilterLabel',
		}, this.textFilterDiv );
    },

    //----------------------------------------------------------------

    /** @private */
    createTrackList: function( renderTo, trackConfigs ) {
    	var that = this;

        var trackPane = this.pane = dojo.create(
            'div',
            { id: 'trackPane' },
            renderTo
        );

        // create text filter input
        this._createTextFilter(trackConfigs, trackPane);
        this._updateTextFilterControl();

        // splitter on right side
        var trackWidget = new ContentPane({region: "right", splitter: true}, trackPane);

        this.div = dojo.create(
            'div',
            { id: 'tracksAvail' },
            trackPane
        );

        // create a DnD source for sequence
        this.trackListWidgets = [];
        this.trackListWidgets.push(
        	this._createDnDSource(1).insertNodes(
            	false,
            	trackConfigs.filter( function(e) {
            		var type = e.coge.type;
            		return ( !type ||
            				  type == 'sequence' ||
            				  type == 'gc_content' );
            	})
            )
        );

        // create a DnD source for each feature group
        var feature_groups = trackConfigs.filter( function(fg) {
    		return (fg.coge.type && fg.coge.type == 'feature_group');
    	});
    	var features = trackConfigs.filter( function(f) {
    		return (f.coge.type && f.coge.type == 'features');
    	});
    	feature_groups.forEach( function(fg) {
    		that.trackListWidgets.push(
    			that._createDnDSource().insertNodes(
                	false,
                	[fg].concat(features.filter( function(f) {
                		return f.coge.dataset_id == fg.coge.id;
                	}))
                )
            );
    	});

        // create a DnD source for each notebook
        var notebooks = trackConfigs.filter( function(e) {
    		return (e.coge.type && e.coge.type == 'notebook');
    	});
    	var experiments = trackConfigs.filter( function(e) {
    		return (e.coge.type && e.coge.type == 'experiment');
    	});
        notebooks.forEach( function(n) {
    		that.trackListWidgets.push(
    			that._createDnDSource().insertNodes(
                	false,
                	[n].concat(experiments.filter( function(e) {
                		return e.coge.notebooks && dojo.indexOf(e.coge.notebooks, n.coge.id) != -1;
                	}))
                )
            );
    	});

        // Initialize text filter label
        this._textFilter();

        return this.div;
    },

    //----------------------------------------------------------------

    _delete: function(trackConfig, type, id, div) {
		// Determine if node is inside a notebook
		if (type == 'experiment') {
			var notebookName;
        	var notebookId;
    		dojo.query('.coge-notebook', div).forEach( function(n) {
    			notebookName = n.id;
    			notebookId = notebookName.match(/\d+/);
    		});
    		if (notebookId) { // it's inside a notebook
    			dojo.xhrPut({ // FIXME: make webservice for this
					    url: "NotebookView.pl",
				    putData: {
				    	fname: 'remove_list_item',
				    	lid: notebookId,
				    	item_type: '3', // FIXME: hardcoded type to experiment
				    	item_id: id
				    },
				    handleAs: "json",
				    load: dojo.hitch(this, function(data) {
				    	if (!data) {
				    		new InfoDialog({
			                            title: 'Permission denied',
			                            content: "You don't have permission to do that."
			                        }).show();
				    		return;
				    	}

				    	// Remove node from tracklist
				    	div.removeChild(container);
				    	// Reload track in browser
	                	this.browser.view.tracks.forEach( function(track) {
	                		if (track.config.label == notebookName) {
	                			track.changed();
//	                			track.config.coge.count--;
//	                			this._setTrackTitle(track.config,);
	                		}
	                	});
				    })
			    });
    			return;
    		}
		}

		// Else, it's not inside a notebook
		new ConfirmDialog({
				title: 'Delete ' + type + '?',
				message: 'Really delete this ' + type + '?  Deleting it will move it to the trash.'
			})
            .show( dojo.hitch(this, function( confirmed ) {
                 if( confirmed ) {
				    	if (type == 'experiment') {
             			// Update database
             			dojo.xhrPut({ // FIXME: make webservice for this
	      					    url: "Experiments.pl",
	    					    putData: {
	    					    	fname: 'delete_experiment',
	    					    	eid: id
	    					    },
	    					    handleAs: "json",
	    					    load: dojo.hitch(this, function(data) {
	    					    	if (!data) {
	    					    		new InfoDialog({
    				                            title: 'Permission denied',
    				                            content: "You don't have permission to do that."
    				                        }).show();
    					    		return;
    					    	}
	    					    	// Remove node from tracklist
	    					    	dojo.query( '.coge-tracklist-label', this.div )
	    					        	.forEach( function( labelNode ) {
	    					        		if (labelNode.id == node.id) {
	    					        			labelNode.parentNode.parentNode.removeChild(labelNode.parentNode);
	    					        		}
	    					        	});
	    					    	// Remove track in browser
	    					    	this.browser.publish( '/jbrowse/v1/v/tracks/hide', [trackConfig] );
	    					   })
						    });
				    	}
				    	else if (type == 'notebook') {
				    		// Update database
				    		dojo.xhrPut({ // FIXME: make webservice for this
	      					    url: "NotebookView.pl",
	    					    putData: {
	    					    	fname: 'delete_list',
	    					    	lid: id
	    					    },
	    					    handleAs: "json",
	    					    load: dojo.hitch(this, function(data) {
	    					    	if (!data) {
	    					    		new InfoDialog({
    				                            title: 'Permission denied',
    				                            content: "You don't have permission to do that."
    				                        }).show();
    					    		return;
    					    	}
	    					    	// Remove all notebook nodes from tracklist
                     			div.parentNode.removeChild(div);
	    					    	// Remove track in browser
	    					    	this.browser.publish( '/jbrowse/v1/v/tracks/hide', [trackConfig] );
	    					   })
						    });
				    	}
                 }
             }));
    },

    //----------------------------------------------------------------

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

    //----------------------------------------------------------------

    _getVisibleConfigs: function (div, trackConfigs) {
    	var visibleConfigs = [];
    	dojo.query( '.coge-experiment, .coge-notebook', div ) // TODO: optimize this
    	    .forEach( function( labelNode ) {
    	    	if (!dojo.hasClass(labelNode.parentNode, 'collapsed')) {
    	        	for (var i = 0;  i < trackConfigs.length;  i++) {
    	        		var conf = trackConfigs[i];
    	        		var trackId = conf.coge.type + conf.coge.id;
    	        		if (labelNode.id == trackId) {
    	        			visibleConfigs.push(conf);
    	    				return;
    	        		}
    	        	}
    	    	}
    	    });
    	return visibleConfigs;
    },

    //----------------------------------------------------------------

    _hasLabelNode: function (div, nodes) {
    	var container;

    	nodes.forEach( function(node) {
    		dojo.query( '#'+node.id+' > .coge-tracklist-label', div )
    		 	.forEach( function( labelNode ) {
    		 		if (labelNode.id == node.id) {
    		 			container = labelNode.parentNode;
    		 			return;
    		 		}
    		 	});
    	});

    	return container;
    },

    //----------------------------------------------------------------

    moveTracks: function( source, nodes, copy, target ) {
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

            // return if no confs; whatever was dragged here probably wasn't a track from browser
            if( confs.length ) {
                // un-highlight track to show it is disabled
                this.dndDrop = true;
                this.browser.publish( '/jbrowse/v1/v/tracks/hide', confs ); // mdb: why not just call setTrackInactive directly?
                this.dndDrop = false;
            }
            else { // dragged-in from track selector
            	var notebookName;
            	var notebookId;
            	target.getAllNodes().forEach( function(node) { // FIXME: kludge
            		dojo.query('.coge-notebook', node).forEach( function(n) {
            			notebookName = n.id;
            			notebookId = notebookName.match(/\d+/);
            		});
            	});

            	if (notebookId && notebookId > 0) { // target is notebook
            		// Add experiment to notebook
            		nodes.forEach( dojo.hitch(this, function(node) {
            			var nodeId;
            			dojo.query('.coge-experiment', node).forEach( function(n) {
                			nodeId = n.id;
                		});
            			var conf = this.browser.trackConfigsByName[nodeId];
            			if (conf) {
            				dojo.xhrPut({
	      					    url: "NotebookView.pl",
	    					    putData: {
	    					    	fname: 'add_item_to_list',
	    					    	lid: notebookId,
	    					    	item_spec: '3:'+conf.coge.id, // FIXME: hardcoded type to experiment
	    					    },
	    					    handleAs: "json",
	    					    load: dojo.hitch(this, function(data) {
	    					    	if (!data) {
	    					    		new InfoDialog({
	    				                            title: 'Permission denied',
	    				                            content: "You don't have permission to do that."
	    				                        }).show();
	    					    		return;
	    					    	}
	    					    	// Reload track in browser
	    		                	this.browser.view.tracks.forEach( function(track) {
	    		                		if (track.config.label == notebookName) {
	    		                			track.changed();
	    		                		}
	    		                	});
	    					    })
						    });
            			}
            		}));
            	}
            }
		}
    },

    //----------------------------------------------------------------

    replaceTracks: function( trackConfigs ) { // mdb: unused now
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

    //----------------------------------------------------------------

    _setTrackTitle: function(config, node) {
//    	var coge = config.coge;
    	var name = config.key;
//    	var html;
//    	if (coge.type == 'notebook')
//    		html = '<img src="picts/notebook-icon-small.png"/>' + ' ';
//    	else if (coge.type == 'experiment')
//    		html = '<img src="picts/testtube-icon-small.png"/>' + ' ';
//    	html += '<img height="19" width="0" style="visibility:hidden;"/>'; // force min height
//    	html += '<span>' + name + '</span>';
    	node.innerHTML = name;
    },

    //----------------------------------------------------------------

    _getFeatureColor: function (id) { //FIXME: dup'ed in MultiXYPlot.js
    	return '#' + ((((id * 1234321) % 0x1000000) | 0x444444) & 0xe7e7e7 ).toString(16);
    },

    //----------------------------------------------------------------

    _info: function (trackConfig, type) {
        var dialog = new dijit.Dialog( { title: this._capitalize(type) + ' View' } );

        var iframe = dojo.create(
            'iframe', {
                tabindex: "0",
                width: $(window).width() * 0.8,
                height: $(window).height() * 0.8,
                style: { border: 'none' },
                src: trackConfig.coge.onClick
            });

        dialog.set( 'content', iframe );

        var updateIframeSize = function() {
            // hitch a ride on the dialog box's
            // layout function, which is called on
            // initial display, and when the window
            // is resized, to keep the iframe
            // sized to fit exactly in it.
            var cDims = domGeom.position( dialog.containerNode );
            var width  = cDims.w;
            var height = cDims.h - domGeom.position(dialog.titleBar).h;
            iframe.width = width;
            iframe.height = height;
        };
        aspect.after( dialog, 'layout', updateIframeSize );
        aspect.after( dialog, 'show', updateIframeSize );

        dialog.show();
    },

    //----------------------------------------------------------------

    _menu_popped_up: function () {
    	var m = $('.dijitMenuPopup');
    	if (m.length == 0)
    		return false;
    	if (m.css('display') == 'none')
    		return false;
    	return true;
    },

    //----------------------------------------------------------------

    /**
     * Given an array of track configs, update the track list to show
     * that they are turned on.
     */
    setTracksActive: function( /**Array[Object]*/ trackConfigs ) {
        var browser = this.browser;
    	console.log('setTracksActive ');
        dojo.query( '.coge-tracklist-label', this.div )
	        .forEach( function( labelNode, i ) {
	        	trackConfigs.forEach( function (trackConfig) {
	        		var trackId = trackConfig.coge.type + trackConfig.coge.id;
	        		if (labelNode.id == trackId) {
	    				dojo.addClass(labelNode, 'selected');
	        			if (dojo.hasClass(labelNode, 'coge-experiment')) {
	        				var id = trackConfig.coge.id;
                            var color;
                            var style = trackConfig.style;
                            var cookie = browser.cookie('track-' + trackConfig.track);

                            if (cookie) {
                                style = dojo.fromJson(cookie);
                            }

                            if (style.featureColor &&
                                style.featureColor[id]) {
	        				    color = style.featureColor[id];
                            } else {
	        				    color = this._getFeatureColor(id);
                            }
	        				dojo.style(labelNode, 'background', color);
	        			}
	        			else {
	        				dojo.style(labelNode, 'background', 'lightgray');
	        			}
	        		}
	        	});
	        });
    },

    //----------------------------------------------------------------

    _textFilter: function( text ) {
    	// Filter tracks
        if( text && /\S/.test(text) ) { // filter on text
            text = text.toLowerCase();
        	var n = dojo.byId('tracksAvail').firstChild;
        	var c = n.firstChild; // sequence
        	dojo.addClass(c, 'collapsed');
        	c = c.nextSibling; // gc
        	dojo.addClass(c, 'collapsed');
            this._traverseTracks(function(container) {
             	if (container.lastChild.title.toLowerCase().indexOf(text) != -1)
            		dojo.removeClass(container, 'collapsed');
            	else
            		dojo.addClass(container, 'collapsed');
            });
        }
        else { // empty string, show all
        	var n = dojo.byId('tracksAvail').firstChild;
        	var c = n.firstChild; // sequence
        	dojo.removeClass(c, 'collapsed');
        	c = c.nextSibling; // gc
        	dojo.removeClass(c, 'collapsed');
        	this._traverseTracks(function(container) {
        		if (dojo.hasClass(container.lastChild, 'coge-tracklist-collapsible'))
                	dojo.removeClass(container, 'collapsed');
            	else
                    dojo.addClass(container, 'collapsed');
            });
        }

        // Update filter label
        var count = dojo.query( '.coge-tracklist-container:not(.collapsed)', this.div ).length;
        this.textFilterLabel.innerHTML = count + ' track' + (count == 1 ? '' : 's') + ' shown';
    },
    
    //----------------------------------------------------------------

    _traverseTracks: function(f) {
    	var n = dojo.byId('tracksAvail').firstChild.nextSibling; // skip sequence & gc
    	while (n) {
    		var c = n.firstChild; // container
    		while (c) {
    			f(c);
    			c = c.nextSibling;
    		}
    		n = n.nextSibling;
    	}
    },

    //----------------------------------------------------------------

    /**
     * Update the display of the text filter control based on whether
     * it has any text in it.
     * @private
     */
    _updateTextFilterControl: function() {
        if( this.textFilterInput.value.length )
            dojo.setStyle( this.textFilterCancel, 'display', 'block' );
        else
            dojo.setStyle( this.textFilterCancel, 'display', 'none' );
    },

    //----------------------------------------------------------------

    /**
     * Given an array of track configs, update the track list to show
     * that they are turned off.
     */
    setTracksInactive: function( /**Array[Object]*/ trackConfigs ) {
    	console.log('setTracksInactive');
        dojo.query( '.coge-tracklist-label', this.div )
	        .forEach( function( labelNode, i ) {
	        	trackConfigs.forEach(function (trackConfig) {
	        		var trackId = trackConfig.coge.type + trackConfig.coge.id;
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

    //----------------------------------------------------------------

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

    //----------------------------------------------------------------

    /**
     * Make the track selector visible.
     * This does nothing for the Simple track selector, since it is always visible.
     */
    show: function() {
    },

    //----------------------------------------------------------------

    /**
     * Make the track selector invisible.
     * This does nothing for the Simple track selector, since it is always visible.
     */
    hide: function() {
    },

    //----------------------------------------------------------------

    /**
     * Toggle visibility of this track selector.
     * This does nothing for the Simple track selector, since it is always visible.
     */
    toggle: function() {
    }

});
});

//----------------------------------------------------------------

function create_notebook() {
	var name = $('#notebook_name').val();
	var description = $('#notebook_description').val();
	var restricted = $('#notebook_restricted').is(':checked');
	var coge_api = api_base_url.substring(0, api_base_url.length - 8);
	$.ajax(coge_api + '/notebooks?username='+un, {
		contentType: 'application/json',
		data: JSON.stringify({
			metadata: {
				name: name,
				description: description,
				restricted: restricted,
				type: 'experiment'
			}
		}),
		dataType: 'json',
		success: function(data) {
			if (data.error)
				alert(JSON.stringify(data.error));
			else {
				coge_track_list.trackListWidgets.push(coge_track_list._createDnDSource().insertNodes(
				    false,
				    [new_notebook(data.id, name, description, restricted)]));
				coge_track_list._textFilter();
				$("#tracksAvail").animate({ scrollTop: $('#tracksAvail').prop("scrollHeight")}, 1000);
				create_notebook_dialog.hide();
			}
		},
		type: 'PUT'
	});
}

//----------------------------------------------------------------

function new_notebook(id, name, description, restricted) {
	return {
        key: (restricted ? '&reg; ' : '') + name,
        baseUrl: api_base_url + '/experiment/notebook/' + id + '/',
        autocomplete: 'all',
        track: 'notebook' + id,
        label: 'notebook' + id,
        type: 'CoGe/View/Track/Wiggle/MultiXYPlot',
        storeClass: 'JBrowse/Store/SeqFeature/REST',
        style: { featureScale: 0.001 },
        showHoverScores: 1,
        coge: {
            id: id,
            type: 'notebook',
            classes: [
                'coge-tracklist-collapsible',
                'coge-tracklist-deletable',
                'coge-tracklist-info'
            ],
            collapsed: false,
            name: name,
            description: description,
            editable: true,
            experiments: null,
            count: 0,
            onClick: 'NotebookView.pl?embed=1&lid=' + id,
            menuOptions: [
                {
                    label: 'NotebookView',
                    action: "function() { window.open( 'NotebookView.pl?lid=" + id + "' ); }"
                }
            ]
        }
    };
}
