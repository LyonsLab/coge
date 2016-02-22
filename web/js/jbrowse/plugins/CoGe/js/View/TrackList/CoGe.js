var coge_track_list;

//----------------------------------------------------------------

define(['dojo/_base/declare',
        'dojo/_base/array',
        'dojo/query',
        'dojo/dom-attr',
        'dojo/dom-construct',
        'dojo/dom-geometry',
        'dojo/dom-style',
        'dojo/aspect',
        'dijit/layout/ContentPane',
        'dojo/dnd/Source',
        'dojo/mouse',
        'dijit/form/Button',
        'dijit/form/DropDownButton',
        'dijit/Menu',
        'dijit/MenuItem',
        'dijit/MenuSeparator',
        'dijit/Dialog',
        'JBrowse/View/ConfirmDialog',
        'JBrowse/View/InfoDialog'
       ],
       function( declare, array, query, attr, dom, domGeom, style, aspect, ContentPane, dndSource, mouse, Button, DropDownButton, Menu, MenuItem, MenuSeparator, Dialog, ConfirmDialog, InfoDialog ) {
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
        this._track_configs = args.trackConfigs;

        // make the track list DOM nodes and widgets
        this._create_track_list( args.browser.container );

        // maximum tracks that can be added via "+" button
        this.maxTracksToAdd = 20;
        
        this._create_search_button();

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
    			className: (collapsed ? 'coge-tracklist-expand' : 'coge-tracklist-collapse'),
    			src: (collapsed ? 'js/jbrowse/plugins/CoGe/img/arrow-right-icon.png' : 'js/jbrowse/plugins/CoGe/img/arrow-down-icon.png')
        	},
        	container
        );
    	dojo.connect(button, 'click', dojo.hitch(this, function() {
    		dojo.toggleClass(button, 'coge-tracklist-expand coge-tracklist-collapse');
            if (dojo.hasClass(button, 'coge-tracklist-expand'))
            	this._collapse(container);
            else
            	this._expand(container);
        	this._update_tracks_shown();
        }));
    },

    //----------------------------------------------------------------

    _add_to_notebook: function(items, notebook_id, create) {
	  	var coge_api = api_base_url.substring(0, api_base_url.length - 8);
	  	dojo.xhrPost({
	  		url: coge_api + '/notebooks/' + notebook_id + '/items/add?username='+un,
	  		postData: JSON.stringify({
	  			items: items
	  		}),
	  		handleAs: 'json',
	  		load: dojo.hitch(this, function(data) {
	  			if (data.error)
	  				this._error('Add Notebook Item', data);
	  			else {
	  				var n = dojo.byId('notebook' + notebook_id);
	  				if (create)
	  					items.forEach(function(item) {
	  						var e = dojo.byId(item.type + item.id);
	  						dojo.place(coge_track_list._new_track(e.config), n.parentNode);
	  					});
	  				this._expand(n);
	  			}
	  		}),
	  		error: dojo.hitch(this, function(data) {
	  			this._error('Add Notebook Item', data);
	  		})
	  	});
	  	if (this._add_dialog) {
	  		this._add_dialog.hide();
	  		this.add_dialog = null;
	  	}
    },

    //----------------------------------------------------------------

    _add_to_notebook_dialog: function(type, id, name) {
    	var content = '<table><tr><td><label>Notebook:</label></td><td><select id="coge_notebook">';
    	var notebooks = this._track_configs.filter(function(e) {
    		return (e.coge.type && e.coge.type == 'notebook' && e.coge.id != 0);
    	});
    	notebooks.sort();
    	for (var i=0; i<notebooks.length; i++)
    		content += '<option value="' + notebooks[i].coge.id + '">' + notebooks[i].coge.name + '</option>';
    	content += '</select></td></tr></table><div class="dijitDialogPaneActionBar"><button data-dojo-type="dijit/form/Button" type="button" onClick="var s=dojo.byId(\'coge_notebook\');coge_track_list._add_to_notebook([{type:\'' + type + '\',id:' + id + '}],s.options[s.selectedIndex].value,true);">OK</button><button data-dojo-type="dijit/form/Button" type="button" onClick="coge_track_list._add_dialog.hide()">Cancel</button></div>';
    	this._add_dialog = new Dialog({
            title: 'Add to Notebook ' + name,
            onHide: function(){this.destroyRecursive()},
            content: content
        });
    	this._add_dialog.show();
    },

    //----------------------------------------------------------------

    addTracks: function( track_configs ) {
//        // note that new tracks are, by default, hidden, so we just put them in the list
//        this.trackListWidget.insertNodes(
//            false,
//            track_configs
//        );
    },

    //----------------------------------------------------------------

    _build_title: function(track_config) {
    	var coge = track_config.coge;
    	var title = this._capitalize(coge.type) + " id" + coge.id;
    	if (coge.name)
    		title += "\nName: " + coge.name;
    	if (coge.description)
    		title += "\nDescription: " + coge.description;
    	if (coge.annotations)
    		title += coge.annotations;
    	return title;
    },

    //----------------------------------------------------------------

    _capitalize: function (string) {
        return string.charAt(0).toUpperCase() + string.substring(1);
    },

    //----------------------------------------------------------------

    _check_all: function(element, value) {
    	var cb = element.getElementsByTagName('INPUT');
    	for (var i=0; i<cb.length; i++)
    		cb[i].checked = value;
    },

    //----------------------------------------------------------------

    _collapse: function(container) {
    	container.firstChild.src = 'js/jbrowse/plugins/CoGe/img/arrow-right-icon.png';
    	var n = container.nextSibling;
    	while (n) {
    		dojo.addClass(n, 'collapsed');
    		n = n.nextSibling;
    	}
    },

    //----------------------------------------------------------------

    _create_notebook: function() {
	  	var name = dojo.getAttr('notebook_name', 'value');
	  	var description = dojo.getAttr('notebook_description', 'value');
	  	var restricted = dojo.getAttr('notebook_restricted', 'checked');
	  	var coge_api = api_base_url.substring(0, api_base_url.length - 8);
	  	dojo.xhrPut({
	  		url: coge_api + '/notebooks?username='+un,
	  		postData: JSON.stringify({
	  			metadata: {
	  				name: name,
	  				description: description,
	  				restricted: restricted,
	  				type: 'experiment'
	  			}
	  		}),
	  		handleAs: 'json',
	  		load: dojo.hitch(this, function(data) {
	  			if (data.error)
	  				this._error('Create Notebook', data);
	  			else {
	  				var config = this._new_notebook_config(data.id, name, description, restricted);
	  				this._track_configs.push(config);
	  	        	this._tracks.push(this._new_notebook().insertNodes(false, [config]));
	  				this._filter_tracks();
	  				this.div.scrollTop = this.div.scrollHeight;
	  				this._create_notebook_dialog.hide();
	  			}
	  		}),
	  		error: dojo.hitch(this, function(data) {
	  			this._error('Create Notebook', data);
	  		})
	  	});
    },

    //----------------------------------------------------------------

    _create_search_button: function() {
    	var content = '<div id="coge-search-dialog"><table><tr><td>Name:</td><td><input id="coge_search_text"></td></tr><tr><td>Chromosome:</td><td><select id="coge_ref_seq"><option>Any</option>';
    	var browser = this.browser;
    	browser.refSeqOrder.forEach(function(rs){
    		content += '<option>' + rs + '</option>';
    	})
    	content += '</select></td></tr><tr><td style="vertical-align:top;">Features:</td><td id="coge_search_features">';
    	var features = this._track_configs.filter(function(f) {
    		return (f.coge.type && f.coge.type == 'features');
    	});
    	features.forEach(function(t) {
    		content += '<div><input type="checkbox"> <label>' + t.coge.id + '</label></div>';
    	});
    	content += '</td></tr><tr><td></td><td><button onClick="coge_track_list._check_all(this.parentNode.parentNode.parentNode,true)">check all</button> <button onClick="coge_track_list._check_all(this.parentNode.parentNode.parentNode,false)">uncheck all</button></td></tr></table>';
    	content += '<div class="dijitDialogPaneActionBar"><button data-dojo-type="dijit/form/Button" type="button" onClick="coge_track_list._search_features()">OK</button><button data-dojo-type="dijit/form/Button" type="button" onClick="search_dialog.hide()">Cancel</button></div></div>';
        new Button({
	        label: 'Search',
	        onClick: function(event) {
	        	search_dialog = new Dialog({
                    title: "Search",
                    content: content,
                    onHide: function(){this.destroyRecursive()},
                    style: "width: 300px"
                });
//	        	for (var i=0; i<browser.refSeqOrder.length; i++)
//	        		if (browser.refSeqOrder[i] == browser.refSeq.name) {
//	        			dojo.byId('coge_ref_seq').selectedIndex = i + 1;
//	        			break;
//	        		}
	        	search_dialog.show();
	        	dojo.stopEvent(event);
	        },
	    }, dojo.create('button', null, this.browser.navbox));
    },

    //----------------------------------------------------------------

    _create_text_filter: function() {
        var div = dom.create( 'div', { className: 'coge-textfilter' }, this.pane);

        var d = dom.create('div',{ style: 'display:inline;overflow:show;position:relative' }, div );
		this.text_filter_input = dom.create('input', {
			className: 'coge-filter-input',
			type: 'text',
			placeholder: 'filter by text',
			onkeypress: dojo.hitch(this, function( evt ) {
				if (this._text_filter_timeout)
					window.clearTimeout( this._text_filter_timeout );
				this._text_filter_timeout = window.setTimeout(
					dojo.hitch( this, function() {
					      this._update_text_filter_control();
					      this._filter_tracks(this.text_filter_input.value);
					  }),
					500
				);
				this._update_text_filter_control();
				evt.stopPropagation();
			})
		}, d);
		this.text_filter_cancel = dom.create('div', {
			className: 'jbrowseIconCancel',
			id: 'text_filter_cancel',
			onclick: dojo.hitch( this, function() {
		        this.text_filter_input.value = '';
		        this._update_text_filter_control();
				this._filter_tracks();
			})
		}, d );

		var menu_button = dom.create('div', {id: 'coge_menu_button'}, div);
        var menu = new Menu();
        menu.addChild(new MenuItem({
            label: "Add All Tracks Shown",
            onClick: dojo.hitch(this, function() {
            	var visible_configs = [];
            	dojo.query('.coge-tracklist-container', this.div).forEach(function(container) {
           	    	if (!dojo.hasClass(container, 'collapsed'))
           	    		visible_configs.push(container.config);
            	});
				if (visible_configs.length) {
					if (visible_configs.length > this.maxTracksToAdd) {
					    var myDialog = new dijit.Dialog({
					        title: "Warning",
					        content: "There are too many tracks to add (>" + this.maxTracksToAdd + "), please filter them further.",
					        style: "width: 300px"
					    });
					    myDialog.show();
					}
					else {
						this.browser.publish('/jbrowse/v1/v/tracks/show', visible_configs);
					}
				}
			})
        }));
        menu.addChild(new MenuItem({
            label: "Clear All Tracks",
            onClick: dojo.hitch(this, function() {
            	var visible_configs = [];
	        	for (var i=0;  i<this._track_configs.length; i++)
	        		if (this._track_configs[i].label)
	        			visible_configs.push(this._track_configs[i]);
	        	if (visible_configs.length)
        			this.browser.publish('/jbrowse/v1/v/tracks/hide', visible_configs);
			})
        }));
        menu.addChild(new MenuSeparator());
        menu.addChild(new MenuItem({
            label: "Create New Notebook",
            onClick: dojo.hitch(this, function() {
            	this._create_notebook_dialog = new Dialog({
                    title: 'Create New Notebook',
                    content: '<table><tr><td><label>Name:</label></td><td><input id="notebook_name"></td></tr><tr><td><label>Description:</label></td><td><input id="notebook_description"></td></tr><tr><td><label>Restricted:</label></td><td><input type="checkbox" checked="checked" id="notebook_restricted"></td></tr></table><div class="dijitDialogPaneActionBar"><button data-dojo-type="dijit/form/Button" type="button" onClick="coge_track_list._create_notebook()">OK</button><button data-dojo-type="dijit/form/Button" type="button" onClick="coge_track_list._create_notebook_dialog.hide()">Cancel</button></div>',
                    onHide: function(){this.destroyRecursive()},
                    style: 'width: 300px'
                });
            	this._create_notebook_dialog.show();
            })
        }));
        menu.addChild(new MenuSeparator());
        menu.addChild(new MenuItem({
            label: 'Move Track Selector to ' + (dijit.byId('track_pane').get('region') == 'right' ? 'Left' : 'Right') + ' Side',
            onClick: function() {
            	var bc = dijit.byId('jbrowse');
            	var pane = dijit.byId('track_pane');
            	bc.removeChild(pane);
            	var region = pane.get('region');
            	this.set('label', 'Move Track Selector to ' + (region == 'right' ? 'Right' : 'Left') + ' Side');
            	localStorage.setItem("track-selector-side", region == 'right' ? 'left' : 'right');
            	pane.set('region', region == 'right' ? 'left' : 'right');
            	bc.addChild(pane);
            	
            }
        }));
        var btn = new DropDownButton({ dropDown: menu }, menu_button);
        menu.startup();
        btn.startup();

		this._tracks_shown = dom.create('div', { className: 'coge-tracks-shown' }, div );
    },

    //----------------------------------------------------------------

    _create_track_list: function(parent) {
        this.pane = dojo.create('div', { id: 'track_pane' }, parent);
        var side = 'right';
        if (localStorage) {
        	var s = localStorage.getItem('track-selector-side');
        	if (s)
        		side = s;
        }
        new ContentPane({region: side, splitter: true}, this.pane);
        this._create_text_filter(this.pane);
        this._update_text_filter_control();
        this.div = dojo.create('div', { id: 'coge-tracks' }, this.pane);

        this._tracks = [];
        this._track_configs.filter(function(tc) {
    		var type = tc.coge.type;
    		return (!type || type == 'sequence' || type == 'gc_content');
    	}).forEach(function(t){
//    		this._tracks.push(t);
    		this.div.appendChild(this._new_track(t));
    	}, this);

        var feature_groups = this._track_configs.filter(function(fg) {
    		return (fg.coge.type && fg.coge.type == 'feature_group');
    	});
    	var features = this._track_configs.filter(function(f) {
    		return (f.coge.type && f.coge.type == 'features');
    	});
    	feature_groups.forEach(function(fg) {
        	[fg].concat(features.filter(function(f) {
        		return f.coge.dataset_id == fg.coge.id;
        	})).forEach(function(t){
//        		this._tracks.push(t);
        		this.div.appendChild(this._new_track(t));
        	}, this);
    	}, this);

        // create a DnD source for each notebook
        var notebooks = this._track_configs.filter(function(e) {
    		return (e.coge.type && e.coge.type == 'notebook');
    	});
    	var experiments = this._track_configs.filter(function(e) {
    		return (e.coge.type && e.coge.type == 'experiment');
    	});
        notebooks.forEach(function(n) {
        	this._tracks.push(this._new_notebook().insertNodes(
            	false,
            	[n].concat(experiments.filter(function(e) {
            		return e.coge.notebooks && dojo.indexOf(e.coge.notebooks, n.coge.id) != -1;
            	}))
            ));
    	}, this);

        // show all tracks
        this._filter_tracks();
    },

    //----------------------------------------------------------------

    _delete: function(track_config, type, id, container) {
    	var Type = this._capitalize(type);
    	var message = 'Delete this ' + Type + '?  Deleting it will move it to the trash.';
    	if (type == 'notebook')
    		message += '<br>Note: Experiments in this notebook will NOT be deleted.'
		new ConfirmDialog({
				title: 'Delete ' + Type,
				message: message,
                onHide: function(){this.destroyRecursive()}
			})
            .show( dojo.hitch(this, function( confirmed ) {
                 if ( confirmed ) {
			    	if (type == 'experiment') {
             			dojo.xhrPut({ // FIXME: make webservice for this
      					    url: "Experiments.pl",
    					    putData: {
    					    	fname: 'delete_experiment',
    					    	eid: id
    					    },
    					    handleAs: "json",
    					    load: dojo.hitch(this, function(data) {
    					    	if (!data) {
    					    		this._error('Permission denied', "You don't have permission to do that.");
    					    		return;
    					    	}
    					    	dojo.destroy(container);
    					    	this.browser.publish( '/jbrowse/v1/v/tracks/hide', [track_config] );
    					   })
					    });
			    	} else if (type == 'notebook') {
			    		dojo.xhrPut({ // FIXME: make webservice for this
      					    url: "NotebookView.pl",
    					    putData: {
    					    	fname: 'delete_list',
    					    	lid: id
    					    },
    					    handleAs: "json",
    					    load: dojo.hitch(this, function(data) {
    					    	if (!data) {
    					    		this._error('Permission denied', "You don't have permission to do that.");
    					    		return;
    					    	}
    					    	dojo.destroy(container.parentNode);
    					    	this.browser.publish( '/jbrowse/v1/v/tracks/hide', [track_config] );
    					   })
					    });
			    	}
                 }
             }));
    },

    //----------------------------------------------------------------

    deleteTracks: function( /**Array[Object]*/ track_configs ) { // mdb: unused now ...?
//    	console.log('deleteTracks');
//        // remove any tracks in our track list that are being set as visible
//        array.forEach( track_configs || [], function( conf ) {
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

    //----------------------------------------------------------------

    _expand: function(container) {
	    container.firstChild.src = 'js/jbrowse/plugins/CoGe/img/arrow-down-icon.png';
		var n = container.nextSibling;
		while (n) {
			dojo.removeClass(n, 'collapsed');
			n = n.nextSibling;
		}
    },

    //----------------------------------------------------------------

    _filter_tracks: function( text ) {
        if (text && /\S/.test(text)) { // filter on text
            text = text.toLowerCase();
        	var n = this.div.firstChild;
        	while (n && dojo.hasClass(n, 'coge-tracklist-container')) {
        		dojo.addClass(n, 'collapsed');
        		n = n.nextSibling;
        	}
        	var matching_tracks = [];
            this._traverse_tracks(function(container) {
            	var t = container.lastChild.title ? container.lastChild.title : container.lastChild.innerHTML;
             	if (t.toLowerCase().indexOf(text) != -1) {
             		var i = 0;
             		for (; i<matching_tracks.length; i++)
             			if (container.id === matching_tracks[i].id)
             				break;
             		if (i === matching_tracks.length) {
             			dojo.removeClass(container, 'collapsed');
             			matching_tracks.push(container);
             		} else
                		dojo.addClass(container, 'collapsed');
             	} else
            		dojo.addClass(container, 'collapsed');
            });
        } else { // empty string, show all
        	var n = this.div.firstChild;
        	while (n && dojo.hasClass(n, 'coge-tracklist-container') && !dojo.hasClass(n.lastChild, 'coge-tracklist-collapsible')) {
        		dojo.removeClass(n, 'collapsed');
        		n = n.nextSibling;
        	}
        	while (n && dojo.hasClass(n, 'coge-tracklist-container')) {
        		if (dojo.hasClass(n.lastChild, 'coge-tracklist-collapsible'))
        			dojo.removeClass(n, 'collapsed');
        		n = n.nextSibling;
        	}
        	this._traverse_tracks(function(container) {
        		if (dojo.hasClass(container.lastChild, 'coge-tracklist-collapsible'))
                	dojo.removeClass(container, 'collapsed');
            	else
                    dojo.addClass(container, 'collapsed');
            });
        }
        this._update_tracks_shown();
    },

    //----------------------------------------------------------------

    _get_feature_color: function(id) { //FIXME: dup'ed in MultiXYPlot.js
    	return '#' + ((((id * 1234321) % 0x1000000) | 0x444444) & 0xe7e7e7 ).toString(16);
    },

    //----------------------------------------------------------------

    /**
     * Make the track selector invisible.
     * This does nothing for the Simple track selector, since it is always visible.
     */
    hide: function() {
    },

    //----------------------------------------------------------------

    _info: function (track_config, type) {
        var dialog = new dijit.Dialog( { title: this._capitalize(type) + ' View' } );

        var iframe = dojo.create(
            'iframe', {
                tabindex: "0",
                width: $(window).width() * 0.8,
                height: $(window).height() * 0.8,
                style: { border: 'none' },
                src: track_config.coge.onClick
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

    _in_notebook(container)
    {
    	return container.parentNode.firstChild.lastChild.innerHTML != 'All Experiments';    	
    },

    //----------------------------------------------------------------

    _menu_popped_up: function () {
    	var m = dojo.query('.dijitMenuPopup');
    	if (m.length == 0)
    		return false;
    	var display = m.style('display');
    	if (display instanceof Array) {
    		if (display.some(function(d){return d != 'none';}))
    			return true;
    		return false;
    	} else if (display == 'none')
    		return false;
    	return true;
    },

    //----------------------------------------------------------------

    _mouse_enter: function(track_config, node, container) {
    	if (this._menu_popped_up())
    		return;
    	var b = dijit.byId('coge_track_menu_button');
    	if (b)
    		b.destroy();
 		this._menu_node = node;
		this._menu_track_config = track_config;
		var id = track_config.coge.id;
		var type = track_config.coge.type;
		var Type = this._capitalize(type);
		var menu_button = dom.create('div', {id: 'coge_track_menu_button'}, container);
        var menu = new Menu();
        if (dojo.hasClass(node, 'coge-tracklist-info'))
            menu.addChild(new MenuItem({
                label: 'Info',
                onClick: dojo.hitch( this, function() {
                	this._info(track_config, type);
    			})
            }));
        if (dojo.hasClass(node, 'coge-tracklist-editable')) {
            menu.addChild(new MenuItem({
                label: 'Rename',
                onClick: dojo.hitch( this, function() {
                	this._rename_dialog(type, id, track_config.coge.name);
    			})
            }));
            if (type == 'experiment') {
                menu.addChild(new MenuItem({
                    label: 'Add to Notebook',
                    onClick: dojo.hitch( this, function() {
                    	this._add_to_notebook_dialog(type, id, track_config.coge.name);
        			})
                }));
            }
            if (track_config.coge.notebooks && this._in_notebook(container)) {
        		var notebook_node = container.parentNode.firstChild;                	
        		var notebook_name = notebook_node.lastChild.innerHTML;
                menu.addChild(new MenuItem({
                    label: 'Remove from Notebook',
                    onClick: dojo.hitch( this, function() {
                    	if (notebook_name.substring(0, 2) == '® ')
                    		notebook_name = notebook_name.substring(2);
                    	this._remove_from_notebook_dialog(type, id, track_config.coge.name, notebook_node.id.substring(8), notebook_name);
        			})
                }));
        	}
            menu.addChild(new MenuItem({
                label: 'Delete',
                onClick: dojo.hitch( this, function() {
                	this._delete(track_config, type, id, container);
    			})
            }));
        }
        var btn = new DropDownButton({ dropDown: menu }, menu_button);
        menu.startup();
        btn.startup();
    },

    //----------------------------------------------------------------

    _mouse_leave: function() {
		if (!this._menu_popped_up()) {
			var b = dijit.byId('coge_track_menu_button');
			if (b)
				b.destroy();
		}
	},

    //----------------------------------------------------------------

    moveTracks: function( source, nodes, copy, target ) {
    	if (source == target) // dropping in same place dragged from. should only happen in jbrowse track container
    		return;

    	var source_is_in_selector = source.node.firstChild.config;
    	var target_is_in_selector = target.node.firstChild.config;

    	if (source_is_in_selector && !target_is_in_selector) { // dragging a track from the selector onto jbrowse's track container
        	// get the configs from the tracks being dragged in
            var confs = dojo.filter(
                dojo.map( nodes, function(n) { return source.map[n.id].data; }),
                function(c) {return c;}
            );
            this.browser.publish('/jbrowse/v1/v/tracks/show', confs);
            return;
        }

        if (target_is_in_selector) {
            // get the configs from the tracks being dragged in
            var confs = dojo.filter(
                dojo.map( nodes, function(n) { return n.track && n.track.config; }),
                function(c) {return c;}
            );

            // return if no confs; whatever was dragged here probably wasn't a track from browser
            if (confs.length) {
                // un-highlight track to show it is disabled
                this.browser.publish( '/jbrowse/v1/v/tracks/hide', confs ); // mdb: why not just call setTrackInactive directly?
            } else { // dragged-in from track selector
            	var notebook_label;
            	target.getAllNodes().forEach( function(node) { // FIXME: kludge
            		dojo.query('.coge-notebook', node).forEach( function(n) {
            			notebook_label = n;
            		});
            	});

            	if (notebook_label) {
            		var items = [];
            		nodes.forEach(function(node){items.push({type: node.config.coge.type, id: node.config.coge.id});});
            		this._add_to_notebook(items, target.node.firstChild.config.coge.id);
            	}
            }
		}
    },

    //----------------------------------------------------------------

    _new_notebook: function() {
        var div = dojo.create( 'div', null, this.div );
    	return new dndSource(div, {
            accept: ["track"],
            checkAcceptance: function(source, nodes) {
            	for (var i=0; i<nodes.length; i++)
            		if (!nodes[i].config) // only accept experiments from the track selector (not from jbrowse's track container)
            			return false;
            	if (nodes[0].id.substring(0, 8) == 'notebook') {
            		dojo.publish("/dnd/cancel");
            		dojo.dnd.manager().stopDrag();
            		return false;
            	}
            	if (this.node.firstChild.id == 'notebook0') // "All Experiments"
            		return false;
            	if (this.node.firstChild.id == source.node.firstChild.id) // same notebook
            		return false;
            	for (var i=0; i<this.node.children.length; i++)
            		if (this.node.children[i].id == nodes[0].id) // already has experiment
            			return false;
            	return true;
            },
            copyOnly: true,
            creator: dojo.hitch(this, function(track_config, hint) {
                return {node: this._new_track(track_config), data: track_config, type: ["track", track_config.coge.type]};
            }),
            delay: 2,
            selfAccept: false
        });
    },

    //----------------------------------------------------------------

    _new_notebook_config: function(id, name, description, restricted) {
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
    			classes: [ 'coge-tracklist-collapsible', 'coge-tracklist-editable', 'coge-tracklist-info' ],
    			collapsed: false,
    			name: name,
    			description: description,
    			editable: true,
    			experiments: null,
    			count: 0,
    			onClick: 'NotebookView.pl?embed=1&lid=' + id,
    			menuOptions: [{
                    label: 'NotebookView',
                    action: "function() { window.open( 'NotebookView.pl?lid=" + id + "' ); }"
                }]
    		}
    	};
    },

    //----------------------------------------------------------------

    _new_track: function(track_config) {
    	var coge = track_config.coge;
        var container = dojo.create( 'div', {
        	className: 'coge-tracklist-container',
    		id: coge.type + coge.id,
        });
        container.config = track_config;
    	var label = dojo.create('div', {
			className: 'coge-tracklist-label coge-' + coge.type,
			title: this._build_title(track_config)
        });
        if (coge.classes)
        	dojo.addClass(label, coge.classes.join(' '));

        this._set_track_title(track_config, label);

        dojo.connect(label, "click", dojo.hitch(this, function(e) {
        	if (dojo.hasClass(label, 'selected'))
        		this.browser.publish( '/jbrowse/v1/v/tracks/hide', [track_config] );
        	else
        		this.browser.publish( '/jbrowse/v1/v/tracks/show', [track_config] );
        }));

        if (dojo.hasClass(label, 'coge-tracklist-editable') || dojo.hasClass(label, 'coge-tracklist-info')) {
        	dojo.connect(container, "onmouseenter", dojo.hitch(this, function(){this._mouse_enter(track_config, label, container)}));
        	dojo.connect(container, "onmouseleave", dojo.hitch(this, this._mouse_leave));
        }

        if (dojo.hasClass(label, 'coge-tracklist-collapsible'))
        	this._add_expander(container, coge.collapsed);
        else if (coge.collapsed)
        	dojo.addClass(container, 'collapsed');

        container.appendChild(label);
        return container;
    },

    //----------------------------------------------------------------

    _remove_from_notebook: function(type, id, notebook_id) {
	  	var coge_api = api_base_url.substring(0, api_base_url.length - 8);
	  	dojo.xhrPost({
	  		url: coge_api + '/notebooks/' + notebook_id + '/items/remove?username='+un,
	  		postData: JSON.stringify({
	  			items: [{
	  				type: type,
	  				id, id
	  			}]
	  		}),
	  		handleAs: 'json',
	  		load: dojo.hitch(this, function(data) {
	  			if (data.error)
	  				this._error('Remove Notebook Item', data);
	  			else
	  				dojo.destroy(this._menu_node.parentNode);
	  		}),
	  		error: dojo.hitch(this, function(data) {
	  			this._error('Remove Notebook Item', data);
	  		})
	  	});
    },

    //----------------------------------------------------------------

    _remove_from_notebook_dialog: function(type, id, name, notebook_id, notebook_name) {
		new ConfirmDialog({
			title: 'Remove ' + this._capitalize(type),
			message: 'Remove ' + type + ' "' + name + '" from notebook "' + notebook_name + '"?',
            onHide: function(){this.destroyRecursive()}
		}).show(dojo.hitch(this, function(confirmed) {
             if (confirmed)
            	 this._remove_from_notebook(type, id, notebook_id);
         }));
    },

    //----------------------------------------------------------------

    _rename: function (type, id, old_name) {
	  	var name = dojo.getAttr('name', 'value');
	  	if (name == old_name)
	  		return;
	  	var coge_api = api_base_url.substring(0, api_base_url.length - 8);
	  	dojo.xhrPost({
	  		url: coge_api + '/' + type + 's/' + id + '?username='+un,
	  		postData: JSON.stringify({
	  			metadata: {
	  				name: name
	  			}
	  		}),
	  		handleAs: 'json',
	  		load: dojo.hitch(this, function(data) {
	  			if (data.error)
	  				this._error('Rename ' + type, data);
	  			else {
	  				this._menu_track_config.coge.name = name;
	  				var key = this._menu_track_config.key;
	  				if (key.length > 6 && key.substring(0, 6) == '&reg; ') {
	  					key = '&reg; ' + name;
	  					old_name = '&reg; ' + old_name;
	  				} else
	  					key = name;
	  				this._menu_track_config.key = key;
	  				this._set_track_title(this._menu_track_config, this._menu_node);
	  				this._track_changed(old_name, key);
	  				this._rename_dialog.hide();
	  			}
	  		}),
	  		error: dojo.hitch(this, function(data) {
	  			this._error('Rename ' + type, data);
	  		})
	  	});
    },

    //----------------------------------------------------------------

    _rename_dialog: function(type, id, name) {
    	this._rename_dialog = new Dialog({
            title: 'Rename ' + type,
            content: '<table><tr><td><label>Name:</label></td><td><input id="name" value=' + JSON.stringify(name) + '></td></tr></table><div class="dijitDialogPaneActionBar"><button data-dojo-type="dijit/form/Button" type="button" onClick="coge_track_list._rename(\'' + type + '\',' + id + ',\'' + name.replace(/'/g, "\\'") + '\')">OK</button><button data-dojo-type="dijit/form/Button" type="button" onClick="coge_track_list._rename_dialog.hide()">Cancel</button></div>',
            onHide: function(){this.destroyRecursive()},
            style: "width: 300px"
        });
    	this._rename_dialog.show();
    	var i = dojo.byId("name");
    	i.setSelectionRange(0, i.value.length);
    },

    //----------------------------------------------------------------

    replaceTracks: function( track_configs ) { // mdb: unused now
    	console.log('replaceTracks');
//        // for each one
//        array.forEach( track_configs, function( conf ) {
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

	_search_features: function() {
    	var features = document.getElementById('coge_search_features').getElementsByTagName('INPUT');
    	var types = [];
    	for (var i=0; i<features.length; i++)
    		if (features[i].checked)
    			types.push("'" + features[i].nextElementSibling.innerText + "'");
    	if (!types.length) {
    		this._error('Search', 'Please select one or more feature types to search.');
    		return;
    	}
    	var url = api_base_url + '/genome/' + gid + '/features?name=' + encodeURIComponent(dojo.byId('coge_search_text').value) + '&features=' + (types.length == features.length ? 'all' : types.join());
    	var ref_seq = dojo.byId('coge_ref_seq');
    	if (ref_seq.selectedIndex > 0)
    		url += '&chr=' + ref_seq.options[ref_seq.selectedIndex].innerHTML;
    	dojo.xhrGet({
    		url: url,
    		handleAs: 'json',
	  		load: dojo.hitch(this, function(data) {
	  			if (data.error) {
	  				this._error('Search', data);
	  				return;
	  			}
	  			if (data.length == 0) {
	  				this._error('Search', 'no features found');
	  				return;
	  			}
	  			dojo.query('.dijitDialogUnderlay').style('opacity', 0);
  				var div = dojo.byId('coge-search-dialog');
  				div.style.maxHeight = '500px';
  				div.style.overflow = 'auto';
    			dojo.empty(div);
    			data.forEach(function(hit){
    				dojo.create('a', {
    					innerHTML: hit.name,
    					onclick: dojo.hitch(hit, function() {
	    					coge_track_list.browser.navigateToLocation(this.location);
	    					return false;
	    				})
    				}, div);
    				dojo.create('br', null, div);
    			});
    		}),
    		error: dojo.hitch(this, function(data) {
    			this._error('Search', data);
    		})
    	})
//    	search_dialog.hide();
    },

    //----------------------------------------------------------------

    /**
     * Given an array of track configs, update the track list to show
     * that they are turned on.
     */
    setTracksActive: function( /**Array[Object]*/ track_configs ) {
        var browser = this.browser;
        dojo.query( '.coge-tracklist-label', this.div )
	        .forEach( function( labelNode, i ) {
	        	track_configs.forEach( function (trackConfig) {
	        		var trackId = trackConfig.coge.type + trackConfig.coge.id;
	        		if (labelNode.parentNode.id == trackId) {
	    				dojo.addClass(labelNode, 'selected');
	        			if (dojo.hasClass(labelNode, 'coge-experiment')) {
	        				var id = trackConfig.coge.id;
                            var color;
                            var style = trackConfig.style;
                            var cookie = browser.cookie('track-' + trackConfig.track);
                            if (cookie)
                                style = dojo.fromJson(cookie);
                            if (style.featureColor && style.featureColor[id])
	        				    color = style.featureColor[id];
                            else
	        				    color = coge_track_list._get_feature_color(id);
	        				dojo.style(labelNode, 'background', color);
	        			} else
	        				dojo.style(labelNode, 'background', 'lightgray');
	        		}
	        	});
	        });
    },

    //----------------------------------------------------------------

    /**
     * Given an array of track configs, update the track list to show
     * that they are turned off.
     */
    setTracksInactive: function( /**Array[Object]*/ track_configs ) {
        dojo.query( '.coge-tracklist-label', this.div )
	        .forEach( function( labelNode, i ) {
	        	track_configs.forEach(function (trackConfig) {
	        		var trackId = trackConfig.coge.type + trackConfig.coge.id;
	        		if (labelNode.parentNode.id == trackId) {
	    				dojo.style(labelNode, 'background', '');
	        			dojo.removeClass(labelNode, 'selected');
	        		}
	        	});
	        });
    },

    //----------------------------------------------------------------

    _set_track_title: function(track_config, node) {
//    	var coge = track_config.coge;
    	var name = track_config.key;
//    	var html;
//    	if (coge.type == 'notebook')
//    		html = '<img src="picts/notebook-icon-small.png"/>' + ' ';
//    	else if (coge.type == 'experiment')
//    		html = '<img src="picts/testtube-icon-small.png"/>' + ' ';
//    	html += '<img height="19" width="0" style="visibility:hidden;"/>'; // force min height
//    	html += '<span>' + name + '</span>';
    	node.innerHTML = name;
    	node.title = this._build_title(track_config);
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
     * Toggle visibility of this track selector.
     * This does nothing for the Simple track selector, since it is always visible.
     */
    toggle: function() {
    },
    
    //----------------------------------------------------------------

    _track_changed: function(key, new_key) {
    	this.browser.view.tracks.forEach( function(track) {
    		if (track.config.key == key) {
    			if (new_key) {
    				track.config.key = new_key;
    				track.label.children[1].innerHTML = new_key;
    			}
    			track.changed();
    		}
    	});
    },

    //----------------------------------------------------------------

    _traverse_tracks: function(f) {
    	var n = this.div.firstChild;
    	while (n && dojo.hasClass(n, 'coge-tracklist-container'))
    		n = n.nextSibling;
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
    // show cancel button only when text filter has text

    _update_text_filter_control: function() {
        if (this.text_filter_input.value.length)
            dojo.setStyle(this.text_filter_cancel, 'display', 'block');
        else
            dojo.setStyle(this.text_filter_cancel, 'display', 'none');
    },

    //----------------------------------------------------------------

    _update_tracks_shown: function() {
        var count = dojo.query('.coge-tracklist-container:not(.collapsed)').length;
        this._tracks_shown.innerHTML = count + ' track' + (count == 1 ? '' : 's') + ' shown';
    }
});
});
