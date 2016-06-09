var coge_track_list;

// ----------------------------------------------------------------

define(['dojo/_base/declare',
		'dojo/dom-construct',
		'dojo/dom-geometry',
		'dojo/aspect',
		'dijit/layout/ContentPane',
		'dojo/dnd/Source',
		'dijit/form/DropDownButton',
		'dijit/Menu',
		'dijit/MenuItem',
		'dijit/MenuSeparator',
		'dojo/Deferred',
		'dijit/Dialog'
	   ],
	   function( declare, dom, domGeom, aspect, ContentPane, DndSource, DropDownButton, Menu, MenuItem, MenuSeparator, Deferred, Dialog ) {
	function natural_sort (a, b) {
	    var re = /(^-?[0-9]+(\.?[0-9]*)[df]?e?[0-9]?$|^0x[0-9a-f]+$|[0-9]+)/gi,
	        sre = /(^[ ]*|[ ]*$)/g,
	        dre = /(^([\w ]+,?[\w ]+)?[\w ]+,?[\w ]+\d+:\d+(:\d+)?[\w ]?|^\d{1,4}[\/\-]\d{1,4}[\/\-]\d{1,4}|^\w+, \w+ \d+, \d{4})/,
	        hre = /^0x[0-9a-f]+$/i,
	        ore = /^0/,
	        i = function(s) { return (''+s).toLowerCase() || ''+s },
	        // convert all to strings strip whitespace
	        x = i(a).replace(sre, '') || '',
	        y = i(b).replace(sre, '') || '',
	        // chunk/tokenize
	        xN = x.replace(re, '\0$1\0').replace(/\0$/,'').replace(/^\0/,'').split('\0'),
	        yN = y.replace(re, '\0$1\0').replace(/\0$/,'').replace(/^\0/,'').split('\0'),
	        // numeric, hex or date detection
	        xD = parseInt(x.match(hre)) || (xN.length != 1 && x.match(dre) && Date.parse(x)),
	        yD = parseInt(y.match(hre)) || xD && y.match(dre) && Date.parse(y) || null,
	        oFxNcL, oFyNcL;
	    // first try and sort Hex codes or Dates
	    if (yD)
	        if ( xD < yD ) return -1;
	        else if ( xD > yD ) return 1;
	    // natural sorting through split numeric strings and default strings
	    for(var cLoc=0, numS=Math.max(xN.length, yN.length); cLoc < numS; cLoc++) {
	        // find floats not starting with '0', string or 0 if not defined (Clint Priest)
	        oFxNcL = !(xN[cLoc] || '').match(ore) && parseFloat(xN[cLoc]) || xN[cLoc] || 0;
	        oFyNcL = !(yN[cLoc] || '').match(ore) && parseFloat(yN[cLoc]) || yN[cLoc] || 0;
	        // handle numeric vs string comparison - number < string - (Kyle Adams)
	        if (isNaN(oFxNcL) !== isNaN(oFyNcL)) { return (isNaN(oFxNcL)) ? 1 : -1; }
	        // rely on string comparison if different types - i.e. '02' < 2 != '02' < '2'
	        else if (typeof oFxNcL !== typeof oFyNcL) {
	            oFxNcL += '';
	            oFyNcL += '';
	        }
	        if (oFxNcL < oFyNcL) return -1;
	        if (oFxNcL > oFyNcL) return 1;
	    }
	    return 0;
	}

	// ----------------------------------------------------------------

    return declare( 'JBrowse.View.TrackList.CoGe', null,

	/** @lends JBrowse.View.TrackList.CoGe.prototype */
	{
	/**
	 * CoGe drag-and-drop track selector.
	 * 
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

	// ----------------------------------------------------------------

	_add_expander: function (container) {
		var button = dom.create('img', {
				className: 'coge-tracklist-expand',
				src: 'js/jbrowse/plugins/CoGe/img/arrow-right-icon.png'
			},
			container
		);
		dojo.connect(button, 'click', dojo.hitch(this, function() {
			if (button.src.indexOf('down') != -1)
				this._collapse(container);
			else
				this._expand(container);
			this._update_tracks_shown();
		}));
	},

	// ----------------------------------------------------------------
	// add experiments to notebook in CoGe

	_add_to_notebook: function(track_configs, notebook_id, create) {
		var on_error = function() {
			if (!create) {
				var notebook_node = dojo.byId('notebook' + notebook_id);
				track_configs.forEach(function(item) {
					dojo.destroy(dojo.query('#' + item.type + item.id, notebook_node)[0]);
				});
			}
		}
		var coge_api = api_base_url.substring(0, api_base_url.length - 8);
		var items = [];
		track_configs.forEach(function(track_config){
			items.push({type: track_config.coge.type, id: track_config.coge.id});
		});
		dojo.xhrPost({
			url: coge_api + '/notebooks/' + notebook_id + '/items/add?username='+un,
			postData: JSON.stringify({ items: items }),
			handleAs: 'json',
			load: dojo.hitch(this, function(data) {
				if (data.error) {
					on_error();
					coge_plugin.error('Add Notebook Item', data);
				} else {
					var n = dojo.byId('notebook' + notebook_id);
					if (create)
						track_configs.forEach(function(track_config) {
							var e = dojo.byId(track_config.coge.type + track_config.coge.id);
							var track = this._new_track(e.config);
							track.style.display = 'none';
							this._add_track_to_notebook(track, n);
						}, this);
					if (n.style.display != 'none')
						this._expand(n);
				}
			}),
			error: function(data) {
				on_error();
				coge_plugin.error('Add Notebook Item', data);
			}
		});
		if (this._add_dialog)
			this._add_dialog.hide();
	},

	// ----------------------------------------------------------------

	_add_to_notebook_dialog: function(track_config) {
		var content = '<table><tr><td><label>Notebook:</label></td><td><select id="coge_notebook">';
		var notebooks = this._track_configs.filter(function(e) {
			return (e.coge.type && e.coge.type == 'notebook' && e.coge.id != 0);
		});
		notebooks.sort(function(a, b) { return natural_sort(a.coge.name, b.coge.name); });
		for (var i=0; i<notebooks.length; i++)
			content += '<option value="' + notebooks[i].coge.id + '">' + notebooks[i].coge.name + '</option>';
		content += '</select></td></tr></table><div class="dijitDialogPaneActionBar"><button data-dojo-type="dijit/form/Button" type="button" onClick="coge_track_list._add_dialog.onExecute()">OK</button><button data-dojo-type="dijit/form/Button" type="button" onClick="coge_track_list._add_dialog.hide()">Cancel</button></div>';
		this._add_dialog = new Dialog({
			title: 'Add experiment ' + track_config.coge.name + ' to Notebook',
			onExecute: function() {
				var s=dojo.byId('coge_notebook');
				coge_track_list._add_to_notebook([track_config],s.options[s.selectedIndex].value,true);
				this.onHide();
			},
			onHide: function() {
				this.destroyRecursive();
				coge_track_list._add_dialog = null;
			},
			content: content
		});
		this._add_dialog.show();
	},

	// ----------------------------------------------------------------

	_add_track_to_notebook: function(track, notebook) {
		if (!notebook.config.coge.experiments)
			notebook.config.coge.experiments = [];
		notebook.config.coge.experiments.push({id:track.config.coge.id, name:track.config.coge.name, type:track.config.coge.data_type})
		if (!notebook.nextSibling) {
			dojo.place(track, notebook, 'after');
			return;
		}
		var n = notebook.nextSibling;
		if (natural_sort(n.config.coge.name, track.config.coge.name) > 0) {
			dojo.place(track, n, 'before');
			return;
		}
		while (n.nextSibling && natural_sort(track.config.coge.name, n.nextSibling.config.coge.name) > 0)
			n = n.nextSibling;
		dojo.place(track, n, 'after');
		if (dojo.byId('track_notebook' + notebook.config.coge.id))
			this._traverse_tracks(function(container){
				if (container.config.coge.type == 'experiment' && container.config.coge.id == track.config.coge.id)
					dojo.create('div', { className: 'coge-circle', style: { backgroundColor: coge_track_list._get_track_color(container) } }, container, 'first');
			})
		this._reload_notebook(notebook.config.coge.id);
	},

	// ----------------------------------------------------------------

	addTracks: function(track_configs) {
		track_configs.forEach(function(track_config) {
			if (track_config.coge)
				if (track_config.coge.search_track)
					this.tracks_div.insertBefore(this._new_track(track_config), this.tracks_div.firstChild); // insert before Sequence track at top
				else if (track_config.coge.type != 'notebook') {
					this._add_to_notebook([track_config], 0, true);
				}
		}, this);
	},

	// ----------------------------------------------------------------

	_build_title: function(track_config) {
		var coge = track_config.coge;
		var title = this._capitalize(coge.type) + ' id: ' + coge.id;
		if (coge.type == 'experiment')
			title += "\nData Type: " + (coge.data_type == 4 ? 'Markers' : coge.data_type == 3 ? 'Alignments' : coge.data_type == 2 ? 'Polymorphism Data' : 'Quantitative Data');
		if (coge.name)
			title += "\nName: " + coge.name;
		if (coge.description)
			title += "\nDescription: " + coge.description;
		if (coge.annotations)
			title += coge.annotations;
		return title;
	},

	// ----------------------------------------------------------------

	_capitalize: function (string) {
		return string.charAt(0).toUpperCase() + string.substring(1);
	},

	// ----------------------------------------------------------------

	_collapse: function(container) {
		container.config.coge.expanded = false;
		container.firstChild.src = 'js/jbrowse/plugins/CoGe/img/arrow-right-icon.png';
		var n = container.nextSibling;
		while (n) {
			n.style.display = 'none';
			n = n.nextSibling;
		}
	},

	// ----------------------------------------------------------------

	_create_notebook: function() {
		var self = this;
		var browser = this.browser;
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
			load: function(data) {
				if (data.error)
					coge_plugin.error('Create Notebook', data);
				else {
			        var d = new Deferred();
					var config = self._new_notebook_config(data.id, name, description, restricted);
					var store_config = {
						browser: browser,
						refSeq: browser.refSeq,
						type: 'JBrowse/Store/SeqFeature/REST',
						baseUrl: config.baseUrl
					};
					var store_name = browser.addStoreConfig(undefined, store_config);
					store_config.name = store_name;
					browser.getStore(store_name, function(store) {
			           d.resolve(true);
			       	});
			       	d.promise.then(function() {
						config.store = store_name;
						self._track_configs.push(config);
						self._new_notebook_source().insertNodes(false, [config]);
						self._filter_tracks();
						dojo.byId('notebook' + config.coge.id).parentNode.scrollIntoView();
						self._create_notebook_dialog.hide();
						self.browser.publish('/jbrowse/v1/v/tracks/new', [config]);
					});
				}
			},
			error: function(data) {
				coge_plugin.error('Create Notebook', data);
			}
		});
	},

	// ----------------------------------------------------------------

	_create_notebook_and_experiment_tracks: function(notebook_config, experiments) {
		var source = this._new_notebook_source();
		source.insertNodes(false, [notebook_config]);
		if (notebook_config.coge.id != 0)
			experiments = experiments.filter(function(e) {
				return e.coge.notebooks && dojo.indexOf(e.coge.notebooks, notebook_config.coge.id) != -1;
			});
		source.insertNodes(false, experiments.sort(function(a, b) { return natural_sort(a.coge.name, b.coge.name); }));
	},

	// ----------------------------------------------------------------

	_create_text_filter: function(parent) {
		var div = dom.create( 'div', { className: 'coge-textfilter' }, parent);

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
		this.text_filter_cancel = dom.create('span', {
			className: 'glyphicon glyphicon-remove',
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
				dojo.query('.coge-track', this.tracks_div).forEach(function(container) {
					if (container.style.display != 'none')
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
			label: "Show Quantitative Tracks",
			onClick: dojo.hitch(this, function() {
				this._type_filter = 1;
				this._filter_tracks(this.text_filter_input.value);
			})
		}));
		menu.addChild(new MenuItem({
			label: "Show Polymorphism Tracks",
			onClick: dojo.hitch(this, function() {
				this._type_filter = 2;
				this._filter_tracks(this.text_filter_input.value);
			})
		}));
		menu.addChild(new MenuItem({
			label: "Show Alignment Tracks",
			onClick: dojo.hitch(this, function() {
				this._type_filter = 3;
				this._filter_tracks(this.text_filter_input.value);
			})
		}));
		menu.addChild(new MenuItem({
			label: "Show Marker Tracks",
			onClick: dojo.hitch(this, function() {
				this._type_filter = 4;
				this._filter_tracks(this.text_filter_input.value);
			})
		}));
		if (un != 'public') {
			menu.addChild(new MenuSeparator());
			menu.addChild(new MenuItem({
				label: "Create New Notebook",
				onClick: dojo.hitch(this, function() {
					this._create_notebook_dialog = new Dialog({
						title: 'Create New Notebook',
						content: '<table><tr><td><label>Name:</label></td><td><input id="notebook_name"></td></tr><tr><td><label>Description:</label></td><td><input id="notebook_description"></td></tr><tr><td><label>Restricted:</label></td><td><input type="checkbox" checked="checked" id="notebook_restricted"></td></tr></table><div class="dijitDialogPaneActionBar"><button data-dojo-type="dijit/form/Button" type="button" onClick="coge_track_list._create_notebook()">OK</button><button data-dojo-type="dijit/form/Button" type="button" onClick="coge_track_list._create_notebook_dialog.hide()">Cancel</button></div>',
						onHide: function() {
							this.destroyRecursive();
							coge_track_list._create_notebook_dialog = null;
						},
						style: 'width: 300px'
					});
					this._create_notebook_dialog.show();
				})
			}));
		}
		menu.addChild(new MenuSeparator());
		menu.addChild(new MenuItem({
			label: 'Move Track Selector to ' + (dijit.byId('coge').get('region') == 'right' ? 'Left' : 'Right') + ' Side',
			onClick: function() {
				var parent = dijit.byId('jbrowse');
				var root = dijit.byId('coge');
				parent.removeChild(root);
				var region = root.get('region');
				this.set('label', 'Move Track Selector to ' + (region == 'right' ? 'Right' : 'Left') + ' Side');
				localStorage.setItem("track-selector-side", region == 'right' ? 'left' : 'right');
				root.set('region', region == 'right' ? 'left' : 'right');
				parent.addChild(root);
				
			}
		}));
		var btn = new DropDownButton({ dropDown: menu }, menu_button);
		menu.startup();
		btn.startup();

		this._tracks_shown = dom.create('div', { className: 'coge-tracks-shown' }, div );
	},

	// ----------------------------------------------------------------

	_create_track_list: function(parent) {
		var root = dojo.create('div', { id: 'coge', style: 'display:flex;overflow:hidden;' }, parent);
		var side = 'right';
		if (localStorage) {
			var s = localStorage.getItem('track-selector-side');
			if (s)
				side = s;
		}
		var cp = new ContentPane({region: side, splitter: true}, root);
		aspect.after(cp, "resize", function(size) { dojo.byId('track_pane').style.width=size.w - dojo.position('feature_hits').w + 'px'; }, true);

		dojo.create('div', { id: 'feature_hits' }, root);
	
		var pane = dojo.create('div', { id: 'track_pane' }, root);
		this._create_text_filter(pane);
		this._update_text_filter_control();
		var scroller = dojo.create('div', { style: 'height:100%;overflow-x:hidden;overflow-y:auto;' }, pane);
		this.tracks_div = dojo.create('div', null, scroller);

		this._track_configs.filter(function(tc) {
			var type = tc.coge.type;
			return (!type || type == 'sequence' || type == 'gc_content');
		}).forEach(function(t) {
			this.tracks_div.appendChild(this._new_track(t));
		}, this);

		var feature_groups = this._track_configs.filter(function(fg) {
			return (fg.coge.type && fg.coge.type == 'feature_group');
		});
		var features = this._track_configs.filter(function(f) {
			return (f.coge.type && f.coge.type == 'features');
		});
		feature_groups.forEach(function(fg) {
			var d = dojo.create('div', null, this.tracks_div);
			d.appendChild(this._new_track(fg));
			features.filter(function(f) {
				return f.coge.dataset_id == fg.coge.id;
			}).forEach(function(t) {
				d.appendChild(this._new_track(t, true));
			}, this);
		}, this);

		// create a DnD source for each notebook
		var notebooks = this._track_configs.filter(function(e) {
			return (e.coge.type && e.coge.type == 'notebook');
		});
		var experiments = this._track_configs.filter(function(e) {
			return (e.coge.type && e.coge.type == 'experiment');
		});
		this._create_notebook_and_experiment_tracks(notebooks.shift(), experiments);
		notebooks.sort(function(a, b) { return natural_sort(a.coge.name, b.coge.name); });
		notebooks.forEach(function(n) {
			this._create_notebook_and_experiment_tracks(n, experiments);
		}, this);

		// show all tracks
		this._filter_tracks();
	},

	// ----------------------------------------------------------------

	_delete: function(track_config, type, id, container) {
		var Type = this._capitalize(type);
		var message = 'Delete this ' + Type + '?  Deleting it will move it to the trash.';
		if (type == 'notebook')
			message += '<br>Note: Experiments in this notebook will NOT be deleted.'
		coge_plugin.confirm('Delete ' + Type, message, dojo.hitch(this, function() {
			var coge_api = api_base_url.substring(0, api_base_url.length - 8);
			dojo.xhrDelete({
				url: coge_api + '/' + type + 's/' + id,
				handleAs: "json",
				load: dojo.hitch(this, function(data) {
					if (data.error)
						coge_plugin.error('Delete ' + Type, data);
					else {
					   dojo.destroy(type == 'notebook' ? container.parentNode : container);
					   this.browser.publish('/jbrowse/v1/v/tracks/hide', [track_config]);
					}
				}),
				error: function(data) {
					coge_plugin.error('Delete ' + Type, data);
				}
			});
		}));
	},

	// ----------------------------------------------------------------
	// are we supposed to delete the experiment(s) from the database? for now just acknowledge their removal from the view

	deleteTracks: function(track_configs) {
		this.setTracksInactive(track_configs);
	},

	// ----------------------------------------------------------------

	_expand: function(container) {
		container.config.coge.expanded = true;
		container.firstChild.src = 'js/jbrowse/plugins/CoGe/img/arrow-down-icon.png';
		var n = container.nextSibling;
		while (n) {
			n.style.display = '';
			n = n.nextSibling;
		}
	},

	// ----------------------------------------------------------------

	_filter_tracks: function(text) {
		var type_filter = this._type_filter;
		if (text && /\S/.test(text)) { // filter on text
			text = text.toLowerCase();
			var already_shown = {};
			this._traverse_tracks(function(container) {
				var t = container.lastChild.title ? container.lastChild.title : container.lastChild.innerHTML;
				if (t.toLowerCase().indexOf(text) != -1) {
					if (!already_shown[container.id]) {
						if (!type_filter || type_filter === container.config.coge.data_type) {
							container.style.display = '';
							already_shown[container.id] = true;
						}
					} else
						container.style.display = 'none';
				} else
					container.style.display = 'none';
			});
		} else if (type_filter) {
			var already_shown = {};
			this._traverse_tracks(function(container) {
				if (!already_shown[container.id] && type_filter === container.config.coge.data_type) {
					container.style.display = '';
					already_shown[container.id] = true;
				} else
					container.style.display = 'none';
			});
		} else { // empty string, show all
			var expanded = true;
			this._traverse_tracks(function(container) {
				if (container.config.coge.collapsible) {
					expanded = container.config.coge.expanded;
					container.style.display = '';
				} else
					container.style.display = expanded ? '' : 'none';
			});
		}
		this._update_tracks_shown();
	},

	// ----------------------------------------------------------------

	get_search_config: function(eid) {
		var n = this.tracks_div.firstChild;
		while (n && dojo.hasClass(n, 'coge-track')) {
			if (n.config.coge.id == eid)
				return n.config;
			n = n.nextSibling;
		}
		return null;
	},

	// ----------------------------------------------------------------

	_get_track_color: function(container) {
		var id = container.config.coge.id;
		var style = container.config.style;
		var cookie = this.browser.cookie('track-' + container.config.track);
		if (cookie)
			style = dojo.fromJson(cookie);
		if (style.featureColor && style.featureColor[id])
			return style.featureColor[id];
		return coge_plugin.calc_color(id);
	},

	// ----------------------------------------------------------------
	/**
	 * Make the track selector invisible
	 */

	hide: function() {
	},

	// ----------------------------------------------------------------

	_info: function(track_config, type) {
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

	// ----------------------------------------------------------------

	_in_notebook(container)
	{
		return container.parentNode.firstChild.lastChild.innerHTML != 'All Experiments';    	
	},

	// ----------------------------------------------------------------

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

	// ----------------------------------------------------------------

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
		menu.addChild(new MenuItem({
			label: this._capitalize(type) + ' View',
			onClick: function() {
				coge_track_list._info(track_config, type);
			}
		}));
		if (track_config.coge.editable) {
			menu.addChild(new MenuItem({
				label: 'Rename',
				onClick: function() {
					coge_track_list.rename_dialog(type, id, track_config.coge.name);
				}
			}));
			if (type == 'experiment') {
				menu.addChild(new MenuItem({
					label: 'Add to Notebook',
					onClick: function() {
						coge_track_list._add_to_notebook_dialog(track_config);
					}
				}));
			}
			if (track_config.coge.notebooks && this._in_notebook(container)) {
				var notebook_node = container.parentNode.firstChild;                	
				var notebook_name = notebook_node.lastChild.innerHTML;
				menu.addChild(new MenuItem({
					label: 'Remove from Notebook',
					onClick: function() {
						if (notebook_name.substring(0, 2) == '® ')
							notebook_name = notebook_name.substring(2);
						coge_track_list._remove_from_notebook_dialog(type, id, track_config.coge.name, notebook_node.id.substring(8), notebook_name);
					}
				}));
			}
			menu.addChild(new MenuItem({
				label: 'Delete',
				onClick: function() {
					coge_track_list._delete(track_config, type, id, container);
				}
			}));
		}
		var btn = new DropDownButton({ dropDown: menu }, menu_button);
		menu.startup();
		btn.startup();
	},

	// ----------------------------------------------------------------

	_mouse_leave: function() {
		if (!this._menu_popped_up()) {
			var b = dijit.byId('coge_track_menu_button');
			if (b)
				b.destroy();
		}
	},

	// ----------------------------------------------------------------

	moveTracks: function( source, nodes, copy, target ) {
		if (source == target) // dropping in same place dragged from. should only happen in jbrowse track container
			return;
		var target_is_in_selector = target.node.firstChild.config;
		if (!target_is_in_selector) // dragging a track from the selector onto jbrowse's track container
			this.browser.publish('/jbrowse/v1/v/tracks/show', nodes.map(function(n){ return source.map[n.id].data; }));
		else {
			var items = [];
			nodes.forEach(function(node) {
				var n = target.node.firstChild;
				while (n.id != node.id)
					n = n.nextSibling;
				target.node.removeChild(n);
				items.push(node.config);
			}, this);
			this._add_to_notebook(items, target.node.firstChild.config.coge.id, true);
		}
	},

	// ----------------------------------------------------------------

	_new_notebook_config: function(id, name, description, restricted) {
		return {
			key: (restricted ? '&reg; ' : '') + name,
			baseUrl: api_base_url + '/experiment/notebook/' + id + '/',
			autocomplete: 'all',
			track: 'notebook' + id,
			label: 'notebook' + id,
			type: 'CoGe/View/Track/Wiggle/XYPlot',
			storeClass: 'JBrowse/Store/SeqFeature/REST',
			style: { featureScale: 0.001 },
			showHoverScores: 1,
			coge: {
				id: id,
				type: 'notebook',
				collapsible: true,
				name: name,
				description: description,
				editable: true,
				experiments: null,
				onClick: 'NotebookView.pl?embed=1&lid=' + id,
				menuOptions: [{
					label: 'NotebookView',
					action: "function() { window.open( 'NotebookView.pl?lid=" + id + "' ); }"
				}]
			}
		};
	},

	// ----------------------------------------------------------------

	_new_notebook_source: function() {
		var div = dojo.create( 'div', null, this.tracks_div );
		return new DndSource(div, {
			accept: ["track"],
			checkAcceptance: function(source, nodes) {
				for (var i=0; i<nodes.length; i++) {
					if (!nodes[i].config) // only accept experiments from the track selector (not from jbrowse's track container)
						return false;
					if (nodes[i].config.coge.type == 'notebook') {
						dojo.dnd.manager().stopDrag();
						return false;
					}
				}
				var container = this.node.firstChild;
				if (container.id == 'notebook0') // "All Experiments"
					return false;
				if (container.id == source.node.firstChild.id) // same notebook
					return false;
				if (!container.config.coge.editable)
					return false;
				for (var i=0; i<this.node.children.length; i++)
					if (this.node.children[i].id == nodes[0].id) // already has experiment
						return false;
				return true;
			},
			copyOnly: true,
			creator: dojo.hitch(this, function(track_config, hint) {
				return {node: this._new_track(track_config, track_config.coge.type == 'experiment'), data: track_config, type: ["track", track_config.coge.type]};
			}),
			delay: 2,
			// onDrop: function(source, nodes, copy) {
			// 	var experiments = Array.prototype.slice.call(div.childNodes, 0);
			// 	experiments.shift();
			// 	experiments.sort(
			// 	    function( a,b ) {  
			// 	        return natural_sort(a.config.coge.name, b.config.coge.name);
			// 	    }
			// 	).forEach(
			// 	    function(a, idx) { 
			// 	        div.insertBefore(a, div.childNodes[idx]);  
			// 	});
			// },
			selfAccept: false
		});
	},

	// ----------------------------------------------------------------

	_new_track: function(track_config, hide) {
		var coge = track_config.coge;
		var container = dojo.create( 'div', {
			className: 'coge-track',
			id: coge.type + coge.id,
		});
		container.config = track_config;
		var label = dojo.create('div', {
			className: 'coge-track-label coge-' + coge.type
		});
		if (coge.type == 'experiment' || coge.type == 'features')
			dojo.addClass(label, 'coge-track-indented');

		this._set_track_label(track_config, label);

		dojo.connect(label, "click", dojo.hitch(this, function() {
			if (track_config.coge.selected)
				this.browser.publish( '/jbrowse/v1/v/tracks/hide', [track_config] );
			else
				this.browser.publish( '/jbrowse/v1/v/tracks/show', [track_config] );
		}));

		if (coge.type == 'experiment' || coge.type == 'notebook') {
			dojo.connect(container, "onmouseenter", dojo.hitch(this, function(){this._mouse_enter(track_config, label, container)}));
			dojo.connect(container, "onmouseleave", dojo.hitch(this, this._mouse_leave));
		}
		if (coge.collapsible)
			this._add_expander(container);
		else if (hide)
			container.style.display = 'none';

		container.appendChild(label);
		return container;
	},

	// ----------------------------------------------------------------

	_reload_notebook: function(notebook_id) {
		var track = dojo.byId('track_notebook' + notebook_id);
		if (track) {
			var config = dojo.byId('notebook' + notebook_id).config;
			var d = new Deferred();
			var store;
			this.browser.getStore(config.store, function(s) {
				store = s;
				d.resolve(true);
       		});
       		d.promise.then(function() {
       			var next = track.nextSibling;
				store.browser.publish('/jbrowse/v1/v/tracks/hide', [config]);
				var cache = store._getCache();
				cache._prune(cache.maxSize);
				store.config.noCache = true;
				store.config.feature_range_cache = false;
				store.browser.publish('/jbrowse/v1/v/tracks/show', [config]);
				store.config.noCache = false;
				store.config.feature_range_cache = true;
				if (next) {
					dojo.place(dojo.byId('track_notebook' + notebook_id), next, 'before');
					store.browser.view.updateTrackList();
				}
			});
		}
	},

	// ----------------------------------------------------------------

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
					coge_plugin.error('Remove Notebook Item', data);
				else {
					dojo.destroy(this._menu_node.parentNode);
					this._reload_notebook(notebook_id);
					this._traverse_tracks(function(container){
						if (container.config.coge.type == 'experiment' && container.config.coge.id == id && dojo.hasClass(container.firstChild, 'coge-circle'))
							dojo.destroy(container.firstChild);
					})
				}
			}),
			error: function(data) {
				coge_plugin.error('Remove Notebook Item', data);
			}
		});
	},

	// ----------------------------------------------------------------

	_remove_from_notebook_dialog: function(type, id, name, notebook_id, notebook_name) {
		coge_plugin.confirm(
			'Remove ' + this._capitalize(type),
			'Remove ' + type + ' "' + name + '" from notebook "' + notebook_name + '"?',
			dojo.hitch(this, function(confirmed) {
				this._remove_from_notebook(type, id, notebook_id);
			})
		);
	},

	// ----------------------------------------------------------------

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
					coge_plugin.error('Rename ' + type, data);
				else {
					this._menu_track_config.coge.name = name;
					var key = this._menu_track_config.key;
					if (key.length > 6 && key.substring(0, 6) == '&reg; ') {
						key = '&reg; ' + name;
						old_name = '&reg; ' + old_name;
					} else
						key = name;
					this._menu_track_config.key = key;
					this._set_track_label(this._menu_track_config, this._menu_node);
					this._track_changed(old_name, key);
					this._rename_dialog.hide();
				}
			}),
			error: function(data) {
				coge_plugin.error('Rename ' + type, data);
			}
		});
	},

	// ----------------------------------------------------------------

	rename_dialog: function(type, id, name) {
		this._rename_dialog = new Dialog({
			title: 'Rename ' + type,
			content: '<table><tr><td><label>Name:</label></td><td><input id="name" value=' + JSON.stringify(name) + '></td></tr></table><div class="dijitDialogPaneActionBar"><button data-dojo-type="dijit/form/Button" type="button" onClick="coge_track_list._rename(\'' + type + '\',' + id + ',\'' + name.replace(/'/g, "\\'") + '\')">OK</button><button data-dojo-type="dijit/form/Button" type="button" onClick="coge_track_list._rename_dialog.hide()">Cancel</button></div>',
			onHide: function() {
				this.destroyRecursive();
				coge_track_list._rename_dialog = null;
			},
			style: "width: 300px"
		});
		this._rename_dialog.show();
		var i = dojo.byId("name");
		i.setSelectionRange(0, i.value.length);
	},

	// ----------------------------------------------------------------

	replaceTracks: function( track_configs ) { // mdb: unused now
		console.log('replaceTracks');
	},

	// ----------------------------------------------------------------

	_set_track_label: function(track_config, node) {
		var name = track_config.key;
		node.innerHTML = name;
		node.title = this._build_title(track_config);
	},

	// ----------------------------------------------------------------

	_set_tracks_active: function(container, combined) {
		var label_node = container.firstChild;
		while (!dojo.hasClass(label_node, 'coge-track-label'))
			label_node = label_node.nextSibling;
		if (!combined)
			container.config.coge.selected = true;
		if (container.config.coge.type == 'experiment') {
			if (combined)
				dojo.create('div', { className: 'coge-circle', style: { backgroundColor: coge_track_list._get_track_color(container) } }, container, 'first');
			else
				dojo.style(label_node, 'background', coge_track_list._get_track_color(container));
		} else if (container.config.coge.type == 'notebook') {
			dojo.style(label_node, 'background', 'lightgray');
			var experiment = container.nextSibling;
			while (experiment) {
				this.setTracksActive([experiment.config], true);
				experiment = experiment.nextSibling;
			}
		} else
			dojo.style(label_node, 'background', 'lightgray');
	},

	// ----------------------------------------------------------------

	_set_tracks_inactive: function(container, combined) {
		if (combined) {
			if (dojo.hasClass(container.firstChild, 'coge-circle'))
				dojo.destroy(container.firstChild);
			return;
		}
		var label_node = container.firstChild;
		while (!dojo.hasClass(label_node, 'coge-track-label'))
			label_node = label_node.nextSibling;
		container.config.coge.selected = false;
		dojo.style(label_node, 'background', '');
		if (container.config.coge.type == 'notebook') {
			var experiment = container.nextSibling;
			while (experiment) {
				this.setTracksInactive([experiment.config], true);
				experiment = experiment.nextSibling;
			}
		}
	},

	// ----------------------------------------------------------------
	/**
	 * Given an array of track configs, update the track list to show that they
	 * are turned on.
	 */

	setTracksActive: function(track_configs, combined) {
		track_configs.forEach(function(track_config) {
			coge_track_list._traverse_tracks(function(container) {
				if (track_config.coge && container.id == track_config.coge.type + track_config.coge.id)
					coge_track_list._set_tracks_active(container, combined);
			});
		});
	},

	// ----------------------------------------------------------------
	/**
	 * Given an array of track configs, update the track list to show that they
	 * are turned off.
	 */

	setTracksInactive: function(track_configs, combined) {
		var search_tracks = [];
		track_configs.forEach(function(track_config) {
			coge_track_list._traverse_tracks(function(container) {
				if (track_config.coge && container.id == track_config.coge.type + track_config.coge.id)
					if (container.config.coge.search_track)
						search_tracks.push(container);
					else
						coge_track_list._set_tracks_inactive(container, combined);
			});
		});
		search_tracks.forEach(function(container) {
			coge_track_list.browser._deleteTrackConfigs([container.config]);
			dojo.destroy(container);
		});
	},

	// ----------------------------------------------------------------

	/**
	 * Make the track selector visible. This does nothing for the Simple track
	 * selector, since it is always visible.
	 */
	show: function() {
	},

	// ----------------------------------------------------------------

	/**
	 * Toggle visibility of this track selector. This does nothing for the
	 * Simple track selector, since it is always visible.
	 */
	toggle: function() {
	},
	
	// ----------------------------------------------------------------

	_track_changed: function(key, new_key) {
		this.browser.view.tracks.forEach(function(track) {
			if (track.config.key == key) {
				if (new_key) {
					track.config.key = new_key;
					track.label.children[1].innerHTML = new_key;
				}
				track.changed();
			}
		});
	},

	// ----------------------------------------------------------------

	_traverse_tracks: function(f) {
		var n = this.tracks_div.firstChild;
		while (n && dojo.hasClass(n, 'coge-track')) {
			f(n);
			n = n.nextSibling;
		}
		while (n) {
			var c = n.firstChild;
			while (c) {
				f(c);
				c = c.nextSibling;
			}
			n = n.nextSibling;
		}
	},

	// ----------------------------------------------------------------
	// show cancel button only when text filter has text

	_update_text_filter_control: function() {
		if (this.text_filter_input.value.length)
			dojo.setStyle(this.text_filter_cancel, 'display', 'flex');
		else
			dojo.setStyle(this.text_filter_cancel, 'display', 'none');
	},

	// ----------------------------------------------------------------

	_update_tracks_shown: function() {
		var num_experiments = 0;
		var num_notebooks = 0;
		var total_experiments = 0;
		var total_notebooks = 0;
		var all = dojo.byId('notebook0');
		if (all) {
			total_experiments = all.parentNode.children.length - 1;
			var n = all.parentNode;
			while (n) {
				var c = n.firstChild;
				while (c) {
					if (c.config.coge.type == 'experiment') {
						if (c.style.display != 'none')
							++num_experiments;
					} else {
						if (c.style.display != 'none')
							++num_notebooks;
						++total_notebooks;
					}
					c = c.nextSibling;
				}
				n = n.nextSibling;
			}
		}
		var html = Math.min(num_experiments, total_experiments) + ' of ' + total_experiments + ' experiment' + (total_experiments == 1 ? '' : 's') + ' shown<br>';
		if (this._type_filter) {
			html += 'showing ' + (this._type_filter === 1 ? 'quantitative' : this._type_filter === 2 ? 'polymorphism' : this._type_filter === 3 ? 'alignment' : 'marker') + ' tracks&nbsp;&nbsp;&nbsp;&nbsp;';
			html += '<span class="glyphicon glyphicon-remove" style="color:black;cursor:pointer" onclick="coge_track_list._type_filter=null;coge_track_list._filter_tracks(coge_track_list.text_filter_input.value);"></span>';
		}
		else
			html += num_notebooks + ' of ' + total_notebooks + ' notebook' + (total_notebooks == 1 ? '' : 's') + ' shown';
		this._tracks_shown.innerHTML = html;
	}
});
});
