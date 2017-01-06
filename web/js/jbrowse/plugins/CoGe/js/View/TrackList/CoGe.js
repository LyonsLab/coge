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
		'dijit/Dialog',
		'dijit/Tooltip'
	   ],
	   function( declare, dom, domGeom, aspect, ContentPane, DndSource, DropDownButton, Menu, MenuItem, MenuSeparator, Deferred, Dialog, Tooltip ) {
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
		this.browser.subscribe( '/dnd/drop', dojo.hitch( this, 'moveTracks' ));

		// subscribe to commands coming from the the controller
		this.browser.subscribe( '/jbrowse/v1/c/tracks/show', 	dojo.hitch( this, 'setTracksActive' ));
		this.browser.subscribe( '/jbrowse/v1/c/tracks/hide', 	dojo.hitch( this, 'setTracksInactive' ));
		this.browser.subscribe( '/jbrowse/v1/c/tracks/new', 	dojo.hitch( this, 'addTracks' ));
		this.browser.subscribe( '/jbrowse/v1/c/tracks/replace', dojo.hitch( this, 'replaceTracks' ));
		this.browser.subscribe( '/jbrowse/v1/c/tracks/delete',  dojo.hitch( this, 'deleteTracks' ));

		this.browser.addGlobalMenuItem('view', new MenuItem({
			label: 'Hide Track Rulers',
			onClick: function() {
				var ss = document.styleSheets;
				for (var i=0; i<ss.length; i++)
					if (ss[i].href && ss[i].href.endsWith('plugins/CoGe/css/main.css'))
						for (var j=0; j<ss[i].cssRules.length; j++) {
							var rule = ss[i].cssRules[j];
							if (rule.cssText.startsWith('.vertical_rule')) {
								rule.style.visibility = rule.style.visibility == 'hidden' ? '' : 'hidden';
								this.set('label', this.get('label').startsWith('Hide') ? 'Show Track Rulers' : 'Hide Track Rulers')
								return;
							}
						}
			}
		}));
		this.browser.addGlobalMenuItem('view', new MenuItem({
			label: 'Hide Track Gridlines',
			onClick: function() {
				dojo.byId('gridtrack').style.display = dojo.byId('gridtrack').style.display == 'none' ? '' : 'none';
				this.set('label', this.get('label').startsWith('Hide') ? 'Show Track Gridlines' : 'Hide Track Gridlines')
			}
		}));

		this.browser.addGlobalMenuItem('help', new MenuSeparator());
		this._add_help('Create experiment from search results', 'EPIC_CoGe_Reference#Create_New_Experiment_from_Search_Tracks');
		this._add_help('Managing experiments', 'EPIC_CoGe_Reference#Managing_Experiments');
		this._add_help('Managing notebooks', 'EPIC_CoGe_Reference#Managing_Notebooks');
		this._add_help('Search for features by name', 'EPIC_CoGe_Reference#Feature_Search');
		this._add_help('Search tracks', 'EPIC_CoGe_Reference#Search_Tracks');
		this._add_help('Settings', 'EPIC_CoGe_Reference#Settings');
		this._add_help('Track overlaps', 'EPIC_CoGe_Reference#Track_Overlaps');
		this._add_help('Videos', 'EPIC-CoGe_Videos');

		var splash = localStorage.getItem('EPIC-CoGe splash');
		if (!splash || 'yes' == splash)
			this._help('Welcome to EPIC-CoGe', 'splash.html', function(){localStorage.setItem('EPIC-CoGe splash', dojo.byId('splash').checked ? 'yes' : 'no')});

		setTimeout(function(){dijit.byId('jbrowse').resize();}, 100);
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

	_add_help: function(label, page) {
		this.browser.addGlobalMenuItem('help', new MenuItem({
			label: label,
			onClick: function(){ window.open('https://genomevolution.org/wiki/index.php/' + page, '_blank'); }
		}));
	},

	// ----------------------------------------------------------------

	add_to_notebook: function(track_configs, notebook_id) {
		var coge_api = api_base_url.substring(0, api_base_url.length - 8);
		var items = [];
		track_configs.forEach(function(track_config){
			items.push({type: track_config.coge.type, id: track_config.coge.id});
		});
		dojo.xhrPost({
			url: coge_api + '/notebooks/' + notebook_id + '/items/add?username=' + USER_NAME,
			postData: JSON.stringify({ items: items }),
			handleAs: 'json',
			load: dojo.hitch(this, function(data) {
				if (data.error) {
					on_error();
					coge_plugin.error('Add Notebook Item', data);
				} else {
					var n = dojo.byId('notebook' + notebook_id);
					track_configs.forEach(function(track_config) {
						var e = dojo.byId(track_config.track);
						this._add_track_to_notebook(e.config, n);
					}, this);
					if (n.style.display != 'none')
						this._expand(n);
				}
			}),
			error: function(data) {
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
		notebooks.sort(function(a, b) { return coge_plugin.natural_sort(a.coge.name, b.coge.name); });
		for (var i=0; i<notebooks.length; i++)
			content += '<option value="' + notebooks[i].coge.id + '">' + notebooks[i].coge.name + '</option>';
		content += '</select></td></tr></table>';
		content += coge_plugin.build_buttons('coge_track_list._add_dialog.onExecute()', 'coge_track_list._add_dialog.hide()');
		this._add_dialog = new Dialog({
			title: 'Add experiment ' + track_config.coge.name + ' to Notebook',
			onExecute: function() {
				var s=dojo.byId('coge_notebook');
				coge_track_list.add_to_notebook([track_config], s.options[s.selectedIndex].value);
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

	_add_track_to_notebook: function(track_config, notebook) {
		var track = this._new_track(track_config);
		if (notebook.config.coge.id != 0) {
			if (!notebook.config.coge.experiments)
				notebook.config.coge.experiments = [];
			notebook.config.coge.experiments.push({value:track.config.coge.id, label:track.config.coge.name})
		}
		if (!notebook.expanded)
			track.style.display = 'none';
		if (!notebook.nextSibling)
			dojo.place(track, notebook, 'after');
		else {
			var n = notebook.nextSibling;
			if (coge_plugin.natural_sort(n.config.coge.name, track.config.coge.name) > 0)
				dojo.place(track, n, 'before');
			else {
				while (n.nextSibling && coge_plugin.natural_sort(track.config.coge.name, n.nextSibling.config.coge.name) > 0)
					n = n.nextSibling;
				dojo.place(track, n, 'after');
			}
		}
		if (dojo.byId('track_experiment' + track_config.coge.id)) {
			var label_node = this._get_label_node(track);
			track.config.coge.selected = true;
			dojo.style(label_node, 'background', this._get_track_color(track.config, track_config.coge.id));
		}
		if (dojo.byId('track_notebook' + notebook.config.coge.id))
			dojo.create('div', { className: 'coge-circle', style: { backgroundColor: this._get_track_color(notebook.config, track_config.coge.id) } }, track, 'first');
		if (notebook.expanded)
			track.scrollIntoView();
		this._reload_notebook(notebook.config.coge.id);
	},

	// ----------------------------------------------------------------

	addTracks: function(track_configs) {
		track_configs.forEach(function(track_config) {
			if (track_config.coge)
				if (track_config.coge.type == 'merge' || track_config.coge.type == 'search')
					this.tracks_div.insertBefore(this._new_track(track_config), this.tracks_div.firstChild); // insert before Sequence track at top
				else if (track_config.coge.type == 'notebook') {
					this._track_configs.push(track_config);
					this._new_notebook_source().insertNodes(false, [track_config]);
					this._filter_tracks();
					this._move_track(dojo.byId('notebook' + track_config.coge.id).parentNode);
				} else
					this._add_track_to_notebook(track_config, dojo.byId('notebook0'));
		}, this);
	},

	// ----------------------------------------------------------------

	_build_tooltip: function(track_config) {
		var coge = track_config.coge;
		var html = coge.type == 'merge' ? 'Merge of Experiments' : coge.type == 'search' ? 'Search of Experiment' : this._capitalize(coge.type);
		html += ': <b>' + coge.name + '</b>';
		if (coge.id != 0)
			html += ' (id ' + (coge.type == 'merge' ? coge.eids : coge.type == 'search' ? coge.eid : coge.id) + ')';
		if (coge.description)
			html += '<br>' + coge.description;
		if (coge.type == 'experiment')
			html += '<br><i>' + (coge.data_type == 4 ? 'Markers' : coge.data_type == 3 ? 'Alignments' : coge.data_type == 2 ? 'Variant Data' : 'Quantitative Data') + '</i>';
		if (coge.annotations) {
			var a = coge.annotations.split('\n');
			html += '<hr><table style="max-width:500px;">';
			a.forEach(function(line, i){
				html += '<tr';
				if (i % 2)
					html += ' style="background:#eee"';
				html += '><td>';
				html += line.replace(':', '</td><td>');
				html += '</td></tr>';
			});
			html += '</table>';
		}
		return html;
	},

	// ----------------------------------------------------------------

	_capitalize: function (string) {
		return string.charAt(0).toUpperCase() + string.substring(1);
	},

	// ----------------------------------------------------------------

	_collapse: function(container) {
		container.expanded = false;
		container.firstChild.src = 'js/jbrowse/plugins/CoGe/img/arrow-right-icon.png';
		var n = container.nextSibling;
		while (n) {
			n.style.display = 'none';
			n = n.nextSibling;
		}
	},

	// ----------------------------------------------------------------

	_create_notebook_and_experiment_tracks: function(notebook_config, experiments) {
		var source = this._new_notebook_source();
		source.insertNodes(false, [notebook_config]);
		if (notebook_config.coge.id != 0)
			experiments = experiments.filter(function(e) {
				return e.coge.notebooks && dojo.indexOf(e.coge.notebooks, notebook_config.coge.id) != -1;
			});
		source.insertNodes(false, experiments.sort(function(a, b) { return coge_plugin.natural_sort(a.coge.name, b.coge.name); }));
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
		menu.addChild(new MenuItem({
			label: "Show Quantitative Tracks",
			onClick: dojo.hitch(this, function() {
				this._type_filter = 1;
				this._filter_tracks(this.text_filter_input.value);
			})
		}));
		menu.addChild(new MenuItem({
			label: "Show Variant Tracks",
			onClick: dojo.hitch(this, function() {
				this._type_filter = 2;
				this._filter_tracks(this.text_filter_input.value);
			})
		}));
		if (USER_NAME != 'public') {
			menu.addChild(new MenuSeparator());
			menu.addChild(new MenuItem({
				label: "Show My Tracks",
				onClick: dojo.hitch(this, function() {
					this._show = 'my';
					this._filter_tracks(this.text_filter_input.value);
				})
			}));
			menu.addChild(new MenuItem({
				label: "Show Shared Tracks",
				onClick: dojo.hitch(this, function() {
					this._show = 'shared';
					this._filter_tracks(this.text_filter_input.value);
				})
			}));
			menu.addChild(new MenuItem({
				label: "Show Public Tracks",
				onClick: dojo.hitch(this, function() {
					this._show = 'public';
					this._filter_tracks(this.text_filter_input.value);
				})
			}));
			menu.addChild(new MenuSeparator());
			menu.addChild(new MenuItem({
				label: "Create New Notebook",
				onClick: coge_plugin.create_notebook_dialog.bind(coge_plugin)
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

		var pane = dojo.create('div', { id: 'track_pane', style: 'width:250px;' }, root);
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
			return fg.coge.type && fg.coge.type == 'feature_group';
		});
		var features = this._track_configs.filter(function(f) {
			return f.coge.type && f.coge.type == 'features';
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
		this._notebooks = {}; // hash of notebook configs by coge id
		if (notebooks.length) {
			this._create_notebook_and_experiment_tracks(notebooks.shift(), experiments);
			notebooks.sort(function(a, b) { return coge_plugin.natural_sort(a.coge.name, b.coge.name); });
			notebooks.forEach(function(n) {
				this._notebooks[n.coge.id] = n;
				this._create_notebook_and_experiment_tracks(n, experiments);
			}, this);
		}

		// show all tracks
		this._filter_tracks();

		var root = dojo.create('div', { id: 'coge tooltip', style: 'background:#eee;display:none;padding:10px;position:absolute;white-space:pre-line;z-index:100;' }, dojo.byId('jbrowse'));
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
		container.expanded = true;
		container.firstChild.src = 'js/jbrowse/plugins/CoGe/img/arrow-down-icon.png';
		var n = container.nextSibling;
		while (n) {
			n.style.display = '';
			n = n.nextSibling;
		}
	},

	// ----------------------------------------------------------------

	_filter_track: function(container, text) {
		if (text) {
			var t = container.lastChild.title ? container.lastChild.title : container.lastChild.innerHTML;
			if (t.toLowerCase().indexOf(text) == -1)
				return false;
		}
		if (this._type_filter && this._type_filter != container.config.coge.data_type)
			return false;
		if (this._show) {
			if (container.config.coge.type != 'experiment')
				return false;
			if (this._show == 'my' && container.config.coge.role != 2)
				return false;
			if (this._show == 'shared' && !(container.config.coge.role == 3 || container.config.coge.role == 4))
				return false;
			if (this._show == 'public' && (container.config.coge.restricted || container.config.coge.role))
				return false;
		}
		return true;
	},

	// ----------------------------------------------------------------

	_filter_tracks: function(text) {
		if (text && /\S/.test(text))
			text = text.toLowerCase();
		else
			text = null;
		if (text || this._type_filter || this._show) {
			var already_shown = {};
			this._traverse_tracks(function(container) {
				if (coge_track_list._filter_track(container, text)) {
					if (!already_shown[container.id]) {
						container.style.display = '';
						already_shown[container.id] = true;
					} else
						container.style.display = 'none';
				} else
					container.style.display = 'none';
			});
		} else { // show all
			var expanded = true;
			this._traverse_tracks(function(container) {
				if (container.config.coge.collapsible) {
					expanded = container.expanded;
					container.style.display = '';
				} else
					container.style.display = expanded ? '' : 'none';
			});
		}
		this._update_tracks_shown();
	},

	// ----------------------------------------------------------------

	_get_label_node: function(container) {
		var label_node = container.firstChild;
		while (!dojo.hasClass(label_node, 'coge-track-label'))
			label_node = label_node.nextSibling;
		return label_node;
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

	_get_track_color: function(config, id) {
		if (config.coge.type == 'merge')
			return 'lightgray';
		var cookie = this.browser.cookie('track-style-' + config.track);
		if (cookie)
			config.style = dojo.fromJson(cookie);
		if (config.coge.type == 'search')
			id = config.coge.eid;
		if (config.style && config.style.featureColor && config.style.featureColor[id])
			return config.style.featureColor[id];
		return coge_plugin.calc_color(id);
	},

	// ----------------------------------------------------------------

	_help: function(title, file, hide) {
		dojo.xhrGet({
			url: 'js/jbrowse/plugins/CoGe/' + file,
			load: function(data) { coge_plugin.info(title, data, null, hide); }
		});
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
		Tooltip.show(this._build_tooltip(track_config), container);
	},

	// ----------------------------------------------------------------

	_mouse_leave: function(container) {
		if (!this._menu_popped_up()) {
			var b = dijit.byId('coge_track_menu_button');
			if (b)
				b.destroy();
		}
		Tooltip.hide(container);
	},

	// ----------------------------------------------------------------

	_move_track: function(div) {
		var type = div.firstChild.config.coge.type;
		var prev = div.previousSibling;
		while (prev.firstChild.id != 'notebook0' && (type == 'notebook' || prev.firstChild.config.coge.type != 'notebook') && coge_plugin.natural_sort(div.firstChild.config.coge.name, prev.firstChild.config.coge.name) < 0)
			prev = prev.previousSibling;
		if (prev != div.previousSibling) {
			dojo.place(div, prev, 'after');
			div.scrollIntoView();
			return;
		}
		var next = div.nextSibling;
		while (next && coge_plugin.natural_sort(div.firstChild.config.coge.name, next.firstChild.config.coge.name) > 0)
			next = next.nextSibling;
		if (next != div.nextSibling) {
			if (next)
				dojo.place(div, next, 'before');
			else if (div.nextSibling)
				dojo.place(div, div.parentNode, 'last');
		}
		div.scrollIntoView();
	},

	// ----------------------------------------------------------------

	moveTracks: function( source, nodes, copy, target ) {
		if (source == target) { // dropping in same place dragged from. should only happen in jbrowse track container
			if (!source.current || !copy)
				return;
			var track1 = $('#' + source.anchor.id.substring(6))[0];
			var track2 = $('#' + source.current.id.substring(6))[0];
			coge_plugin.dnd_dialog(track1, track2);
			return;
		}
		var target_is_in_selector = target.node.firstChild.config;
		if (!target_is_in_selector) { // dragging a track from the selector onto jbrowse's track container
			var data = nodes.map(function(n){ return source.map[n.id].data; });
			dojo.destroy(dojo.byId('track_' + data[0].coge.type + data[0].coge.id));
			this.browser.publish('/jbrowse/v1/c/tracks/show', data);
		} else {
			var items = [];
			nodes.forEach(function(node) {
				var n = target.node.firstChild;
				while (n.id != node.id)
					n = n.nextSibling;
				target.node.removeChild(n);
				items.push(node.config);
			}, this);
			this.add_to_notebook(items, target.node.firstChild.config.coge.id);
		}
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
			selfAccept: false
		});
	},

	// ----------------------------------------------------------------

	_new_track: function(track_config, hide) {
		var coge = track_config.coge;
		var container = dojo.create( 'div', {
			className: 'coge-track',
			id: track_config.track,
		});
		container.config = track_config;
		var label = dojo.create('div', {
			className: 'coge-track-label coge-' + coge.type
		});
		if (coge.type == 'experiment' || coge.type == 'features')
			dojo.addClass(label, 'coge-track-indented');

		label.innerHTML = track_config.key;

		dojo.connect(label, "click", dojo.hitch(this, function() {
			if (track_config.coge.selected)
				this.browser.publish( '/jbrowse/v1/v/tracks/hide', [track_config] );
			else
				this.browser.publish( '/jbrowse/v1/v/tracks/show', [track_config] );
		}));

		if (coge.type == 'experiment' || coge.type == 'notebook') {
			dojo.connect(container, "onmouseenter", dojo.hitch(this, function(){this._mouse_enter(track_config, label, container)}));
			dojo.connect(container, "onmouseleave", dojo.hitch(this, function(){this._mouse_leave(container)}));
		} else if (coge.type == 'merge') {
			dojo.connect(container, "onmouseenter", function(){Tooltip.show(track_config.coge.keys.join('<br>'), container)});
			dojo.connect(container, "onmouseleave", function(){Tooltip.hide(container)});
		} else if (coge.type == 'search') {
			dojo.connect(container, "onmouseenter", dojo.hitch(this, function(){Tooltip.show(this._build_tooltip(track_config), container)}));
			dojo.connect(container, "onmouseleave", function(){Tooltip.hide(container)});
		}
		if (coge.collapsible)
			this._add_expander(container);
		else if (hide)
			container.style.display = 'none';

		container.appendChild(label);
		return container;
	},

	// ----------------------------------------------------------------

	notebook_is_editable: function(notebook_id) {
		return this._notebooks[notebook_id].coge.editable;
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
			url: coge_api + '/notebooks/' + notebook_id + '/items/remove?username=' + USER_NAME,
			postData: JSON.stringify({
				items: [{
					type: type,
					id: id
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
			url: coge_api + '/' + type + 's/' + id + '?username=' + USER_NAME,
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
					if (key.length > 2 && key.substring(0, 2) == '🔒 ') {
						key = '🔒 ' + name;
						old_name = '🔒 ' + old_name;
					} else
						key = name;
					this._menu_track_config.key = key;
					this._menu_node.innerHTML = this._menu_track_config.key;
					this._track_changed(old_name, key);
					this._move_track(this._menu_node.parentNode.parentNode);
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
			content: '<table style="width:100%"><tr><td><label>Name:</label></td><td><input id="name" value=' + JSON.stringify(name) + ' style="width:100%"></td></tr></table>' +
				coge_plugin.build_buttons('coge_track_list._rename(\'' + type + '\',' + id + ',\'' + name.replace(/'/g, "\\'") + '\')', 'coge_track_list._rename_dialog.hide()'),
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

	replaceTracks: function(track_configs) { // mdb: unused now, but required to be here
		console.log('replaceTracks');
	},

	// ----------------------------------------------------------------

	setTrackColor: function(container_id, experiment_id, color) {
		var container = dojo.byId(container_id);
		if (!container.config.style.featureColor)
			container.config.style.featureColor = {};
		container.config.style.featureColor[experiment_id] = color;
		this.setTracksActive([container.config]);
	},

	// ----------------------------------------------------------------

	_set_tracks_active: function(container) {
		var label_node = this._get_label_node(container);
		container.config.coge.selected = true;
		var type = container.config.coge.type;
		if (type == 'experiment') {
			dojo.style(label_node, 'background', this._get_track_color(container.config, container.config.coge.id));
		} else if (type == 'notebook') {
			dojo.style(label_node, 'background', 'lightgray');
			var experiment = container.nextSibling;
			while (experiment) {
				if (dojo.hasClass(experiment.firstChild, 'coge-circle'))
					dojo.destroy(experiment.firstChild);
				dojo.create('div', { className: 'coge-circle', style: { backgroundColor: this._get_track_color(container.config, experiment.config.coge.id) } }, experiment, 'first');
				experiment = experiment.nextSibling;
			}
		} else if (type == 'search')
			dojo.style(label_node, 'background', this._get_track_color(container.config, container.config.coge.id));
		else
			dojo.style(label_node, 'background', 'lightgray');
	},

	// ----------------------------------------------------------------

	_set_tracks_inactive: function(container, combined) {
		if (combined) {
			if (dojo.hasClass(container.firstChild, 'coge-circle'))
				dojo.destroy(container.firstChild);
			return;
		}
		var label_node = this._get_label_node(container);
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

	setTracksActive: function(track_configs) {
		track_configs.forEach(function(track_config) {
			coge_track_list._traverse_tracks(function(container) {
				if (container.id == track_config.track)
					coge_track_list._set_tracks_active(container);
			});
		});
	},

	// ----------------------------------------------------------------
	/**
	 * Given an array of track configs, update the track list to show that they
	 * are turned off.
	 */

	setTracksInactive: function(track_configs, combined) {
		var tmp_tracks = [];
		track_configs.forEach(function(track_config) {
			coge_track_list._traverse_tracks(function(container) {
				if (container.id == track_config.track)
					if (container.config.coge.type == 'merge' || container.config.coge.type == 'search')
						tmp_tracks.push(container);
					else
						coge_track_list._set_tracks_inactive(container, combined);
			});
		});
		tmp_tracks.forEach(function(container) {
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
		if (this._type_filter || this._show) {
			html += 'showing ';
			if (this._show)
				html += this._show + ' ';
			if (this._type_filter)
				html += (this._type_filter === 1 ? 'quantitative' : this._type_filter === 2 ? 'variant' : this._type_filter === 3 ? 'alignment' : 'marker');
			html += ' tracks&nbsp;&nbsp;&nbsp;&nbsp;<span class="glyphicon glyphicon-remove" style="color:black;cursor:pointer" onclick="coge_track_list._type_filter=null;coge_track_list._show=null;coge_track_list._filter_tracks(coge_track_list.text_filter_input.value);"></span>';
		}
		else
			html += num_notebooks + ' of ' + total_notebooks + ' notebook' + (total_notebooks == 1 ? '' : 's') + ' shown';
		this._tracks_shown.innerHTML = html;
	}
});
});
