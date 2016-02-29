var coge;
define([
           'dojo/_base/declare',
           'dojo/Deferred',
           'dijit/Dialog',
           'dijit/form/Button',
           'JBrowse/View/ConfirmDialog',
           'JBrowse/View/InfoDialog',
           'JBrowse/Plugin'
       ],
       function(
           declare,
           Deferred,
           Dialog,
           Button,
           ConfirmDialog,
           InfoDialog,
           JBrowsePlugin
       ) {
return declare( JBrowsePlugin,
{
    constructor: function( args ) {
    	coge = this;
    	this.browser = args.browser;
    	JBrowse.afterMilestone('initView', function() {
    		coge.create_search_button();
    	});
    },

    // ----------------------------------------------------------------

    adjust_nav: function(eid) {
     	var l = dojo.byId('label_search_' + eid);
     	if (l) {
     		var nav = dojo.byId('nav_' + eid);
     		if (nav) {
    	 		var track = dojo.byId('track_search_' + eid);
    	     	dojo.style(nav, 'left', dojo.style(l, 'left') + 10);
    	     	dojo.style(nav, 'top', dojo.style(track, 'top') + 26);
     		}
     	}
    },

    // ----------------------------------------------------------------

    check_all: function(element, value) {
    	var cb = element.getElementsByTagName('INPUT');
    	for (var i=0; i<cb.length; i++)
    		cb[i].checked = value;
    },

    // ----------------------------------------------------------------

    confirm: function(title, message, on_confirmed) {
		new ConfirmDialog({
			title: title,
			message: message,
            onHide: function(){this.destroyRecursive()}
		}).show(function(confirmed) {
             if (confirmed)
            	 on_confirmed();
        });
    },

    // ----------------------------------------------------------------

    create_search_button: function() {
    	var content = '<div id="coge-search-dialog"><table><tr><td>Name:</td><td><input id="coge_search_text"></td></tr><tr><td>Chromosome:</td><td><select id="coge_ref_seq"><option>Any</option>';
    	this.browser.refSeqOrder.forEach(function(rs){
    		content += '<option>' + rs + '</option>';
    	})
    	content += '</select></td></tr><tr><td style="vertical-align:top;">Features:</td><td id="coge_search_features">';
    	var features = this.browser.config.tracks.filter(function(f) {
    		return (f.coge.type && f.coge.type == 'features');
    	});
    	features.forEach(function(t) {
    		content += '<div><input type="checkbox"> <label>' + t.coge.id + '</label></div>';
    	});
    	content += '</td></tr><tr><td></td><td><button onClick="coge.check_all(this.parentNode.parentNode.parentNode,true)">check all</button> <button onClick="coge.check_all(this.parentNode.parentNode.parentNode,false)">uncheck all</button></td></tr></table>';
    	content += '<div class="dijitDialogPaneActionBar"><button data-dojo-type="dijit/form/Button" type="button" onClick="coge.search_features()">OK</button><button data-dojo-type="dijit/form/Button" type="button" onClick="coge.search_dialog.hide()">Cancel</button></div></div>';
        new Button({
	        label: 'Search',
	        onClick: function(event) {
	        	coge.search_dialog = new Dialog({
                    title: "Search",
                    content: content,
                    onHide: function(){
                    	this.destroyRecursive();
                    	coge.search_dialog = null;
                    },
                    style: "width: 300px"
                });
	        	coge.search_dialog.show();
	        	dojo.stopEvent(event);
	        },
	    }, dojo.create('button', null, this.browser.navbox));
    },

    // ----------------------------------------------------------------

    error: function(title, content) {
    	if (content.responseText)
    		content = content.responseText;
    	else if (content.error)
    		content = JSON.stringify(content.error);
    	this.info(title, content);
    },

    // ----------------------------------------------------------------

    calc_color: function(id) {
    	return '#' + ((((id * 1234321) % 0x1000000) | 0x444444) & 0xe7e7e7 ).toString(16);
    },

    // ----------------------------------------------------------------

    go_to_hit: function(nav, hit) {
    	nav.hit = hit % nav.hits.length;
    	nav.num_span.innerHTML = nav.hit + 1;
    	nav.num_span.title = nav.hits[nav.hit];
    	var loc = JSON.parse('[' + nav.hits[nav.hit] + ']');
    	if (loc[0] != this.browser.refSeq.name)
    		this.browser.navigateToLocation({
    			ref: loc[0],
    			start: loc[1],
    			end: loc[loc.length > 2 ? 2 : 1]
    		});
    	else
    		this.browser.view.centerAtBase(loc.length > 2 ? (loc[1] + loc[2]) / 2 : loc[1], true);
    },

    // ----------------------------------------------------------------

    info: function(title, content) {
    	new InfoDialog({
    		title: title,
    		content: content,
            onHide: function(){this.destroyRecursive()}
    	}).show();	
    },

    // ----------------------------------------------------------------

    new_search_track: function(track, data, search) {
        var browser = this.browser;
    	var eid = track.config.coge.id;
        var d = new Deferred();
        var store_config = {
            browser: browser,
            data: data,
            eid: eid,
            refSeq: browser.refSeq,
            type: 'CoGe/Store/SeqFeature/Search'
        };
        var store_name = browser.addStoreConfig(undefined, store_config);
        store_config.name = store_name;
        browser.getStore(store_name, function(store) {
            d.resolve(true);
        });
        d.promise.then(function(){
        	var config = dojo.clone(track.config);
        	config.key = 'Search: ' + track.config.key + ' (' + coge.search_to_string(search) + ')';
        	config.track = 'search_' + eid;
        	config.label = 'search_' + eid;
            config.metadata = {Description: 'Track to show results of searching a track.'};
            config.store = store_name;
            config.coge.collapsed = false;
            config.coge.search_track = true;
            browser.publish( '/jbrowse/v1/v/tracks/new', [config] );
            browser.publish( '/jbrowse/v1/v/tracks/show', [config] );
        });
		var nav = dojo.byId('nav_' + eid);
		if (nav)
			dojo.destroy(nav);
		nav = dojo.create('div', { id: 'nav_' + eid, style: { background: 'white', opacity: 0.7, position: 'absolute' } }, dojo.byId('container'));
		coge.adjust_nav(eid);
		nav.hits = data;
		nav.hit = 0;
		dojo.create('span', { className: 'glyphicon glyphicon-step-backward', onclick: function() { coge.go_to_hit(nav, 0) }, style: { cursor: 'pointer' } }, nav);
		dojo.create('span', { className: 'glyphicon glyphicon-chevron-left', onclick: function() { if (nav.hit > 0) coge.go_to_hit(nav, nav.hit - 1) }, style: { cursor: 'pointer' } }, nav);
		nav.num_span = dojo.create('span', { innerHTML: '1', style: { cursor: 'default' } }, nav);
		dojo.create('span', { innerHTML: ' of ' + data.length + ' hit' + (data.length != 1 ? 's ' : ' '), style: { cursor: 'default', marginRight: 5 } }, nav);
		dojo.create('span', { className: 'glyphicon glyphicon-chevron-right', onclick: function() { if (nav.hit < nav.hits.length - 1) coge.go_to_hit(nav, nav.hit + 1) }, style: { cursor: 'pointer' } }, nav);
		dojo.create('span', { className: 'glyphicon glyphicon-step-forward', onclick: function() { coge.go_to_hit(nav, nav.hits.length - 1) }, style: { cursor: 'pointer' } }, nav);
        browser.subscribe('/jbrowse/v1/v/tracks/hide', function(configs) {
        	for (var i=0; i<configs.length; i++)
        		if (configs[i].coge.id == eid) {
        			dojo.destroy(dojo.byId('nav_' + eid));
        			return;
        		}
        });
        coge.go_to_hit(nav, 0);
    },

    // ----------------------------------------------------------------

	search_features: function() {
    	var features = document.getElementById('coge_search_features').getElementsByTagName('INPUT');
    	var types = [];
    	for (var i=0; i<features.length; i++)
    		if (features[i].checked)
    			types.push("'" + features[i].nextElementSibling.innerText + "'");
    	if (!types.length) {
    		coge.error('Search', 'Please select one or more feature types to search.');
    		return;
    	}

    	var name = encodeURIComponent(dojo.byId('coge_search_text').value);
		var url = api_base_url + '/genome/' + gid + '/features?name=' + name + '&features=' + (types.length == features.length ? 'all' : types.join());
    	var ref_seq = dojo.byId('coge_ref_seq');
    	if (ref_seq.selectedIndex > 0)
    		url += '&chr=' + ref_seq.options[ref_seq.selectedIndex].innerHTML;

    	var div = dojo.byId('coge-search-dialog');
		dojo.empty(div);
		div.innerHTML = '<img src="picts/ajax-loader.gif">';

    	dojo.xhrGet({
    		url: url,
    		handleAs: 'json',
	  		load: function(data) {
	  			if (data.error) {
	  				coge.error('Search', data);
	  				return;
	  			}
	  			if (data.length == 0) {
	  				coge.error('Search', 'no features found');
	  				return;
	  			}
	  			dojo.query('.dijitDialogUnderlayWrapper')[0].style.display = 'none';
  				var div = dojo.byId('coge-search-dialog');
  				div.style.maxHeight = '500px';
  				div.style.overflow = 'auto';
    			dojo.empty(div);
    			data.forEach(function(hit){
    				dojo.create('a', {
    					innerHTML: hit.name,
    					onclick: dojo.hitch(hit, function() {
	    					coge.browser.navigateToLocation(this.location);
	    					return false;
	    				})
    				}, div);
    				dojo.create('br', null, div);
    			});
    		},
    		error: function(data) {
    			coge.error('Search', data);
    		}
    	})
    },

    // ----------------------------------------------------------------

	search_to_params: function(search) {
		var params;
		if (search.type == 'range')
			params = 'type=range&' + search.gte + '&' + search.lte;
		else
			params = 'type=' + search.type;
		if (search.chr && search.chr != 'Any')
			params += '&chr=' + search.chr;
		return params;		
	},

    // ----------------------------------------------------------------

	search_to_string: function(search) {
		var string;
		if (search.type == 'range')
			string = 'range: ' + search.gte + ' .. ' + search.lte;
		else if (search.type == 'SNPs') {
			string = search.type;
			if (search.features != 'all')
				string += ' in ' + search.features;
		} else
			string = search.type;
		if (search.chr && search.chr != 'Any')
			string += ', chr=' + search.chr;
		return string;
	}
});
});
