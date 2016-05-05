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
	var SearchNav = declare(null, {
		constructor: function(eid, results, browser) {
			this.results = results;
			this.browser = browser;
			this.hit = 0;
			this.div = dojo.create('div', { id: 'nav_' + eid, style: { background: 'white', opacity: 0.7, position: 'absolute' } }, dojo.byId('container'));
			coge.adjust_nav(eid);
			dojo.create('span', { className: 'glyphicon glyphicon-step-backward', onclick: dojo.hitch(this, function() { this.go_to(0) }), style: { cursor: 'pointer' } }, this.div);
			dojo.create('span', { className: 'glyphicon glyphicon-chevron-left', onclick: dojo.hitch(this, function() { if (this.hit > 0) this.go_to(this.hit - 1) }), style: { cursor: 'pointer' } }, this.div);
			this.num_span = dojo.create('span', { innerHTML: '1', style: { cursor: 'default' } }, this.div);
			dojo.create('span', { innerHTML: ' of ' + results.hits.length + ' hit' + (results.hits.length != 1 ? 's ' : ' '), style: { cursor: 'default', marginRight: 5 } }, this.div);
			dojo.create('span', { className: 'glyphicon glyphicon-chevron-right', onclick: dojo.hitch(this, function() { if (this.hit < this.results.hits.length - 1) this.go_to(this.hit + 1) }), style: { cursor: 'pointer' } }, this.div);
			dojo.create('span', { className: 'glyphicon glyphicon-step-forward', onclick: dojo.hitch(this, function() { this.go_to(this.results.hits.length - 1) }), style: { cursor: 'pointer' } }, this.div);
	        browser.subscribe('/jbrowse/v1/v/tracks/hide', function(configs) {
	        	for (var i=0; i<configs.length; i++)
	        		if (configs[i].coge.search_track && configs[i].coge.id == eid) {
	        			dojo.destroy(dojo.byId('nav_' + eid));
	        			return;
	        		}
	        });
		},
	    go_to: function(index) {
	    	this.hit = index % this.results.hits.length;
	    	var hit = this.results.hits[this.hit];
	    	this.num_span.innerHTML = this.hit + 1;
	    	this.num_span.title = JSON.stringify(hit);
	    	var chr = this.results.chr_at(this.hit);
	    	if (chr != this.browser.refSeq.name)
	    		this.browser.navigateToLocation({
	    			ref: chr,
	    			start: hit[0],
	    			end: hit[1]
	    		});
	    	else
	    		this.browser.view.centerAtBase((hit[0] + hit[1]) / 2, true);
	    }
	});

	// ----------------------------------------------------------------

	var SearchResults = declare(null, {
		constructor: function(data, stranded) {
			this.hits = data;
			this.stranded = stranded;
			this.chr = [];
	    	var current_chr;
	    	for (var i=0; i<data.length; i++) {
	    		var index = data[i].indexOf('"', 1);
	    		var chr = data[i].substring(1, index);
	    		if (chr != current_chr) {
	    			if (this.chr.length > 0)
	    				this.chr[this.chr.length - 1][1] = i;
	    			this.chr.push([chr, 0]);
	    			current_chr = chr;
	    		}
	    		this.hits[i] = JSON.parse('[' + data[i].substring(index + 2) + ']');
	    	}
			this.chr[this.chr.length - 1][1] = data.length;
		},
	    boundaries: function(chr) {
	    	var l = 0;
	    	for (var i=0; i<this.chr.length; i++) {
	    		if (chr == this.chr[i][0])
	    			return [l, l + this.chr[i][1]];
	    		l += this.chr[i][1];
	    	}
	    },
	    chr_at: function(index) {
	    	var l = 0;
	    	for (var i=0; i<this.chr.length; i++) {
	    		l += this.chr[i][1];
	    		if (l > index)
	    			return this.chr[i][0];
	    	}
	    },
	    get_hits: function(chr, start, end) {
	    	var b = this.boundaries(chr);
	    	if (!b)
	    		return null;
	    	var i = b[0];
	    	var j = b[1];
	    	if (this.hits[i][0] > end)
	    		return null;
	    	while (i < j && this.hits[i][1] < start)
	    		++i;
	    	--j;
	    	while (j >= i && this.hits[j][0] > end)
	    		--j;
	    	return [this.hits, i, j];
	    }
	});

	// ----------------------------------------------------------------

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

    build_chromosome_select: function(first, onchange) {
        var chr = this.browser.refSeq.name;
        var html = '<select id="coge_ref_seq"';
        if (onchange)
            html += ' onchange="' + onchange + '"';
        html += '>';
        if (first)
            html += '<option>' + first + '</option>';
        this.browser.refSeqOrder.forEach(function(rs) {
            html += '<option';
            if (rs == chr)
                html += ' selected';
            html += '>';
            html += rs;
            html += '</option>';
        });
        html += '</select>';
        return html;
    },

    // ----------------------------------------------------------------

    build_features_checkboxes: function() {
        var html = '';
        var features = this.browser.config.tracks.reduce(function(accum, current) {
            if (current.coge.type && current.coge.type == 'features' && accum.indexOf(current.coge.id) < 0)
                accum.push(current.coge.id);
            return accum;
        }, []);
        features.forEach(function(f) {
            html += '<div><input type="checkbox"';
            if (f == 'gene')
                html += ' checked';
            html += '> <label>' + f + '</label></div>';
        });
        html += '<div><button onClick="coge.check_all(this.parentNode.parentNode.parentNode,true)">check all</button> <button onClick="coge.check_all(this.parentNode.parentNode.parentNode,false)">uncheck all</button></div>';
        return html;
    },

    // ----------------------------------------------------------------

    calc_color: function(id) {
        return '#' + ((((id * 1234321) % 0x1000000) | 0x444444) & 0xe7e7e7 ).toString(16);
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

    create_download_dialog: function(track) {
        this._track = track;
        var content = '<div id="coge-track-download"><table align="center"><tr><td>Chromosome:</td><td>';
        content += this.build_chromosome_select('All');
        content += '</td></tr>';
        if (track.config.coge.transform)
            content += '<tr><td>Transform:</td><td style="white-space:nowrap"><input type="radio" name="transform" checked="checked"> None <input id="transform" type="radio" name="transform"> ' + track.config.coge.transform + '</td></tr>';
        if (track.config.coge.search)
            content += '<tr><td>Search:</td><td style="white-space:nowrap"><input type="radio" name="search" checked="checked"> None <input id="search" type="radio" name="search"> ' + coge.search_to_string(track.config.coge.search, true) + '</td></tr>';
        content += '<tr><td></td><td></td></tr></table><div class="dijitDialogPaneActionBar"><button data-dojo-type="dijit/form/Button" type="button" onClick="coge.download_track()">OK</button><button data-dojo-type="dijit/form/Button" type="button" onClick="coge._track_download_dialog.hide()">Cancel</button></div></div>';
        this._track_download_dialog = new Dialog({
            title: 'Download Track',
            content: content,
            onHide: function() {
                this.destroyRecursive();
                coge._track_download_dialog = null;
            },
            style: "width: 350px"
        });
        this._track_download_dialog.show();
    },

    // ----------------------------------------------------------------

    create_search_button: function() {
    	var content = '<div id="coge-search-dialog"><table><tr><td>Name:</td><td><input id="coge_search_text"></td></tr><tr><td>Chromosome:</td><td>';
        content += this.build_chromosome_select('Any');
    	content += '</td></tr><tr><td style="vertical-align:top;">Features:</td><td id="coge_search_features">';
        content += this.build_features_checkboxes();
    	content += '</td></tr></table>';
    	content += '<div class="dijitDialogPaneActionBar"><button data-dojo-type="dijit/form/Button" type="button" onClick="coge.search_features()">OK</button><button data-dojo-type="dijit/form/Button" type="button" onClick="coge.search_dialog.hide()">Cancel</button></div></div>';
        new Button({
	        label: 'Find Features',
	        onClick: function(event) {
	        	coge.search_dialog = new Dialog({
                    title: "Search",
                    content: content,
                    onHide: function() {
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

    download_track: function() {
        var ref_seq = dojo.byId('coge_ref_seq');
        var url = api_base_url + '/experiment/' + this._track.config.coge.id + '/data/' + ref_seq.options[ref_seq.selectedIndex].innerHTML + '?username='+un;
        if (dojo.byId('search') && dojo.byId('search').checked)
            url += '&' + coge.search_to_params(this._track.config.coge.search, true);
        if (dojo.byId('transform') && dojo.byId('transform').checked)
            url += '&transform=' + this._track.config.coge.transform;
        document.location = url;
        coge._track_download_dialog.hide();
    },

    // ----------------------------------------------------------------

    error: function(title, content) {
    	if (content.responseText) {
    		var error = JSON.parse(content.responseText);
    		if (error.error)
        		if (error.error.Error)
        			content = error.error.Error;
        		else
        			content = JSON.stringify(error.error);
    		else
    			content = content.responseText;
    	} else if (content.error)
    		if (content.error.Error)
    			content = content.error.Error;
    		else
    			content = JSON.stringify(content.error);
    	this.info(title, content);
    },

    // ----------------------------------------------------------------

    get_checked_values: function(id, description, quote) {
        var checkboxes = document.getElementById(id).getElementsByTagName('INPUT');
        var values = [];
        for (var i=0; i<checkboxes.length; i++)
            if (checkboxes[i].checked)
                if (quote)
                    values.push("'" + checkboxes[i].nextElementSibling.innerText + "'");
                else
                    values.push(checkboxes[i].nextElementSibling.innerText);
        if (!values.length) {
            coge.error('Search', 'Please select one or more ' + description + ' to search.');
            return null;
        }
        return values.length == checkboxes.length ? 'all' : values.join();
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

    new_search_track: function(track, data) {
    	var browser = this.browser;
    	var config = track.config;
    	var eid = config.coge.id;
    	if (dojo.byId('nav_' + eid))
    		browser.publish('/jbrowse/v1/v/tracks/hide', [coge_track_list.get_search_config(eid)]);
    	var results = new SearchResults(data);
//        var d = new Deferred();
        var store_config = {
            browser: browser,
            config: config,
            refSeq: browser.refSeq,
            results: results,
            type: 'CoGe/Store/SeqFeature/Search'
        };
        var store_name = browser.addStoreConfig(undefined, store_config);
        store_config.name = store_name;
        browser.getStore(store_name, function(store) {
//            d.resolve(true);
//        });
//        d.promise.then(function() {
        	config = dojo.clone(config);
        	config.key = 'Search: ' + config.key + ' (' + coge.search_to_string(track.config.coge.search) + ')';
        	config.track = 'search_' + eid;
        	config.label = 'search_' + eid;
            config.metadata = {Description: 'Track to show results of searching a track.'};
            config.store = store_name;
            config.coge.collapsed = false;
            config.coge.search_track = true;
            browser.publish('/jbrowse/v1/v/tracks/new', [config]);
            browser.publish('/jbrowse/v1/v/tracks/show', [config]);
   			dojo.place(dojo.byId('track_search_' + eid), dojo.byId('track_experiment' + eid), 'after');
   			browser.view.updateTrackList();
    		new SearchNav(eid, results, browser).go_to(0);
        });
    },

    // ----------------------------------------------------------------

	search_features: function() {
        var types = this.get_checked_values('coge_search_features', 'feature types', true);
        if (!types)
            return;

    	var name = encodeURIComponent(dojo.byId('coge_search_text').value);
		var url = api_base_url + '/genome/' + gid + '/features?name=' + name + '&features=' + types;
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
                coge.search_dialog.hide();
	  			if (data.error) {
	  				coge.error('Search', data);
	  				return;
	  			}
	  			if (data.length == 0) {
	  				coge.error('Search', 'no features found');
	  				return;
	  			}
	  			//dojo.query('.dijitDialogUnderlayWrapper')[0].style.display = 'none';
  				//var div = dojo.byId('coge-search-dialog');
  				//div.style.maxHeight = '500px';
  				//div.style.overflow = 'auto';
    			//dojo.empty(div);
                var div = dojo.byId('feature_hits')
                dojo.create('div', { innerHTML: 'Features <span class="glyphicon glyphicon-remove" onclick="dojo.byId(\'feature_hits\').empty();dijit.byId(\'jbrowse\').resize()"></span>' }, div);
                div = dojo.create('div', { 'class': 'feature_hits' }, div);
    			data.forEach(function(hit) {
    				dojo.create('a', {
    					innerHTML: hit.name,
    					onclick: dojo.hitch(hit, function() {
	    					coge.browser.navigateToLocation(this.location);
	    					return false;
	    				})
    				}, div);
    				dojo.create('br', null, div);
    			});
                dijit.byId('jbrowse').resize();
    		},
    		error: function(data) {
    			coge.error('Search', data);
    		}
    	})
    },

    // ----------------------------------------------------------------

	search_to_params: function(search, without_chr) {
		var params;
        if (search.type == 'SNPs')
            if (search.snp_type)
                params = 'snp_type=' + search.snp_type;
            else
                params = 'features=' + search.features;
        else if (search.type == 'range')
			params = 'type=range&gte=' + search.gte + '&lte=' + search.lte;
		else
			params = 'type=' + search.type;
		if (!without_chr && search.chr && search.chr != 'Any')
			params += '&chr=' + search.chr;
		return params;		
	},

    // ----------------------------------------------------------------

	search_to_string: function(search, without_chr) {
		var string;
		if (search.type == 'range')
			string = 'range: ' + search.gte + ' .. ' + search.lte;
		else if (search.type == 'SNPs') {
            if (search.snp_type)
                string = search.snp_type;
            else {
    			string = 'SNPs';
    			if (search.features != 'all')
    				string += ' in ' + search.features;
            }
		} else
			string = search.type;
		if (!without_chr && search.chr && search.chr != 'Any')
			string += ', chr=' + search.chr;
		return string;
	}
});
});
