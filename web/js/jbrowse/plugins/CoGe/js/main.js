var coge;
define([
           'dojo/_base/declare',
           'JBrowse/View/ConfirmDialog',
           'JBrowse/View/InfoDialog',
           'JBrowse/Plugin'
       ],
       function(
           declare,
           ConfirmDialog,
           InfoDialog,
           JBrowsePlugin
       ) {
return declare( JBrowsePlugin,
{
    constructor: function( args ) {
    	coge = this;
    	this.browser = args.browser;
    },

    // ----------------------------------------------------------------

    adjust_nav: function(eid) {
     	var l = dojo.byId('label_experiment' + eid);
     	if (l) {
     		var nav = dojo.byId('nav_' + eid);
     		if (nav) {
    	 		var track = dojo.byId('track_experiment' + eid);
    	     	dojo.style(nav, 'left', dojo.style(l, 'left') + 10);
    	     	dojo.style(nav, 'top', dojo.style(track, 'top') + 32);
     		}
     	}
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

    new_nav: function(eid, data) {
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
        this.browser.subscribe('/jbrowse/v1/v/tracks/hide', function(configs) {
        	for (var i=0; i<configs.length; i++)
        		if (configs[i].coge.id == eid)
        			dojo.destroy(dojo.byId('nav_' + eid));
        });
        coge.go_to_hit(nav, 0);
    }
});
});
