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

    info: function(tite, content) {
    	new InfoDialog({
    		title: title,
    		content: content,
            onHide: function(){this.destroyRecursive()}
    	}).show();	
    }
});
});
