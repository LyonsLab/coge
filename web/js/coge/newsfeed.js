/* global window, document, coge*/
var coge = window.coge = (function(ns) {
    ns.newsfeed = {
    	init: function(opts) {
    		if (!opts.element_id) {
    			console.error("Missing element_id param");
    		}
    		this.element = $('#'+opts.element_id);
    	},
    	
    	fetch_logs: function(id, type) {
			coge.services.fetch_logs(id, type)
				.done(function(result) { // success
					if (result && result.entries) {
					}
				})
				.fail(function() { // error
					//TODO
				});
		}
    };

    return ns;
})(coge || {});
