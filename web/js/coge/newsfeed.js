/* global window, document, coge*/
var coge = window.coge = (function(ns) {
    ns.newsfeed = {
    	init: function(opts) {
    		if (!opts.element_id) {
    			console.error("Missing required element_id param");
    			return;
    		}
    		if (!opts.type || !opts.id) {
    			console.error("Missing required type/id param");
    			return;
    		}
    		this.element = $('#'+opts.element_id);
    		this.type = opts.type;
    		this.id = opts.id;
    		
    		this.fetch_logs();
    	},
    	
    	fetch_logs: function() {
    		var self = this;
			coge.services.fetch_logs(self.id, self.type)
				.done(function(result) { // success
					//console.log(result);
					if (result && result.entries) {
						result.entries.forEach(function(entry) {
							self.element.append(
								'<tr>' + 
								'<td>' + entry.time + '</td>' +
								'<td>' + '<img width="18" height="18" src="image.pl?id=' + entry.user.image_id + '"/>' + entry.user.name + '</td>' +
								'<td>' + entry.description + '</td>' +
								'</tr>'
							);
						});
					}
				})
				.fail(function() { // error
					//TODO
				});
		}
    };

    return ns;
})(coge || {});
