/*
 * infopanel.js
 *
 *
 *
 * Requires: jQuery
 *
 */

class InfoPanel {
    constructor(params) {
        var element = $('#'+params.elementId);
        if (!params.elementID && !element) {
            console.error('InfoPanel: Missing required target element');
            return;
        }
        this.element = element;
        this.defaultInfo = params.defaultInfo;
	}

    busy() {
    	this.element.html('<img src="picts/ajax-loader.gif"/>');
    	return this;
    }

    update(items) {
    	console.log('InfoPanel.update');
    	var self = this;

    	if (!items || !items.length) {
    		if (self.defaultInfo)
    		    self.element.html( self.defaultInfo.call(self) );
    		return;
    	}

    	var numItems = items.length;
    	if (numItems > 0) {
    		if (numItems == 1) {
    			var item = items[0];
    			if (item) {
	    			item.getInfo().pipe(function(info) {
	    				self.element.html(info);
	    			});
    			}
    		}
    		else {
    			self.element.html(numItems + ' item' + (numItems > 1 ? 's' : '') +
    				' selected.<br><br>Click an action icon at the top to share, organize, delete, or analyze.');
    				//TODO add action links for sharing, adding to notebook, deleting, etc...
    		}
    	}
    }

    scheduleUpdate(items) {
    	if (this.timer)
    		window.clearTimeout(this.timer);

    	this.timer = window.setTimeout(
    		function() {
    			infoPanel.busy().update(items);
    		},
    		500
    	);
    }
}