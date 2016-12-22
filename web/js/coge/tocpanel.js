/*
 * tocpanel.js
 *
 *
 *
 * Requires: jQuery
 *
 */

class TocPanel {
    constructor(params) {
        this.element = $('#'+params.elementId);
        this.selection = params.selection;
        this.initialize();
    }

	initialize() {
		var self = this;

		// Style TOC
		this.element.addClass('coge-side-menu');

		// Add click handler
		this.element.find('span').each(function(index, value) {
			$(value).on('click', function () {
				self.clearSelection().selectItem(this);
		    });
		});

		//TODO dynamically generate html from types here instead of statically in .tmpl

		this.element.show();
    }

    setCount(itemType, count) {
        var item = this._getItem(itemType);
        var span = $('<span></span>').html(count).addClass('coge-item-count').css('display', 'inline');
        item.append(span);
    }

    clearSelection() {
    	this.element.find('span').removeClass('selected');
    	return this;
    }

    selectItemType(itemType) {
    	var item = this._getItem(itemType);
    	this.selectItem(item);
    }

    selectItem(item) {
    	this.element.find('span').removeClass('selected'); // disable all
    	$(item).addClass('selected'); // enable this one
    	var itemType = $(item).data('type');
    	console.log('TocPanel.selectItem ' + itemType);

    	if (this.selectedTypeId && itemType === this.selectedTypeId) // already selected
    		return;
    	this.selectedTypeId = itemType;

    	// Call user handler
    	if (this.selection)
			this.selection(itemType);

    	return this;
    }

    _getItem(itemType) {
        return $('span[data-type="'+itemType+'"]');
    }
}