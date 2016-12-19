/*
 * contentpanel.js
 *
 *
 *
 * Requires: jQuery
 *
 */

class ContentPanel {
    constructor(params) {
        this.element = $('#'+params.elementId);
        this.views = params.views;
        this.grid = params.grid;
        if (!this.element || !this.views || !this.grid) {
            console.error("Missing required element/view/grid inputs");
            return;
        }

        this.cache = new Array();
        this.selectedView = null;
    }

    getRow(dataTypeId, id) {
    	var row = null;
    	this.grid.dataTable.api().rows().every( function () {
    	    var d = this.data();
    	    if (d.id == id) {
    	    	row = this;
    	    }
    	});
    	return row;
    }

    getRowData(dataTypeId, id) {
    	var data = this.getData(dataTypeId);
    	var rowData = null;
    	if (data) {
    		data.some(function(d) {
    			if (d.id == id) {
    				rowData = d;
    				return true;
    			}
    			return false;
    		});
    	}
    	return rowData;
    }

    setRowData(dataTypeId, id, newData) {
    	this.grid.dataTable.api().rows().every( function () {
    	    var d = this.data();
    	    if (d.id == id) {
	    	    for (key in newData) {
	    	    	d[key] = newData[key];
	    	    }
	    	    this.data(d);
    	    }
    	});
    }

    getData(dataTypeId) {
    	var self = this;
    	var cachedData = new Array();

    	if (dataTypeId instanceof Array)
    		dataTypeId.forEach(function(i) {
    			cachedData = cachedData.concat(self.cache[i]);
    		});
    	else
    		cachedData = self.cache[dataTypeId];

    	return cachedData;
    }

    setData(typeId, data) {
    	console.log("ContentPanel.setData " + typeId);
    	var typeDef = this.views[typeId];

    	if (typeDef.displayType == 'grid') {
			this.cache[typeId] = data.map(function(obj) {
				return new DataGridRow(obj, typeId);
			});
		}
		else {
			this.cache[typeId] = data;
		}

    	return this;
    }

    setView(viewId) {
    	console.log('ContentPanel.update: ' + viewId + ' ');
    	this.selectedView = viewId;
    	return this;
    }

    getView() {
        return this.selectedView;
    }

    render(refresh) {
    	console.log('ContentPanel.render ' + this.selectedView);
    	if (!this.selectedView)
    		return;

        var view = this.views[this.selectedView];
        var isGrid = (view.displayType == 'grid');

        // Disable search bar if specified
        if (view.hasOwnProperty('search') && !view.search)
        	$('#search_input').hide();
        else
        	$('#search_input').show();

        // Disable flags column if specified
        if (view.hasOwnProperty('flagsColumn') && !view.flagsColumn)
        	this.grid.dataTable.api().column(0).visible(false);
        else
        	this.grid.dataTable.api().column(0).visible(true);

    	// Render contents
    	if (isGrid) {
    		// Save selection and scroll position
    		var items = this.grid.getSelectedItems();
    		var scrollPos = this.element.find(".dataTables_scrollBody").scrollTop();

    		// Swap in grid and update contents
    		this.element.children('.html').hide();
    		this.element.children('.grid').show();
    		this.grid.update(this.getData(view.dataTypes));
    		if (!refresh)
    			this.grid.dataTable.api().order(view.defaultSort ? view.defaultSort : [1, 'asc']);
    		this.grid.redraw(); // needed to display column widths properly

    		// Restore selection and scroll position
    		if (items)
    			this.grid.setSelectedItems(items);
    		this.element.find(".dataTables_scrollBody").scrollTop(scrollPos);
    	}
    	else {
    		this.element.children('.grid').hide();
    		this.element.children('.html').html(this.getData(this.selectedView)).show();
    	}

        // Update title with row number
        this.renderTitle();

        // Show/hide action icons based on type of data
    	$('.item-button').hide(); // hide all icons
    	if (view.operations) {
    		view.operations.forEach(function(op) {
    			$('.'+op).show();
    		});
    	}

    	// Icons are initially set to invisible on load to prevent flickering
    	$('.item-button').removeClass('invisible');

    	// Update browser url
    	var params = coge.utils.getURLParameters();
    	params['p'] = this.selectedView;
    	var queryString = coge.utils.toQueryString(params);
    	window.history.pushState({}, "", PAGE_NAME + "?" + queryString);
    }

    renderTitle() {
    	var view = this.views[this.selectedView];
    	var title = view.title;
    	var isGrid = (view.displayType == 'grid');
        if (isGrid)
        	title += '&nbsp;&nbsp;<span class="small info">' + this.grid.getNumRowsDisplayed() + '</span>';
        $('#contents_title').html(title);
    }

    busy() {
    	var spinner = '<div class="spinner" style="display:flex;justify-content:center;align-items:center;margin-top:40%;"></div>';
    	this.element.children('.grid').hide();
    	this.element.children('.html').html(spinner).show();
    }
}