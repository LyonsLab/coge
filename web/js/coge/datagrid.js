/* 
 * datagrid.js
 *
 * A wrapper around DataTables that provides CoGe-specific features.
 *
 * Requires: DataTables, jQuery
 *
 */

class DataGrid {
    constructor(params) {
        if (params.element)
            this.element = params.element;
        else if (params.elementId)
            this.element = $('#'+params.elementId);
        else
            console.warn('DataGrid: please specify target element');

        this.filter       = params.filter;
        this.selectable   = (typeof params.selectable !== 'undefined' ? params.selectable : true);
        this.selectionCallback = params.selectionCallback;
        this.mouseOver    = params.mouseOver;
        this.mouseOut     = params.mouseOut;
        this.dateSortAsc  = params.dateSortAsc;
        this.dateSortDesc = params.dateSortDesc;
        this.height       = params.height;
        this.class        = (typeof params.class !== 'undefined' ? params.class : "dt-cell hover compact row-border");
        this.columns      = params.columns;
        this.options      = params.options || {};
        if (!this.columns) {
            console.error('Missing required "columns" parameter');
            return;
        }

        this.initialize();
    }

	initialize() {
		var self = this;
		this.element.html('<img id="busy" src="picts/ajax-loader.gif" /><table cellpadding="0" cellspacing="0" border="0" class="' + self.class + '" style="cursor:pointer;"></table>');

		// Instantiate grid
		var dataTable = this.dataTable = this.element.children('table').dataTable($.extend(true, {}, {
			paging:    false,
			info:      false,
			searching: true,
			dom:       'lrt', // remove unused elements (like search box)
			scrollY:   self.height,
			language: { zeroRecords: "", emptyTable: "" },
			columns:   self.columns
		}, self.options));

		var dataTableBody = dataTable.children('tbody');

		// Handle row selection event
		if (self.selectable) {
            dataTableBody.on('click', 'tr', function(event) {
                var tr = this;
                var row = dataTable.api().row(tr).data();
                if (!row)
                    return;

                if ( $(tr).hasClass('selected') ) { // unselect
                    $(tr).removeClass('selected');
                }
                else { // select
                    if (event.ctrlKey || event.metaKey) // multi select
                        ; // no action required
                    else if (event.shiftKey) { // block select
                        var oSettings = dataTable.fnSettings(),
                            fromPos = dataTable.fnGetPosition(self.lastRowSelected),
                            toPos = dataTable.fnGetPosition(tr),
                            fromIndex = $.inArray(fromPos, oSettings.aiDisplay),
                            toIndex = $.inArray(toPos, oSettings.aiDisplay),
                            start = (fromIndex < toIndex ? fromIndex : toIndex),
                            end = (fromIndex < toIndex ? toIndex : fromIndex);

                        for (var i = start; i <= end; i++) {
                            var tr2 = dataTable.api().row(oSettings.aiDisplay[i]).node();
                            $(tr2).addClass('selected'); // select item
                        }
                    }
                    else
                        dataTable.$('tr.selected').removeClass('selected'); // unselect all

                    $(tr).addClass('selected'); // select item
                    self.lastRowSelected = tr;
                }

                self.selectItem(row);
            });

            // Handle row double-click event
            dataTableBody.on('dblclick', 'tr', function() {
                var tr = this;
                var row = dataTable.api().row(tr).data();

                if (row) {
                    self.dataTable.$('tr.selected').removeClass('selected'); // unselect all
                    $(tr).addClass('selected'); // select item
                    self.lastRowSelected = tr;
                    self.selectItem(row);
                    self.openItem(row);
                }
            });

            // Handle row hover events
            dataTableBody.on('mouseover', 'tr', function () {
                if (self.getSelectedItems()) // Do nothing if row(s) currently selected
                    return;

                var tr = $(this).closest('tr');
                var row = dataTable.api().row(tr).data();
                if (row && self.mouseOver)
                    self.mouseOver.call(self, row);
            });

            dataTableBody.on('mouseout', 'tr', function () {
                if (self.getSelectedItems()) // Do nothing if row(s) currently selected
                    return;

                if (self.mouseOut)
                    self.mouseOut.call(self);
            });
        }

		// Add custom filter
		if (self.filter) {
            $.fn.dataTable.ext.search.push(
                function(settings, data, dataIndex) {
                    var data = self.dataTable.api().row(dataIndex).data();
                    return self.filter(data);
                }
            );
		}

		// Add custom date sort functions
		if (self.dateSortAsc)
            $.fn.dataTable.ext.oSort['relative-date-asc']  = self.dateSortAsc.bind(this);
        if (self.dateSortDesc)
            $.fn.dataTable.ext.oSort['relative-date-desc'] = self.dateSortDesc.bind(this);
    }

    reset() {
    	// TODO
    	return this;
    }

    setData(data) {
        var busy = document.getElementById('busy');
        if (busy)
            busy.parentNode.removeChild(busy);

    	if (data) {
	    	this.dataTable.api()
				.clear()
				.rows.add(data);
    	}
    	else {
    	    this.dataTable.api()
    	        .clear()
    	}
    }

    update(data) {
        var busy = document.getElementById('busy');
        if (busy)
            busy.parentNode.removeChild(busy);

    	if (data) {
    	    this.setData(data);
    	    this.redraw();
    	}

        return this;
    }

    search(search_term) {
		this.dataTable.api()
			.search(search_term)
			.draw();
    }

    redraw() {
    	this.dataTable.api().draw();
    }

    getNumRows() {
    	return this.dataTable.api().page.info().recordsTotal;
    }

    getNumRowsDisplayed() {
    	return this.dataTable.api().page.info().recordsDisplay;
    }

    getSelectedRows() {
    	var rows = this.dataTable.api().rows('.selected');
    	return rows;
    }

    getSelectedItems() {
    	//console.log('getSelectedItems');
    	var items = this.dataTable.api().rows('.selected').data();
    	if (!items || !items.length)
    		return;
    	return items;
    }

    getSelectedItemList() {
    	var items = this.getSelectedItems();
    	var item_list;
    	if (items && items.length)
    		item_list = $.map(items, function(item) {
					return item.id + '_' + item.type;
				}).join(',');
    	return item_list;
    }

    setSelectedItems(items) {
    	this.dataTable.api().rows().every( function () {
    		var row = this;
    	    var d = row.data();
    	    items.each(function(item) {
	    	    if (d.id == item.id) {
	    	    	var tr = row.node();
	    	    	$(tr).addClass('selected'); // select item
	    	    }
    	    });
    	});
    }

    clearSelection() {
    	this.dataTable.$('tr.selected').removeClass('selected'); // unselect all
    }

    selectItem(item) {
    	console.log('DataGrid.selectItem');
    	var selectedItems = this.getSelectedItems();
    	if (this.selectionCallback)
    		this.selectionCallback.call(this, selectedItems);
    }

    openItem(row) {
        console.log('DataGrid.openItem');
        if (row && row.open)
            row.open.call(row);
    }
}
