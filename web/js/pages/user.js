const POLL_TIME = 30*1000, // polling rate when user is not idle
	  IDLE_TIME = 30*1000; // stop polling after this lapse, then poll on next mousemove

var grid;
var infoPanel;
var tocPanel;
var timestamps = new Array();
var timers = new Array();

$(function() {
	// Initialize AJAX
	$.ajaxSetup({
		type: "GET",
		url: PAGE_NAME,
		dataType: "html",
		cache: false,
	});
	
	// Initialize dialog boxes
	$(".dialog_box").dialog({autoOpen: false, resizable: false});

	// Initialize dropdown menus
	$("#create_menu").menu()
		.position({
			my: "left top",
			at: "left bottom",
			of: "#create_button"
		})
		.css({ position: 'absolute' });
	
	$("#send_menu").menu().css({ position: 'absolute', width: '90px' });

	// Define views in the Content Panel
	var views = {
		mine: {
			title: 'My Data',
			displayType: 'grid',
			dataTypes: ['genome', 'experiment', 'notebook'],
			operations: ['share', 'organize', 'favorite', 'delete', 'sendto']
		},
		genome: {
			title: 'Genomes',
			displayType: 'grid',
			dataTypes: ['genome'],
			operations: ['share', 'organize', 'favorite', 'delete', 'sendto']
		},
		experiment: {
			title: 'Experiments',
			displayType: 'grid',
			dataTypes: ['experiment'],
			operations: ['share', 'organize', 'favorite', 'delete', 'sendto']
		},
		notebook: {
			title: 'Notebooks',
			displayType: 'grid',
			dataTypes: ['notebook'],
			operations: ['share', 'favorite', 'delete', 'sendto', 'add']
		},
		shared: {
			title: 'Shared with me',
			displayType: 'grid',
			dataTypes: ['genome', 'experiment', 'notebook'],
			operations: ['share', 'organize', 'favorite'],
			shared: true
		},
		favorite: {
			title: 'Favorites',
			displayType: 'grid',
			dataTypes: ['genome', 'experiment', 'notebook', 'favorite'],
			operations: ['share', 'organize', 'favorite', 'sendto'],
			favorite: true
		},	
		metadata: {
			title: 'Metadata',
			displayType: 'html',
			dataTypes: ['metadata'],
			search: false
		},
		group: {
			title: 'User Groups',
			displayType: 'grid',
			dataTypes: ['group'],
			operations: ['edit', 'delete', 'add'],
			shared: true,
			flagsColumn: false
		},
		activity: {
			title: 'Activity',
			displayType: 'html',
			dataTypes: ['activity'],
			search: false
		},
		analyses: {
			title: 'Analyses',
			displayType: 'grid',
			dataTypes: ['analyses'],
			noFilter: true,
			defaultSort: [ 2, 'desc' ],
			flagsColumn: false
		},
		loads: {
			title: 'Data loading',
			displayType: 'grid',
			dataTypes: ['loads'],
			noFilter: true,
			defaultSort: [ 2, 'desc' ],
			flagsColumn: false
		},
		graph: {
			title: 'Graph',
			displayType: 'html',
			dataTypes: ['graph'],
			search: false,
			refresh: false
		},
		trash: {
			title: 'Trash',
			displayType: 'grid',
			dataTypes: ['genome', 'experiment', 'notebook', 'group'],
			operations: ['undelete'],
			deleted: true
		}
	};
	
	// Initialize the main panels
	infoPanel = new InfoPanel({
		elementId: 'info_panel'
	});
	
	contentPanel = new ContentPanel({
		elementId: 'contents_panel',
		views: views
	});
	
	tocPanel = new TocPanel({
		elementId: 'toc_panel',
		selection: function(typeId) {
			cancel_poll();
			contentPanel
				.update(typeId)
				.done(function() { 
					contentPanel.grid.search(''); // clear search filter
					contentPanel.render();
					schedule_poll();
				});
			infoPanel.update(null);
			update_icons(null);
			$('#search_input').val(''); //FIXME move into ContentPanel
		}
	});
	
	$('#search_input').on('keyup search', function() { //FIXME move into ContentPanel
		contentPanel.grid.search( $(this).val() );
		contentPanel.renderTitle();
	});
	
	// Get starting page from URL and initialize TOC panel
	var toc_id = getURLParameter('p');
	if (!toc_id || toc_id == 'null')
		toc_id = 'mine';
	tocPanel.selectItemType(toc_id);
	
	// Setup idle timer
	init_timestamp('idle');
	$(document).mousemove(function() {
		var currentTime = new Date().getTime();
		var idleTime = currentTime - timestamps['idle'];
		timestamps['idle'] = currentTime;

		if (idleTime > IDLE_TIME) {
			// User was idle for a while, refresh page immediately
			schedule_poll(0);
		}
	});

	// Initiate refresh loop
	schedule_poll();

	// Initialize confirm cancel job dialog
    $("#cancel_dialog").dialog({
    	modal: true,
    	buttons: {
    		No: function() {
    			$(this).dialog("close");
    		},
    		Yes: function() {
    			cancel_job( $(this).data("log_id") );
    			$(this).dialog("close");
    		}
    	}
    });
    
	// Initialize comment dialog
    $("#comment_dialog").dialog({
    	width: 400,
    	modal: true,
    	buttons: {
    		OK: function() {
    			var log_id = $(this).data("log_id");
    			var comment = $(this).find("input").first().val();
    			contentPanel.setRowData('analyses', log_id, {comment: comment});
    			contentPanel.grid.redraw();
    			comment_job( log_id, comment );
    			$(this).dialog("close");
    		},
    		Cancel: function() {
    			$(this).dialog("close");
    		}
    	}
    });
});

function pad(string, size) {
    while (string.length < size) string = "0" + string;
    return string;
}

function getURLParameter(name) {
    return decodeURI(
        (RegExp(name + '=' + '(.+?)(&|$)').exec(location.search)||[,null])[1]
    );
}

function poll(sync) {
	console.log('poll');
	
	// Refresh contents
	//get_contents(sync, pageObj.content_type);
	contentPanel
		.refresh()
		.done(function() { 
			contentPanel.render(true);
			schedule_poll(); // schedule next poll
		});
	
	// Refresh item info cache
	//grid.reset();
}

function schedule_poll(when) { 
	cancel_poll();

	if (when !== undefined) {
		timers['poll'] = window.setTimeout(
			function() { poll(1); },
			when
		);
		return;
	}

	// Quit polling if page idle for too long
	var idleTime = new Date().getTime() - timestamps['idle'];
	if (idleTime < IDLE_TIME) {
		timers['poll'] = window.setTimeout(
			function() { poll(1); },
			POLL_TIME
		);
	}
}

function cancel_poll() {
	clearTimeout(timers['poll']);
}

function default_info() {
	switch(contentPanel.selectedView) {
		case 'activity':
			return "Here is a summary of all analyses you have performed.";
		case 'analyses':
			return "These are the analyses you have performed or started.<br><br>" + 
				"Select an analysis to open the current progress or finished result in a new tab.<br><br>" +
				"Use the icons to the left of each analysis to 'Favorite' it, add comments, or cancel (if running).";
		case 'loads':
			return "These are the data loading workflows you have performed or started.<br><br>" +
				"Select an item to open the current progress or finished result in a new tab.";
		case 'graph':
			return "This is a graphical representation of the analyses you've run.<br><br>Click a node to see the individual analyses of that type.";
		case 'trash':
			return "These are items you deleted.<br><br>" +
				"Hover over an item to view additional info. Select one or more items to undelete.";
		case 'shared':
			return "These are data items that your collaborators shared with you.<br><br>" +
				"Hover over an item to view additional info. Select one or more items to share with others or add to a notebook.";
		case 'favorite':
			return "These are data items that you have marked as your favorites.<br><br>" +
				   "Hover over an item to view additional info.";
		case 'metadata':
			return "Here is a summary of the metadata for your experiments, genomes and notebooks.";
		case 'mine':
		case 'notebook':
		case 'genome':
		case 'experiment':
			return "<p>These are data items that you added to the system.</p>"
				+ "<p><b>Hover over</b> an item to view additional info.</p>"
                + "<p><b>Single-click</b> to select one or more items to share, organize, delete, or send them to one of CoGe's tools. Use <b>Ctrl-click</b> to select multiple items.</p>"
                + "<p><b>Double-click</b> an item for a detailed view of the item.</p>";
		case 'group':
			return "You are a member of these user groups.<br><br>" +
				"Hover over a group to view additional info. Select one or more groups to edit or delete."; 
	}
}

/*
 * Content Panel
 */
function ContentPanel(params) {
	this.element = $('#'+params.elementId);
	this.views = params.views;
	this.cache = new Array();
	this.selectedView = null;
	this.initialize();
}

$.extend(ContentPanel.prototype, {
	initialize: function() {
		var self = this;
		
		// Create grid
		self.grid = new DataGrid({
			element: self.element.children('.grid'),
			filter: function(data) { 
				// Filter rows based on view
				var view = self.views[self.selectedView];
				if (!view.noFilter) {
					if (view.deleted && data.deleted == '0')
						return false;
					if (!view.deleted && data.deleted == '1')
						return false;
					if (view.shared && data.role_id == '2')
						return false;
					if (!view.shared && data.role_id != '2' && !view.deleted)
						return false;
					if (view.favorite && data.favorite == '0')
						return false;
				}
				return true;
			},
			selection: function(items) {
				// Update icons
				update_icons(items);
			}
		});
	},
    
    getRow: function(dataTypeId, id) {
    	var row = null;
    	this.grid.dataTable.api().rows().every( function () {
    	    var d = this.data();
    	    if (d.id == id) {
    	    	row = this;
    	    }
    	});
    	return row;
    },
    
    getRowData: function(dataTypeId, id) {
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
    },
    
    setRowData: function(dataTypeId, id, newData) {
    	this.grid.dataTable.api().rows().every( function () {
    	    var d = this.data();
    	    if (d.id == id) {
	    	    for (key in newData) {
	    	    	d[key] = newData[key];
	    	    }
	    	    this.data(d);
    	    }
    	});
    },    
    
    getData: function(dataTypeId) {
    	var self = this;
    	var cachedData = new Array();
    	
    	if (dataTypeId instanceof Array)
    		dataTypeId.forEach(function(i) {
    			cachedData = cachedData.concat(self.cache[i]);
    		});
    	else
    		cachedData = self.cache[dataTypeId];
    	
    	return cachedData;
    },
    
    setData: function(typeId, data) {
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
    },
    
    render: function(refresh) {
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
    },
    
    renderTitle: function() {
    	var view = this.views[this.selectedView];
    	var title = view.title;
    	var isGrid = (view.displayType == 'grid');
        if (isGrid)
        	title += '&nbsp;&nbsp;<span class="small info">' + this.grid.getNumRowsDisplayed() + '</span>';
        $('#contents_title').html(title);
    },
    
    busy: function() {
    	var spinner = '<div class="spinner" style="display:flex;justify-content:center;align-items:center;margin-top:40%;"></div>';
    	this.element.children('.grid').hide();
    	this.element.children('.html').html(spinner).show();
    },
    
    update: function(viewId) {
    	console.log('ContentPanel.update: ' + viewId + ' ');
    	var self = this;
    	this.selectedView = viewId;
    	var view = this.views[viewId];
    	
    	// Setup requests for content data
    	var promises = new Array();
    	view.dataTypes.forEach(function (dataType) {
    		var deferred = $.Deferred();
    		var cachedData = self.getData(dataType);
    	
	        if (cachedData) {
		    	deferred.resolve();
	    	}
	    	else {
	    		self.busy();
	    		deferred = self.fetch(false, dataType);
		    	setTimeout(deferred.resolve, 10);
	    	}
	        promises.push(deferred);
    	});
        
        return $.when.apply($, promises).then(function(schemas) {
	            console.log("ContentPanel.update: DONE");
	            self.grid.clearSelection();
	        }, function(e) {
	            console.log("ContentPanel.update: FAILED");
	        });
    },
    
    refresh: function() {
    	console.log('ContentPanel.refresh');
    	var self = this;
    	if (!this.selectedView)
    		return;
    	
    	var view = this.views[this.selectedView];
    	
    	// Skip refresh if specified
    	if (view.hasOwnProperty('refresh') && !view.refresh)
    		return;
    	
    	$('#refresh_label').fadeIn(); //FIXME move into ContentPanel
    	
       	var promises = new Array();
    	view.dataTypes.forEach(function (dataType) {
    		console.log('refresh ' + dataType);
    		var deferred = self.fetch(false, dataType);
	    	setTimeout(deferred.resolve, 10);
	        promises.push(deferred);
    	});
        
        return $.when.apply($, promises).then(function(schemas) {
	            console.log("ContentPanel.refresh: DONE");
	            $('#refresh_label').fadeOut(); //FIXME move into ContentPanel
	        }, function(e) {
	            console.log("ContentPanel.refresh: ajax failed");
	        });
    },    
    
    fetch: function(sync, typeId) {
    	var self = this;
    	console.log('ContentPanel.fetch ' + typeId + ' ' + self.selectedTypeId);
    	if (!typeId)
    		typeId = self.selectedTypeId;
    	
    	var lastUpdate = (sync ? timestamps['lastUpdate'] : 0);

    	return $.ajax({
    		dataType: 'text',
    		data: {
    			fname: 'get_contents',
    			user_id: USER_ID,
    			item_type: typeId,
    			last_update: lastUpdate,
    			timestamp: init_timestamp('get_contents')
    		},
    		success : function(data) {
    			if (!data) {
    				console.warn('get_contents: null data');
    				return;
    			}
    			//console.log(data);
    			if (self.views[typeId].displayType == 'grid') {
    				data = JSON.parse(data);
    			}
    			self.setData(typeId, data);
    		},
    		complete : function() {
    			
    		}
    	});
    }
});

/*
 * Data Grid
 */

function DataGrid(params) {
	if (params.element)
		this.element = params.element;
	else if (params.elementId)
		this.element = $('#'+params.elementId);
	else 
		console.warn('DataGrid: please specify target element');
	
	this.filter = params.filter;
	this.selection = params.selection;
	
	this.initialize();
}

$.extend(DataGrid.prototype, {
	initialize: function() {
		var self = this;
		this.element.html('<table cellpadding="0" cellspacing="0" border="0" class="dt-cell hover compact row-border" style="cursor:pointer;"></table>');
		
		// Instantiate grid
		var dataTable = this.dataTable = this.element.children('table').dataTable({
			paging:    false,
			info:      false,
			searching: true,
			dom:       'lrt', // remove unused elements (like search box)
			sScrollY:  $(window).height() - 210, // this depends on the height of the header/footer and should be passed in as an argument
			oLanguage: { "sZeroRecords": "", "sEmptyTable": "" },
			columns: [
	            { 	title: "", 
	            	targets: 0,
	            	orderable: false,
	            	type: "string",
	            	width: 15,
	            	data: null, // use full data object
	            	render: function(data, type, row, meta) {
	            		return data.getFlags();
	            	}
	            },			          
	            { 	title: "Name", 
	            	targets: 1,
	            	type: "html",
	            	data: null, // use full data object
	            	render: function(data, type, row, meta) {
	            		return data.getDescription();
	            	}
	            },
	            { 	title: "Date added", 
	            	targets: 2, 
	            	type: "relative-date", // this is our own custom type
	            	data: null, // use full data object
					orderSequence: [ 'desc', 'asc' ],
	            	width: "100px",
	            	render: function(data, type, row, meta) {
	            		return data.getDate();
	            	}
	            }
			]
		});
		
		var dataTableBody = dataTable.children('tbody');
		
		// Handle row selection event
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
	        if (row)
	        	infoPanel.busy().scheduleUpdate([row]);
	    });
		
		dataTableBody.on('mouseout', 'tr', function () {
	    	if (self.getSelectedItems()) // Do nothing if row(s) currently selected
	    		return;
	    	
	    	infoPanel.scheduleUpdate();
	    });
		
		// Add custom filter
		$.fn.dataTable.ext.search.push(
			function(settings, data, dataIndex) { 
				var data = self.dataTable.api().row(dataIndex).data();
				return self.filter(data); 
			}
		);
		
		// Add custom date sort functions
		$.fn.dataTable.ext.oSort['relative-date-asc'] = function(x,y) {
			x = this.time_diff(x);
			y = this.time_diff(y);
		    return ((x < y) ? 1 : ((x > y) ? -1 : 0));
		}.bind(this);
		$.fn.dataTable.ext.oSort['relative-date-desc'] = function(x,y) {
			x = this.time_diff(x);
			y = this.time_diff(y);
			return ((x < y) ? -1 : ((x > y) ? 1 : 0));
		}.bind(this);
    },

	time_diff: function(s) {
		var matches = s.match(/^\<\!\-\-(\d+)\-\-\>/); // the time difference is hidden in an html comment
		if (!matches || !matches[1])
			return 999999; // sort last
		var diff = matches[1];
		return parseInt(diff);
	},
    
    reset: function() {
    	// TODO
    	return this;
    },
    
    update: function(data) {
    	console.log('DataGrid.update');
    	
    	if (data) {
	    	this.dataTable.api()
				.clear()
				.rows.add(data)
				.draw();
    	}
		
        return this;
    },
    
    search: function(search_term) {
		this.dataTable.api()
			.search(search_term)
			.draw();
    },
    
    redraw: function() {
    	this.dataTable.api().draw();
    },
    
    getNumRows: function() {
    	return this.dataTable.api().page.info().recordsTotal;
    },    
    
    getNumRowsDisplayed: function() {
    	return this.dataTable.api().page.info().recordsDisplay;
    },
    
    getSelectedRows: function() {
    	var rows = this.dataTable.api().rows('.selected');
    	return rows;
    },
    
    getSelectedItems: function() {
    	//console.log('getSelectedItems');
    	var items = this.dataTable.api().rows('.selected').data();
    	if (!items || !items.length)
    		return;
    	return items;
    },
    
    getSelectedItemList: function() {
    	var items = this.getSelectedItems();
    	var item_list;
    	if (items && items.length)
    		item_list = $.map(items, function(item) {
					return item.id + '_' + item.type;
				}).join(',');
    	return item_list;
    },
    
    setSelectedItems: function(items) {
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
    },
    
    clearSelection: function() {
    	this.dataTable.$('tr.selected').removeClass('selected'); // unselect all
    },
    
    selectItem: function(item) {
    	console.log('DataGrid.selectItem');
    	
    	var selectedItems = this.getSelectedItems();
    	infoPanel.busy().update(selectedItems); //FIXME move into selection handler
    	
    	if (this.selection)
    		this.selection(selectedItems);
    },

    openItem: function(row) {
    	if (row.type == 'group')
    		group_dialog();
    	else if (row.type == 'analyses' || row.type == 'loads')
    		window.open(row.link, '_blank');
    	else {
    		var title = row.getDescription();
    		var link = row.getLink();
    		var flags = row.getFlags({noSpaces: 1});
    		title = flags + ' ' + title + "<br><a class='xsmall' style='color:#eeeeee;' href='" + link + "' target='_blank'>[Open in new tab]</a> ";
    		link = link + "&embed=1";
    		console.log('DataGrid.openItem: ' + link);
    		var height = $(window).height() * 0.8;
    		var d = $('<div class="dialog_box"><iframe src="'+link+'" height="100%" width="100%" style="border:none;"/></div>')
    			.dialog({
    				//title: title,
    				width: '80%',
    				height: height,
                    open: function() { // mdb added 10/16/16 -- fix html in dialog title bar for jQuery 3.1.1 update
                        $(this).prev().find("span.ui-dialog-title").append('<span>'+title+'</span>');
                    }
    			})
    			.dialog('open');
    	}
    }
});

/* 
 * Data Grid Row
 */
		
function DataGridRow(data, type) {
	$.extend(this, data);
	if (!this.type) // mdb added condition 9/16/18 COGE-388 -- prevent native type (genome,etc) from being set to "favorite"
		this.type = type;
}

$.extend(DataGridRow.prototype, { // TODO consider extending this into separate classes for each type (genome, experiment, etc...)
    getFlags: function(opts) {
    	if (this.type == 'genome' || 
    		this.type == 'experiment' || 
    		this.type == 'notebook' || 
    		this.type == 'favorite') 
    	{
    		var noSpaces = (opts && opts.noSpaces);
    		
    		var flags = '';
    		if (!noSpaces || this.favorite == '1')
    			flags = '<span style="color:goldenrod;visibility:' + 
	    			(this.favorite == '1' ? 'visible' : 'hidden') + 
	    			'">&#9733;</span>&nbsp;';
    		if (this.restricted == '1')
	    		flags += '&#x1f512;' + '&nbsp;';
    		
    		return flags;
    	}
    	return '';
    },
    
    getDescription: function() {
    	if (this.type == 'genome' || this.type == 'favorite')
    		return this._formatGenome();
    	if (this.type == 'experiment')
    		return this._formatExperiment();
    	if (this.type == 'notebook')
    		return this._formatNotebook();
    	if (this.type == 'group')
    		return this._formatGroup();
    	if (this.type == 'analyses')
    		return this._formatAnalysis();
    	if (this.type == 'loads')
    		return this._formatLoad();
    },
    
    _formatGenome: function() {
    	var icon = '<img src="picts/dna-icon.png" width="15" height="15" style="vertical-align:middle;"/> ';
    	var certified = '<span class="glyphicon glyphicon-ok coge-certified-icon"></span> <span class="coge-small-text">Certified Genome<span>';
    	var descStr = 
    		icon +
    	   	(this.organism ? this.organism : '') + 
    	   	(this.name ? ' (' + this.name + ')' : '') +
    	   	(this.description ? ': ' + this.description : '') +
    	   	' (v' + this.version + ', id' + this.id + ')' +
    	   	(this.certified == '1' ? '&nbsp;&nbsp;' + certified : '');
    	return descStr;
    },
    
    _formatExperiment: function() {
    	var descStr = 
    		'<img src="picts/testtube-icon.png" width="15" height="15" style="vertical-align:middle;"/> ' +
    	   	this.name +
    	   	(this.description ? ': ' + this.description : '') +
    	   	' (v' + this.version + ', id' + this.id + ')';
    	return descStr;
    },
    
    _formatNotebook: function() {
    	var descStr =
    		'<img src="picts/notebook-icon.png" width="15" height="15" style="vertical-align:middle;"/> ' +
    		this.name +
    		(this.description ? ': ' + this.description : '') +
    		(this.type_name ? ' (' + this.type_name + ')' : '');
    	return descStr;
    },
    
    _formatGroup: function() {
    	var descStr =
    		'<img src="picts/group-icon.png" width="15" height="15" style="vertical-align:middle;"/> ' +
    		this.name +
    		(this.description ? ': ' + this.description : '');;
    	return descStr;
    },
    
    _formatWorkflowStatus: function(status) {
    	status = status.toLowerCase();
        var color;
        
        if (status == 'terminated') status = 'cancelled';
        
        switch (status) {
        	case 'running':   color = 'yellowgreen'; 	break;
        	case 'completed': color = 'cornflowerblue'; break;
        	case 'scheduled': color = 'goldenrod'; 		break;
            default:          color = 'salmon';
        }
        
        return '<span style="padding-bottom:1px;padding-right:5px;padding-left:5px;border-radius:15px;color:white;background-color:' + color + ';">' + coge.utils.ucfirst(status) + '</span>';
    },
    
    _formatAnalysis: function() {
        var isRunning   = (this.status.toLowerCase() == 'running');
        var isCancelled = (this.status.toLowerCase() == 'cancelled');
        var star_icon    = '<img title="Favorite this analysis" src="picts/star-' + (this.is_important ? 'full' : 'hollow') + '.png" width="15" height="15" class="link" style="vertical-align:middle;" onclick="toggle_star(this, '+this.id+');" />';
        var cancel_icon  = '<img title="Cancel this analysis" class="link" height="15" style="vertical-align:middle;" src="picts/cancel.png" width="15" onclick="cancel_job_dialog('+this.id+');"/>';
        var restart_icon = '<img title="Restart this analysis" class="link" height="15" style="vertical-align:middle;" src="picts/refresh-icon.png" width="15" onclick="restart_job('+this.id+');"/>';
        var comment_icon = '<img title=' + (this.comment && this.comment.length ? JSON.stringify(this.comment) : '"Add comment"') + ' class="link" height="15" style="vertical-align:middle;" src="picts/' + (this.comment && this.comment.length ? 'comment' : 'no-comment') + '-icon.png" width="15" onclick="comment_dialog('+this.id+');" />';
        var icons = star_icon + ' ' + comment_icon + ' ' + (isCancelled ? restart_icon : '') + ' ' + (isRunning ? cancel_icon : '');
    	var descStr = icons + ' ' + this._formatWorkflowStatus(this.status) + ' ' + this.page + ' | ' + this.description + (this.comment ? ' | ' + '<span class="highlighted">' + this.comment : '') + '</span>' + ' | ' + this.elapsed + (this.workflow_id ? ' | id' + this.workflow_id : '');
    	return descStr;
    },
    
    _formatLoad: function() {
    	var descStr =
    		this._formatWorkflowStatus(this.status) + ' ' + this.page + ' | ' + this.description + ' | ' + this.elapsed + (this.workflow_id ? ' | id' + this.workflow_id : '');
    	return descStr;
    },

    getInfo: function() {
    	console.log('DataGridRow.getInfo');
    	var self = this;
    	
    	return $.ajax({
    		dataType: 'json',
    		data: {
    			fname: 'get_item_info',
    			item_id: self.id,
    			item_type: self.type,
    			timestamp: init_timestamp('get_item_info')
    		}
    	}).pipe(function(data) {
    		if (data && test_timestamp('get_item_info', data.timestamp))
				return data.html;
    		return;
    	});
    },
    
    getLink: function() {
    	var type = this.type;
    	
    	if (type == 'genome')
    		return 'GenomeInfo.pl?gid=' + this.id;
    	else if (type == 'experiment')
    		return 'ExperimentView.pl?eid=' + this.id;
    	else if (type == 'notebook')
    		return 'NotebookView.pl?nid=' + this.id;
    	else
    		return this.link;
    },
    
    getDate: function() {
    	var dateStr = this.date;
    	
    	// Handle null date
    	if (!dateStr || dateStr.indexOf('0000') == 0)
    		dateStr = this.dataset_date;
    	if (!dateStr || dateStr.indexOf('0000') == 0)
    		return '&nbsp;'; // return blank (needed for cell to render properly in Safari)
    	
    	// Convert from database time zone (Tucson) to user's time zone -- mdb added 9/22/16 COGE-196
    	dateStr = coge.utils.timeToLocal(dateStr);
    	
        // Convert date from aboslute into relative form ("yesterday", etc)
    	const MONTHS = [ 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec' ];
    	dateStr = dateStr.replace(/-/g, '/'); // needed for Firefox & Safari
    	var date = new Date(dateStr);
    	var today = new Date();
    	var diffMS = Math.abs(today.getTime() - date.getTime());
    	var diffDays = Math.floor(diffMS/(24*60*60*1000));
    	
    	if (diffDays == 0 && today.getDate() == date.getDate()) { // same day as today
    		var hour = (date.getHours() + 24) % 12 || 12;
    		dateStr = hour + ':' + pad(''+date.getMinutes(), 2) + ' ' + (date.getHours() < 12 ? 'am' : 'pm');
    	}
    	else if (diffDays <= 1) // yesterday
    		dateStr = 'Yesterday';
    	else if (diffDays <= 4) // last several days
    		dateStr = diffDays + ' days ago';
    	else if (date.getFullYear() == today.getFullYear()) // same year 
    		dateStr = MONTHS[date.getMonth()] + ' ' + date.getDate()
    	else // last year or older
    		dateStr = MONTHS[date.getMonth()] + ' ' + date.getDate() + ', ' + date.getFullYear();

    	return '<!--' + diffMS + '-->' + dateStr; // embed the time difference in a hidden html comment for sorting
    }
});

/* 
 * Info Panel
 */

function InfoPanel(params) {
	this.element = $('#'+params.elementId);
    this.initialize();
}

$.extend(InfoPanel.prototype, {
	initialize: function() {
    },
    
    busy: function() {
    	this.element.html('<img src="picts/ajax-loader.gif"/>');
    	return this;
    },
    
    update: function(items) {
    	console.log('InfoPanel.update');
    	var self = this;
    	
    	if (!items || !items.length) {
    		self.element.html( default_info() );
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
    },
    
    scheduleUpdate: function(items) {
    	if (this.timer)
    		window.clearTimeout(this.timer);

    	this.timer = window.setTimeout( 
    		function() { 
    			infoPanel.busy().update(items);
    		},
    		500
    	);
    }
});

/* 
 * Table of Contents Panel
 */

function TocPanel(params) {
	this.element = $('#'+params.elementId);
	this.selection = params.selection;
    this.initialize();
}

$.extend(TocPanel.prototype, {
	initialize: function() {
		var self = this;
		
		// Style TOC
		this.element.addClass('coge-side-menu');
		
		// Add click handler
		this.element.find('span').each(function(index, value) {
			$(value).on('click', function () {
				self.clearSelection().selectItem(this);
		    });
		});
		
		//TODO dynamically generate html from types here instead of statically in User.tmpl
		
		this.element.show();
    },
    
    clearSelection: function() {
    	this.element.find('span').removeClass('selected');
    	return this;
    },
    
    selectItemType: function(itemType) {
    	var item = $('span[data-type="'+itemType+'"]');
    	this.selectItem(item);
    },
    
    selectItem: function(item) {
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
});

function update_icons(items) { //TODO move into ContentPanel
	if ( items && items.length > 0) 
		$('.item-button:not(#add_button)').removeClass('coge-disabled');
	else
		$('.item-button:not(#add_button)').addClass('coge-disabled');
}

function favorite_items() {
	var selected_rows = contentPanel.grid.getSelectedRows();
	var item_list = contentPanel.grid.getSelectedItemList();
	if (item_list) {
		// Optimistically toggle favorite in UI
		selected_rows.every(function() {
			var d = this.data();
			d.favorite = (d.favorite == '0' ? '1' : '0');
			this.data(d);
		});
		contentPanel.grid.redraw();
		infoPanel.update(null);
		
		// Toggle favorite on server
		$.ajax({
			data: {
				fname: 'favorite_items',
				item_list: item_list
			},
			success : function(data) {
				// TODO handle error
			}
		});
	}
}

function delete_items() {
	var selected_rows = contentPanel.grid.getSelectedRows();
	var item_list = contentPanel.grid.getSelectedItemList();
	if (item_list) {
		$.ajax({
			data: {
				fname: 'delete_items',
				item_list: item_list
			},
			success : function(data) {
				selected_rows.every(function() {
					var d = this.data();
					d.deleted = '1';
					this.data(d);
				});
				contentPanel.grid.redraw();
				infoPanel.update(null);
			}
		});
	}
}

function undelete_items() {
	var selected_rows = contentPanel.grid.getSelectedRows();
	var item_list = contentPanel.grid.getSelectedItemList();
	if (item_list) {
		$.ajax({
			data: {
				fname: 'undelete_items',
				item_list: item_list
			},
			success : function(data) {
				selected_rows.every(function() {
					var d = this.data();
					d.deleted = '0';
					this.data(d);
				});
				contentPanel.grid.redraw();
				infoPanel.update(null);
			}
		});
	}
}

function cancel_job_dialog(id) {
	if (id) {
		$('#cancel_dialog')
			.data("log_id", id)
			.dialog('open');
	}
}

function cancel_job(id) {
	if (id) {
		var row = contentPanel.getRow('analyses', id);
		if (row) {
			var data = row.data();
			$.ajax({
				data: {
					fname: 'cancel_job',
					workflow_id: data.workflow_id
				},
				success : function(rc) {
					if (rc) {
						//poll(0);//schedule_poll(0); // FIXME mdb changed 10/7/14
						
						// Update status to cancelled in displayed row
						var data = row.data();
						data.status = 'Cancelled';
						row.data(data);
		            }
				}
			});
		}
	}
}

function restart_job(id) {
	if (id) {
		var row = contentPanel.getRow('analyses', id);
		var data = row.data();
		if (data.link) {
			// Change status to running in displayed row
			data.status = 'Running';
			row.data(data);
			
			// Open the link in a "hidden" window.
			// Workaround for restarting a workflow until JEX implements a 
			// restart command.
			var w = window.open(data.link,'_blank', 'toolbar=no,status=no,menubar=no,scrollbars=no,resizable=no,left=10000,top=10000,width=1,height=1,visible=none', ''); 
			setTimeout(function() {
					//poll(0);//schedule_poll(0); // FIXME mdb changed 10/7/14
					w.close();
				},
				5*1000
			);
		}
	}
}

function comment_dialog(id) {
	if (id) {
		var data = contentPanel.getRowData('analyses', id);
		$('#comment_dialog').find("input").first().val(data.comment);
		$('#comment_dialog')
			.data("log_id", id)
			.dialog('open');
	}
}

function comment_job(id, comment) {
	$.ajax({
		data: {
			fname: 'comment_job',
			log_id: id,
			comment: comment
		},
		success : function(rc) {
			if (rc) {
				//schedule_poll(0);
				contentPanel.grid.redraw();
            }
		}
	});
}

function add_to_notebook_dialog() {
	var selected = contentPanel.grid.getSelectedItems();
	if (selected.length) {
		// Initialize add-to-notebook dialog
		window.setTimeout(search_notebooks, 1000);
		
		// Open dialog window
		$('#add_to_notebook_dialog').dialog({width:500}).dialog('open');
	}
}

function init_timestamp(name) { // TODO move into class
	timestamps[name] = new Date().getTime()
	return timestamps[name];
}
function test_timestamp(name, time) {
	return time >= timestamps[name];
}

function wait_to_search (search_func, search_term) {
	if (!search_term || search_term.length > 2) {
		pageObj.search_term = search_term;
		if (timers['search']) {
			clearTimeout(timers['search']);
		}

		timers['search'] = window.setTimeout(
			function() {
				search_func(pageObj.search_term);
			},
			500
		);
	}
}

function search_notebooks () {
	var search_term = $('#notebook_search_input').val();

	$("#wait_notebook").animate({opacity:1});
	$("#notebook_select").html("<option disabled='disabled'>Searching...</option>");

	$.ajax({
		data: {
			fname: 'search_notebooks',
			search_term: search_term,
			timestamp: init_timestamp('search_notebooks')
		},
		success : function(data) {
			if (data) {
				var obj = jQuery.parseJSON(data);
				if (test_timestamp('search_notebooks', obj.timestamp)) {
					$("#notebook_select").html(obj.html);
					$("#wait_notebook").animate({opacity:0});
				}
			}
		},
	});
}

function add_items_to_notebook() {
	var item_list = contentPanel.grid.getSelectedItemList();
	var nid = $('#notebook_select').find('option:selected').val();
	if (nid && item_list) {
		$.ajax({
			data: {
				fname: 'add_items_to_notebook',
				nid: nid,
				item_list: item_list,
			},
			success : function(data) {
				$('#add_to_notebook_dialog').dialog('close');
			}
		});
	}
	else {
		$('#add_to_notebook_dialog').dialog('close');
	}
}

function share_dialog() {
	var item_list = contentPanel.grid.getSelectedItemList();
	if (item_list) {
		$.ajax({
			data: {
				fname: 'get_share_dialog',
				item_list: item_list,
			},
			success : function(data) {
				if (data)
					$('#share_dialog').html(data).dialog({width:'31em'}).dialog('open');
			}
		});
	}
}

function make_items_public(make_public) {
	var item_list = contentPanel.grid.getSelectedItemList();
	if (item_list) {
		$.ajax({
			data: {
				fname: 'make_items_public',
				item_list: item_list,
				make_public: make_public
			},
			success : function(data) {
				if (data)
					$('#share_dialog').html(data);
			}
		});
	}
}

function remove_items_from_user_or_group(target_item) {
	var item_list = contentPanel.grid.getSelectedItemList();
	if (target_item && item_list) {
		$.ajax({
			data: {
				fname: 'remove_items_from_user_or_group',
				target_item: target_item,
				item_list: item_list,
			},
			success : function(data) {
				if (data)
					$('#share_dialog').html(data);
			}
		});
	}
}

function add_items_to_user_or_group() {
	var item_list = contentPanel.grid.getSelectedItemList();
	var target_item = $('#share_input').data('select_id');
	var role_id = $('#share_role_select').val();
	if (target_item && item_list) {
		$.ajax({
			data: {
				fname: 'add_items_to_user_or_group',
				target_item: target_item,
				role_id: role_id,
				item_list: item_list,
			},
			success : function(data) {
				if (data) 
					$('#share_dialog').html(data);
			}
		});
	}
}

function search_share () {
	var search_term = $('#share_input').val();

	//$("#wait_notebook").animate({opacity:1});

	$.ajax({
		data: {
			fname: 'search_share',
			search_term: search_term,
			timestamp: init_timestamp('search_share')
		},
		success : function(data) {
			var obj = jQuery.parseJSON(data);
			if (obj && test_timestamp('search_share', obj.timestamp) && obj.items) {
				$("#share_input").autocomplete({source: obj.items}).autocomplete("search");
				//$("#wait_notebook").animate({opacity:0});
			}
		},
	});
}

function search_group () { // FIXME dup of above routine but for group dialog
	var search_term = $('#group_input').val();

	$.ajax({
		data: {
			fname: 'search_share',
			search_term: search_term,
			timestamp: init_timestamp('search_group')
		},
		success : function(data) {
			var obj = jQuery.parseJSON(data);
			if (obj && test_timestamp('search_group', obj.timestamp) && obj.items) {
				$("#group_input").autocomplete({source: obj.items}).autocomplete("search");
			}
		},
	});
}

function group_dialog() {
	var item_list = contentPanel.grid.getSelectedItemList();
	if (item_list) {
		$.ajax({
			data: {
				fname: 'get_group_dialog',
				item_list: item_list,
			},
			success : function(data) {
				if (data)
					$('#group_dialog').html(data).dialog({width:500}).dialog('open');
			}
		});
	}
}

function change_group_role() {
	var item_list = contentPanel.grid.getSelectedItemList();
	var role_id = $('#group_role_select').val();
	if (role_id && item_list) {
		$.ajax({
			data: {
				fname: 'change_group_role',
				target_items: item_list,
				role_id: role_id,
			},
			success : function(data) {
				if (data)
					$('#group_dialog').html(data);
			}
		});
	}
}

function add_users_to_group() {
	var item_list = contentPanel.grid.getSelectedItemList();
	var new_item = $('#group_input').data('select_id');
	if (new_item && item_list) {
		$.ajax({
			data: {
				fname: 'add_users_to_group',
				target_items: item_list,
				new_item: new_item,
			},
			success : function(data) {
				if (data) 
					$('#group_dialog').html(data);
			}
		});
	}
}

function remove_user_from_group(user_id) {
	var item_list = contentPanel.grid.getSelectedItemList();
	if (user_id && item_list) {
		$.ajax({
			data: {
				fname: 'remove_user_from_group',
				target_items: item_list,
				user_id: user_id,
			},
			success : function(data) {
				if (data) 
					$('#group_dialog').html(data);
			}
		});
	}
}

function edit_dialog() {
	if (contentPanel.selectedView == 'group') {
		group_dialog();
	}
//	else if (contentPanel.selectedView == 'notebook') {
//		add_to_notebook_dialog();
//	}
}

//function hide_top_panel() {
//	top_panel_height = $('#top_panel').height();
//	$('#top_panel').slideUp('slow',
//		function() {
//			$('#show_panel_button').show();
//		}
//	);
//	var contents_panel_height = $('#contents_table').height() + top_panel_height;
//	$('#contents_table').animate({height: contents_panel_height}, 'slow');
//}

//function show_top_panel() {
//	$('#show_panel_button').hide();
//	$('#top_panel').slideDown('slow');
//	var contents_panel_height = $('#contents_table').height() - top_panel_height;
//	$('#contents_table').animate({ height: contents_panel_height}, 'slow');
//}

function create_menu() {
	var menu = $("#create_menu");

	if (menu.is(":visible")) {
		menu.hide();
	}
	else {
		menu.show();
		menu.one("mouseleave", function() { menu.hide(); } );
	}
}

function create_group_dialog() {
	$('#edit_group_name,#edit_group_desc').val('');
	$('#create_group_dialog').dialog({width:'30em'}).dialog('open');
	$('#create_menu').hide();
}

function create_notebook_dialog() {
	$('#edit_notebook_name,#edit_notebook_desc').val('');
	$('#create_notebook_dialog').dialog({width:'30em'}).dialog('open');
	$('#create_menu').hide();
}

function add_dialog() {
	if (contentPanel.selectedView == 'group') {
		create_group_dialog();
	}
	else if (contentPanel.selectedView == 'notebook') {
		create_notebook_dialog();
    }
}

function create_new_group() {
	var name = $('#edit_group_name').val();
	if (!name) {
		alert('Please enter a group name.');
		return;
	}
	var desc = $('#edit_group_desc').val();
	var role_id = $('#edit_group_role').val();

    $.ajax({
        data: {
            fname: 'create_new_group',
            name: name,
            desc: desc,
            role_id: role_id
        },
        success: function(rc) {
            if (rc) {
            	schedule_poll(0);
            	tocPanel.selectItemType('group');
            }
        },
        complete: function() {
        	$('#create_group_dialog').dialog('close');
        }
    });
}

function create_new_notebook() {
	var name = $('#edit_notebook_name').val();
	if (!name) {
		alert('Please enter a notebook name.');
		return;
	}
	var desc = $('#edit_notebook_desc').val();
	var type_id = $('#edit_notebook_type').val();

	var item_list = contentPanel.grid.getSelectedItemList(); // optional
	
    $.ajax({
        data: {
            fname: 'create_new_notebook',
            name: name,
            desc: desc,
            type_id: type_id,
            item_list: item_list,
        },
        success: function(rc) {
            if (rc) {
            	schedule_poll(0);
            	tocPanel.selectItemType('notebook');
            }
        },
        complete: function() {
        	$('#create_notebook_dialog').dialog('close');
        }
    });
}

function upload_metadata_dialog(type) {
	$('#upload_metadata_dialog').dialog({
		width:'28em'
	}).dialog('open');
	$('#metadata_file').fileupload({
		submit: function(e,data) {
			if ('text/plain' != data.files[0].type) {
				alert('The file must be a plain text file.');
				return false;
			}
			return true;
		},
		done: function(e,data) {
			alert(data.result);
			$('#upload_metadata_dialog').dialog('close');
		},
		formData: {
			fname: 'upload_metadata',
			type: type
		}
	});
}

function send_menu() {
	var menu = $("#send_menu");

	// Positioning is done here instead of onload to work-around misplacement problem
	// due to contents title missing on page load.
	if (!pageObj.positionMenu) {
		menu.position({
			my: "left top",
			at: "left bottom",
			of: "#send_button"
		});
		pageObj.positionMenu = 1;
	}

	if (menu.is(":visible")) {
		menu.hide();
	}
	else {
		menu.show();
		menu.one("mouseleave", function() { menu.hide(); } );
	}
}

function send_items_to(page_name, format) {
	var item_list = contentPanel.grid.getSelectedItemList();
	if (item_list) {
		$.ajax({
			data: {
				fname: 'send_items_to',
				page_name: page_name,
				format: format,
				item_list: item_list,
			},
			success : function(url) {
				if (url)
					window.open(url);
			}
		});
	}
}

function toggle_star(img, id) {
	$.ajax({
		data: {
			fname: 'toggle_star',
			log_id: id,
		},
		success :  function(val) {
			$(img).attr({ src: (val == 0 ? "picts/star-hollow.png" : "picts/star-full.png") });
		}
	});
}

// For "Create New Genome" and "Create New Experiment" //FIXME merge with ContentPanel.openItem ...?
function open_item(item_type, title, link) {
	title = title + "<br><a class='xsmall' style='color:#eeeeee;' href='"+link+"' target='_blank'>[Open in new tab]</a> ";
	link = link + "&embed=1";
	console.log(link);
	var height = $(window).height() * 0.8;
	var d = $('<div class="dialog_box"><iframe src="'+link+'" height="100%" width="100%" style="border:none;"/></div>')
		.dialog({
			//title: title,
			width: '80%',
			height: height, //'80%',
			open: function() { // mdb added 10/16/16 -- fix html in dialog title bar for jQuery 3.1.1 update
                $(this).prev().find("span.ui-dialog-title").append('<span>'+title+'</span>');
            }
		})
		.dialog('open');
}

function search_metadata(type, key) {
	window.open('SearchResults.pl?s=type::' + type + ' "metadata_key::' + key + '" role::owner');
}