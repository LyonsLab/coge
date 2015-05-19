const POLL_TIME = 30*1000, // polling rate when user is not idle
	  IDLE_TIME = 30*1000; // stop polling after this lapse, then poll on next mousemove

var grid;
var infoPanel;
var timestamps = new Array();
var timers = new Array();
	
$(function() {
	// Initialize globals
	timers['item'] = null;
	init_timestamp('idle');
	
	// Initialize AJAX
	$.ajaxSetup({
		type: "GET",
		url: PAGE_NAME,
		dataType: "html",
		cache: false,
	});
	
	// Initialize dialog boxes
	$(".dialog_box").dialog({autoOpen: false, resizable: false});

	// Initialize fileupload plugin
	$('#input_upload_file').fileupload({
    	dataType: 'json',
    	formData: {
    		fname: 'upload_image_file',
    	},
       	add:
    		function(e, data) {
				if ( verify_image_file(data.files[0]) ) {
					$('#user_image').attr('src', 'picts/ajax-loader-large.gif');
					data.submit();
				}
    		},
		done:
			function(e, data) {
				if (data.result && data.result.link) {
					$('#user_image').attr('src', data.result.link);
				}
			},
	});

	// Initialize dropdown menus
	$("#create_menu").menu()
		.position({
			my: "left top",
			at: "left bottom",
			of: "#create_button"
		})
		.css({ position: 'absolute' });
	
	$("#send_menu").menu().css({ position: 'absolute', width: '90px' });

	// Get starting page from URL if set
	var toc_id = parseInt( getURLParameter('p') );
	if (!toc_id || toc_id == 'null') {
		toc_id = ITEM_TYPES.mine;
	}

	// Initialize main panels
	grid = new DataGrid({
		elementId: 'contents_table'
	});
	
	infoPanel = new InfoPanel({
		elementId: 'info_panel',
		grid: grid
	});
	
	//exit_item(); // initialize info panel
	toc_select(toc_id); // initialize toc panel
	
	// Setup idle timer
//	$(document).mousemove(function() {
//		var currentTime = new Date().getTime();
//		var idleTime = currentTime - timestamps['idle'];
//		timestamps['idle'] = currentTime;
//
//		if (idleTime > IDLE_TIME) {
//			// User was idle for a while, refresh page
//			schedule_poll(0);
//		}
//	});

	// Initiate refresh loop
	//schedule_poll();

	// Initialize add-to-notebook dialog
	window.setTimeout(search_notebooks, 1000);
	
	// Initialize confirm cancel job dialog
    $("#cancel_dialog").dialog({
    	modal: true,
    	buttons: {
    		No: function() {
    			$(this).dialog("close");
    		},
    		Yes: function() {
    			cancel_job( $(this).data("job_id") );
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

function formatDate(dateStr) {
	const MONTHS = [ 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec' ];
	if (!dateStr || dateStr === '0000-00-00 00:00:00')
		return '';
	
	dateStr = dateStr.replace(/-/g, '/'); // needed for Firefox & Safari
	var date = new Date(dateStr);
	var today = new Date();
	var diffDays = Math.round(Math.abs((today.getTime() - date.getTime())/(24*60*60*1000)));
	var dateStr;
	if (diffDays == 0) // same day as today
		dateStr = (date.getHours() % 12) + ':' + pad(date.getMinutes(), 2) + ' ' + (date.getHours() < 12 ? 'am' : 'pm');
	else if (diffDays == 1) // yesterday
		dateStr = 'Yesterday';
	else if (diffDays <= 4) // last several days
		dateStr = diffDays + ' days ago';
	else if (date.getFullYear() == today.getFullYear()) // same year 
		dateStr = MONTHS[date.getMonth()] + ' ' + date.getDate()
	else // last year or older
		dateStr = MONTHS[date.getMonth()] + ' ' + date.getDate() + ', ' + date.getFullYear();

	return dateStr;
}

function formatGenome(genome) {
	var descStr = 
	   	(genome.restricted ? '&reg; '  : '') +
	   	(genome.organism ? genome.organism : '') + 
	   	(genome.name ? ' (' + genome.name + ')' : '') +
	   	(genome.description ? ': ' + genome.description : '') +
	   	' (v' + genome.version + ', id' + genome.id + ')';
	return descStr;
}

function getURLParameter(name) {
    return decodeURI(
        (RegExp(name + '=' + '(.+?)(&|$)').exec(location.search)||[,null])[1]
    );
}

function poll(sync) {
	// Refresh contents
	get_contents(sync, pageObj.content_type);
	
	// Refresh item info cache
	grid.reset();
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

function toc_toggle_children(toc_id, num_children) {
	$('#toc_'+toc_id)
	.next()
	.attr('src', function(idx, oldSrc) {
        return (oldSrc == 'picts/arrow-down-icon.png' ? 'picts/arrow-right-icon.png' : 'picts/arrow-down-icon.png');
    });
	$('#toc_'+toc_id)
		.parent()
		.nextAll().slice(0, num_children)
		.toggle();
}

function default_info() {
	console.log('default_info');
	var text = "";
	switch(pageObj.content_type) {
		case ITEM_TYPES.activity_summary:
			text = "Here is a summary of all analyses you have performed.";
			break;
		case ITEM_TYPES.activity_analyses:
			text = "These are the analyses you have performed or started.<br><br>" + 
				"Select an analysis to open the current progress or finished result in a new tab.<br><br>" +
				"Use the icons to the left of each analysis to 'Favorite' it, add comments, or cancel (if running).";
			break;
		case ITEM_TYPES.activity_loads:
			text = "These are the data loading workflows you have performed or started.<br><br>" +
				"Select an item to open the current progress or finished result in a new tab.";
			break;
		case ITEM_TYPES.activity_viz:
			text = "Woah, cool!";
			break;
		case ITEM_TYPES.trash:
			text = "These are items you deleted.<br><br>" +
				"Hover over an item to view additional info. Select one or more items to undelete.";
			break;
		case ITEM_TYPES.shared:
			text = "These are data items that your collaborators shared with you.<br><br>" +
				"Hover over an item to view additional info. Select one or more items to share with others or add to a notebook.";
			break;
		case ITEM_TYPES.mine:
		case ITEM_TYPES.notebook:
		case ITEM_TYPES.genome:
		case ITEM_TYPES.experiment:
			text = "<p>These are data items that you added to the system.</p>"
				+ "<p><b>Hover over</b> an item to view additional info.</p>"
                + "<p><b>Single-click</b> to select one or more items to share, organize, delete, or send them to one of CoGe's tools.</p>"
                + "<p><b>Double-click</b> an item for a detailed view of the item.</p>";
			break;
		case ITEM_TYPES.group:
			text = "You are a member of these user groups.<br><br>" +
				"Hover over a group to view additional info. Select one or more groups to edit or delete."; 
			break;
	}
	$('#info_panel').html(text);
}

/*
 * Data Grid
 */

function DataGrid(params) {
	this.element = $('#'+params.elementId);
	this.selectedItemId = null;
	this.initialize();
}

$.extend(DataGrid.prototype, {
	initialize: function() {
		var self = this;
		this.element.html('<table cellpadding="0" cellspacing="0" border="0" class="dt-cell hover compact row-border"></table>');
		
		var dataTable = this.dataTable = this.element.children('table').dataTable({
			"paging": false,
			"info" : false,
			"searching": false,
			"sScrollY": $(window).height() - 245, // this depends on the height of the header/footer
			"columns": [
	            { 	title: "Name", 
	            	targets: 0,
	            	type: "string",
	            	data: null, // use full data object
	            	render: function(data, type, row, meta) {
	            		return formatGenome(data);
	            	}
	            },
	            { 	title: "Date added", 
	            	targets: 1, 
	            	type: "date",
	            	data: "date",
	            	width: "100px",
	            	render: function(data, type, row, meta) {
	            		return formatDate(data);
	            	}
	            }
			]
		});
		
		// Handle row selection event
		var dataTableBody = dataTable.children('tbody');
		dataTableBody.on('click', 'tr', function() {
			var tr = this;
			var row = dataTable.api().row(tr);
			
	        if ( $(tr).hasClass('selected') ) { // unselect
	            $(tr).removeClass('selected');
	        }
	        else { // select
	        	self.dataTable.$('tr.selected').removeClass('selected'); // unselect all
	            $(tr).addClass('selected'); // select item
	        }
	        
	        self.selectItem(row.data());
		});
		
		// Handle row hover events
		dataTableBody.on('mouseover', 'tr', function () {
			var tr = $(this).closest('tr');
	        var row = dataTable.api().row(tr).data();
	        row.enter.apply(row);
	    });
		dataTableBody.on('mouseout', 'tr', function () {
			var tr = $(this).closest('tr');
	        var row = dataTable.api().row(tr).data();
			row.exit.apply(row);
	    });	

		$.ajax({
			dataType: 'json',
			data: {
				fname: 'get_contents2',
			},
			success : function(data) {
				//console.log(data);
				var rows = data.map(function(obj) {
					return new DataGridRow(obj, self);
				});
				dataTable.api()
					.clear()
					.rows.add(rows).draw();
			}
		});		
    },
    
    reset: function() {
    	
    },
    
    getSelectedItems: function() {
    	console.log('getSelectedItems');
    	//return $('#contents_table input[type=checkbox]').filter(':checked').filter(':not(:hidden)');
    	return this.dataTable.api().rows('.selected').data();//.filter(':not(:hidden)')
    },

    clearSelection: function() {
    	$('.coge-list-item,.coge-selected').removeClass('coge-selected');
    	$('#contents_table input[type=checkbox]').filter(':checked').prop('checked', false);
    	update_icons();
        selection_hint();
        this.selectedItemId = null;
    },
    
    clickItem: function(item) {
    	$('#'+item.id).toggleClass('coge-selected')
    	infoPanel.update();
    },

    selectItem: function(item) {
    	console.log('selectItem');
    	//this.clearSelection();
    	
//    	$('#'+item.id)
//    		.addClass('coge-selected')
//    		.find('input[type=checkbox]')
//    		.prop('checked', true);
    	this.selectedItemId = item.id;
    	
    	infoPanel.busy();
    	infoPanel.update();
    },

    openItem: function(item_type, title, link) {
    	if (item_type == ITEM_TYPES.group) // FIXME this is a kludge
    		group_dialog();
    	else {
            if (!link) {
                return alert("The following link could not be generated");
            }

    		title = title + "<br><a class='xsmall' href='"+link+"' target='_blank'>[Open in new tab]</a> ";
    		link = link + "&embed=1";
    		console.log(link);
    		var height = $(window).height() * 0.8;
    		var d = $('<div class="dialog_box"><iframe src="'+link+'" height="100%" width="100%" style="border:none;"/></div>')
    			.dialog({
    				title: title,
    				width: '80%',
    				height: height//'80%'
    			})
    			.dialog('open');
    	}
    }
});

/* 
 * Data Grid Row
 */
		
function DataGridRow(data, grid) {
	$.extend(this, data);
	if (grid)
		this.grid = grid; // parent DataGrid object
    this.initialize();
}

$.extend(DataGridRow.prototype, {
	initialize: function() {
    },

    getInfo: function() {
    	var self = this;
    	
//    	if (self.info)
//    		return $.Deferred({done: function() { return self.info; }}).promise();
    	
    	return $.ajax({
    		dataType: 'json',
    		data: {
    			fname: 'get_item_info',
    			item_id: self.id,
    			item_type: self.type,
    			timestamp: init_timestamp('get_item_info')
    		}
//    		success : function(data) {
//    			if (data && test_timestamp('get_item_info', data.timestamp))
//    				self.info = data.html;
//    		}
    	}).pipe(function(data) {
    		if (data && test_timestamp('get_item_info', data.timestamp))
				return data.html;
    		return;
    	});
    },

    enter: function() {
    	console.log('DataGridRow.enter');
    	var self = this;
    	
    	if (self.grid && self.grid.selectedItemId) // Do nothing if row currently selected
    		return;

    	if (timers['item'])
    		window.clearTimeout(timers['item']);

    	timers['item'] = window.setTimeout( 
    		function() { 
    			infoPanel.busy().update();
    		},
    		500
    	);
    },
    
    exit: function() {
    	console.log('DataGridRow.exit');
    	if (self.grid && self.grid.selectedItemId) // Do nothing if row currently selected
    		return;
    	
    	if (timers['item'])
    		window.clearTimeout(timers['item']);
    	
    	timers['item'] = window.setTimeout( default_info, 500 );
    }
});

/* 
 * Info Panel
 */

function InfoPanel(params) {
	this.element = $('#'+params.elementId);
	this.grid = params.grid; // parent DataGrid object
    this.initialize();
}

$.extend(InfoPanel.prototype, {
	initialize: function() {
    },
    
    busy: function() {
    	this.element.html('<img src="picts/ajax-loader.gif"/>');
    	return this;
    },
    
    update: function() {
    	console.log('InfoPanel.update');
    	var self = this;
    	var selected = grid.getSelectedItems();
    	var num_items = selected.length;
    	console.log(num_items);
    	console.log(selected);

    	if (num_items > 0) {
    		if (num_items == 1) {
    			var item = selected[0];//.parentNode;
    			item.getInfo().pipe(function(info) {
    				//console.log(info);
    				self.element.html(info);
    			});
    		}
    		else
    			self.element.html(num_items + ' item' + (num_items > 1 ? 's' : '') + ' selected.<br><br>Click an action icon at the top to share, organize, delete, or analyze.');
    	}
    	else {
    		//default_info();
    		parent.selectedItemId = null;
    	}
    	
    	update_icons();
        selection_hint();    	
    }
});

function get_item_type(obj) {
	return obj.id.match(/content_\w+_(\w+)/)[1];
}

function sync_items(html) {
	var content1 = $('#contents_table .coge-list-item');
	var content2 = $(html).filter('.coge-list-item'); // FIXME: this is slow

	var insertIndex = 0;
	content2.each(
		function() {
			var match = document.getElementById(this.id);
			if (!match) // item doesn't exist
				$(this).insertBefore( content1.get(insertIndex) );
			else { // item exists
				var src_info = $(this).find('span[name="info"]').html();
				var dest = $(match).find('span[name="info"]');
				if (dest.html() !== src_info)
					dest.html(src_info);
				insertIndex++;
			}
		}
	);
}

function get_contents(sync, item_type) {
	$('#refresh_label').show();
	
	var lastUpdate = (sync ? timestamps['lastUpdate'] : 0);

	// Clear pending refresh
	cancel_poll();

	$.ajax({
		data: {
			fname: 'get_contents',
			item_type: item_type,
			last_update: lastUpdate,
			timestamp: init_timestamp('get_contents')
		},
		success : function(data) {
			var obj = jQuery.parseJSON(data);
			if (obj) {
				if (test_timestamp('get_contents', obj.timestamp)) {
					if (sync) { // merge with existing contents
						sync_items(obj.html);
					}
					else { // replace existing contents
						$('#contents_table').html(obj.html);
					}
					timestamps['lastUpdate'] = obj.lastUpdate;
					filter_contents();
				}
			}

			// Setup next refresh
			schedule_poll();
		},
		complete : function() {
			$('#refresh_label').hide();
		}
	});
}

function filter_contents() {
	var search_term = $('#search_input').val();
	search_term = search_term.toLowerCase();

	$('#contents_table div.coge-list-item').each(
		function() {
			var item_type = get_item_type(this);
			var show;

			if (pageObj.content_type == ITEM_TYPES.group) {
				show = item_type == ITEM_TYPES.group
						&& !$(this).hasClass('deleted');
			}
			else if (pageObj.content_type == ITEM_TYPES.mine) {
				show = !$(this).hasClass('deleted')
						&& !$(this).hasClass('shared')
						&& item_type != ITEM_TYPES.activity_summary 
						&& item_type != ITEM_TYPES.activity_analyses
						&& item_type != ITEM_TYPES.activity_loads
						&& item_type != ITEM_TYPES.activity_viz
						&& item_type != ITEM_TYPES.group
						&& item_type != ITEM_TYPES.notebook;
			}
			else if (pageObj.content_type == ITEM_TYPES.shared) {
				show = $(this).hasClass('shared')
						&& !$(this).hasClass('deleted')
						&& item_type != ITEM_TYPES.activity_summary 
						&& item_type != ITEM_TYPES.activity_analyses
						&& item_type != ITEM_TYPES.activity_loads
						&& item_type != ITEM_TYPES.activity_viz;
			}
			else if (pageObj.content_type == ITEM_TYPES.trash) {
				show = $(this).hasClass('deleted');
			}
			else {
				show = (item_type == pageObj.content_type
						&& !$(this).hasClass('deleted')
						&& !$(this).hasClass('shared'));
			}

			if (show && search_term) {
				show = (this.innerHTML.toLowerCase().indexOf(search_term) >= 0);
			}

			if (show) { $(this).show(); }
			else { $(this).hide(); }
		}
	);
}

function update_icons(on) {
	var selected = grid.getSelectedItems();
	if ( selected.length > 0) 
		$('.item-button:not(#add_button)').removeClass('coge-disabled');
	else
		$('.item-button:not(#add_button)').addClass('coge-disabled');
}

function selection_hint() {
	var selected = grid.getSelectedItems();
	if ( selected.length == 1)
		$('#selected-hint').removeClass('hidden');
	else
		$('#selected-hint').addClass('hidden');
}

function toc_select(toc_id) {
	if (toc_id == pageObj.content_type) 
		return;
	
	pageObj.content_type = toc_id;

	// Clear search bar
	$('#search_input').show().val('');
	
	// Show/hide action icons based on type of data
	switch(toc_id) {
		case ITEM_TYPES.mine:
		case ITEM_TYPES.genome:
		case ITEM_TYPES.experiment:
			$('#undelete_button,#edit_button,#add_button').hide();
			$('#share_button,#notebook_button,#delete_button,#send_button').show();
			break;
		case ITEM_TYPES.notebook:
			$('#undelete_button,#edit_button,#notebook_button').hide();
			$('#share_button,#delete_button,#send_button,#add_button').show();
			break;
		case ITEM_TYPES.trash:
			$('#share_button,#notebook_button,#edit_button,#delete_button,#send_button,#add_button').hide();
			$('#undelete_button').show();
			break;
		case ITEM_TYPES.group:
			$('#notebook_button,#send_button,#share_button,#undelete_button').hide();
			$('#edit_button,#delete_button,#add_button').show();
			break;
		case ITEM_TYPES.shared:
			$('#edit_button,#delete_button,#send_button,#undelete_button,#add_button').hide();
			$('#share_button,#notebook_button').show();
			break;
		case ITEM_TYPES.activity_analyses:
		case ITEM_TYPES.activity_loads:
			$('#share_button,#notebook_button,#edit_button,#delete_button,#send_button,#undelete_button,#add_button').hide(); // hide all
			break;
		case ITEM_TYPES.activity_summary:
		case ITEM_TYPES.activity_viz:
			$('#share_button,#notebook_button,#edit_button,#delete_button,#send_button,#undelete_button,#add_button').hide(); // hide all
			$('#search_input').hide();
			break;
		default:
			$('#share_button,#notebook_button,#edit_button,#delete_button,#send_button,#undelete_button,#add_button').hide(); // hide all
	}
	
	// Icons are set to invisible on load to prevent flickering
	$('.item-button').removeClass('invisible');

	// Refilter based on content selection
	filter_contents();
	
	// Clear selected items
	grid.clearSelection();
	default_info();
	
	// Show TOC item as selected
	$('#toc_'+toc_id)
		.parent()
		.css({'font-weight':'bold'})
		.siblings()
		.css({'font-weight':'normal'});

	// Update title
	var title = $('#toc_label_'+toc_id).html();
	$('#contents_title').html(title);

	// Update browser url
	window.history.pushState({}, "", PAGE_NAME + "?p="+toc_id);
}

function delete_items() {
	var selected = grid.getSelectedItems();
	if (selected.length) {
		selected.parent().fadeOut('fast',
			function() {
				selected.parent().addClass('deleted');
				selected.prop('checked', false);
				infoPanel.update(); // reset button states
			}
		);

		var item_list = selected.map(function(){return this.parentNode.id;}).get().join(',');
		$.ajax({
			data: {
				fname: 'delete_items',
				item_list: item_list
			},
			success : function(data) {
			}
		});
	}
}

function undelete_items() {
	var selected = grid.getSelectedItems();
	if (selected.length) {
		selected.parent().fadeOut('fast',
			function() {
				selected.parent().removeClass('deleted');
				selected.prop('checked', false);
				infoPanel.update(); // reset button states
			}
		);

		var item_list = selected.map(function(){return this.parentNode.id;}).get().join(',');
		$.ajax({
			data: {
				fname: 'undelete_items',
				item_list: item_list
			},
			success : function(data) {
			}
		});
	}
}

/*
function cancel_jobs() {
	var selected = get_selected_items();
	if (selected.length) {
		var item_list = selected.map(function(){return this.parentNode.id;}).get().join(',');
		$.ajax({
			data: {
				fname: 'cancel_jobs',
				item_list: item_list
			},
			success : function(rc) {
				if (rc) {
	            	poll();
	            }
			}
		});
	}
}
*/

function cancel_job_dialog(id) {
	if (id) {
		$('#cancel_dialog')
			.data("job_id", id)
			.dialog('open');
	}
}

function cancel_job(id) {
	if (id) {
		$.ajax({
			data: {
				fname: 'cancel_job',
				workflow_id: id
			},
			success : function(rc) {
				if (rc) {
					poll(0);//schedule_poll(0); // FIXME mdb changed 10/7/14, reevaluate someday
	            }
			}
		});
	}
}

function restart_job(link) {
	if (link) {
		// Doesn't work for http vs. https:
//		var element = document.createElement("iframe"); 
//		element.setAttribute('id', 'myframe');
//		element.setAttribute('src', link);
//		document.body.appendChild(element);

		// Open the link in a "hidden" window.
		// Workaround for restarting a workflow until JEX implements a 
		// restart command.
		var w = window.open(link,'_blank', 'toolbar=no,status=no,menubar=no,scrollbars=no,resizable=no,left=10000,top=10000,width=1,height=1,visible=none', ''); 
		setTimeout(function() {
				poll(0);//schedule_poll(0); // FIXME mdb changed 10/7/14, reevaluate someday
				w.close();
			},
			5*1000
		);
	}
}

function comment_dialog(id, val) {
	if (id) {
		$('#comment_dialog').find("input").first().val(val);
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
				schedule_poll(0);
            }
		}
	});
}

function add_to_notebook_dialog() {
	var selected = grid.getSelectedItems();
	if (selected.length) {
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
	var search_term = $('#notebook_search_input').attr('value');

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
	var selected = grid.getSelectedItems();
	var nid = $('#notebook_select').find('option:selected').val();
	if (nid && selected.length) {
		var item_list = selected.map(function(){return this.parentNode.id;}).get().join(',');
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
	var selected = grid.getSelectedItems();
	if (selected.length) {
		var item_list = selected.map(function(){return this.parentNode.id;}).get().join(',');
		$.ajax({
			data: {
				fname: 'get_share_dialog',
				item_list: item_list,
			},
			success : function(data) {
				$('#share_dialog').html(data).dialog({width:500}).dialog('open');
			}
		});
	}
}

function remove_items_from_user_or_group(target_item) {
	var selected = grid.getSelectedItems();
	if (target_item && selected.length) {
		var item_list = selected.map(function(){return this.parentNode.id;}).get().join(',');
		$.ajax({
			data: {
				fname: 'remove_items_from_user_or_group',
				target_item: target_item,
				item_list: item_list,
			},
			success : function(data) {
				if (data) {
					$('#share_dialog').html(data);
//					infoCache.refresh();
				}
			}
		});
	}
}

function add_items_to_user_or_group() {
	var selected = grid.getSelectedItems();
	var target_item = $('#share_input').data('select_id');
	var role_id = $('#share_role_select').val();
	if (target_item && selected.length) {
		var item_list = selected.map(function(){return this.parentNode.id;}).get().join(',');
		$.ajax({
			data: {
				fname: 'add_items_to_user_or_group',
				target_item: target_item,
				role_id: role_id,
				item_list: item_list,
			},
			success : function(data) {
				if (data) {
					$('#share_dialog').html(data);
//					infoCache.refresh();
				}
			}
		});
	}
}

function search_share () {
	var search_term = $('#share_input').attr('value');

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
	var search_term = $('#group_input').attr('value');

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
	var selected = grid.getSelectedItems();
	if (selected.length) {
		var item_list = selected.map(function(){return this.parentNode.id;}).get().join(',');
		$.ajax({
			data: {
				fname: 'get_group_dialog',
				item_list: item_list,
			},
			success : function(data) {
				$('#group_dialog').html(data).dialog({width:500}).dialog('open');
			}
		});
	}
}

function change_group_role() {
	var selected = grid.getSelectedItems();
	var role_id = $('#group_role_select').val();
	if (role_id && selected.length) {
		var target_items = selected.map(function(){return this.parentNode.id;}).get().join(',');
		$.ajax({
			data: {
				fname: 'change_group_role',
				target_items: target_items,
				role_id: role_id,
			},
			success : function(data) {
				if (data) {
					$('#group_dialog').html(data);
				}
			}
		});
	}
}

function add_users_to_group() {
	var selected = grid.getSelectedItems();
	var new_item = $('#group_input').data('select_id');
	if (new_item && selected.length) {
		var target_items = selected.map(function(){return this.parentNode.id;}).get().join(',');
		$.ajax({
			data: {
				fname: 'add_users_to_group',
				target_items: target_items,
				new_item: new_item,
			},
			success : function(data) {
				if (data) {
					$('#group_dialog').html(data);
				}
			}
		});
	}
}

function remove_user_from_group(user_id) {
	var selected = grid.getSelectedItems();
	if (user_id && selected.length) {
		var target_items = selected.map(function(){return this.parentNode.id;}).get().join(',');
		$.ajax({
			data: {
				fname: 'remove_user_from_group',
				target_items: target_items,
				user_id: user_id,
			},
			success : function(data) {
				if (data) {
					$('#group_dialog').html(data);
				}
			}
		});
	}
}

function edit_dialog() {
	if (pageObj.content_type == ITEM_TYPES.group) {
		group_dialog();
	}
//	else if (pageObj.content_type == ITEM_TYPES.notebook) {
//		add_to_notebook_dialog();
//	}
}

function hide_top_panel() {
	top_panel_height = $('#top_panel').height();
	$('#top_panel').slideUp('slow',
		function() {
			$('#show_panel_button').show();
		}
	);
	var contents_panel_height = $('#contents_table').height() + top_panel_height;
	$('#contents_table').animate({height: contents_panel_height}, 'slow');
}

function show_top_panel() {
	$('#show_panel_button').hide();
	$('#top_panel').slideDown('slow');
	var contents_panel_height = $('#contents_table').height() - top_panel_height;
	$('#contents_table').animate({ height: contents_panel_height}, 'slow');
}

function show_recent_activity() {
	$.ajax({
		data: {
			fname: 'get_logs',
			type: 'recent'
		},
		success : function(data) {
			$('#logs').html(data);
			$('#recent').css('font-weight', 'bold');
			$('#important').css('font-weight', 'normal');
		}
	});
}

function show_important_activity() {
	$.ajax({
		data: {
			fname: 'get_logs',
			type: 'important'
		},
		success : function(data) {
			if (data) {
				$('#logs').html(data);
			}
			else {
				$('#logs').html("<span style='font-style:italic;color:gray;'>None ... click the <img src='picts/star-hollow.png'> icon  in your <a href='History.pl' target='_blank'>History</a> to mark items as important.</span>");
			}
			$('#recent').css('font-weight', 'normal');
			$('#important').css('font-weight', 'bold');
		}
	});
}

function select_image_file() {
	$('#input_upload_file').click();
}

function verify_image_file(file) {
	var ext = file.name.split('.').pop();
	if (ext != 'jpg' && ext != 'gif' && ext != 'png') {
		alert('Error: specified file is not an image');
		return 0;
	}

	if (file.size > 2*1024*1024) {
		alert('Error: image file is too large (>2MB)');
		return 0;
	}

	return 1;
}

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
	$('#create_group_dialog').dialog({width:'27em'}).dialog('open');
	$('#create_menu').hide();
}

function create_notebook_dialog() {
	$('#edit_notebook_name,#edit_notebook_desc').val('');
	$('#create_notebook_dialog').dialog({width:'27em'}).dialog('open');
	$('#create_menu').hide();
}

function add_dialog() {
	if (pageObj.content_type == ITEM_TYPES.group) {
		create_group_dialog();
	} else if (pageObj.content_type == ITEM_TYPES.notebook) {
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
            	toc_select(ITEM_TYPES.group);
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

	var item_list;
	var selected = grid.getSelectedItems();
	if (selected.length) {
		item_list = selected.map(function(){return this.parentNode.id;}).get().join(',');
	}

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
            	toc_select(ITEM_TYPES.notebook);
            }
        },
        complete: function() {
        	$('#create_notebook_dialog').dialog('close');
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
	var selected = grid.getSelectedItems();
	if (selected.length) {
		var item_list = selected.map(function(){return this.parentNode.id;}).get().join(',');
		$.ajax({
			data: {
				fname: 'send_items_to',
				page_name: page_name,
				format: format,
				item_list: item_list,
			},
			success : function(url) {
				if (url) {
					window.open(url);
				}
			}
		});
	}
}

function toggle_star(img) {
	$.ajax({
		data: {
			fname: 'toggle_star',
			log_id: img.name,
		},
		success :  function(val) {
			$(img).attr({ src: (val == 0 ? "picts/star-hollow.png" : "picts/star-full.png") });
		}
	});
}
