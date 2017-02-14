const POLL_TIME = 30*1000, // polling rate when user is not idle
	  IDLE_TIME = 30*1000; // stop polling after this lapse, then poll on next mousemove

var grid;
var infoPanel;
var tocPanel;
var timestamps = new Array();
var timers = new Array();

class UserContentPanel extends ContentPanel {
    update(viewId) {
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
    }

    refresh() {
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
    }

    fetch(sync, typeId) {
    	var self = this;

    	if (self.fetchHandler)
    	    self.fetchHandler.call(self);

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
}

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
			operations: ['share', 'organize', 'favorite', 'sendto'],
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
		elementId: 'info_panel',
		defaultInfo: default_info
	});
	
	contentPanel = new UserContentPanel({
		elementId: 'contents_panel',
		views: views,
		grid: new DataGrid({
			element: $('#contents_panel').children('.grid'),
			height: $(window).height() - 210, // this depends on the height of the header/footer and should be passed in as an argument
			filter: function(data) {
				// Filter rows based on view
				var view = views[contentPanel.getView()];
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
			selectionCallback: function(items) {
			    infoPanel.busy().update(items);
				contentPanel.renderButtons(items && items.length);
			},
			mouseOver: function(row) {
			    infoPanel.busy().scheduleUpdate([row]);
			},
			mouseOut: function() {
			    infoPanel.scheduleUpdate();
			},
			dateSortAsc:  dateSortAscending,
			dateSortDesc: dateSortDescending,
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
		})
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
			$('#search_input').val(''); //FIXME move into ContentPanel
		}
	});
	
	$('#search_input').on('keyup search', function() { //FIXME move into ContentPanel
		contentPanel.grid.search( $(this).val() );
		contentPanel.renderTitle();
	});
	
	// Get starting page from URL and initialize TOC panel
	var toc_id = coge.utils.getURLParameter('p');
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
	switch(contentPanel.getView()) {
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
 * Data Grid Row
 */

class DataGridRow { //FIXME duplicated in search-results.js
    constructor(data, type) {
        $.extend(this, data);
        if (!this.type) // mdb added condition 9/16/18 COGE-388 -- prevent native type (genome,etc) from being set to "favorite"
            this.type = type;
    }

    getFlags(opts) {
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
    }

    getDescription() {
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
    }

    _formatGenome() {
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
    }

    _formatExperiment() {
    	var descStr =
    		'<img src="picts/testtube-icon.png" width="15" height="15" style="vertical-align:middle;"/> ' +
    	   	this.name +
    	   	(this.description ? ': ' + this.description : '') +
    	   	' (v' + this.version + ', id' + this.id + ')';
    	return descStr;
    }

    _formatNotebook() {
    	var descStr =
    		'<img src="picts/notebook-icon.png" width="15" height="15" style="vertical-align:middle;"/> ' +
    		this.name +
    		(this.description ? ': ' + this.description : '');
    	return descStr;
    }

    _formatGroup() {
    	var descStr =
    		'<img src="picts/group-icon.png" width="15" height="15" style="vertical-align:middle;"/> ' +
    		this.name +
    		(this.description ? ': ' + this.description : '');;
    	return descStr;
    }

    _formatWorkflowStatus(status) {
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
    }

    _formatAnalysis() {
        var isRunning   = (this.status.toLowerCase() == 'running');
        var isCancelled = (this.status.toLowerCase() == 'cancelled');
        var star_icon    = '<img title="Favorite this analysis" src="picts/star-' + (this.is_important ? 'full' : 'hollow') + '.png" width="15" height="15" class="link" style="vertical-align:middle;" onclick="toggle_star(this, '+this.id+');" />';
        var cancel_icon  = '<img title="Cancel this analysis" class="link" height="15" style="vertical-align:middle;" src="picts/cancel.png" width="15" onclick="cancel_job_dialog('+this.id+');"/>';
        var restart_icon = '<img title="Restart this analysis" class="link" height="15" style="vertical-align:middle;" src="picts/refresh-icon.png" width="15" onclick="restart_job('+this.id+');"/>';
        var comment_icon = '<img title=' + (this.comment && this.comment.length ? JSON.stringify(this.comment) : '"Add comment"') + ' class="link" height="15" style="vertical-align:middle;" src="picts/' + (this.comment && this.comment.length ? 'comment' : 'no-comment') + '-icon.png" width="15" onclick="comment_dialog('+this.id+');" />';
        var icons = star_icon + ' ' + comment_icon + ' ' + (isCancelled ? restart_icon : '') + ' ' + (isRunning ? cancel_icon : '');
    	var descStr = icons + ' ' + this._formatWorkflowStatus(this.status) + ' ' + this.page + ' | ' + this.description + (this.comment ? ' | ' + '<span class="highlighted">' + this.comment : '') + '</span>' + ' | ' + this.elapsed + (this.workflow_id ? ' | id' + this.workflow_id : '');
    	return descStr;
    }

    _formatLoad() {
    	var descStr =
    		this._formatWorkflowStatus(this.status) + ' ' + this.page + ' | ' + this.description + ' | ' + this.elapsed + (this.workflow_id ? ' | id' + this.workflow_id : '');
    	return descStr;
    }

    getInfo() {
    	console.log('DataGridRow.getInfo');
    	var self = this;

    	return coge.utils.ajaxWithTimestamp({
    		dataType: 'json',
    		data: {
    			fname: 'get_item_info',
    			item_id: self.id,
    			item_type: self.type,
    		}
    	}).pipe(function(data) {
    		if (data)
				return data.html;
    		return;
    	});
    }

    getLink() {
    	var type = this.type;

    	if (type == 'genome')
    		return 'GenomeInfo.pl?gid=' + this.id;
    	else if (type == 'experiment')
    		return 'ExperimentView.pl?eid=' + this.id;
    	else if (type == 'notebook')
    		return 'NotebookView.pl?nid=' + this.id;
    	else
    		return this.link;
    }

    getDate() {
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

    open(url, title) {
        if (this.type == 'group')
            group_dialog();
        else if (this.type == 'analyses' || this.type == 'loads')
            window.open(this.link, '_blank');
        else {
            var title = title || this.getDescription();
            var link  = url || this.getLink();
            var flags = this.getFlags({noSpaces: 1});
            title = flags + ' ' + title + "<br><a class='xsmall' style='color:#eeeeee;' href='" + link + "' target='_blank'>[Open in new tab]</a> ";
            link = link + "&embed=1";
            console.log('DataGrid.openItem: ' + link);
            var height = $(window).height() * 0.8;
            var d = $('<div class="dialog_box"><iframe src="'+link+'" height="100%" width="100%" style="border:none;"/></div>')
                .dialog({
					autoOpen: true,
                    //title: title,
                    width: '80%',
                    height: height,
                    open: function() { // mdb added 10/16/16 -- fix html in dialog title bar for jQuery 3.1.1 update
                        $(this).prev().find("span.ui-dialog-title").append('<span>'+title+'</span>');
                    },
					close: function() {
						schedule_poll(0);
					}
                });
        }
    }
}

function dateSortAscending(x,y) {
    x = time_diff(x);
    y = time_diff(y);
    return ((x < y) ? 1 : ((x > y) ? -1 : 0));
}

function dateSortDescending(x,y) {
    x = time_diff(x);
    y = time_diff(y);
    return ((x < y) ? -1 : ((x > y) ? 1 : 0));
}

function time_diff(s) {
    var matches = s.match(/^\<\!\-\-(\d+)\-\-\>/); // the time difference is hidden in an html comment
    if (!matches || !matches[1])
        return 999999; // sort last
    var diff = matches[1];
    return parseInt(diff);
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

function init_timestamp(name) { // TODO move into class
	timestamps[name] = new Date().getTime()
	return timestamps[name];
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
	if (contentPanel.getView() == 'group') {
		create_group_dialog();
	}
	else if (contentPanel.getView() == 'notebook') {
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
	//var type_id = $('#edit_notebook_type').val(); // mdb removed 12/14/16 COGE-800

	var item_list = contentPanel.grid.getSelectedItemList(); // optional
	
    $.ajax({
        data: {
            fname: 'create_new_notebook',
            name: name,
            desc: desc,
            //type_id: type_id, // mdb removed 12/14/16 COGE-800
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

function search_metadata(type, key) {
	window.open('SearchResults.pl?s=type::' + type + ' "metadata_key::' + key + '" role::owner');
}