var ITEM_TYPE_USER = 5; //TODO: This is duplicated elsewhere, move to a common location
var timestamps = new Array();
var jobs_timers = new Array();
var hist_timers = new Array();
var current_tab = 0;
var previous_search = ""; //indicates the previous search term, used to refresh after a delete
var jobs_updating = true;
var jobs_init = false;
var running_only = 1;
var hist_init = false;
var user_graph_init = false;
var group_graph_init = false;
var hist_updating = true;
var hist_entries = 0;
var last_hist_update;
var IDLE_TIME = 30*1000; // stop polling after this lapse, then poll on next mousemove
var reports_grid;
var tree;

$(function () {
	$( "#tabs" ).tabs({
		select: function(event, ui) {
            var theSelectedTab = ui.index;
            change_tab(theSelectedTab);
            schedule_update(5000);
        },
		show: function(event, ui) {
			if(current_tab == 1 && !jobs_init) {
				init_jobs_grid();
			}
			if(current_tab == 2 && !hist_init) {
				init_hist_grid();
			}
			if(current_tab == 4 && !reports_grid) {
				init_reports();
			}
			if(current_tab == 5) {
				init_taxon_tree("");
			}
		}
    });
	
	$( "#tabs" ).show();
	
	timestamps['idle'] = new Date().getTime();
	
	// Configure dialogs
    $(".dialog_box").dialog({autoOpen: false, width: 500});
    $("#job_search_bar").keyup(function (e) {
        Slick.GlobalEditorLock.cancelCurrentEdit();

        if (e.which == 27) { // Clear on Esc
            this.value = "";
        }
        update_filter();
    });
    $("#show_select,#job_search_type").change(function(e) {
        update_filter();
    });
    $("#job_update_checkbox").change(function(e) {
    	toggle_job_updater();
    });
    $("#running_checkbox").change(function(e) {
    	toggle_running();
    });
    $("#hist_update_checkbox").change(function(e) {
    	toggle_hist_updater();
    });
    
    //Setup idle timer
    $(document).mousemove(function() {
		var currentTime = new Date().getTime();
		var idleTime = currentTime - timestamps['idle'];
		timestamps['idle'] = currentTime;

		if (idleTime > IDLE_TIME) {
			// User was idle for a while, refresh page
			schedule_update(5000);
		}
	});
    
    //var histo_test = new Histogram("histo_dialog");
});

//Initialize the Jobs tab
function init_jobs_grid() {
	var searchFilter = function(item, args) {
        var link = item['link'] ? item['link'].toLowerCase() : '',
            tool = item['tool'] ? item['tool'].toLowerCase() : '',
            status = item['status'] ? item['status'].toLowerCase() : '',
            started = item['started'] ? item['started'].toLowerCase() : '',
            completed = item['completed'] ? item['completed'].toLowerCase() : '',
            user = item['user'] ? item['user'].toLowerCase() : '';

        if (args.searchType == 1) {
            if (args.searchString != "" &&
                link.indexOf(args.searchString) == -1 &&
                tool.indexOf(args.searchString) == -1 &&
                status.indexOf(args.searchString) == -1 &&
                started.indexOf(args.searchString) == -1 &&
                completed.indexOf(args.searchString) == -1
                && user.indexOf(args.searchString) == -1) {
                return false;
            }
        } else {
            if (args.searchString != "" &&
                link.indexOf(args.searchString) != -1 ||
                tool.indexOf(args.searchString) != -1 ||
                status.indexOf(args.searchString) != -1 ||
                started.indexOf(args.searchString) != -1 ||
                completed.indexOf(args.searchString) != -1

                || user.indexOf(args.searchString) != -1 ) {
                return false;
            }
        }

        return true;
    };

    var job_options = {
        editable: true,
        enableCellNavigation: true,
        asyncEditorLoading: true,
        forceFitColumns: true,
        filter: searchFilter,
        comparator: coge.ascending,
    };

    var checkbox = new Slick.CheckboxSelectColumn({
            cssClass: 'slick-cell-checkboxsel'
    });

    var linkformatter = function(row, cell, value, columnDef, dataContext) {
        return '<a href="' + dataContext['link'] + '" target="_blank">'
        + dataContext['link'] + '</a>'
    }
    var job_columns = [
        checkbox.getColumnDefinition(),
        {id: 'id', name: 'Id', field: 'workflow_id', maxWidth: 50, sortable: true},
        {id: 'started', name: 'Started', field: 'started', minWidth: 75,
            sortable: true},
        {id: 'completed', name: 'Completed', field: 'completed', minWidth: 75,
            sortable: true},
        {id: 'elapsed', name: 'Elapsed', field: 'elapsed', minWidth: 55,
            sortable: true},
        {id: 'user', name: 'User', field: 'user', sortable: true, minWidth: 75},
        {id: 'tool', name: 'Tool', field: 'tool', minWidth: 75,
            sortable: true},
        {id: 'link', name: 'Link to Analysis', field: 'link', minWidth: 250,
            sortable: false, formatter: linkformatter },
        {id: 'status', name: 'Status', field: 'status', minWidth: 75,
            sortable: true}
    ];

    window.jobs = new coge.Grid('#jobs', job_options, job_columns);
    jobs.grid.registerPlugin(checkbox);
    get_jobs();
}

//Initialize the History tab
function init_hist_grid() {
	var hist_filter = function(item, args) {
    	var date_time 	= (item['date_time'] ? item['date_time'].toLowerCase() : '');
    	var user_name 	= (item['user'] ? item['user'].toLowerCase() : '');
    	var description = (item['description'] ? item['description'].toLowerCase() : '');
    	var page 		= (item['page'] ? item['page'].toLowerCase() : '');
    	var link 		= (item['link'] ? item['link'].toLowerCase() : '');
    	var comment 	= (item['comment'] ? item['comment'].toLowerCase() : '');

    	var show = 1;
    	if (args.show != 0) {
    		if (args.show == -1) { // Starred
    			show = item['starred'];
    		}
    		else if (args.show == -2) { // Commented
    			show = comment;
    		}
    		else if (args.show == -3) { // Mine
    			show = (user_name == '<TMPL_VAR NAME="USER_NAME">');
    		}
    		else if (args.show > 0) { // Time Range
    			var diff = new Date() - new Date(date_time.replace(/-/g, '/'));
    			show = (diff <= args.show*60*60*1000);
    		}
    	}
    	if (!show) {
    		return false;
    	}

    	if (args.searchString != "") {
    		//FIXME optimize
    		if (args.searchType == 1) { // Contains
    			if (date_time.indexOf(args.searchString) == -1 &&
    				user_name.indexOf(args.searchString) == -1 &&
    				description.indexOf(args.searchString) == -1 &&
    				page.indexOf(args.searchString) == -1 &&
    				link.indexOf(args.searchString) == -1 &&
    				comment.toLowerCase().indexOf(args.searchString) == -1 )
    			{
    				return false;
    			}
    		}
    		else { // Does not contain
    			if (date_time.indexOf(args.searchString) != -1 ||
    				user_name.indexOf(args.searchString) != -1 ||
    				description.indexOf(args.searchString) != -1 ||
    				page.indexOf(args.searchString) != -1 ||
    				link.indexOf(args.searchString) != -1 ||
    				comment.toLowerCase().indexOf(args.searchString) != -1 )
    			{
    				return false;
    			}
    		}
    	}

    	return true;
    };
    
    var hist_options = {
    		editable: true,
    		enableCellNavigation: true,
    		asyncEditorLoading: true,
    		forceFitColumns: true,
    		filter: hist_filter,
            comparator: coge.ascending,
    };
    
    var hist_columns = [
               	{id: "starred", name: "", field: "starred", maxWidth: 25, cssClass: "cell-centered",
               		formatter: function (row, cell, value, columnDef, dataContext) {
               			if (value) {
               				return '<img id="'+dataContext['id']+'" src="picts/star-full.png" onclick="toggle_star(this);">'
               			}
               			return '<img id="'+dataContext['id']+'" src="picts/star-hollow.png" onclick="toggle_star(this);">';
               		}},
               	{id: "date_time", name: "Date/Time", field: "date_time", minWidth: 160, maxWidth: 160, sortable: true, cssClass: "cell-centered"},
               	{id: "user", name: "User", field: "user", minWidth: 30, maxWidth: 80, sortable: true, cssClass: "cell-normal"/*,
               		formatter: function ( row, cell, value, columnDef, dataContext ) {
                           return '<a target="_blank" href="User.pl?name=' + value + '">' + value + '</a>';
                       }*/},
               	{id: "page", name: "Page", field: "page", minWidth: 90, maxWidth: 100, sortable: true, cssClass: "cell-normal"},
               	{id: "description", name: "Description", field: "description", minWidth: 100, sortable: true, cssClass: "cell-normal",
               		formatter: function ( row, cell, value, columnDef, dataContext ) {
                           return '<span>' + value + '</span>';
                       }},
               	{id: "link", name: "Link", field: "link", minWidth: 100, maxWidth: 250, cssClass: "cell-normal",
               		formatter: function ( row, cell, value, columnDef, dataContext ) {
                           return '<a target="_blank" href="' + value + '">' + value + '</a>';
                       }},
               	{id: "comment", name: "Comments (click to edit)", field: "comment", minWidth: 100, sortable: true, cssClass: "cell-normal",
               		editor: Slick.Editors.Text, validator: requiredFieldValidator}
               ];
    
    window.hist = new coge.Grid('#history', hist_options, hist_columns);
    get_history();
    
    hist.grid.onCellChange.subscribe(function (e, args) {
		$.ajax({
			data: {
				jquery_ajax: 1,
				fname: 'update_comment',
				log_id: args.item.id,
				comment: args.item.comment
			},
			success: function() {
				hist.dataView.updateItem(args.item.id, args.item);
			}
		});
	});
    
    hist.grid.onSort.subscribe(function (e, args) {
		sortcol = args.sortCol.field;
		hist.dataView.sort(args.sortAsc);
	});

	// Wire up model events to drive the grid
	hist.dataView.onRowCountChanged.subscribe(function (e, args) {
		hist.grid.updateRowCount();
		hist.grid.render();
		if ($("#history").is(":not(visible)")) {
			$("#history").slideDown();
		}
	});
	hist.dataView.onRowsChanged.subscribe(function (e, args) {
		hist.grid.invalidateRows(args.rows);
		hist.grid.render();
	});

	// Wire up the show selector to apply the filter to the model
	$("#hist_show_select,#hist_search_type").change(function (e) {
		updateHistFilter();
	});

	// Wire up the search textbox to apply the filter to the model
	$("#hist_search_input").keyup(function (e) {
		Slick.GlobalEditorLock.cancelCurrentEdit();

		if (e.which == 27) { // Clear on Esc
			this.value = "";
		}

		updateHistFilter();
	});
}

//initialize Reports tab
function init_reports() {
	reports_grid = new DataGrid({
		elementId: "reports",
		filter: "none",
		selection: "total",
	});
}

function change_tab(tab) {
	current_tab = tab;
	//console.log(tab);
}

//Tab 1 Search
function search_stuff (search_term) {
	if(search_term.length > 2) {
		$("#loading").show();
		timestamps['search_stuff'] = new Date().getTime();
		$.ajax({
			//type: "POST",
			//dataType: "json",
			//contentType: "application/json; charset=utf8",
			data: {
				fname: 'search_stuff',
				search_term: search_term,
				timestamp: timestamps['search_stuff']
			},
			success : function(data) {
				//console.log(data);
				//$('#loading').hide();
				var obj = jQuery.parseJSON(data);
				
				if (obj && obj.items && obj.timestamp != timestamps['search_stuff']) {
					return;
				}

				var userCounter = 0, orgCounter = 0, genCounter = 0, expCounter = 0, noteCounter = 0, usrgroupCounter = 0;
				var userList = "", orgList = "", genList = "", expList = "", noteList = "", usrgroupList = "";

				for (var i = 0; i < obj.items.length; i++) {
					if (obj.items[i].type == "user") {
						userList = userList + "<tr><td><span>";
						userList = userList + (obj.items[i].label) + " (ID: " + (obj.items[i].id) + ") ";
						userList = userList + "<button onclick=\"search_user(" + (obj.items[i].id) + ",'user')\">Show Data</button>";
						userList = userList + "</span></td></tr>";
						userCounter++;
					}

					if (obj.items[i].type == "organism") {
						orgList = orgList + "<tr><td><span title='" + obj.items[i].description + "'>";
						orgList = orgList + (obj.items[i].label) + " (ID: " + (obj.items[i].id) + ")";
						orgList = orgList + " <a href=\"OrganismView.pl?oid=" + (obj.items[i].id) + "\">Info </a>";
						orgList = orgList + "</span></td></tr>";
						orgCounter++;
					}
	
					if (obj.items[i].type == "genome") {
						genList = genList + "<tr><td><span onclick=\"modify_item(" + obj.items[i].id + ", 'Genome', 'restrict');\"";
						if (obj.items[i].restricted == 1) {
							genList = genList + " class=\"link ui-icon ui-icon-locked\"></span>";
						} else {
							genList = genList + " class=\"link ui-icon ui-icon-unlocked\"></span>";
						}
						genList = genList + "<span onclick=\"modify_item(" + obj.items[i].id + ", 'Genome', 'delete');\"";
						genList = genList + " class=\"link ui-icon ui-icon-trash\"></span>";
						if (obj.items[i].deleted == 1) {
							genList = genList + "<span style=\"color: red\">";
						} else {
							genList = genList + "<span>";
						}
						genList = genList + (obj.items[i].label) + " <a href=\"GenomeInfo.pl?gid=" + (obj.items[i].id) + "\">Info </a>";
						genList = genList + "<button onclick='share_dialog(" + obj.items[i].id + ", 2, " + obj.items[i].restricted + ")'>Edit Access</button>";
						genList = genList + "</span></td></tr>";
						genCounter++;
					}
	
					if (obj.items[i].type == "experiment") {
						expList = expList + "<tr><td><span onclick=\"modify_item(" + obj.items[i].id + ", 'Experiment', 'restrict');\"";
						if (obj.items[i].restricted == 1) {
							expList = expList + " class=\"link ui-icon ui-icon-locked\"></span>";
						} else {
							expList = expList + " class=\"link ui-icon ui-icon-unlocked\"></span>";
						}
						expList = expList + "<span onclick=\"modify_item(" + obj.items[i].id + ", 'Experiment', 'delete');\"";
						expList = expList + " class=\"link ui-icon ui-icon-trash\"></span>";
						if (obj.items[i].deleted == 1) {
							expList = expList + "<span style=\"color: red\">";
						} else {
							expList = expList + "<span>";
						}
						expList = expList + (obj.items[i].label) + " (ID: " + (obj.items[i].id) + ") <a href=\"ExperimentView.pl?eid=" + (obj.items[i].id) + "\">Info </a>";
						expList = expList + "<button onclick='share_dialog(" + obj.items[i].id + ", 3, " + obj.items[i].restricted + ")'>Edit Access</button>";
						expList = expList + "</span></td></tr>";
						expCounter++;
					}
	
					if (obj.items[i].type == "notebook") {
						noteList = noteList + "<tr><td><span onclick=\"modify_item(" + obj.items[i].id + ", 'List', 'restrict');\"";
						if (obj.items[i].restricted == 1) {
							noteList = noteList + " class=\"link ui-icon ui-icon-locked\"></span>";
						} else {
							noteList = noteList + " class=\"link ui-icon ui-icon-unlocked\"></span>";
						}
						noteList = noteList + "<span onclick=\"modify_item(" + obj.items[i].id + ", 'List', 'delete');\"";
						noteList = noteList + " class=\"link ui-icon ui-icon-trash\"></span>";
						if (obj.items[i].deleted == 1) {
							noteList = noteList + "<span style=\"color: red\">";
						} else {
							noteList = noteList + "<span>";
						}
						noteList = noteList + (obj.items[i].label) + " (ID: " + (obj.items[i].id) + ") <a href=\"NotebookView.pl?lid=" + (obj.items[i].id) + "\">Info </a>";
						noteList = noteList + "<button onclick='share_dialog(" + obj.items[i].id + ", 1 , " + obj.items[i].restricted + ")'>Edit Access</button>";
						noteList = noteList + "</span></td></tr>";
						noteCounter++;
					}
							
					if (obj.items[i].type == "user_group") {
						usrgroupList = usrgroupList + "<tr><td>";
						if (obj.items[i].deleted == 1) {
							usrgroupList = usrgroupList + "<span style=\"color: red\">";
						} else {
							usrgroupList = usrgroupList + "<span>";
						}
						usrgroupList = usrgroupList + (obj.items[i].label) + " (ID: " + (obj.items[i].id) + ") ";
						usrgroupList = usrgroupList + "<button onclick=\"search_user(" + (obj.items[i].id) + ",'group')\">Show Data</button>";
						usrgroupList = usrgroupList + "</span></td></tr>";
						usrgroupCounter++;
					}
				}
				
	
				//Populate the html with the results
				$("#loading").show();
				$(".result").fadeIn( 'fast');
				
				//user
				if(userCounter > 0) {
					$('#user').show();
					$('#userCount').html("Users: " + userCounter);
					$('#userList').html(userList);
					if(userCounter <= 10) {
						$( "#userList" ).show();
						//$( "#userArrow" ).find('img').toggle();
						$("#userArrow").find('img').attr("src", "picts/arrow-down-icon.png");
					} else {
						$( "#userList" ).hide();
						$("#userArrow").find('img').attr("src", "picts/arrow-right-icon.png");
					}
				} else {
					$('#user').hide();
				}
				
				//organism
				if(orgCounter > 0) {
					$('#organism').show();
					$('#orgCount').html("Organisms: " + orgCounter);
					$('#orgList').html(orgList);
					if(orgCounter <= 10) {
						$( "#orgList" ).show();
						//$( "#orgArrow" ).find('img').toggle();
						$( "#orgArrow" ).find('img').attr("src", "picts/arrow-down-icon.png");
					} else {
						$( "#orgList" ).hide();
						$("#orgArrow").find('img').attr("src", "picts/arrow-right-icon.png");
					}
				} else {
					$('#organism').hide();
				}
				
				//genome
				if(genCounter > 0) {
					$('#genome').show();
					$('#genCount').html("Genomes: " + genCounter);
					$('#genList').html(genList);
					if(genCounter <= 10) {
						$( "#genList" ).show();
						//$( "#genArrow" ).find('img').toggle();
						$( "#genArrow" ).find('img').attr("src", "picts/arrow-down-icon.png");
					} else {
						$( "#genList" ).hide();
						$("#genArrow").find('img').attr("src", "picts/arrow-right-icon.png");
					}
				} else {
					$('#genome').hide();
				}
				
				//experiment
				if(expCounter > 0) {
					$('#experiment').show();
					$('#expCount').html("Experiments: " + expCounter);
					$('#expList').html(expList);
					if(expCounter <= 10) {
						$( "#expList" ).show();
						//$( "#expArrow" ).find('img').toggle();
						$( "#expArrow" ).find('img').attr("src", "picts/arrow-down-icon.png");
					} else {
						$( "#expList" ).hide();
						$("#expArrow").find('img').attr("src", "picts/arrow-right-icon.png");
					}
				} else {
					$('#experiment').hide();
				}
				
				//notebook
				if(noteCounter > 0) {
					$('#notebook').show();
					$('#noteCount').html("Notebooks: " + noteCounter);
					$('#noteList').html(noteList);
					if(noteCounter <= 10) {
						$( "#noteList" ).show();
						//$( "#noteArrow" ).find('img').toggle();
						$( "#noteArrow" ).find('img').attr("src", "picts/arrow-down-icon.png");
					} else {
						$( "#noteList" ).hide();
						$("#noteArrow").find('img').attr("src", "picts/arrow-right-icon.png");
					}
				} else {
					$('#notebook').hide();
				}
				
				//user group
				if(usrgroupCounter > 0) {
					$('#user_group').show();
					$('#usrgroupCount').html("User Groups: " + usrgroupCounter);
					$('#usrgroupList').html(usrgroupList);
					if(usrgroupCounter <= 10) {
						$( "#usrgroupList" ).show();
						//$( "#usrGArrow" ).find('img').toggle();
						$("#usrGArrow").find('img').attr("src", "picts/arrow-down-icon.png");
					} else {
						$( "#usrgroupList" ).hide();
						$("#usrGArrow").find('img').attr("src", "picts/arrow-right-icon.png");
					}
				} else {
					$('#user_group').hide();
				}
				
				$("#loading").hide();
			},
		});
		previous_search = search_term;
	}
}

function show_table(id) {
	$(id).fadeToggle('fast');	
}

function toggle_arrow(id) {
	//$(id).find('img').toggle();
	if( $(id).find('img').attr('src') == "picts/arrow-right-icon.png" ) {
        	$(id).find('img').attr("src", "picts/arrow-down-icon.png");
        } else {
		$(id).find('img').attr("src", "picts/arrow-right-icon.png");
	}
}

function toggle_master() {
	if (current_page == 0) {
		current_page = 1;
		$('#master').fadeToggle('slow', function() {
			$('#userInfo').fadeToggle('slow');
		});
	} else {
		current_page = 0;
		$('#userInfo').fadeToggle('slow', function() {
			$('#master').fadeToggle('slow');
		});
	}

}

function open_dialog() {
	$("#user_dialog").dialog("open");
}

function edit_access(id, type) {
	//console.log("edit_access called " + id + " " + type);

        //var selected = get_selected_items();
	var selected = "content_" + id + "_" + type;

        var target_item = $('#share_input').data('select_id');
        var role_id = $('#share_role_select').val();
        if (target_item && selected.length) {
                var item_list = selected; //.map(function(){return this.parentNode.id;}).get().join(',');
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
                                        //refresh_info_cache();
                                }
                        }
                });
        }
}

function remove_items_from_user_or_group(target_item, id, type) {
	var selected = "content_" + id + "_" + type; //get_selected_items();
	if (target_item && selected.length) {
		var item_list = selected; //.map(function(){return this.parentNode.id;}).get().join(',');
		$.ajax({
			data: {
				fname: 'remove_items_from_user_or_group',
				target_item: target_item,
				item_list: item_list,
			},
			success : function(data) {
				if (data) {
					$('#share_dialog').html(data);
					//refresh_info_cache();
				}
			}
		});
	}
}

//indicates the previous user searched, used by search_user and refresh_data
var previous_user = 0;
var previous_type = 0;

// indicates the part of the Admin page currently being displayed. 0 --> "master", 1 --> "user info"
var current_page = 0;

function search_user(userID, search_type) {
	if (current_page == 0) {	
		toggle_master();
		//current_page=1;
	}
	if(previous_user != userID) {
		//$('#userResults').hide();
		//$('#userResults').html("Loading...");
		user_info(userID, search_type);
	}
	previous_user = userID;
	previous_type = search_type;
}

function refresh_data() {
	if (current_page == 0) {
		//$('#masterTable').html("Loading...");
		search_stuff(previous_search);
	} else {
		//$('#userResults').html("Loading...");
		user_info(previous_user, previous_type);
	}
}

function user_info(userID, search_type) {

	var search_term = userID;
	$("#loading2").show();
	timestamps['user_info'] = new Date().getTime();
	$.ajax({
		data: {
			fname: 'user_info',
			search_term: search_term,
			search_type: search_type,
			timestamp: timestamps['user_info']
		},
		success : function(data) {
			//console.log("Ajax success");
			var obj = jQuery.parseJSON(data);
        	//console.log(obj.items);
        		
        	var htmlBlock = "";

        	//for each user
        	for (var i = 0; i < obj.items.length; i++) {

        		var genList = "", expList = "", noteList = "", userList = "";
        		var genCounter = 0, expCounter = 0, noteCounter = 0, userCounter = 0;
				
        		//for each object belonging to that user, populate tables
        		for (var j = 0; j < obj.items[i].result.length; j++) {
        			if (obj && obj.items && obj.timestamp == timestamps['user_info']) {
						
        				var current = obj.items[i].result[j];
	
        				if (current.type == "genome") {
        					genList = genList + "<tr><td><span onclick=\"modify_item(" + current.id + ", 'Genome', 'restrict');\"";
    						if (current.restricted == 1) {
    							genList = genList + " class=\"link ui-icon ui-icon-locked\"></span>";
    						} else {
    							genList = genList + " class=\"link ui-icon ui-icon-unlocked\"></span>";
    						}
        					genList = genList + "<span onclick=\"modify_item(" + current.id + ", 'Genome', 'delete');\"";
        					genList = genList + " class=\"link ui-icon ui-icon-trash\"></span>";
	                                               	
        					if (current.role == 2) {
        						genList = genList + "<span style='color: green'>Owner: </span>";
        					} else if (current.role == 3) {
        						genList = genList + "<span style='color: red'>Editor: </span>";
        					} else if (current.role == 4) {
        						genList = genList + "<span style='color: blue'>Reader: </span>";
        					}

        					if (current.deleted == 1) {
                               	genList = genList + "<span style=\"color: red\">";
                            } else {
                            	genList = genList + "<span>";
                            }
        					genList = genList + (current.label) + " <a href=\"GenomeInfo.pl?gid=" + (current.id) + "\">Info </a>";
        					genList = genList + "<button onclick='share_dialog(" + current.id + ", 2, " + current.restricted + ")'>Edit Access</button>";
        					genList = genList + "</span></td></tr>";
        					genCounter++;
        				}

        				if (current.type == "experiment") {
        					expList = expList + "<tr><td><span onclick=\"modify_item(" + current.id + ", 'Experiment', 'restrict');\"";
    						if (current.restricted == 1) {
    							expList = expList + " class=\"link ui-icon ui-icon-locked\"></span>";
    						} else {
    							expList = expList + " class=\"link ui-icon ui-icon-unlocked\"></span>";
    						}
        					expList = expList + "<span onclick=\"modify_item(" + current.id + ", 'Experiment', 'delete');\"";
        					expList = expList + " class=\"link ui-icon ui-icon-trash\"></span>";
        					
        					if (current.role == 2) {
        						expList = expList + "<span style='color: green'>Owner: </span>";
        					} else if (current.role == 3) {
        						expList = expList + "<span style='color: red'>Editor: </span>";
        					} else if (current.role == 4) {
        						expList = expList + "<span style='color: blue'>Reader: </span>";
        					}

        					if (current.deleted == 1) {
        						expList = expList + "<span style=\"color: red\">";
        					} else {
        						expList = expList + "<span>";
        					}
        					expList = expList + (current.label) + " (ID: " + (current.id) + ") <a href=\"ExperimentView.pl?eid=" + (current.id) + "\">Info </a>";
        					expList = expList + "<button onclick='share_dialog(" + current.id + ", 3, " + current.restricted + ")'>Edit Access</button>";
        					expList = expList + "</span></td></tr>";
        					expCounter++;
        				}
						
        				if (current.type == "notebook") {
        					noteList = noteList + "<tr><td><span onclick=\"modify_item(" + current.id + ", 'List', 'restrict');\"";
    						if (current.restricted == 1) {
    							noteList = noteList + " class=\"link ui-icon ui-icon-locked\"></span>";
    						} else {
    							noteList = noteList + " class=\"link ui-icon ui-icon-unlocked\"></span>";
    						}
        					noteList = noteList + "<span onclick=\"modify_item(" + current.id + ", 'List', 'delete');\"";
        					noteList = noteList + " class=\"link ui-icon ui-icon-trash\"></span>";

        					if (current.role == 2) {
        						noteList = noteList + "<span style='color: green'>Owner: </span>";
        					} else if (current.role == 3) {
        						noteList = noteList + "<span style='color: red'>Editor: </span>";
        					} else if (current.role == 4) {
        						noteList = noteList + "<span style='color: blue'>Reader: </span>";
        					}

        					if (current.deleted == 1) {
        						noteList = noteList + "<span style=\"color: red\">";
        					} else {
        						noteList = noteList + "<span>";
        					}
        					noteList = noteList + (current.label) + " (ID: " + (current.id) + ") <a href=\"NotebookView.pl?lid=" + (current.id) + "\">Info </a>";
        					noteList = noteList + "<button onclick='share_dialog(" + current.id + ", , " + current.restricted + ")'>Edit Access</button>";
        					noteList = noteList + "</span></td></tr>";
        					noteCounter++;
        				}
						
        				if (current.type == "user") {
        					userList = userList + "<tr><td><span>";
        					userList = userList + (current.label) + " (ID: " + (current.id) + ") ";
        					userList = userList + "<button onclick=\"search_user(" + current.id + ",'user')\">Search</button>";
        					userList = userList + "</span></td></tr>";
        					userCounter++;
        				}

        			}
        		} //end of single user loop

        		var genBlock = "", noteBlock = "", expBlock = "", userBlock = "", nameBlock = "";
				
        		nameBlock = nameBlock + "<div style=\"padding-top:10px;\">"
        		if (i == 0) {
        			nameBlock = nameBlock + "<img src='picts/user-icon.png' width='15' height='15'><span> ";
        		} else {
        			nameBlock = nameBlock + "<img src='picts/group-icon.png' width='15' height='15'><span> ";
        		}
        		nameBlock = nameBlock + obj.items[i].user + " (ID: " + obj.items[i].user_id + "): </span>";
        		
        		if (search_type =='group') {
        			nameBlock = nameBlock + "<button onclick='group_dialog(" + obj.items[i].user_id + ", 6 )'>Edit Group</button></div>";
        		} else {
        			nameBlock = nameBlock + "</div>";
        		}
				

        		if (genCounter > 0) {
        			genBlock = genBlock + "<div style=\"padding-top:10px;padding-left:30px;\">";
        			genBlock = genBlock + "<span id='genCount" + i + "' class='coge-table-header' style='color:119911;' onclick=\"toggle_arrow('#genArrow" + i + "');show_table('#genList" + i + "')\">";	
        			genBlock = genBlock + "Genomes: " + genCounter + " </span>";
        			genBlock = genBlock + "<div id=\"genArrow" + i + "\" onclick=\"toggle_arrow('#genArrow" + i + "');show_table('#genList" + i + "')\" style='display:inline;'>";
        			genBlock = genBlock + "<img src=\"picts/arrow-right-icon.png\" class=\"link\" style=\"width:10px;height:10px;\"/></div>";
        			genBlock = genBlock + "<table cellspacing=\"5\" class=\"hidden\" id='genList" + i + "' style=\"border-top:0px solid green; padding-left:20px; padding-bottom:10px;\">";
        			genBlock = genBlock + genList + "</table></div>";
        		}

        		if (noteCounter > 0) {
        			noteBlock = noteBlock + "<div style=\"padding-top:10px;padding-left:30px;\">";
        			noteBlock = noteBlock + "<span id='noteCount" + i + "' class='coge-table-header' style='color:119911;' onclick=\"toggle_arrow('#noteArrow" + i + "');show_table('#noteList" + i + "')\">";
        			noteBlock = noteBlock +	"Notebooks: " + noteCounter + " </span>";
        			noteBlock = noteBlock + "<div id=\"noteArrow" + i + "\" onclick=\"toggle_arrow('#noteArrow" + i + "');show_table('#noteList" + i + "')\" style='display:inline;'>";
        			noteBlock = noteBlock + "<img src=\"picts/arrow-right-icon.png\" class=\"link\" style=\"width:10px;height:10px;\"/></div>";
        			noteBlock = noteBlock + "<table cellspacing=\"5\" class=\"hidden\" id='noteList" + i + "' style=\"border-top:0px solid green; padding-left:20px; padding-bottom:10px;\">"; 
        			noteBlock = noteBlock + noteList + "</table></div>";
        		}

        		if (expCounter > 0) {
        			expBlock = expBlock + "<div style=\"padding-top:10px;padding-left:30px;\">";
        			expBlock = expBlock + "<span id='expCount" + i + "' class='coge-table-header' style='color:119911;' onclick=\"toggle_arrow('#expArrow" + i + "');show_table('#expList" + i + "')\">";
        			expBlock = expBlock + "Experiments: " + expCounter + " </span>";
        			expBlock = expBlock + "<div id=\"expArrow" + i + "\" onclick=\"toggle_arrow('#expArrow" + i + "');show_table('#expList" + i + "')\" style='display:inline;'>";
        			expBlock = expBlock + "<img src=\"picts/arrow-right-icon.png\" class=\"link\" style=\"width:10px;height:10px;\"/></div>";
        			expBlock = expBlock + "<table cellspacing=\"5\" class=\"hidden\" id='expList" + i + "' style=\"border-top:0px solid green; padding-left:20px; padding-bottom:10px;\">";
        			expBlock = expBlock + expList + "</table></div>";
        		}

        		if (userCounter > 0) {
        			userBlock = userBlock + "<div style=\"padding-top:10px;padding-left:30px;\">";
        			userBlock = userBlock + "<span id='userCount" + i + "' class='coge-table-header' style='color:119911;' onclick=\"toggle_arrow('#userArrow" + i + "');show_table('#userList" + i + "')\">";
        			userBlock = userBlock + "Users: " + userCounter + " </span>";
        			userBlock = userBlock + "<div id=\"userArrow" + i + "\" onclick=\"toggle_arrow('#userArrow" + i + "');show_table('#userList" + i + "')\" style='display:inline;'>";
        			userBlock = userBlock + "<img src=\"picts/arrow-right-icon.png\" class=\"link\" style=\"width:10px;height:10px;\"/></div>";
        			userBlock = userBlock + "<table cellspacing=\"5\" class=\"hidden\" id='userList" + i + "' style=\"border-top:0px solid green; padding-left:20px; padding-bottom:10px;\">"; 
        			userBlock = userBlock + userList + "</table></div>";
        		}
					
        		htmlBlock = htmlBlock + nameBlock + genBlock + expBlock + noteBlock + userBlock;

        		//genome for user_info
				
        	} //end of all users loop

        	$('#userResults').html(htmlBlock);
        	//$('#userResults').show();

        	if(search_type == 'group') {
        		if(expCounter <= 10 && expCounter > 0) {
        			$( "#expList0" ).show();
        			$( "#expArrow0" ).find('img').attr("src", "picts/arrow-down-icon.png");
        		} else {
        			$( "#expList0" ).hide();
        			$("#expArrow0" ).find('img').attr("src", "picts/arrow-right-icon.png");
        		}

        		if(genCounter <= 10 && genCounter > 0) {
        			$( "#genList0" ).show();
        			$( "#genArrow0" ).find('img').attr("src", "picts/arrow-down-icon.png");
        		} else {
        			$( "#genList0" ).hide();
        			$("#genArrow0" ).find('img').attr("src", "picts/arrow-right-icon.png");
        		}

        		if(noteCounter <= 10 && noteCounter > 0) {
        			$( "#noteList0" ).show();
        			$( "#noteArrow0" ).find('img').attr("src", "picts/arrow-down-icon.png");
        		} else {
        			$( "#noteList0" ).hide();
        			$("#noteArrow0" ).find('img').attr("src", "picts/arrow-right-icon.png");
        		}

        		if(userCounter <= 10 && userCounter > 0) {
        			$( "#userList0" ).show();
        			$( "#userArrow0" ).find('img').attr("src", "picts/arrow-down-icon.png");
        		} else {
        			$( "#userList0" ).hide();
        			$("#userArrow0" ).find('img').attr("src", "picts/arrow-right-icon.png");
        		}
        	}
        	$("#loading2").hide();
        }
	});
}

function share_dialog(id, type, restricted) {
	console.log(type);
	var item_list = "content_" + id + "_" + type;  //selected.map(function(){return this.parentNode.id;}).get().join(',');
	var type_string;
	switch(type) {
		case 1: type_string = 'List';
			break;
		case 2: type_string = 'Genome';
			break;
		case 3: type_string = 'Experiment';
			break;
	}
	$.ajax({
		data: {
			fname: 'get_share_dialog',
			item_list: item_list,
		},
		success : function(data) {			
			$('#share_dialog').html(data).dialog({width:500}).dialog('open');
			if(restricted == 0) {
				$('#restrict_checkbox').prop('checked', true);
			} else if(restricted == 1) {
				$('#restrict_checkbox').prop('checked', false);
			}
			$('#restrict_checkbox').change(function(e) {
				modify_item(id, type_string, 'restrict');
			});
		}
	});
}

function search_share () {
	var search_term = $('#share_input').attr('value');

	//$("#wait_notebook").animate({opacity:1});
	timestamps['search_share'] = new Date().getTime();

	$.ajax({
		data: {
			fname: 'search_share',
			search_term: search_term,
			timestamp: timestamps['search_share']
		},
		success : function(data) {
			var obj = jQuery.parseJSON(data);
			if (obj && obj.timestamp == timestamps['search_share'] && obj.items) {
				//console.log(obj.items);
				$("#share_input").autocomplete({source: obj.items}).autocomplete("search");
				//$("#wait_notebook").animate({opacity:0});
			}
		},
	});
}

function search_group () { // FIXME dup of above routine but for group dialog
	var search_term = $('#group_input').attr('value');

	timestamps['search_group'] = new Date().getTime();
	$.ajax({
		data: {
			fname: 'search_share',
			search_term: search_term,
			timestamp: timestamps['search_group']
		},
		success : function(data) {
			var obj = jQuery.parseJSON(data);
			if (obj && obj.timestamp == timestamps['search_group'] && obj.items) {
				$("#group_input").autocomplete({source: obj.items}).autocomplete("search");
			}
		},
	});
}

function group_dialog(id, type) {
	var item_list = "content_" + id + "_" + type;
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

function change_group_role(id, type) {
	var selected = "content_" + id + "_" + type; //get_selected_items();
	var role_id = $('#group_role_select').val();
	if (role_id && selected.length) {
		var target_items = selected; //.map(function(){return this.parentNode.id;}).get().join(',');
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

function add_users_to_group(id, type) {
	var selected = "content_" + id + "_" + type; //get_selected_items();
	var new_item = $('#group_input').data('select_id');
	if (new_item && selected.length) {
		var target_items = selected; //.map(function(){return this.parentNode.id;}).get().join(',');
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

function remove_user_from_group(user_id, id, type) {
	var selected = "content_" + id + "_" + type; //get_selected_items();
	if (user_id && selected.length) {
		var target_items = selected; //.map(function(){return this.parentNode.id;}).get().join(',');
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


function modify_item(id, type, modification) {
	console.log(type);
	$.ajax({
		data: {
			fname: 'modify_item',
			id: id,
			modification: modification,
			type: type,
		},
		success : function(val) {
			refresh_data();
		},
	});
}

function wait_to_search(search_func, search_term) {
	//console.log(search_term);
	pageObj.search_term = search_term;

	if (pageObj.time) {
		clearTimeout(pageObj.time);
	}

	// FIXME: could generalize by passing select id instead of separate search_* functions
	pageObj.time = setTimeout(
		function() {
			search_func(pageObj.search_term);
			//console.log("request sent");
		},
		1000
	);
}

//The following javascript deals with Tab2, the Jobs tab
function get_jobs() {
	cancel_update("jobs");
	$.ajax({
		dataType: 'json',
	    data: {
	        jquery_ajax: 1,
	        fname: 'get_jobs',
	        time_range: 0,
	        running_only: running_only,
	    },
	    success: function(data) {
	    	//console.log(data.jobs);
	        jobs.load(data.jobs);
	        entries = data.jobs.length;
	        $("#filter_busy").hide();
	        update_filter();
	        jobs_init = true;
	    },
	    complete: function(data) {
	    	schedule_update(5000);
	    }
	});
}

function update_filter() {
    jobs.dataView.setFilterArgs({
        show: $('#show_select').val(),
        searchType: $('#job_search_type').val(),
        searchString: $('#job_search_bar').val().toLowerCase()
    });

    jobs.filter();
    $('#job_filter_count').html('Showing ' + jobs.dataView.getLength() + ' of ' + entries + ' results');
}

function toggle_job_updater() {
	jobs_updating = !jobs_updating;
	if (jobs_updating) {
		schedule_update(5000);
	}
}

function toggle_running() {
	if(running_only == 0) {
		running_only = 1;
	} else {
		running_only = 0;
	}
	get_jobs();
}

function schedule_update(delay) {
	var idleTime = new Date().getTime() - timestamps['idle'];
	if (idleTime < IDLE_TIME && delay !== undefined) {
		if(current_tab == 1 && jobs_updating) {
			console.log("Updating jobs");
			cancel_update("jobs");
			jobs_timers['update'] = window.setTimeout(
				function() { get_jobs(); },
				delay
			);
		}
		if(current_tab == 2 && hist_updating && hist_init) {
			console.log("Updating hist");
			cancel_update("hist");
			hist_timers['update'] = window.setTimeout(
				function() { update_history(); },
				delay
			);
		}
		return;
	}	
}

function cancel_update(page) {
	if(page == "jobs") {
		clearTimeout(jobs_timers['update']);
	}
	if(page == "hist") {
		clearTimeout(hist_timers['update']);
	}
}

function cancel_job() {
    submit_task("cancel_job", function(row) {
        return row.status.toLowerCase() === 'running'
    });
}

function restart_job() {
    submit_task("restart_job", function(row) {
        return row.status.toLowerCase() === 'cancelled' ||
               row.status.toLowerCase() === 'stopped';
    });
}

function submit_task(task, predicate) {
    var selectedIndexes = window.jobs.grid.getSelectedRows();

    var selectedRows = selectedIndexes.map(function(item) {
        return window.jobs.dataView.getItem(item);
    });

    var validRows = selectedRows.filter(predicate);

    // No rows were valid
    if (!validRows.length) return;

    jQuery.each(validRows, function(index,row) {
        var argument_list =  {
            fname: task,
            job: row.workflow_id,
        };

        $.ajax({
            type: "GET",
            dataType: "json",
            data: argument_list,
            success: function(data) {
                if (data.status) {
                    row.status = data.status;
                    window.jobs.dataView.updateItem(row.id, row);
                }
            }
        });
    });

    // Deselect all rows
    window.jobs.grid.setSelectedRows([]);
}


//The following Javascript deals with Tab3, the History page 
function get_history() {
	$("#loading3").show();
	$.ajax({
		dataType: 'json',
		data: {
			jquery_ajax: 1,
			fname: 'get_history_for_user',
			time_range: 0,
		},
		success : function(data) {
			//console.log(data[0]);
			hist.load(data);
			hist_entries = data.length;
			last_hist_update = data[0].date_time;
			updateHistFilter();
			
			hist_init = true;
			$("#loading3").hide();
		},
	    complete: function(data) {
	    	schedule_update(5000);
	    }
	});
}

function update_history() {
	//console.log(last_hist_update);
	$.ajax({
		dataType: 'json',
		data: {
			jquery_ajax: 1,
			fname: 'update_history',
			timestamp: last_hist_update,
			time_range: 0,
		},
		success: function(data) {
			console.log(data);
			if(data[0]) {
				hist.insert(data);
				last_hist_update = data[0].date_time;
			}
		},
		complete: function(data) {
			schedule_update(5000);
	    }
	})
}

function toggle_hist_updater() {
	hist_updating = !hist_updating;
	if (hist_updating && hist_init) {
		schedule_update(5000);
	}
}

function requiredFieldValidator(value) {
	return {valid: true, msg: null};
}

function updateHistFilter() {
	hist.dataView.setFilterArgs({
		show: $('#hist_show_select').val(),
		searchType: $('#hist_search_type').val(),
		searchString: $('#hist_search_input').val().toLowerCase()
	});
    hist.filter();
    $('#hist_filter_count').html('Showing ' + hist.dataView.getLength() + ' of ' + hist_entries + ' results');
}

function toggle_star(img) {
	$.ajax({
		data: {
			jquery_ajax: 1,
			fname: 'toggle_star',
			log_id: img.id,
		},
		success :  function(val) {
			if (val == 0) { $(img).attr({src:"picts/star-hollow.png"}); }
			else { $(img).attr({src:"picts/star-full.png"}); }
		}
	});
}


//Tab 4, the User Graph
var colors = [
        	{ name: 'list',       link: 'NotebookView.pl?nid=',   color: 'Tomato',      show: 1 },
        	{ name: 'genome',     link: 'GenomeInfo.pl?gid=',     color: 'YellowGreen', show: 1 },
        	{ name: 'experiment', link: 'ExperimentView.pl?eid=', color: 'Orchid',      show: 1 },
        	{ name: 'feature',    link: '',                       color: 'orange',      show: 1 },
        	{ name: 'user',       link: '',                       color: 'DeepSkyBlue', show: 1 },
        	{ name: 'group',      link: 'GroupView.pl?ugid=',     color: 'Turquoise',   show: 1 },
];

var w = Math.max(800, $(window).width()-200),
	h = Math.max(800, $(window).height()-300),
	user_force,
	group_force;

function init_graph(selectId) {
	//console.log(selectId);
	if(selectId == 1) {
		$('#group_chart').hide();
		$('#group_legend').hide();
		$('#user_chart').show();
		$('#user_legend').show();
		
		if(!user_graph_init) {	
			user_graph_init = true;
			$("#loading4").show();
			d3.json("?fname=get_user_nodes", function(json) {
				console.log(json);
				var root = json;
				var nodes = flatten(root);
				nodes.forEach(function(d) {
					d._children = d.children;
					d.children = null;
				});
				
				svg = d3.select("#user_chart").append("svg")
					.attr("width", w)
					.attr("height", h);
				
				var link = svg.selectAll(".link");
				var node = svg.selectAll(".node");
				
				var user_filters = [
				                    { name: 'deleted', 		color: 'Red', 		show: 0 },
				                    { name: 'restricted', 	color: 'Grey',	 	show: 0 },
				                    { name: '||||||||',		color: 'White',		show: 1 },
				];
				user_filters.forEach(function(element, index) {
					var item =
						$('<div class="link legend selected">'+user_filters[index].name+'</div>')
							.css('color', user_filters[index].color)
							.css('background-color', '')
							.click(function() {
								$(this).toggleClass('selected');
								if ($(this).hasClass('selected')) {
									$(this).css('color', user_filters[index].color);
									$(this).css('background-color', '');
								}
								else {
									$(this).css('color', 'white');
									$(this).css('background-color', user_filters[index].color);
								}
								user_filters[index].show = !user_filters[index].show;
								
								var nodes = flatten(root);
								nodes.forEach(function(d) {
									return color(d, user_filters);
								});
								user_force.update();
							});

					$('#user_legend')
						.append(item);
				});
				
				colors.forEach(function(element, index) {
					var item =
						$('<div class="link legend selected">'+colors[index].name+'</div>')
							.css('color', 'white')
							.css('background-color', colors[index].color);
	
					$('#user_legend')
						.append(item);
				});
				
				var user_force = new Force(root, link, node, user_filters);
				$("#loading4").hide();
				user_force.update();
			});
		} 
	} else if (selectId == 2) {
		$('#user_chart').hide();
		$('#user_legend').hide();
		$('#group_chart').show();
		$('#group_legend').show();
		
		if(!group_graph_init) {
			group_graph_init = true;
			$("#loading4").show();
			d3.json("?fname=get_group_nodes", function(json) {
				console.log(json);
				var root = json;
				var nodes = flatten(root);
				nodes.forEach(function(d) {
					d._children = d.children;
					d.children = null;
				});
				
				svg = d3.select("#group_chart").append("svg")
					.attr("width", w)
					.attr("height", h);
			
				var link = svg.selectAll(".link");
				var node = svg.selectAll(".node");
				
				var group_filters = [
				                     { name: 'deleted', 		color: 'Red', 		show: 0 },
				                     { name: 'restricted', 	color: 'Grey',	 	show: 0 },
				                     { name: '||||||||',		color: 'White',		show: 1 },
				];
				group_filters.forEach(function(element, index) {
					var item =
						$('<div class="link legend selected">'+group_filters[index].name+'</div>')
							.css('color', group_filters[index].color)
							.css('background-color', 'white')
							.click(function() {
								$(this).toggleClass('selected');
								if ($(this).hasClass('selected')) {
									$(this).css('color', group_filters[index].color);
									$(this).css('background-color', '');
								}
								else {
									$(this).css('color', 'white');
									$(this).css('background-color', group_filters[index].color);
								}
								group_filters[index].show = !group_filters[index].show;
								
								var nodes = flatten(root);
								nodes.forEach(function(d) {
									return color(d, group_filters);
								});
								group_force.update();
							});

					$('#group_legend')
						.append(item);
				});
				
				colors.forEach(function(element, index) {
					var item =
						$('<div class="link legend selected">'+colors[index].name+'</div>')
							.css('color', 'white')
							.css('background-color', colors[index].color);
	
					$('#group_legend')
						.append(item);
				});
				
				var group_force = new Force(root, link, node, group_filters);
				$("#loading4").hide();
				group_force.update();
			});
		}
	}
}

var Force = function(root, link, node, filters) {
	var self = this;
	this.link = link;
	this.node = node;
	this.root = root;
	this.filters = filters;
	this.force = d3.layout.force()
		.size([w, h])
		.on("tick", function() {
			self.tick.call(self);
		})
		.gravity(1);
}

Force.prototype.update = function() {
	var self = this;
	var nodes = flatten(this.root),
		links = d3.layout.tree().links(nodes);
	
	// Restart the force layout.
	this.force
 		.nodes(nodes)
 		.links(links)
 		.start();

	// Update the links…
	this.link = this.link.data(links, function(d) { return d.target.id; });

	// Exit any old links.
	this.link.exit().remove();

	// Enter any new links.
	this.link.enter().insert("line", ".node")
		.attr("class", "link")
		.attr("x1", function(d) { return d.source.x; })
		.attr("y1", function(d) { return d.source.y; })
		.attr("x2", function(d) { return d.target.x; })
		.attr("y2", function(d) { return d.target.y; });

	// Update the nodes…
	this.node = this.node.data(nodes, function(d) { return d.id; })
		.style("fill", function(d) {
			return color(d, self.filters);
		})
		//.attr("r", function(d) { return Math.sqrt(d.size) / 10 || 4.5; });
		.attr("r", function(d) { return Math.pow(d.size, 1.0/3.0)/3 || 4.5});

	// Exit any old nodes.
	this.node.exit().remove();

	// Enter any new nodes.
	this.node.enter().append("circle")
		.attr("class", "node")
		.attr("cx", function(d) { return d.x; })
		.attr("cy", function(d) { return d.y; })
		//.attr("r", function(d) { return Math.sqrt(d.size) / 10 || 4.5; })
		.attr("r", function(d) { return Math.pow(d.size, 1.0/3.0)/3 || 4.5})
		.style("fill", function(d) {
			return color(d, self.filters);
		})
		.on("click", function(d) {
			self.click.call(self, d);
		})
		.call(this.force.drag)
		.append("svg:title").text(function(d) { 
			if(d.info) {
				return d.info; 
			} else {
				return d.name;
			}
		});
}

Force.prototype.tick = function() {
	/*this.link.attr("x1", function(d) { return d.source.x; })
		.attr("y1", function(d) { return d.source.y; })
		.attr("x2", function(d) { return d.target.x; })
		.attr("y2", function(d) { return d.target.y; });
	
	this.node.attr("cx", function(d) { return d.x; })
		.attr("cy", function(d) { return d.y; });*/
	
	var self = this;
	self.moveItems.call(self);
}

//Toggle children on click.
Force.prototype.toggle_children = function(d) {
	if (d.children) {
		d._children = d.children;
		d.children = null;
	} else {
		d.children = d._children;
		d._children = null;
	}
	//console.log(d.size);
}

Force.prototype.click = function(d) {
	if (!d3.event.defaultPrevented) {
		//if (d._size) {
		//	var temp = d._size;
		//	d._size = d.size;
		//	d.size = temp;
		//} else {
		//	d._size = d.size;
		//	d.size = 2025;
		//}
		
		this.force.charge(
			function(d, i) {
				if(d.size) {
					//return Math.pow(d.size, 1.0/3.0)*(-12);
					//console.log(d.size + '!');
					return Math.sqrt(d.size)*(-2.25);
				}
			}
		);
		//force.getAttribute("charge");
		this.force.start();
		
		this.toggle_children(d);
		this.update();
	}
}

//Returns a list of all nodes under the root.
function flatten(root) {
	var nodes = [], i = 0;

	function recurse(node) {
		if (node.children) node.children.forEach(recurse);
		if (!node.id) node.id = ++i;
		nodes.push(node);
	}

	recurse(root);
	return nodes;
}

function color(d, new_filters) {
	//console.log(new_filters);
	if(!new_filters[0].show && !new_filters[1].show) {
		//normal filter
		if((d.children || d._children) || (d.type == 2 || d.type == 3 || d.type == 4)) {
			if (d.type) {
				return colors[d.type-1].color;
			}
			return 'white';
		} else {
			return 'black';
		}
	} else {
		if(d.deleted == 1 && d.restricted == 1 && new_filters[0].show && new_filters[1].show) {
			return 'yellow';
		} else if (d.deleted == 1 && new_filters[0].show) {
			return new_filters[0].color;
		} else if (d.restricted == 1 && new_filters[1].show) {
			return new_filters[1].color;
		} else {
			return 'yellowGreen';
		}
	}
}

//Optimization of the Tick function.
Force.prototype.moveItems = (function(){
	var self = this;
	var node = self.node;
	var link = self.link;
    var todoNode = 0;
    var todoLink = 0;
    var MAX_NODES = 240;
    var MAX_LINKS = MAX_NODES/2;
      
    var restart = false;
       
    function moveSomeNodes(){
    	var self = this;
        var n;
        var goal = Math.min(todoNode+MAX_NODES, self.node[0].length);
          
        for(var i=todoNode ; i < goal ; i++){
            n = self.node[0][i];
            n.setAttribute('cx', n.__data__.x);
            n.setAttribute('cy', n.__data__.y);
        }
            
        todoNode = goal;
        requestAnimationFrame(function() {
        	moveSome.call(self);
        });
    }
      
    function moveSomeLinks(){
    	var self = this;
        var l;
        var goal = Math.min(todoLink+MAX_LINKS, self.link[0].length);
           
        for(var i=todoLink ; i < goal ; i++){
            l = self.link[0][i];
            //console.log(l);
            l.setAttribute('x1', l.__data__.source.x);
            l.setAttribute('y1', l.__data__.source.y);
            l.setAttribute('x2', l.__data__.target.x);
            l.setAttribute('y2', l.__data__.target.y);
        }
          
        todoLink = goal;
        requestAnimationFrame(function() {
        	moveSome.call(self);
        });
    }
        
    function moveSome(){
    	var self = this;
        //console.time('moveSome')
        if(todoNode < self.node[0].length) // some more nodes to do
            moveSomeNodes.call(self)
        else{ // nodes are done
            if(todoLink < self.link[0].length) // some more links to do
                moveSomeLinks.call(self)
            else{ // both nodes and links are done
                if(restart){
                    restart = false;
                    todoNode = 0;
                    todoLink = 0;
                    requestAnimationFrame(function() {
                    	moveSome.call(self);
                    });
                }
            }
        }
        //console.timeEnd('moveSome')
    }
        
        
    return function moveItems(){
    	var self = this;
        if(!restart){
            restart = true;
            requestAnimationFrame(function() {
            	moveSome.call(self);
            });
        }
    };
 
})();


//Tab 5, Summary
function DataGrid(params) {
	this.elementId = params.elementId;
	this.filter = params.filter;
	this.selection = params.selection;
	this.data;
	
	this.initialize();
}

$.extend(DataGrid.prototype, {
	initialize: function() {
		var self = this;
		var element = this.elementId;
		var fname;
		$('#' + element).hide();
		$('#' + element + '_loading').show();
		switch (this.selection) {
			case "user":
				fname = 'get_user_table';
				$('#' + element).html('<table id="' + element + '_table" cellpadding="0" cellspacing="0" border="0" class="dt-cell hover compact row-border">'
					+ '<thead><tr><th>User Name</th><th>Notebooks</th><th>Genomes</th><th>Experiments</th><th>Groups</th></tr></thead></table>');
				$('#histogram_button').prop('disabled', false);
				break;
			case "group":
				fname = 'get_group_table';
				$('#' + element).html('<table id="' + element + '_table" cellpadding="0" cellspacing="0" border="0" class="dt-cell hover compact row-border">'
					+ '<thead><tr><th>Group Name</th><th>Notebooks</th><th>Genomes</th><th>Experiments</th><th>Users</th></tr></thead></table>');
				$('#histogram_button').prop('disabled', false);
				break;
			case "total":
				fname = 'get_total_table';
				if (this.filter == "none") {
					$('#' + element).html('<table id="' + element + '_table" cellpadding="0" cellspacing="0" border="0" class="dt-cell hover compact row-border">'
						+ '<thead><tr><th>Notebooks</th><th>Genomes</th><th>Experiments</th><th>Users</th><th>Groups</th></tr></thead></table>');
				} else {
					$('#' + element).html('<table id="' + element + '_table" cellpadding="0" cellspacing="0" border="0" class="dt-cell hover compact row-border">'
						+ '<thead><tr><th>Notebooks</th><th>Genomes</th><th>Experiments</th></tr></thead></table>');
				}
				$('#histogram_button').prop('disabled', true);
				break;
		}
		$.ajax({
			data: {
				fname: fname,
				filter: this.filter,
			},
			success: function(data) {
				$('#' + element + '_loading').hide();
				self.data = JSON.parse(data);
				//console.log(self.data);
				$('#' + element + '_table').dataTable(JSON.parse(data));
				$('#' + element).show();
			}
		});
		
	}
});

function change_report(index) {
	switch (index) {
		case 0: reports_grid.selection = "total";
			break;
		case 1: reports_grid.selection = "user";
			break;
		case 2: reports_grid.selection = "group";
			break;
	}
	reports_grid.initialize();
}

function filter_report(index) {
	switch (index) {
		case 0: reports_grid.filter = "none";
			break;
		case 1: reports_grid.filter = "restricted";
			break;
		case 2: reports_grid.filter = "deleted";
			break;
		case 3: reports_grid.filter = "public";
			break;
	}
	reports_grid.initialize();
}

var Histogram = function(element, json, parent) {
	this.element = element;
	this.json = json;
	this.parent = parent;
	//this.selection = selection;
	
	/*this.values = [];
	this.max_value = 0;
	this.min_value = Number.MAX_SAFE_INTEGER;
	for(var i = 0; i < json.data.length; i++) {
		var total_items = 0;
		for(var j = 0; j < json.data[i].length; j++) {
			if(selection == "total" || j != 0) {
				var num = parseInt(json.data[i][j]);
				total_items += num;
				//console.log(json.data[i][j]);
			}
		}
		this.values.push(total_items);
		if(total_items > this.max_value) {
			this.max_value = total_items;
		}
		if(total_items < this.min_value) {
			this.min_value = total_items;
		}
	}*/
	this.values = this.json.values;
	this.min_value = this.json.min;
	this.max_value = this.json.max;
	
	console.log(this.json);
	
	this.formatCount;
	this.margin = {top: 10, right: 30, bottom: 30, left: 50};
	this.width = 500 - this.margin.left - this.margin.right;
	this.height = 500 - this.margin.top - this.margin.bottom;
	this.x;
	this.data;
	this.y;
	this.xAxis;
	this.yAxis;
	this.svg;
	this.bar;
	this.brush;
	
	this.initialize();
}

$.extend(Histogram.prototype, {
	initialize: function() {
		var self = this;
		$("#" + this.element).html('');
		//Generate a Bates distribution of 10 random variables.
		//this.values = d3.range(1000).map(d3.random.bates(10));
		//console.log(this.values);
	
		// A formatter for counts.
		this.formatCount = d3.format(",.0f");
	
		/*this.margin = {top: 10, right: 30, bottom: 30, left: 50},
		this.width = 500 - this.margin.left - this.margin.right,
		this.height = 500 - this.margin.top - this.margin.bottom;*/
	
		this.x = d3.scale.linear()
		    .domain([this.min_value, this.max_value])
		    .range([0, this.width]);
	
		// Generate a histogram using twenty uniformly-spaced bins.
		this.data = d3.layout.histogram()
		    .bins(this.x.ticks(this.max_value))
		    (this.values);
	
		this.y = d3.scale.log()
		    .domain([1, d3.max(this.data, function(d) { return d.y; })])
		    .range([this.height, 0]);
	
		this.xAxis = d3.svg.axis()
		    .scale(this.x)
		    .orient("bottom");
		
		this.yAxis = d3.svg.axis()
	    	.scale(this.y)
	    	.orient("left");
		
		this.brush = d3.svg.brush()
	    	.x(this.x)
	    	.on("brushend", function() {
	    		self.brushed.call(self);
	    	});
	
		this.svg = d3.select("#" + this.element).append("svg")
	    	.attr("width", this.width + this.margin.left + this.margin.right)
	    	.attr("height", this.height + this.margin.top + this.margin.bottom)
	    	.append("g")
	    	.attr("transform", "translate(" + this.margin.left + "," + this.margin.top + ")");
	
		this.bar = this.svg.selectAll(".bar")
		    .data(self.data)
			.enter().append("g")
		    .attr("class", "bar")
		    .attr("transform", function(d) { 
		    	if(d.y == 0) {
		    		return "translate(" + self.x(d.x) + ", 0)";
		    	}
		    	return "translate(" + self.x(d.x) + "," + self.y(d.y) + ")"; 
		    });
		
		this.bar.append("rect")
		    .attr("x", 1)
		    .attr("width", self.x((self.data[0].dx * 2) - 1))
		    .attr("height", function(d) { 
		    	if(d.y == 0) {
		    		return 0;
		    	}
		    	return self.height - self.y(d.y); 
		    });
	
		/*this.bar.append("text")
		    .attr("dy", ".75em")
		    .attr("y", 6)
		    .attr("x", this.x(this.data[0].dx) / 2)
		    .attr("text-anchor", "top")
		    .text(function(d) { return self.formatCount(d.y); });*/
	
		this.svg.append("g")
		    .attr("class", "x axis")
		    .attr("transform", "translate(0," + this.height + ")")
		    .call(this.xAxis);
		
		this.svg.append("g")
	    	.attr("class", "y axis")
	    	.call(this.yAxis);
		
		this.svg.append("g")
	    	.attr("class", "x brush")
	    	.call(self.brush)
	    	.selectAll("rect")
	    	.attr("y", -6)
	    	.attr("height", this.height + 7);
		
		$('#' + this.element)
			.dialog({
				resizeStop: function( event, ui ) { 
					//self.svg.attr("transform", "scale(" + ui.size.width/ui.originalSize.width + " " + ui.size.height/ui.originalSize.height + ") translate(" + self.margin.left + "," + self.margin.top + ")");

					var widthScale = ui.size.width/ui.originalSize.width;
					var heightScale = ui.size.height/ui.originalSize.height;
					self.width = (self.width + self.margin.left + self.margin.right)*widthScale - self.margin.left - self.margin.right;
					self.height = (self.height + self.margin.top + self.margin.bottom + 30)*heightScale - self.margin.top - self.margin.bottom - 30;
					self.initialize();
				},
			})
			.dialog('open');
	},
	brushed: function() {
		var self = this;
		var extent = this.brush.extent();
		var rangeExtent = [this.x( extent[0] ), this.x( extent[1] ) ];
		var rangeWidth  = rangeExtent[1] - rangeExtent[0];
		console.log(extent[0]);
		console.log(rangeExtent[0]);
		
		var newValues = [];
		for (var i = 0; i < this.values.length; i++) {
			if (this.values[i] >= extent[0] && this.values[i] <= extent[1]) {
				newValues.push(this.values[i]);
			}
		}
		var newJson = {
				min: Math.floor(extent[0]),
				max: Math.ceil(extent[1]),
				values: newValues,
		};
		//console.log(newJson);
		var childHistogram = new Histogram(this.element, newJson, this);
	}
});

function init_histogram(element) {
	var json = reports_grid.data;
	var values = [];
	var max_value = 0;
	var min_value = Number.MAX_SAFE_INTEGER;
	for(var i = 0; i < json.data.length; i++) {
		var total_items = 0;
		for(var j = 0; j < json.data[i].length; j++) {
			if(reports_grid.selection == "total" || j != 0) {
				var num = parseInt(json.data[i][j]);
				total_items += num;
				//console.log(json.data[i][j]);
			}
		}
		values.push(total_items);
		if(total_items > max_value) {
			max_value = total_items;
		}
		if(total_items < min_value) {
			min_value = total_items;
		}
	}
	var histogram = new Histogram(element, {min: min_value, max: max_value, values: values}, null);
	return histogram;
}

function init_taxon_tree(element) {
	if (!tree) {
		$.ajax({
			data: {
				fname: "gen_tree_json",
			},
			success: function(data) {
				console.log(JSON.parse(data));
				tree = new Taxon_tree(JSON.parse(data), "taxonomic_tree");
			}
		});
	}
}

var Taxon_tree = function(json, element) {
	this.margin = {top: 20, right: 120, bottom: 20, left: 120};
	this.width = 1960 - this.margin.right - this.margin.left;
	this.height = 3200 - this.margin.top - this.margin.bottom;
	this.duration = 750;
	this.root = json;
	this.current_root = this.root;
	this.node;
	this.link;
	this.svg;
	this.i = 0;
	
	this.tree = d3.layout.tree()
		.size([this.height, this.width]);
	
	this.diagonal = d3.svg.diagonal()
		.projection(function(d) { return [d.y, d.x]; });
	

	d3.select(self.frameElement).style("height", "800px");
	
	this.initialize(element);
}

$.extend(Taxon_tree.prototype, {
	initialize: function(element) {
		var self = this;
		self.svg = d3.select("#" + element).append("svg")
			.attr("width", self.width + self.margin.right + self.margin.left)
			.attr("height", self.height + self.margin.top + self.margin.bottom)
			.append("g")
			.attr("transform", "translate(" + self.margin.left + "," + self.margin.top + ")");
		
		// Prepopulate the parent nodes of the tree.
		var nodes = self.tree.nodes(self.root).reverse();
		
		self.fix.call(self);
	},
	fix: function() {
		var self = this;
		self.root.x0 = self.height / 2;
		self.root.y0 = 0;

		self.root.children.forEach(function(d) {
			self.collapse.call(self, d);
		});
		self.update.call(self, self.root, self.current_root);
	},
	update: function(source, new_root) {
		var self = this;
		self.current_root = new_root;
		console.log(source.depth);
		
		//Adjust height based on the number of children
		if (self.current_root.children) {
			self.tree.size([(self.current_root.children.length*20) * Math.pow(1.1,source.depth), 1]);
		} else if (self.current_root._children && self.current_root._children.length != 0) {
			self.tree.size([self.current_root._children.length*20, 1]);
		} else {
			self.tree.size([1, 1]);
		}
		//set minimum height to 800
		if (self.tree.size()[0] < 800) {
			self.tree.size([800, 1]);
		}
		console.log(self.tree.size())
		
		// Compute the new tree layout.
		var nodes = self.tree.nodes(self.current_root); //.reverse();
		var links = self.tree.links(nodes);
		
		// deal with filtering complications
		var add_depth = 0;
		var add_nodes = [];
		while (new_root.parent) {
			add_depth++;
			if (new_root.parent._children) {
				new_root.parent.children = new_root.parent._children;
				new_root.parent._children = null;
			}
			new_root.parent.x = new_root.x;
			add_nodes.push(new_root.parent);
			links.push({
				source: new_root.parent,
				target: new_root,
			});
			new_root = new_root.parent;
		}
		
		// fix depth
		var max_depth = 0;
		for(var i = 0; i < nodes.length; i++) {
			nodes[i].depth += add_depth;
			max_depth = Math.max(max_depth, nodes[i].depth);
		}
		for(var i = 0; i < add_nodes.length; i++) {
			add_depth--;
			add_nodes[i].depth = add_depth;
			nodes.push(add_nodes[i]);
		}
		
		// Normalize for fixed-depth.
		nodes.forEach(function(d) { d.y = d.depth * 1260/max_depth; });
		
		// Update the nodes…
		self.node = self.svg.selectAll("g.node")
			.data(nodes, function(d) { return d.id || (d.id = ++self.i); });

		// Enter any new nodes at the parent's previous position.
		var nodeEnter = self.node.enter().append("g")
			.attr("class", "node")
			.attr("transform", function(d) { return "translate(" + source.y0 + "," + source.x0 + ")"; })
			.on("click", function(d) {
				self.click.call(self, d);
			})
			.on("dblclick", function(d) {
				self.double_click.call(self, d);
			});

		nodeEnter.append("circle")
			.attr("class", "tree_node noselect")
			.attr("r", 1e-6)
			.style("fill", self.color);

		nodeEnter.append("text")
			.attr("x", function(d) { return d.children || d._children ? -10 : 10; })
			.attr("dy", ".35em")
			.attr("text-anchor", function(d) { return d.children || d._children ? "end" : "start"; })
			.text(function(d) { return d.name; })
			.style("fill-opacity", 1e-6);

		// Transition nodes to their new position.
		var nodeUpdate = self.node.transition()
			.duration(self.duration)
			.attr("transform", function(d) { return "translate(" + d.y + "," + d.x + ")"; });

		nodeUpdate.select("circle")
			.attr("class", "tree_node noselect")
			.attr("r", 6)
			.style("fill", self.color);
		
		nodeUpdate.select("text")
			.attr("class", "noselect")
			.style("fill-opacity", 1);

		// Transition exiting nodes to the parent's new position.
		var nodeExit = self.node.exit().transition()
			.duration(self.duration)
			.attr("transform", function(d) { return "translate(" + source.y + "," + source.x + ")"; })
			.remove();

		nodeExit.select("circle")
			.attr("class", "tree_node")
			.attr("r", 1e-6);

		nodeExit.select("text")
			.style("fill-opacity", 1e-6);

		// Update the links…
		self.link = self.svg.selectAll("path.link")
			.data(links, function(d) { return d.target.id; });

		// Enter any new links at the parent's previous position.
		self.link.enter().insert("path", "g")
			.attr("class", "link")
			.attr("d", function(d) {
				var o = {x: source.x0, y: source.y0};
				return self.diagonal({source: o, target: o});
			})
			.attr("fill", "none")
			.attr("stroke", "#ccc")
			.attr("stroke-width", "1.5px");

		// Transition links to their new position.
		self.link.transition()
			.duration(self.duration)
			.attr("d", self.diagonal)
			.attr("fill", "none")
			.attr("stroke", "#ccc")
			.attr("stroke-width", "1.5px");

		// Transition exiting nodes to the parent's new position.
		self.link.exit().transition()
			.duration(self.duration)
			.attr("d", function(d) {
				var o = {x: source.x, y: source.y};
				return self.diagonal({source: o, target: o});
			})
			.attr("fill", "none")
			.attr("stroke", "#ccc")
			.attr("stroke-width", "1.5px")
			.remove();

		// Stash the old positions for transition.
		nodes.forEach(function(d) {
			d.x0 = d.x;
			d.y0 = d.y;
		});
	},
	filter: function(search_text) {
		var self = this;
		if (search_text) {
			var find = self.find.call(self, search_text, [self.root]);
			if(find) {
				self.update.call(self, self.root, find);
			} else {
				// Not found = no change
				//self.update.call(self, self.root, current.root);
			}
		} else {
			// Empty search text = reset filter
			self.update.call(self, self.root, self.root);
		}
	},
	find: function(search_text, nodes) {	//nodes is expected to be an array
		var self = this;
		//console.log(nodes);
		var new_nodes = [];
		if (nodes.length == 0) {
			return null;
		}
		for (var i = 0; i < nodes.length; i++) {
			var node = nodes[i];
			if (node.name.toUpperCase() == search_text.toUpperCase()) {
				return node;
			}
			if(node.children) {
				for (var j = 0; j < node.children.length; j++) {
					new_nodes.push(node.children[j]);
				}
			} else if (node._children) {
				for (var j = 0; j < node._children.length; j++) {
					new_nodes.push(node._children[j]);
				}
			}
		}
		return self.find.call(self, search_text, new_nodes);
	},
	color: function(d) {
		if (d.children) {
			return "white";
		} else if (d._children && d._children.length != 0) {
			return "rgb(154, 205, 50)";
		} else {
			//no children, hidden or otherwise
			return "gray";
		}
	},
	click: function(d) {
		var self = this;
		if (d.children) {
			d._children = d.children;
			d.children = null;
		} else {
			d.children = d._children;
			d._children = null;
		}
		self.update.call(self, d, self.current_root);
	},
	double_click: function(d) {
		var self = this;
		if (d.children) {
			self.collapse.call(self, d);
		} else {
			self.open.call(self, d);
		}
		tree.update.call(self, d, tree.current_root);
	},
	collapse: function(d) {
		var self = this;
		if (d.children) {
			d._children = d.children;
			d._children.forEach(function(d) {
				self.collapse.call(self, d);
			});
			d.children = null;
		}
	},
	open: function(d) {
		var self = this;
		if (d._children) {
			d.children = d._children;
			d.children.forEach(function(d) {
				self.open.call(self, d);
			});
			d._children = null;
		}
	},
});

function filter_tree() {
	var search_text = $('#tree_filter').val();
    tree.filter(search_text);
}