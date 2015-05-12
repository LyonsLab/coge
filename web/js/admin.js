var ITEM_TYPE_USER = 5; //TODO: This is duplicated elsewhere, move to a common location
var timestamps = new Array();
var jobs_timers = new Array();
var hist_timers = new Array();
var previous_search = ""; //indicates the previous search term, used to refresh after a delete
var jobs_updating = true;
var hist_updating = true;
var hist_entries = 0;
var last_hist_update;

$(function () {
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
    $("#update_checkbox").change(function(e) {
    	toggle_job_updater();
    });
    
    //Initialize Jobs tab
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
    //document.getElementById("job_search_bar").value = "running";
    
    
    //Initialize History tab
    
    
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
		hist.dataView.sort(comparer, args.sortAsc);
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
});

function search_stuff (search_term) {
	if(search_term.length > 2) {
		$("#loading_gears").show();
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
				//$('#loading_gears').hide();
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
						orgList = orgList + (obj.items[i].label) + " (ID: " + (obj.items[i].id) + ")" + "</span></td></tr>";
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
						genList = genList + "<button onclick='share_dialog(" + obj.items[i].id + ", 2 )'>Edit Access</button>";
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
						expList = expList + "<button onclick='share_dialog(" + obj.items[i].id + ", 3 )'>Edit Access</button>";
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
						noteList = noteList + "<button onclick='share_dialog(" + obj.items[i].id + ", 1 )'>Edit Access</button>";
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
				$("#loading_gears").show();
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
				
				$("#loading_gears").hide();
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
	} else {
		current_page = 0;
	}

	$('#master').fadeToggle('slow');
	$('#userInfo').fadeToggle('slow');
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
	$("#loading_gears2").show();
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
        					genList = genList + "<button onclick='share_dialog(" + current.id + ", 2 )'>Edit Access</button>";
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
        					expList = expList + "<button onclick='share_dialog(" + current.id + ", 3 )'>Edit Access</button>";
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
        					noteList = noteList + "<button onclick='share_dialog(" + current.id + ", 1 )'>Edit Access</button>";
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
        			genBlock = genBlock + "<span id='genCount" + i + "' class='coge-table-header' style='color:119911;'>Genomes: " + genCounter + " </span>";
        			genBlock = genBlock + "<div id=\"genArrow" + i + "\" onclick=\"toggle_arrow('#genArrow" + i + "');show_table('#genList" + i + "')\" style='display:inline;'>";
        			genBlock = genBlock + "<img src=\"picts/arrow-right-icon.png\" class=\"link\" style=\"width:10px;height:10px;\"/></div>";
        			genBlock = genBlock + "<table cellspacing=\"5\" class=\"hidden\" id='genList" + i + "' style=\"border-top:0px solid green; padding-left:20px; padding-bottom:10px;\">";
        			genBlock = genBlock + genList + "</table></div>";
        		}

        		if (noteCounter > 0) {
        			noteBlock = noteBlock + "<div style=\"padding-top:10px;padding-left:30px;\">";
        			noteBlock = noteBlock + "<span id='noteCount" + i + "' class='coge-table-header' style='color:119911;'>Notebooks: " + noteCounter + " </span>";
        			noteBlock = noteBlock + "<div id=\"noteArrow" + i + "\" onclick=\"toggle_arrow('#noteArrow" + i + "');show_table('#noteList" + i + "')\" style='display:inline;'>";
        			noteBlock = noteBlock + "<img src=\"picts/arrow-right-icon.png\" class=\"link\" style=\"width:10px;height:10px;\"/></div>";
        			noteBlock = noteBlock + "<table cellspacing=\"5\" class=\"hidden\" id='noteList" + i + "' style=\"border-top:0px solid green; padding-left:20px; padding-bottom:10px;\">"; 
        			noteBlock = noteBlock + noteList + "</table></div>";
        		}

        		if (expCounter > 0) {
        			expBlock = expBlock + "<div style=\"padding-top:10px;padding-left:30px;\">";
        			expBlock = expBlock + "<span id='expCount" + i + "' class='coge-table-header' style='color:119911;'>Experiments: " + expCounter + " </span>";
        			expBlock = expBlock + "<div id=\"expArrow" + i + "\" onclick=\"toggle_arrow('#expArrow" + i + "');show_table('#expList" + i + "')\" style='display:inline;'>";
        			expBlock = expBlock + "<img src=\"picts/arrow-right-icon.png\" class=\"link\" style=\"width:10px;height:10px;\"/></div>";
        			expBlock = expBlock + "<table cellspacing=\"5\" class=\"hidden\" id='expList" + i + "' style=\"border-top:0px solid green; padding-left:20px; padding-bottom:10px;\">";
        			expBlock = expBlock + expList + "</table></div>";
        		}

        		if (userCounter > 0) {
        			userBlock = userBlock + "<div style=\"padding-top:10px;padding-left:30px;\">";
        			userBlock = userBlock + "<span id='userCount" + i + "' class='coge-table-header' style='color:119911;'>Users: " + userCounter + " </span>";
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
        	$("#loading_gears2").hide();
        }
	});
}

function share_dialog(id, type) {
	var item_list = "content_" + id + "_" + type;  //selected.map(function(){return this.parentNode.id;}).get().join(',');
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


function modify_item (id, type, modification) {
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

function wait_to_search (search_func, search_term) {
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
	$.ajax({
		dataType: 'json',
	    data: {
	        jquery_ajax: 1,
	        fname: 'get_jobs',
	        time_range: 0,
	    },
	    success: function(data) {
	    	//console.log(data.jobs);
	        jobs.load(data.jobs);
	        entries = data.jobs.length;
	        $("#filter_busy").hide();
	        update_filter();
	    },
	    complete: function(data) {
	    	if (jobs_updating) {
	        	schedule_update("jobs", 1000);
	        }
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
		schedule_update("jobs", 1000);
	}
}

function schedule_update(page, delay) {
	console.log("Updating " + page);
	cancel_update(page);
	
	if (delay !== undefined) {
		if(page == "jobs") {
			jobs_timers['update'] = window.setTimeout(
					function() { get_jobs(); },
					delay
			);
		}
		if(page == "hist") {
			hist_timers['update'] = window.setTimeout(
					function() { update_history(); },
					delay
			);
		}
		return;
	}	
}

function cancel_update(page) {
	if(page == "job") {
		clearTimeout(job_timers['update']);
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
			//console.log(last_hist_update);
			//console.log(data[1].date_time);
			//console.log(last_hist_update > data[1].date_time);
			//console.log(last_hist_update < data[1].date_time);
			
			//dataView.beginUpdate();
			//dataView.setItems(data);
			//dataView.setFilterArgs({
			//	show: 0,
			//	searchType: 1,
			//	searchString: ''
			//});
			//dataView.setFilter(myFilter);
			//dataView.endUpdate();
			updateHistFilter();
		},
	    complete: function(data) {
	    	if (hist_updating && last_hist_update) {
	        	schedule_update("hist", 1000);
	        }
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
			//console.log(data);
	    	if (hist_updating) {
	        	schedule_update("hist", 1000);
	        }
	    }
	})
}

function toggle_hist_updater() {
	hist_updating = !hist_updating;
	if (hist_updating) {
		schedule_update("hist", 1000);
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
