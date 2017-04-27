var user_is_admin = false;
var ITEM_TYPE_USER = 5; //TODO: This is duplicated elsewhere, move to a common location
var timestamps = new Array();
//var jobs_timers = new Array();
var hist_timers = new Array();
var current_tab = 0;
var previous_search = ""; //indicates the previous search term, used to refresh after a delete
var jobs_grid;
var hist_grid;
var user_graph_init = false;
var group_graph_init = false;
//var hist_entries = 0;
var last_hist_update;
var IDLE_TIME = 30*1000; // stop polling after this lapse, then poll on next mousemove
var reports_grid;
var histogram;
var tree;
var system_graph;
var system_graph2;
var query_counter;
var tab7_init = false;

$(function () {
	// Initialize CoGe web services
    coge.services.init({
    	baseUrl: API_BASE_URL,
    	userName: USER_NAME
    });

	// Initialize tabs
	$( "#tabs" ).tabs({
		activate: function(event, ui) {
            current_tab = ui.newTab.index();
            schedule_update(5000);
			//Remember that current_tab is 0-indexed but #tabs is 1-indexed.
			if (current_tab == 1 && !jobs_grid) {
				init_jobs_grid();
			}
			if (current_tab == 2 && !hist_grid) {
				init_hist_grid();
			}
			if (current_tab == 4 && !reports_grid) {
				init_reports();
			}
			if (current_tab == 5) {
				init_taxon_tree("taxonomic_tree");
			}
			if (current_tab == 6) {
				if (!system_graph) {
					init_system_load();
				}
			}
			if (current_tab == 7) {
				if (!tab7_init) {
					var height = Math.max(400, $(window).height() - 250);
					$("#tabs-8").html('<iframe src="https://genomevolution.org/greentea/" height="'+height+'" width="100%"></iframe>');
					tab7_init = true;
				}
			}
			if (current_tab == 8) {
				if (!query_counter) {
					init_database_tab();
				}
			}
		}
    }).show();
	
	// Configure dialogs
    $(".dialog_box").dialog({autoOpen: false, width: 500});
    $("#job_search_bar").keyup(function (e) {
        Slick.GlobalEditorLock.cancelCurrentEdit();

        if (e.which == 27) // Clear on Esc
            this.value = "";
        update_filter();
    });
    $("#show_select,#job_search_type").change(function(e) {
        update_filter();
    });
    $("#histogram").dialog({autoOpen: false, width: 500, height: 500});
    
    //Setup idle timer
    timestamps['idle'] = new Date().getTime();
    $(document).mousemove(function() {
		var currentTime = new Date().getTime();
		var idleTime = currentTime - timestamps['idle'];
		timestamps['idle'] = currentTime;

		if (idleTime > IDLE_TIME) {
			// User was idle for a while, refresh page
			schedule_update(5000);
		}
	});
    
    //See if the current user is an admin
    $.ajax({
		data: {
			fname: 'user_is_admin',
		},
		success: function(data) {
			user_is_admin = (data == 1);
		}
	});
});

//Initialize the Jobs tab -- TODO move this code into an object similar to code for other tabs (mdb 11/9/16)
function init_jobs_grid() {
	jobs_grid = new JobGrid({
		elementId: "jobs",
		//filter: "none",
		running_only: 1,
		schedule_update: schedule_update,
		cancel_update: cancel_update,
	});
}

//Initialize the History tab
function init_hist_grid() {
	hist_grid = new HistGrid({
		elementId: "history",
		schedule_update: schedule_update,
		cancel_update: cancel_update,
	});
}

//initialize Reports tab
function init_reports() {
	reports_grid = new ReportGrid({
		elementId: "reports",
		filter: "none",
		selection: "total",
	});
}

//Initialize Database tab
function init_database_tab() {
	query_counter = new Query_Counter({
		elementId: "database",
		schedule_update: schedule_update,
		cancel_update: cancel_update,
	});
}

//Tab 1 Search
function search_stuff (search_term) {
	if(search_term.length > 2) {
		$("#loading").show();
		
		coge.services.search_global(search_term)
			.done(function(response) {
				console.log(response);
				var obj = response;
				
				if (!obj || !obj.results) {
					return;
				}
				
				var userCounter = 0, orgCounter = 0, genCounter = 0, expCounter = 0, noteCounter = 0, usrgroupCounter = 0;
				var userList = "", orgList = "", genList = "", expList = "", noteList = "", usrgroupList = "";

				for (var i = 0; i < obj.results.length; i++) {
					if (obj.results[i].type == "user") {
						userList = userList + "<tr><td><span>";
						userList = userList + obj.results[i].first + " " + obj.results[i].last + ", ";
						if (user_is_admin) {
							userList = userList + obj.results[i].email;
						}
						userList = userList + " (" + obj.results[i].username + ", id" + obj.results[i].id + ") ";
						userList = userList + "<button onclick=\"search_user(" + obj.results[i].id + ",'user')\">Show Data</button>";
						userList = userList + "</span></td></tr>";
						userCounter++;
					}

					if (obj.results[i].type == "organism") {
						orgList = orgList + "<tr><td><span title='" + obj.results[i].description + "'>";
						orgList = orgList + obj.results[i].name + " (id" + obj.results[i].id + ")";
						orgList = orgList + " <a href=\"OrganismView.pl?oid=" + obj.results[i].id + "\">Info </a>";
						orgList = orgList + "</span></td></tr>";
						orgCounter++;
					}
	
					if (obj.results[i].type == "genome") {
						genList = genList + "<tr><td><span onclick=\"modify_item(" + obj.results[i].id + ", 'Genome', 'restrict');\"";
						if (obj.results[i].restricted == 1) {
							genList = genList + " class=\"link ui-icon ui-icon-locked\"></span>";
						} else {
							genList = genList + " class=\"link ui-icon ui-icon-unlocked\"></span>";
						}
						genList = genList + "<span onclick=\"modify_item(" + obj.results[i].id + ", 'Genome', 'delete');\"";
						genList = genList + " class=\"link ui-icon ui-icon-trash\"></span>";
						if (obj.results[i].deleted == 1) {
							genList = genList + "<span style=\"color: red\">";
						} else {
							genList = genList + "<span>";
						}
						genList = genList + obj.results[i].name + " <a href=\"GenomeInfo.pl?gid=" + obj.results[i].id + "\">Info </a>";
						genList = genList + "<button class='access' onclick='share_dialog(" + obj.results[i].id + ", 2, " + obj.results[i].restricted + ")'>Edit Access</button>";
						genList = genList + "</span></td></tr>";
						genCounter++;
					}
	
					if (obj.results[i].type == "experiment") {
						expList = expList + "<tr><td><span onclick=\"modify_item(" + obj.results[i].id + ", 'Experiment', 'restrict');\"";
						if (obj.results[i].restricted == 1) {
							expList = expList + " class=\"link ui-icon ui-icon-locked\"></span>";
						} else {
							expList = expList + " class=\"link ui-icon ui-icon-unlocked\"></span>";
						}
						expList = expList + "<span onclick=\"modify_item(" + obj.results[i].id + ", 'Experiment', 'delete');\"";
						expList = expList + " class=\"link ui-icon ui-icon-trash\"></span>";
						if (obj.results[i].deleted == 1) {
							expList = expList + "<span style=\"color: red\">";
						} else {
							expList = expList + "<span>";
						}
						expList = expList + obj.results[i].name + " (id" + obj.results[i].id + ") <a href=\"ExperimentView.pl?eid=" + obj.results[i].id + "\">Info </a>";
						expList = expList + "<button class='access' onclick='share_dialog(" + obj.results[i].id + ", 3, " + obj.results[i].restricted + ")'>Edit Access</button>";
						expList = expList + "</span></td></tr>";
						expCounter++;
					}
	
					if (obj.results[i].type == "notebook") {
						noteList = noteList + "<tr><td><span onclick=\"modify_item(" + obj.results[i].id + ", 'List', 'restrict');\"";
						if (obj.results[i].restricted == 1) {
							noteList = noteList + " class=\"link ui-icon ui-icon-locked\"></span>";
						} else {
							noteList = noteList + " class=\"link ui-icon ui-icon-unlocked\"></span>";
						}
						noteList = noteList + "<span onclick=\"modify_item(" + obj.results[i].id + ", 'List', 'delete');\"";
						noteList = noteList + " class=\"link ui-icon ui-icon-trash\"></span>";
						if (obj.results[i].deleted == 1) {
							noteList = noteList + "<span style=\"color: red\">";
						} else {
							noteList = noteList + "<span>";
						}
						noteList = noteList + obj.results[i].name + " <a href=\"NotebookView.pl?lid=" + obj.results[i].id + "\">Info </a>";
						noteList = noteList + "<button class='access' onclick='share_dialog(" + obj.results[i].id + ", 1 , " + obj.results[i].restricted + ")'>Edit Access</button>";
						noteList = noteList + "</span></td></tr>";
						noteCounter++;
					}
							
					if (obj.results[i].type == "user_group") {
						usrgroupList = usrgroupList + "<tr><td>";
						if (obj.results[i].deleted == 1) {
							usrgroupList = usrgroupList + "<span style=\"color: red\">";
						} else {
							usrgroupList = usrgroupList + "<span>";
						}
						usrgroupList = usrgroupList + obj.results[i].name + " (id" + obj.results[i].id + ") ";
						usrgroupList = usrgroupList + "<button onclick=\"search_user(" + obj.results[i].id + ",'group')\">Show Data</button>";
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
				if (!user_is_admin) {
					$(".access").hide();
				} else {
					$(".access").show();
				}

			})
			.fail(function() {
				//TODO: Find some way to test this without breaking anything.
				$("#loading").hide();
				$("#masterTable").html("An error occured. Please reload the page and try again.");
			});
		
		previous_search = search_term;
	}
}

function show_table(id) {
	$(id).fadeToggle('fast');	
}

function toggle_arrow(id) {
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

// indicates the part of the search tab currently being displayed. 0 --> "master", 1 --> "user info"
var current_page = 0;

function search_user(userID, search_type) {
	if (current_page == 0) {	
		toggle_master();
		//current_page=1;
	}
	if(previous_user != userID) {
		$('#userResults').hide();
		//$('#userResults').html("Loading...");
		user_info(userID, search_type);
		$('#userResults').show();
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
			var obj = jQuery.parseJSON(data);
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
        					genList = genList + current.label + " <a href=\"GenomeInfo.pl?gid=" + current.id + "\">Info </a>";
        					genList = genList + "<button class='access' onclick='share_dialog(" + current.id + ", 2, " + current.restricted + ")'>Edit Access</button>";
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
        					expList = expList + current.label + " (id" + current.id + ") <a href=\"ExperimentView.pl?eid=" + current.id + "\">Info </a>";
        					expList = expList + "<button class='access' onclick='share_dialog(" + current.id + ", 3, " + current.restricted + ")'>Edit Access</button>";
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
        					noteList = noteList + current.label + " <a href=\"NotebookView.pl?lid=" + current.id + "\">Info </a>";
        					noteList = noteList + "<button class='access' onclick='share_dialog(" + current.id + ", , " + current.restricted + ")'>Edit Access</button>";
        					noteList = noteList + "</span></td></tr>";
        					noteCounter++;
        				}
						
        				if (current.type == "user") {
        					userList = userList + "<tr><td><span>";
        					userList = userList + current.first + " " + current.last + ", " + current.email + " (id" + current.id + ") ";
        					userList = userList + "<button onclick=\"search_user(" + current.id + ",'user')\">Search</button>";
        					userList = userList + "</span></td></tr>";
        					userCounter++;
        				}

        			}
        		} //end of single user loop

        		var genBlock = "", noteBlock = "", expBlock = "", userBlock = "", nameBlock = "";
				
        		nameBlock = nameBlock + "<div style=\"padding-top:10px;\">"
        		if (search_type == "user" && i == 0) {
        			nameBlock = nameBlock + "<img src='picts/user-icon.png' width='15' height='15'><span> ";
        			nameBlock = nameBlock + obj.items[i].first + " " + obj.items[i].last + ", " + obj.items[i].email + " (" + obj.items[i].username + ", id" + obj.items[i].user_id + "): </span>";
        		} else {
        			nameBlock = nameBlock + "<img src='picts/group-icon.png' width='15' height='15'><span> ";
        			nameBlock = nameBlock + obj.items[i].user + " (id" + obj.items[i].user_id + "): </span>";
        		}
        		
        		
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
        	if (!user_is_admin) {
				$(".access").hide();
			} else {
				$(".access").show();
			}
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
	pageObj.search_term = search_term;

	if (pageObj.time) {
		clearTimeout(pageObj.time);
	}

	// FIXME: could generalize by passing select id instead of separate search_* functions
	pageObj.time = setTimeout(
		function() {
			search_func(pageObj.search_term);
		},
		500
	);
}

/*
 * Jobs Tab (tab #2)
 */

function JobGrid(params) {
	this.elementId = params.elementId;
	//this.filter = params.filter;
	this.running_only = params.running_only;
	this.schedule_update = params.schedule_update;
	this.cancel_update = params.cancel_update;
	this.timers = new Array();
	this.updating = true;
	this.flag = false;
	this.data;
	this.table;
	
	var element = $('#' + this.elementId);
	this.width = element.outerWidth();
	this.height = Math.max(400, $(window).height() - 450); //element.outerWidth(); // mdb changed 2/8/16

	this.initialize();
}

$.extend(JobGrid.prototype, {
	initialize: function() {
		var self = this;
		$('#' + self.elementId).hide();
		$('#' + self.elementId + '_loading').show();
		//$('#report_filter').prop('disabled', true);
		//$('#histogram_button').prop('disabled', true);
		$('#' + this.elementId).html('<table id="' + this.elementId + '_table" cellpadding="0" cellspacing="0" border="0" class="dt-cell hover compact row-border">'
				+ '<thead><tr><th>ID</th><th>Started</th><th>Completed</th><th>Elapsed</th><th>User</th><th>Tool</th><th>Link to Analysis</th><th>Status</th></tr></thead></table>');
		$('#' + self.elementId + '_update_checkbox').change(function(e) {
    		self.toggle_updating.call(self);
    	});
		$('#' + self.elementId + '_running_checkbox').change(function(e) {
			self.toggle_running.call(self);
		});
		$('#' + self.elementId + '_cancel').click( function () {
			self.cancel_job.call(self);
		});
		$('#' + self.elementId + '_restart').click( function () {
			self.restart_job.call(self);
		});
		
		//Setup table formatting
		self.table = $('#' + self.elementId + '_table').DataTable({
    		columnDefs : [
				{ 
					orderSequence : [ "desc", "asc" ], 
					targets : [0, 1, 2] 
				},
				{
					"render": function ( data, type, full, meta ) {
						if (meta.col == 0)
							return '<a href="'+'jex/status/'+data+'" target="_blank">'+data+'</a>';
						if (meta.col == 6)
							return '<a href="'+data+'" target="_blank">'+data+'</a>';
					},
					targets : [ 0, 6 ]
				}
    		],
    		iDisplayLength: Math.floor(self.height/24), //Each row is 24 pixels tall
    		order: [[1, "desc"]],
    		scrollY: self.height,
    		lengthChange: false,
	    });
		
		self.get_data.call(self);
	},
	get_data: function() {
		var self = this;
		if(!self.flag) {
			self.flag = true;
			self.cancel_update();
			//$("#" + self.elementId + "_update_checkbox").prop('disabled', true);
		    //$("#" + self.elementId + "_running_checkbox").prop('disabled', true);
			$.ajax({
				dataType: 'json',
			    data: {
			        fname: 'get_jobs',
			        time_range: 0,
			        status: self.running_only ? 'running' : null,
			    },
			    success: function(data) {
			    	//console.log(data)
			    	self.data = data.data;
			    	self.table
				    	.clear()
				    	.rows.add(self.data)
				    	.draw();
					
			    	$('#' + self.elementId + '_loading').hide();
					$('#' + self.elementId).show();
			    },
			    complete: function(data) {
			    	self.flag = false;

			    	$('#' + self.elementId + '_table tbody')
			    	    .off( 'click' )
					    .on( 'click', 'tr', function () {
                            if ( $(this).hasClass('selected') ) {
                                $(this).removeClass('selected');
                            }
                            else {
                                self.table.$('tr.selected').removeClass('selected');
                                $(this).addClass('selected');
                            }
                        } );
			    	
			    	self.schedule_update(10*1000);
			    }
			});
		}
	},
	update: function(delay) {
		var self = this;
		if (self.updating) {
			console.log("Updating jobs");
			self.cancel_update();
			self.timers['update'] = window.setTimeout(
				function() { self.get_data.call(self); },
				delay
			);
		}
	},
	toggle_updating: function() {
		var self = this;
		self.updating = !self.updating;
		if (self.updating) {
			self.schedule_update(5000);
		}
	},
	toggle_running: function() {
		var self = this;
		self.running_only = (self.running_only ? 0 : 1);
		self.get_data();
	},
	cancel_job: function() {
		var self = this;
		var parent_id = self.table.row('.selected').data()[0];
	    self.submit_task("cancel_job", parent_id);
		
		self.get_data.call(self);
	},
	restart_job: function() {
		var self = this;
		var parent_id = self.table.row('.selected').data()[0];
	    self.submit_task("restart_job", parent_id);

	    self.get_data.call(self);
	},
	submit_task: function(task, parent_id) {
		var self = this;
		var argument_list =  {
			fname: task,
			job: parent_id,
		};

		$.ajax({
			type: "GET",
			dataType: "json",
			data: argument_list,
			success: function(data) {
				//DO NOTHING (atm)
			}
		});
	}
});

function schedule_update(delay) {
	var idleTime = new Date().getTime() - timestamps['idle'];
	if (idleTime < IDLE_TIME && delay !== undefined) {
		if (current_tab == 1 && jobs_grid)
			jobs_grid.update(delay);
		if (current_tab == 2 && hist_grid)
			hist_grid.update(delay);
		if (current_tab == 8 && query_counter)
			query_counter.update(delay);
		return;
	}	
}

function cancel_update() {
	if (jobs_grid)
		clearTimeout(jobs_grid.timers['update']);
	if (hist_grid)
		clearTimeout(hist_grid.timers['update']);
	if (query_counter)
		clearTimeout(query_counter.timers['update']);
}

function HistGrid(params) {
	this.elementId = params.elementId;
	this.schedule_update = params.schedule_update;
	this.cancel_update = params.cancel_update;
	this.timers = new Array();
	this.updating = true;
	this.flag = false;
	this.data;
	this.table;
	this.last_update;
	this.oldest_timestamp;
	
	this.width = $('#' + this.elementId).outerWidth();
	this.height = Math.max(400, $(window).height() - 370); //this.height = $('#' + this.elementId).outerHeight() - 100; // mdb changed 11/9/16
	
	this.initialize();
}

$.extend(HistGrid.prototype, {
	initialize: function() {
		var self = this;
		$('#' + self.elementId).hide();
		$('#' + self.elementId + '_loading').show();
		$('#' + this.elementId).html('<table id="' + this.elementId + '_table"'
		      + ' cellpadding="0" cellspacing="0" border="0" class="dt-cell hover compact row-border">'
			  + '<thead><tr><th>Date/Time</th><th>User</th><th>Page</th><th>Description</th><th>Link</th><th>Comments</th></tr></thead></table>');

		$('#' + self.elementId + '_update_checkbox').change(function(e) {
    		self.toggle_updating.call(self);
    	});
		
		//Setup table format
		self.table = $('#' + self.elementId + '_table').DataTable({
    		columnDefs: [
    		             { 
    		            	 orderSequence : [ "desc", "asc" ], 
    		            	 targets : 0 
    		             },
    		             {
    		            	 "render": function ( data, type, full, meta ) {
    		            	      return '<a href="'+data+'" target="_blank">'+data+'</a>';
    		            	 },
    		            	 targets : 4
    		             }
    		            ],
    		columns:	[
    		        	 { data: 0 },
    		        	 { data: 1 },
    		        	 { data: 2 },
    		        	 { data: 3 },
    		        	 { data: 4 },
    		        	 { data: 5 }
    		        	],
    		iDisplayLength: Math.floor(self.height/24), //Each row is 24 pixels tall
    		order: [[0, "desc"]],
    		scrollY: self.height,
    		lengthChange: false,
	    });
		
		self.get_data.call(self);
	},
	get_data: function() {
		var self = this;
		if (!self.flag) {
			self.flag = true;
			self.cancel_update();
			
			$.ajax({
				dataType: 'json',
			    data: {
			        fname: 'get_history',
			        time_range: 0,
			    },
			    success: function(data) {
			    	//console.log(data)
			    	self.data = data.data;
			    	self.table
				    	.clear()
				    	.rows.add(self.data)
				    	.draw();

			    	//Record most recent timestamp
			    	if (!self.last_update)
			    		self.last_update = self.data[0][0];

			    	//Record oldest timestamp
			    	self.oldest_timestamp = self.data[self.data.length - 1][0];
			    	//console.log(self.oldest_timestamp);
					
			    	$('#' + self.elementId + '_loading').hide();
					$('#' + self.elementId).show();

					//Populate remaining pages
					if (self.data.length > 0)
						self.get_more_data();
			    },
			    complete: function(data) {
			    	self.flag = false;

			    	$('#' + self.elementId + '_table tbody')
			    	    .off('click')
					    .on('click', 'tr', function () {
                            if ( $(this).hasClass('selected') ) {
                                $(this).removeClass('selected');
                            }
                            else {
                                self.table.$('tr.selected').removeClass('selected');
                                $(this).addClass('selected');
                            }
                        } );
			    }
			});
		}
	},
	get_more_data: function() { // TODO combine this with get_data() above -- mdb 11/9/16
		var self = this;
		$.ajax({
			dataType: 'json',
		    data: {
		        fname: 'get_history',
		        time_range: 0,
		        timestamp: self.oldest_timestamp,
		    },
		    success: function(data) {
		    	//console.log(data)
		    	self.data = data.data;
		    	self.table
			    	//.clear()
			    	.rows.add(self.data)
			    	.draw(false);
			   	
			   	//Record oldest timestamp && recurse
		    	if(self.data.length > 0) {
		    		self.oldest_timestamp = self.data[self.data.length - 1][0];
					self.get_more_data();
				} else {
					self.schedule_update(5000);
				}
		    }
		});
	},
	update: function(delay) {
		var self = this;
		if (self.updating) {
			console.log("Updating history");
			self.cancel_update();
			self.timers['update'] = window.setTimeout(
				function() {
					$.ajax({
						dataType: 'json',
						data: {
							fname: 'update_history',
							timestamp: self.last_update,
							time_range: 0,
						},
						success: function(data) {
							//console.log(data);
							if(data.new_rows[0]) {
								self.table.rows.add(data.new_rows).draw(false);
								self.last_update = data.new_rows[0][0];
							}
						},
						complete: function(data) {
							self.schedule_update(5000);
					    }
					})
				},
				delay
			);
		}
	},
	toggle_updating: function() {
		var self = this;
		self.updating = !self.updating;
		if (self.updating)
			self.schedule_update(5000);
	}
});

//Tab 4, the User Graph
var colors = [
        	{ name: 'list',       link: 'NotebookView.pl?nid=',   color: 'Tomato',      show: 1 },
        	{ name: 'genome',     link: 'GenomeInfo.pl?gid=',     color: 'YellowGreen', show: 1 },
        	{ name: 'experiment', link: 'ExperimentView.pl?eid=', color: 'Orchid',      show: 1 },
        	{ name: 'feature',    link: '',                       color: 'orange',      show: 1 },
        	{ name: 'user',       link: '',                       color: 'DeepSkyBlue', show: 1 },
        	{ name: 'group',      link: 'GroupView.pl?ugid=',     color: 'Turquoise',   show: 1 },
];
var user_force,
	group_force;

function init_graph(selectId) {
	if(selectId == 1) {
		$('#group_chart').hide();
		$('#group_legend').hide();
		$('#user_chart').show();
		$('#user_legend').show();
		
		if(!user_graph_init) {	
			user_graph_init = true;
			$("#loading4").show();
			d3.json("?fname=get_user_nodes", function(json) {
				var user_filters = [
				                    { name: 'deleted', 		color: 'Red', 		show: 0 },
				                    { name: 'restricted', 	color: 'Grey',	 	show: 0 },
				                    { name: '||||||||',		color: 'White',		show: 1 },
				];
				
				var user_force = new Force(json, user_filters, "user", colors);
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
				var group_filters = [
				                     { name: 'deleted', 		color: 'Red', 		show: 0 },
				                     { name: 'restricted', 		color: 'Grey',	 	show: 0 },
				                     { name: '||||||||',		color: 'White',		show: 1 },
				];
				
				var group_force = new Force(json, group_filters, "group", colors);
				$("#loading4").hide();
				group_force.update();
			});
		}
	}
}

var Force = function(root, filters, element, colors) {
	var self = this;
	this.width = Math.max(800, $(window).width()-200);
	this.height = Math.max(800, $(window).height()-300);
	this.root = root;
	this.colors = colors;
	this.filters = filters;
	this.element = element;
	
	this.svg = d3.select('#' + self.element + '_chart').append("svg")
		.attr("width", self.width)
		.attr("height", self.height);
	
	this.link = self.svg.selectAll(".link");
	this.node = self.svg.selectAll(".node");
	
	var nodes = self.flatten(root);
	nodes.forEach(function(d) {
		d._children = d.children;
		d.children = null;
	});
	
	self.filters.forEach(function(element, index) {
		var item =
			$('<div class="link legend selected">'+self.filters[index].name+'</div>')
				.css('color', self.filters[index].color)
				.css('background-color', '')
				.click(function() {
					$(this).toggleClass('selected');
					if ($(this).hasClass('selected')) {
						$(this).css('color', self.filters[index].color);
						$(this).css('background-color', '');
					}
					else {
						$(this).css('color', 'white');
						$(this).css('background-color', self.filters[index].color);
					}
					self.filters[index].show = !self.filters[index].show;
					
					var nodes = self.flatten(root);
					nodes.forEach(function(d) {
						return self.color.call(self, d);
					});
					self.update();
				});

		$('#' + self.element + '_legend')
			.append(item);
	});
	
	self.colors.forEach(function(element, index) {
		var item =
			$('<div class="link legend selected">'+self.colors[index].name+'</div>')
				.css('color', 'white')
				.css('background-color', self.colors[index].color);

		$('#' + self.element + '_legend')
			.append(item);
	});
	
	this.force = d3.layout.force()
		.size([self.width, self.height])
		.on("tick", function() {
			self.tick.call(self);
		})
		.gravity(1);
}

$.extend(Force.prototype, {
	update: function() {
		var self = this;
		var nodes = self.flatten(self.root),
			links = d3.layout.tree().links(nodes);
		
		// Restart the force layout.
		this.force
	 		.nodes(nodes)
	 		.links(links)
	 		.start();

		// Update the links�
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

		// Update the nodes�
		this.node = this.node.data(nodes, function(d) { return d.id; })
			.style("fill", function(d) {
				return self.color(d, self.filters);
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
				return self.color.call(self, d);
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
	},
	tick: function() {
		/*this.link.attr("x1", function(d) { return d.source.x; })
			.attr("y1", function(d) { return d.source.y; })
			.attr("x2", function(d) { return d.target.x; })
			.attr("y2", function(d) { return d.target.y; });
		
		this.node.attr("cx", function(d) { return d.x; })
			.attr("cy", function(d) { return d.y; });*/
		
		var self = this;
		self.moveItems.call(self);
	},
	toggle_children: function(d) {
		if (d.children) {
			d._children = d.children;
			d.children = null;
		} else {
			d.children = d._children;
			d._children = null;
		}
	},
	click: function(d) {
		if (!d3.event.defaultPrevented) {
			this.force.charge(
				function(d, i) {
					if(d.size) {
						return Math.sqrt(d.size)*(-2.25);
					}
				}
			);
			this.force.start();
			
			this.toggle_children(d);
			this.update();
		}
	},	
	flatten: function(root) {	//Returns a list of all nodes under the root.
		var nodes = [], i = 0;

		function recurse(node) {
			if (node.children) node.children.forEach(recurse);
			if (!node.id) node.id = ++i;
			nodes.push(node);
		}

		recurse(root);
		return nodes;
	},
	color: function(d) {
		var self = this;
		if(!self.filters[0].show && !self.filters[1].show) {
			//normal filter
			if((d.children || d._children) || (d.type == 2 || d.type == 3 || d.type == 4)) {
				if (d.type) {
					return self.colors[d.type-1].color;
				}
				return 'white';
			} else {
				return 'black';
			}
		} else {
			if(d.deleted == 1 && d.restricted == 1 && self.filters[0].show && self.filters[1].show) {
				return 'yellow';
			} else if (d.deleted == 1 && self.filters[0].show) {
				return self.filters[0].color;
			} else if (d.restricted == 1 && self.filters[1].show) {
				return self.filters[1].color;
			} else {
				return 'yellowGreen';
			}
		}
	},
	moveItems: (function() {		//Optimization of the Tick function.
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
	 
	})(),
});

//Tab 5, Reports
function ReportGrid(params) {
	this.elementId = params.elementId;
	this.filter = params.filter;
	this.selection = params.selection;
	this.data;
	
	this.initialize();
}

$.extend(ReportGrid.prototype, {
	initialize: function() {
		var self = this;
		var element = this.elementId;
		var fname;
		let jobs = ['API','CoGeBlast','GEvo','LoadAnnotation','LoadExperiment','LoadGenome','SynFind','SynMap','SynMap2','SynMap3D'];
		$('#' + element).hide();
		$('#reports_loading').css('visibility', 'visible');
		$('#report_type').prop('disabled', true);
		$('#report_filter').prop('disabled', true);
		$('#histogram_button').prop('disabled', true);
		switch (this.selection) {
			case "group":
				fname = 'get_group_table';
				$('#' + element).html('<table id="' + element + '_table" cellpadding="0" cellspacing="0" border="0" class="dt-cell hover compact row-border stripe">'
					+ '<thead><tr><th>Group Name</th><th>Genomes</th><th>Experiments</th><th>Notebooks</th><th>Total</th><th>Users</th></tr></thead>'
					+ '<tfoot><tr style="font-weight:bold"><td style="text-align:right">Totals:</td><td></td><td></td><td></td><td></td><td></td></tr></tfoot></table>');
				break;
			case "jobs":
				fname = 'get_jobs_table';
				$('#' + element).html('<table id="' + element + '_table" cellpadding="0" cellspacing="0" border="0" class="dt-cell hover compact row-border stripe">'
					+ '<thead><tr><th>Job</th><th>Times Run</th></tr></thead>'
					+ '<tfoot><tr style="font-weight:bold"><td style="text-align:right">Total:</td><td></td></tr></tfoot></table>');
				break;
			case "user":
				fname = 'get_user_table';
				$('#' + element).html('<table id="' + element + '_table" cellpadding="0" cellspacing="0" border="0" class="dt-cell hover compact row-border stripe">'
					+ '<thead><tr><th>User Name</th><th>Genomes</th><th>Experiments</th><th>Notebooks</th><th>Total</th><th>Groups</th></tr></thead>'
					+ '<tfoot><tr style="font-weight:bold"><td style="text-align:right">Totals:</td><td></td><td></td><td></td><td></td><td></td></tr></tfoot></table>');
				break;
			case "user jobs":
				fname = 'get_user_jobs_table';
				let html = '<table id="' + element + '_table" cellpadding="0" cellspacing="0" border="0" class="dt-cell hover compact row-border stripe"><thead><tr><th>User Name</th>';
				jobs.forEach(function(job){
					html += '<th>' + job + '</th>';
				});
				html += '<th>Total</th></tr></thead><tfoot><tr style="font-weight:bold"><td style="text-align:right">Total:</td>';
				jobs.forEach(function(job){
					html += '<td></td>';
				});
				html += '<td></td></tr></tfoot></table>';
				$('#' + element).html(html);
				break;
			case "total":
				fname = 'get_total_table';
					$('#' + element).html('<table id="' + element + '_table" cellpadding="0" cellspacing="0" border="0" class="dt-cell hover compact row-border stripe">'
						+ '<thead><tr><th></th><th>Genomes</th><th>Experiments</th><th>Notebooks</th><th>Total</th><th>Users</th><th>Groups</th></tr></thead></table>');
				break;
		}
		$.ajax({
			data: {
				fname: fname,
				filter: this.filter,
			},
			success: function(data) {
				$('#reports_loading').css('visibility', 'hidden');
				var json = JSON.parse(data);
				if (self.selection != 'jobs')
					for (var i=0; i<json.data.length; i++) {
						let num_cols = self.selection == 'user jobs' ? jobs.length : 3;
						var total_items = 0;
						for(var j=1; j<=num_cols; j++)
							if (json.data[i][j]) {
								var num = parseInt(json.data[i][j].replace(',', ''));
								total_items += num;
							}
						if (self.selection == 'total') {
							json.data[i].push(json.data[i][5]);
							json.data[i][5] = json.data[i][4];
						} else
							json.data[i].push(json.data[i][num_cols + 1]);
						json.data[i][num_cols + 1] = total_items > 0 ? total_items : null;
					}
				self.data = json;
				if (self.selection != 'total')
					self.data.footerCallback = function ( row, data, start, end, display ) {
						var api = this.api();
						var add_total = function(col) {
							var total = api.column(col).data().reduce(function(a, b) { return b ? typeof b === 'number' ? a + b : a + parseInt(b.replace(',', '')) : a; }, 0 );
							if (total > 0)
								$(api.column(col).footer()).html(total.toLocaleString());
						}
						let num_cols = self.selection == 'jobs' ? 1 : self.selection == 'user jobs' ? jobs.length + 1 : 5;
						for (let i=1; i<=num_cols; i++)
							add_total(i);
					};
				self.data.scrollY = window.innerHeight - 350;
				self.data.scrollCollapse = true;
				$('#' + element).show();
				$('#' + element + '_table').dataTable(self.data);
				if (self.selection == 'user' || self.selection == 'group') {
					$('#histogram_button').prop('disabled', false);
					$('#report_filter').prop('disabled', false);
				}
				$('#report_type').prop('disabled', false);
			}
		});
		
	}
});

function user_job_plot(user_id, user_name) {
	var $dialog = $('<div id="plot" style="overflow:hidden"><img src="picts/ajax-loader.gif" style="width:16px;"/></div>')
    .dialog({
        height: 500,
        width: 1000,
        title: 'Jobs run by ' + user_name,
		close: function() {
			$dialog.dialog('destroy').remove();
		}
	});
	$.ajax({
		data: {
			fname: 'get_user_jobs',
			user_id: user_id
		},
		success: function(data) {
			$('#plot').empty();
			Plotly.newPlot('plot', JSON.parse(data));
		}
	});
}

function change_report(report) {
	if (report == 'Data Summary')
		reports_grid.selection = 'total';
	else if (report == 'Jobs Summary')
		reports_grid.selection = 'jobs';
	else if (report == 'Groups')
		reports_grid.selection = 'group';
	else if (report == 'Users Data')
		reports_grid.selection = 'user';
	else if (report == 'Users Jobs')
		reports_grid.selection = 'user jobs';
	reports_grid.initialize();
}

function filter_report(index) {
	var s = $('#report_filter')[0];
	reports_grid.filter = s.options[s.selectedIndex].value;
	reports_grid.initialize();
}

var Histogram = function(element, json, parent) {
	var self = this;
	this.element = element;
	this.json = json;
	this.parent = parent;
	this.values = this.json.values;
	this.min_value = this.json.min;
	this.max_value = this.json.max;
	
	$('#' + this.element).dialog('open');
	
	this.margin = {top: 10, right: 30, bottom: 80, left: 50};
	this.width = $('#' + this.element).outerWidth() - this.margin.left - this.margin.right;
	this.height = $('#' + this.element).outerHeight() - this.margin.top - this.margin.bottom;
	
	this.initialize();
}

$.extend(Histogram.prototype, {
	initialize: function() {
		var self = this;
		$("#" + this.element).html(
			'<div><div id="' + self.element + '_back_button" class="coge-button" style="margin-right:' + (this.width/2 - 150) + 'px;">Zoom Out</div>' +
			'<span style="margin-right:40px;">Data: ' + $('#report_type').val() + ', Filter: ' + $('#report_filter').val() + '</span></div>'
		);
		$('#' + this.element + '_back_button').on("click", function() {
			self.zoom_out.call(self);
		});

		if (self.parent) {
			$('#' + self.element + '_back_button').prop('disabled', false);
		} else {
			$('#' + self.element + '_back_button').prop('disabled', true);
		}
	
		// A formatter for counts.
		this.formatCount = d3.format(",.0f");
	
		this.x = d3.scale.linear()
		    .domain([this.min_value, this.max_value])
		    .range([0, this.width]);
	
		// Generate a histogram using twenty uniformly-spaced bins.
		this.data = d3.layout.histogram()
		    .bins(this.x.ticks(this.max_value - this.min_value))
		    (this.values);
	
		this.y = d3.scale.log()
		    .domain([1, Math.max(
		    		d3.max(this.data, function(d) { return d.y; }),
		    		2
		    )])
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
		    	if (d.y == 1) {
		    		return "translate(" + self.x(d.x) + ", " + self.y(1.1) + ")";
		    	} else if (d.y == 0) {
		    		return "translate(" + self.x(d.x) + ", " + self.y(1.1) + ")";
		    	}
		    	return "translate(" + self.x(d.x) + "," + self.y(d.y) + ")"; 
		    });
		
		this.bar.append("rect")
		    .attr("x", 1)
		    .attr("width", (self.width / self.data.length)*0.9)
		    .attr("height", function(d) {
		    	if(d.y == 1) {
		    		return self.height - self.y(1.1);
		    	} else if (d.y == 0) {
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

		var x_label = this.svg.append("text")
			.attr("font-size", "1.25em")
			.attr("transform", "translate(" + (this.width/2 - 150) + "," + (this.height + this.margin.bottom/2.5) + ")")
			.text("Total Data Sets (Notebooks + Genomes + Experiments)");
		
		this.svg.append("g")
	    	.attr("class", "y axis")
	    	.call(this.yAxis);
		
		var y_label = this.svg.append("text")
			.attr("font-size", "1.25em")
			.attr("transform", "translate(" + this.margin.left/(-1.25) + "," + this.height/2 + ") rotate(270)")
			.text("Frequency");
		
		this.svg.append("g")
	    	.attr("class", "x brush")
	    	.call(self.brush)
	    	.selectAll("rect")
	    	.attr("y", -6)
	    	.attr("height", this.height + 7);
		
		$('#' + this.element)
			.dialog({
				resizeStop: function( event, ui ) {
					self.width = $(this).outerWidth() - self.margin.left - self.margin.right;
					self.height = $(this).outerHeight() - self.margin.top - self.margin.bottom;
					self.initialize();
				},
			})
			.dialog('open');
	},
	brushed: function() {
		var self = this;
		var extent = this.brush.extent();
		extent[0] = Math.floor(extent[0]);
		extent[1] = Math.ceil(extent[1]);
		
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
		var parent = self;
		self = new Histogram(self.element, newJson, parent);
	},
	zoom_out: function() {
		var self = this;
		if (self.parent) {
			self = self.parent;
			self.child = null;
		}
		self.initialize();
	}
});

function init_histogram(element) {
	var json = reports_grid.data;
	var values = [];
	var max_value = 0;
	var min_value = Number.MAX_SAFE_INTEGER;
	for(var i = 0; i < json.data.length; i++) {
		var total_items = json.data[i][json.data[i].length - 1];
		values.push(total_items);
		if(total_items > max_value) {
			max_value = total_items;
		}
		if(total_items < min_value) {
			min_value = total_items;
		}
	}
	histogram = new Histogram(element, {min: min_value, max: max_value, values: values}, null);
}

function init_taxon_tree(element) {
	$("#loading6").show();
	if (!tree) {
		$.ajax({
			data: {
				fname: "gen_tree_json",
			},
			success: function(data) {
				tree = new Taxon_tree(JSON.parse(data), element);
				$("#loading6").hide();
			}
		});
	}
}

var Taxon_tree = function(json, element) {
	var container_width = $("#tabs").width();
	this.margin = {top: 20, right: 120, bottom: 20, left: 120};
	this.width = container_width - this.margin.right - this.margin.left;
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
	
	this.shifted;
	var object = this;
	$(document).on('keyup keydown', function(e){object.shifted = e.shiftKey;} );
	
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
		//show inital tree
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
		//set maximum height to 3100
		if (self.tree.size()[0] > 3100) {
			self.tree.size([3100, 1]);
		}
		
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
		for (var i = 0; i < nodes.length; i++) {
			nodes[i].depth += add_depth;
			max_depth = Math.max(max_depth, nodes[i].depth);
		}
		for (var i = 0; i < add_nodes.length; i++) {
			add_depth--;
			add_nodes[i].depth = add_depth;
			nodes.push(add_nodes[i]);
		}
		
		// Normalize for fixed-depth.
		nodes.forEach(function(d) { 
			if (max_depth == 0) {
				d.y = 0;
			} else {
				d.y = d.depth * (self.width-100)/max_depth;
			}
		});
		
		// Update the nodes�
		self.node = self.svg.selectAll("g.node")
			.data(nodes, function(d) { return d.id || (d.id = ++self.i); });

		// Enter any new nodes at the parent's previous position.
		var nodeEnter = self.node.enter().append("g")
			.attr("class", "node")
			.attr("transform", function(d) { return "translate(" + source.y0 + "," + source.x0 + ")"; })
			.on("click", function(d) {
				if (self.shifted) {
					self.double_click.call(self, d);
				} else {
					self.click.call(self, d);
				}
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

		// Update the links�
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
		$("#loading6").show();
		if (search_text) {
			var find = self.find.call(self, search_text, [self.root]);
			if(find) {
				self.update.call(self, self.root, find);
				$("#not").hide();
				//does not show $("#loading6").show();
			} else {
				// Not found = reset filter
				//does not show $("#loading6").show();
				self.update.call(self, self.root, self.root);
				$("#not").show();
			}
			setTimeout(function(){ $("#loading6").hide(); }, 500);
		} else {
			// Empty search text = reset filter
			self.update.call(self, self.root, self.root);
			//hides the loading6 image
			$("#loading6").hide();
			$("#not").hide();
		}
	},
	find: function(search_text, nodes) {	//nodes is expected to be an array
		var self = this;
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
	//loads the found tree
	if (tree) {
		var search_text = $('#tree_filter').val();
		tree.filter(search_text);
	}
}

function init_line_graph(index) {
	if (index == 0) {
		if (system_graph) {
			system_graph.initialize();
		} else {
			init_system_load();
		}
		$("#system_graph2").hide();
		$("#system_graph").show();
	} else {
		if (system_graph2) {
			system_graph2.initialize();
		} else {
			init_system_load2();	
		}
		$("#system_graph").hide();
		$("#system_graph2").show();
	}
}

function init_system_load() {
	$.ajax({
		url: "https://genomevolution.org/coge/data/system_load.txt",
		success: function (data){
			var strings =  data.split("\n");
			var json = [];
			for (var i = 0; i < strings.length - 1; i++) {
				strings[i] = strings[i].split("\t");
				json.push({
					"time": strings[i][0],
					"load": strings[i][3],
					"memory": parseInt(strings[i][4], 10)/1000,
				});
			}
			system_graph = new System_graph(json, "system_graph", null);
		}
	});
}

function init_system_load2() {
	$.ajax({
		url: "https://geco.iplantcollaborative.org/coge/data/system_load.txt",
		success: function (data){
			var strings =  data.split("\n");
			var json = [];
			for (var i = 0; i < strings.length - 1; i++) {
				strings[i] = strings[i].split("\t");
				json.push({
					"time": strings[i][0],
					"load": strings[i][3],
					"memory": parseInt(strings[i][4], 10)/1000,
				});
			}
			system_graph2 = new System_graph(json, "system_graph2", null);
		}
	});
}

var System_graph = function(data, element, parent) {
	var self = this;
	this.parent = parent;
	this.child = null;
	this.margin = {top: 30, right: 80, bottom: 30, left: 50},
	this.width = $(window).width()*.75 - this.margin.left - this.margin.right,
	this.height = $(window).height()*.5 - this.margin.top - this.margin.bottom;
	this.data = data;
	this.element = element;
	this.initialize();
}

$.extend(System_graph.prototype, {
	initialize: function() {
		var self = this;
		
		// Clear the element, add the zoom out button and svg container
		$("#" + this.element).html(
		        '<div class="inline">' +
				'<div id="' + self.element + '_back_button" class="coge-button" style="width:7em;">Zoom Out</div>' +
				'<div id="' + self.element + '_zoom_month_button" class="coge-button" style="width:7em;">Zoom Month</div>' +
				'<div id="' + self.element + '_zoom_week_button" class="coge-button" style="width:7em;">Zoom Week</div>' +
				'<div id="' + self.element + '_zoom_button" class="coge-button" style="width:7em;">Zoom Day</div>' +
				'</div>' +
				'<div><br><br></div>'+
				'<div id="' + self.element + '_container" style="height:' + (self.height+100) + 'px;">' +
				'<div id="' + self.element + '_graph" style="float:left;"></div></div>'
		);

		// Register button events
		if (self.parent) {
			$('#' + self.element + '_back_button')
				.on("click", function() {
					self.zoom_out.call(self);
				});
		}
        $('#' + self.element + '_zoom_button')
            .on("click", function() {
                self.zoom_hours.call(self, 24); // past day
        });
        $('#' + self.element + '_zoom_week_button')
            .on("click", function() {
                self.zoom_hours.call(self, 7*24); // past week
        });
        $('#' + self.element + '_zoom_month_button')
            .on("click", function() {
                self.zoom_hours.call(self, 30*24); // past month
        });

		// Parse the date / time
		self.timeFormat = d3.time.format("%H:%M:%S %m/%d/%Y");
		
		// Format the data
		self.data.forEach(function(d) {
			if (Object.prototype.toString.call(d.time) !== "[object Date]")
				d.time = self.timeFormat.parse(d.time);
			if (d.load)
				d.load = +d.load;
	        if (d.memory)
	        	d.memory = +d.memory;
	    });

		// Set the ranges
		self.x = d3.time.scale().range([0, self.width]);
		self.yLoad = d3.scale.linear().range([self.height, 0]);
		self.yMemory = d3.scale.linear().range([self.height, 0]);

		// Define the axes
		self.xAxis = d3.svg.axis().scale(self.x)
			.orient("bottom").ticks(10);

		self.yLoadAxis = d3.svg.axis().scale(self.yLoad)
		    .orient("left").ticks(10);

		self.yMemAxis = d3.svg.axis().scale(self.yMemory)
	    	.orient("right").ticks(10);
		
		// Define the lines
		self.valueline = d3.svg.line()
		    .x(function(d) { return self.x(d.time); })
		    .y(function(d) { return self.yLoad(d.load); });
		
		self.valueline2 = d3.svg.line()
	    	.x(function(d) { return self.x(d.time); })
	    	.y(function(d) { 
	    		if (d.memory)
	    			return self.yMemory(d.memory); 
	    		return self.yMemory(0);
	    	});
		    
		// Adds the svg canvas
		self.svg = d3.select("#" + self.element + "_graph")
		    .append("svg")
		        .attr("width", self.width + self.margin.left + self.margin.right)
		        .attr("height", self.height + self.margin.top + self.margin.bottom)
		    .append("g")
		        .attr("transform", 
		              "translate(" + self.margin.left + "," + self.margin.top + ")");

	    // Scale the range of the data
	    self.x.domain(d3.extent(self.data, function(d) { return d.time; }));
	    self.yLoad.domain([0, d3.max(self.data, function(d) { return d.load; })]);
	    self.yMemory.domain([0, d3.max(self.data, function(d) { return d.memory; })]);

	    // Add the valueline paths.
	    self.svg.append("path")
	        .attr("class", "line")
	        .attr("d", self.valueline(self.data))
	    	.attr("fill", "none")
	    	.attr("stroke", "#119911");
	    
	    self.svg.append("path")
        	.attr("class", "line")
        	.attr("d", self.valueline2(self.data))
        	.attr("fill", "none")
        	.attr("stroke", "blue");

	    // Add the X Axis
	    self.svg.append("g")
	        .attr("class", "x axis")
	        .attr("transform", "translate(0," + self.height + ")")
	        .call(self.xAxis);

	    // Add the Y Axes
	    self.svg.append("g")
	        .attr("class", "y axis")
	        .style("fill", "#119911")
	        .call(self.yLoadAxis);
	    
	    var yLoadLabel = self.svg.append("text")
			.attr("font-size", "1.25em")
			.attr("transform", "translate(" + self.margin.left/(-1.25) + "," + self.height/2 + ") rotate(270)")
			.text("Load");
	    
	    self.svg.append("g")
        	.attr("class", "y axis")
        	.attr("transform", "translate(" + self.width + ", 0)")
        	.style("fill", "blue")
        	.call(self.yMemAxis);
	    
	    var yMemLabel = self.svg.append("text")
			.attr("font-size", "1.25em")
			.attr("transform", "translate(" + (self.width + self.margin.right/2) + "," + (self.height/2 - 20) + ") rotate(90)")
			.text("Memory (GB)");
	    
	    // Add a line representing the number of cores the server has.
	    var numCores = 44;
	    self.svg.append("svg:line")
	    	.attr("x1", 0)
	    	.attr("x2", self.width)
	    	.attr("y1", self.yLoad(numCores))
	    	.attr("y2", self.yLoad(numCores))
	    	.attr("stroke", "green")
	    	.attr("opacity", "0.5");
	    
	    // mdb added 7/14/16 -- Add a line representing the size of RAM the server has.
	    var memSize = 500; // in GB
	    self.svg.append("svg:line")
	    	.attr("x1", 0)
	    	.attr("x2", self.width)
	    	.attr("y1", self.yMemory(memSize))
	    	.attr("y2", self.yMemory(memSize))
	    	.attr("stroke", "blue")
	    	.attr("opacity", "0.5");
	    
	    // Add brush selection
	    self.brush = d3.svg.brush()
			.x(self.x)
			.on("brushend", function() {
				self.brushed.call(self);
			});
	    
	    self.svg.append("g")
    		.attr("class", "x brush")
    		.call(self.brush)
    		.selectAll("rect")
    		.attr("y", -6)
    		.attr("height", this.height + 7);
	    
	    self.init_data_table();
	},
	init_data_table: function () {
		var self = this;
		var loadTotal = 0;
		var loadCount = 0;
		var minLoad = Number.MAX_SAFE_INTEGER;
		var maxLoad = 0;
		
		var memTotal = 0;
		var memCount = 0;
		var minMem = Number.MAX_SAFE_INTEGER;
		var maxMem = 0;
		
		for (var i = 0; i < self.data.length; i++) {
			if (self.data[i].load) {
				if (self.data[i].load < minLoad)
					minLoad = self.data[i].load;
				if (self.data[i].load > maxLoad)
					maxLoad = self.data[i].load;
				loadTotal += self.data[i].load;
				loadCount++;
			}
			
			if (self.data[i].memory) {
				if (self.data[i].memory < minMem)
					minMem = self.data[i].memory;
				if (self.data[i].memory > maxMem)
					maxMem = self.data[i].memory;
				memTotal += self.data[i].memory;
				memCount++;
			}
		}
		var avgLoad = Math.round(loadTotal/loadCount * 100)/100;
		var avgMem = Math.round(memTotal/memCount * 100)/100;
		$("#" + self.element + "_container").append("<div style='width:100%'>" +
				"<div>Average Load: " + avgLoad + "</div>" +
				"<div>Min Load: " + minLoad + "</div>" +
				"<div>Max Load: " + maxLoad + "</div>" +
				"<div><hr></div>" +
				"<div>Average Memory: " + avgMem + "</div>" +
				"<div>Min Memory: " + minMem + "</div>" +
				"<div>Max Memory: " + maxMem + "</div>" +
				"<div><hr></div>" +
				"<div>Total Data Points: " + self.data.length + "</div>" +
				"</div>"
		);
	},
	root_node: function() { // return the root node (zoomed out all the way)
	    var self = this;
        while (self.parent) {
			self = self.parent;
		}
		return self;
	},
	brushed: function() {
		var self = this;
		var extent = self.brush.extent();
		
		var newJson = this.extract_data(extent[0], extent[1]);
		if (newJson.length > 1) {
			var parent = self;
			self = new System_graph(newJson, self.element, parent);
			parent.child = self;
		}
	},
	zoom_out: function() { // zoom out as far as possible
		var self = this.root_node();
		self.initialize();
	},
	zoom_hours: function(hours) { // zoom to past number of hours
		var self = this.root_node();

		var now = new Date();
		var start = new Date();
		start.setHours(start.getHours() - hours);
		var format = d3.time.format("%H:%M:%S %m/%d/%Y");

		var newHours = self.extract_data(start, now);
		if (newHours.length > 1) {
			var parent = self;
			self = new System_graph(newHours, self.element, parent);
			parent.child = self;
		}
	},
	extract_data: function(start, end) { // extract range of data values
		var self = this;
		var newData = [];
		for (var i = 0; i < self.data.length; i++) {
			if (start <= self.data[i].time && end >= self.data[i].time)
				newData.push(self.data[i]);
		}
		return newData;
	}	
});

var Query_Counter = function(params) {
	this.elementId = params.elementId;
	this.schedule_update = params.schedule_update;
	this.cancel_update = params.cancel_update;
	this.timers = new Array();
	this.total_queries;
	this.queries_per_second = 0;
	this.updating = true;
	this.get_data();
}

$.extend(Query_Counter.prototype, {
	update: function(delay) {
		var self = this;
		if (self.updating) {
			console.log("Updating query count");
			self.cancel_update();
			self.timers['update'] = window.setTimeout(
				function() { self.get_data.call(self); },
				delay
			);
		}
	},
	get_data: function() {
		var self = this;
		$.ajax({
	        data: {
	                fname: 'get_total_queries',
	        },
	        success : function(data) {
	                if (data) {
	                		var data = JSON.parse(data);
	                		if (self.total_queries) {
	                			var delta = data.Queries - self.total_queries;
	                			//console.log(delta);
	                			self.queries_per_second = delta/5;
	                		}
	                		self.total_queries = data.Queries;
	                        
	                		$("#" + self.elementId + "_total").html("<span>Total Database Queries: " + coge.utils.numberWithCommas(self.total_queries) + "</span>");
	                		$("#" + self.elementId + "_per_second").html("<span>Queries per Second: " + self.queries_per_second + "</span>");
	                		$("#" + self.elementId + "_uptime").html("<span>Uptime: " + Math.ceil((data.Uptime)/86400) + " days</span>");
	                }
	        },
	
	        complete: function() {
	        	self.schedule_update(5000);
		}
		});
	}
});

function group_dialog(group_id) {
	$.ajax({
		url: 'User.pl',
		data: {
			fname: 'get_group_dialog',
			item_list: group_id + '_group',
		},
		success : function(data) {
			if (data)
				$('#group_dialog').html(data).dialog({width:500}).dialog('open');
		}
	});
}
