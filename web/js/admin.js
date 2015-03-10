var ITEM_TYPE_USER = 5; //TODO: This is duplicated elsewhere, move to a common location

$(function () {
	    // Configure dialogs
    $(".dialog_box").dialog({autoOpen: false, width: 500});
});

var timestamps = new Array();

function search_stuff (search_term) {
	if(search_term.length > 2) {
	//if(search_term[0].search_term.length > 2) {
		//console.log("Submit");
		//for(var i = 0; i < search_term.length; i++) {
			//console.log(search_term[i]);
		//}
		//console.log(search_term);

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
				//console.log("data received");
                                //console.log(data);
                                var obj = jQuery.parseJSON(data);
				//console.log(obj.items[0]);

				var userCounter = 0, orgCounter = 0, genCounter = 0, expCounter = 0, noteCounter = 0, usrgroupCounter = 0;
				var userList = "", orgList = "", genList = "", expList = "", noteList = "", usrgroupList = "";

				for (var i = 0; i < obj.items.length; i++) {

					if (obj && obj.items && obj.timestamp == timestamps['search_stuff']) {

		                                if (obj.items[i].type == "user") {
	                                	        userList = userList + "<tr><td><span>";
							userList = userList + (obj.items[i].label) + " (ID: " + (obj.items[i].id) + ") ";
							userList = userList + "<button onclick=\"search_user(" + (obj.items[i].id) + ",'user')\">Show Data</button>";
							userList = userList + "</span></td></tr>";
	                                        	userCounter++;
		              	                }

						if (obj.items[i].type == "organism") {
                                           		orgList = orgList + "<tr><td><span>" + (obj.items[i].label) + " (ID: " + (obj.items[i].id) + ")" + "</span></td></tr>";
                                           		orgCounter++;
						}

						if (obj.items[i].type == "genome") {
							genList = genList + "<tr><td>";
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
                                                        expList = expList + "<tr><td>";
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
                                                        noteList = noteList + "<tr><td>";
							if (obj.items[i].deleted == 1) {
                                                                noteList = noteList + "<span style=\"color: red\">";
                                                        } else {
                                                                noteList = noteList + "<span>";
                                                        }
							noteList = noteList + (obj.items[i].label) + " (ID: " + (obj.items[i].id) + ") <a href=\"NotebookView.pl?lid=" + (obj.items[i].id) + "\">Info </a>";
							noteList = noteList + "<button onclick='share_dialog(" + obj.items[i].id + ", 1 )'>Edit Access</button>";
							NoteList = noteList + "</span></td></tr>";
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
				}

				//Populate the html with the results
				$(".result").fadeIn( 'fast');

				//user
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

				//organism
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

				//genome
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

				//experiment
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

				//notebook
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

				//user group
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

                       	},
                });
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
	console.log("Open_dialog called");
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
	console.log("" + target_item + " " + id + " " + type);
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
var previous = 0;
var previous_type = 0;

// indicates the part of the ADmin page currently being displayed. 0 --> "master", 1 --> "user info"
var current_page = 0;

function search_user(userID, search_type) {
	if (current_page == 0) {	
		toggle_master();
		//current_page=1;
	}
	if(previous != userID) {
		//$('#userResults').hide();
		$('#userResults').html("Loading...");
		console.log("Before search.");
		user_info(userID, search_type);
	}
	previous = userID;
	previous_type = search_type;
}

function refresh_data() {
	$('#userResults').html("Loading...");
	user_info(previous, previous_type);
}

function user_info(userID, search_type) {

	var search_term = userID;
	//console.log(search_term);
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

				//for (var j = 0; j < obj.items[i].result.length; j++) {
				//	console.log(current);
				//}				

				
				//for each object belonging to that user, populate tables
				for (var j = 0; j < obj.items[i].result.length; j++) {
	                        	if (obj && obj.items && obj.timestamp == timestamps['user_info']) {
						
						var current = obj.items[i].result[j];
	
	                                        if (current.type == "genome") {
	                                                genList = genList + "<tr><td>";
	                                               	
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
                                                	expList = expList + "<tr><td>";

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
                                                	noteList = noteList + "<tr><td>";

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
				nameBlock = nameBlock + obj.items[i].user + " (ID: " + obj.items[i].user_id + "):</span></div>";
				

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
			
			console.log("After search.");
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
		}
	});
}

/*function search_users (search_term) {
        if (search_term.length > 2) {
		search_term = search_term + " type:user"
                timestamps['search_organisms'] = new Date().getTime();
                $.ajax({
                        data: {
                                fname: 'search_stuff',
                                search_term: search_term,
                                timestamp: timestamps['search_stuff']
                        },
                        success : function(data) {
                                var obj = jQuery.parseJSON(data);
                                if (obj && obj.items && obj.timestamp == timestamps['search_stuff']) {
                                        $("#user_field").autocomplete({source: obj.items}).autocomplete("search");
                                }
                        },
                });
        }
}*/

function share_dialog(id, type) {
	var item_list = "content_" + id + "_" + type;  //selected.map(function(){return this.parentNode.id;}).get().join(',');
	console.log(item_list);
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
	timestamps['search_share'] = new Date().getTime()

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

/*function update_dialog(request, user, identifier, formatter) {
    var get_status = function () {
        $.ajax({
            type: 'GET',
            url: request,
            dataType: 'json',
            data: {
                username: user
            },
            success: update_callback,
            error: update_callback,
            xhrFields: {
                withCredentials: true
            }
        });
    };

    var update_callback = function(json) {
        var dialog = $(identifier);
        var workflow_status = $("<p></p>");
        var data = $("<ul></ul>");
        var results = [];
        var current_status;
        var timeout = 2000;

        var callback = function() {
            update_dialog(request, user, identifier, formatter);
        }

        if (json.error) {
            pageObj.error++;
            if (pageObj.error > 3) {
                workflow_status.html('<span class=\"alert\">The job engine has failed.</span>');
                load_failed();
                return;
            }
        } else {
            pageObj.error = 0;
        }

        if (json.status) {
            current_status = json.status.toLowerCase();
            workflow_status.html("Workflow status: ");
            workflow_status.append($('<span></span>').html(json.status));
            workflow_status.addClass('bold');
        } else {
            setTimeout(callback, timeout);
            return;
        }

        if (json.tasks) {
            var jobs = json.tasks;
            for (var index = 0; index < jobs.length; index++) {
                var item = formatter(jobs[index]);
                if (item) {
                    results.push(item);
                }
            }
        }

        if (!dialog.dialog('isOpen')) {
            return;
        }

        //FIXME Update when a workflow supports elapsed time
        if (current_status == "completed") {
            var total = json.tasks.reduce(function(a, b) {
                if (!b.elapsed) return a;

                return a + b.elapsed;
            }, 0);

            var duration = coge.utils.toPrettyDuration(total);

            workflow_status.append("<br>Finished in " + duration);
            workflow_status.find('span').addClass('completed');
            get_load_log(function(result) {
                load_succeeded(result);
            });

        }
        else if (current_status == "failed"
                || current_status == "error"
                || current_status == "terminated"
                || current_status == "cancelled")
        {
            workflow_status.find('span').addClass('alert');
            get_load_log(function(result) {
                load_failed(result);
            });
        }
        else if (current_status == "notfound") {
            setTimeout(callback, timeout);
            return;
        }
        else {
            workflow_status.find('span').addClass('running');
            setTimeout(callback, timeout);
        }

        results.push(workflow_status);
        data.append(results);
        dialog.find('#load_log').html(data);
    };

    get_status();
}*/

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

