//Global Variables
var user_is_admin = false;

$(function () {
	// Initialize CoGe web services
    coge.services.init({
    	baseUrl: API_BASE_URL,
    	userName: USER_NAME
    });
    
    // See if the current user is an admin
    $.ajax({
		data: {
			fname: 'user_is_admin',
		},
		success: function(data) {
			if (data == 1)
				user_is_admin = true;
			else
				user_is_admin = false;
		}
	}).done(function() {
		search_stuff(SEARCH_TERM);
	});
});

function search_stuff(search_term) {
	if (!search_term || search_term.length <= 2) {
		$("#noresult").html('Please specify a search term longer than 2 characters').show();
		$("#loading").hide();
		return;
	}
	
	$("#noresult").hide();
	$("#loading").show();
	
	coge.services.search_global(search_term)
		.done(function(response) {
			//console.log(response);
			var obj = response;
			
			if (!obj || !obj.results)
				return;
			
			var userCounter = 0, orgCounter = 0, genCounter = 0, expCounter = 0, noteCounter = 0, usrgroupCounter = 0;
			var userList = "", orgList = "", genList = "", expList = "", noteList = "", usrgroupList = "";

			for (var i = 0; i < obj.results.length; i++) {
				/*if (obj.results[i].type == "user") {
					userList = userList + "<tr><td><span>";
					userList = userList + obj.results[i].first + " " + obj.results[i].last + ", ";
					if (user_is_admin) {
						userList = userList + obj.results[i].email;
					}
					userList = userList + " (" + obj.results[i].username + ", id" + obj.results[i].id + ") ";
					userList = userList + "<button onclick=\"search_user(" + obj.results[i].id + ",'user')\">Show Data</button>";
					userList = userList + "</span></td></tr>";
					userCounter++;
				}*/

				if (obj.results[i].type == "organism") {
					orgList = orgList + "<tr><td><span title='" + obj.results[i].description + "'>";
					orgList = orgList + obj.results[i].name + " (id" + obj.results[i].id + ")";
					orgList = orgList + " <a href=\"OrganismView.pl?oid=" + obj.results[i].id + "\">Info </a>";
					orgList = orgList + "</span></td></tr>";
					orgCounter++;
				}

				if (obj.results[i].type == "genome") {
					if(obj.results[i].deleted != 1 || user_is_admin) {
						//genList = genList + "<tr><td><span onclick=\"modify_item(" + obj.results[i].id + ", 'Genome', 'restrict');\"";
						genList = genList + "<tr><td><span";
						/*if (obj.results[i].restricted == 1) {
							genList = genList + " class=\"link ui-icon ui-icon-locked\"></span>";
						} else {
							genList = genList + " class=\"link ui-icon ui-icon-unlocked\"></span>";
						}
						genList = genList + "<span onclick=\"modify_item(" + obj.results[i].id + ", 'Genome', 'delete');\"";
						genList = genList + " class=\"link ui-icon ui-icon-trash\"></span>";*/
						if (obj.results[i].deleted == 1) {
							genList = genList + "<span style=\"color: red\">";
						} else {
							genList = genList + "<span>";
						}
						genList = genList + obj.results[i].name + " <a href=\"GenomeInfo.pl?gid=" + obj.results[i].id + "\">Info </a>";
						//genList = genList + "<button class='access' onclick='share_dialog(" + obj.results[i].id + ", 2, " + obj.results[i].restricted + ")'>Edit Access</button>";
						genList = genList + "</span></td></tr>";
						genCounter++;
					}
				}

				if (obj.results[i].type == "experiment") {
					if(obj.results[i].deleted != 1 || user_is_admin) {
						//expList = expList + "<tr><td><span onclick=\"modify_item(" + obj.results[i].id + ", 'Experiment', 'restrict');\"";
						expList = expList + "<tr><td><span";
						/*if (obj.results[i].restricted == 1) {
							expList = expList + " class=\"link ui-icon ui-icon-locked\"></span>";
						} else {
							expList = expList + " class=\"link ui-icon ui-icon-unlocked\"></span>";
						}
						expList = expList + "<span onclick=\"modify_item(" + obj.results[i].id + ", 'Experiment', 'delete');\"";
						expList = expList + " class=\"link ui-icon ui-icon-trash\"></span>";*/
						if (obj.results[i].deleted == 1) {
							expList = expList + "<span style=\"color: red\">";
						} else {
							expList = expList + "<span>";
						}
						expList = expList + obj.results[i].name + " (id" + obj.results[i].id + ") <a href=\"ExperimentView.pl?eid=" + obj.results[i].id + "\">Info </a>";
						//expList = expList + "<button class='access' onclick='share_dialog(" + obj.results[i].id + ", 3, " + obj.results[i].restricted + ")'>Edit Access</button>";
						expList = expList + "</span></td></tr>";
						expCounter++;
					}
				}

				if (obj.results[i].type == "notebook") {
					if(obj.results[i].deleted != 1 || user_is_admin) {
						//noteList = noteList + "<tr><td><span onclick=\"modify_item(" + obj.results[i].id + ", 'List', 'restrict');\"";
						noteList = noteList + "<tr><td><span";
						/*if (obj.results[i].restricted == 1) {
							noteList = noteList + " class=\"link ui-icon ui-icon-locked\"></span>";
						} else {
							noteList = noteList + " class=\"link ui-icon ui-icon-unlocked\"></span>";
						}
						noteList = noteList + "<span onclick=\"modify_item(" + obj.results[i].id + ", 'List', 'delete');\"";
						noteList = noteList + " class=\"link ui-icon ui-icon-trash\"></span>";*/
						if (obj.results[i].deleted == 1) {
							noteList = noteList + "<span style=\"color: red\">";
						} else {
							noteList = noteList + "<span>";
						}
						noteList = noteList + obj.results[i].name + " <a href=\"NotebookView.pl?lid=" + obj.results[i].id + "\">Info </a>";
						//noteList = noteList + "<button class='access' onclick='share_dialog(" + obj.results[i].id + ", 1 , " + obj.results[i].restricted + ")'>Edit Access</button>";
						noteList = noteList + "</span></td></tr>";
						noteCounter++;
					}
				}
						
				if (obj.results[i].type == "user_group") {
					if(obj.results[i].deleted != 1 || user_is_admin) {
						usrgroupList = usrgroupList + "<tr><td>";
						if (obj.results[i].deleted == 1) {
							usrgroupList = usrgroupList + "<span style=\"color: red\">";
						} else {
							usrgroupList = usrgroupList + "<span>";
						}
						usrgroupList = usrgroupList + obj.results[i].name + " (id" + obj.results[i].id + ") ";
						//usrgroupList = usrgroupList + "<button onclick=\"search_user(" + obj.results[i].id + ",'group')\">Show Data</button>";
						usrgroupList = usrgroupList + "</span></td></tr>";
						usrgroupCounter++;
					}
				}
			}
			

			//Populate the html with the results
			$("#loading").show();
			$(".result").fadeIn( 'fast');
			
			if (userCounter + orgCounter + genCounter + expCounter + noteCounter + usrgroupCounter == 0)
				$('#noresult').html('No matching results found').show();
			
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
				if (usrgroupCounter <= 10) {
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
			$("#masterTable").show();
			if (!user_is_admin)
				$(".access").hide();
			else
				$(".access").show();
		})
		.fail(function() {
			$("#loading").hide();
			$("#masterTable").html("An error occured. Please reload the page and try again.");
		});
	
	//previous_search = search_term;
}

function toggle_arrow(id) {
	if ( $(id).find('img').attr('src') == "picts/arrow-right-icon.png" ) {
        $(id).find('img').attr("src", "picts/arrow-down-icon.png");
    }
	else {
		$(id).find('img').attr("src", "picts/arrow-right-icon.png");
	}
}

function show_table(id) {
	$(id).fadeToggle('fast');	
}