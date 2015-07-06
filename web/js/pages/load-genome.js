function file_selected(filename, url) {
	$('#no_files').hide();
	$('#file_select_dialog').dialog('close');
	$('#files_clear').removeClass('ui-state-disabled');
}

function get_sequence_types(id) {
	$.ajax({
		data: {
			fname: 'get_sequence_types'
		},
		success : function(html) {
			$('#select_type').html(html);
			if (id) {
				$('#select_type').val(id);
			}
		}
	});
}

function create_sequence_type() {
	var name = $('#edit_type_name').val();
	var desc = $('#edit_type_desc').val();
	$.ajax({
		data: {
			fname: 'create_sequence_type',
			name: name,
			desc: desc,
		},
		success : function(id) {
			$('#create_new_type_dialog').dialog('close');
			if (id) {
				get_sequence_types(id);
			}
		}
	});
}

function create_source() {
	var name = $('#edit_source_name').val();
	var desc = $('#edit_source_desc').val();
	var link = $('#edit_source_link').val();
	$.ajax({
		data: {
			fname: 'create_source',
			name: name,
			desc: desc,
			link: link,
		},
		success : function(name) {
			$('#create_new_source_dialog').dialog('close');
			if (name) {
				$('#edit_source').val(name);
			}
		}
	});
}

function create_organism() {
	if ($('#create_organism_button').hasClass('ui-state-disabled')) {
		return;
	}

	var name = $('#edit_organism_name').val();
	if (!name) {
		alert('Please specify the organism name.');
		return;
	}

	var desc = $('#edit_organism_desc').val();
	if (!desc) {
		alert('Please specify the organism description.');
		return;
	}
	if (desc.split(';').length < 2) {
		alert('Please specify an NCBI taxonomy of classes separated\nby semicolons for the organism description.');
		return;
	}

	$.ajax({
		data: {
			fname: 'create_organism',
			name: name,
			desc: desc,
		},
		success : function(organism_id) {
			if (organism_id) {
				$('#create_new_organism_dialog').dialog('close');
				$('#edit_organism').val(name).data("organism_id", organism_id);
			}
		}
	});
}

function reset_log() {
	$('#load_log').html('');
	$('#loading_msg').show();
	$('#finished_msg,#error_msg,#ok_button,#logfile').hide();
	$('#finish_actions,#cancel_button').hide();
}

function check_login() {
	var logged_in = false;

	$.ajax({
		async: false,
		data: {
			fname: 'check_login',
		},
		success : function(rc) {
			logged_in = rc;
		}
	});

	return logged_in;
}

function error_help(s) {
	$('#error_help_text')
		.html(s)
		.show()
		.delay(10*1000)
		.fadeOut(1500);
}

function load_genome() {
	// Validate data items
	var items = get_selected_files();
	if (items == null) {
		error_help('Files are still being transferred, please wait.');
		return;
	}
	else if (items.length == 0) {
		error_help('Please select some sequence files by clicking <b>Add Data</b>.');
		return;
	}

	// Prevent mix of NCBI and file data types
	var types = {};
	items.forEach(function(item) { types[item.type] = 1; });
	var isNCBI = 'ncbi' in types;
	if (Object.keys(types).length > 1 && isNCBI) {
		error_help('Cannot mix NCBI data with other types.');
		return;
	}

	// Validate other input fields
	var organism_name = $('#edit_organism').val();
	var organism_id = $('#edit_organism').data("organism_id");
	if (!isNCBI && (!organism_id || !organism_name)) {
		error_help('Organism not found.');
		return;
	}

	var version = $('#edit_version').val();
	if (!isNCBI && !version) {
		error_help('Please specify a genome version.');
		return;
	}

	var source = $('#edit_source').val();
	if (!isNCBI && !source) {
		error_help('Please specify a data source.');
		return;
	}

	var user_name = $('#edit_user').val(); // only exists if admin
	var keep_headers = $('#keep_headers').is(':checked'); // only exists if admin
	var name = $('#edit_name').val();
	var description = $('#edit_description').val();
	var link = $('#edit_link').val();
	var type_id = $('#select_type').val();
	var restricted = $('#restricted').is(':checked');
	var json = JSON.stringify(items);

	// Prevent concurrent executions - issue 101
	if ( $("#load_dialog").dialog( "isOpen" ) )
		return;

	// Make sure user is still logged-in - issue 206
	if (!check_login()) {
		alert('Your session has expired, please log in again.');
		location.reload(true)
		return;
	}

	// Open status dialog right away - issue 101
	reset_log();
	$('#load_dialog').dialog('open');
	$('#load_log').html('Initializing ...');
	newLoad = true;

	$.ajax({
		data: {
			fname: 'load_genome',
			load_id: load_id,
			name: name,
			description: description,
			link: link,
			version: version,
			type_id: type_id,
			restricted: restricted,
			organism_id: organism_id,
			source_name: source,
			user_name: user_name,
			keep_headers: keep_headers,
			items: json,
			timestamp: new Date().getTime()
		},
		success : function(data) {
			var obj = jQuery.parseJSON(data);
			if (!obj || obj.error) {
				if (!obj)
					alert("Error: load_genome: invalid response from server");
				else
					alert(obj.error);
				reset_load();
				return;
			}

			// Set link in status dialog
			$('#loading_msg span a').attr('href', obj.link).html(obj.link);

			// Start status update
            if (obj.job_id) { // JEX status for load FASTQ
                job_id = obj.job_id;
                window.history.pushState({}, "Title", PAGE_NAME + "?job_id=" + obj.job_id); // Add job_id to browser URL
                update_dialog(STATUS_URL + obj.job_id, pageObj.user, "#load_dialog", progress_formatter);
            }
		}
		// TODO: handle error, show in status dialog
	});
}

function get_load_log(callback) {
    $.ajax({
        data: {
            dataType:    'text',
            fname:       'get_load_log',
            workflow_id: job_id,
            timestamp:   new Date().getTime()
        },
        success : function(data) {
            if (callback) {
            	var obj = jQuery.parseJSON(data);
                callback(obj);
                return;
            }
        }
    });
}

function load_failed(obj) {
	// Handle special case of genbank load of existing genome
	if ( obj && obj.links && obj.links.length ) {
		var link_text = obj.links.reduce(function(prev, cur) {
				return prev + '<a href="' + cur + '" target=_new>' + cur + '</a>' + '<br>';
			},
			'<b>Cannot load because the data already exist in the system:</b><br>');
		var log = $('#load_log');
		log.html( log.html() + link_text );
	}
	else { // mdb added 6/24/14 - temporary message until JEX logging is improved
		var msg =
			'<div class="alert">' +
			'The CoGe Support Team has been notified of this error but please ' + 
			'feel free to contact us at <a href="mailto:<TMPL_VAR NAME=SUPPORT_EMAIL>"><TMPL_VAR NAME=SUPPORT_EMAIL></a> ' +
			'and we can help to determine the cause.' +
			'</div>';
		var log = $('#load_log');

        if (obj) {
            $("#logfile a").attr("href", obj);
            $('#logfile').fadeIn();
        }

		log.html( log.html() + msg );
	}

    // Update dialog
    $('#loading_msg').hide();
    $('#error_msg').fadeIn();
    $('#cancel_button').fadeIn();

    if (newLoad) { // mdb added check to prevent redundant emails, 8/14/14 issue 458
	    $.ajax({
	        data: {
	            fname: "send_error_report",
	            load_id: load_id,
	            job_id: job_id
	        }
	    });
    }
}

function load_succeeded(obj) {
    // Update globals
    genome_id = obj.genome_id; // for continuing to GenomeInfo

    // Update dialog
    $('#loading_msg').hide();
    $('#finished_msg,#finish_actions,#ok_button').fadeIn();
}

function reset_load() {
	clear_list();
    window.history.pushState({}, "Title", PAGE_NAME);
    $('#load_dialog').dialog('close');
}

function handle_action() {
    var action = $("#finish_actions select").val() || "genome";

    if (action === "genome") {
	    window.location.href = "GenomeInfo.pl?embed=" + embed + "&gid=" + genome_id;
    } else if (action === "annotation") {
	    window.location.href = "LoadAnnotation.pl?embed=" + embed + "&gid=" + genome_id;
    } else if (action === "new") {
        reset_load();
    }
}

function wait_to_search (search_func, search_term) {
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

function search_organisms (search_term) {
	if (search_term.length > 2) {
		timestamps['search_organisms'] = new Date().getTime();
		$.ajax({
			data: {
				fname: 'search_organisms',
				search_term: search_term,
				timestamp: timestamps['search_organisms']
			},
			success : function(data) {
				var obj = jQuery.parseJSON(data);
				if (obj && obj.items && obj.timestamp == timestamps['search_organisms']) {
					$("#edit_organism").autocomplete({source: obj.items}).autocomplete("search");
				}
			},
		});
	}
}

function search_users (search_term) {
	timestamps['search_users'] = new Date().getTime();
	$.ajax({
		data: {
			fname: 'search_users',
			search_term: search_term,
			timestamp: timestamps['search_users']
		},
		success : function(data) {
			var obj = jQuery.parseJSON(data);
			if (obj && obj.items && obj.timestamp == timestamps['search_users']) {
				$("#edit_user").autocomplete({source: obj.items});
				$("#edit_user").autocomplete("search");
			}
		},
	});
}

function get_sources () {
	$.ajax({
		data: {
			fname: 'get_sources',
		},
		success : function(data) {
			var obj = jQuery.parseJSON(data);
			if (obj) {
				$("#edit_source").autocomplete({source: obj});
			}
		},
	});
}

function more() {
	$('#more').hide();
	$('#edit_name,#edit_description,#edit_link').closest('tr').fadeIn();
}

function build_taxonomy_tree(items) {
	if (items) {
		var tree = $("#tax_tree");
		if (items.length == 0) {
			tree.slideUp('fast',
				function() {
					tree.empty();
					$('#tax_empty').fadeIn();
					$('#edit_organism_name').val('');
					$('#edit_organism_desc').val('');
					activate_on_input(['edit_organism_name', 'edit_organism_desc'], 'create_organism_button');
				}
			);
		}
		else {
			var list = $('<ul></ul>');
			items.sort(sort_name).forEach(
				function(e) {
					var id = 'taxon' + e.id;
					$(list)
						.append('<li id="' + id + '" name="' + e.name + '">' +
								'<a href="#">' + e.name + '</a></li>');
					//tree.jstree("create_node", null, "last", {attr: {id: id}, data: e.name})
				}
			);
			tree.empty()
				.append(list)
				.jstree()
					.bind("select_node.jstree",
				      	function (event, data) {
							var id = data.rslt.obj.attr("id");
							taxonomy_get_node(id);
						}
					)
				.slideDown();
		}
	}
	else {
		$('#ncbi_result').html('<i>No result</i>');
		$('#create_organism_button').addClass('ui-state-disabled');
	}
}

//TODO break this out into a widget/plugin
function search_ncbi_taxonomy(search_term) {
	if (!search_term || search_term.length < 3) {
		return;
	}

	$("#wait_ncbi").animate({opacity:1});
	$('#tax_empty').fadeOut();
	$("#edit_organism_info").slideUp();

	$.get(
		ENTREZ_URL + "esearch.fcgi?db=taxonomy&term=" + search_term + "*",
		function(xml) {
			var ids;
			$(xml).find("Id").each(
				function() {
					ids += $(this).text() + ',';
				}
			);

			if (ids) {
				$.get(
					ENTREZ_URL + "efetch.fcgi?db=taxonomy&id=" + ids,
					function(xml) {
						var results = new Array();

						$(xml).children("Taxon").each(
							function() {
								var id = $(this).children('TaxId').text();
								var name = $(this).children('ScientificName').text();
								var lineage = $(this).children('Lineage').text();
								results.push({id: id, name: name, lineage: lineage});
							}
						);

						build_taxonomy_tree(results);

						$("#wait_ncbi").animate({opacity:0});
					}
				);
			}
			else {
				var tree = $("#tax_tree");
				tree.slideUp('fast',
					function() {
						tree.empty();
						$('#tax_empty').fadeIn();
						$('#edit_organism_name').val('');
						$('#edit_organism_desc').val('');
						activate_on_input(['edit_organism_name', 'edit_organism_desc'], 'create_organism_button');
					}
				);

				$("#wait_ncbi").animate({opacity:0});
			}
		}
	);
}

function taxonomy_get_node(id) { // FIXME: cleanup and merge common stuff with search_ncbi_taxonomy
	var tree = $("#tax_tree");
	var node = $("#"+id);
	var name = $(node).attr('name');

	// If the node is open then it has children so just
	// close the node and return.
	var is_open = tree.jstree("is_open", node);
	if (is_open) {
		tree.jstree("close_node", node);
		return;
	}

	// If not a leaf node then we already retreived the children so
	// just open the node and return.
	var is_leaf = tree.jstree("is_leaf", node);
	if (!is_leaf) {
		tree.jstree("open_node", node);
		return;
	}

	// Retrieve the children for this taxon.
	$("#wait_ncbi").animate({opacity:1});

	$.get(
		ENTREZ_URL + "esearch.fcgi?db=taxonomy&term=" + name + "&field=nxlv",
		function(xml) {
			var ids;
			$(xml).find("Id").each(
				function() {
					ids += $(this).text() + ',';
				}
			);

			if (ids) {
				$.get(
					ENTREZ_URL + "efetch.fcgi?db=taxonomy&id=" + ids,
					function(xml) {
						var results = new Array();

						$(xml).children("Taxon").each(
							function() {
								var id = $(this).children('TaxId').text();
								var name = $(this).children('ScientificName').text();
								var lineage = $(this).children('Lineage').text();
								results.push({id: id, name: name, lineage: lineage});
							}
						);

						if (results.length) {
							results.sort(sort_name).forEach(
								function(e) {
									tree.jstree("create_node", node, "last", {attr: {id: 'taxon'+e.id, name: e.name}, data: e.name})
										.data("lineage", e.lineage);
								}
							);
							tree.jstree("open_node", node);
						}

						$("#wait_ncbi").animate({opacity:0});
					}
				);
			}
			else {
				var lineage = node.data("lineage");
				$("#edit_organism_name").val(name);
				$("#edit_organism_desc").val(lineage);
				activate_on_input(['edit_organism_name', 'edit_organism_desc'], 'create_organism_button');
				$("#edit_organism_info").slideDown();

				$("#wait_ncbi").animate({opacity:0});
			}
		}
	);

}

function sort_name (a,b) {
	var nameA=a.name.toLowerCase(), nameB=b.name.toLowerCase()
	if (nameA < nameB) //sort string ascending
		return -1
	if (nameA > nameB)
		return 1
	return 0 //default return value (no sorting)
}

function progress_formatter(item) {
    var msg;
    var row = $('<li>'+ item.description + ' </li>');

    var job_status = $('<span></span>');

    if (item.status == 'scheduled')
        job_status.append(item.status).addClass('down bold');
    else if (item.status == 'completed')
        job_status.append(item.status).addClass('completed bold');
    else if (item.status == 'running')
        job_status.append(item.status).addClass('running bold');
    else if (item.status == 'skipped')
        job_status.append("already generated").addClass('skipped bold');
    else if (item.status == 'cancelled')
        job_status.append(item.status).addClass('alert bold');
    else if (item.status == 'failed')
        job_status.append(item.status).addClass('alert bold');
    else
        return;

    row.append(job_status);

    if (item.elapsed)  {
        row.append(" in " + coge.utils.toPrettyDuration(item.elapsed));
    }

    if (item.log) {
        var p = item.log.split("\n");

        var pElements = p.map(function(item) {
            var norm = item.replace(/\\t/g, " ").replace(/\\'/g, "'");
            return $("<div></div>").append(norm);
        });

        var log = $("<div></div>").html(pElements).addClass("padded");
        row.append(log);
    }

    return row;
}

function update_dialog(request, user, identifier, formatter) {
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
                var logfile;

                if (json.results && json.results.length) {
                    logfile = json.results[0].path;
                }
                load_failed(logfile);
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
            var logfile;

            if (json.results && json.results.length) {
                logfile = json.results[0].path;
            }
            load_failed(logfile);
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
}