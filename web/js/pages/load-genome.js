// Global genome data
var current_genome = {};

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
				$('#edit_organism')
					.val(name)
					.data("organism_id", organism_id);
			}
		}
	});
}

function load_genome(genome) {
	var keep_headers = $('#keep_headers').val();
	
	// Prevent mix of NCBI and file data types
//	var types = {};
//	items.forEach(function(item) { types[item.type] = 1; });
//	var isNCBI = 'ncbi' in types;
//	if (Object.keys(types).length > 1 && isNCBI) {
//		error_help('Cannot mix NCBI data with other types.');
//		return;
//	}
	
	// Build file list
	var items = genome.data.map(function(item) {
		return { path: item.path, type: item.type };
	});

	// Open progress window
	coge.progress.begin();
	newLoad = true;
	
	console.log(genome);

	$.ajax({ // TODO migrate to web services like load-experiment.js
		dataType: 'json',
		data: {
			fname: 'load_genome',
			load_id: load_id,
			name: genome.metadata.name,
			description: genome.metadata.description,
			link: genome.metadata.link,
			version: genome.metadata.version,
			type_id: genome.metadata.type,
			restricted: genome.metadata.restricted,
			organism_id: genome.organism_id,
			source_name: genome.metadata.source,
			user_name: USER_NAME,
			keep_headers: keep_headers,
			items: JSON.stringify(items),
			timestamp: new Date().getTime()
		},
		success : function(response) {
			// Handle reponse error
			if (!response || response.error) {
				if (!response)
					coge.progress.failed("Error: load_genome: invalid response from server");
				else
					coge.progress.failed(response.error);
				return;
			}

			// Start status update
            if (response.job_id) { // JEX status for load FASTQ
                window.history.pushState({}, "Title", PAGE_NAME + "?job_id=" + response.job_id); // Add job_id to browser URL
                coge.progress.update(response.job_id, response.link);
            }
		},
		error: function(jqXHR, textStatus, errorThrown) { // transaction error callback
	    	coge.progress.failed("Couldn't talk to the server: " + textStatus + ': ' + errorThrown);
	    }
	});
}

//function get_load_log(callback) {
//    $.ajax({
//        data: {
//            dataType:    'text',
//            fname:       'get_load_log',
//            workflow_id: job_id,
//            timestamp:   new Date().getTime()
//        },
//        success : function(data) {
//            if (callback) {
//            	var obj = jQuery.parseJSON(data);
//                callback(obj);
//                return;
//            }
//        }
//    });
//}
//
//function load_failed(obj) {
//	// Handle special case of genbank load of existing genome
//	if ( obj && obj.links && obj.links.length ) {
//		var link_text = obj.links.reduce(function(prev, cur) {
//				return prev + '<a href="' + cur + '" target=_new>' + cur + '</a>' + '<br>';
//			},
//			'<b>Cannot load because the data already exist in the system:</b><br>');
//		var log = $('#load_log');
//		log.html( log.html() + link_text );
//	}
//	else { // mdb added 6/24/14 - temporary message until JEX logging is improved
//		var msg =
//			'<div class="alert">' +
//			'The CoGe Support Team has been notified of this error but please ' + 
//			'feel free to contact us at <a href="mailto:<TMPL_VAR NAME=SUPPORT_EMAIL>"><TMPL_VAR NAME=SUPPORT_EMAIL></a> ' +
//			'and we can help to determine the cause.' +
//			'</div>';
//		var log = $('#load_log');
//
//        if (obj) {
//            $("#logfile a").attr("href", obj);
//            $('#logfile').fadeIn();
//        }
//
//		log.html( log.html() + msg );
//	}
//
//    // Update dialog
//    $('#loading_msg').hide();
//    $('#error_msg').fadeIn();
//    $('#cancel_button').fadeIn();
//
//    if (newLoad) { // mdb added check to prevent redundant emails, 8/14/14 issue 458
//	    $.ajax({
//	        data: {
//	            fname: "send_error_report",
//	            load_id: load_id,
//	            job_id: job_id
//	        }
//	    });
//    }
//}
//
//function load_succeeded(obj) {
//    // Update globals
//    genome_id = obj.genome_id; // for continuing to GenomeInfo
//
//    // Update dialog
//    $('#loading_msg').hide();
//    $('#finished_msg,#finish_actions,#ok_button').fadeIn();
//}
//
//function reset_load() {
//	clear_list();
//    window.history.pushState({}, "Title", PAGE_NAME);
//    $('#load_dialog').dialog('close');
//}

function handle_action() {
    var action = $("#finish_actions select").val() || "genome";

    if (action === "genome") {
	    window.location.href = "GenomeInfo.pl?embed=" + EMBED + "&gid=" + genome_id;
    } 
    else if (action === "annotation") {
	    window.location.href = "LoadAnnotation.pl?embed=" + EMBED + "&gid=" + genome_id;
    } 
    else if (action === "new") {
        reset_load();
    }
}

function wait_to_search (search_func, search_obj) {
	if (pageObj.time)
		clearTimeout(pageObj.time);

	// FIXME: could generalize by passing select id instead of separate search_* functions
	pageObj.time = setTimeout(
		function() {
			search_func(search_obj.value);
		},
		500
	);
}

function search_organisms (search_term) {
	if (search_term && search_term.length > 2) {
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
					$("#edit_organism")
						.autocomplete({source: obj.items})
						.autocomplete("search");
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
				$("#edit_user")
					.autocomplete({source: obj.items})
					.autocomplete("search");
			}
		},
	});
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

function GenomeDescriptionView(opts) {
    this.genome = opts.genome;
    this.metadata = opts.metadata;
    this.organism_id = opts.organism_id;
    this.onError = opts.onError;
    this.sources = undefined;
    this.title = "Describe your genome";
    this.initialize();
}

$.extend(GenomeDescriptionView.prototype, {
    initialize: function() {
        this.el = $($("#description-template").html());
        this.edit_source = this.el.find("#edit_source");
        this.edit_organism = this.el.find("#edit_organism");
        this.select_type = this.el.find('#select_type');
        this.get_sources();
        this.get_sequence_types();

        if (this.metadata) {
            this.el.find('#edit_name').val(this.metadata.name);
            this.el.find('#edit_description').val(this.metadata.description);
            this.el.find('#edit_version').val(this.metadata.version);
            this.select_type.val(this.metadata.type);
            this.edit_source.val(this.metadata.source);

            if (!this.metadata.restricted)
                this.el.find('#restricted').removeAttr('checked');

            this.edit_organism.val(this.metadata.organism);
        }
    },

    get_sources: function() { // TODO move into web services
        var self = this;

        return $.ajax({
            dataType: "json",
            data: {
                fname: 'get_sources',
            },
            success : function(sources) {
                self.sources = sources;
                if (sources)
                    self.edit_source.autocomplete({source: sources});
            },
        });
    },
    
    get_sequence_types: function(id) { // TODO move into web services
    	var self = this;
    	
    	$.ajax({
    		dataType: "html",
    		data: {
    			fname: 'get_sequence_types'
    		},
    		success : function(html) {
    			self.select_type.html(html);
    			if (id)
    				self.select_type.val(id);
    		}
    	});
    },

    render: function() {
        var self = this;

        // jQuery Events
        this.edit_organism.unbind().change(function() {
            // Reset organism_id when item has changed
            self.organism_id = undefined;
        });

        // jQuery UI
        this.edit_organism.autocomplete({
            source:[],
            select: function(event, ui) {
                $(this).val(ui.item.label);
                self.organism_id = ui.item.value;
                return false; // Prevent the widget from inserting the value.
            },

            focus: function(event, ui) {
                return false; // Prevent the widget from inserting the value.
            }
        });

        this.edit_source.autocomplete({source: this.sources});
    },

    is_valid: function() {
        var name = this.el.find('#edit_name').val();
        var description = this.el.find('#edit_description').val();
        var link = this.el.find('#edit_link').val();
        var version = this.el.find('#edit_version').val();
        var type = this.el.find('#select_type').val();
        var restricted = this.el.find('#restricted').is(':checked');
        var organism = this.el.find('#edit_organism').val();

        if (!organism || organism === 'Search' || !this.organism_id) {
        	if (this.onError)
            	this.onError('Please specify an organism.');
            return false;
        }        
        
        if (!version) {
        	if (this.onError)
            	this.onError('Please specify a genome version.');
            return false;
        }

        var source = $('#edit_source').val();
        if (!source || source === 'Search') {
        	if (this.onError)
            	this.onError('Please specify a data source.');
            return false;
        }

       $.extend(this.genome, {
            metadata: {
                name: name,
                description: description,
                link: link,
                version: version,
                type: type,
                restricted: restricted,
                source: source,
                organism: organism
            },

            organism_id: this.organism_id
        });

        return true;
    },
});

function reset_load() {
    window.history.pushState({}, "Title", PAGE_NAME);

    // Reset the wizard and set description step
    initialize_wizard({
        admin: IS_ADMIN,
        metadata: current_genome.metadata,
        organism_id: current_genome.organism_id
    });

    $('#wizard-container').hide().fadeIn();
}

function initialize_wizard(opts) {
    current_genome = {};
    var root = $("#wizard-container");
    var wizard = new Wizard({ 
    	onCompleted: load_genome.bind(current_genome), 
    	data: current_genome, 
    	helpUrl: opts.helpUrl 
    });
    wizard.addStep(new GenomeDescriptionView({
        genome: current_genome,
        metadata: opts.metadata,
        organism_id: opts.organism_id,
        onError: wizard.error_help.bind(wizard)
    }));
    wizard.addStep(new DataView(current_genome, { supportedFileTypes: ['fasta', 'faa', 'fa'], onError: wizard.error_help.bind(wizard) }));
    //wizard.addStep(new OptionsView({genome: current_genome, admin: opts.admin}));
    wizard.addStep(new ConfirmationView(current_genome));

    // Render all the wizard sub views
    wizard.render();

    // Add the wizard to the document
    root.html(wizard.el);
    return wizard;
}
