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

	coge.services.add_organism({name: name, description: desc})
		.done(function(response) {
			if (response.id) {
				$('#create_new_organism_dialog').dialog('close');
				genomeDescriptionView.set_organism(response.id, response.name);
			}
		})
		.fail(function() {
			//TODO
		});
}

function load(genome) {
	coge.progress.begin();
    newLoad = true;

	// Convert request into format for job service
	var request = {
		type: 'load_genome',
		requester: {
			page: PAGE_NAME
		},
		parameters: {
			organism_id: genome.organism_id,
			metadata:    genome.metadata,
			//email:       genome.options.email,
			source_data: genome.data,
			load_id:     load_id
		}
	};
    
    coge.services.submit_job(request)
    	.done(function(response) {
    		if (!response) {
    			coge.progress.failed("Error: empty response from server");
    			return;
    		}
    		else if (!response.success || !response.id) {
    			coge.progress.failed("Error: failed to start workflow", response.error);
    			return;
    		}
    		
	        // Start status update
            window.history.pushState({}, "Title", PAGE_NAME + "?wid=" + response.id); // Add workflow id to browser URL
            coge.progress.update(response.id, response.site_url);
	    })
	    .fail(function(jqXHR, textStatus, errorThrown) {
	    	coge.progress.failed("Couldn't talk to the server: " + textStatus + ': ' + errorThrown);
	    });
}

function handle_action(action) {
	var w = window;
	if (w.parent)
		w = w.parent;
    if (action === "genome")
	    w.location.href = "GenomeInfo.pl?embed=" + EMBED + "&gid=" + genome_id;
    else if (action === "annotation")
	    w.location.href = "LoadAnnotation.pl?embed=" + EMBED + "&gid=" + genome_id;
    else //if (action === "new") 
        coge.progress.reset();
}

function search_organisms (search_term) {
	if (search_term && search_term.length > 2) {
		var edit_organism = $("#edit_organism");
		edit_organism.autocomplete("close");
		var spinner = $('#edit_organism_busy');
		spinner.css('visibility', 'visible');
		
		coge.services.search_organisms(search_term)
			.done(function(response) {
				var transformed = response.organisms.map(function(obj) {
					return { label: obj.name, value: obj.id };
				});
				edit_organism
					.autocomplete({source: transformed})
					.autocomplete("search");
				spinner.css('visibility', 'hidden');
			})
			.fail(function() {
				//TODO
			});
	}
}

function search_users (search_term) {
	coge.services.search_users(search_term)
		.done(function(result) {
			var transformed = result.users.map(function(obj) {
				return obj.user_name;
			});
			$("#edit_user")
				.autocomplete({source: transformed})
				.autocomplete("search");
		})
		.fail(function() {
			//TODO
		});
}

/*
 * NCBI Taxonomy Browser //TODO break this out into separate widget
 */

function search_ncbi_taxonomy(search_term) {
	if (!search_term || search_term.length < 3) {
		return;
	}

	$("#wait_ncbi").animate({opacity:1});
	$('#tax_empty').fadeOut();
	$("#edit_organism_info").slideUp();

	$.get( //FIXME use promises (mdb 10/25/16)
		ENTREZ_URL + "esearch.fcgi?db=taxonomy&term=" + search_term + "*",
		function(xml) {
			var ids = '';
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
								results.push({
								    id:      'taxon' + $(this).children('TaxId').text(),
								    text:    $(this).children('ScientificName').text(),
								    data: {
								        lineage: $(this).children('Lineage').text()
								    }
                                });
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

function build_taxonomy_tree(items) {
	if (items) {
		var tree = $("#tax_tree");

		if (!items || items.length == 0) {
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
		    tree.hide().jstree("destroy").empty();
		    tree.jstree({
                    core: {
                        check_callback : true,
                        data: items
                    }
                })
                .bind("select_node.jstree",
                    function (event, data) {
                        var id = data.selected[0];
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

function taxonomy_get_node(id) { // FIXME: merge common stuff with search_ncbi_taxonomy
	var tree = $("#tax_tree");
    var node = tree.jstree().get_node(id);

	// If the node is open then it has children so just close it and return
	var is_open = tree.jstree("is_open", node);
	if (is_open) {
		tree.jstree("close_node", node);
		return;
	}

	// If not a leaf node then we already retrieved the children so just open it and return.
	var is_leaf = tree.jstree("is_leaf", node);
	if (!is_leaf) {
		tree.jstree("open_node", node);
		return;
	}

	// Retrieve the children for this taxon.
	$("#wait_ncbi").animate({opacity:1});

	$.get(
		ENTREZ_URL + "esearch.fcgi?db=taxonomy&term=" + node.text + "&field=nxlv",
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
								results.push({
								    id:   'taxon' + $(this).children('TaxId').text(),
								    text: $(this).children('ScientificName').text(),
								    data: {
								        lineage: $(this).children('Lineage').text()
								    }
								});
							}
						);

						if (results.length) {
							results.sort(sort_name).forEach(
								function(newNode) {
									tree.jstree().create_node(node, newNode, "last");
								}
							);
							tree.jstree().open_node(node);
						}

						$("#wait_ncbi").animate({opacity:0});
					}
				);
			}
			else {
				$("#edit_organism_name").val(node.text);
				$("#edit_organism_desc").val(node.data.lineage);
				activate_on_input(['edit_organism_name', 'edit_organism_desc'], 'create_organism_button');
				$("#edit_organism_info").slideDown();
				$("#wait_ncbi").animate({opacity:0});
			}
		}
	);
}

function sort_name (a,b) {
	var nameA=a.text.toLowerCase(), nameB=b.text.toLowerCase()
	if (nameA < nameB) //sort string ascending
		return -1
	if (nameA > nameB)
		return 1
	return 0 //default return value (no sorting)
}

function activate_on_input(elements, target) {
	if (elements) {
		elements.forEach(function(e) {
			if ($('#'+e).val())
				$('#'+target).removeClass('ui-state-disabled');
		});
	}
}

/*
 * Wizard
 */

function GenomeDescriptionView(opts) {
    this.genome = opts.genome;
    this.metadata = opts.metadata;
    this.organism_id = opts.organism_id;
    this.onError = opts.onError;
    this.sources = undefined;
    this.title = "Describe Genome";
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
            this.select_type.val(this.metadata.type_id);
            this.edit_source.val(this.metadata.source_name);

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
    
    set_organism: function(id, name) {
    	this.edit_organism.val(name);
    	this.organism_id = id;
    },

    render: function() {
        var self = this;
        var edit_organism = this.edit_organism;

// mdb removed 10/17/16 -- broken autocomplete when updated to jQuery 3.1.1
//        edit_organism.unbind().change(function() {
//            // Reset organism_id when item has changed
//            self.organism_id = undefined;
//        });

        edit_organism.autocomplete({
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
        
        edit_organism.keyup(function() {
        	coge.utils.wait_to_search(search_organisms, self.edit_organism.get(0));
        });
        edit_organism.click(function() {
        	$(this).autocomplete('search');
        });

        this.edit_source.autocomplete({source: this.sources});
    },

    is_valid: function() {
        var name        = this.el.find('#edit_name').val();
        var description = this.el.find('#edit_description').val();
        var link        = this.el.find('#edit_link').val();
        var version     = this.el.find('#edit_version').val();
        var type_id     = this.el.find('#select_type').val();
        var restricted  = this.el.find('#restricted').is(':checked');
        var organism    = this.el.find('#edit_organism').val();

        if (!organism || !this.organism_id) {
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
                name:             coge.utils.removeSpecialChars(name),
                description:      coge.utils.removeSpecialChars(description),
                link:             link,
                version:          coge.utils.removeSpecialChars(version),
                sequence_type_id: type_id,
                restricted:       restricted,
                source_name:      coge.utils.removeSpecialChars(source),
                organism:         organism
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

var genomeDescriptionView; // needed by taxonomy browser
function initialize_wizard(opts) {
    current_genome = {};
    var root = $("#wizard-container");
    var wizard = new Wizard({ 
    	onCompleted: load.bind(current_genome), 
    	data: current_genome, 
    	helpUrl: opts.helpUrl 
    });
    genomeDescriptionView = new GenomeDescriptionView({
        genome: current_genome,
        metadata: opts.metadata,
        organism_id: opts.organism_id,
        onError: wizard.error_help.bind(wizard)
    });
    wizard.addStep(genomeDescriptionView);
    wizard.addStep(new DataView(current_genome, { supportedFileTypes: ['fasta', 'faa', 'fa'], onError: wizard.error_help.bind(wizard) }));
    //wizard.addStep(new OptionsView({genome: current_genome, admin: opts.admin}));
    wizard.addStep(new ConfirmationView(current_genome));

    // Render all the wizard sub views
    wizard.render();

    // Add the wizard to the document
    root.html(wizard.el);
    return wizard;
}
