/*global document,$,alert,window,JSON */

var concat = Array.prototype.concat;

// Global experiment data
var current_experiment = {};
var focus = undefined;

//AKB - Need to work on this "search-genomes"
function search_genomes (search_term) {
	coge.services.search_genomes(search_term, { fast: true })
		.done(function(result) { // success
			if (result && result.genomes) {
				var transformed = result.genomes.map(function(obj) {
					var label = obj.info.replace(/&reg;/g, "\u00ae"); // (R) symbol
					return { label: label, value: obj.id };
				});
				if (focus == 'x') {
					$('#edit_xgenome')
						.autocomplete({source: transformed})
						.autocomplete("search");
				} else if (focus == 'y') {
					$('#edit_ygenome')
						.autocomplete({source: transformed})
						.autocomplete("search");
				} else if (focus == 'z') {
					$('#edit_zgenome')
						.autocomplete({source: transformed})
						.autocomplete("search");
				}				
			}
		})
		.fail(function() { // error
			//TODO
		});
}

/* AKB removed - probably wont need
function search_users (search_term) {
	coge.services.search_users(search_term)
		.done(function(result) {
			if (result && result.users) {
				var transformed = result.users.map(function(obj) {
					return obj.user_name;
				});
				$("#edit_user")
					.autocomplete({source: transformed})
					.autocomplete("search");
			}
		})
		.fail(function() {
			//TODO
		});
}
*/

/* AKB removed - probably won't need
function search_notebooks (search_term) {
	coge.services.search_notebooks(search_term)
		.done(function(data) { 
			if (data.notebooks) {
				var items = [];
				data.notebooks.forEach(function(notebook) {
					items.push({label:notebook.name,value:notebook.id});
			    });
			    $("#edit_notebook")
			    	.autocomplete({source: items})
			    	.autocomplete("search");
			}
		})
		.fail(function() {
			//TODO
		});
}
*/

function render_template(template, container) {
    container.empty()
    .hide()
    .append(template)
    .show();//.slideDown();
}

function ExperimentDescriptionView(opts) {
    this.experiment = opts.experiment;
    this.metadata = opts.metadata;
    //this.gid = opts.gid; AKB - replaced with x_gid, y_gid, z_gid
    this.x_gid = opts.x_gid;
    this.y_gid = opts.y_gid;
    this.z_gid = opts.z_gid;
    this.onError = opts.onError;
    this.sources = undefined;
    this.title = "Select Genomes";
    this.initialize();
}

$.extend(ExperimentDescriptionView.prototype, {
    initialize: function() {
        this.el = $($("#genomes-template").html());
        this.edit_xgenome = this.el.find("#edit_xgenome");
	this.edit_ygenome = this.el.find("#edit_ygenome");
	this.edit_zgenome = this.el.find("#edit_zgenome")
        if (this.metadata) {
            //this.el.find('#edit_name').val(this.metadata.name);
            //this.el.find('#edit_description').val(this.metadata.description);
            //this.el.find('#edit_version').val(this.metadata.version);
            //this.edit_source.val(this.metadata.source_name);

            //if (!this.metadata.restricted)
            //    this.el.find('#restricted').removeAttr('checked');

            //this.el.find('#edit_genome').val(this.metadata.genome); AKB REMOVED 9/15
        }
    },

    render: function() {
        var self = this;
        
        var edit_xgenome = this.edit_xgenome;
	var edit_ygenome = this.edit_ygenome;
	var edit_zgenome = this.edit_zgenome;

        edit_xgenome.unbind().change(function() {
            // Reset gid when item has changed
            self.x_gid = undefined;
        });
        edit_ygenome.unbind().change(function() {
            // Reset gid when item has changed
            self.y_gid = undefined;
        });
        edit_zgenome.unbind().change(function() {
            // Reset gid when item has changed
            self.z_gid = undefined;
        });

        edit_xgenome.autocomplete({
            source:[],
            select: function(event, ui) {
                $(this).val(ui.item.label);
                self.x_gid = ui.item.value;
                return false; // Prevent the widget from inserting the value.
            },

            focus: function(event, ui) {
                //$("#edit_genome").val(ui.item.label);
                return false; // Prevent the widget from inserting the value.
            }
        });
        edit_ygenome.autocomplete({
            source:[],
            select: function(event, ui) {
                $(this).val(ui.item.label);
                self.y_gid = ui.item.value;
                return false; // Prevent the widget from inserting the value.
            },

            focus: function(event, ui) {
                //$("#edit_genome").val(ui.item.label);
                return false; // Prevent the widget from inserting the value.
            }
        });
        edit_zgenome.autocomplete({
            source:[],
            select: function(event, ui) {
                $(this).val(ui.item.label);
                self.z_gid = ui.item.value;
                return false; // Prevent the widget from inserting the value.
            },

            focus: function(event, ui) {
                //$("#edit_genome").val(ui.item.label);
                return false; // Prevent the widget from inserting the value.
            }
        });
       
 
        edit_xgenome.keyup(function() {
		focus = 'x';
        	coge.utils.wait_to_search(search_genomes, self.edit_xgenome.get(0));
        });
        edit_xgenome.click(function() {
		focus = 'x';
        	$(this).autocomplete('search');
        });
        edit_ygenome.keyup(function() {
        	focus = 'y';
		coge.utils.wait_to_search(search_genomes, self.edit_ygenome.get(0));
        });
        edit_ygenome.click(function() {
        	focus = 'y';
		$(this).autocomplete('search');
        });
        edit_zgenome.keyup(function() {
        	focus = 'z';
		coge.utils.wait_to_search(search_genomes, self.edit_zgenome.get(0));
        });
        edit_zgenome.click(function() {
        	focus = 'z';
		$(this).autocomplete('search');
        });

        //this.edit_source.autocomplete({source: this.sources});
    },

    is_valid: function() {
        var xgenome = this.el.find('#edit_xGenome').val();
	var ygenome = this.el.find('#edit_yGenome').val();
	var zgenome = this.el.find('#edit_zGenome').val();
        
	if (!xgenome || xgenome === 'Search' || !this.x_gid) {
        	if (this.onError)
            	this.onError('Please specify an X-axis genome.');
            return false;
        }

	if (!ygenome || ygenome === 'Search' || !this.y_gid) {
        	if (this.onError)
            	this.onError('Please specify a Y-axis genome.');
            return false;
        }

	if (!zgenome || zgenome === 'Search' || !this.z_gid) {
        	if (this.onError)
            	this.onError('Please specify a Z-axis genome.');
            return false;
        }

       $.extend(this.experiment, {
            metadata: {
		xgenome: xgenome,
		ygenome: ygenome,
		zgenome: zgenome
                //name: name,
                //description: description,
                //version: version,
                //restricted: restricted,
                //source_name: source,
                //genome: genome,
            },
	    
	    //gid: this.gid,
            x_gid: this.x_gid,
	    y_gid: this.y_gid,
	    z_gid: this.z_gid
        });

        return true;
    },
});

function GeneralOptionsView(opts) {
    this.data = {};
    this.initialize();
    this.onError = opts.onError;
}

$.extend(GeneralOptionsView.prototype, {
    initialize: function() {
        this.el = $($("#general-options-template").html());
        //this.edit_notebook = this.el.find("#edit_notebook");
        //this.notebook_container = this.el.find("#notebook-container");
    },

    is_valid: function() {
        return true;
    },

    get_options: function() {
        return this.data;
    }
});

function AdminOptionsView() {
    this.data = {};
    this.initialize();
}

$.extend(AdminOptionsView.prototype, {
    initialize: function() {
        this.el = $($("#admin-options-template").html());
    },

    is_valid: function() {
        return true;
    },

    get_options: function() {
        return this.data;
    },
});

function OptionsView(opts) {
    this.experiment = opts.experiment;
    this.admin = opts.admin;
    this.onError = opts.onError;
    this.title = "Options";
    this.initialize();
}

$.extend(OptionsView.prototype, {
    initialize: function() {
        this.admin_view = new AdminOptionsView();
        this.general_view = new GeneralOptionsView({onError: this.onError});
        this.layout_view = new LayoutView({
            template: "#options-layout-template",
            layout: {
                "#general-options": this.general_view
            }
        });

        if (this.admin)
            this.layout_view.updateLayout({"#admin-options": this.admin_view});

        this.el = this.layout_view.el;
        
        this.experiment.options = {};
    },

    // Validate and add all options to the experiment
    is_valid: function() {
        if (!this.advanced_view.is_valid() || !this.general_view.is_valid())
            return false;

        var options = $.extend({}, this.general_view.get_options(), this.advanced_view.get_options());
        if (this.admin) {
            if (!this.admin_view.is_valid())
                return false;
            $.extend(options, this.admin_view.get_options());
        }

        $.extend(this.experiment.options, options);
        return true;
    },

    render: function() {
	this.advanced_view = new AdvancedOptionView();
        this.layout_view.updateLayout(
            {"#advanced-options": this.advanced_view}
        );

        // Render the views added to the layout view
        this.layout_view.renderLayout();
    },
});

function AdvancedOptionView() {
    this.initialize();
}

$.extend(AdvancedOptionView.prototype, {
    initialize: function() {
	this.el = $($("#advanced-options-template").html());
    },

    is_valid: function() {
	return true;
    },

    get_options: function() {
	return {};
    }
});

/* TODO: AKB-Need to make launch(experiment), but use this as a template.
function load(experiment) {
	coge.progress.begin();
    newLoad = true;

	// Convert request into format for job service
	var request = {
		type: 'load_experiment',
		requester: {
			page:      PAGE_NAME,
			user_name: USER_NAME
		},
		parameters: {
			genome_id:         experiment.gid,
			metadata:          experiment.metadata,
			alignment_params:  experiment.options.alignment_params,
			trimming_params:   experiment.options.trimming_params,
			expression_params: experiment.options.expression_params,
			snp_params:        experiment.options.snp_params,
			normalize:         experiment.options.normalize,
			normalize_method:  experiment.options.normalize_method,
			email:             experiment.options.email,
			notebook:          experiment.options.notebook,
			notebook_name:     experiment.options.notebook_name,
			notebook_id:       experiment.options.notebook_id,
			source_data:       experiment.data,
			load_id:           load_id
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
            window.history.pushState({}, "Title", "LoadExperiment.pl" + "?wid=" + response.id); // Add workflow id to browser URL
            coge.progress.update(response.id, response.site_url);
	    })
	    .fail(function(jqXHR, textStatus, errorThrown) {
	    	coge.progress.failed("Couldn't talk to the server: " + textStatus + ': ' + errorThrown);
	    });
}
*/

function launch(experiment) {};

function reset_launch() { //AKB - Renamed from reset_load()
    window.history.pushState({}, "Title", PAGE_NAME);

    // Reset the wizard and set description step
    initialize_wizard({
        admin: IS_ADMIN,
        metadata: current_experiment.metadata,
        //gid: current_experiment.gid, //AKB - replaced with x_gid, y_gid, z_gid
	x_gid: current_experiment.x_gid,
	y_gid: current_experiment.y_gid,
	z_gid: current_experiment.z_gid
    });

    $('#wizard-container').hide().fadeIn();
}

function initialize_wizard(opts) {
    current_experiment = {};
    var root = $("#wizard-container");
    var wizard = new Wizard({ 
    	onCompleted: launch, //TODO - AKB, changed from 'load' 
    	data: current_experiment, 
    	helpUrl: opts.helpUrl 
    });
    wizard.addStep(new ExperimentDescriptionView({
        experiment: current_experiment,
        metadata: opts.metadata,
        //gid: opts.gid, AKB - replaced with x_gid, y_gid, z_gid
        x_gid: opts.x_gid,
        y_gid: opts.y_gid,
        z_gid: opts.z_gid,
        onError: wizard.error_help.bind(wizard)
    }));
    //wizard.addStep(new DataView(current_experiment, { supportedFileTypes: SUPPORTED_FILE_TYPES, onError: wizard.error_help.bind(wizard) })); //AKB - removed
    wizard.addStep(new OptionsView({experiment: current_experiment, admin: opts.admin, onError: wizard.error_help }));
    wizard.addStep(new ConfirmationView(current_experiment));

    // Render all the wizard sub views
    wizard.render();

    // Add the wizard to the document
    root.html(wizard.el);
    return wizard;
}
