function BatchDescriptionView(opts) {
    this.batch = opts.batch;
    this.metadata = opts.metadata;
    this.gid = opts.gid;
    this.onError = opts.onError;
    this.sources = undefined;
    this.title = "Describe your experiments";
    this.initialize();
}

$.extend(BatchDescriptionView.prototype, {
    initialize: function() {
        this.el = $($("#description-template").html());
        this.edit_genome = this.el.find("#edit_genome");

        if (this.metadata) {
            this.el.find('#edit_genome').val(this.metadata.genome);
        }
    },

    render: function() {
        var self = this;

        // jQuery Events
        this.edit_genome.unbind().change(function() {
            // Reset gid when item has changed
            self.gid = undefined;
        });

        // jQuery UI
        this.edit_genome.autocomplete({
            source:[],
            select: function(event, ui) {
                $(this).val(ui.item.label);
                self.gid = ui.item.value;
                return false; // Prevent the widget from inserting the value.
            },

            focus: function(event, ui) {
                //$("#edit_genome").val(ui.item.label);
                return false; // Prevent the widget from inserting the value.
            }
        });
    },

    is_valid: function() {
        var genome = this.el.find('#edit_genome').val();

        if (!genome || genome === 'Search' || !this.gid) {
        	if (this.onError)
            	this.onError('Please specify a genome.');
            return false;
        }

       $.extend(this.batch, {
            metadata: {
                genome: genome,
            },
            gid: this.gid
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
        this.edit_notebook = this.el.find("#edit_notebook");
        this.notebook_container = this.el.find("#notebook-container");
    },

    is_valid: function() {
        var notebook = this.edit_notebook.val();

        this.data.notebook = this.el.find("#notebook").is(":checked");
        this.data.notebook_type = this.el.find("[name=notebook] :checked").val();
        this.data.notebook_name = notebook;
        this.data.notebook_id = this.notebook_id;
        this.data.email = this.el.find("#email").is(":checked");

        if (this.data.notebook && this.data.notebook_type === "existing" && 
        		(!notebook || notebook === 'Search' || !this.notebook_id)) 
        {
        	if (this.onError)
            	this.onError('Please specify a notebook.');
            return false;
        }

        return true;
    },

    get_options: function() {
        return this.data;
    },
    
    render: function() {
        var self = this;

        // jQuery Events
        this.el.find("#notebook").unbind().change(this.toggleNotebook.bind(this));
        this.el.find("[name=notebook]").unbind().click(function() {
        	var option = $(this).val();
        	self.edit_notebook.prop("disabled", (option === 'new' ? true : false));
        });
        this.edit_notebook.unbind().change(function() {
            // Reset notebook_id when item has changed
            self.notebook_id = undefined;
        });

        // jQuery UI
        this.edit_notebook.autocomplete({
            source:[],
            select: function(event, ui) {
                $(this).val(ui.item.label);
                self.notebook_id = ui.item.value;
                return false; // Prevent the widget from inserting the value.
            },

            focus: function(event, ui) {
                return false; // Prevent the widget from inserting the value.
            }
        });
    },
    
    toggleNotebook: function() {
        this.notebook_enabled = this.el.find("#notebook").is(":checked");

        if (this.notebook_enabled) 
            this.notebook_container.slideDown();
        else 
            this.notebook_container.slideUp();
    }
});

function AdminOptionsView() {
    this.data = {};
    this.initialize();
}

$.extend(AdminOptionsView.prototype, {
    initialize: function() {
        this.el = $($("#admin-options-template").html());
        this.edit_user = this.el.find("#edit_user");
    },

    render: function() {
        // jQuery UI
        this.edit_user.unbind().autocomplete({
            source:[],
            focus: function() { return false; },
        });
    },

    is_valid: function() {
        //var ignore_cb = this.el.find('#ignore_missing_chrs');
        //this.ignore_missing_chrs = ignore_cb.is(':checked');
        if (this.edit_user.val()) {
            this.data.user = this.edit_user.val();
        }

        return true;
    },

    get_options: function() {
        return this.data;
    }
});

function OptionsView(opts) {
    this.batch = opts.batch;
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
        
        this.batch.options = {};
    },

    // Validate and add all options to the batch
    is_valid: function() {
        if (!this.general_view.is_valid())
            return false;

        var options = $.extend({}, this.general_view.get_options());
        if (this.admin) {
            if (!this.admin_view.is_valid())
                return false;
            $.extend(options, this.admin_view.get_options());
        }

        $.extend(this.batch.options, options);
        return true;
    },

    render: function() {
        var file_type = this.batch.data[0].file_type;
        if (!file_type) {
            if (this.onError)
            	this.onError("Please set the file type.");
            return;
        }

        this.layout_view.updateLayout(
            {"#analysis-options": this.analysis_view}
        );

        // Render the views added to the layout view
        this.layout_view.renderLayout();
    }
});

function render_template(template, container) {
    container.empty()
        .hide()
        .append(template)
        .show();//.slideDown();
}

function search_genomes (search_term) {
	coge.services.search_genomes(search_term, { fast: true })
		.done(function(result) { // success
			var transformed = result.genomes.map(function(obj) {
				var label = obj.info.replace(/&reg;/g, "\u00ae"); // (R) symbol
				return { label: label, value: obj.id };
			});
			$("#edit_genome")
				.autocomplete({source: transformed})
				.autocomplete("search");
		})
		.fail(function() { // error
			//TODO
		});
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

function search_notebooks (search_term) {
	coge.services.search_notebooks(search_term)
		.done(function(obj) {
			var notebooks = obj.notebooks;
			if (notebooks && notebooks.length > 0) {
				var items = [];
				notebooks.forEach(function(n) {
					var label = n.name + (n.description ? ': ' + n.description : '');
					//TODO //if (n.restricted) label = "\u00ae" + label; // (R)// Add (R) html symbol
	                items.push({
	                	label: label,
	                	value: n.id
	                });
	            });
	            $("#edit_notebook")
	            	.autocomplete({source: items})
	            	.autocomplete("search");
	        }
		});
}

function load(batch) {
	coge.progress.begin();
    newLoad = true;

	// Convert request into format for job service
	var request = {
		type: 'load_batch',
		requester: {
			page:      PAGE_NAME,
			user_name: USER_NAME
		},
		parameters: {
			genome_id:         batch.gid,
			metadata:          batch.metadata,
			email:             batch.options.email,
			notebook:          batch.options.notebook,
			notebook_name:     batch.options.notebook_name,
			notebook_id:       batch.options.notebook_id,
			source_data:       batch.data,
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
            window.history.pushState({}, "Title", "LoadBatch.pl" + "?wid=" + response.id); // Add workflow id to browser URL
            coge.progress.update(response.id, response.site_url);
	    })
	    .fail(function(jqXHR, textStatus, errorThrown) {
	    	coge.progress.failed("Couldn't talk to the server: " + textStatus + ': ' + errorThrown);
	    });
}

function reset_load() {
    window.history.pushState({}, "Title", PAGE_NAME);

    // Reset the wizard and set description step
    initialize_wizard({
        admin: IS_ADMIN,
        metadata: current_batch.metadata,
        gid: current_batch.gid
    });

    $('#wizard-container').hide().fadeIn();
}

function initialize_wizard(opts) {
    current_batch = {};
    var root = $("#wizard-container");
    var wizard = new Wizard({
    	onCompleted: load, 
    	data: current_batch, 
    	helpUrl: opts.helpUrl 
    });
    wizard.addStep(new BatchDescriptionView({
        batch: current_batch,
        metadata: opts.metadata,
        gid: opts.gid,
        onError: wizard.error_help.bind(wizard)
    }));
    wizard.addStep(new DataView(current_batch, { supportedFileTypes: ['gz'], onError: wizard.error_help.bind(wizard) }));
    wizard.addStep(new OptionsView({ batch: current_batch, admin: opts.admin, onError: wizard.error_help.bind(wizard) }));
    wizard.addStep(new ConfirmationView(current_batch));

    // Render all the wizard sub views
    wizard.render();

    // Add the wizard to the document
    root.html(wizard.el);
    return wizard;
}