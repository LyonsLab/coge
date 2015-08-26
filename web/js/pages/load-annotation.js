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

function search_genomes (search_term) {
	coge.services.search_genomes(search_term)
		.done(function(result) {
			var transformed = result.genomes.map(function(obj) {
				return { label: obj.info, value: obj.id };
			});
			$("#edit_genome")
				.autocomplete({source: transformed})
				.autocomplete("search");
		})
		.fail(function() {
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

function load(annotation) {
	coge.progress.begin();
    newLoad = true;

	// Convert request into format for job service
	var request = {
		type: 'load_annotation',
		requester: {
			page:      PAGE_NAME,
			user_name: USER_NAME
		},
		parameters: {
			genome_id:         annotation.gid,
			metadata:          annotation.metadata,
			//email:             annotation.options.email,
			source_data:       annotation.data,
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
            window.history.pushState({}, "Title", "LoadAnnotation.pl" + "?wid=" + response.id); // Add workflow id to browser URL
            coge.progress.update(response.id, response.site_url);
	    })
	    .fail(function(jqXHR, textStatus, errorThrown) {
	    	coge.progress.failed("Couldn't talk to the server: " + textStatus + ': ' + errorThrown);
	    });
}

function AnnotationDescriptionView(opts) {
    this.annotation = opts.annotation;
    this.metadata = opts.metadata;
    this.gid = opts.gid;
    this.onError = opts.onError;
    this.sources = undefined;
    this.title = "Describe your annotation";
    this.initialize();
}

$.extend(AnnotationDescriptionView.prototype, {
    initialize: function() {
        this.el = $($("#description-template").html());
        this.edit_source = this.el.find("#edit_source");
        this.edit_genome = this.el.find("#edit_genome");
        this.get_sources();

        if (this.metadata) {
            this.el.find('#edit_name').val(this.metadata.name);
            this.el.find('#edit_description').val(this.metadata.description);
            this.el.find('#edit_version').val(this.metadata.version);
            this.edit_source.val(this.metadata.source);
            this.el.find('#edit_genome').val(this.metadata.genome);
        }
    },

    get_sources: function() {
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

        this.edit_source.autocomplete({source: this.sources});
    },

    is_valid: function() {
        var name = this.el.find('#edit_name').val();
        var description = this.el.find('#edit_description').val();
        var link = this.el.find('#edit_link').val();
        var version = this.el.find('#edit_version').val();
        var genome = this.el.find('#edit_genome').val();

        if (!genome || genome === 'Search' || !this.gid) {
        	if (this.onError)
            	this.onError('Please specify a genome.');
            return false;
        }
        
        if (!version) {
        	if (this.onError)
            	this.onError('Please specify an annotation version.');
            return false;
        }

        var source = $('#edit_source').val();
        if (!source || source === 'Search') {
        	if (this.onError)
            	this.onError('Please specify a data source.');
            return false;
        }

       $.extend(this.annotation, {
            metadata: {
                name: name,
                description: description,
                link: link,
                version: version,
                source: source,
                genome: genome,
            },

            gid: this.gid
        });

        return true;
    },
});

function reset_load() {
    window.history.pushState({}, "Title", PAGE_NAME);

    // Reset the wizard and set description step
    initialize_wizard({
        admin: IS_ADMIN,
        metadata: current_annotation.metadata,
        gid: current_annotation.gid
    });

    $('#wizard-container').hide().fadeIn();
}

function initialize_wizard(opts) {
    current_annotation = {};
    var root = $("#wizard-container");
    var wizard = new Wizard({ 
    	onCompleted: load, 
    	data: current_annotation, 
    	helpUrl: opts.helpUrl 
    });
    wizard.addStep(new AnnotationDescriptionView({
        annotation: current_annotation,
        metadata: opts.metadata,
        gid: opts.gid,
        onError: wizard.error_help.bind(wizard)
    }));
    wizard.addStep(new DataView(current_annotation, { supportedFileTypes: ['gff', 'gtf'], onError: wizard.error_help.bind(wizard) }));
    //wizard.addStep(new OptionsView({annotation: current_annotation, admin: opts.admin, onError: wizard.error_help }));
    wizard.addStep(new ConfirmationView(current_annotation));

    // Render all the wizard sub views
    wizard.render();

    // Add the wizard to the document
    root.html(wizard.el);
    return wizard;
}