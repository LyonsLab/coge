/*global document,$,alert,window,JSON */

var concat = Array.prototype.concat;

// Global experiment data
var current_experiment = {};

// Supported file types
var POLY_FILES = [
    "vcf"
];

var ALIGN_FILES = [
    "bam"
];

var SEQ_FILES = [
    "fastq", "fq"
];

var QUANT_FILES = [
    "csv", "tsv", "bed", "gff", "gtf"
];

var SUPPORTED_FILE_TYPES = concat.call(QUANT_FILES, ALIGN_FILES, SEQ_FILES, POLY_FILES);

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

function render_template(template, container) {
    container.empty()
    .hide()
    .append(template)
    .show();//.slideDown();
}

function ExperimentDescriptionView(opts) {
    this.experiment = opts.experiment;
    this.metadata = opts.metadata;
    this.gid = opts.gid;
    this.onError = opts.onError;
    this.sources = undefined;
    this.title = "Describe your experiment";
    this.initialize();
}

$.extend(ExperimentDescriptionView.prototype, {
    initialize: function() {
        this.el = $($("#description-template").html());
        this.edit_source = this.el.find("#edit_source");
        this.edit_genome = this.el.find("#edit_genome");
        this.get_sources();

        if (this.metadata) {
            this.el.find('#edit_name').val(this.metadata.name);
            this.el.find('#edit_description').val(this.metadata.description);
            this.el.find('#edit_version').val(this.metadata.version);
            this.edit_source.val(this.metadata.source_name);

            if (!this.metadata.restricted)
                this.el.find('#restricted').removeAttr('checked');

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
        
        var edit_genome = this.edit_genome;

        edit_genome.unbind().change(function() {
            // Reset gid when item has changed
            self.gid = undefined;
        });

        edit_genome.autocomplete({
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
        
        edit_genome.keyup(function() {
        	coge.utils.wait_to_search(search_genomes, self.edit_genome.get(0));
        });
        edit_genome.click(function() {
        	$(this).autocomplete('search');
        });

        this.edit_source.autocomplete({source: this.sources});
    },

    is_valid: function() {
        var name = this.el.find('#edit_name').val();
        var description = this.el.find('#edit_description').val();
        var version = this.el.find('#edit_version').val();
        var restricted = this.el.find('#restricted').is(':checked');
        var genome = this.el.find('#edit_genome').val();

        if (!name) {
            if (this.onError)
            	this.onError('Please specify an experiment name.');
            return false;
        }
        
        if (!version) {
        	if (this.onError)
            	this.onError('Please specify an experiment version.');
            return false;
        }

        var source = $('#edit_source').val();
        if (!source || source === 'Search') {
        	if (this.onError)
            	this.onError('Please specify a data source.');
            return false;
        }

        if (!genome || genome === 'Search' || !this.gid) {
        	if (this.onError)
            	this.onError('Please specify a genome.');
            return false;
        }

       $.extend(this.experiment, {
            metadata: {
                name: name,
                description: description,
                version: version,
                restricted: restricted,
                source_name: source,
                genome: genome,
            },

            gid: this.gid
        });

        return true;
    },
});

function FindSNPView() {
    this.data = {};
    this.initialize();
};

$.extend(FindSNPView.prototype, {
    initialize: function() {
        this.el = $($("#snp-template").html());
        this.snp_container = this.el.find("#snp-container");

        this.snp_templates = {
            coge: $($("#coge-snp-template").html()),
            samtools: $($("#samtools-snp-template").html()),
            platypus: $($("#platypus-snp-template").html()),
            gatk: $($("#gatk-snp-template").html())
        };
    },

    render: function() {
        var self = this;

        var method = this.el.find("#snp-method");
        this.el.find("#snps").unbind().change(this.update_snp.bind(this));

        // Events to rebind when the view is added to the dom
        method.unbind().change(function() {
            var selected = $(this).val();
            render_template(self.snp_templates[selected], self.snp_container);
        });

        if (this.data.snp_params) {
            method.val(this.data.snp_params.method);
        }
    },

    // Callback to display the selected snp pipeline
    update_snp: function (ev) {
        var enabled = $(ev.target).is(":checked"),
            method = this.el.find("#snp-method");

        var el = $(document.getElementById(method.val()));

        if (enabled) {
            this.data.snp_params = $.extend({}, this.data.snp_params, { method: method.val() });
            el.show();
            method.removeAttr("disabled");
            this.snp_container.slideDown();
            var selected = $("#snp-method").val();
            render_template(this.snp_templates[selected], this.snp_container);
        } else {
            this.data.snp_params = undefined;
            method.attr("disabled", 1);
            this.snp_container.slideUp();
        }
    },

    is_valid: function() {
        // SNP pipeline
        var enabled = this.el.find("#snps").is(":checked");
        var method = this.el.find("#snp-method").val();

        if (enabled) {
        	//TODO this can be automated
            if (method === "coge") {
                this.data.snp_params = {
                    method: method,
                    'min-read-depth':   this.el.find("#min-read-depth").val(),
                    'min-base-quality': this.el.find("#min-base-quality").val(),
                    'min-allele-count': this.el.find("#min-allele-count").val(),
                    'min-allele-freq':  this.el.find("#min-allele-freq").val(),
                    scale: this.el.find("#scale").val()
                };
            } else if (method === "samtools") {
                this.data.snp_params = {
                    method: method,
                    'min-read-depth': this.el.find("#min-read-depth").val(),
                    'max-read-depth': this.el.find("#max-read-depth").val(),
                };
            } else if (method === "platypus") {
                this.data.snp_params = {
                    method: method
                };
            } else if (method === "gatk") {
                this.data.snp_params = {
                    method: method
                };
            }
        }
        return true;
    },

    get_options: function() {
        return this.data;
    },
});

function AlignmentView() {
    this.data = {};
    this.initialize();
}

$.extend(AlignmentView.prototype, {
    initialize:function() {
        this.el = $($("#align-template").html());
        this.align_container = this.el.find("#align-container");
        this.align_templates = {
            gsnap: $($("#gsnap-template").html()),
            tophat: $($("#tophat-template").html())
        };
    },

    render: function() {
        // jQuery events
        this.el.find("[name=aligner]").unbind().click(this.update_aligner.bind(this));
        this.update_aligner();
    },

    // Callback to display the selected aligner
    update_aligner: function() {
        var selected = this.el.find("#alignment :checked").val();
        render_template(this.align_templates[selected], this.align_container);
    },

    is_valid: function() {
        var aligner = this.el.find("#alignment :checked").val();

        // Pick the aligner and set the options
        if (aligner === "gsnap") {
            this.data = {
                alignment_params: { //TODO is there a way to automate this parameter passing?
                    tool: "gsnap",
                    '-n': this.el.find("[id='-n']").val(),
                    '-Q': this.el.find("[id='-Q']").is(":checked"),
                    '--gap-mode': this.el.find("[id='--gap-mode']").val(),
                    '--nofails': this.el.find("[id='--nofails']").is(":checked"),
                    read_type: this.el.find("#read_type :checked").val()
                },
                trimming_params: {
                    '-q': this.el.find("[id='-q']").val(),
                    '-m': this.el.find("[id='-m']").val(),
                    '--quality-base': this.el.find("[id='--quality-base']").val()
                }
            };
        } else {
            this.data = {
                alignment_params: {
                    tool: "tophat",
                    '-g': this.el.find("[id='-g']").val(),
                    read_type: this.el.find("#read_type :checked").val()
                }
            }
        }

        return true;
    },

    get_options: function() {
        return this.data;
    }
});

function PolymorphismView() {
    this.initialize();
}

$.extend(PolymorphismView.prototype, {
    initialize: function() {
        this.el = $($("#poly-template").html());
    },

    is_valid: function() {
        return true;
    },

    get_options: function() {
        return {};
    },
});

function AlignmentOptionView() {
    this.initialize();
}

$.extend(AlignmentOptionView.prototype, {
    initialize: function() {
        this.snp_view = new FindSNPView();
        this.expression_view = new ExpressionView();

        this.layout_view = new LayoutView({
            template: "#align-option-template",

            layout: {
                "#expression-view": this.expression_view,
                "#snp-view": this.snp_view
            }
        });

        this.el = this.layout_view.el;
    },

    render: function() {
        this.layout_view.renderLayout();
    },

    is_valid: function() {
        return this.snp_view.is_valid();
    },

    get_options: function() {
        return $.extend(this.snp_view.get_options(),
                        this.expression_view.get_options());
    },
});

function QuantativeView(){
    this.initialize();
    this.data = {};
}

$.extend(QuantativeView.prototype, {
    initialize: function() {
        this.el = $($("#quant-template").html());
        this.container = this.el.find("#normalize_method");
    },

    render: function() {
        this.el.find("#normalize").unbind().change(this.toggleAnalysis.bind(this));
    },

    toggleAnalysis: function() {
        this.enabled = this.el.find("#normalize").is(":checked");

        if (this.enabled) {
            this.container.slideDown();
        } else {
            this.container.slideUp();
        }
    },

    is_valid: function() {
        this.data.normalize = this.el.find("#normalize").is(":checked");
        return true;
    },

    get_options: function() {
        if (this.enabled)
            this.data.normalize_method = this.el.find("#percentage").is(":checked") ? 'percentage' : this.el.find("#log10").is(":checked") ? 'log10' : this.el.find("#loge").is(":checked") ? 'loge' : null;
        return this.data;
    },
});

function ExpressionView() {
    this.initialize();
    this.data = {};
}

$.extend(ExpressionView.prototype, {
    initialize: function() {
        this.el = $($("#expression-template").html());
        this.enabled = false;
        this.container = this.el.find("#expression-container");
    },

    render: function() {
        this.el.find("#expression").unbind().change(this.toggleAnalysis.bind(this));
    },

    toggleAnalysis: function() {
        this.enabled = this.el.find("#expression").is(":checked");

        if (this.enabled) {
            this.container.slideDown();
        } else {
            this.container.slideUp();
        }
    },

    is_valid: function() {
        return true;
    },

    get_options: function() {
        if (this.enabled) {
            this.data.expression_params = {
                '-Q': this.el.find("[id='-Q']").val()
            };
        }

        return this.data;
    },

})

function FastqView() {
    this.initialize();
}

$.extend(FastqView.prototype, {
    initialize: function() {
        this.expression_view = new ExpressionView();
        this.snp_view = new FindSNPView();
        this.align_view = new AlignmentView();

        this.layout_view = new LayoutView({
            template: "#fastq-template",

            layout: {
                "#expression-view": this.expression_view,
                "#snp-view": this.snp_view,
                "#align-view": this.align_view
            }
        });

        // pass through to the layout
        this.el = this.layout_view.el;
    },

    render: function() {
        this.layout_view.renderLayout();
    },

    is_valid: function() {
        if (!this.snp_view.is_valid())
            return false;

        if (!this.align_view.is_valid())
            return false;

        if (!this.expression_view.is_valid())
            return false;

        return true;
    },

    get_options: function() {
        return $.extend(this.expression_view.get_options(),
                        this.snp_view.get_options(),
                        this.align_view.get_options());
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
        if (!this.analysis_view.is_valid() || !this.general_view.is_valid())
            return false;

        var options = $.extend({}, this.general_view.get_options(), this.analysis_view.get_options());
        if (this.admin) {
            if (!this.admin_view.is_valid())
                return false;
            $.extend(options, this.admin_view.get_options());
        }

        $.extend(this.experiment.options, options);
        return true;
    },

    render: function() {
        var file_type = this.experiment.data[0].file_type;
        if (!file_type) {
            if (this.onError)
            	this.onError("Please set the file type.");
            return;
        }

        //FIXME: An aggregate view should add analysis options for multiple file types
        if ($.inArray(file_type, POLY_FILES) > -1)
            this.analysis_view = new PolymorphismView();
        else if ($.inArray(file_type, SEQ_FILES) > -1)
            this.analysis_view = new FastqView();
        else if ($.inArray(file_type, QUANT_FILES) > -1)
            this.analysis_view = new QuantativeView();
        else if ($.inArray(file_type, ALIGN_FILES) > -1)
            this.analysis_view = new AlignmentOptionView();

        this.layout_view.updateLayout(
            {"#analysis-options": this.analysis_view}
        );

        // Render the views added to the layout view
        this.layout_view.renderLayout();
    },
});

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

function reset_load() {
    window.history.pushState({}, "Title", PAGE_NAME);

    // Reset the wizard and set description step
    initialize_wizard({
        admin: IS_ADMIN,
        metadata: current_experiment.metadata,
        gid: current_experiment.gid
    });

    $('#wizard-container').hide().fadeIn();
}

function initialize_wizard(opts) {
    current_experiment = {};
    var root = $("#wizard-container");
    var wizard = new Wizard({ 
    	onCompleted: load, 
    	data: current_experiment, 
    	helpUrl: opts.helpUrl 
    });
    wizard.addStep(new ExperimentDescriptionView({
        experiment: current_experiment,
        metadata: opts.metadata,
        gid: opts.gid,
        onError: wizard.error_help.bind(wizard)
    }));
    wizard.addStep(new DataView(current_experiment, { supportedFileTypes: SUPPORTED_FILE_TYPES, onError: wizard.error_help.bind(wizard) }));
    wizard.addStep(new OptionsView({experiment: current_experiment, admin: opts.admin, onError: wizard.error_help }));
    wizard.addStep(new ConfirmationView(current_experiment));

    // Render all the wizard sub views
    wizard.render();

    // Add the wizard to the document
    root.html(wizard.el);
    return wizard;
}
