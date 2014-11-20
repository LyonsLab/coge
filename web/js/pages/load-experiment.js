/*global document,$,alert,window,JSON */

var concat = Array.prototype.concat;

// Global experiment data
var current_experiment = {};

// Support file types
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

var SUPPORTED_FILES = concat.call(QUANT_FILES, ALIGN_FILES, SEQ_FILES, POLY_FILES);
var FILE_TYPE_PATTERNS = new RegExp("(:?" + SUPPORTED_FILES.join("|") + ")$");

// Template support
var snp_templates  = {};
var align_templates = {};


function autodetect_file_type(file) {
    var stripped_file = file.replace(/.gz$/, '');

    if (FILE_TYPE_PATTERNS.test(stripped_file)) {
        return stripped_file.match(FILE_TYPE_PATTERNS)[0];
    }
}

function error_help(s) {
    $('#error_help_text')
        .html(s)
        .show()
        .delay(10*1000)
        .fadeOut(1500);
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

function LayoutView(options) {
    this.template = $(options.template);
    this.layout = options.layout;
    this.initialize();
}

$.extend(LayoutView.prototype, {
    initialize: function() {
        this.el = $(this.template.html());
    },

    renderLayout: function() {
        var elementId, section, view;

        for (elementId in this.layout) {
            if (this.layout.hasOwnProperty(elementId)) {
                section = this.el.find(elementId);
                view = this.layout[elementId];

                section.empty();
                section.html(view.el);

                if (view.render) {
                    view.render();
                }
            }
        }
    },

    updateLayout: function(layout) {
        this.layout = $.extend({}, this.layout, layout);
    }
});

//
// Requires a done callback and a data object to pass to the callback
//
function Wizard(options) {
    this.completed = options.completed;
    this.data = options.data;
    this.steps = [];
    this.currentIndex = 0;
    this.initialize();
}

$.extend(Wizard.prototype, {
    initialize: function() {
        this.el = $($("#wizard-template").html());
        this.tabs = this.el.find(".sections");
        this.next = this.el.find(".next");
        this.prev = this.el.find(".prev");
        this.done = this.el.find(".done");
        this.viewer = this.el.find("#step-container");
        this.notifications = this.el.find("#error_help_text");

        // jQuery events
        this.prev.unbind().click(this.movePrevious.bind(this));
        this.next.unbind().click(this.moveNext.bind(this));
        this.done.unbind().click(this.submit.bind(this));
    },

    at_first: function() {
        return this.currentIndex === 0;
    },

    at_last: function() {
        return this.currentIndex >= (this.steps.length - 1);
    },

    render: function() {
        var titles = this.steps.map(function(step) {
            return $("<div></div>", { text:  step.title });
        });

        this.tabs.html(titles);
        titles[this.currentIndex].addClass("active");

        var step = this.steps[this.currentIndex];
        if (step.render) {
            step.render();
        }
        this.viewer.html(step.el);

        if (this.at_first()) {
            this.prev.attr("disabled", 1);
        } else {
            this.prev.removeAttr("disabled");
        }

        if (this.at_last()) {
            this.next.attr("disabled", 1);
            this.done.removeAttr("disabled");
        } else {
            this.next.removeAttr("disabled");
            this.done.attr("disabled", 1);
        }
    },

    movePrevious: function() {
        if (!this.at_first()) {
            this.currentIndex--;
            this.render();
            this.notifications.stop(true, true).fadeOut(1500);
        }
    },

    moveNext: function() {
        var step = this.steps[this.currentIndex];

        if (!this.at_last() && step.is_valid()) {
            this.currentIndex++;
            this.render();
            this.notifications.stop(true, true).fadeOut(1500);
        }
    },

    message: function(message) {
        this.notifications.html(message)
            .show()
            .stop(true, true)
            .delay(10*1000)
            .fadeOut(1500);
    },

    submit: function() {
        if (!check_login()) {
            this.message('Your session has expired, please log in again.');
            return;
        }

        if (this.at_last()) {
            this.completed(this.data);
        }
    },

    // Expects a view with render, is_valid methods and a element property el
    addStep: function(step, index) {
        if (index !== undefined && index < this.steps.length) {
            this.steps.slice(index, 0, step);
        } else {
            this.steps.push(step);
        }
    }
});

function DataView(experiment) {
    this.experiment = experiment || {};
    this.title = "Data";
    this.files = [];
    this.initialize();
}

$.extend(DataView.prototype, {
    initialize: function() {
        this.el = $($("#data-template").html());
        this.file_selector = $($("#selector-template").html());
        this.selector_container = this.el.find("#selector_container");
    },

    render: function() {
        //FIXME: This selector should be in another view
        var selector = this.file_selector.clone();
        this.selector_container.empty();
        selector.appendTo(this.selector_container);
        selector.tabs();

        //FIXME: selector view should track the current path
        if (pageObj.current_path) {
            irods_get_path(pageObj.current_path);
        } else {
            irods_get_path();
        }

        selector.find('#input_url').bind('keyup focus click', function() {
            var button = selector.find("#ftp_get_button"),
                disabled = !selector.find("#input_url").val();

            button.toggleClass("ui-state-disabled", disabled);
        });

        selector.find('#input_accn').bind('keyup focus click', function() {
            var button = selector.find("#ncbi_get_button"),
                disabled = !selector.find("#input_accn").val();

            button.toggleClass("ui-state-disabled", disabled);
        });

        selector.find('#input_upload_file').fileupload({
            dataType: 'json',
            add: this.add.bind(this),
            done: this.uploaded.bind(this)
        });
    },

    add: function(e, data) {
        var filename = data.files[0].name;

        if ( !add_file_to_list(filename, 'file://'+filename) ) {
            alert('File already exists.');
        } else {
            // mdb 10/29/13 - moved from above to prevent stale load_id value, issue 236
            $('#input_upload_file').fileupload('option', { formData: {
                fname: 'upload_file',
                load_id: load_id
            }});

            data.submit();
        }
    },

    uploaded: function(e, data) {
        this.files.push({
            name: 'file://'+data.result.filename,
            path: data.result.path,
            size: units(data.result.size)
        });

        finish_file_in_list('file', 'file://'+data.result.filename, data.result.path, data.result.size);
    },

    is_valid: function() {
        var items = get_selected_files();

        if (!items.length) {
            return false;
        }

        var file_type = items[0].file_type = this.el.find("#select_file_type option:selected").val();

        if (items === null) {
            error_help('Files are still being transferred, please wait.');
            return false;
        }

        if (items.length === 0) {
            error_help('Please select a data file.');
            return false;
        }

        if (!file_type) {
            error_help("Please select the file type to continue");
            return false;
        }

        $.extend(current_experiment, {
            data: items
        });

        return true;
    }
});

function DescriptionView(experiment) {
    this.experiment = experiment;
    this.sources = undefined;
    this.title = "Describe your experiment";
    this.initialize();
}

$.extend(DescriptionView.prototype, {
    initialize: function() {
        this.el = $($("#description-template").html());
        this.edit_source = this.el.find("#edit_source");
        this.edit_genome = this.el.find("#edit_genome");

        this.get_sources();
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

                if (sources) {
                    self.edit_source.autocomplete({source: sources});
                }
            },
        });
    },

    render: function() {
        var self = this;

        // jQuery Events
        this.edit_genome.unbind().change(function() {
            // Reset gid when item has changed
            self.experiment.gid = undefined;
        });

        // jQuery UI
        this.edit_genome.autocomplete({
            source:[],
            select: function(event, ui) {
                $(this).val(ui.item.label);
                self.experiment.gid = ui.item.value;
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
        var version = this.el.find('#edit_version').val();
        var restricted = $('#restricted').is(':checked');
        var genome = this.el.find('#edit_genome').val();
        var aligner = this.el.find("#alignment").find(":checked").val();
        var ignore_cb = this.el.find('#ignore_missing_chrs');
        var ignore_missing_chrs = ignore_cb.is(':checked');

        if (!name) {
            error_help('Please specify an experiment name.');
            return false;
        }
        if (!version) {
            error_help('Please specify an experiment version.');
            return false;
        }

        var source = $('#edit_source').val();
        if (!source || source === 'Search') {
            error_help('Please specify a data source.');
            return false;
        }

        if (!genome || genome === 'Search' || !this.experiment.gid) {
            error_help('Please specify a genome.');
            return false;
        }

       $.extend(this.experiment, {
            description: {
                name: name,
                description: description,
                version: version,
                restricted: restricted,
                genome: genome,
            }
        });

        return true;
    },
});

function render_template(template, container) {
    container.empty()
        .hide()
        .append(template)
        .slideDown();
}

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
});

function AlignmentView() {
    this.initialize();
}

$.extend(AlignmentView.prototype, {
    initialize: function() {
        this.el = $($("#align-template").html());
    },

    is_valid: function() {
        return true;
    },
});

function QuantativeView(){
    this.initialize();
}

$.extend(QuantativeView.prototype, {
    initialize: function() {
        this.el = $($("#quant-template").html());
    },

    is_valid: function() {
        return true;
    },
});

function FastqView() {
    this.initialize();
}

$.extend(FastqView.prototype, {
    initialize: function() {
        this.el = $($("#fastq-template").html());
        this.snp_container = this.el.find("#snp-container");
        this.align_container = this.el.find("#align-container");
    },

    render: function() {
        var self = this;

        // jQuery events
        this.el.find("#snps").unbind().change(this.update_snp.bind(this));
        this.el.find("[name=aligner]").unbind().change(this.update_aligner.bind(this));
        this.el.find("#snp-method").unbind().change(function() {
            var selected = $(this).val();
            render_template(snp_templates[selected], self.snp_container);
        });
    },

    update_snp: function (ev) {
        var enabled = $(ev.target).is(":checked"),
            method = this.el.find("#snp-method");

        var el = $(document.getElementById(method.val()));

        if (enabled) {
            el.show();
            method.removeAttr("disabled");
            this.snp_container.slideDown();
            var initial = $("#snp-method").val();
            render_template(snp_templates[initial], this.snp_container);
        } else {
            method.attr("disabled", 1);
            this.snp_container.slideUp();
        }
    },

    update_aligner: function() {
        var selected = $("#alignment").find(":checked").val();
        render_template(align_templates[selected], this.align_container);
    },

    is_valid: function() {
        return true;
    },
});

function GeneralOptionsView() {
    this.initialize();
}

$.extend(GeneralOptionsView.prototype, {
    initialize: function() {
        this.el = $($("#general-options-template").html());
    },

    is_valid: function() {
        return true;
    },
});

function AdminOptionsView() {
    this.initialize();
}

$.extend(AdminOptionsView.prototype, {
    initialize: function() {
        this.el = $($("#admin-options-template").html());
    },

    render: function() {
        // jQuery UI
        this.el.find("#edit_user").unbind().autocomplete({
            source:[],
            focus: function() { return false; },
        });
    },

    is_valid: function() {
        return true;
    },
});

function OptionsView(experiment) {
    this.experiment = experiment;
    this.title = "Options";
    this.initialize();
}

$.extend(OptionsView.prototype, {
    initialize: function() {
        this.fastq_view = new FastqView();
        this.admin_view = new AdminOptionsView();
        this.general_view = new GeneralOptionsView();

        this.layout_view = new LayoutView({
            template: "#options-layout-template",

            layout: {
                "#admin-options": this.admin_view,
                "#analysis-options": this.fastq_view,
                "#general-options": this.general_view
            }
        });

        this.el = this.layout_view.el;
    },

    is_valid: function() {
        //if (!this.fastq_view.is_valid()) {
        //    return false;
        //}

        //if (!this.general_view.is_valid()) {
        //    return false;
        //}

        //if (!this.admin_view.is_valid()) {
        //    return false;
        //}

        return true;
    },

    render: function() {
        this.layout_view.renderLayout();
    },
});

function ConfirmationView(experiment) {
    this.experiment = experiment;
    this.initialize();
    this.title = "Review and Load";
}

$.extend(ConfirmationView.prototype, {
    initialize: function() {
        this.el = $($("#confirm-template").html());
        this.description = this.el.find(".confirm-description");
        this.data = this.el.find(".confirm-data");
        this.options = this.el.find(".confirm-options");
        this.pair_template = $($("#summary-pair-template").html());
    },

    render: function() {
        this.renderDescription(this.experiment.description);
        this.renderData(this.experiment.data);
        this.renderOptions(this.experiment.options);
    },

    renderDescription: function(description) {
        var key, newpair;
        this.description.empty();

        // Description Confirmation
        for(key in description) {
            if (description.hasOwnProperty(key)) {
                newpair = this.pair_template.clone();
                newpair.find(".name").html(key);
                newpair.find(".data").html(description[key]);
                this.description.append(newpair);
            }
        }
    },

    renderData: function(data) {
        var index, newpair;

        this.data.empty();

        for(index = 0; index < data.length; index++) {
            newpair = this.pair_template.clone();
            newpair.find(".name").html("File");
            newpair.find(".data").html(data[index].path);
            this.data.append(newpair);
        }
    },

    renderOptions: function(options) {
        this.options.empty();
    },

    is_valid: function() {
        return true;
    }
});

function load(experiment) {
    // Open status dialog right away - issue 101
    reset_log();
    $('#load_dialog').dialog('open');
    $('#load_log').html('Initializing ...');
    newLoad = true;

    var payload = $.extend({fname: "load_experiment"}, experiment);

    $.ajax({
        dataType: "json",
        type: "POST",
        data: JSON.stringify(payload),
        contentType: "application/json",
        success: function(obj) {
            if (obj && obj.error) {
                if (obj.error.PAYLOAD)  {
                    alert(obj.error.PAYLOAD);
                }
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
    });
}

function initialize_wizard() {
    var root = $("#wizard-container");
    var wizard = new Wizard({ completed: load, data: current_experiment });
    wizard.addStep(new DescriptionView(current_experiment));
    wizard.addStep(new DataView(current_experiment));
    wizard.addStep(new OptionsView(current_experiment));
    wizard.addStep(new ConfirmationView(current_experiment));
    wizard.render();

    // Create psudeo templates
    snp_templates = {
        coge: $($("#coge-snp-template").html()),
        gatk: $($("#samtools-snp-template").html()),
        platypus: $($("#platypus-snp-template").html()),
        samtools: $($("#gatk-snp-template").html())
    };
    align_templates = {
        gsnap: $($("#gsnap-template").html()),
        tophat: $($("#tophat-template").html())
    };

    root.append(wizard.el);
    return wizard;
}
