/*global document,$,alert,window,JSON */

// Global experiment data
var current_experiment = {};

// Support file types
var QUANT_FILE_TYPES = /(:?csv|tsv|bed|gff|gtf)$/;
var POLY_FILE_TYPES = /(:?vcf)$/;
var ALIGN_FILE_TYPES =/(:?bam)$/;
var SEQ_FILE_TYPES = /(:?fastq|fq)$/;

// Template support
var snp_templates  = {};
var align_templates = {};

function autodetect_file_type(file) {
    var stripped_file = file.replace(/.gz$/, '');

    switch (true) {
        case QUANT_FILE_TYPES.test(stripped_file): return "quant";
        case POLY_FILE_TYPES.test(stripped_file): return "poly";
        case ALIGN_FILE_TYPES.test(stripped_file): return "align";
        case SEQ_FILE_TYPES.test(stripped_file): return "seq";
        default: return;
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

function description() {
    var name = $('#edit_name').val();
    var description = $('#edit_description').val();
    var version = $('#edit_version').val();
    var restricted = $('#restricted').is(':checked');
    var genome = $('#edit_genome').val();
    var gid = $('#gid').val();
    var aligner = $("#alignment").find(":checked").val();
    var ignore_cb = $('#ignore_missing_chrs');
    var ignore_missing_chrs = ignore_cb.is(':checked');

    if (!name) {
        error_help('Please specify an experiment name.');
        return;
    }
    if (!version) {
        error_help('Please specify an experiment version.');
        return;
    }

    var source = $('#edit_source').val();
    if (!source || source === 'Search') {
        error_help('Please specify a data source.');
        return;
    }

    if (!genome || genome === 'Search') {
        error_help('Please specify a genome.');
        return;
    }

    // Prevent concurrent executions - issue 101
    if ( $("#load_dialog").dialog( "isOpen" ) ) {
        return;
    }

    // Make sure user is still logged-in - issue 206
    if (!check_login()) {
        alert('Your session has expired, please log in again.');
        window.location.reload(true);
        return;
    }

    $.extend(current_experiment, {
        description: {
            name: name,
            description: description,
            version: version,
            restricted: restricted,
            genome: genome,
            gid: gid
        }
    });

    return true;
}

//FIXME: Add support for multiple files
function data() {
    var items = get_selected_files();
    var file_type = items[0].file_type = $("#select_file_type option:selected").val();

    if (items === null) {
        error_help('Files are still being transferred, please wait.');
        return;
    }

    if (items.length === 0) {
        error_help('Please select a data file.');
        return;
    }

    if (file_type === "autodetect") {
        file_type = items[0].file_type = autodetect_file_type(items[0].path);
    }

    if (!file_type) {
        return error_help("The file type could not be auto detected please select the filetype");
    }

    $('#fastq,#poly,#quant,#align').addClass('hidden');

    if (QUANT_FILE_TYPES.test(file_type)) {
        $("#quant").removeClass("hidden");
    }

    if (POLY_FILE_TYPES.test(file_type)) {
        $("#poly").removeClass("hidden");
    }

    if (SEQ_FILE_TYPES.test(file_type)) {
        $("#fastq").removeClass("hidden");
    }

    if (ALIGN_FILE_TYPES.test(file_type)) {
        $("#align").removeClass("hidden");
    }

    $.extend(current_experiment, {
        data: items
    });

    return true;
}

function options() {
    error_help("Test options");
    return true;
}

var setup_wizard = function () {
    var el = $("#wizard"),
        $next = el.find(".next"),
        $prev = el.find(".prev"),
        $done = el.find(".done"),
        header = el.find(".sections"),
        steps = [],
        currentIndex = 0;

    function noop() { return true; }

    el.find(".step").each(function(el) {
        var step = {
            el: $(this),
            title: $("<div></div>", { text:  $(this).attr('data-title') }),
            validateFn: window[$(this).attr('data-validate')] || noop,
            validated: false
        };

        step.el.hide();
        header.append(step.title);
        steps.push(step);
    });

    steps[currentIndex].el.show();
    steps[currentIndex].title.addClass("active");

    function my (options) {
        this.options = options || {};
        this.steps = [];

        $prev.click(this.prev.bind(this));
        $next.click(this.next.bind(this));
        $done.click(this.done.bind(this));

        this.initialize();
    }

    my.prototype = {
        initialize: function() {
            this.el = $($("#wizard-template").html());
            this.tabs = this.el.find(".sections");
            this.next = this.el.find(".next");
            this.prev = this.el.find(".prev");
            this.done = this.el.find(".done");
            this.viewer = this.el.find("#step-container");
            this.notifications = this.el.find("#error_help_text");
        },

        render: function() {
            console.log(this.tabs);
            var titles = this.steps.map(function(step) {
                return $("<div></div>", { text:  step.title });
            });

            this.tabs.html(titles);
            titles[0].addClass("active");
            this.viewer.html(this.steps[0].el);
        },

        prev: function() {
            var cur = steps[currentIndex],
                prev = steps[currentIndex - 1];

            if ((currentIndex - 1) >= 0) {
                cur.el.slideUp();
                cur.title.removeClass("active");
                prev.el.slideDown();
                prev.title.addClass("active");

                currentIndex--;
            } else {
                $prev.attr("disabled", 1);
            }

            if (currentIndex == 0) {
                $prev.attr("disabled", 1);
            }

            $next.removeAttr("disabled")
            $done.attr("disabled", 1);
        },

function Wizard(options) {
    this.options = options || {};
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
        }
    },

    moveNext: function() {
        var step = this.steps[this.currentIndex];

        if (!this.at_last() && step.is_valid()) {
            this.currentIndex++;
            this.render();
        }
    },

    message: function(message) {
        this.notifications.html(message)
            .show()
            .delay(10*1000)
            .fadeOut(1500);
    },

    submit: function() {
        if (this.at_last()) {
            this.options.success();
        }
    },

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
    this.initialize();
}

$.extend(DataView.prototype, {
    initialize: function() {
        this.el = $($("#data-template").html());
        this.file_selector = $($("#selector-template").html());
        this.selector_container = this.el.find("#selector_container");
    },

    render: function() {
        //FIXME: This selector should be another view
        var selector = this.file_selector.clone();
        this.selector_container.empty();
        selector.appendTo(this.selector_container);
        selector.tabs();
    },

    is_valid: function() {
        return true;
    }
});

function DescriptionView(experiment) {
    this.experiment = experiment || {};
    this.title = "Describe your experiment";
    this.initialize();
}

$.extend(DescriptionView.prototype, {
    initialize: function() {
        this.el = $($("#description-template").html());
    },

    is_valid: function() {
        return true;
        var name = this.el.find('#edit_name').val();
        var description = this.el.find('#edit_description').val();
        var version = this.el.find('#edit_version').val();
        var restricted = $('#restricted').is(':checked');
        var genome = this.el.find('#edit_genome').val();
        var gid = this.el.find('#gid').val();
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

        if (!genome || genome === 'Search') {
            error_help('Please specify a genome.');
            return false;
        }

        // Prevent concurrent executions - issue 101
        if ( $("#load_dialog").dialog( "isOpen" ) ) {
            return false;
        }

        // Make sure user is still logged-in - issue 206
        if (!check_login()) {
            alert('Your session has expired, please log in again.');
            window.location.reload(true);
            return false;
        }

        $.extend(this.experiment, {
            description: {
                name: name,
                description: description,
                version: version,
                restricted: restricted,
                genome: genome,
                gid: gid
            }
        });

        return true;
    }
});

function OptionsView() {
    this.initialize();
    this.title = "Options";
}

$.extend(OptionsView.prototype, {
    initialize: function() {
        this.el = $($("options-template").html());
    },
    is_valid: function() {
        return true;
    }
});

function ConfirmationView() {
    this.initialize();
    this.title = "Review and Load";
}

$.extend(ConfirmationView.prototype, {
    initialize: function() {
        this.el = $($("#confirm-template").html());
    },
    is_valid: function() {
        return true;
    }
});

function render_template(template, container) {
    container.empty()
        .hide()
        .append(template)
        .slideDown();
}

function update_snp(ev) {
    var enabled = $(ev.target).is(":checked"),
        method = $("#snp-method"),
        container = $("#snp-container");

    var el = $(document.getElementById(method.val()));

    if (enabled) {
        el.show();
        method.removeAttr("disabled");
        container.slideDown();

        method.unbind().change(function() {
            var selected = $("#snp-method").val();

            render_template(snp_templates[selected], container);
        });
    } else {
        method.attr("disabled", 1);
        container.slideUp();
    }
}

function update_aligner() {
    var selected = $("#alignment").find(":checked").val();
    var container = $("#align-container");

    render_template(align_templates[selected], container);
}

function initialize_wizard() {
    var root = $("#wizard-container");
    var wizard = new Wizard({ success: load_experiment });
    wizard.addStep(new DescriptionView());
    wizard.addStep(new DataView());
    wizard.addStep(new OptionsView());
    wizard.addStep(new ConfirmationView());
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

//    root.append(wizard.el);
    return wizard;
}
