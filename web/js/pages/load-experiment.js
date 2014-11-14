/*global document,$,alert,window,JSON */

// Support file types
var QUANT_FILE_TYPES = /(:?csv|tsv|bed|gff|gtf)$/;
var POLY_FILE_TYPES = /(:?vcf)$/;
var ALIGN_FILE_TYPES =/(:?bam)$/;
var SEQ_FILE_TYPES = /(:?fastq|fq)$/;

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

    return true;
}

//FIXME: Add support for multiple files
function data() {
    var items = get_selected_files();
    var file_type = $("#select_file_type option:selected").val();

    if (items === null) {
        error_help('Files are still being transferred, please wait.');
        return;
    }

    if (items.length === 0) {
        error_help('Please select a data file.');
        return;
    }

    file_type = autodetect_file_type(items[0].path);

    if (!file_type) {
        return error_help("The file type could not be auto detected please select the filetype");
    }

    var json = JSON.stringify(items);

    return true;
}

function options() {
    error_help("Test options");
    return true;
}

function setup_wizard(options) {
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

    function my () {}

    my.prototype = (function() {
        return {
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

            next: function(options) {
                options = (options || {});

                var cur = steps[currentIndex],
                    next = steps[currentIndex + 1];

                if (options.force || (cur.validated = cur.validateFn())) {
                } else {
                    return;
                }

                if ((currentIndex + 1) < steps.length) {
                    cur.el.slideUp();
                    cur.title.removeClass("active");
                    next.el.slideDown();
                    next.title.addClass("active");

                    currentIndex++;
                } else {
                    $next.attr("disabled", 1);
                    $done.removeAttr("disabled");
                }

                if (currentIndex == (steps.length - 1)) {
                    $next.attr("disabled", 1);
                    $done.removeAttr("disabled");
                }

                $prev.removeAttr("disabled")
            },

            done: function() {
                if ((currentIndex + 1) === steps.length) {
                    options.success();
                }
            }
        };
    }());

    $prev.click(my.prototype.prev);
    $next.click(my.prototype.next);
    $done.click(my.prototype.done);

    return my;
}

function update_snp(ev) {
    var enabled = $(ev.target).is(":checked"),
        method = $("#snp-method"),
        wrapper = $("#snp-wrapper");

    var el = $(document.getElementById(method.val()));

    if (enabled) {
        $("#coge,#platypus,#samtools,#gatk").hide();
        el.show();
        method.removeAttr("disabled");
        wrapper.slideDown();

        method.unbind().change(function() {
            var selected = $("#snp-method").val();

            $("#coge,#platypus,#samtools,#gatk").hide();
            var el = $(document.getElementById(selected));
            el.slideDown();
        });
    } else {
        method.attr("disabled", 1);
        wrapper.slideUp();
    }
}

function update_aligner() {
    var aligner = $("#alignment").find(":checked").val();
    $("#gsnap,#tophat").hide();
    $(document.getElementById(aligner)).slideDown();
}

function initialize() {
    var Wizard = setup_wizard({ success: load_experiment });
    wizard = new Wizard();
}
