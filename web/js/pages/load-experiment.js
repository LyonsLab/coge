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

var SUPPORTED_FILES = concat.call(QUANT_FILES, ALIGN_FILES, SEQ_FILES, POLY_FILES);
var FILE_TYPE_PATTERNS = new RegExp("(:?" + SUPPORTED_FILES.join("|") + ")$");

// Returns the file extension detected or undefined
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

function check_login() { //TODO move to services.js
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

function file_selected(filename, url) {
    $('#select_file_button').hide();
    $('#select_file_type').show();
    $('#files').show();
}

function file_finished(size, url) {
    var files = get_selected_files();
    if (!files || !files.length)
    	return;
    
    var paths = files.map(function (item) { return item.path; });
    var file_type = autodetect_file_type(paths[0])

    if (file_type)
        $("#file_type_selector").val(file_type);
}

function file_canceled() {
    if (!get_selected_files()) {
// mdb removed 2/20/15 -- causes js error in Firefox
//    	$('#select_file_type')
//    		.hide()
//	        .find("option[value=autodetect")
//	        .prop("selected", true)
//	        .change();
    	
    	$('#files').hide();
    	
    	$("#select_file_type option:first").attr("selected", "selected");
    }
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

function wait_to_search (search_func, search_obj) {
    var search_term = search_obj.value;
    if (search_term && search_term.length >= 3) {
        if (pageObj.time) {
            clearTimeout(pageObj.time);
        }

        pageObj.time = setTimeout(
            function() {
                search_func(search_obj.value);
            },
            250
        );
    }
}

function search_genomes (search_term) {
    $.ajax({
        data: {
            fname: 'search_genomes',
            search_term: search_term,
            timestamp: new Date().getTime()
        },
        success : function(data) {
            var obj = jQuery.parseJSON(data);
            if (obj.items) {
                obj.items.forEach(function(element) {
                    element.label = element.label.replace(/&reg;/g, "\u00ae"); // (R)
                });
                $("#edit_genome").autocomplete({source: obj.items});
                $("#edit_genome").autocomplete("search");
            }
        }
    });
}

function search_users (search_term) {
    $.ajax({
        data: {
            fname: 'search_users',
            search_term: search_term,
            timestamp: new Date().getTime()
        },
        success : function(data) {
            var obj = jQuery.parseJSON(data);
            if (obj && obj.items) {
                $("#edit_user").autocomplete({source: obj.items});
                $("#edit_user").autocomplete("search");
            }
        }
    });
}

function search_notebooks (search_term) {
	coge.services.search_notebooks(search_term,USER_NAME,function(data){
		if (data.notebooks) {
			var items = [];
			data.notebooks.forEach(function(notebook) {
				items.push({label:notebook.name,value:notebook.id});
		    });
		    $("#edit_notebook").autocomplete({source: items});
		    $("#edit_notebook").autocomplete("search");
		}
	},function(){console.log('error');});
}

//function load_begin() {
//    // Open status dialog right away - issue 101
//    reset_log();
//    $('#load_dialog').dialog('open');
//    $('#load_log').html('Initializing ...');
//}

//function load_failed(logfile) {
//    // mdb added 6/24/14 - temporary message until JEX logging is improved
//    var msg =
//        '<div class="alert">' +
//        'The CoGe Support Team has been notified of this error but please ' +
//        'feel free to contact us at <a href="mailto:' + SUPPORT_EMAIL + '">' +
//        SUPPORT_EMAIL + '</a> ' +
//        'and we can help to determine the cause.' +
//        '</div>';
//    var log = $('#load_log');
//    log.html( log.html() + msg );
//
//    if (logfile) {
//        $("#logfile a").attr("href", logfile);
//        $('#logfile').fadeIn();
//    }
//
//    // Update dialog
//    $('#loading_msg').hide();
//    $('#error_msg,#cancel_load_experiment_button').fadeIn();
//
//    if (newLoad) { // mdb added check to prevent redundant emails, 8/14/14 issue 458
//        $.ajax({
//            data: {
//                fname: "send_error_report",
//                load_id: load_id,
//                job_id: WORKFLOW_ID
//            }
//        });
//    }
//}

//function load_succeeded() {
//    // Update dialog
//    $('#loading_msg').hide();
//    $('#finished_msg,#ok_button').fadeIn();
//}

//function progress_formatter(item) {
//    var msg;
//    var row = $('<li>'+ item.description + ' </li>');
//
//    var job_status = $('<span></span>');
//
//    if (item.status == 'scheduled')
//        job_status.append(item.status).addClass('down bold');
//    else if (item.status == 'completed')
//        job_status.append(item.status).addClass('completed bold');
//    else if (item.status == 'running')
//        job_status.append(item.status).addClass('running bold');
//    else if (item.status == 'skipped')
//        job_status.append("already generated").addClass('skipped bold');
//    else if (item.status == 'cancelled')
//        job_status.append(item.status).addClass('alert bold');
//    else if (item.status == 'failed')
//        job_status.append(item.status).addClass('alert bold');
//    else
//        return;
//
//    row.append(job_status);
//
//    if (item.elapsed)  {
//        row.append(" in " + coge.utils.toPrettyDuration(item.elapsed));
//    }
//
//    if (item.log) {
//        var p = item.log.split("\n");
//
//        var pElements = p.map(function(item) {
//            var norm = item.replace(/\\t/g, " ").replace(/\\'/g, "'");
//            return $("<div></div>").html(norm);
//        });
//
//        var log = $("<div></div>").html(pElements).addClass("padded");
//        row.append(log);
//    }
//
//    return row;
//}

//var workflow_start_time;
//function update_dialog(job_id, identifier, formatter) { //update_dialog(request, user, identifier, formatter) {
////    var get_status = function () { 
////        $.ajax({
////            type: 'GET',
////            url: request,
////            dataType: 'json',
////            data: {
////                username: user
////            },
////            success: update_callback,
////            error: update_callback,
////            xhrFields: {
////                withCredentials: true
////            }
////        });
////    };
//    
//    var update_callback = function(json) {
//        var dialog = $(identifier);
//        var workflow_status = $("<p></p>");
//        var log_content = $("<ul></ul>");
//        var results = [];
//        var current_status;
//        var refresh_interval = 2000;
//        var retry_interval = 5*1000;
//        
//        // Set refresh rate based on elapsed time
//        if (!workflow_start_time)
//        	workflow_start_time = new Date().getTime();
//        var run_time = new Date().getTime() - workflow_start_time;
//        if (run_time > 10*60*1000)
//        	refresh_interval = 60*1000;
//        else if (run_time > 5*60*1000)
//        	refresh_interval = 30*1000;
//        else if (run_time > 60*1000)
//        	refresh_interval = 15*1000;
//        //console.log('Refresh run_time=' + run_time + ' refresh_interval=' + refresh_interval);
//
//        var callback = function() {
//            update_dialog(job_id, identifier, formatter);
//        }
//
//        if (!json || json.error) {
//            pageObj.error++;
//            if ('Auth' in json.error) {
//            	dialog.find('#loading_msg').html('Login required to continue');
//            	dialog.find('#load_log')
//            		.css({'font-size': '1em'})
//            		.html("<br>Your session has expired.<br><br>" + 
//            			"Please log in again by clicking " +
//            			"<a onclick='login_cas();' style='font-weight:bold'>here</a>.");
//            	return;
//            }
//            else {
//            	//mdb removed 3/5/15 -- retry indefinitely
////	            if (pageObj.error > 3) {
////	                workflow_status.html('<span class=\"alert\">The job engine has failed.</span>');
////	                var logfile;
////		                if (json.results.length) {
////	                    logfile = json.results[0].path;
////	                }
////	                load_failed(logfile);
////	                return;
////	            }
//	            $("#load_error").html('Server not responding ('+pageObj.error+')').show();
//	            setTimeout(callback, retry_interval);
//	            return;
//            }
//        }
//        
//        pageObj.error = 0;
//        $("#load_error").hide();
//
//        if (json.status) {
//            current_status = json.status.toLowerCase();
//            workflow_status
//                .html("Workflow status: ")
//                .append( $('<span></span>').html(json.status) )
//                .addClass('bold');
//        } 
//        else {
//        	console.log('Error: missing status');
//            setTimeout(callback, refresh_interval);
//            return;
//        }
//
//        if (json.tasks) {
//            var jobs = json.tasks;
//            for (var index = 0; index < jobs.length; index++) {
//                var item = formatter(jobs[index]);
//                if (item)
//                    results.push(item);
//            }
//        }
//        else {
//        	console.log('Error: missing tasks');
//        }
//
//        // Sanity check -- progress dialog should be open
//        if (!dialog.dialog('isOpen')) {
//        	console.log('Error: progress dialog is closed');
//            return;
//        }
//
//        //FIXME Update when a workflow supports elapsed time
//        if (current_status == "completed") {
//            var total = json.tasks.reduce(function(a, b) {
//                if (!b.elapsed) return a;
//
//                return a + b.elapsed;
//            }, 0);
//
//            var duration = coge.utils.toPrettyDuration(total);
//
//            workflow_status.append("<br>Finished in " + duration);
//            workflow_status.find('span').addClass('completed');
//            if (json.results && json.results.length) 
//            	load_succeeded(json.results);
//        }
//        else if (current_status == "failed"
//                || current_status == "error"
//                || current_status == "terminated"
//                || current_status == "cancelled")
//        {
//            workflow_status.find('span').addClass('alert');
//
//            if (json.results && json.results.length)
//                logfile = json.results[0].path;
//            load_failed(logfile);
//        }
//        else if (current_status == "notfound") {
//        	console.log('Error: status is "notfound"');
//            setTimeout(callback, refresh_interval);
//            return;
//        }
//        else {
//            workflow_status.find('span').addClass('running');
//            setTimeout(callback, refresh_interval);
//        }
//
//        results.push(workflow_status);
//        log_content.append(results);
//        
//        if (json.results && json.results.length) {
//        	log_content.append("<div class='bold'>Here are the results:</div>");
//    	    json.results.forEach(function(result) {
//    	    	var html = format_result(result);
//    	    	log_content.append(html);
//    	    });
//        }
//        
//        dialog.find('#load_log').html(log_content);
//    };
//
//    //get_status();
//    coge.services.fetch_job(job_id, update_callback, update_callback);
//}

//function reset_log() {
//    $('#load_log').html('');
//    $('#load_link').html('');
//    $('#loading_msg').show();
//    $('#ok_button,#error_msg,#finished_msg,#finish_load_experiment_button,#cancel_load_experiment_button,#logfile').hide();
//}

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
    this.onCompleted = options.onCompleted;
    this.data = options.data;
    this.steps = [];
    this.currentIndex = 0;
    this.initialize(options);
}

$.extend(Wizard.prototype, {
    initialize: function(options) {
        this.el = $($("#wizard-template").html());
        this.tabs = this.el.find(".sections");
        this.next = this.el.find(".next");
        this.prev = this.el.find(".prev");
        this.done = this.el.find(".done");
        this.viewer = this.el.find("#step-container");
        this.notifications = this.el.find("#error_help_text");
        this.help = this.el.find(".link");

        // jQuery events
        this.prev.unbind().click(this.movePrevious.bind(this));
        this.next.unbind().click(this.moveNext.bind(this));
        this.done.unbind().click(this.submit.bind(this));
        this.help.unbind().click(function() {
        	window.open(options.helpUrl);
        });
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
//        for (var i=0;i<=this.currentIndex;i++)
//        	titles[i].addClass("active");
//        for (var i=0;i<this.currentIndex;i++) {
//        	var self = this;
//        	var index = i;
//        	titles[i].click(function(){self.move(index)});
//        	titles[i].css('cursor','pointer');
//        }

        var step = this.steps[this.currentIndex];
        if (step.render)
            step.render();
        this.viewer.html(step.el);

        if (this.at_first())
            this.prev.hide();//this.prev.attr("disabled", 1);
        else
            this.prev.removeAttr("disabled").show();

        if (this.at_last()) {
            this.next.hide();
            this.done.show();
        } else {
            this.next.show();
            this.done.hide();
        }
    },

    move: function(index) {
    	this.currentIndex = index;
    	this.render();
        this.notifications.stop(true, true).hide();
    },

    movePrevious: function() {
        if (!this.at_first())
        	this.move(this.currentIndex - 1);
    },

    moveNext: function() {
        var step = this.steps[this.currentIndex];

        if (!this.at_last() && step.is_valid())
        	this.move(this.currentIndex + 1);
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

        if (this.at_last())
            this.onCompleted(this.data);
    },

    // Expects a view with render, is_valid methods and a element property el
    addStep: function(step, index) {
        if (index !== undefined && index < this.steps.length)
            this.steps.slice(index, 0, step);
        else
            this.steps.push(step);
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
        if (pageObj.current_path)
            irods_get_path(pageObj.current_path);
        else
            irods_get_path();

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

    //FIXME: Add files to file list view and mark them as being transferred
    add: function(e, data) {
        var filename = data.files[0].name;

        if ( !add_file_to_list(filename, 'file://'+filename) ) {
            error_help('File already exists.');
        } 
        else {
            // mdb added 10/29/13 - prevent stale load_id value, COGE-236
            $('#input_upload_file').fileupload('option', { formData: {
                fname: 'upload_file',
                load_id: load_id
            }});

            data.submit();
        }
    },

    //FIXME: Update files to a completed state after transfer
    uploaded: function(e, data) {
        finish_file_in_list('file', 'file://'+data.result.filename, data.result.path, data.result.size);
    },

    //FIXME: Add multiple file support
    is_valid: function() {
        var items = get_selected_files();
        if (!items || items.length === 0) {
            error_help('Please select a valid data file.');
            return false;
        }

        items[0].file_type = this.el.find("#select_file_type option:selected").val();
        if (!items[0].file_type) {
            error_help("Please select the file type to continue");
            return false;
        }

        this.experiment.data = items;
        return true;
    }
});

function DescriptionView(opts) {
    this.experiment = opts.experiment;
    this.metadata = opts.metadata;
    this.gid = opts.gid;
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

        if (this.metadata) {
            this.el.find('#edit_name').val(this.metadata.name);
            this.el.find('#edit_description').val(this.metadata.description);
            this.el.find('#edit_version').val(this.metadata.version);
            this.edit_source.val(this.metadata.source);

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
        var version = this.el.find('#edit_version').val();
        var restricted = this.el.find('#restricted').is(':checked');
        var genome = this.el.find('#edit_genome').val();

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

        if (!genome || genome === 'Search' || !this.gid) {
            error_help('Please specify a genome.');
            return false;
        }

       $.extend(this.experiment, {
            metadata: {
                name: name,
                description: description,
                version: version,
                restricted: restricted,
                source: source,
                genome: genome,
            },

            gid: this.gid
        });

        return true;
    },
});

function render_template(template, container) {
    container.empty()
        .hide()
        .append(template)
        .show();//.slideDown();
}

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

function GeneralOptionsView() {
    this.data = {};
    this.initialize();
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

        console.log(this.data.notebook_type);
        if (this.data.notebook && this.data.notebook_type === "existing" && 
        		(!notebook || notebook === 'Search' || !this.notebook_id)) 
        {
            error_help('Please specify a notebook.');
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
    this.title = "Options";
    this.initialize();
}

$.extend(OptionsView.prototype, {
    initialize: function() {
        this.admin_view = new AdminOptionsView();
        this.general_view = new GeneralOptionsView();
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
            error_help("Please set the file type.");
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

function ConfirmationView(experiment) {
    this.experiment = experiment;
    this.initialize();
    this.title = "Review and Submit";
}

$.extend(ConfirmationView.prototype, {
    initialize: function() {
        this.el = $($("#confirm-template").html());
        this.description = this.el.find(".confirm-description");
        this.options = this.el.find(".confirm-options");
        this.data = this.el.find(".confirm-data");
        this.pair_template = $($("#summary-pair-template").html());
    },

    render: function() {
        this.renderDescription(this.experiment.metadata);
        this.renderData(this.experiment.data);
        this.renderOptions(this.experiment.options);
    },

    // Render description summary
    renderDescription: function(description) {
        this.description.empty();

        var key, newpair;
        for(key in description) {
            if (description.hasOwnProperty(key)) {
            	var value = description[key];
            	if (typeof value === "boolean")
            		value = (value ? "yes" : "no");
                newpair = this.pair_template.clone();
                newpair.find(".name").html(coge.utils.ucfirst(key));
                newpair.find(".data").html(value);
                this.description.append(newpair);
            }
        }
    },

    // Render data files summary
    renderData: function(data) {
        this.data.empty();

        var index, newpair;
        for(index = 0; index < data.length; index++) {
            newpair = this.pair_template.clone();
            newpair.find(".name").html("File");
            var filename = data[index].path.replace(/^.*[\\\/]/, '')
            newpair.find(".data").html(filename);
            this.data.append(newpair);
        }
    },

    // Render options summary
    renderOptions: function(options) {
        this.options.empty();

        var key, newpair;
        for(key in options) {
            if (options.hasOwnProperty(key)) {
            	var val;
            	if (typeof options[key] === 'object')
            		val = objToString(options[key]);
            	else if (typeof options[key] === 'boolean')
            		val = options[key] ? 'yes' : 'no';
            	else
                	val = String(options[key]);
                newpair = this.pair_template.clone();
                newpair.find(".name").html(coge.utils.ucfirst(key.replace('_', ' ')));
                newpair.find(".data").html(val);
                this.options.append(newpair);
            }
        }
    },

    // Validates the confirmation view (nothing to do here)
    is_valid: function() {
        return true;
    }
});

function objToString(obj) {
    var str = '<br>';
    for (var p in obj) {
        if (obj.hasOwnProperty(p)) {
            str += '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;' + p + ': ' + obj[p] + '<br>';
        }
    }
    return str;
}

function load(experiment) {
	// TODO make sure still logged in
	
	coge.progress.begin();
    newLoad = true;

	// Convert request into format for job service
	var request = {
		type: 'load_experiment',
		requester: {
			page: PAGE_NAME,
			user_name: USER_NAME
		},
		parameters: {
			gid: experiment.gid,
			metadata: experiment.metadata,
			alignment_params: experiment.options.alignment_params,
			trimming_params: experiment.options.trimming_params,
			expression_params: experiment.options.expression_params,
			snp_params: experiment.options.snp_params
		},
		options: {
			load_id: load_id,
			email: experiment.options.email,
			normalize: experiment.options.normalize,
			normalize_method: experiment.options.normalize_method,
			notebook: experiment.options.notebook,
			notebook_name: experiment.options.notebook_name,
			notebook_id: experiment.options.notebook_id,
			source_data: experiment.data
		}
	};
    
    coge.services.submit_job(request, 
    	function(response) { // success callback
    		if (!response) {
    			coge.progress.failed("Error: empty response from server");
    			return;
    		}
    		else if (!response.success || !response.id) {
    			coge.progress.failed("Error: failed to start workflow");
    			return;
    		}
    		
	        // Start status update
            window.history.pushState({}, "Title", "LoadExperiment.pl" + "?wid=" + response.id); // Add workflow id to browser URL
            coge.progress.update(response.id, response.site_url);
	    },
	    function(jqXHR, textStatus, errorThrown) { // error callback
	    	coge.progress.failed("Couldn't talk to the server: " + textStatus + ': ' + errorThrown);
	    }
	);
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
    wizard.addStep(new DescriptionView({
        experiment: current_experiment,
        metadata: opts.metadata,
        gid: opts.gid
    }));
    wizard.addStep(new DataView(current_experiment));
    wizard.addStep(new OptionsView({experiment: current_experiment, admin: opts.admin}));
    wizard.addStep(new ConfirmationView(current_experiment));

    // Render all the wizard sub views
    wizard.render();

    // Add the wizard to the document
    root.html(wizard.el);
    return wizard;
}
