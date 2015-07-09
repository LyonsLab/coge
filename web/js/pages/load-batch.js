//function load_batch() {
//	// Validate data items
//	var items = get_selected_files();
//	if (items == null) {
//		error_help('Files are still being transferred, please wait.');
//		return;
//	}
//	else if (items.length == 0) {
//		error_help('Please select some sequence files by clicking <b>Add Data</b>.');
//		return;
//	}
//	var json = JSON.stringify(items);
//	
//	// Validate other fields
//    var name = $('#edit_name').val();
//    var description = $('#edit_description').val();
//
//    if (!genome_id) {
//    	error_help('Please specify a genome.');
//        return;
//    }
//    
//    if (!name) {
//    	error_help('Please specify a name.');
//        return;
//    }
//
//    var items = get_selected_files();
//    if (items == null) {
//    	error_help('Files are still being transferred, please wait.');
//        return;
//    }
//    else if (items.length == 0) {
//    	error_help('Please select a data file.');
//        return;
//    }
//
//    var assignee_user_name = $('#edit_user').val(); // input only exists if admin
//
//    // if Notebook field was cleared then ignore it and create new one based on name
//    if ( !$('#edit_notebook').val() )
//    	notebook_id = '';
//    
//    // Prevent concurrent executions - issue 101
//    if ( $("#load_dialog").dialog( "isOpen" ) )
//        return;
//
//    // Make sure user is still logged-in - issue 206
//    if (!check_login()) {
//        alert('Your session has expired, please log in again.');
//        location.reload(true)
//        return;
//    }
//
//    // Open status dialog right away - issue 101
//    reset_log();
//    $('#load_dialog').dialog('open');
//    $('#load_log').html('Initializing ...');
//    newLoad = true;
//
//    $.ajax({
//        data: {
//            fname: 'load_batch',
//            load_id: load_id,
//            name: name,
//            description: description,
//            gid: genome_id,
//            nid: notebook_id,
//            assignee_user_name: assignee_user_name,
//            items: json,
//            timestamp: new Date().getTime()
//        },
//        success : function(data) {
//            var obj = jQuery.parseJSON(data);
//            if (obj && obj.error) {
//                alert(obj.error);
//                return;
//            }
//
//            // Set link in status dialog
//            $('#loading_msg span a').attr('href', obj.link).html(obj.link);
//
//            // Start status update
//            if (obj.job_id) { // JEX status for load FASTQ
//                job_id = obj.job_id;
//                window.history.pushState({}, "Title", PAGE_NAME + "?job_id=" + obj.job_id); // Add job_id to browser URL
//                update_dialog(STATUS_URL + obj.job_id,  user_name, "#load_dialog", progress_formatter);
//            }
//        }
//        // TODO: handle error, show in status dialog
//    });
//}
//
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
//function load_failed(logfile){
//	// mdb added 6/24/14 - temporary message until JEX logging is improved
//	var msg =
//		'<div class="alert">' +
//		'The CoGe Support Team has been notified of this error but please ' + 
//		'feel free to contact us at <a href="mailto:<TMPL_VAR NAME=SUPPORT_EMAIL>"><TMPL_VAR NAME=SUPPORT_EMAIL></a> ' +
//		'and we can help to determine the cause.' +
//		'</div>';
//	var log = $('#load_log');
//	log.html( log.html() + msg );
//
//    if (logfile) {
//        $("#logfile a").attr("href", logfile);
//        $('#logfile').fadeIn();
//    }
//
//    // Update dialog
//    $('#loading_msg').hide();
//    $('#error_msg').fadeIn();
//    $('#cancel_load_button').fadeIn();
//
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
//	batch_id = obj.batch_id;  // for continuing to ExperimentView
//    notebook_id = obj.notebook_id;      // for continuing to NotebookView
//
//    // Update dialog
//
//    $('#loading_msg').hide();
//    $('#finished_msg').fadeIn();
//    $('#ok_button').fadeIn();
//    if (notebook_id) { // qTeller pipeline experiment load
//        $('#finish_load_button')
//            .html('NotebookView').fadeIn()
//            .unbind().on('click', function() {
//                window.location.href = "NotebookView.pl?nid=" + notebook_id;
//        });
//    }
//    else { // normal experiment load
//        $('#finish_load_button')
//            .html('ExperimentView').fadeIn()
//            .unbind().on('click', function() {
//                window.location.href = "ExperimentView.pl?eid=" + experiment_id;
//        });
//    }
//}

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
            this.el.find('#edit_name').val(this.metadata.name);
            this.el.find('#edit_description').val(this.metadata.description);
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
        var name = this.el.find('#edit_name').val();
        var description = this.el.find('#edit_description').val();
        var genome = this.el.find('#edit_genome').val();

        if (!name) {
            if (this.onError)
            	this.onError('Please specify an experiment name.');
            return false;
        }
        
        if (!genome || genome === 'Search' || !this.gid) {
        	if (this.onError)
            	this.onError('Please specify a genome.');
            return false;
        }

       $.extend(this.batch, {
            metadata: {
                name: name,
                description: description,
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

function wait_to_search (search_func, search_obj) {
    var search_term = search_obj.value;
    if (!search_term || search_term.length >= 2) {
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
        },
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
        },
    });
}

function search_notebooks (search_term) {
	coge.services.search_notebooks(
		search_term, 
		user_name, 
		function(obj) {
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
		}
	);
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
    
    coge.services.submit_job(request, 
    	function(response) { // success callback
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
    wizard.addStep(new DataView(current_batch, { supportedFileTypes: ['gz'], onError: wizard.error_help.bind(wizard) }));
    //wizard.addStep(new OptionsView({batch: current_batch, admin: opts.admin, onError: wizard.error_help }));
    wizard.addStep(new ConfirmationView(current_batch));

    // Render all the wizard sub views
    wizard.render();

    // Add the wizard to the document
    root.html(wizard.el);
    return wizard;
}