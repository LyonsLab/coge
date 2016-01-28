var coge = coge || {};
coge.GenomeInfo = {};

// EXPORT DIALOG GLOBALS
var irods_home = $("<p>Sending to: " + irods_home_path + "</p>");
var export_error = $("<p></p>")
    .text("Failed to export to the CyVerse Data Store.")
    .addClass("alert");
var download_error = $("<p></p>")
    .text("The file could not be fetched")
    .addClass("alert");
var spinner = $("<img></img>").attr("src", "picts/ajax-loader.gif");
var message = $("<div></div>").text("Please wait...  ").append(spinner.clone());
var note = $("<p></p>").addClass("small")
        .text("(This may take several hours)");

coge.GenomeInfo.init = function(options) {
    $.ajaxSetup({
        type: "GET",
        url: options.PAGE_NAME,
        timeout: 86400,
        dataType: "html",
        cache: false,
    });

    coge.GenomeInfo.pageObj = options;

    $(".dialog_box").dialog({ autoOpen: false, width: 500 });

    $("#status_dialog").dialog({modal: true});

    $("#export_dialog").dialog({modal: true});

    set_annotation_table();

    init_annotation_dialog(coge.GenomeInfo.pageObj.genome_id, default_type);

    // Open status dialog
    if (coge.GenomeInfo.pageObj.job_id) {
        $('#loading_msg').hide();
        $('#load_dialog').dialog('open');
        update_dialog("jex/status/" + coge.GenomeInfo.pageObj.job_id, "#load_dialog", progress_formatter);
    }

    $("#edit_user").autocomplete({
        source:[],
        focus: function() { return false; },
    });
};

function get_gc_content(id, fname) {
    var elem = $(id);
    elem.hide().html(spinner.clone()).fadeIn();
    elem.attr("onclick", '')
        .unbind()
        .removeClass("link")
        .addClass("data5");

    $.ajax({
        dataType: "html",
        data: {
            fname: fname,
            dsgid: coge.GenomeInfo.pageObj.genome_id,
            gstid: $("#gstid").val()
        },
        success: function(data) {
            elem.fadeOut(function() {;
                elem.html(data).fadeIn();
            });
        }
    });
}

function get_feat_gc(opts){
    var elem = $('#gc_histogram');
    elem.dialog("option", "position", {
        my: "center",
        at: "top",
        of: $(".box")
    });
    elem.dialog("open");

    if (!opts){opts={};}
    chr = opts.chr;
    typeid = opts.typeid;
    min = $('#feat_gc_min').val();
    max = $('#feat_gc_max').val();
    hist_type = $('#feat_hist_type').val();
    dsid = opts.dsid;
    dsgid = $('#dsg_id').val();
    elem.html('loading. . .');

    $.ajax({
        data: {
                dsgid: coge.GenomeInfo.pageObj.genome_id,
                jquery_ajax: 1,
                    fname: 'get_gc_for_feature_type',
                dsid: dsid,
                typeid: typeid,
                chr: chr,
                min: min,
                max: max,
                hist_type: hist_type,
            },
            success: function(data) {
                elem.html(data);
            }
        });
}

function chr_hist (dsgid) {
    $('#chromosome_hist').dialog({
        autoOpen: true,
        position: {
           my: "center",
           at: "top",
           of: $(".box")
        }
    });

    $('#chromosome_hist').html('loading. . .');
    $.ajax({
        data: {
            jquery_ajax: 1,
            fname: 'get_chr_length_hist',
            dsgid: dsgid,
        },
        success: function (data) {$('#chromosome_hist').html(data);}
    });
}

function export_features_to_irods(id, feature, isGenome, type) {
    // Prevent concurrent executions
    if ( $("#export_dialog").dialog( "isOpen" ) )
        return;

    // Make sure user is still logged-in
    if (!check_and_report_login())
        return;

    // RESET DIALOG
    reset_dialog();

    // Open status dialog right away
    $('#export_dialog')
        .unbind()
        .dialog('open');

    var data = {
        fname: "export_features",
        protein: type,
        fid: feature
    };

    if (data.protein) {
        $('#export_log').html(
                "<br>Copying protein sequences to <br><br>" +
                '<a class="bold" target="_blank" href="http://data.iplantcollaborative.org/">' +
                irods_home_path + '</a>'
        );
    } else {
        $('#export_log').html(
                    "<br>Copying nucleotide sequences to <br><br>" +
                '<a class="bold" target="_blank" href="http://data.iplantcollaborative.org/">' +
                irods_home_path + '</a>'
        );

    }

    if (isGenome) {
        data.gid = id;
    } else {
        data.dsid = id;
    }

    $.ajax({
        data: data,
        dataType: "json",
        success: function(json) {
            if (json.error) {
                show_error();
            } else {
                $('#export_log').append('<br><br>File: ' + json.file);
                show_success();
            }
        }
    });
}

function wait_to_search (search_func, search_term) {
    coge.GenomeInfo.pageObj.search_term = search_term;

    if (coge.GenomeInfo.pageObj.time) {
        clearTimeout(coge.GenomeInfo.pageObj.time);
    }

    // FIXME: could generalize by passing select id instead of separate search_* functions
    coge.GenomeInfo.pageObj.time = setTimeout(
        function() {
            search_func(coge.GenomeInfo.pageObj.search_term);
        },
        500
    );
}

function search_organisms (search_term) {
    if (search_term.length > 2) {
        $.ajax({
            data: {
                jquery_ajax: 1,
                fname: 'search_organisms',
                search_term: search_term,
                timestamp: new Date().getTime()
            },
            success : function(data) {
                var obj = jQuery.parseJSON(data);
                if (obj && obj.items) {
                    $("#edit_organism").autocomplete({source: obj.items});
                    $("#edit_organism").autocomplete("search");
                }
            }
        });
    }
}

function get_features(selector) {
    var elem = $(selector);
    elem.html(spinner.clone())
        .attr("onclick", "")
        .unbind()
        .fadeIn();

    $.ajax({
        data: {
            fname: "get_features",
            dsgid: coge.GenomeInfo.pageObj.genome_id,
            gstid: $("#gstid").val()
        },
        success: function(data) {
            elem.fadeOut(function() {;
                elem.html(data).removeClass("small link padded ui-widget-content ui-corner-all").slideDown();
            });
        }
    })
};

function search_users (search_term) {
    $.ajax({
        data: {
            jquery_ajax: 1,
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

function get_genome_info () {
    $.ajax({
        data: {
            jquery_ajax: 1,
            fname: 'get_genome_info',
            gid: coge.GenomeInfo.pageObj.genome_id
        },
        success : function(data) {
            if (data) {
                $("#genome_info").html(data);
            }
        }
    });
}

function edit_genome_info () {
    $.ajax({
        data: {
            jquery_ajax: 1,
            fname: 'edit_genome_info',
            gid: coge.GenomeInfo.pageObj.genome_id
        },
        success : function(data) {
            if (data) {
                $("#genome_info_edit_box").html(data).dialog('open');
            }
        }
    });
}

function update_genome_info (){
    var org_name = $('#edit_organism').val();
    if (!org_name) {
        alert('Please select an organism.');
        return;
    }

    var version = $('#edit_version').val();
    if (!version) {
        alert('Please specify a genome version.');
        return;
    }

    var source_name = $('#edit_source').val();
    if (!source_name) {
        alert('Please specify a data source.');
        return;
    }

    var link = $('#edit_link').val();
    var name = $('#edit_name').val();
    var description = $('#edit_description').val();
    var type_id = $('#select_type').val();
    var restricted = $('#restricted').is(':checked');

    $.ajax({
        data: {
            jquery_ajax: 1,
            fname: 'update_genome_info',
            gid: coge.GenomeInfo.pageObj.genome_id,
            name: name,
            description: description,
            version: version,
            type_id: type_id,
            restricted: restricted,
            org_name: org_name,
            source_name: source_name,
            link: link,
            timestamp: new Date().getTime()
        },
        success : function(val) {
            if (val) {
                alert(val);
            }
            else {
                get_genome_info();
                $("#genome_info_edit_box").dialog('close');
            }
        }
    });
}

function delete_genome () {
    $.ajax({
        data: {
            fname: 'delete_genome',
            gid: coge.GenomeInfo.pageObj.genome_id,
        },
        success : function(rc) {
            if (rc) {
                location.reload();
            }
            else {
                alert('Error occurred!');
            }
        },
    });
}

function show_error() {
    $('#export_loading_msg').hide();
    $('#export_finished_msg').hide();
    $('#export_error_msg').fadeIn();
}

function show_success() {
    $('#export_ok_button').fadeIn();
    $('#export_loading_msg').hide();
    $('#export_finished_msg').fadeIn();
    $('#export_ok_button').fadeIn();
}

function get_gff() {
    var id = coge.GenomeInfo.pageObj.genome_id;
    $('#gff_export').dialog('option', 'width', 400).dialog('open');

    $('#gff_export').unbind().on("dialogclose", function(event, ui) {
        var id_type = $('#gff_id_type').val(),
            cds = +$('#cds_only')[0].checked,
            annos = +$('#annos')[0].checked,
            nu = +$('#name_unique')[0].checked,
            upa = +$('#upa')[0].checked;

        var height = $('#status_log')[0].scrollHeight;
        $("#status_log").animate({ scrollTop: height}, 500);

        var title = $("<h4>Generating gff file</h4>");

        // RESET DIALOG
        reset_dialog();

        $("#export_log")
            .html(title)
            .append(note)

        $("#export_dialog")
            .unbind()
            .dialog({ title: "Preparing download", autoOpen: true});

        $.ajax({
            dataType: "json",
            data: {
                fname: "get_gff",
                gid: id,
                id_type: id_type,
                cds: cds,
                annos: annos,
                nu: nu,
                upa: upa
            },
            success: function(json) {
                if (json.error) {
                    show_error();
                    $("#export_log").html(download_error);
                } else {
                    $('#export_log').append('<br><br>File: ' + json.file);
                    show_success();
                    $("#export_dialog").on("dialogclose", function(event, ui){
                        json.files.forEach(coge.utils.open);
                    });
                }
            }
        });
    });
}

function get_tbl() {
    var title = $("<h4>Generating tbl file</h4>");

    reset_dialog();

    $("#export_log")
        .html(title)
        .append(note);

    $("#export_dialog")
        .unbind()
        .dialog({ title: "Preparing download", autoOpen: true});

    $.ajax({
        dataType: "json",
        data: {
            fname: "get_tbl",
            gid: coge.GenomeInfo.pageObj.genome_id
        },
        success: function(json) {
            if (json.error) {
                show_error();
                $("#export_log").html(download_error);
            } else {
                $('#export_log').append('<br><br>File: ' + json.file);
                show_success();
                $("#export_dialog").on("dialogclose", function(event, ui){
                    json.files.forEach(coge.utils.open);
                });
            }
        }
    });
}

function get_bed() {
    var title = $("<h4>Generating bed file</h4>");

    reset_dialog();

    $("#export_log")
        .html(title)
        .append(note);

    $('#export_dialog')
        .unbind()
        .dialog({ title: "Preparing download", autoOpen: true});

    $.ajax({
        dataType: "json",
        data: {
            fname: "get_bed",
            gid: coge.GenomeInfo.pageObj.genome_id
        },
        success: function(json) {
            if (json.error) {
                $('#export_loading_msg').hide();
                $('#export_finished_msg').hide();
                $('#export_error_msg').fadeIn();
                $("#export_log").html(download_error);
            } else {
                $('#export_log').append('<br><br>File: ' + json.file);
                $('#export_ok_button').fadeIn();
                $('#export_loading_msg').hide();
                $('#export_finished_msg').fadeIn();
                $('#export_ok_button').fadeIn();
                $("#export_dialog").on("dialogclose", function(event, ui){
                    json.files.forEach(coge.utils.open);
                });
            }
        }
    });
}

function reset_dialog() {
    $('#export_finished_msg').hide();
    $('#export_ok_button').hide();
    $('#export_error_msg').hide();
    $('#export_loading_msg').show();
}

function export_gff() {
    $('#gff_export').dialog('option', 'width', 400).dialog('open');
    $('#gff_export').unbind().on("dialogclose", function() {

        // ANALYSIS OPTIONS
        var id = coge.GenomeInfo.pageObj.genome_id;
        var id_type = $('#gff_id_type').val(),
            cds = +$('#cds_only')[0].checked,
            annos = +$('#annos')[0].checked,
            nu = +$('#name_unique')[0].checked,
            upa = +$('#upa')[0].checked;

        // Prevent concurrent executions
        if ( $("#export_dialog").dialog( "isOpen" ) )
            return;

        // Make sure user is still logged-in
        if (!check_and_report_login())
            return;

        // RESET DIALOG
        reset_dialog();

        var title = $("<h4>Generating gff file</h4>");

        $("#export_log")
            .html(title)
            .append(note)

        $("#export_dialog")
            .unbind()
            .dialog({ title: "Exporting file", autoOpen: true});

        $.ajax({
            dataType: "json",
            data: {
                fname: "export_gff",
                gid: id,
                id_type: id_type,
                cds: cds,
                annos: annos,
                nu: nu,
                upa: upa
            },
            success: function(json) {
                $('#export_log').append('<br><br>File: ' + json.file);
                if (json.error) {
                    show_error();
                } else {
                    show_success();
                }
            }
        });
    });
}

function export_tbl() {
    // Prevent concurrent executions
    if ( $("#export_dialog").dialog( "isOpen" ) )
        return;

    // Make sure user is still logged-in
    if (!check_and_report_login())
        return;

    // RESET DIALOG
    reset_dialog();

    // Open status dialog right away
    $('#export_dialog')
        .unbind()
        .dialog('open');

    $('#export_log').html(
            "<br>Copying this genomes's TBL file to <br><br>" +
            '<a class="bold" target="_blank" href="http://data.iplantcollaborative.org/">' +
            irods_home_path + '</a>'
    );

    $.ajax({
        dataType: "json",
        data: {
            fname: "export_tbl",
            gid: coge.GenomeInfo.pageObj.genome_id
        },
        success: function(json) {
            $('#export_log').append('<br><br>File: ' + json.file);

            if (json.error) {
                show_error();
            } else {
                show_success();
            }
        }
    });
}

function export_bed() {
    // Prevent concurrent executions
    if ( $("#export_dialog").dialog( "isOpen" ) )
        return;

    // Make sure user is still logged-in
    if (!check_and_report_login())
        return;

    // RESET DIALOG
    reset_dialog();

    // Open status dialog right away
    $('#export_dialog')
        .unbind()
        .dialog('open');

    $('#export_log').html(
            "<br>Generating/copying this genomes's BED file to <br><br>" +
            '<a class="bold" target="_blank" href="http://data.iplantcollaborative.org/">' +
            irods_home_path + '</a>'
    );

    $.ajax({
        dataType: "json",
        data: {
            fname: "export_bed",
            gid: coge.GenomeInfo.pageObj.genome_id
        },
        success: function(json) {
            if (json.error) {
                show_error();
            } else {
                $('#export_log').append('<br><br>File: ' + json.file);
                show_success();
            }
        }
    });
}

function export_fasta_irods() {
    // Prevent concurrent executions
    if ( $("#export_dialog").dialog( "isOpen" ) )
        return;

    // Make sure user is still logged-in
    if (!check_and_report_login())
        return;

    // RESET DIALOG
    reset_dialog();

    // Open status dialog right away
    $('#export_dialog')
        .unbind()
        .dialog('open');

    $('#export_log').html(
            "<br>Copying this genomes's FASTA file to <br><br>" +
            '<a class="bold" target="_blank" href="http://data.iplantcollaborative.org/">' +
            irods_home_path + '</a>'
    );

    $.ajax({
        data: {
            fname: 'export_fasta_irods',
            gid: coge.GenomeInfo.pageObj.genome_id
        },
        success : function(filename) {
            if (filename) {  // finished successfully
                show_success();
                $('#export_log').append('<br><br>File: ' + filename);
            }
            else { // error occurred
                show_error();
            }
        }
    });
}

function reset_log() {
    $('#status_log').html('');
    $('#status_msg').show();
    $('#finished_msg,#error_msg,#ok_button').hide();
    $('#finish_button,#cancel_button').hide();
}

function reset_load() {
    window.history.pushState({}, "Title", PAGE_NAME + "?gid=" + coge.GenomeInfo.pageObj.genome_id);
    $('#load_dialog').dialog('close');
}

function load_failed(obj) {
    // Handle special case of genbank load of existing genome - FIXME: should this be here or miscopy from LoadGenome.tmpl?
    if ( obj && obj.links && obj.links.length ) {
        var link_text = obj.links.reduce(function(prev, cur) {
                return prev + '<a href="' + cur + '" target=_new>' + cur + '</a>' + '<br>';
            },
            '<b>Cannot load because the data already exist in the system:</b><br>');
        var log = $('#load_log');
        log.html( log.html() + link_text );
    }
	else { // mdb added 6/24/14 - temporary message until JEX logging is improved
		var msg =
			'<div class="alert">' +
			'We recently made some changes that limit our ability to report errors (more details <a href="https://genomevolution.org/wiki/index.php/Jex_Logging" target=_blank>here</a>). ' +
			'While we fix this please feel free to contact us at <a href="mailto:coge.genome@gmail.com">coge.genome@gmail.com</a> regarding this failure to load data ' +
			'and we can help to determine the cause.  Thanks!' +
			'</div>';
		var log = $('#load_log');
		log.html( log.html() + msg );
	}

    // Update dialog
    $('#loading_msg').hide();
    $('#error_msg').fadeIn();
    $('#cancel_button').fadeIn();

    $.ajax({
        data: {
            fname: "send_error_report",
            load_id: load_id,
            job_id: coge.GenomeInfo.pageObj.job_id
        }
    });
}

function load_succeeded(obj) {
    // Update globals
    coge.GenomeInfo.pageObj.genome_id = obj.coge.GenomeInfo.pageObj.genome_id; // for continuing to GenomeInfo

    // Update dialog
    $('#loading_msg').hide();
    $('#finished_msg,#finish_button,#ok_button').fadeIn();
}

function check_login() {
    var logged_in = false;

    $.ajax({
        async: false,
        data: {
            fname: 'check_login'
        },
        success : function(rc) {
            logged_in = rc;
        }
    });

    return logged_in;
}

function check_and_report_login() {
    if (!check_login()) {
        alert('Your session has expired, please log in again.');
        location.reload(true)
        return false;
    }
    return true;
}

function copy_genome(mask, seq_only) {
    // Prevent concurrent executions - issue 101
    if ( $("#status_dialog").dialog( "isOpen" ) )
        return;

    // Make sure user is still logged-in - issue 206
    check_and_report_login();

    // Open status dialog right away - issue 101
    reset_log();
    if (mask) {
        $('#load_dialog').dialog('option', 'title', 'Copying & Masking Genome');
    }
    else if (seq_only) {
        $('#load_dialog').dialog('option', 'title', 'Copying Genome (No Annotations)');
    }
    else {
        $('#load_dialog').dialog('option', 'title', 'Copying Genome');
    }
    $('#load_dialog').dialog('open');
//   $('#status_log').html('Initializing ...');

    $.ajax({
        dataType: "json",
        data: {
            fname: 'copy_genome',
            load_id: load_id,
            gid: coge.GenomeInfo.pageObj.genome_id,
            mask: mask,
            seq_only: seq_only,
            timestamp: new Date().getTime()
        },
        success : function(data) {
            if (data.error) {
                alert(data.error);
                return;
            }

            coge.GenomeInfo.pageObj.job_id = data.job_id;

            // Set link in status dialog
            $('#loading_msg span a').attr('href', data.link).html(data.link);

            // Add load_id to browser URL
            window.history.pushState({}, "Title", PAGE_NAME + "?gid=" + coge.GenomeInfo.pageObj.genome_id + "&job_id=" + data.job_id); // Add job_id to browser URL
            update_dialog("jex/status/" + data.job_id, "#load_dialog", progress_formatter);
        }
        // TODO: handle error, show in status dialog
    });
}

function progress_formatter(item) {
    var msg;
    var row = $('<li>'+ item.description + ' </li>');
    row.addClass('small');

    var job_status = $('<span></span>');

    if (item.status == 'scheduled')
        job_status.append(item.status).addClass('down bold');
    else if (item.status == 'completed')
        job_status.append(item.status).addClass('completed bold');
    else if (item.status == 'running')
        job_status.append(item.status).addClass('running bold');
    else if (item.status == 'skipped')
        job_status.append("already generated").addClass('skipped bold');
    else if (item.status == 'cancelled')
        job_status.append(item.status).addClass('alert bold');
    else if (item.status == 'failed')
        job_status.append(item.status).addClass('alert bold');
    else
        return;

    row.append(job_status);

    if (item.elapsed)  {
        row.append(" in " + coge.utils.toPrettyDuration(item.elapsed));
    }

    return row;
}

function update_dialog(request, identifier, formatter) {
    var get_status = function () {
        $.ajax({
            type: 'GET',
            url: request,
            dataType: 'json',
            success: update_callback,
            error: update_callback,
        });
    };

    var update_callback = function(json) {
        var dialog = $(identifier);
        var workflow_status = $("<p></p>");
        var data = $("<ul></ul>");
        var results = [];
        var current_status;
        var timeout = 2000;

        var callback = function() {
            update_dialog(request, identifier, formatter);
        }

        if (json.error) {
            coge.GenomeInfo.pageObj.error++;
            if (coge.GenomeInfo.pageObj.error > 3) {
                workflow_status.html('<span class=\"alert\">The job engine has failed.</span>');
                load_failed();
                return;
            }
        } else {
            coge.GenomeInfo.pageObj.error = 0;
        }

        if (json.status) {
            current_status = json.status.toLowerCase();
            workflow_status.html("Workflow status: ");
            workflow_status.append($('<span></span>').html(json.status));
            workflow_status.addClass('bold');
        } else {
            setTimeout(callback, timeout);
            return;
        }

        if (json.jobs) {
            var jobs = json.jobs;
            for (var index = 0; index < jobs.length; index++) {
                var item = formatter(jobs[index]);
                if (item) {
                    results.push(item);
                }
            }
        }

        if (!dialog.dialog('isOpen')) {
            return;
        }

        //FIXME Update when a workflow supports elapsed time
        if (current_status == "completed") {
            var total = json.jobs.reduce(function(a, b) {
                if (!b.elapsed) return a;

                return a + b.elapsed;
            }, 0);

            var duration = coge.utils.toPrettyDuration(total);

            workflow_status.append("<br>Finished in " + duration);
            workflow_status.find('span').addClass('completed');
            get_load_log(function(result) {
                load_succeeded(result);
            });

        }
        else if (current_status == "failed"
                || current_status == "error"
                || current_status == "terminated"
                || current_status == "cancelled")
        {
            workflow_status.find('span').addClass('alert');
            get_load_log(function(result) {
                load_failed(result);
            });
        }
        else if (current_status == "notfound") {
            setTimeout(callback, timeout);
            return;
        }
        else {
            workflow_status.find('span').addClass('running');
            setTimeout(callback, timeout);
        }

        results.push(workflow_status);
        data.append(results);
        dialog.find('#load_log').html(data);
    };

    get_status();
}

function get_load_log(callback) {
    $.ajax({
        dataType:    'json',
        data: {
            fname:       'get_load_log',
            workflow_id: coge.GenomeInfo.pageObj.job_id,
            timestamp:   new Date().getTime()
        },
        success : function(data) {
            if (callback) {
                callback(data);
                return;
            }
        }
    });
}

function continue_to_view() {
    window.location.href = "GenomeInfo.pl?gid=" + coge.GenomeInfo.pageObj.genome_id;
}

function set_annotation_table() {
    $('#genome_annotation_table').tablesorter({widgets: ['zebra']});
}

function get_annotations() {
    $.ajax({
        data: {
            fname: 'get_annotations',
            gid: coge.GenomeInfo.pageObj.genome_id,
        },
        success : function(data) {
            $('#genome_annotations').html(data);
            set_annotation_table();
        }
    });
}

function remove_annotation (gaid) {
    $.ajax({
        data: {
            fname: 'remove_annotation',
            gid: coge.GenomeInfo.pageObj.genome_id,
            gaid: gaid,
        },
        success : function() {
            get_annotations();
        },
    });
}

function get_datasets() {
    $.ajax({
        data: {
            fname: 'get_datasets',
            gid: coge.GenomeInfo.pageObj.genome_id,
        },
        success : function(data) {
            $('#datasets').html(data);
        }
    });
}

function delete_dataset (dsid) {
    $.ajax({
        data: {
            fname: 'delete_dataset',
            gid: coge.GenomeInfo.pageObj.genome_id,
            dsid: dsid,
        },
        success : function() {
            get_datasets();
        },
    });
}

function open_aa_usage_table(chromosome)
{

    $("#aa_usage_table")
        .html("Please wait... ")
        .append(spinner.clone())
        .dialog({
            autoOpen: true,
            position: {
                my: "center",
                at: "top",
                of: ".box"
            }
        });

    $.ajax({
        data: {
            fname: "get_aa_usage",
            dsgid: coge.GenomeInfo.pageObj.genome_id,
            chr: chromosome
        },
        dataType: "html",
        success: function(html) {
            $('#aa_usage_table').html(html);
            $('#aa_table').tablesorter();
        }
    });
}

function get_content_dialog(id, request, chromosome) {
    var elem = $(id);

    chromosome = (chromosome === "") ? undefined : chromosome;

    elem.dialog('option', 'position', {
        my: "center",
        at: "top",
        of: ".box"
    });

    elem.html("Please wait... ")
        .append(spinner.clone())
        .dialog('open')

    $.ajax({
        data: {
            fname: request,
            dsgid: coge.GenomeInfo.pageObj.genome_id,
            chr: chromosome,
        },
        dataType: "html",
        success: function(html) {
            elem.html(html).slideDown();
        }
    })
}

function update_owner () {
    var user_name = $('#edit_user').val();
    if (!user_name) {
        alert('Please specify a user.');
        return;
    }

    $.ajax({
        data: {
            jquery_ajax: 1,
            fname: 'update_owner',
            gid: coge.GenomeInfo.pageObj.genome_id,
            user_name: user_name,
            timestamp: new Date().getTime()
        },
        success : function(data) {
            if (data) {
                alert(data);
            }
        }
    });
}
