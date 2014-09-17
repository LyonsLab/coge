/*global $,jQuery,pageObj,document */

var concat = String.prototype.concat;
var cache = {};

function update_gapcost(pro) {
    $('.gapcosts').hide();
    if(pro)
    {
      var val = $('#matrix').val();
      $('#gapcosts_'+val).toggle();
    }
    else
    {
      var str = $('#match_score').val();
      var num1 = str.substr(0,1);
      var num2 = str.substr(2);
      $('#gapcosts_'+num1+num2).toggle();
    }
}

function update_gapcost_ncbi(pro) {
    $('.ncbi_gapcosts').hide();
    if(pro)
    {
      var val = $('#ncbi_matrix').val();
      $('#ncbi_gapcosts_'+val).toggle();
    }
    else
    {
      var str = $('#ncbi_match_score').val();
      var num1 = str.substr(0,1);
      var num2 = str.substr(2);
       $('#ncbi_gapcosts_'+num1+num2).toggle();
    }
}


function animate_params (html,version,pro){
    if(version === "coge_radio") {
        $('#pro_or_nu_param').hide(0).html(html).toggle();
        update_gapcost(pro);
    }
    else{
        $('#ncbi_pro_or_nu_param').hide(0).html(html).toggle();
        update_gapcost_ncbi(pro);
    }
}

function blast_param(blast_type, translate, version) {
    var cache_id,
        deferred = $.Deferred();

    //FIXME: This is hack to cover an issue with blast_param being called twice
    try {
        cache_id = concat.call(blast_type, translate, version);
    } catch(e) {
        console.error(e);
        deferred.reject('invalid arguments passed');
        return deferred.promise();
    }

    if (cache[cache_id]) {
        var entry = cache[cache_id];
        deferred.resolve(entry);
        animate_params(entry.html, entry.version, entry.pro);
        return deferred.promise();
    }

    return $.ajax({
        data: {
            fname: 'blast_param',
            blast_type: blast_type,
            translate: translate,
            version: version
        },
        success : function(data) {
            var obj = jQuery.parseJSON(data);
            cache[cache_id] = obj;
            animate_params(obj.html, obj.version, obj.pro);
        },
    });
}

function database_param(program) {
    var cache_id = program,
        deferred = $.Deferred(),
        entry;

    if (cache[cache_id]) {
        entry = cache[cache_id];
        deferred.resolve(entry);

        $('#database').html(entry);

        return deferred.promise();
    }

    return $.ajax({
        data: {
            fname: 'database_param',
            program: program,
        },
        success : function(data) {
            cache[cache_id] = data;
            $('#database').html(data);
        },
    });
}

function generate_seq_obj(dsid,dsgid, upstream,downstream,seqview,chr,rc,featid) {
    this.featid=featid;
    this.dsid=dsid;
    this.dsgid=dsgid;
    this.upstream=upstream;
    this.downstream=downstream;
    this.seqview=seqview;
    this.chr=chr;
    this.rc=rc;
    this.gstid= pageObj.gstid;
}

function update_info_box(featid) {
    //generate_feat_info(['args__'+featid],['feature_info_popup']);
    $.ajax({
        data: {
            fname: 'generate_feat_info',
            featid: featid,
        },
        success : function(data) {
            $('#feature_info_popup').html(data);
            $('#feature_info_popup').dialog('open');
        },
    });
}

function loading(id,msg) {
    var message = '<font class="loading">Loading '+msg+' . . .</font>';
    $('#'+id).html(message);
}

function update_hsp_info (featid) {
    loading('image_info','Information');
    loading('query_image','Image');
    loading('subject_image','Image');
    //get_hsp_info(['args__blastfile','args__'+pageObj.basename,'args__num','args__'+featid],['image_info','query_image','subject_image']);
    $.ajax({
        data: {
            fname: 'get_hsp_info',
            hsp_id: featid,
            blastfile: pageObj.basename
        },
        success : function(data) {
            var obj = jQuery.parseJSON(data);
            $('#image_info').html(obj.html);
            $('#query_image').html(obj.query_link);
            $('#subject_image').html(obj.subject_link);
            $('#result_visual_popup').dialog('open');
        },
    });
}

function update_checkbox(name, dist, hspid, id){
    $('#feat'+hspid).html(name);
    $('#dist'+hspid).html(dist);
    var id_array = [];
    if (id) {
        id_array = id.split(',');

        // mdb removed 5/19/14 issue 382 - jquery not working with html injection
        //$('#checkbox'+id_array[0]).attr("id","checkbox"+id_array[1]);
        //$('#'+id_array[0]).attr("value",id_array[1]).attr("id",id_array[1]);

        // mdb added 5/19/14 issue 382
        var old_id = id_array[0];
        var new_id = id_array[1];
        var e = document.getElementById('checkbox' + old_id);
        e.id = "checkbox" + new_id;
        e = document.getElementById(old_id);
        e.id = e.value = new_id;
    }
}

function init_table () {
    $.tablesorter.addParser({
        id: 'percent',
        is: function() { return false; },
        format: function(s) { return s.replace(/%/,''); },
        type: 'numeric'
    });

    $(function(){
        $("#hsp_result_table").tablesorter({
            //sortColumn: 'HSP#',               // Integer or String of the name of the column to sort by.
            sortClassAsc: 'headerSortUp',       // Class name for ascending sorting action to header
            sortClassDesc: 'headerSortDown',    // Class name for descending sorting action to header
            headerClass: 'header',              // Class name for headers (th's)
            widgets: ['zebra'],
            textExtraction: 'complex',
            headers: {0: {sorter: false},
                      4: {sorter: 'digit'},
                      5: {sorter: 'digit'},
                      6: {sorter: 'digit'},
                      7: {sorter: 'percent'},
                      8: {sorter: 'digit'},
                      9: {sorter: 'percent'},
                      10: {sorter: 'digit'},
                      11: {sorter: 'percent'}
            },
            sortList: [[2,0],[5,0]]
        });
    });
}

//Need to instantiate this seperately from other dialog boxes, need to do this AFTER results are generated
function init_table_opts() {
    //substaniate dialog box
    $("#table_opts_dialog").dialog({ height: 240,
                            width: 746,
                            autoOpen: false,
    });

    //button effects on events
    $('#table_opts').click(function() {
        $('#table_opts_dialog').dialog('open');
    });
}


// FIXME mdb 3/8/13 - instead of separate ajax requests, all HSP's should be handled in one
function fill_nearby_feats(id_array) { // mdb rewritten 3/8/13 issue 47
    //get_nearby_feats(['args__basefile','args__'+pageObj.basename,'args__num','args__'+id],[update_checkbox]);
    if (id_array.length) {
        var id = id_array.shift();
        $.ajax({
            data: {
                fname: 'get_nearby_feats',
                basefile: pageObj.basename,
                num: id
            },
            success: function(data) {
                var obj = jQuery.parseJSON(data);
                update_checkbox(obj.name, obj.distance, obj.hsp_id, obj.new_checkbox_info);
            },
            complete : function() {
                setTimeout(
                    function() {
                        fill_nearby_feats(id_array);
                    },
                    0
                );
            }
        });
    }
    else { // all done, init table sorting
        init_table();
        init_table_opts();
        //$("#hsp_result_table").trigger("update");
    }
}

function click_all_feat_links(feature_links) { // mdb rewritten 3/8/13 issue 47
    var link_array = feature_links.split(',');
    var id_array = [];

    link_array.forEach(
        function(element, index, array) {
            if (element) {
                id_array.push(element);
            }
        }
    );

    fill_nearby_feats(id_array);
}

function overlap_checkboxes() {
    var accn = "";
    $('#hsp_result_table :checkbox').each(function() {
        if (this.checked)
            accn = accn + this.id+",";
    });
    if (!accn || accn == ",") {
        alert("Please select one or more features.");
        return;
    }
    
    var action = $('#overlap_action').val();
    if (action == "gevo")
        overlap_feats_parse(accn);
    else if (action == "fasta")
        export_fasta_file(accn);
    else if (action == "seqview") {
        var locations = [];
        $('#hsp_result_table :checkbox').each(function(){
            if (this.checked) {
                var loc = $(this).parents('tr').find('.location').html();
                locations.push(loc);
            }
        });
        window.open("SeqView.pl?locations=" + locations.join(','));
    }
    else if (action == "phylo")
        export_fasta_file(accn);
    else if (action == "xls")
        export_to_excel(accn, pageObj.basename);
    else if (action == "tophits")
        export_top_hits(accn, pageObj.basename);
    else if(action == "tab")
        generate_tab_deliminated(accn, pageObj.basename);
    else if(action == "list")
        generate_feat_list(accn, pageObj.basename);
    else if(action == "blast")
        generate_blast(accn, pageObj.basename);
    else if(action == "hsp")
        export_hsp_info();
    else if (action == "CodeOn")
        export_CodeOn(accn);
    else
        alert("Internal Error: missing '"+action+"'");
}

function generate_blast(accn, filename) {
    $.ajax({
    	type: "POST",
        data: {
            fname: 'generate_blast',
            accn: accn,
            filename: filename,
        },
        success : function(data) {
            window.open(data);
        },
    });
}

function generate_feat_list(accn, filename) {
    $.ajax({
    	type: "POST",
        data: {
            fname: 'generate_feat_list',
            accn: accn,
            filename: filename,
        },
        success : function(data) {
            window.open(data);
        },
    });
}

function generate_tab_deliminated(accn, filename) {
    $.ajax({
    	type: "POST",
        data: {
            fname: 'generate_tab_deliminated',
            accn: accn,
            filename: filename,
        },
        success : function(data) {
            window.open(data);
        },
    });
}

function export_to_excel(accn, filename) {
    $.ajax({
    	type: "POST",
        data: {
            fname: 'export_to_excel',
            accn: accn,
            filename: filename,
        },
        success : function(data) {
            window.open(data);
        },
    });
}

function export_top_hits(accn, filename) {
    $.ajax({
    	type: "POST",
        data: {
            fname: 'export_top_hits',
            accn: accn,
            filename: filename,
        },
        success : function(data) {
            window.open(data);
        },
    });
}

function export_fasta_file(accn) {
    $.ajax({
    	type: "POST",
        data: {
            fname: 'export_fasta_file',
            accn: accn
        },
        success : function(data) {
            window.open(data);
        },
    });
}

function export_CodeOn(accn) {
    $.ajax({
    	type: "POST",
        data: {
            fname: 'export_CodeOn',
            accn: accn
        },
        success : function(data) {
            window.open(data);
        },
    });
}

function overlap_feats_parse(accn) {
    $.ajax({
    	type: "POST",
        data: {
            fname: 'overlap_feats_parse',
            accn: accn
        },
        success : function(data) {
            var obj = jQuery.parseJSON(data);
            if (obj) {
                if (obj.count > 10) {
                    var remove = obj.count - 10;
                    alert("You have exceeded the number of features you can send to GEvo ( 10 Max ). You currently have "+obj.count+" selected. Please uncheck "+remove+" of your checked item(s).");
                }
                else {
                    window.open(obj.url);
                }
            }
        },
    });
}

function export_hsp_info() {
    $.ajax({
        data: {
            fname: 'export_hsp_info',
            filename: pageObj.basename
        },
        success : function(data) {
            window.open(data);
        },
    });
}

function export_hsp_query_fasta() {
    $.ajax({
        data: {
            fname: 'export_hsp_query_fasta',
            filename: pageObj.basename
        },
        success : function(data) {
            window.open(data);
        },
    });
}

function export_hsp_subject_fasta(dna) {
    $.ajax({
        data: {
            fname: 'export_hsp_subject_fasta',
            filename: pageObj.basename,
            dna: dna,
        },
        success : function(data) {
            window.open(data);
        },
    });
}

function export_alignment_file() {
    $.ajax({
        data: {
            fname: 'export_alignment_file',
            filename: pageObj.basename,
        },
        success : function(data) {
            window.open(data);
        },
    });
}

function get_all_hsp_data() {
    export_hsp_info(); //export_hsp_info(['args__all','args__'+pageObj.basename],[export_file]);
}

function get_query_fasta() {
    export_hsp_query_fasta(); //export_hsp_query_fasta(['args__'+pageObj.basename],[export_file]);
}

function get_subject_fasta(val) {
    if (val) {
        export_hsp_subject_fasta(val); //export_hsp_subject_fasta(['args__'+pageObj.basename,'args__'+val],[export_file]);
    }
    else {
        export_hsp_subject_fasta(); //export_hsp_subject_fasta(['args__'+pageObj.basename],[export_file]);
    }
}

function get_alignment_file() {
    export_alignment_file(); //export_alignment_file(['args__'+pageObj.basename],[export_file]);
}

function enlarge_picture_window(id) {
    var html = $('#large_db'+id).html();
    html = html+'<br/>'+$('#large_pic'+id).html();
    $('#big_picture').html(html).dialog('open');
}

//TODO - clean this function up a bit
function show_seq(seq,name,num,dsid,chr,start,stop, rc) {
    if (num == 1) which_seq = "Query";
    else if (num == 2) which_seq = "Subject";
    else which_seq = "Alignment";
    if (seqObj){
        if (dsid == "seqObj") dsid = seqObj.dsid;
        if (chr == "seqObj") chr = seqObj.chr;
    }

    var html ="";
    hmtl =  html + "<TABLE width=100%><TR><TD valign='top' align='right'>"+which_seq+" Name: </TD><TD align=justify>"+name+"</TD></TR>";

    //*1 to make sure js sees integer
    if (stop*1 > start*1) {
        html = html + "<TR><TD valign='top' align='right'>"+which_seq+" Position: </TD><TD align=justify><pre>"+start+"-"+stop+"</pre></TD></TR>";
    }

    html = html + "<TR><TD valign='top' align='right'>"+which_seq+" Sequence: </TD><TD align=justify>"+seq+"</TD></TR></TABLE>";

    //check to see if SeqView is an option, then display link
    var regex = /^\d+$/;
    if (regex.test(dsid) && dsid != 0 && which_seq == "Subject") {
        var seqview = "SeqView.pl?dsid="+dsid+"&chr="+chr+"&start="+start+"&stop="+stop+"&rc="+rc;
        $('#sequence_popup').dialog('option','buttons',{"View Sequence in SeqView": function() { window.open(seqview); }});
    }
    $('#sequence_popup').html(html).dialog('open');
}

function get_params(){
    radio = get_radio("coge_radio","coge");
    var word_size = $('#word_size').val();
    var expect = $('#e_value').val();
    var match_mismatch = $('#match_score').val();
    var matrix = $('#matrix').val();
    var program = $('#'+radio).val();

    var gapcost;
    if (program == "blastn" || program == 'mega' || program == 'dcmega')
    {
        var num1 = match_mismatch.substr(0,1);
        var num2 = match_mismatch.substr(2);
        gapcost = $('#gapcosts_'+num1+num2).val();
    }
    else
    {
        gapcost = $('#gapcosts_'+matrix).val();
    }

    var filter_query = 0;
    if ($('#filter_query')[0].checked) {filter_query=1;}
    var reslimit = $('#resultslimit').val();
    //blastz parameters
    var zwordsize = $('#blastz_wordsize').val();
    var zgap_start = $('#blastz_gap_start').val();
    var zgap_extension = $('#blastz_gap_extension').val();
    var zchaining = $('#blastz_chaining').val();
    var zthreshold = $('#blastz_threshold').val();
    var zmask = $('#blastz_mask').val();

    var comp = $('#comp_adj').val();
    var seq = $('#seq_box').val();
    return {w : word_size, e : expect,g : gapcost,p : program,mm: match_mismatch,m : matrix,c : comp,s : seq, zw : zwordsize, zgs : zgap_start, zge : zgap_extension, zc : zchaining, zt : zthreshold, zm : zmask, fq : filter_query, rl : reslimit};
}

function get_ncbi_params(){
    radio = get_radio("ncbi_radio","ncbi");
    var word_size = $('#ncbi_word_size').val();
    var expect = $('#ncbi_e_value').val();
    var db = $('#ncbi_db').val();
    var match_mismatch = $('#ncbi_match_score').val();
    if (match_mismatch)
    {
        var num1 = match_mismatch.substr(0,1);
        var num2 = match_mismatch.substr(2);
    }
    var gapcost = $('#ncbi_gapcosts_'+num1+num2).val();
    var job_title = $('#job_title').val();
    var program = $('#'+radio).val();

    var matrix = $('#ncbi_matrix').val();
    var comp = $('#ncbi_comp_adj').val();
    var seq = $('#seq_box').val();

    return {w : word_size, e : expect, db : db,g : gapcost,j : job_title,p : program,mm: match_mismatch,m : matrix,c : comp,s : seq};
}

function reset_basename(){
    if(pageObj.basename) pageObj.basename=0;
}

function generate_basefile() {
    return $.ajax({
        data: {
            fname: 'generate_basefile',
        }
    });
}

function run_coge_blast() {
    reset_basename();

    generate_basefile().then(function(basename) {
        pageObj.basename = basename;
        blastOff("#status_dialog", "#results");
        ga('send', 'event', 'cogeblast', 'run', 'coge');
    });
}


function blastOff(dialog, results, basename) {
    var validator = $('#validator').hide();

    var check = $('#genome_choice').getLength();
    if (($('#blank').val())||(check==0)) {
        var msg = 'You have not selected any genomes to blast.';
        validator.html(msg).fadeIn();
        $('#coge_blast').toggle();

        return;
    }

    var params = get_params();

    program = params.p;
    expect = params.e;
    job_title = params.j;
    word_size = params.w;
    comp = params.c;
    matrix = params.m
    match_score = params.mm;
    gapcost = params.g;
    seq = params.s;
    seq = seq.replace(/__/g,/_/);

    if (!seq.length) {
        var msg = 'You have not specified a query sequence to blast.';
        validator.html(msg).fadeIn();
        return;
    }

    $('#results').slideUp();
    //$('#log_text').html('').slideDown();
    pageObj.nolog=1;

    filter_query = params.fq;
    resultslimit = params.rl;

    //blastz params
    var zwordsize = params.zw;
    var zgap_start = params.zgs;
    var zgap_extension = params.zge;
    var zchaining = params.zc;
    var zthreshold = params.zt;
    var zmask = params.zm;

    var page_width = ($('#co').width())*(95/100);
    pageObj.image_width = page_width;
    var blastable_db = $('#genome_choice').getLength(1);

    var _results = $(results);

    var status_dialog = $(dialog).dialog({
        title: 'Running CoGeBlast ...',
        modal: true,
        width: 500,
        autoOpen: false
    });

    // Close results then display dialog
    _results.slideUp(function() {
        status_dialog.dialog("open");
    });

    var options = {
        fname:          'blast_search',
        program:        program,
        expect:         expect,
        job_title:      job_title,
        wordsize:       word_size,
        comp:           comp,
        matrix:         matrix,
        matchscore:     match_score,
        gapcost:        gapcost,
        filter_query:   filter_query,
        resultslimit:   resultslimit,
        zwordsize:      zwordsize,
        zgap_start:     zgap_start,
        zgap_extension: zgap_extension,
        zchaining:      zchaining,
        zthreshold:     zthreshold,
        zmask:          zmask,
        basename:       pageObj.basename,
        seq:            seq,
        blastable:      blastable_db,
        fid:            pageObj.fid,
        width:          page_width,
        color_hsps:     $('#color_by').val(),
    };

    $.ajax({
        type: "POST",
        dataType: 'json',
        data: options,
        success : function(response) {
            //console.log(data);
            //if (data.error) {
            //    validator.html(data.error).fadeIn();
            //    $('#log_text').slideUp();
            //} else {
            //    //blastresults(data.html, data.click_all_links);
            //}

            status_dialog.unbind().on("dialogclose", function() {
                _results.removeClass('hidden').slideDown();

                // reset dialog
                status_dialog.find(".dialog-error,.dialog-complete").hide();
                status_dialog.find(".dialog-running").show();
                status_dialog.find(".dialog-text").html("");
                status_dialog.find(".dialog-link").html("");
                status_dialog.find(".dialog-log").html("");
            });

            if(response.success) {
                pageObj.runtime = 0;
                pageObj.fetch_error = 0;
                pageObj.error = 0;
                pageObj.engine = "<span class=\"alert\">The job engine has failed.</span><br>Please use the link below to use the previous version of SynMap.";

                var link = $("<a></a>")
                    .attr("href", response.link)
                    .html(response.link);

                var link_message = $("<span></span>")
                    .html("Return to this analysis: ")
                    .append(link);

                var logfile = $("<a></a>")
                    .attr("href", response.logfile)
                    .html("Logfile");

                status_dialog.find(".dialog-link").html(link_message);
                status_dialog.find(".dialog-log").html(logfile);

                options.fname = "get_results";
                options.logfile = response.logfile;
                delete options.seq;

                update_dialog("jex/status/" + response.id, status_dialog, results, formatter, options);
            } else {

                var error = $("<div></div>")
                    .addClass("alert")
                    .html(response.message);

                _results.append(results, error);
                status_dialog.find(".dialog-error").slideDown();
                status_dialog.find(".dialog-running").hide();
            }
        },
    });
}

function formatter(item) {
    var msg;
    var row = $('<li>'+ item.description + ' </li>');

    var job_status = $('<span></span>');

    if (item.status == 'scheduled') {
        job_status.append(item.status);
        job_status.addClass('down');
        job_status.addClass('bold');
    } else if (item.status == 'completed') {
        job_status.append(item.status);
        job_status.addClass('completed');
        job_status.addClass('bold');
    } else if (item.status == 'running') {
        job_status.append(item.status);
        job_status.addClass('running');
        job_status.addClass('bold');
    } else if (item.status == 'skipped') {
        job_status.append("already generated");
        job_status.addClass('skipped');
        job_status.addClass('bold');
    } else if (item.status == 'cancelled') {
        job_status.append(item.status);
        job_status.addClass('alert');
        job_status.addClass('bold');
    } else if (item.status == 'failed') {
        job_status.append(item.status);
        job_status.addClass('alert');
        job_status.addClass('bold');
    } else {
        return;
    }

    row.append(job_status);

    if (item.elapsed)  {
        row.append(" in " + coge.utils.toPrettyDuration(item.elapsed));
    }

    /*
    if (item.status == "skipped") {
        row.append("<p>The analyses previously was generated</p>");
    }
    */

    return row;
}

function update_dialog(request, identifier, result, formatter, args) {
    var get_status = function () {
        $.ajax({
            type: 'GET',
            url: request,
            dataType: 'json',
            success: update_callback,
            error: update_callback,
        });
    };

    var get_poll_rate = function() {
        pageObj.runtime += 1;

        if (pageObj.runtime <= 5) {
            return 1000;
        } else if (pageObj.runtime <= 60) {
            return 2000;
        } else if (pageObj.runtime <= 300) {
            return 5000;
        } else if (pageObj.runtime <= 1800) {
            return 30000;
        } else if (pageObj.runtime <= 10800) {
            return 60000;
        } else {
            return 300000;
        }
    };

    var showImmediately = function (started) {
        return (Date.now() - started) < 3000;
    };

    var fetch_results = function(completed, attempts) {
        dialog = $(identifier);

        var started = Date.now();

        $.ajax({
            data: args,
            dataType: "json",
            success: function(data) {
                if (completed && data.success) {
                    handle_results(result, data);

                    if (showImmediately(started)) {
                        dialog.dialog('close');
                    } else {
                        dialog.find('.dialog-running').hide();
                        dialog.find('.dialog-complete').slideDown();
                    }
                } else {
                    handle_results(result, data);
                    dialog.find('.dialog-running').hide();
                    dialog.find('.dialog-error').slideDown();
                }
            },
            error: function(data) {
                if (attempts >= 3) {
                    dialog.find('.dialog-running').hide();
                    dialog.find('.dialog-error').slideDown();
                } else {
                    var callback = function() {fetch_results(completed, attempts + 1)};
                    setTimeout(callback, 100);
                }
            }
        });
    }

    var update_callback = function(json) {
        var dialog = $(identifier);
        var workflow_status = $("<p></p>");
        var data = $("<ul></ul>");
        var results = [];
        var current_status;
        var timeout = get_poll_rate();

        var callback = function() {
            update_dialog(request, identifier, result, formatter, args);
        }

        if (json.error) {
            pageObj.error++;
            if (pageObj.error > 3) {
                workflow_status.html(pageObj.engine);
                dialog.find('.dialog-text').html(workflow_status);
                dialog.find('.dialog-running').hide();
                dialog.find('.dialog-error').slideDown();
                return;
            }
        } else {
            pageObj.error = 0;
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
            fetch_results(true, 1);
        } else if (current_status == "failed" || current_status == "error"
                || current_status == "terminated"
                || current_status == "cancelled") {
            workflow_status.find('span').addClass('alert');
            fetch_results(false, 1);
        } else if (current_status == "notfound") {
            setTimeout(callback, timeout);
            return;
        } else {
            workflow_status.find('span') .addClass('running');
            setTimeout(callback, timeout);
        }

        results.push(workflow_status);
        data.append(results);
        dialog.find('.dialog-text').html(data);
    };

    get_status();
}

function handle_results(selector, data) {
    var feature_links = data.click_all_links || "",
        results = $(selector).html(data.html);

    if ( $('#null').html() == "null" ) {
        $('#result_toggle').hide();
        $('#result_features').hide();
        return;
    }
    else {
        var table_height = $('#hsp_result_table_body').height();
        var max_height = 400;
        if (table_height > max_height) $('#hsp_result_table_body').height(max_height);
    }

    setup_button_states();
    click_all_feat_links(feature_links);
    check_display();
}

function ncbi_blast(url) {
    if (url == 1) {
        alert('You have not selected a BLAST program! Please select a program to run.');
    }
    else {
        var params = get_ncbi_params();
        seq = params.s;
        seq = seq.replace(/\n/,'%0D');
        db = params.db;
        expect = params.e;
        job_title = params.j;
        word_size = params.w;
        comp = params.c;
        matrix = params.m
        match_score = params.mm;

        url = url+'&DATABASE='+db+'&EXPECT='+expect+'&QUERY='+seq+'&JOB_TITLE='+job_title+'&WORD_SIZE='+word_size;

        if (program == 'blastn') {
            window.open(url+'&MATCH_SCORES='+match_score);
        }
        else {
            if ($('#comp_adj').is(':hidden')) {window.open(url+'&MATRIX_NAME='+matrix);}
            else {window.open(url+'&MATRIX_NAME='+matrix+'&COMPOSITION_BASED_STATISTICS='+comp);}
        }
    }
}

$.fn.getLength = function(val){
  var opt_length;
  var blastable;
  var blanked=0; //otherwise get math problems later...boo javascript
  this.each(
    function()
    {
        var opts = this.options;
        opt_length = opts.length;
        if (opt_length == 0) {return opt_length;}
        blastable = opts[0].id;
        if (blastable == 'blank') {blanked++;} //Counts the number of instances of blank
        if (val){
          for(var i=1;i<opts.length;i++) //>
          {
            blastable += ","+opts[i].id;
            //need to chopoff last comma
          }
        }
    }
  );
  if(val) return blastable;
  if (blanked) {opt_length-=blanked;} //subtract elements that are classified as blank above
  return opt_length;
};

$.fn.sortSelect = function(){
  this.each(
      function()
      {
        if(this.nodeName.toLowerCase() != "select"){
          return;}
        var opts = this.options;
        var sortArray = [];
        for(var i=0;i<opts.length;i++) //>
        {
            sortArray[i] = {
                    v: opts[i].value,
                    t: opts[i].text,
                    d: opts[i].id,
                    }
        }
        sortArray.sort(
                function(obj1,obj2)
                {
                    obj1t = obj1.t.toLowerCase(),
                    obj2t = obj2.t.toLowerCase();
                    if(obj1t == obj2t){
                      return 0;}
                    return obj1t < obj2t ? -1 : 1; //>
                }
        );
        for(var i=0;i<opts.length;i++) //>
        {
            opts[i].id = sortArray[i].d;
            opts[i].text = sortArray[i].t;
            opts[i].value = sortArray[i].v;
        }
      }
    );
    return this;
};

function select_blast() {
    var radio = get_radio('ncbi_radio','ncbi');
    //get_url([radio],[ncbi_blast]);
    $.ajax({
        data: {
            fname: 'get_url',
            program: $('#'+radio).val()
        },
        success : function(data) {
            ncbi_blast(data);
        },
    });
}

function get_radio(which_type,val){
    if ($('#'+which_type)[0].checked) { return val+"_blast_type_n"; }
    else { return val+"_blast_type_p"; }
}

function get_seq(which_type) {
    var cache_id,
        dsid = seqObj.dsid,
        dsgid = seqObj.dsgid,
        featid = seqObj.featid,
        chr = seqObj.chr,
        deferred = $.Deferred(),
        program = get_radio(which_type,"coge");

    if (featid) {
        cache_id = concat.call(featid, program, seqObj.upstream,
                               seqObj.downstream, seqObj.rc, seqObj.gstid);

        if (cache[cache_id]) {
            deferred.resolve(cache[cache_id]);
            $('#seq_box').val(cache[cache_id]);
            return deferred.promise();
        }

        $('#seq_box').val('Loading ...');
        return $.ajax({
            data: {
                fname: 'get_sequence',
                fid: featid,
                blast_type: program,
                upstream: seqObj.upstream,
                downstream: seqObj.downstream,
                rc: seqObj.rc,
                gstid: seqObj.gstid
            },
            success : function(html) {
                cache[cache_id] = html;
                $('#seq_box').val(html);
            },
        });
    }
    else if (chr) {
        $('#seq_box').val('Loading ...');

        cache_id = concat.call(dsid, dsgid, program, seqObj.upstream,
                               seqObj.downstream, seqObj.gstid);

        if (cache[cache_id]) {
            deferred.resolve(cache[cache_id]);
            $('#seq_box').val(cache[cache_id]);
            return deferred.promise();
        }

        return $.ajax({
            data: {
                fname: 'get_sequence',
                dsid: dsid,
                dsgid: dsgid,
                blast_type: program,
                start: seqObj.upstream,
                stop: seqObj.downstream,
                gstid: seqObj.gstid
            },
            success : function(html) {
                cache[cache_id] = html;
                $('#seq_box').val(html);
            },
        });
    }
    else if (pageObj.locations) {
        $('#seq_box').val('Loading ...');
        return $.ajax({
            data: {
                fname: 'get_sequence',
                blast_type: program,
                locations: pageObj.locations
            },
            success : function(html) {
                $('#seq_box').val(html);
            },
        });
    }
    else {
        deferred.resolve("");
        return deferred.promise();
        //$('#seq_box').val('');
        //return;
    }
}

function blast_param_on_select(which_type, val) {
    radio = get_radio(which_type, val);
    program = $('#'+radio).val();
    database_param(program); //database_param([radio], ['database']);

    $('#blast_parameters').hide();
    $('#blastz_parameters').hide();

    if (program == 'lastz') {
        $('#blastz_parameters').toggle();
    }
    else {
        $('#blast_parameters').toggle();
        if ((program == 'blastx') || (program == 'tblastx')) {
            blast_param("blast_type_p", 1, which_type); //blast_param(['args__blast_type','args__'+"blast_type_p",'args__translate','args__1','args__version','args__'+which_type],[animate_params]);
            $('#word_size').val(3);
        }
        else if ((program == 'blastp') || (program == 'tblastn')) {
            blast_param("blast_type_p", 0, which_type); //blast_param(['args__blast_type','args__'+"blast_type_p",'args__version','args__'+which_type],[animate_params]);
            $('#word_size').val(3);
        }
        else {
            blast_param('', 0, which_type); //blast_param(['args__version','args__'+which_type],[animate_params]);
            if (program == "dcmega") {
                $('#word_size').val(11);
            }
            else {
                $('#word_size').val(8);
            }
        }
    }
}

function org_search(desc_search){
    if (pageObj.time) {
        clearTimeout(pageObj.time);
    }

    name_desc = $('#org_name_desc').val();
    if (name_desc.length < 3) { return; } //>

    pageObj.time = setTimeout(
        function() {
            $("#wait_indicator").css({opacity:1});
            $.ajax({
                data: {
                    fname: 'get_orgs',
                    name_desc: name_desc,
                    timestamp: new Date().getTime()
                },
                success : function(data) {
                    var items = jQuery.parseJSON(data);
                    if (!pageObj.timestamp_get_orgs || items.timestamp > pageObj.timestamp_get_orgs) {
                        pageObj.timestamp_get_orgs = items.timestamp;
                        $("#wait_indicator").animate({opacity:0});
                        $('#org_id').html(items.html);
                        seq_type_search();
                        update_buttons();
                        count_organisms();
                    }
                }
            });
        },
        500
    );
}


function select_by_value($elements, property, value) {
    return $elements.filter(function() {
        return this.value === String(value);
    }).prop(property, true);
}

function Blast(params) {
    this.params = params || {};

    this.defaults = {
        type: 'coge_blast_type_n',
        match_score: null,
        evalue: 1e-3,
        program: null
    };

    this.params = $.extend(this.defaults, this.params);
};

$.extend(Blast.prototype, {
    _select_type: function () {
        var elements = $('input[name="cogeblast"]');
        select_by_value(elements, 'checked', this.params['type']);
    },

    _select_program: function () {
        $("#" + this.params['type']).val(this.params['program']);
    },

    _select_color_by: function () {
        var elements = $('input[name="color_by"]');
        select_by_value(elements, 'checked', this.params['color']);
    },

    _select_match_score: function () {
        var elements = $('#match_score option');
        select_by_value(elements, 'selected', this.params['match_score']);
    },

    _select_evalue: function () {
        var elements = $('#e_value option');
        select_by_value(elements, 'selected', this.params['evalue']);
    },

    update_display: function () {
        // Select the blast type (nucleotide vs protein)
        this._select_type();

        // Set the blast tool being used
        this._select_program();

        // Select the blast hit coloring scheme
        this._select_color_by();

        // Set the match score to be used
        this._select_match_score();

        // Set the e-value parameter
        this._select_evalue();
    }
});

function adjust_blast_types(val){
    if(val == 1){
        if($('#ncbi_blast_type').is(":hidden")){
            $('#coge_blast_type').hide(0);
            $('#ncbi_blast_type').show(0);
            get_seq('ncbi_radio');
        }
        else return;
    }
    else{
        if($('#coge_blast_type').is(":hidden")){
            $('#coge_blast_type').show(0);
            $('#ncbi_blast_type').hide(0);
            get_seq('coge_radio');
        }
        else return;
    }
}

function matrix_view (){
    var matrix = $('#matrix').val();
    window.open('MatrixView.pl?matrix='+matrix);
}

function toggle_hsp_column(index) {
    show = 0;
    //console.log ($('#show_columns td:eq('+(1*index-1)+')').html());
    if ($('#show_columns td:eq('+(1*index-1)+')').children()[0] && $('#show_columns td:eq('+(1*index-1)+')').children()[0].checked) { show=1; }
    if (show) {
        $('#hsp_result_table td:nth-child('+(1*index+1)+')').show();
        $('#hsp_result_table th:nth-child('+(1*index+1)+')').show();
    }
    else {
        $('#hsp_result_table td:nth-child('+(1*index+1)+')').hide();
        $('#hsp_result_table th:nth-child('+(1*index+1)+')').hide();
    }
}

function check_display() {
    i=1;
    $('#show_columns td').each(function() {
        if (!$(this).children()[0].checked){ toggle_hsp_column(i);}
        i++;
    });
}

function save_display_settings() {
    var i=1;
    var index;
    $(".hsp_display").each(function(){
        if($(this)[0].checked) {
            if (index) {index= index +","+i;}
            else {index = i;}
        }
        i++;
    });
    //save_settings_cogeblast(['args__display','args__'+index],[]);
    $.ajax({
        data: {
            fname: 'save_settings',
            display: index,
        },
        success: function() {
            $('#table_opts_dialog').dialog('close');
        }
    });
}

function toggle_display(button_obj, target_id) {
    var label = $(button_obj).html();

    if (label == 'hide') {
        $(button_obj).html('show');
        $('#'+target_id).slideUp('fast');
    }
    else {
        $(button_obj).html('hide');
        $('#'+target_id).slideDown('fast');
    }
}
