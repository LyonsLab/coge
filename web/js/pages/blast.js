/*global $,jQuery,pageObj,document */

var concat = String.prototype.concat;
var cache = {};

function update_gapcost(pro) {
    var root = $("#coge-params");
    root.find('.gapcosts').hide();
    if(pro) {
      var val = root.find('#matrix').val();
      root.find('#gapcosts_'+val).toggle();
    } else {
      var str = root.find('#match_score').val();
      var num1 = str.substr(0,1);
      var num2 = str.substr(2);
      root.find('#gapcosts_'+num1+num2).toggle();
    }
}

function update_gapcost_ncbi(pro) {
    var root = $("#ncbi-params");
    root.find('.gapcosts').hide();

    if(pro) {
      var val = root.find('#matrix').val();
      root.find('#gapcosts_'+val).toggle();
    } else {
      var str = root.find('#match_score').val();
      var num1 = str.substr(0,1);
      var num2 = str.substr(2);
      root.find('#gapcosts_'+num1+num2).toggle();
    }
}

function animate_params (html,version,pro){
    if(version === "coge_radio") {
        $("#coge-params").find('#pro_or_nu_param').hide(0).html(html).toggle();
        update_gapcost(pro);
    } else {
        $("#ncbi-params").find('#pro_or_nu_param').hide(0).html(html).toggle();
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
        database = $("#database"),
        deferred = $.Deferred(),
        entry;

    if (cache[cache_id]) {
        entry = cache[cache_id];
        deferred.resolve(entry);

        database.html(entry);

        return deferred.promise();
    }

    return $.ajax({
        data: {
            fname: 'database_param',
            program: program,
        },
        success : function(data) {
            cache[cache_id] = data;
            database.html(data);
        },
    });
}

function generate_seq_obj(dsid, dsgid, upstream, downstream, seqview, chr, rc, featid) {
    this.featid = featid;
    this.dsid = dsid;
    this.dsgid = dsgid;
    this.upstream = upstream;
    this.downstream = downstream;
    this.seqview = seqview;
    this.chr = chr;
    this.rc = rc;
    this.gstid = pageObj.gstid;
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
    if (id) {
        var id_array = id.split(',');

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

function popup_blocker_check(windowObject) {
    if (!windowObject) {
        alert("Unable to open a new window. Check your popup blocker settings.");
    }
}

function overlap_checkboxes() {
    var action = $('#overlap_action').val();
    var accn = "";
    $('#hsp_result_table :checkbox').each(function(i) {
        if (this.checked) {
            if (accn)
                accn += ',';
            accn += this.id;
            if (action == 'xls') {
                var t = this.id.split('_');
                accn += '_' + $('#dist' + t[1] + '_' + t[2]).html();
            }
        }
    });
    if (!accn) {
        alert("Please select one or more features.");
        return;
    }

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

        popup_blocker_check(window.open("SeqView.pl?locations=" + locations.join(',')));
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
            popup_blocker_check(window.open(data));
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
            popup_blocker_check(window.open(data));
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
            popup_blocker_check(window.open(data));
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
            popup_blocker_check(window.open(data));
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
            popup_blocker_check(window.open(data));
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
            popup_blocker_check(window.open(data));
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
            popup_blocker_check(window.open(data));
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

            if (obj.error) {
                return alert(obj.error);
            }

            if (obj) {
                if (obj.count > 10) {
                    var remove = obj.count - 10;
                    alert("You have exceeded the number of features you can send to GEvo ( 10 Max ). You currently have "+obj.count+" selected. Please uncheck "+remove+" of your checked item(s).");
                }
                else {
                    popup_blocker_check(window.open(obj.url));
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
            popup_blocker_check(window.open(data));
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
            popup_blocker_check(window.open(data));
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
            popup_blocker_check(window.open(data));
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
            popup_blocker_check(window.open(data));
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
        $('#sequence_popup').dialog('option','buttons',{"View Sequence in SeqView": function() { popup_blocker_check(window.open(seqview)); }});
    }
    $('#sequence_popup').html(html).dialog('open');
}

function get_params(){
    var root = $("#coge-params");
    radio = get_radio("coge_radio","coge");
    var word_size = root.find('#word_size').val();
    var expect = root.find('#e_value').val();
    var match_mismatch = root.find('#match_score').val();
    var matrix = root.find('#matrix').val();
    var program = root.find('#'+radio).val();

    var gapcost;
    if (program == "blastn" || program == 'mega' || program == 'dcmega') {
        var num1 = match_mismatch.substr(0,1);
        var num2 = match_mismatch.substr(2);
        gapcost = root.find('#gapcosts_'+num1+num2).val();
    }
    else {
        gapcost = root.find('#gapcosts_'+matrix).val();
    }

    var filter_query = 0;
    if (root.find('#filter_query')[0].checked) {filter_query=1;}
    var reslimit = root.find('#resultslimit').val();
    //blastz parameters
    var zwordsize = root.find('#blastz_wordsize').val();
    var zgap_start = root.find('#blastz_gap_start').val();
    var zgap_extension = root.find('#blastz_gap_extension').val();
    var zchaining = root.find('#blastz_chaining').val();
    var zthreshold = root.find('#blastz_threshold').val();
    var zmask = root.find('#blastz_mask').val();

    var comp = root.find('#comp_adj').val();
    var seq = $('#seq_box').val();
    return {w : word_size, e : expect,g : gapcost,p : program,mm: match_mismatch,m : matrix,c : comp,s : seq, zw : zwordsize, zgs : zgap_start, zge : zgap_extension, zc : zchaining, zt : zthreshold, zm : zmask, fq : filter_query, rl : reslimit, type: radio };
}

function get_ncbi_params(){
    var gapcost,
        num1,
        num2;

    var root = $("#ncbi-params");
    radio = get_radio("ncbi_radio","ncbi");
    var word_size = root.find('#word_size').val();
    var expect = root.find('#e_value').val();
    var db = root.find('#db').val();
    var match_mismatch = root.find('#match_score').val();
    var matrix = root.find('#matrix').val();

    if (match_mismatch) {
        num1 = match_mismatch.substr(0,1);
        num2 = match_mismatch.substr(2);
        gapcost = root.find('#gapcosts_'+num1+num2).val();
    } else {
        gapcost = root.find('#gapcosts_'+matrix).val();
    }

    var job_title = escape(root.find('#job_title').val());
    var program = root.find('#'+radio).val();

    var comp = root.find('#comp_adj').val();
    var seq = $('#seq_box').val();
    var filter = $("#complexity").val();

    return {w : word_size, e : expect, db : db,g : gapcost,j : job_title,p : program,mm: match_mismatch,m : matrix,c : comp,s : seq, f: filter};
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
    var genomes = $('#genome_choice').getLength(1);

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
        basename:       pageObj.basename,
        color_hsps:     $('#color_by:checked').val(),
        comp:           comp,
        e_value:        expect,
        fid:            pageObj.fid,
        filter_query:   filter_query,
        gapcost:        gapcost ? gapcost.split(' ') : null,
        genomes:        genomes ? genomes.split(',') : null,
        job_title:      job_title,
        match_score:    match_score ? match_score.split(',') : null,
        matrix:         matrix,
        max_results:    resultslimit,
        program:        program,
        query_seq:      seq,
        type:           params.type,
        width:          page_width,
        wordsize:       word_size,
        zchaining:      zchaining,
        zgap_extension: zgap_extension,
        zgap_start:     zgap_start,
        zmask:          zmask,
        zthreshold:     zthreshold,
        zwordsize:      zwordsize
    };

    $.ajax({
        type: "PUT",
        url: 'api/v1/jobs',
        dataType: 'json',
        contentType: "application/json",
        data: JSON.stringify({
            type: 'blast',
            requester: {
                page: PAGE_NAME
            },
            parameters: options
        }),
        success : function(response) {
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
                    .attr("href", response.site_url)
                    .html(response.site_url);

                var link_message = $("<span></span>")
                    .html("Return to this analysis: ")
                    .append(link);

                var logfile = pageObj.tempdir.substring(pageObj.tempdir.indexOf('web') + 4) + '/' + pageObj.basename + '.log';

                status_dialog.find(".dialog-link").html(link_message);
                status_dialog.find(".dialog-log").html($("<a></a>")
                    .attr("href", logfile)
                    .html("Logfile"));
 
                options.fname = "get_results";
                options.logfile = logfile;
                options.gapcost = gapcost;
                options.genomes = genomes;
                options.match_score = match_score;

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
            type: 'POST',
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

    if (typeof setup_button_states !== 'undefined') // mdb added condition 5/10/16 - kludge for EMBED support: this routine is in header.tmpl. Not calling it prevents button highlight on hover.
    	setup_button_states();
    click_all_feat_links(feature_links);
    check_display();
}

//FIXME: separate ncbi and coge parameters
function ncbi_blast(url) {
    var params = get_ncbi_params(),
        pairs,
        coge_pairs,
        options,
        coge_options,
        request;

    if (url == 1) {
        return alert('You have not selected a BLAST program! Please select a program to run.');
    }

    options = {
        "DATABASE": params.db,
        "EXPECT": params.e,
        "QUERY": params.s.replace(/\n/,'%0D'),
        "WORD_SIZE": params.w,
        "GAP_COSTS": params.g,
        "JOB_TITLE": params.j,
        "FILTER": params.f
    };

    var radio = get_radio('ncbi_radio','ncbi');

    //FIXME: CoGe specific options should not be included in ncbi-blast url
    coge_options = {
        program: $("#" + radio).val(),
        expect: params.e,
        database: params.db,
        word_size: params.w,
        gapcost: params.g,
        job: params.j,
        type: radio,
        filter: params.f
    };

    if(seqObj.featid) {
        options["fid"] = seqObj.featid;
    }

    if (program == 'blastn') {
        options["MATCH_SCORES"] = params.mm;
        coge_options["match_score"] = params.mm;
    } else {
        options["MATRIX_NAME"] = params.m;
        coge_options["matrix"] = params.m;

        if (!$('#comp_adj').is(':hidden')) {
            options["COMPOSITION_BASED_STATISTICS"] = params.c;
            coge_options["comp"] = params.c;
        }
    }

    var parameterify = function(value, key) {
        return concat.call(key, "=", value);
    };

    //FIXME: Replace with underscore ie: _.map if library is included
    pairs = map(options, parameterify)
    coge_pairs = map(coge_options, parameterify);

    // Find the selected tab
    var hash = $(".ui-tabs-selected:first a").attr("href");
    var parts = location.href.split(/[#?]/g);

    history.pushState(null, null, concat.call(parts[0], "?", coge_pairs.join("&"), hash));

    request = concat.call(url, "&", pairs.join("&"));
    popup_blocker_check(window.open(request));
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

//FIXME: Replace with a simplified jquery selector
function get_radio(which_type,val){
    if ($('#'+which_type)[0].checked) {
        return val+"_blast_type_n";
    } else {
        return val+"_blast_type_p";
    }
}

function get_seq(which_type) {
    var cache_id,
        dsid = seqObj.dsid,
        dsgid = seqObj.dsgid,
        featid = seqObj.featid,
        chr = seqObj.chr,
        deferred = $.Deferred(),
        program = get_radio(which_type,"coge")
        seqbox = $("#seq_box");

    if (featid) {
        cache_id = concat.call(featid, program, seqObj.upstream,
                               seqObj.downstream, seqObj.rc, seqObj.gstid);

        if (cache[cache_id]) {
            deferred.resolve(cache[cache_id]);
            seqbox.val(cache[cache_id]);
            return deferred.promise();
        }

        seqbox.val('Loading ...');

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
                seqbox.val(html);
            },
        });
    } else if (chr) {
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
    } else if (pageObj.locations) {
        seqbox.val('Loading ...');

        return $.ajax({
            data: {
                fname: 'get_sequence',
                blast_type: program,
                locations: pageObj.locations
            },
            success : function(html) {
                seqbox.val(html);
            },
        });
    } else {
        deferred.resolve("");
        return deferred.promise();
    }
}

function blast_param_on_select(which_type, val) {
    var promise, wordsize = $("#word_size");

    radio = get_radio(which_type, val);
    program = $('#'+radio).val();
    database_param(program);

    // mdb added 11/4/16 -- kludge: default to blastn not working after update to jQuery 3.1.1
    if (!program)
        $('#coge_blast_type_n').val('blastn');

    $('#blast_parameters').hide();
    $('#blastz_parameters').hide();

    if (program == 'lastz') {
        $('#blastz_parameters').toggle();
    } else {
        $('#blast_parameters').toggle();

        if ((program == 'blastx') || (program == 'tblastx')) {
            promise = blast_param("blast_type_p", 1, which_type);
            promise.then(function() {
                wordsize.val(3);
            });
        } else if ((program == 'blastp') || (program == 'tblastn')) {
            promise = blast_param("blast_type_p", 0, which_type);
            promise.then(function() {
                wordsize.val(3);
            })
        } else {
            promise = blast_param('', 0, which_type);
            promise.then(function() {
                if (program == "dcmega") {
                    wordsize.val(11);
                } else {
                    wordsize.val(8);
                }
            })
        }
        (program == 'tblastn' || program == 'tblastx') ? $('#filter_query_row').show() : $('#filter_query_row').hide();
    }

    return promise;
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

//FIXME: remove if underscore library is included
function map(object, func) {
    var key,
        result = [];

    for(key in object) {
        if (object.hasOwnProperty(key)) {
            result.push(func(object[key], key, object));
        }
    }

    return result;
}

//FIXME: remove if underscore library is included
function toObject(pairs) {
    return pairs.reduce(function(a, b) {
        a[b[0]] = b[1];
        return a;
    }, {});
}

function unescapify(value, key) {
    return [key, unescape(value)];
}

function select_by_value($elements, property, value) {
    return $elements.filter(function() {
        return this.value === String(value);
    }).prop(property, true);
}

var TypeSelectorMixin = {
    _select_type: function ($elements) {
        select_by_value($elements, 'checked', this.params['type']);
    },

    _select_program: function () {
        var program = this.params['program'];
        this.root.find("#" + this.params['type']).val(program);
    }
};

var ProteinMixin = {
    _select_composition: function () {
        this.root.find("#comp_adj").val(this.params['composition']);
    }
};

var ScoringMixin = {
    // This belongs in NucleotideMixin
    _select_match_score: function () {
        var elements = this.root.find('#match_score option');
        select_by_value(elements, 'selected', this.params['match_score']);
    },

    // This belongs in ProteinMixin
    _select_matrix_score: function () {
        var elements = this.root.find('#matrix option');
        select_by_value(elements, 'selected', this.params['matrix_score']);
    },

    _select_evalue: function () {
        var elements = this.root.find('#e_value option');
        select_by_value(elements, 'selected', this.params['evalue']);
    },

    // Tightly coupled to matrix/match scoring for picking the gapcost select
    _select_gapcost: function($element) {
        if (!$element || !$element.val()) {
            console.warn('blast:_select_gapcost: null input');
            return;
        }

        var matchPattern = /[,]/g;
        var val = $element.val().replace(matchPattern, "");

        //FIXME: This should really only be one gap cost element
        this.root.find('.gapcosts').hide()
        var gapcost = this.root.find('#gapcosts_' + val).toggle();

        // Requires a space between characters
        if (this.params['gapcost']) {
            cost = this.params['gapcost'].split(/[\s+,]/).join(" ");
            gapcost.val(cost);
        }
    },

    _select_word_size: function () {
        this.root.find("#word_size").val(this.params['wordsize']);
    }
};

function Ncbi(selector, params) {
    this.params = params || {};

    this.defaults = {
        type: 'coge_blast_type_n',
        match_score: null,
        matrix_score: null,
        evalue: 1e-3,
        wordsize: 8,
        limit: 100,
        gapcost: null,
        filter: 0,
        composition: 1,
        database: null,
        program: null
    };

    //FIXME: Replace with underscore ie: _.map if library is included
    this.params = toObject(map($.extend(this.defaults, this.params), unescapify));
    this.root = $(selector);
}

$.extend(Ncbi.prototype, TypeSelectorMixin, ScoringMixin, ProteinMixin, {
    _select_database: function () {
        var elements = this.root.find("#db option");
        select_by_value(elements, 'selected', this.params['database']);
    },

    _select_filter: function () {
        var elements = this.root.find("#complexity option");
        select_by_value(elements, 'selected', this.params['filter']);
    },

    _select_job_title: function () {
        if (this.params['job']) {
            this.root.find("#job_title").val(this.params['job']);
        }
    },

    update_nucleotide: function () {
        // Set the match score to be used
        this._select_match_score();

        // Set the e-value parameter
        this._select_evalue();

        // Set the word size parameter
        this._select_word_size();

        // Select the database to be searched
        this._select_database();

        // Set the filter to be used
        this._select_filter();

        // Selects the gap cost select and option
        this._select_gapcost(this.root.find("#match_score"));
    },

    update_default: function () {
        // Set the matrix score to be used
        this._select_matrix_score();

        // Set the e-value parameter
        this._select_evalue();

        // Set the word size parameter
        this._select_word_size();

        // Select the database to be searched
        this._select_database();

        // Set the filter to be used
        this._select_filter();

        // Selects the gap cost select and option
        this._select_gapcost(this.root.find("#matrix"));

        // Select the composition adjustments
        this._select_composition();
    },

    update_display: function () {
        var self = this;

        // Select the blast type (nucleotide vs protein)
        var elements = this.root.find('input[name="ncbiblast"]');
        this._select_type(elements);

        // Set the blast tool being used (depends on type)
        this._select_program();


        // Set the title of the job
        this._select_job_title();

        // dispatch fetch the blast parameters'
        var promise = blast_param_on_select('ncbi_radio', 'ncbi');

        // Set the options after the parameters have been returned
        if (promise)
            promise.always(function() {
                switch (self.params['program']) {
                    case 'blastn': self.update_nucleotide(); break;
                    default: self.update_default(); break;
                }
            })
        }
});

function Blast(selector, params) {
    this.params = params || {};

    this.defaults = {
        type: 'coge_blast_type_n',
        match_score: null,
        matrix_score: null,
        evalue: 1e-5,
        wordsize: 8,
        limit: 20,
        gapcost: "5 2",
        filtered: 1,
        composition: 1,
        blastz_wordsize: 8,
        blastz_gap_start: 400,
        blastz_gap_extension: 30,
        blastz_chaining: 0,
        blastz_threshold: 3000,
        blastz_mask: 0,
        program: null
    };

    //FIXME: Replace with underscore ie: _.map if library is included
    this.params = toObject(map($.extend(this.defaults, this.params), unescapify));
    this.root = $(selector);
};

$.extend(Blast.prototype, TypeSelectorMixin, ScoringMixin, ProteinMixin, {
    _select_color_by: function () {
        var elements = $('input[name="color_by"]');
        select_by_value(elements, 'checked', this.params['color']);
    },

    _select_limit: function () {
        $('#resultslimit').val(this.params['limit']);
    },

    _select_filtered: function () {
        var elements = $('input[name="filter_query"]');
        select_by_value(elements, 'checked', this.params['filtered'])
    },

    _select_blastz_options: function () {
        $('#blastz_wordsize').val(this.params['blastz_wordsize']);
        $('#blastz_gap_start').val(this.params['blastz_gap_start']);
        $('#blastz_gap_extension').val(this.params['blastz_gap_extension']);
        $('#blastz_threshold').val(this.params['blastz_threshold']);
        $('#blastz_mask').val(this.params['blastz_mask']);

        var elements =$('#blastz_chaining options');
        select_by_value(elements, 'selected', this.params['blastz_chaining'])
    },

    update_default: function () {
        // Set the match score to be used
        this._select_match_score();

        // Set the e-value parameter
        this._select_evalue();

        // Set the word size parameter
        this._select_word_size();

        // Set the result limits
        this._select_limit();

        // Set whether the query sequence will be filtered
        this._select_filtered();

        // Selects the gap cost select and option
        this._select_gapcost(this.root.find("#match_score"));
    },

    update_blastz: function () {
        // Set blastz specific options
        this._select_blastz_options();
    },

    update_protein: function () {
        // Set the matrix score to be used
        this._select_matrix_score();

        // Set the e-value parameter
        this._select_evalue();

        // Set the word size parameter
        this._select_word_size();

        // Set the result limits
        this._select_limit();

        // Set whether the query sequence will be filtered
        this._select_filtered();

        // Selects the gap cost select and option
        this._select_gapcost(this.root.find("#matrix"));

        // Select the composition adjustments
        this._select_composition();
    },

    update_display: function () {
        var self = this;

        // Select the blast type (nucleotide vs protein)
        var elements = this.root.find('input[name="cogeblast"]');
        this._select_type(elements);

        // Set the blast tool being used (depends on type)
        this._select_program();

        // Select the blast hit coloring scheme
        this._select_color_by();

        // dispatch fetch the blast parameters'
        var promise = blast_param_on_select('coge_radio', 'coge');

        // Set the options after the parameters have been returned
        if (promise)
            promise.always(function() {
                switch (self.params['program']) {
                    case 'lastz': self.update_blastz(); break;
                    case 'tblastx': self.update_protein(); break;
                    case 'tblastn': self.update_protein(); break;
                    default: self.update_default(); break;
                }
            })
        }
});

function getParamsFromUrl() {
    var query = location.search.substr(1),
        data = query.split(/[&;]/),
        params = {},
        pair, i;

    for(i = 0; i < data.length; i++) {
        pair = data[i].split("=");
        params[pair[0]] = pair[1];
    }

    return params;
}

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
    var matrix = $("#coge-params").find('#matrix').val();
    popup_blocker_check(window.open('MatrixView.pl?matrix='+matrix));
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

var initialized = false;
function select_tab(event, ui) {
    var button = $("#run_blast");
    var index = (ui.newTab ? ui.newTab.index() : 0);
    switch(index) {
        case 0:
            adjust_blast_types();
            button.unbind().click(run_coge_blast).html("Run CoGe BLAST");

            if (!initialized) {
                var queryParams = getParamsFromUrl();

                var analysis = new Blast("#coge-params", {
                    type: queryParams['type'],
                    color: queryParams['color_hsps'],
                    match_score: queryParams['match_score'],
                    matrix_score: queryParams['matrix'],
                    evalue: queryParams['expect'],
                    wordsize: queryParams['wordsize'],
                    limit: queryParams['resultslimit'],
                    gapcost: queryParams['gapcost'],
                    filtered: queryParams['filter_query'],
                    composition: queryParams['comp'],
                    blastz_wordsize: queryParams['zwordsize'],
                    blastz_gap_start: queryParams['zgap_start'],
                    blastz_gap_extension: queryParams['zgap_exten'],
                    blastz_chaining: queryParams['zchaining'],
                    blastz_threshold: queryParams['zthreshold'],
                    blastz_mask: queryParams['zmask'],
                    program: queryParams['program']
                });

                analysis.update_display();
                initialized = true;
            }
            break;
        case 1:
            adjust_blast_types(1);
            blast_param_on_select('ncbi_radio','ncbi');
            button.unbind().click(function() {
                select_blast();
                ga('send', 'event', 'cogeblast', 'run', 'ncbi');
            }).html("Run NCBI BLAST");

            if (!initialized) {
                var queryParams = getParamsFromUrl();

                var analysis = new Ncbi("#ncbi-params", {
                    type: queryParams['type'],
                    match_score: queryParams['match_score'],
                    matrix_score: queryParams['matrix'],
                    evalue: queryParams['expect'],
                    wordsize: queryParams['word_size'],
                    limit: queryParams['resultslimit'],
                    gapcost: queryParams['gapcost'],
                    filter: queryParams['filter'],
                    composition: queryParams['comp'],
                    database: queryParams['database'],
                    job: queryParams['job'],
                    program: queryParams['program']
                });

                analysis.update_display();
                initialized = true;
            }
            break;
    }
}
