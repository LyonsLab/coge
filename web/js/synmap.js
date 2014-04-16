function init_table() {
    $("#prev_table").tablesorter({
        sortClassAsc: 'headerSortUp',       // Class name for ascending sorting action to header
        sortClassDesc: 'headerSortDown',    // Class name for descending sorting action to header
        headerClass: 'header',      // Class name for headers (th's)
        widgets: ['zebra'],
        textExtraction: 'complex'
    });
 }

function has_organisms () {
    return ($('#org_id1').val() != "") &&
        ($('#org_id2').val() != "");
}

$(function() {$("#pair_info").draggable();});
$(function() {$("#tabs").tabs({selected:0});});
$(function() {$(".resizable").resizable();});

function get_gc (dsgid, divid)
{
    $('#'+divid).removeClass('link').html('loading...');
    get_dsg_gc(['args__dsgid','args__'+dsgid,'args__text','args__1'],[divid]);
}

function get_organism_chain(type,val,i)
{
    $('#org_list').html('<input type=hidden id = "org_id"+i><font class="loading"></font>');
    if (type == 'name') {get_orgs(['args__name','args__'+val,'args__i','args__'+i], ['org_list'+i]);}
    else if (type == 'desc') {get_orgs(['args__desc','args__'+val,'args__i','args__'+i], ['org_list'+i]);}
    $('#dsg_info'+i).html('<div class="loading dna_small small">loading. . .</div>');
    ajax_wait("gen_dsg_menu(['args__oid','org_id"+i+"', 'args__num','args__"+i+"'],['dsg_menu"+i+"', 'genome_message"+i+"']);");
    ajax_wait("check_previous_analyses();");
    ajax_wait("get_genome_info(['args__dsgid','dsgid"+i+"','args__org_num','args__"+i+"'],[handle_dsg_info]);");
}

function get_genome_info_chain(i) {
    $('#dsg_info'+i).html('<div class=dna_small class=loading class=small>loading. . .</div>');
    // ajax_wait("gen_dsg_menu(['args__oid','org_id"+i+"', 'args__num','args__"+i+"'],['dsg_menu"+i+"','genome_message"+i+"']);");
    gen_dsg_menu(['args__oid','org_id'+i, 'args__num','args__'+i],['dsg_menu'+i, 'genome_message'+i]);
    $('#depth_org_1').html($('#org_id1 option:selected').html());
    $('#depth_org_2').html($('#org_id2 option:selected').html());

    ajax_wait("check_previous_analyses();");
    ajax_wait("get_genome_info(['args__dsgid','dsgid"+i+"','args__org_num','args__"+i+"'],[handle_dsg_info]);");
}

function rand () {
    return ( Math.floor ( Math.random ( ) * 99999999 + 1 ) );
}

function populate_page_obj(basefile) {
    if (!basefile) {
        basefile = "SynMap_"+rand();
    }
    pageObj.basename = basefile;
    pageObj.nolog = 0;
    pageObj.waittime = 1000;
    pageObj.runtime = 0;
    pageObj.fetch_error = 0;
    pageObj.error = 0;
    pageObj.engine = "<span class=\"alert\">The job engine has failed.</span><br>Please use the link below to use the previous version of SynMap.";
}

function run_synmap(scheduled, regenerate){
    populate_page_obj();

    var org_name1 = pageObj.org_name1;
    var org_name2 = pageObj.org_name2;
    //feat_type 1 == CDS, 2 == genomic
    var feat_type1 = $('#feat_type1').val();
    var feat_type2 = $('#feat_type2').val();
    var org_length1 = pageObj.org_length1;
    var org_length2 = pageObj.org_length2;
    //seq_type == 1 is unmasked
    var seq_type1 = pageObj.seq_type1;
    var seq_type2 = pageObj.seq_type2;
    //check to see if we will allow this run
    var max_size = 100000000;
    // console.log (org_name1, org_length1, seq_type1);
    // console.log (org_name2, org_length2, seq_type2);
    if (( org_length1 > max_size && feat_type1 == 2 && seq_type1 == 1) &&
        ( org_length2 > max_size && feat_type2 == 2 && seq_type2 == 1) ) {
        var message = "You are trying to compare unmasked genomic sequences that are large!  This is a bad idea.  Chances are there will be many repeat sequences that will cause the entire pipeline to take a long time to complete.  This usually means that the analyses will use a lot of RAM and other resources.  As such, these jobs are usually killed before they can complete.  Please contact coge.genome@gmail.com for assistance with your analysis.";
            alert(message);
            $('#log_text').hide(0);
                $('#results').show(0).html("<span class=alert>Analysis Blocked:</span>"+message);
                return;
    }


    if (!has_organisms())
        return;

    if ($('#blast').val() == 5 && (feat_type1 != 1 || feat_type2 != 1) ) {
        alert('BlastP only works if both genomes have protein coding sequences (CDS) AND CDS is selected for both!');
        return;
    }

    var overlay = $("#overlay").show();

    //TODO Scale polling time linearly with long running jobs
    var duration = pageObj.waittime;

    $('#results').hide();

    if (regenerate) {
        return schedule(get_params("go", regenerate));
    }

    return $.ajax({
        type: 'GET',
        data: get_params("get_results"),
        dataType: "json",
        success: function(data) {
            if (!data.error) {
                overlay.hide();
                $("#synmap_zoom_box").draggable();
                $('#results').html(data.html).slideDown();
            } else {
                schedule(get_params("go"))
            }
        },
        error: schedule.bind(this, get_params("go"))
    });
}

function schedule(params) {
    pageObj.nolog=1;
    var status_dialog = $('#synmap_dialog');

    return $.ajax({
        data: params,
        dataType: 'json',
        success: function(data) {
            $("#overlay").hide();
            status_dialog.dialog('open');
            if (data.status == 'Attached' || data.status == 'Scheduled') {
                var link = "Return to this analysis: <a href="
                + data.link + " onclick=window.open('tiny')"
                + "target = _new>" + data.link + "</a><br>"
                + "To run this analysis on the previous version of SynMap <a href="
                + data.old_link + " onclick=window.open('tiny')"
                + "target = _new>click here</a>";

                var logfile = '<a href="tmp/SynMap/'
                + pageObj.basename + '.log">Logfile</a>';

                $('#dialog_log').html(logfile);
                $('#synmap_link').html(link);

                update_dialog(data.request, "#synmap_dialog", synmap_formatter,
                        get_params("get_results"));
            } else {
                status_dialog.find('#text').html(pageObj.engine);
                status_dialog.find('#progress').hide();
                status_dialog.find('#dialog_error').slideDown();
            }
        },
        error: function(err) {
            $("#overlay").hide();
            status_dialog.dialog('open');
            status_dialog.find('#progress').hide();
            status_dialog.find('#dialog_error').slideDown();
        }
    })
}

function get_params(name, regenerate) {
    return {
        fname: name,
        tdd: $('#tdd').val(),
        D: $('#D').val(),
        A: $('#A').val(),
        beta: pageObj.beta,
        gm: $('#gm').val(),
        Dm: $('#Dm').val(),
        blast: $('#blast').val(),
        feat_type1: $('#feat_type1').val(),
        feat_type2: $('#feat_type2').val(),
        dsgid1: $('#dsgid1').val(),
        dsgid2: $('#dsgid2').val(),
        jobtitle: $('#jobtitle').val(),
        basename: pageObj.basename,
        email: $('#email').val(),
        regen_images: regenerate,
        width: $('#master_width').val(),
        dagchainer_type: $('#dagchainer_type').filter(':checked').val(),
        ks_type: $('#ks_type').val(),
        assemble: $('#assemble')[0].checked,
        axis_metric: $('#axis_metric').val(),
        axis_relationship: $('#axis_relationship').val(),
        min_chr_size: $('#min_chr_size').val(),
        spa_ref_genome: $('#spa_ref_genome').val(),
        show_non_syn: $('#show_non_syn')[0].checked,
        color_type: $('#color_type').val(),
        box_diags: $('#box_diags')[0].checked,
        merge_algo: $('#merge_algo').val(),
        depth_algo: $('#depth_algo').val(),
        depth_org_1_ratio: $('#depth_org_1_ratio').val(),
        depth_org_2_ratio: $('#depth_org_2_ratio').val(),
        depth_overlap: $('#depth_overlap').val(),
        fid1: pageObj.fid1,
        fid2: pageObj.fid2,
        show_non_syn_dots: $('#show_non_syn_dots')[0].checked,
        flip: $('#flip')[0].checked,
        clabel: $('#clabel')[0].checked,
        skip_rand: $('#skiprand')[0].checked,
        color_scheme: $('#color_scheme').val(),
        chr_sort_order: $('#chr_sort_order').val(),
        codeml_min: $('#codeml_min').val(),
        codeml_max: $('#codeml_max').val(),
        logks: $('#logks')[0].checked,
        csco: $('#csco').val(),
        jquery_ajax: 1,
    };
}

function close_dialog() {
    var dialog_window = $('#synmap_dialog');
    if(dialog_window.dialog('isOpen')) {
        dialog_window.dialog('close');
    }

    dialog_window.find('#text').html('');
    dialog_window.find('#progress').show();
    dialog_window.find('#dialog_error').hide();
    dialog_window.find('#dialog_success').hide();
}

function load_results() {
    $('#intro').hide();
    $('#log_text').hide();
    $('#results').fadeIn();
}

function handle_results(val){
    $('#results').html(val);
    $(function() {$("#synmap_zoom_box").draggable();});
    setup_button_states();
    ajax_wait("check_previous_analyses();");
}

function check_previous_analyses(){
    var gid1 = $('#org_id1').val();
    var gid2 = $('#org_id2').val();

    if (gid1 && gid2) {
        $.ajax({
            data: {
                jquery_ajax: 1,
                oid1: gid1,
                oid2: gid2,
                fname: 'get_previous_analyses',
            },
            success: function(data)  {
                load_previous_analyses(data);
            }
        });
    }
//get_previous_analyses(['args__oid1','org_id1', 'args__oid2','org_id2'],[load_previous_analyses]);
}

function load_previous_analyses (stuff) {
    $('#previous_analyses').html(stuff);
    init_table();
}

function update_params(val) {
    var cmd;
    var params;
    var type;

    if (val) {
        params = val.split('_');
    } else {
        params = $('#prev_params').val()[0].split('_');;
    }

    if ($('#org_id1').val() == params[3]) {
        cmd = "$('#feat_type1').attr('value', '"+params[5]+"');$('#feat_type2').attr('value', '"+params[8]+"');";
        $('#dsgid1').attr('value',params[4]);
        $('#dsgid2').attr('value',params[7]);
    } else {
        cmd = "$('#feat_type2').attr('value', '"+params[5]+"');$('#feat_type1').attr('value', '"+params[8]+"');";
        $('#dsgid2').attr('value',params[4]);
        $('#dsgid1').attr('value',params[7]);
    }

    get_genome_info(['args__dsgid','dsgid1','args__org_num','args__1'],[handle_dsg_info]);
    get_genome_info(['args__dsgid','dsgid2','args__org_num','args__2'],[handle_dsg_info]);
    ajax_wait(cmd);

    $('#blast').attr('value',params[9]);

    if (params[10] == 'Distance') {
        $("input[name='dagchainer_type']:nth(1)").attr("checked","checked");
        type=" bp";
    } else {
        $("input[name='dagchainer_type']:nth(0)").attr("checked","checked");
        type= " genes";
    }

    display_dagchainer_settings([params[1],params[2]],type);
    $('#c').val(params[11]);
    merge_select_check();
    depth_algo_check();
}

function handle_dsg_info(dsg_html, feat_menu, genome_message, length, org_num, org_name, seq_id) {
    $('#dsg_info'+org_num).html(dsg_html);
    $('#feattype_menu'+org_num).html(feat_menu);
    $('#genome_message'+org_num).html(genome_message);

    if (org_num == 1) {
        pageObj.org_length1 = length;
        pageObj.org_name1 = org_name;
        pageObj.seq_type1 = seq_id;
    } else {
        pageObj.org_length2 = length;
        pageObj.org_name2 = org_name;
        pageObj.seq_type2 = seq_id;
    }
}


function set_dagchainer_defaults(params, type) {
    var settings = $('#dagchainer_default').val();

    if (!(params && type)) {
        if ($('#dagchainer_type')[0].checked) {
            params = [20,5,0,0];
            type = " genes";
        } else {
            if (settings == 1) { // for plant genomes
                params = [120000, 5, 96000,480000];
            } else if (settings == 2) { // for microbe genomes
                params = [2000, 5, 4000, 8000];
            }

            type=" bp";
        }
    }

    if (!params) {
        return;
    }

    $('#D').val(params[0]);
    $('#A').val(params[1]);

    if (typeof(params[2]) == 'undefined') {
        params[2] = 4*params[0];
    }

    if (typeof(params[3]) == 'undefined') {
        params[3] = 4*params[1];
    }

    $('#gm').val(params[2]);
    $('#Dm').val(params[3]);
    $('.distance_type').html(type);
}

function ajax_wait (val){
    if (ajax.length) {
        setTimeout("ajax_wait("+'"'+val+'"'+")",100);
        return;
    }

    eval(val);
}

function timing(val, val2){
    var searchterm;
    namere = /name/;
    descre = /desc/;

    if (namere.exec(val)) {
        searchterm = $('#'+val).val();
    } else if (descre.exec(val)) {
        searchterm = $('#'+val).val();
    }

    if (!searchterm) {
        val=0;
    }

    if(searchterm == "Search") {
        searchterm = "";
    }

    pageobjsearch = "search"+val;
    pageobjtime = "time"+val;

    if (pageObj.pageobjsearch && pageObj.pageobjsearch == searchterm+val) {
    //    return;
    }

    pageObj.pageobjsearch=searchterm+val;

    if (pageObj.pageobjtime){
        clearTimeout(pageObj.pageobjtime);
    }

    re = /(\d+)/;
    i = re.exec(val);

    if (namere.exec(val)) {
        if (val2) {
            get_organism_chain('name',$('#'+val).val(),i[0])
        } else {
            pageObj.pageobjtime = setTimeout("get_organism_chain('name',$('#"+val+"').val(),i[0])",500);
        }
    } else if (descre.exec(val)) {
        if (val2) {
            get_organism_chain('desc',$('#'+val).val(),i[0])
        } else {
            pageObj.pageobjtime = setTimeout("get_organism_chain('desc',$('#"+val+"').val(),i[0])",200);
        }
    }
}


function display_dagchainer_settings(params,type) {

    if ($('#dagchainer_type')[0].checked) {
        $('#dagchainer_distance').hide(0);
    } else {
        $('#dagchainer_distance').show(0);
    }

    set_dagchainer_defaults(params, type);
}

//These two functions are involved with emailing results
function toggle_email_field() {
    if ($('#check_email')[0].checked) {
        $('.email_box').show(0);
    } else {
        $('.email_box').hide(0);
        $('#email').val('');
        $('#email_error').hide(0);
    }
}

function address_validity_check(validity) {
    if (validity) {
        if(validity == 'invalid') {
            $('#email_error').show(0);
        } else {
            $('#email_error').hide(0);
        }
    } else {
        check_address_validity(['email'],[address_validity_check]);
    }
}

function fill_jobtitle(){
    var title;
    var org1 = $('#org_id1 option:selected').html() || 0;
    var org2 = $('#org_id2 option:selected').html() || 0;

    if (org1 != 0) {
        org1 = org1.replace(/\s+\(id\d+\)$/,"");
    }

    if (org2 != 0) {
        org2 = org2.replace(/\s+\(id\d+\)$/,"");
    }

    if (org1 != 0 && org2 != 0) {
        title = org1 + " v. " + org2;
    } else if (org1 != 0) {
        title = org1;
    } else if (org2 != 0) {
        title = org2;
    } else {
        return;
    }

    $('#jobtitle').val(title);
}

function update_basename(basename){
    pageObj.basename=basename;
}

function reset_basename(){
    if(pageObj.basename) pageObj.basename=0;
}

function synmap_formatter(item) {
    var msg;
    var row = $('<li>'+ item.description + ' </li>');
    row.addClass('small');

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


function update_dialog(request, identifier, formatter, args) {
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

    var fetch_results = function(completed) {
        var request = window.location.href.split('?')[0];
        args.fname = 'get_results';
        dialog = $(identifier);

        $.ajax({
            type: 'GET',
            url: request,
            data: args,
            dataType: "json",
            success: function(data) {
                if (completed && !data.error) {
                    $("#synmap_zoom_box").draggable();
                    $('#results').html(data.html);
                    dialog.find('#progress').hide();
                    dialog.find('#dialog_success').slideDown();
                } else {
                    $('#results').html(data.error);
                    dialog.find('#progress').hide();
                    dialog.find('#dialog_error').slideDown();
                }
            },
            error: function(data) {
                if (pageObj.fetch_error >= 3) {
                    dialog.find('#progress').hide();
                    dialog.find('#dialog_error').slideDown();
                } else {
                    pageObj.fetch_error += 1;
                    console.log("error");
                    var callback = function() {fetch_results(completed)};
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
            update_dialog(request, identifier, formatter, args);
        }

        if (json.error) {
            pageObj.error++;
            if (pageObj.error > 3) {
                workflow_status.html(pageObj.engine);
                dialog.find('#text').html(workflow_status);
                dialog.find('#progress').hide();
                dialog.find('#dialog_error').slideDown();
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
            fetch_results(true);
        } else if (current_status == "failed" || current_status == "error"
                || current_status == "terminated"
                || current_status == "cancelled") {
            workflow_status.find('span').addClass('alert');
            fetch_results(false);
        } else if (current_status == "notfound") {
            setTimeout(callback, timeout);
            return;
        } else {
            workflow_status.find('span') .addClass('running');
            setTimeout(callback, timeout);
        }

        results.push(workflow_status);
        data.append(results);
        dialog.find('#text').html(data);
    };

    get_status();
}

function synteny_zoom(dsgid1, dsgid2, basename, chr1, chr2, ksdb) {
    var url = 'dsg1='+dsgid1+';dsg2='+dsgid2+';chr1='+chr1+';chr2='+chr2+';base='+basename;
    var loc = $('#map_loc').val();
    var width = $('#zoom_width').val();
    var min = $('#zoom_min').val();
    var max = $('#zoom_max').val();
    var am = $('#axis_metric').val();
    var fid1=0;
    if (pageObj.fid1) {fid1 = pageObj.fid1;}
    var fid2=0;
    if (pageObj.fid2) {fid2 = pageObj.fid2;}
    var ct = $('#color_type').val();
    var loc = pageObj.loc;

    if (!loc) {loc=1;}

    loc++;
    pageObj.loc=loc;

    win = window.open ('DisplayMessage.pl', 'win'+loc,'width=400,height=200,scrollbars=1');
    win.focus();

    get_dotplot(
        ['args__url','args__'+url, 'args__loc','args__'+loc, 'args__flip','args__'+$('#flip')[0].checked,'args__regen_images','args__'+$('#regen_images')[0].checked, 'args__width', 'args__'+width, 'args__ksdb','args__'+ksdb,'args__kstype','ks_type','args__min', 'args__'+min,'args__max', 'args__'+max, 'args__am', 'args__'+am, 'args__ct','args__'+ct, 'args__bd', 'args__'+$('#box_diags')[0].checked, 'args__color_scheme','color_scheme', 'args__am','axis_metric', 'args__ar','axis_relationship', 'args__fid1','args__'+fid1, 'args__fid2', 'args__'+fid2],[open_window]);
}

function open_window (url, loc, width, height) {
    if (!loc) {
        loc = pageObj.loc;
    }

    if (!loc) {
        loc=1;
    }
    my_window = window.open(url,"win"+loc,'"width='+width+',height='+height+', scrollbars=1"');
    my_window.resizeTo(width,height);
}

function merge_select_check () {
    var merge_algo = $('#merge_algo').val();

    if (merge_algo == 0) {
        $('#merge_algo_options').hide();
    } else if (merge_algo == 1) {
        $('#merge_algo_options').show();
        $('#max_dist_merge').hide();
    } else {
        $('#merge_algo_options').show();
        $('#max_dist_merge').show();
    }
}

function depth_algo_check() {
   var depth_algo = $('#depth_algo').val();

    if (depth_algo == 0) {
        $('#depth_options').hide();
    } else if (depth_algo == 1) {
        $('#depth_options').show();
    }
}

function post_to_grimm(seq1, seq2) {
    var url = "http://nbcr.sdsc.edu/GRIMM/grimm.cgi#report";
    var query_form = document.createElement("form");
    var input1 = document.createElement("textarea");
    var input2 = document.createElement("textarea");

    seq1 = seq1.replace(/\|\|/g,"\n");
    seq2 = seq2.replace(/\|\|/g,"\n");

    input1.name="genome1";
    input1.value=seq1;

    input2.name="genome2";
    input2.value=seq2;

    query_form.method="post" ;
    query_form.action=url;
    query_form.setAttribute("target", "_blank");
    query_form.setAttribute("name", "genomeForm");
    query_form.appendChild(input1);
    query_form.appendChild(input2);
    query_form.submit("action");
}

function add(a, b) {
    return a + b;
};

function sum(iterable) {
    return iterable.reduce(add);
}

function pick(attribute) {
    return function(object, index, list) {
        return object[attribute];
    };
};

function zip(coll1, coll2) {
    var length = Math.min(coll1.length, coll2.length),
        index,
        zipped = [];

    for(index = 0; index < length; index++) {
        zipped.push([coll1[index], coll2[index]]);
    }

    return zipped;
}

function scan(func, initial, iterable) {
    var values = [initial],
        index;

    for(index = 0; index < iterable.length; index++) {
        values.push(func(values[index], iterable[index]));
    }

    return values;
}

function toObject(pairs) {
    return pairs.reduce(function(object, value, index) {
        object[value[0]] = value[1];
        return object
    }, {});
}

function compareAlphaNumeric(a, b) {
    var regAlpha = /[^a-zA-Z_]/g,
        regNumeric = /[^0-9]/g,
        aString = a.toString(),
        bString = b.toString(),
        aAlpha,
        aNum,
        bAlpha,
        bNum;

    aAlpha = aString.replace(regAlpha, "");
    bAlpha = bString.replace(regAlpha, "");

    if(aAlpha === bAlpha) {
        aNum = parseInt(aString.replace(regNumeric, ""), 10);
        bNum = parseInt(bString.replace(regNumeric, ""), 10);

        return aNum === bNum ? 0 : aNum > bNum ? 1 : -1;
    }

    return aAlpha.localeCompare(bAlpha);
}

function sortBy(attribute, cmp) {
    return function (a, b) {
        if (cmp) {
            return cmp(a[attribute], b[attribute]);
        } else {
            return a[attribute] - b[attribute];
        }
    };
}

function identity(x) {
    return x;
}

function inverse(func) {
    return function(a, b) {
        return func(b, a);
    }
}

function getValues(object, keys) {
    var index,
        coll = [];

    for(index = 0; index < keys.length; index++) {
        coll.push(object[keys[index]]);
    }

    return coll;
}

function checkRequestSize(url) {
    var size,
        request;

    request = $.ajax({
        type: "HEAD",
        async: false,
        url: url
    });

    if (request.status === 200) {
        return +request.getResponseHeader("Content-Length");
    }
}

var synmap = function(element, metric, sort) {
    var genomes,
        pairs = [],
        plots = [];

    my = {};

    my.loadPlots = function(json) {
        var keys,
            type,
            layer,
            source,
            reference;

        genomes = json.genomes;

        for(layer in json.layers) {
            type = Object.keys(json.layers[layer].data)[0];
            for (reference in json.layers[layer].data[type]) {
                for (source in json.layers[layer].data[type][reference]) {
                    pairs.push([reference, source]);
                }
            }
        }

        // generate plots
        keys = pairs.pop();
        data = generatePlot(keys[0], keys[1], json.layers, metric, sort);

        var handler =  function() {
            var active = [],
                ids = Object.keys(data.layers).reverse();

            for (var index = 0; index < ids.length; index++) {
                if (data.layers[ids[index]].enabled) {
                    active.push(data.layers[ids[index]]);
                }
            }

            return active;
        };

        var plot = new DotPlot(element, {
            size: { width: 1000, height: 800 },
            genomes: [{
                name: data.xtitle,
                length: data.xtotal,
                chromosomes: data.xlabels,
                fetchDataHandler: handler
            },
            {
                name: data.ytitle,
                length: data.ytotal,
                chromosomes: data.ylabels,
                fetchDataHandler: handler
            }],
            fetchDataHandler: handler,
            style: {
                position: "relative",
                minHeight: "800px"
            }
        });
        plot.redraw();

        plots.length = 0;
        plots.push({
            "data": data,
            "plot": plot
        });
    }

    my.toggleLayer = function(layerId, enabled) {
        var index;

        for(index = 0; index < plots.length; index++) {
            if (plots[index].data.layers[layerId]) {
                plots[index].data.layers[layerId].enabled = enabled;
                plots[index].plot.redraw();
            }
        }
    }

    my.reset = function() {
        var index;

        for(index = 0; index < plots.length; index++) {
            plots[index].plot.reset();
        }
    }

    function generatePlot(genome1, genome2, layers, metric, by) {
        var xlabels,
            ylabels,
            xtotal,
            ytotal,
            layers;

        // Construct labels with property axis and offsets
        // x is flipped to sort by maximum value
        // y is not flipped because increasing value is down in canvas
        xlabels = generateLabels(genomes[genome1].chromosomes, metric, inverse(by));
        ylabels = generateLabels(genomes[genome2].chromosomes, metric, by);

        // Generate the layers
        layers = createAllLayers(layers, genome1, genome2, xlabels, ylabels);

        return {
            "xtitle": genomes[genome1].name,
            "ytitle" : genomes[genome2].name,
            "xlabels": xlabels,
            "ylabels": ylabels,
            "xtotal": sum(xlabels.map(pick("length"))),
            "ytotal": sum(ylabels.map(pick("length"))),
            "layers": layers
        };
    }

    function generateLabels(chromosomes, metric, by) {
        return orderBy(getLengths(chromosomes, metric), by);
    }

    function getLengths(chromosomes, prop) {
        return chromosomes.map(function(x) {
            return { name : x.name, length: x[prop] };
        });
    }

    function orderBy(chromosomes, by) {
        return chromosomes.sort(by);
    };

    function offsets(chromosomes) {
        return scan(add, 0, chromosomes.map(pick("length")));
    };

    function generateIndexBy(collection, by) {
        var index = {};

        for(i = 0; i < collection.length; i++) {
            index[by(collection[i])] = i;
        };

        return index;
    };

    function createAllLayers(layers, genome1, genome2, xlabels, ylabels) {
        var xindex = generateIndexBy(xlabels, pick("name")),
            yindex = generateIndexBy(ylabels, pick("name")),
            xoffsets = offsets(xlabels),
            // For Canvas coordinates the chromosome should be offset
            // starting from the end of its length
            yoffsets = offsets(ylabels).slice(1),
            output = {},
            create,
            layerId,
            rawLayers,
            transformedLayers;

        create = createLayer.bind(null, genome1, genome2, xoffsets, yoffsets,
                    xindex, yindex);

        layerIds = Object.keys(layers);
        rawLayers = getValues(layers, layerIds);
        transformedLayers = rawLayers.map(create);

        return toObject(zip(layerIds, transformedLayers));
    }

    function createLayer(genome1, genome2, xoffset, yoffset, xindex, yindex, layer) {
        var style = layer.style || {},
            enabled = layer.enabled || false,
            points,
            rects,
            lines,
            raw;

        if (layer.data.points) {
            raw = layer.data.points[genome1][genome2];
            points = transformBy(raw, xoffset, yoffset, xindex, yindex, transformPoint);
        }

        if (layer.data.lines) {
            raw = layer.data.lines[genome1][genome2]
            lines = transformBy(raw, xoffset, yoffset, xindex, yindex, transformLine);
        }

        if (layer.data.rects) {
            raw = layer.data.rects[genome1][genome2];
            rects = transformBy(raw, xoffset, yoffset, xindex, yindex, transformRect);
        }

        return {
            "enabled": enabled,
            "style": style,
            "points": points,
            "lines": lines,
            "rects": rects
        };
    }

    function transformBy(gridLayer, rowOffset, colOffset, rowIndex, colIndex, transform) {
        var row,
            col,
            r,
            c,
            data,
            transformed,
            collection = [];

        for(row in gridLayer) {
            for (col in gridLayer[row]) {
                r = rowOffset[rowIndex[row]];
                c = colOffset[colIndex[col]];
                data = gridLayer[row][col];

                if (data !== undefined)  {
                    transformed = data.map(transform.bind(undefined, r, c));
                    collection = [].concat.call(collection, transformed);
                }
            }
        }

        return collection;
    }

    function toCanvasRow(offset, position) {
        return offset - position;
    }

    function transformPoint(rowoffset, colOffset, point) {
        return {
            x: point[0] + rowOffset,
            y: toCanvasRow(colOffset, point[1]),
        };
    }

    function transformLine(rowOffset, colOffset, line) {
        return {
            x1: line[0] + rowOffset,
            x2: line[1] + rowOffset,
            y1: toCanvasRow(colOffset, line[2]),
            y2: toCanvasRow(colOffset, line[3])
        };
    }

    function transformRect(rowOffset, colOffset, rect) {
        var width = Math.abs(rect[0] - rect[1]),
            height = Math.abs(rect[2] - rect[3]),
            x = Math.min(rect[0], rect[1]);
            y = Math.max(rect[2], rect[3]);

        return [x + rowOffset, colOffset - y, width, height];
    }

    return my;
}

synmap.dropdown = function(element, datasets) {
    var my = {},
        listeners = [],
        el = $(element);

    el.on("change", function(e) {
        return my.handleSelected(datasets[$(this).val()].dataset);
    });

    my.selected = function(callback) {
        listeners.push(callback);
    };

    my.handleSelected = function(selected) {
        listeners.forEach(function(callback) {
            callback(selected);
        });
    };

    my.datasets = function() {
        return datasets;
    }

    var options = datasets.map(function(dataset, index) {
        return $("<option>").attr("value", index)
            .html(dataset.title);
    });

    el.append(options);

    return my;
}

function Histogram(dataset) {
    var my = {};

    my.get = function() {
        return dataset;
    };

    return my;
}

function HistogramController(element, model, config) {
    var histogram = bioplot.histogram(),
        scheme,
        my = {};

    my.setModel = function(newModel) {
        model = newModel;
        setTimeout(update, 200);
    };

    my.setColorScheme = function(colors) {
        histogram.configure({colors: colors});
        setTimeout(update, 200);
    };

    my.select = function(selection) {
        model.get({selected: selection});
    }

    my.render = function() {
        if(!model) return;
        update();
    }

    function update() {
        var bins = histogram.bin(model.get());
        histogram(element, bins);
    }

    return my;
}
