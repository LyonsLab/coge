function init_table() {
    $("#prev_table").tablesorter({
		sortClassAsc: 'headerSortUp',		// Class name for ascending sorting action to header
		sortClassDesc: 'headerSortDown',	// Class name for descending sorting action to header
		headerClass: 'header',		// Class name for headers (th's)
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
}

function run_synmap(scheduled){
    //  generate_basefile([],[update_basename]);
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

    if(!pageObj.basename) {
        var run_callback = function() {
            run_synmap(scheduled);
        };

        setTimeout(run_callback, 100);
        return;
    }

    if (!has_organisms())
        return;

    if ($('#blast').val() == 5 && (feat_type1 != 1 || feat_type2 != 1) ) {
        alert('BlastP only works if both genomes have protein coding sequences (CDS) AND CDS is selected for both!');
        return;
    }

    var dagchainer = $('#dagchainer_type').filter(':checked');
    var argument_list = {
        tdd: $('#tdd').val(),
        D: $('#D').val(),
        A: $('#A').val(),
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
        regen_images: $('#regen_images')[0].checked,
        width: $('#master_width').val(),
        dagchainer_type: dagchainer.val(),
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

    // TODO: Scale polling time linearly with long running jobs
    var duration = pageObj.waittime;
    var request = window.location.href.split('?')[0];

    var readlog_callback = function () {
        read_log(pageObj.basename, pageObj.tempdir);
    };

    var start_callback = function(tiny_link, status_request) {
        pageObj.nolog=1;

        close_dialog();
        $('#synmap_dialog').dialog('open');

        update_dialog_callback = function() {
            update_dialog(status_request, "#synmap_dialog", synmap_formatter);
        }

        setTimeout(update_dialog_callback, duration);

        argument_list.fname = 'go';
        argument_list.tiny_link = tiny_link;

        $.ajax({
            url: request,
            data: argument_list,
            success: function(data) {
                handle_results(data);
            }
        });

    };

    argument_list.fname = 'get_query_link';

    $.ajax({
        url: request,
        dataType: 'json',
        data: argument_list,
        success: function(data) {
            var link = "Return to this analysis: <a href=" + data.link + " onclick=window.open('tiny')"
            + "target = _new>" + data.link + "</a>";
            var logfile = '<a href="tmp/SynMap/'
            + pageObj.basename + '.log">Logfile</a>';

            $('#results').hide(0);
            $('#dialog_log').html(logfile);
            $('#synmap_link').html(link);

            start_callback(data.link, data.request);
        }
    });
}

function read_log(name, dir, callback) {
    $.ajax({
        data: {
            jquery_ajax: 1,
            fname: 'read_log',
            logfile: name,
            tempdir: dir,
        },
        success : function(data) {
            monitor_log(data);
        },
    });
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

function search_bar(div_id){
    if ($('#'+div_id).val() == "Search")
        $('#'+div_id).val("");

    if ($('#'+div_id).val() != "Search")
        $('#'+div_id).css({fontStyle: "normal"});
}

function restore_search_bar(div_id) {
  if(! $('#'+div_id).val()) {
    $('#'+div_id).val("Search").css({fontStyle: "italic"});
  }
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
    var row = $('<li>'+ item.description + '</li>');
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
        job_status.addClass('down');
        job_status.addClass('bold');
    } else if (item.status == 'skipped') {
        job_status.append(item.status);
        job_status.addClass('down');
        job_status.addClass('bold');
    } else if (item.status == 'failed') {
        job_status.append(item.status);
        job_status.addClass('alert');
        job_status.addClass('bold');
    } else {
        return;
    }

    row.append(job_status);

    /*
    if (item.status == "skipped") {
        row.append("<p>The analyses previously was generated</p>");
    }
    */

    return row;
}

function update_dialog(request, identifier, formatter) {
    var timeout = pageObj.waittime;

    var get_status = function () {
        $.ajax({
            type: 'GET',
            url: request,
            dataType: 'json',
            success: update_callback
        });
    };

    var update_callback = function(json) {
        var dialog = $(identifier);
        var workflow_status = $("<p></p>");
        var data = $("<ul></ul>");
        var results = [];
        var current_status;

        var callback = function() {
            update_dialog(request, identifier, formatter);
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

        if (json.status) {
            current_status = json.status.toLowerCase();
            workflow_status.html("Workflow status: " + json.status);
        }

        if (current_status == "completed") {
            dialog.find('#progress').hide();
            dialog.find('#dialog_success').slideDown();
        } else if (current_status == "failed" || current_status == "error"
                || current_status == "terminated") {
            dialog.find('#progress').hide();
            dialog.find('#dialog_error').slideDown();
        } else if (current_status == "notfound") {
            setTimeout(callback, timeout);
            return;
        } else {
            setTimeout(callback, timeout);
        }

        results.push(workflow_status);
        data.append(results);
        dialog.find('#text').html(data);
    };

    get_status();
}

function monitor_log(log)
{
    var waittime = pageObj.waittime;
    var fasta = 0;
    var blast = 0;
    var blastdb = 0;
    var bed = 0;
    var tandem = 0;
    var converting1 =0;
    var evalue = 0;
    var dag = 0;
    var converting2=0;
    var ks = 0;
    var results = 0;
    var match;
    pageObj.finished = 0;

    if (waittime < 60000) {
        pageObj.waittime = pageObj.waittime * 1.25;
    }

    if (log) {
        if(log.match(/fasta sequence/i))
            fasta="Generating fasta files . . . ";
        if(log.match(/Fasta creation passed/i)) {
            fasta += "done!<br/>";
            blastdb = "Generating blastable databases . . . ";
        }
        if(log.match(/BlastDB creation passed/i)) {
            blastdb += "done!<br/>";
            blast = "Running genome comparison . . . ";
        }
        if(log.match(/Completed blast run/i)) {
            blast += "done!<br/>";
            bed = "Creating .bed files . . .";
        }
        if(log.match(/Filtering results of tandem/i)) {
            bed += "done!<br/>";
            tandem = "Filtering tandem dups . . . ";
        }
        if(log.match(/Converting blast file to dagchainer input/i)) {
            tandem += "done!<br/>";
            converting1 = "Formatting for DagChainer . . . ";
        }
        if(log.match(/Adjusting evalue/i)) {
            converting1 += "done!<br/>";
            evalue = "Adjusting e-values . . . ";
        }

        if(log.match(/Running DagChainer/i)) {
            evalue += "done!<br/>";
            dag = "Running DAGChainer . . . ";
        }

        if(log.match(/Completed dagchainer run/i)) {
            dag += "done!<br/>";
            results = "Generating images . . . ";
            pageObj.finished = 0;
        }
        if(log.match(/Generating ks data/i)) {
            ks = "Calculating synonymous changes (slow) . . . ";
                results = 0;
            pageObj.finished = 0;
        }
        if(ks && log.match(/Completed generating ks data/i)) {
            ks += "done!<br/>";
            results = "Generating images . . . ";
            pageObj.finished = 0;
        }

        if(log.match(/#finished/i)) {
            pageObj.finished = 1;
        }
    } else {
        pageObj.nolog += 1;
    }

    var message = "Initializing search . . . ";

    if (fasta) message += "done!<br/>"+fasta;
    if (blastdb) message += blastdb;
    if (blast) message += blast;
    if (bed) message += bed;
    if (tandem) message += tandem;
    if (converting1) message += converting1;
    if (evalue) message += evalue;
    if (dag) message += dag;
    if (ks) message += ks;
    if (results) message += results;

    if (!pageObj.finished && pageObj.nolog<8) {
        var waittime = pageObj.waittime;
        message+="<br/><br/>Next progress check in "+Math.floor(waittime/1000)+" seconds.";

        // TODO: Scale polling time linearly with long running jobs
        var duration = pageObj.waittime;
        var readlog_callback = function() {
            read_log(pageObj.basename, pageObj.tempdir);
        };

        if ($('#synmap_dialog').dialog('isOpen')) {
            setTimeout(readlog_callback, duration);
        }
    }

    if (pageObj.finished == 0 && pageObj.waittime > 3) {
    }

    if (message) $('#text').html(message);
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

