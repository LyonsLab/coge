/*global $,window */

//
// Initialize jQuery AJAX
//
$.ajaxSetup({
        type: "GET",
        url: "OrganismView.pl",
        dataType: "html",
        cache: false,
});

var pageObj = {};

/*
 * Accessory functions
 */

//FIXME: Wrapper around window.open due to perl cgi ajax limitation
function send_to_url(url){
    window.open(url);
}

function comparator(a, b) {
    return (a > b) ? 1 : (a < b) ? -1 : 0;
}

/*
 * JQuery extensions
 */
$.fn.sortSelect = function(){
    this.each(function() {
        var i,
            opts = this.options,
            sortArray = [];

        if(this.nodeName.toLowerCase() !== "select") {
            return;
        }

        for(i = 0; i < opts.length; i++) {
            sortArray[i] = {
                v: opts[i].value,
                t: opts[i].text,
                d: opts[i].id,
            };
        }

        sortArray.sort(function(obj1, obj2) {
            var obj1t = obj1.t.toLowerCase(),
                obj2t = obj2.t.toLowerCase();

            return comparator(obj1t, obj2t);
        });

        for(i = 0; i < opts.length; i++) {
            opts[i].id = sortArray[i].d;
            opts[i].text = sortArray[i].t;
            opts[i].value = sortArray[i].v;
        }
    });
    return this;
};

$.fn.getLength = function(val){
    var opt_length;
    var blastable;

    this.each(function() {
        var i,
            opts = this.options;

        opt_length = opts.length;

        if (opt_length === 0) {
            return opt_length;
        }

        blastable = opts[0].id;

        if (val) {
            for(i = 1; i <opts.length; i++) {
                blastable += ","+opts[i].id;
                //need to chopoff last comma
            }
        }
    });

    if (val) {
        return blastable;
    }

    return opt_length;
};

/*
 * Genome List
 */

function counting() {
    var count=0;
    if($('#blank').val()){ count =0;}
    else { count = $('#genomelist_choice')[0].length;}
        if(count==0) { $('#genomelist_choice').html('<option id="blank" value=null>No genome selected</option>') ; }
    $('#count').html('Genome Count:'+count);
}


function add_to_genomelist_nowarn(genome_name,gstid)
{
    var check = $('#'+gstid).val();
    if (check){

        return;
    }else{
        var html = '<option id='+gstid+' value='+gstid+ ' >'+genome_name+'</option>';
        //alert(html);
        $('#blank').remove();
        $('#genomelist_choice').append(html);
        counting();
    }
}

function add_to_genomelist(genome_name,gstid) {
    var check = $('#'+gstid).val();
    if (check){
        alert('You have already added '+genome_name+'.');
        return;
    }
    var html = '<option id='+gstid+' value='+gstid+ ' >'+genome_name+'</option>';
    //alert(html);
    $('#blank').remove();
    $('#genomelist_choice').append(html);
    counting();

}

function add_all_genomes(){
    $('#dsg_id option').each(function(){
            add_to_genomelist_nowarn($(this).text(),$(this).attr("value"));
    });
    $('#geno_list').dialog('option', 'width', 500).dialog('open');
}

function remove_selected_genomes(){
    $('#genomelist_choice option:selected').each(function() {
            //$('#'+$(this).val()).remove();
            $(this).remove();
    });
    counting();
}


function clear_genome_list(){
    var i,
        listlength = $('#genomelist_choice')[0].length;

    for(i=0; i < listlength; i++) {
        $('#'+$('#genomelist_choice')[0][0].id).remove();
    }
    counting();
}


function add_all_org(){
    $("#org_id option").each(function(){
            get_genome_list_for_org(['args__oid','args__'+$(this).attr('value')],[add_all_genomes_from_org]);
    });
}

function add_all_genomes_from_org(val){
    if (val.length>0) {
        var genomes = val.split("&&");
        var i=0;
        for(i=0; i < genomes.length;i++) {
            var genome = genomes[i].split('%%');
            if (genome.length>1) {
                add_to_genomelist_nowarn(genome[1],genome[0]);
            }
        }
    }
}

function send_to_GenoList(){
    var check = $('#genomelist_choice').getLength();
    if(($('#blank').val())||(check==0)){
        alert('You have not selected any Genomes to examine. You must select at least one.');
        return;
    }
    var genolist = $('#genomelist_choice').getLength(1);
    parse_for_GenoList(['args__'+genolist],[send_to_url]);
}


/*
 * Searching
 */

function get_organism_chain(val) {
    get_orgs(['args__name','args__'+val, 'args__desc','args__'+val, 'args__dsgid','dsg_id'], [get_org_info_chain]);
}

function get_org_info_chain(val, val2) {
    if (val) $('#org_list').html(val);
    if (val2*1) {
        $('#org_count').html(val2);
        $('#allorg').show();
    }
    else $('#allorg').hide();
    get_org_info(['args__oid','org_id'],[genome_chain]);
}

function get_dataset_name_chain(val) {
    if (val) {
        get_dataset(['args__dsname','args__'+val], [dataset_info_chain]);
        $('#org_list').html('<hidden id="org_id">');
        $('#dsg_list').html('<hidden id="dsg_id">');
        $('#org_count,#dsg_info,#dsg_count').html(' ');
    }
    else get_orgs([],[dataset_chain]);
}

function genome_chain(val) {
    if (val) $('#org_info').html(val);
    get_genomes(['args__oid','org_id','args__dsgid','dsg_id'],[genome_info_chain]);
    $('#dsg_list').html('<hidden id="dsg_id">');
    $('#dsg_count,#dsg_info,#ds_count,#ds_info,#chr_list,#chr_count,#chr_info,#viewer,#get_seq').empty();
    $('#ds_list').html('<hidden id="ds_id">');
}

function genome_info_chain(val, val2) {
    if (val) $('#dsg_list').html(val);
    if (val2) $('#dsg_count').html(val2);
    $('#dsg_info').html('<span class="small alert">loading...</span>');
    get_genome_info(['args__dsgid','dsg_id'],[dataset_chain]);
}

function dataset_chain(val) {
    if (val) {$('#dsg_info').html(val);}
    get_dataset(['args__dsgid','dsg_id'],[dataset_info_chain]);
}

function recent_dataset_chain() {
    get_dataset(['args__oid','recent_org_id'],[dataset_info_chain]);
}

function dataset_info_chain(val, val2) {
    if (val) $('#ds_list').html(val);
    if (val2) $('#ds_count').html(val2);
    $('#ds_info').html('<span class="small alert">loading...</span>');
    $('#chr_list,#chr_count,#chr_info,#viewer,#get_seq').empty();
    get_dataset_info(['ds_id'],[dataset_chr_info_chain]);
}

function dataset_chr_info_chain(val, val2, val3) {
    if (val) $('#ds_info').html(val);
    if (val2) $('#chr_list').html(val2);
    if (val3) $('#chr_count').html(val3);
    $('#viewer,#get_seq,#feature_count_data').empty();
    get_dataset_chr_info(['ds_id', 'chr', 'dsg_id'],[populate_dataset_chr_info]);
}

function populate_dataset_chr_info (chr_info, viewer, get_seq){
    $('#chr_info').html(chr_info);
    $('#viewer').html(viewer);
    $('#get_seq').html(get_seq);
    setup_button_states();
    $('#busy').animate({opacity:0});
    $("._orgviewresult").fadeIn();
}

function timing(val){
    var searchterm = $('#org_search').val();

    if (pageObj.time) {
        clearTimeout(pageObj.time);
    }

    if (!searchterm || searchterm.length < 3 || !val) {
        $("._orgviewresult").hide();
        $('#busy').animate({opacity:0});
    }
    else {
        $('#busy').animate({opacity:1});
        pageObj.time = setTimeout(function() {
            get_organism_chain($('#org_search').val())
        }, 500);
    }
}

/*
 * Services and Resources
 */

function get_feat_gc(opts) {
    opts = opts || {};
    chr = opts.chr;
    typeid = opts.typeid;
    min = $('#feat_gc_min').val();
    max = $('#feat_gc_max').val();
    hist_type = $('#feat_hist_type').val();
    dsid = opts.dsid;
    dsgid = $('#dsg_id').val();
    $('#gc_histogram').html('loading...');

    $.ajax({
        data: {
            dsgid: dsgid[0],
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
            $('#gc_histogram').html(data);
        }
    });
}

function edit_genome_info (dsgid) {
    $.ajax({
        data: {
            jquery_ajax: 1,
            fname: 'edit_genome_info',
            dsgid: dsgid,
        },
        success : function(val) {
            $('#edit_genome_info').html(val);
            $('#edit_genome_info').dialog('open');
            get_genome_info(['args__dsgid','dsg_id'],[dataset_chain]);
        },
    });
}

function update_genome_info (dsgid) {
    dsg_name = $('#dsg_name').val();
    dsg_desc = $('#dsg_desc').val();
    dsg_ver = $('#dsg_ver').val();
    dsg_message = $('#dsg_message').val();
    dsg_link = $('#dsg_link').val();
    $.ajax({
            data: {
                jquery_ajax: 1,
                fname: 'update_genome_info',
                dsgid: dsgid,
                name: dsg_name,
                desc: dsg_desc,
                ver: dsg_ver,
                message: dsg_message,
                link: dsg_link
            },
            success : function(val) {
                if (val != 1) {alert(val);}
                $('#edit_genome_info').html('');
                $('#edit_genome_info').dialog('close');
                get_genome_info(['args__dsgid','dsg_id'], [dataset_chain]);
            },
    });
}

function export_gff () {
    var link = "bin/export/coge_gff.pl?dsgid=";
    link += $('#dsg_id').val();
    if ($('#cds_only')[0].checked) {link += ";cds=1";}
    if ($('#annos')[0].checked) {link += ";annos=1";}
    if ($('#name_unique')[0].checked) {link += ";nu=1";}
    if ($('#upa')[0].checked) {link += ";upa=1";}
    link += ";id_type="+$('#gff_id_type').val();
    $('#export_gff_link').html("<a href="+link+" target=_new>Link: "+ link+"<\a>");
    window.open(link);
}

function export_tbl () {
    var link = "bin/export/export_NCBI_TBL.pl?dsgid=";
    link += $('#dsg_id').val();
    window.open(link);
}

function export_bed () {
    var link = "bin/export/coge2bed.pl?gid=";
    link += $('#dsg_id').val();
    window.open(link);
}

function launch_seqview(dsgid, chr, dsid) {
    start = $('#start').val();
    stop = $('#stop').val();
    window.open("SeqView.pl?dsgid="+dsgid+"&dsid="+dsid+"&chr="+chr+"&start="+start+"&stop="+stop);
}

function launch_viewer (dsgid, chr) {
    x = $('#x').val();
    z = $('#z').val();
    if (z > 12) z=12;
    link = "GenomeView.pl?z="+z+"&x="+x+"&gid="+dsgid+"&chr="+chr;
    window.open(link);
}

/*
 * DOM manipulation functions
 */

function open_aa_usage_table (html) {
    $('#aa_usage_table').html(html);
    $('#aa_table').tablesorter();
}
