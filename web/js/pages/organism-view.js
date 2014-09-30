/*global $,window */

/*
 * Initialization
 */

var page;
var SPINNER = '<img src="picts/ajax-loader.gif"/>';

function init(params) {
	page = params;
	
	// Setup jQuery AJAX
	$.ajaxSetup({
        type: "GET",
        url: "OrganismView.pl",
        dataType: "json",
        cache: false,
	});
	
	// # Setup dialog popups
	$(".dialog").dialog({ autoOpen: false});
    
    // Auto-show results
    if (page.show_results) {
	    //search_organisms();
	    $('#results').fadeIn();
    } else {
	    $('#results').hide();
	    $('#no_results').html('Enter a search term').show();
    }
	
    // Setup event handlers
	$('#org_list > select').on('change', '', org_list_event);
	$('#genome_list > select').on('change', '', genome_list_event);
	$('#ds_list > select').on('change', '', dataset_list_event);	
	$('#chr_list > select').on('change', '', chr_list_event);
}

/*
 * Event Handlers
 */

function org_list_event() {
	console.log('org_list_event');
	var oid = $('#org_list > select').find(':selected').val();
	if (oid) {
		page.oid = oid;
		get_org_info_chain();
	}
}

function genome_list_event() {
	console.log('genome_list_event');
	var gid = $('#genome_list > select').find(':selected').val();
	if (gid) {
		page.gid = gid;
		genome_info_chain();
	}
}

function dataset_list_event() {
	console.log('dataset_list_event');
	var gid = $('#dataset_list > select').find(':selected').val();
	if (gid) {
		page.gid = gid;
		dataset_info_chain();
	}
}

function chr_list_event() {
	console.log('chr_list_event');
	var chr = $('#chr_list > select').find(':selected').val();
	if (chr) {
		page.chr = chr;
		chr_info_chain();
	}
}

/*
 * Searching
 */

//FIXME: Remove name and desc for search
function get_organism_chain(search_term) {
	console.log("get_organism_chain");

	$('#org_count,org_list').html("");
    
    return $.ajax({
        data: {
            fname: "get_orgs",
            jquery_ajax: 1,
            name: search_term,
            desc: search_term,
            //gid: page.gid
        },
        success: function(response) {
            if (response.organisms)
            	$('#org_list > select').html(response.organisms);
            if (response.count * 1) {
                $('#org_count').html(response.count);
                $('#allorg').show();
            } else {
                $('#allorg').hide();
            }
        },
        error: function () {
            $('#org_list > select').html('<span class="small alert">The results could not be loaded.</span>');
        }
    }).always(function() {
    	var org_count = $('#org_count').html();
    	console.log('orgs: '+org_count);
    	if (org_count) {
    		var oid = $('#org_list > select').find(':selected').val();
    		if (!oid)
    			oid = $('#org_list > select').first().val();
    		page.oid = oid;
    		get_org_info_chain();
    		$('#no_results').hide();
    		$('#results').show();
    	}
    	else {
    		$('#results').hide();
    		$('#no_results').html('No matching organisms were found').show();
    	}
    });
}

function get_org_info_chain() {
	console.log('get_org_info_chain '+page.oid);
    $.ajax({
        data: {
            fname: "get_org_info",
            jquery_ajax: 1,
            oid: page.oid
        },
        success: function(response) {
            if (response.organism)
                $('#org_info').html(response.organism);
            if (response.error)
                $('#org_info').html(response.error);
            $('#genome_list > select').html('');
            $('#genome_count,#genome_info,#ds_count,#ds_info,#chr_list > select,#chr_count,#chr_info,#viewer,#get_seq').empty();
            //$('#genome_list > select,#genome_info,#ds_list > select,#ds_info,#chr_info').html(SPINNER); // TODO
        },
        error: function () {
            $('#org_info').html('<span class="small alert">The results could not be loaded.</span>');
        }
    }).always(genome_chain);
}

function genome_chain() {
	console.log('genome_chain');
    $.ajax({
        data: {
            fname: "get_genomes",
            jquery_ajax: 1,
            oid: page.oid,
            //gid: page.gid
        },
        success: function(response) {
        	if (!response) return;
            if (response.error) 
            	$('#genome_list > select').html(response.error);
            if (response.genomes) {
            	$('#genome_list > select').html(response.genomes);
            	page.gid = response.selected_id;
            }
            if (response.count) 
            	$('#genome_count').html(response.count);
            $('#ds_list > select,#genome_info,#ds_info,#chr_info').html(SPINNER);
        },
        error: function () {
            $('#genome_list > select').html('<span class="small alert">The results could not be loaded.</span>');
        }
    }).always(genome_info_chain);
}

function genome_info_chain() {
	console.log('genome_info_chain gid='+page.gid);
    $.ajax({
        data: {
            fname: "get_genome_info",
            jquery_ajax: 1,
            gid: page.gid
        },
        success: function (response) {
            if (response.genome) {
                $('#genome_info').html(response.genome);
            } else {
                $('#genome_info').html(response.error);
            }
        },
        error: function () {
            $('#genome_info').html('<span class="small alert">The results could not be loaded.</span>');
        }
    }).always(dataset_chain);
}

function dataset_chain() {
	console.log('dataset_chain gid='+page.gid);
    $.ajax({
        data: {
            fname: "get_datasets",
            jquery_ajax: 1,
            gid: page.gid
        },
        success: function (response) {
            if (response.error) $('#ds_list > select').html(response.error);
            if (response.datasets) $('#ds_list > select').html(response.datasets);
            if (response.count) $('#ds_count').html(response.count);
            $('#ds_info').html(SPINNER);//'<span class="small alert">loading...</span>');
            $('#chr_info').html(SPINNER);//'<span class="small alert">loading...</span>');
            $('#chr_list > select').html('<span class="small alert">loading...</span>');
            $('#chr_count,#viewer,#get_seq').empty();
            page.dsid = response.selected_id;
        },
        error: function (response) {
        	console.log(response);
            $('#ds_list > select').html('<span class="small alert">Datasets could not be loaded.</span>');
        }
    }).always(dataset_info_chain);
}

function dataset_info_chain() {
	console.log('dataset_info_chain '+page.dsid);
    $.ajax({
        data: {
            jquery_ajax: 1,
            fname: "get_dataset_info",
            jquery_ajax: 1,
            dsid: page.dsid
        },
        success: function (response) {
        	console.log(response);
        	//if (response.error) $('#ds_info').html(response.error);
            if (response.dataset) $('#ds_info').html(response.dataset);
            if (response.chromosomes) $('#chr_list > select').html(response.chromosomes);
            if (response.count) $('#chr_count').html(response.count);
            $('#viewer,#get_seq,#feature_count_data').empty();
            page.chr = response.selected_chr;
        },
        error: function() {
            $('#ds_info').html('<span class="small alert">The results could not be loaded.</span>');
            $('#chr_list > select').html('<span class="small alert">The results could not be loaded.</span>');
        }
    }).always(chr_info_chain);
}

function chr_info_chain() {
	console.log('chr_info_chain '+page.chr);
    $.ajax({
        data: {
            fname: "get_chr_info",
            jquery_ajax: 1,
            dsid: page.dsid,
            chr: page.chr
        },
        success: function (response) {
        	console.log(response);
            if (response.chr_info) {
            	 $('#chr_info').html(response.chr_info).show();
            	 $('#chr_info').parent().show(); // see kludge in tmpl
            	 $('#viewer').html(response.viewer).show();
                 $('#get_seq').html(response.seqview).show();
            }
            else {
            	$('#chr_info').hide();
            	 $('#viewer').hide();
                 $('#get_seq').hide();
            }
            setup_button_states();
        },
        error: function() {
            $('#chr_info').html('<span class="small alert">The results could not be loaded.</span>').show();
        }
    });
}

function search_organisms() {
    var search_term = $('#org_search').val();
    
    if (pageObj.time)
        clearTimeout(pageObj.time);
    
    if (!search_term || search_term.length < 3) {
	        $('#busy').animate({opacity:0});
    }
    else {
        $('#busy').animate({opacity:1});
        pageObj.time = setTimeout(function() {
        	var search_term = $('#org_search').val();
            get_organism_chain(search_term)
            	.always(function() {
            		$('#busy').animate({opacity:0});
            	});
        }, 500);
    }
}

/*
 * Accessory functions
 */

function comparator(a, b) {
    return (a > b) ? 1 : (a < b) ? -1 : 0;
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
    dsgid = page.gid;
    $('#gc_histogram').html('loading...');

    $.ajax({
        dataType: "html",
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

/*
 * Export
 */

function export_gff () {
    var link = "bin/export/coge_gff.pl?dsgid=" + page.gid;
    if ($('#cds_only')[0].checked) {link += ";cds=1";}
    if ($('#annos')[0].checked) {link += ";annos=1";}
    if ($('#name_unique')[0].checked) {link += ";nu=1";}
    if ($('#upa')[0].checked) {link += ";upa=1";}
    link += ";id_type="+$('#gff_id_type').val();
    $('#export_gff_link').html("<a href="+link+" target=_new>Link: "+ link+"</a>");
    window.open(link);
}

function export_tbl () {
    var link = "bin/export/export_NCBI_TBL.pl?dsgid=";
    link += page.gid;
    window.open(link);
}

function export_bed () {
    var link = "bin/export/coge2bed.pl?gid=";
    link += page.gid;
    window.open(link);
}

/*
 * Links
 */

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

function open_aa_usage_table (html) {
    $('#aa_usage_table').html(html);
    $('#aa_table').tablesorter();
}


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


//function add_all_org(){
//    $("#org_id option").each(function(){
//            get_genome_list_for_org(['args__oid','args__'+$(this).attr('value')],[add_all_genomes_from_org]);
//    });
//}

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
    window.open("GenomeList.pl?dsgid=" + genolist, "_self");
}

/*
 * jQuery extensions
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

