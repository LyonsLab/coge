/*global $,window */

/*
 * Initialization
 */

var page;
var SPINNER = '<img class="padded" src="picts/ajax-loader.gif"/>';

function init(params) {
	page = params;
	console.log(page);
	
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
    if (page.init_msg) {
    	$('#results').hide();
  	    $('#no_results').html(page.init_msg).show();
    } else {
    	$('#results').fadeIn();
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
		organism_info_chain();
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
	var dsid = $('#ds_list > select').find(':selected').val();
	if (dsid) {
		page.dsid = dsid;
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
function organism_chain(search_term) {
	console.log("organism_chain");

    return $.ajax({
        data: {
            fname: "get_orgs",
            name: search_term,
            desc: search_term,
            //gid: page.gid
        },
        success: function(response) {
            if (response.error)
            	return;
            
            if (response.organisms) {
            	// We got orgs to show
            	$('#org_list > select').html(response.organisms);
            	$('#org_count').html(response.count);
        		return;
            }
        },
        error: function () {
            $('#no_results').html('There was an error and the results could not be loaded.');
            // TODO add more detailed error message
        }
    }).always(function(response) {
    	if (!response.error) {
    		page.oid = response.selected_id;
    		organism_info_chain();
    		
    		$('#no_results').hide();
    		$('#results').show();
    	}
    	else {
    		$('#results').hide();
    		$('#no_results').html(response.error).show();
    	}
    });
}

function organism_info_chain() {
	console.log('organism_info_chain '+page.oid);
    $.ajax({
        data: {
            fname: "get_org_info",
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
	$('#genome_info,#ds_info,#chr_info').html(SPINNER).show();
	$('#genome_list > select,#ds_list > select,#chr_list > select').html('<option>loading...</option>').show();
	$('#genome_info,#ds_list,#ds_info,#chr_list,#chr_info').parent().show();
	
    $.ajax({
        data: {
            fname: "get_genomes",
            oid: page.oid,
            //gid: page.gid
        },
        success: function(response) {
        	if (!response) return;
            if (response.error) 
            	$('#genome_list > select').html('<option>'+response.error+'</option>');
            else if (response.genomes) {
            	$('#genome_list > select').html(response.genomes);
            	page.gid = response.selected_id;
            }
            if (response.count) 
            	$('#genome_count').html(response.count);
            //$('#ds_list > select').html('<option>loading...</option>');
            $('#ds_list > select').html('<option>loading...</option>');
            $('#genome_info,#ds_info,#chr_info').html(SPINNER);
        },
        error: function () {
            $('#genome_list > select').html('<span class="small alert">The results could not be loaded.</span>');
        }
    }).always(function() {
    	if ($('#genome_count').html()) {
    		$('#genome_info').parent().show();
    		genome_info_chain();
    	}
    	else {
    		$('#genome_info,#ds_list,#ds_info,#chr_list,#chr_info').parent().hide();
    	}
    });
}

function genome_info_chain() {
	console.log('genome_info_chain dsid='+page.gid);
	if (page.gid) {
	    $.ajax({
	        data: {
	            fname: "get_genome_info",
	            gid: page.gid
	        },
	        success: function (response) {
	            if (response.genome) 
	                $('#genome_info').html(response.genome).show();
	            else 
	                $('#genome_info').html(response.error).show();
	            $('#genome_info').parent().show();
	        },
	        error: function () {
	            $('#genome_info').html('<span class="small alert">The results could not be loaded.</span>').show();
	            $('#genome_info').parent().show();
	        }
	    }).always(dataset_chain);
	}
	else {
		$('#genome_info,#ds_list,#ds_info,#chr_list,#chr_info').parent().hide();
	}
}

function dataset_chain() {
	console.log('dataset_chain gid='+page.gid);
    $.ajax({
        data: {
            fname: "get_datasets",
            gid: page.gid
        },
        success: function (response) {
            if (response.error) $('#ds_list > select').html(response.error);
            if (response.datasets) $('#ds_list > select').html(response.datasets);
            if (response.count) $('#ds_count').html(response.count);
            $('#ds_info').html(SPINNER);//'<span class="small alert">loading...</span>');
            $('#chr_info').html(SPINNER);//'<span class="small alert">loading...</span>');
            $('#chr_list > select').html('<option>loading...</option>');
            $('#chr_count,#viewer,#get_seq').empty();
            page.dsid = response.selected_id;
            
            $('#ds_list').parent().show();
        },
        error: function (response) {
            $('#ds_list > select').html('<span class="small alert">Datasets could not be loaded.</span>');
        }
    }).always(dataset_info_chain);
}

function dataset_info_chain() {
	console.log('dataset_info_chain '+page.dsid);
    $.ajax({
        data: {
            fname: "get_dataset_info",
            dsid: page.dsid
        },
        success: function (response) {
        	//if (response.error) $('#ds_info').html(response.error);
            if (response.dataset) $('#ds_info').html(response.dataset);
            if (response.chromosomes) $('#chr_list > select').html(response.chromosomes);
            if (response.count) $('#chr_count').html(response.count);
            $('#viewer,#get_seq,#feature_count_data').empty();
            page.chr = response.selected_chr;
            
            $('#ds_info,#chr_list').parent().show();
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
            gid: page.gid,
            dsid: page.dsid,
            chr: page.chr
        },
        success: function (response) {
            if (response.chr_info) {
            	 $('#chr_info').html(response.chr_info).show();
            	 $('#chr_info').parent().show(); // see kludge in tmpl
            	 $('#viewer').html(response.viewer).show();
                 $('#get_seq').html(response.seqview).show();
            }
            else {
                 $('#chr_info').hide();
                 $('#chr_info').parent().hide();
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
            organism_chain(search_term)
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
    var link = "coge_gff.pl?dsgid=" + page.gid;
    if ($('#cds_only')[0].checked) {link += ";cds=1";}
    if ($('#annos')[0].checked) {link += ";annos=1";}
    if ($('#name_unique')[0].checked) {link += ";nu=1";}
    if ($('#upa')[0].checked) {link += ";upa=1";}
    link += ";id_type="+$('#gff_id_type').val();
    $('#export_gff_link').html("<a href="+link+" target=_new>Link: "+ link+"</a>");
    window.open(link);
}

function export_tbl () {
    var link = "export_NCBI_TBL.pl?dsgid=";
    link += page.gid;
    window.open(link);
}

function export_bed () {
    var link = "coge2bed.pl?gid=";
    link += page.gid;
    window.open(link);
}

/*
 * Links
 */

function launch_seqview() {
    start = $('#start').val();
    stop = $('#stop').val();
    window.open("SeqView.pl?dsgid="+page.gid+"&dsid="+page.dsid+"&chr="+page.chr+"&start="+start+"&stop="+stop);
}

function launch_viewer() {
    var x = $('#x').val();
    var z = $('#z').val();
    if (z > 12) z=12;
    link = "GenomeView.pl?gid="+page.gid+"&loc="+page.chr+':'+x;
    window.open(link);
}

function open_aa_usage_table (html) {
    $('#aa_usage_table').html(html);
    $('#aa_table').tablesorter();
}


/*
 * Genome List
 */

function get_genomelist_count() {
    var count = 0;
    
    if ($('#genomelist_blank').val())
    	count = 0;
    else
    	count = $('#genomelist')[0].length;
    
    if (count == 0) 
    	$('#genomelist').html('<option id="genomelist_blank" value=null>No genome selected</option>');
    
    $('#genomelist_count').html(count);
    
    return count;
}


function add_to_genomelist_nowarn(genome_name,gstid)
{
    var check = $('#'+gstid).val();
    if (check)
        return;
    else {
        $('#genomelist_blank').remove();
        $('#genomelist').append('<option id='+gstid+' value='+gstid+ ' >'+genome_name+'</option>');
        get_genomelist_count();
    }
}

function add_genome_to_list(gid) { // mdb added 10/2/14 COGE-512
	if (!gid) {
		alert('Internal error (add_genome_to_list)');
		return;
	}
	
	if ( $('#genomelist option[value="'+gid+'"]').html() ) {
		alert('You have already added this genome.');
        return;
	}
	
	$.ajax({
        data: {
        	fname: 'get_genome_name',
            gid: gid,
        },
        success: function(response) {
        	if (!response) {
        		alert("Error: couldn't retrieve genome");
        		return;
        	}
        	else if (response.error || !response.info) {
        		alert('Error: '+response);
        		return;
        	}
        	
            $('#genomelist_blank').remove();
            $('#genomelist').append('<option title="'+response.info+'" value="'+gid+'">'+response.info+'</option>');
            get_genomelist_count();
            
            $('#dialog_genomelist').dialog('option', 'width', '32.5em').dialog('open');
        }
    });
}

function add_all_genomes(){
    $('#dsg_id option').each(function(){
        add_to_genomelist_nowarn($(this).text(),$(this).attr("value"));
    });
    $('#dialog_genomelist').dialog('option', 'width', '32.5em').dialog('open');
}

function remove_selected_genomes(){
    $('#genomelist option:selected').each(function() {
            //$('#'+$(this).val()).remove();
            $(this).remove();
    });
    get_genomelist_count();
}


function clear_genome_list(){
    var i, num_genomes = get_genomelist_count();
    
    $('#genomelist').empty();
    
    get_genomelist_count();
}

function add_all_genomes_from_org(val) {
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

function send_to_genomelist() {
    var num_genomes = $('#genomelist').getLength();
    if ( $('#genomelist_blank').val() || (num_genomes == 0) ) {
        alert('You have not selected any genomes to add. You must select at least one.');
        return;
    }
    
    var list = $('#genomelist option');
    var ids = [];
    for (var i = 0;  i < num_genomes;  i++)
    	ids.push(list[i].value);
    gid_list = ids.join(",");
    window.open("GenomeList.pl?dsgid=" + gid_list, "_self");
}

/*
 * jQuery extensions
 */
$.fn.sortSelect = function() {
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

