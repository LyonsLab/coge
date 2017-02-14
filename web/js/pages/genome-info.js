
function wait_to_search (search_func, search_term) {
    pageObj.search_term = search_term;

    if (pageObj.time) {
        clearTimeout(pageObj.time);
    }

    // FIXME: could generalize by passing select id instead of separate search_* functions
    pageObj.time = setTimeout(
        function() {
            search_func(pageObj.search_term);
        },
        500
    );
}

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
            dsgid: GENOME_ID,
            gstid: $("#gstid").val()
        },
        success: function(data) {
            elem.fadeOut(function() {
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
            dsgid: GENOME_ID,
            fname: 'get_gc_for_feature_type',
            dsid: dsid,
            typeid: typeid,
            chr: chr,
            min: min,
            max: max,
            hist_type: hist_type
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
           my: "top",
           at: "top",
           of: window
        }
    });

    $('#chromosome_hist').html('loading. . .');
    $.ajax({
        data: {
            fname: 'get_chr_length_hist',
            dsgid: dsgid
        },
        success: function (data) {$('#chromosome_hist').html(data);}
    });
}

var chr_list_table = null;
function chr_list() {
	$('#chromosome_list').dialog({
		autoOpen: true,
		position: {
			my: "top",
			at: "top",
			of: window
		}
	});

	if (!chr_list_table)
	    $.ajax({
	        data: {
	            fname: 'get_chromosomes',
	            gid: GENOME_ID
	        },
	        dataType: "json",
	        success: function (data) {
				chr_list_table = $('#chr_list_table').DataTable({
			   		columnDefs:[{targets:0,type:'natural'},{targets:1,type:'num'},{orderable:false,targets:[2,3,4]}],
			   		data:data,
					deferRender:true,
					language:{info:'Showing _START_ to _END_ of _TOTAL_ chromosomes'},
					lengthChange:false,
					paging:data.length>10,
					pagingType:'simple',
					searching:data.length>10
				});
				$('#chr_list_loading').hide();
				$('#chr_list').show();
	        }
		});
 }

(function() { //FIXME move into module (mdb 9/27/16)
	 
	/*
	 * Natural Sort algorithm for Javascript - Version 0.7 - Released under MIT license
	 * Author: Jim Palmer (based on chunking idea from Dave Koelle)
	 * Contributors: Mike Grier (mgrier.com), Clint Priest, Kyle Adams, guillermo
	 * See: http://js-naturalsort.googlecode.com/svn/trunk/naturalSort.js
	 */
	function naturalSort (a, b) {
	    var re = /(^-?[0-9]+(\.?[0-9]*)[df]?e?[0-9]?$|^0x[0-9a-f]+$|[0-9]+)/gi,
	        sre = /(^[ ]*|[ ]*$)/g,
	        dre = /(^([\w ]+,?[\w ]+)?[\w ]+,?[\w ]+\d+:\d+(:\d+)?[\w ]?|^\d{1,4}[\/\-]\d{1,4}[\/\-]\d{1,4}|^\w+, \w+ \d+, \d{4})/,
	        hre = /^0x[0-9a-f]+$/i,
	        ore = /^0/,
	        // convert all to strings and trim()
	        x = a.toString().replace(sre, '') || '',
	        y = b.toString().replace(sre, '') || '',
	        // chunk/tokenize
	        xN = x.replace(re, '\0$1\0').replace(/\0$/,'').replace(/^\0/,'').split('\0'),
	        yN = y.replace(re, '\0$1\0').replace(/\0$/,'').replace(/^\0/,'').split('\0'),
	        // numeric, hex or date detection
	        xD = parseInt(x.match(hre), 10) || (xN.length !== 1 && x.match(dre) && Date.parse(x)),
	        yD = parseInt(y.match(hre), 10) || xD && y.match(dre) && Date.parse(y) || null;
	 
	    // first try and sort Hex codes or Dates
	    if (yD) {
	        if ( xD < yD ) {
	            return -1;
	        }
	        else if ( xD > yD ) {
	            return 1;
	        }
	    }
	 
	    // natural sorting through split numeric strings and default strings
	    for(var cLoc=0, numS=Math.max(xN.length, yN.length); cLoc < numS; cLoc++) {
	        // find floats not starting with '0', string or 0 if not defined (Clint Priest)
	        var oFxNcL = !(xN[cLoc] || '').match(ore) && parseFloat(xN[cLoc], 10) || xN[cLoc] || 0;
	        var oFyNcL = !(yN[cLoc] || '').match(ore) && parseFloat(yN[cLoc], 10) || yN[cLoc] || 0;
	        // handle numeric vs string comparison - number < string - (Kyle Adams)
	        if (isNaN(oFxNcL) !== isNaN(oFyNcL)) {
	            return (isNaN(oFxNcL)) ? 1 : -1;
	        }
	        // rely on string comparison if different types - i.e. '02' < 2 != '02' < '2'
	        else if (typeof oFxNcL !== typeof oFyNcL) {
	            oFxNcL += '';
	            oFyNcL += '';
	        }
	        if (oFxNcL < oFyNcL) {
	            return -1;
	        }
	        if (oFxNcL > oFyNcL) {
	            return 1;
	        }
	    }
	    return 0;
	}
	 
	jQuery.extend( jQuery.fn.dataTableExt.oSort, {
	    "natural-asc": function ( a, b ) {
	        return naturalSort(a,b);
	    },
	 
	    "natural-desc": function ( a, b ) {
	        return naturalSort(a,b) * -1;
	    }
	} );
}());

function download_chromosome_sequence(chr) { // FIXME use API genome/sequence (mdb 11/22/16)
    $.ajax({
        data: {
            fname: 'cache_chr_fasta',
            gid: GENOME_ID,
            chr: chr
        },
        success: function(data) {
        	document.location='get_seq_for_chr.pl?gid=' + GENOME_ID + '&chr=' + chr;
        }
    });
}

function percent_gc_at(chr, ws, irods) {
    coge.progress.begin();
    var request = {
        type: 'analyze_nucleotides',
        requester: {
            page: PAGE_NAME
        },
        parameters: {
            gid: GENOME_ID,
            chr: chr,
            ws: ws,
            irods: irods ? 1 : 0
        }
    };
    coge.services.submit_job(request) 
        .done(function(response) {
            if (!response) {
                coge.progress.failed("Error: empty response from server");
                return;
            }
            if (!response.success || !response.id) {
                coge.progress.failed("Error: failed to start workflow", response.error);
                return;
            }
            coge.progress.update(response.id, response.site_url);
        })
        .fail(function(jqXHR, textStatus, errorThrown) {
            coge.progress.failed("Couldn't talk to the server: " + textStatus + ': ' + errorThrown);
        });
}

function download_chr_file() {
	var i = $('input[name=chr]:checked');
	if (!i.length) {
		alert('Please select one of the files to download first.');
		return;
	}
	var id = i.attr('id');
    var type = id.substring(0,1);
	if (type == 'f')
		download_chromosome_sequence(id.substring(1));
	else if (type == 'g')
		get_gff(id.substring(1));
    else {
        var num_windows = chromosome_length - 10000;
        if (num_windows > 2000000)
            num_windows = chromosome_length - 2000000;
        coge.utils.prompt('Window Size:', '%GC/AT Sliding Window', chromosome_length - num_windows, function(val) {
            if (val != '') {
                var chr = id.substring(1);
                coge.progress.init({title: '%GC/AT',
                    onReset: function() {
                        document.location='get_percent_gc_at_for_chr.pl?gid=' + GENOME_ID + '&chr=' + chr + '&ws=' + val;
                    }
                });
                percent_gc_at(chr, val);
            }
        });
    }
}

function export_chr_file() {
	var i = $('input[name=chr]:checked');
	if (!i.length) {
		alert('Please select one of the files to export first.');
		return;
	}
	var id = i.attr('id');
	var type = id.substring(0,1);
	if (type == 'f')
		export_fasta_chr(id.substring(1));
	else if (type == 'g')
		export_gff(id.substring(1));
    else {
        var num_windows = chromosome_length - 10000;
        if (num_windows > 2000000)
            num_windows = chromosome_length - 2000000;
        coge.utils.prompt('Window Size:', '%GC/AT Sliding Window', chromosome_length - num_windows, function(val) {
            if (val != '') {
                coge.progress.init({title: '%GC/AT'});
                percent_gc_at(id.substring(1), val, true);
            }
        });
    }
}

// major hack, rewrite get_chromosomes in GenomeInfo.pl to get json from api, don't send chr_len to this
var chromosome_length;
function update_percent_gc_at_plot_button(chr_len) {
    chromosome_length = chr_len;
	var i = $('input[name=chr]:checked');
	var id = i.attr('id');
	var type = id.substring(0,1);
    var b = $('#percent_gc_at_plot_button');
    if (type == 'n')
        b.removeClass('coge-disabled');
    else
        b.addClass('coge-disabled');
}

function plot_percent_gc_at() {
	var i = $('input[name=chr]:checked');
	if (!i.length)
		return;
	var id = i.attr('id');
    var chr = id.substring(1);
    var num_windows = chromosome_length - 10000;
    if (num_windows > 2000000)
        num_windows = chromosome_length - 2000000;
    coge.utils.prompt('Window Size:', '%GC/AT Sliding Window', chromosome_length - num_windows, function(val) {
        if (val != '') {
            coge.progress.init({title: '%GC/AT',
                onReset: function() {
                    var div = $('<div title="%GC/AT Plot"></div>').appendTo(document.body);
                    $('<div id="percent_gc_at_plot" style="width:100%;height:100%;"><img id="percent_gc_at_busy" src="picts/ajax-loader.gif" /></div>').appendTo(div);
                    div.dialog({
                        modal: true,
                        resizable: true,
                        height: $(window).height() - 40,
                        width: $(window).width() - 40,
                        close: function() { div.remove(); },
                        resizeStop: function() { Plotly.Plots.resize($('#percent_gc_at_plot')[0]); }
                    });
                    $.ajax({
                        url: 'get_percent_gc_at_for_chr_json.pl?gid=' + GENOME_ID + '&chr=' + chr + '&ws=' + val,
                        dataType: "json",
                        success: function(json) {
                            $('#percent_gc_at_plot').empty();
                            Plotly.newPlot('percent_gc_at_plot', [{
                                name: 'AT',
                                y: json.at,
                                mode: 'lines',
                                line: { color: 'rgb(3,141,243)' }
                            },{
                                name: 'GC',
                                y: json.gc,
                                mode: 'lines',
                                line: { color: 'rgb(64,182,77)' }
                            },{
                                name: 'N',
                                y: json.n,
                                mode: 'lines',
                                line: { color: 'rgb(243,145,3)' }
                            },{
                                name: 'X',
                                y: json.x,
                                mode: 'lines',
                                line: { color: 'rgb(171,3,243)' }
                            }], {
                                title: 'Sliding Window for ' + chr,
                                xaxis: { rangeslider: {}, title: 'Window Iteration' },
                                yaxis: { title: 'Percentage (%)' }
                            });
                        },
                        error: function(err) {
                            console.log(err);
                        }
                    });
                }
            });
            percent_gc_at(chr, val);
        }
    });
}

function toggle_load_log() {
	var btn = $('#log_button');
	var log = $('#log_contents');
	var spinner = $('#log_spinner');
	
	var setVisible = function(visible) { // TODO: i want to use jQuery.toggle() instead, oh well
    	if (visible) {
			log.removeClass('hidden');
	    	btn.html('Hide');
	    	spinner.animate({opacity:0});
    	}
    	else {
    		log.addClass('hidden');
    		btn.html('Show');
    		spinner.animate({opacity:0});
    	}
	}
	
	if (log.is(":hidden")) {
		if (log.html() == '') {
			spinner.css({opacity:1});
			$.ajax({
		        data: {
		            fname: 'get_load_log',
		            gid: GENOME_ID
		        },
		        success: function(data) {
		        	if (data)
		        		log.html(data);
		        	else
		        		log.html('<span class="alert">Error: log file not found</span>');
		        	setVisible(true);
		        }
		    });
		}
		else
			setVisible(true);
	}
	else
		setVisible(false);
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

    var link = $("<a></a>")
            .addClass("bold")
            .attr("href", DISCOVERY_ENVIRONMENT.concat(irods_home_path))
            .html(irods_home_path);

    if (data.protein) {
        $('#export_log').html(
            "<br>Copying protein sequences to <br><br>"
        ).append(link);
    } else {
        $('#export_log').html(
            "<br>Copying nucleotide sequences to <br><br>"
        ).append(link);

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

function debounce_search(search_func, search_term) { //FIXME why is this a dup of wait_to_search() above? (mdb 9/27/16)
    pageObj.search_term = search_term;

    if (pageObj.time) {
        clearTimeout(pageObj.time);
    }

    // FIXME: could generalize by passing select id instead of separate search_* functions
    pageObj.time = setTimeout(
        function() {
            search_func(pageObj.search_term);
        },
        500
    );
}

function search_organisms (search_term) {
    if (search_term.length > 2) {
        $.ajax({
            data: {
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
            dsgid: GENOME_ID,
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
            fname: 'get_genome_info',
            gid: GENOME_ID
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
            fname: 'edit_genome_info',
            gid: GENOME_ID
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
            fname: 'update_genome_info',
            gid: GENOME_ID,
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
            gid: GENOME_ID,
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

function get_gff($chr) {
    var id = GENOME_ID;
    $('#gff_export').dialog('option', 'width', 400).dialog('open');

    $('#gff_submit').unbind().click(function(event, ui) {
        var id_type = $('#gff_id_type').val(),
            cds = +$('#cds_only')[0].checked,
            chr = +$chr,
            annos = +$('#annos')[0].checked,
            nu = +$('#name_unique')[0].checked,
            upa = +$('#upa')[0].checked;

        var height = $('#status_log')[0].scrollHeight;
        $("#status_log").animate({ scrollTop: height}, 500);

        var title = $("<br><h4>Generating GFF file</h4>");

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
                upa: upa,
                chr: $chr
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
            gid: GENOME_ID
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
            gid: GENOME_ID
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

function export_gff(chr) {
    $('#gff_export').dialog('option', 'width', 400).dialog('open');
    $('#gff_submit').unbind().click(function() {

        // ANALYSIS OPTIONS
        var id = GENOME_ID;
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

        var title = $("<h4>Generating gff and copying file</h4>");


        var link = $("<a></a>")
                .addClass("bold")
                .attr("href", DISCOVERY_ENVIRONMENT.concat(irods_home_path))
                .html(irods_home_path);

        $("#export_log")
            .html(title)
            .append(link)
            .append(note);

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
                upa: upa,
                chr: chr
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

function export_to_irods(file_description, fname, data) { // FIXME migrate to API (mdb 9/26/16)
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

    var link = $("<a></a>")
            .addClass("bold")
            .attr("href", DISCOVERY_ENVIRONMENT.concat(irods_home_path))
            .html(irods_home_path);

    $('#export_log').html(
        "<br>Copying " + file_description + " to <br><br>"
    ).append(link);

    if (!data)
    	data = {};
    data.fname = fname;
    data.gid = GENOME_ID;
    $.ajax({
        dataType: "json",
        data: data,
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

function export_tbl() {
	export_to_irods("this genomes's TBL file", "export_tbl");
}

function export_bed() {
	export_to_irods("this genomes's BED file", "export_bed");
}

function export_fasta() {
	export_to_irods("this genomes's FASTA file", "export_fasta");
}

function export_fasta_chr(chr) {
    $.ajax({
        data: {
            fname: 'cache_chr_fasta',
            gid: GENOME_ID,
            chr: chr
        },
        success: function(data) {
        	var obj = jQuery.parseJSON(data);
        	export_to_irods("This chromosome's FASTA file", "export_file_chr", { chr: chr, file: obj.file });
        }
    });
}

function annotate() {
    // Prevent concurrent executions
    if ( $("#annotate_dialog").dialog( "isOpen" ) )
        return;

    // Make sure user is still logged-in
    if (!check_and_report_login())
        return;

    // Open status dialog right away
    $('#annotate_log').html("<br><br>Annotating transcriptome using TransDecoder...");
    $('#annotate_dialog')
        .unbind()
        .dialog('open');
    
    $.ajax({
        dataType: "json",
        data: {
        	fname: 'annotate',
        	gid: GENOME_ID
        },
        success: function(json) {
            //$('#annotate_log').append('<br><br>File: ' + json.file);

            if (!json || json.error) {
                $('#annotate_loading_msg,#annotate_finished_msg').hide();
                $('#annotate_error_msg').fadeIn();
            } 
            else {
                $('#annotate_loading_msg').hide();
                $('#annotate_ok_button,#annotate_finished_msg,#annotate_ok_button').fadeIn();
                setTimeout(get_datasets, 0);
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
    window.history.pushState({}, "Title", PAGE_NAME + "?gid=" + GENOME_ID);
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
			'The CoGe Support Team has been notified of this error but please ' + 
			'feel free to contact us at <a href="mailto:<TMPL_VAR NAME=SUPPORT_EMAIL>"><TMPL_VAR NAME=SUPPORT_EMAIL></a> ' +
			'and we can help to determine the cause.' +
			'</div>';
		var log = $('#load_log');
		log.html( log.html() + msg );
	}

    // Update dialog
    $('#loading_msg').hide();
    $('#error_msg').fadeIn();
    $('#cancel_button').fadeIn();

    if (newLoad) { // mdb added check to prevent redundant emails, 8/14/14 issue 458
	    $.ajax({
	        data: {
	            fname: "send_error_report",
	            load_id: load_id,
	            job_id: job_id
	        }
	    });
    }
}

function load_succeeded(obj) {
    // Update globals
	GENOME_ID = obj.genome_id; // for continuing to GenomeInfo

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

    $('#load_dialog').dialog('open');
//   $('#status_log').html('Initializing ...');
    newLoad = true;
    
    $.ajax({
        dataType: "json",
        data: {
            fname: 'copy_genome',
            load_id: load_id,
            gid: GENOME_ID,
            mask: mask,
            seq_only: seq_only,
            timestamp: new Date().getTime()
        },
        success : function(data) {
            if (data.error) {
                alert(data.error);
                return;
            }

            job_id = data.job_id;

            // Set link in status dialog
            $('#loading_msg span a').attr('href', data.link).html(data.link);

            // Add load_id to browser URL
            window.history.pushState({}, "Title", PAGE_NAME + "?gid=" + GENOME_ID + "&job_id=" + data.job_id); // Add job_id to browser URL
            update_dialog("api/v1/jobs/" + data.job_id, pageObj.user, "#load_dialog", progress_formatter);
        }
        // TODO: handle error, show in status dialog
    });
}

function progress_formatter(item) {
    var msg;
    var row = $('<li>'+ item.description + ' </li>');

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

    if (item.log) {
        var p = item.log.split("\n");

        var pElements = p.map(function(item) {
            var norm = item.replace(/\\t/g, " ").replace(/\\'/g, "'");
            return $("<div></div>").append(norm);
        });

        var log = $("<div></div>").html(pElements).addClass("padded");
        row.append(log);
    }

    return row;
}

function update_dialog(request, user, identifier, formatter) {
    var get_status = function () {
        $.ajax({
            type: 'GET',
            url: request,
            dataType: 'json',
            data: {
                username: user
            },
            success: update_callback,
            error: update_callback,
            xhrFields: {
                withCredentials: true
            }
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
            update_dialog(request, user, identifier, formatter);
        }

        if (json.error) {
            pageObj.error++;
            if (pageObj.error > 3) {
                workflow_status.html('<span class=\"alert\">The job engine has failed.</span>');
                load_failed();
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

        if (json.tasks) {
            var jobs = json.tasks;
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
            var total = json.tasks.reduce(function(a, b) {
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
            fname:       'get_progress_log',
            workflow_id: job_id,
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
    window.location.href = "GenomeInfo.pl?gid=" + GENOME_ID;
}

function set_annotation_table() {
    $('#genome_annotation_table').tablesorter({widgets: ['zebra']});
}

function get_annotations() {
    $.ajax({
        data: {
            fname: 'get_annotations',
            gid: GENOME_ID,
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
            gid: GENOME_ID,
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
            gid: GENOME_ID,
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
            gid: GENOME_ID,
            dsid: dsid,
        },
        success : function() {
            get_datasets();
        },
    });
}

function open_aa_usage_table(chromosome) {
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
            dsgid: GENOME_ID,
            chr: chromosome
        },
        dataType: "html",
        success: function(html) {
            $('#aa_usage_table').html(html);
            $('#aa_table').tablesorter();
        }
    });
    
    return false;
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
            dsgid: GENOME_ID,
            chr: chromosome,
        },
        dataType: "html",
        success: function(html) {
            elem.html(html).slideDown();
        }
    })
}

function get_experiments(e) {
    var experiments = $("#experiments");
    experiments.html("Loading experiments... ").append(spinner.clone());

    $.ajax({
        data: {
            fname: 'get_experiments',
            gid: GENOME_ID
        },
        success:function(html) {
            experiments
            	.hide()
            	.html(html)
            	.slideDown();
            var count = experiments.find('span').length;
            $('#exp_count').html('('+count+')').show();
        }
    });
}

function update_owner() {
    var user_name = $('#edit_user').val();
    if (!user_name) {
        alert('Please specify a user.');
        return;
    }

    $.ajax({
        data: {
            fname: 'update_owner',
            gid: GENOME_ID,
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

function update_certified(val) {
    $.ajax({
        data: {
            fname: 'update_certified',
            gid: GENOME_ID,
            certified: val ? 1 : 0,
            timestamp: new Date().getTime()
        },
        success : function(data) {
            if (data) { // error
                alert(data);
            }
            else { // success
	            if (val)
	            	$('#certified_box').show();
	            else
	            	$('#certified_box').hide();
            }
        }
    });
}

function toggle_favorite(img) {
	$.ajax({
		data: {
			fname: 'toggle_favorite',
			gid: GENOME_ID,
		},
		success :  function(val) {
			$(img).attr({ src: (val == '0' ? "picts/star-hollow.png" : "picts/star-full.png") });
		}
	});
}