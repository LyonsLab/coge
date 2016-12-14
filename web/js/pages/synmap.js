//TODO: Create a proper helper function
function get_gc (dsgid, divid)
{
    $('#'+divid).removeClass('link').html('loading...');
    get_dsg_gc(['args__dsgid','args__'+dsgid,'args__text','args__1'],[divid]);
}

//TODO: Replace with proper promise chain
function get_organism_chain(type,val,i)
{
	$('#feattype_menu'+i).hide();
	$('#dsgid'+i).html('');
	
    $('#org_list').html('<input type=hidden id = "org_id"+i><font class="loading"></font>');
  //  if (type == 'name')
    	get_orgs(['args__search','args__'+val,'args__i','args__'+i], ['org_list'+i]);
//    else if (type == 'desc') 
//   	get_orgs(['args__desc','args__'+val,'args__i','args__'+i], ['org_list'+i]);
    //$('#dsg_info'+i).html('<div class="loading dna_small small">loading. . .</div>');
    $('#dsg_info'+i).html('<div class="small note indent">loading... <img src="picts/ajax-loader.gif"/></div>');
    ajax_wait("gen_dsg_menu(['args__oid','org_id"+i+"', 'args__num','args__"+i+"'],['dsg_menu"+i+"', 'genome_message"+i+"']);");
    ajax_wait("get_genome_info(['args__dsgid','dsgid"+i+"','args__org_num','args__"+i+"'],[handle_dsg_info]);");
}

//TODO: Replace with proper promise chain
function get_genome_info_chain(i) {
	//$('#dsg_info'+i).html('<div class=dna_small class=loading class=small>loading. . .</div>');
	$('#dsg_info'+i).html('<div class="small note indent">loading... <img src="picts/ajax-loader.gif"/></div>');
    // ajax_wait("gen_dsg_menu(['args__oid','org_id"+i+"', 'args__num','args__"+i+"'],['dsg_menu"+i+"','genome_message"+i+"']);");
    gen_dsg_menu(['args__oid','org_id'+i, 'args__num','args__'+i],['dsg_menu'+i, 'genome_message'+i]);
    $('#depth_org_1').html($('#org_id1 option:selected').html());
    $('#depth_org_2').html($('#org_id2 option:selected').html());

    ajax_wait("get_genome_info(['args__dsgid','dsgid"+i+"','args__org_num','args__"+i+"'],[handle_dsg_info]);");
}

function close_dialog(dialog) {
    dialog.dialog('close');
    dialog.find('#text').empty();
    dialog.find('#progress').show();
    dialog.find('#dialog_error').hide();
    dialog.find('#dialog_success').hide();
}

function load_results() {
    $('#intro').hide();
    $('#log_text').hide();
    $('#results').fadeIn();
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
    if (dsg_html)
    	$('#feattype_menu'+org_num).show();
    
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
            get_organism_chain('search',$('#'+val).val(),i[0])
        } else {
            pageObj.pageobjtime = setTimeout("get_organism_chain('search',$('#"+val+"').val(),i[0])",500);
        }
    } else if (descre.exec(val)) {
        if (val2) {
            get_organism_chain('search',$('#'+val).val(),i[0])
        } else {
            pageObj.pageobjtime = setTimeout("get_organism_chain('search',$('#"+val+"').val(),i[0])",200);
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

function display_legacy_settings() {  //AKB Added 2016-10-18

    if ($('#visualizer_select')[0].checked) {
        $('#legacy_opts').hide(0);
    } else {
        $('#legacy_opts').show(0);
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

function update_blast_option(val) {
    var l = $('#blast_option');
    if (val == 4) { // lastz
        var c = l.children();
        c[0].innerHTML = '--hspthresh';
        c[1].value = '3000';
        c[2].innerHTML = '(default 3000)';
        l.show();
    } else if (val == 6) // lastal
        l.hide();
    else { // blasts
        var c = l.children();
        c[0].innerHTML = '-evalue';
        c[1].value = '0.0001';
        c[2].innerHTML = '(default 0.0001)';
        l.show();
    }
}

var coge = window.coge = (function(namespace) {
    var ArrayProto = Array.prototype.slice;
    var slice = ArrayProto.slice;
    namespace.synmap = {

        setup: function(options) {
            var root = options.rootElement;

            $.ajaxSetup({
                type: "GET",
                url: options.page,
                dataType: "html",
                cache: false
            });

            $(".dialog_box").dialog({
                autoOpen: false,
                width: 500,
            });

            $("#synmap_dialog").dialog({modal: true});

            if($('#org_name1').val() != "Search") {
                $('#org_name1').css({fontStyle: "normal"});
                timing('org_name1',1);
            }
            if($('#org_desc1').val() != "Search") {
                $('#org_desc1').css({fontStyle: "normal"});
                timing('org_desc1',1);
            }
            if($('#org_name2').val() != "Search") {
                $('#org_name2').css({fontStyle: "normal"});
                timing('org_name2',1);
            }
            if($('#org_desc2').val() != "Search") {
                $('#org_desc2').css({fontStyle: "normal"});
                timing('org_desc2',1);
            }

            if ($('#assemble')[0].checked) {
                $('#assemble_info').toggle();
            }

            merge_select_check();
            depth_algo_check();

            $("#pair_info").draggable();
            $("#tabs").tabs({selected:0});
            $(".resizable").resizable();
            $('#depth_org_1').html($('#org_id1 option:selected').html());
            $('#depth_org_2').html($('#org_id2 option:selected').html());

            if (options.autostart)
                this.get_results(this.go.bind(this));

            $("#tabs").removeClass("invisible");

            // track analysis
            $("#synmap_go").on("click", function() {
                coge.synmap.run_synmap();
                ga('send', 'event', 'synmap', 'run');
            });
        },

        //TODO: Create a proper state object
        populate_page_obj: function(basefile) {
            if (!basefile) {
                basefile = "SynMap_" + Math.floor(Math.random() * 99999999 + 1);
            }
            pageObj.basename = basefile;
            pageObj.nolog = 0;
            pageObj.waittime = 1000;
            pageObj.runtime = 0;
            pageObj.fetch_error = 0;
            pageObj.error = 0;
            pageObj.engine = $("<span></span>", {
                "class": "alert",
                text: "The job engine has failed."
            });
            pageObj.failed = "The workflow could not be submitted.";
        },

        //TODO: Replace helper with state checking
        has_organisms: function() {
            return ($('#org_id1').val() != "") && ($('#org_id2').val() != "");
        },

        run_synmap: function() {
            var self = this;

            this.populate_page_obj();

            //var org_name1 = pageObj.org_name1;
            //var org_name2 = pageObj.org_name2;
            var feat_type1  = $('#feat_type1').val();
            var feat_type2  = $('#feat_type2').val();
            var org_length1 = $('#org_length1').html() || pageObj.org_length1;
            var org_length2 = $('#org_length2').html() || pageObj.org_length2;
            var seq_type1   = $('#seq_type1').html() || pageObj.seq_type1;
            var seq_type2   = $('#seq_type2').html() || pageObj.seq_type2;
            
            // Block this analysis if the genomes are unmasked, unannotated, and too large
            // feat_type 1 == CDS, 2 == genomic, seq_type == 1 is unmasked
            var max_size = 50 * 1000 * 1000;
            console.log('Block check: ' + org_length1 + ' ' + feat_type1 + ' ' + seq_type1 + ' ' + org_length2 + ' ' + feat_type2 + ' ' + seq_type2);
            if (( org_length1 > max_size && feat_type1 == 2 && seq_type1 == 1) &&
                ( org_length2 > max_size && feat_type2 == 2 && seq_type2 == 1))
            {
                var message = "Unfortunately this analysis cannot be performed. "
      			  + "A comparison of two unmasked and unannotated genomes of "
    			  + "these sizes requires many days to weeks to finish. "
    			  + "Please try:  1) select a hard-masked sequence or 2) use at least one annotated genome.";
                $('#log_text').hide();
                $('#results').html("<div class='alert'>Analysis Blocked:</div> " 
                		+ "<div class='text'>" + message + "</div>")
                		.show();
                window.scrollTo(0,0); // scroll to top of window to make error message visible
                return;
            }

            if (!this.has_organisms())
                return;

            if ($('#blast').val() == 5 && (feat_type1 != 1 || feat_type2 != 1) ) {
                alert('BlastP only works if both genomes have protein coding sequences (CDS) AND CDS is selected for both!');
                return;
            }
            
            if ($('#frac_bias')[0].checked && $('#depth_algo').val() == 0) {
            	alert('You can only run Fractination Bias if you select the Quota Align algorithm for Syntenic Depth.');
            	return;
            }

            $('#results').hide();

            // If user checked "Regenerate Images" then force the workflow to be re-executed.
            var regenerate = $('#regen_images')[0].checked;
            if (regenerate) {
                this.go();
                return;
            }

            // If results already exist for this analysis then show them without re-executing the workflow.
            this.get_results().done(function(data) {
                if (data.error) {
                    self.go(); // Results don't exist, re-execute the workflow
                }
            });
        },

        submit_assembly: function(e, input, gid1, gid2,flip) {
        	var self = this;
            e.preventDefault();

            var promise = $.ajax({
                dataType: "json",
                data: {
                    fname: "generate_assembly",
                    jquery_ajax: 1,
                    input: input,
                    gid1: gid1,
                    gid2: gid2,
                    flip: flip
                }
            });

            $("#dialog").dialog({
                autoOpen: true,
                position: {
                    my: "top",
                    at: "top+150",
                }
            });

            promise.then(function(response) { return self.wait_for_assembly.call(self, response); })
                   .then(function(url) { self.download_file.call(self, url); }, self.report_error);
        },

        get_job_status: function(id) {
            return $.getJSON("jex/status/" + id);
        },

        wait_for_assembly: function(response) {
            var deferred = $.Deferred();

            if (response && response.success) {
                this.wait_for_job(response.id, deferred, response.output);
            } else {
            	console.warn('synmap:wait_for_assembly: error response');
                deferred.reject(undefined);
            }

            return deferred;
        },

        wait_for_job: function(id, promise, args) {
            this.get_job_status(id).then(function(response) {
                switch(response.status) {
                    case "Completed":
                    	console.log(args);
                        promise.resolve(args);
                        break;
                    case "Failed":
                    	console.warn(response);
                        promise.reject("The workflow has failed");
                        break;
                    default:
                        setTimeout(function() {
                            coge.synmap.wait_for_job(id, promise, args);
                        }, 3000);
                        break;
                }
            });
        },

        report_error: function() {
            $("#dialog").html("<div class='padded alert'>Error: the pseudo-assembly could not be generated</div>").dialog("open");
        },

        download_file: function(url) {
        	console.log('synmap:download_file: ' + url);
            if (url) {
            	$("#dialog").dialog("close");
            	window.open(url, "_self");
            }
            else
            	this.report_error();
        },

        get_results: function(on_error) {
            var self = this;

            var spinner = $("#overlay").show();
            coge.progress.onReset = coge.progress.saveOnReset;

            var params = this.get_params();
            params['fname'] = 'get_results';
            params['jquery_ajax'] = 1;
            return $.ajax({
                type: 'GET',
                data: params,
                dataType: "json",
                success: function(data) {
                    $('.box').css("float", "left");
                    spinner.hide();
                    if (!data.error) {
                        $("#synmap_zoom_box").draggable();
                        $('#results')
                            .html(data.html) // insert rendered synmap_results.tmpl into page
                            .slideDown(show_results);
                    }
                    else if (on_error) {
                        on_error();
                    }
                },
                error: function() {
                    spinner.hide();
                }
            });
        },

        go: function() {
            coge.progress.begin();
            var params = this.get_params();
            params['genome_id1'] = params['dsgid1'];
            delete params['dsgid1'];
            params['genome_id2'] = params['dsgid2'];
            delete params['dsgid2'];
            var request = {
                type: 'synmap',
                requester: {
                    page: PAGE_NAME
                },
                parameters: params
            };
    
            coge.services.submit_job(request) 
                .done(function(response) {
                    if (!response) {
                        coge.progress.failed("Error: empty response from server");
                        return;
                    }
                    else if (!response.success || !response.id) {
                        coge.progress.failed("Error: failed to start workflow", response.error);
                        return;
                    }
                    
                    // Start status update
//                    window.history.pushState({}, "Title", "SynMap.pl" + "?wid=" + response.id); // Add workflow id to browser URL
                    coge.progress.update(response.id, response.site_url);
                    coge.progress.saveOnReset = coge.progress.onReset;
                    coge.progress.onReset = coge.synmap.get_results.bind(coge.synmap);
                })
                .fail(function(jqXHR, textStatus, errorThrown) {
                    coge.progress.failed("Couldn't talk to the server: " + textStatus + ': ' + errorThrown);
                });
        },

        get_params: function() {
            return {
                tdd: $('#tdd').val(),
                D: $('#D').val(),
                A: $('#A').val(),
                gm: $('#gm').val(),
                Dm: $('#Dm').val(),
                blast: $('#blast').val(),
                blast_option: $('#blast').val() == 6 ? null : $('blast_option').children()[1].value,
                feat_type1: $('#feat_type1').val(),
                feat_type2: $('#feat_type2').val(),
                dsgid1: $('#dsgid1').val(),
                dsgid2: $('#dsgid2').val(),
                jobtitle: $('#jobtitle').val(),
                basename: pageObj.basename,
                email: $('#email').val(),
                regen_images: $('#regen_images')[0].checked,
                width: $('#master_width').val(),
                dagchainer_type: $('#dagchainer_type:checked').val(),
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
                frac_bias: $('#frac_bias')[0].checked,
                fb_window_size: $('#fb_window_size').val(),
                fb_target_genes: $('#fb_target_genes')[0].checked,
                fb_numquerychr: $('#fb_numquerychr').val(),
                fb_numtargetchr: $('#fb_numtargetchr').val(),
                fb_remove_random_unknown: $('#fb_remove_random_unknown')[0].checked,
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
                vis: $('#visualizer_select:checked').val()  //AKB Added 2016-10-18
            };
        },

        onError: function() {
            var l = document.location.pathname;
            l = l.substr(0, l.indexOf('/coge/') + 6);
            var genome_id1 = $('#dsgid1').val();
            var genome_id2 = $('#dsgid2').val();
            if (genome_id2 < genome_id1) {
                var t = genome_id1;
                genome_id1 = genome_id2;
                genome_id2 = t;
            }
            var tiny_url = coge.progress.url;
            tiny_url = tiny_url.substr(tiny_url.lastIndexOf('/') + 1);
            $(".logfile a").attr("href", l + 'data/diags/' + genome_id1 + '/' + genome_id2 + '/' + tiny_url + '.log');
        }
    };
    return namespace;
})(coge || {});
