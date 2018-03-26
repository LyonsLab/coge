/*global document,$,alert,window,JSON */

var concat = Array.prototype.concat;

// Global experiment data
var current_experiment = {};
var options_name;
var final_experiment;
var synmapRenderer = "syn3d-1.0.js";

function search_genomes (search_term) {
    var edit_genome = $(geneSelect[0]);
	edit_genome.autocomplete("close");
    var spinner = $(geneSelect[1]);
	spinner.show();

	coge.services.search_genomes(search_term, { fast: true })
		.done(function(response) { // success
			if (response && response.genomes) {
				var results = response.genomes.map(function(obj) {
					var label = obj.info.replace(/&#x1f512;/g, "\uD83D\uDD12"); // Lock symbol
					//TODO add certified and favorite icons
					return { label: label, value: obj.id };
				});
				edit_genome
					.autocomplete({source: results})
					.autocomplete("search");
			}
			spinner.hide();
		})
		.fail(function() { // error
			spinner.hide();
		});
}

function render_template(template, container) {
    container.empty()
    .hide()
    .append(template)
    .show();//.slideDown();
}

function GenomesView(opts) {
    this.experiment = opts.experiment;
    this.metadata = opts.metadata;
    this.x_gid = opts.x_gid;
    this.y_gid = opts.y_gid;
    this.z_gid = opts.z_gid;
    this.validity = {'x': false, 'y': false, 'z': false};
    this.onError = opts.onError;
    this.title = "Select Genomes";
    this.initialize();
}

$.extend(GenomesView.prototype, {
    initialize: function() {
        this.el = $($("#genomes-template").html());
        this.edit_genome = this.el.find("#edit_genome");
    },

    checkGenome: function(info) {
        // Checks a genome for chromosomes, length > 0, and CDS.
        var self = this;
        if (info.chromosome_count > 0) {
            var len = 0;
            var cds = 0;
            var i = 0;
            var c = info["chromosomes"].length;

            while (i<c && (len < 1 && cds < 1)) {
                len += parseInt(info["chromosomes"][i]["length"]);
                cds += parseInt(info["chromosomes"][i]["CDS_count"]);
                i+=1;
            }
            if (len > 0 && cds > 0) {
                return true
            } else {
                if (self.onError)
                    self.onError('Genome ' + info["organism"]["name"] + ' (id' + info["id"] + ') has zero length or no CDS. Please select a different genome.');
                return false
            }
        } else{
            if (self.onError)
                    self.onError('Genome ' + info["organism"]["name"] + ' (id' + info["id"] + ') has zero length or no CDS. Please select a different genome.');
            return false;
        }
    },

    validateGenome: function(name, gid, wait_indicator_id) {
        // Validates genomes and handles visual elements.
        var self = this;
        var wait = $('#' + wait_indicator_id);
        wait.show();
        var info = $.getJSON(API_BASE_URL + 'genomes/' + gid);
        $.when(info).done(function(data, textStatus, jqXHR) {
            var valid = self.checkGenome(data);
            if (valid) {
                var g = document.getElementById('genomes');
                var div = document.createElement('div');
                div.innerHTML = name;
                div.gid = gid;
                g.appendChild(div);
                var img = document.createElement('img');
                img.src = 'picts/delete.png';
                img.style.cursor = 'pointer';
                img.onclick = function() {var p=this.parentNode;p.removeChild(this.previousElementSibling);p.removeChild(this);}
                g.appendChild(img);
            }
            self.edit_genome.val(null);
            wait.hide();
            return valid;
        });
    },

    validate_on_load: function() {
        // This function is necessary to permit validation after rendering of template components.
        var self = this;
        // if (self.x_gid) {
        //     self.validateGenome('x', self.x_gid, 'edit_xgenome_busy', 'x_status')
        // }
    },

    render: function() {
        var self = this;
        this.edit_genome.autocomplete({
            source:[],
            select: function(event, ui) {
                $(this).val(ui.item.label);
                self.validateGenome(ui.item.label, ui.item.value, 'edit_genome_busy');
                return false; // Prevent the widget from inserting the value.
            },

            focus: function(event, ui) {
                //$("#edit_xgenome").val(ui.item.label);
                return false; // Prevent the widget from inserting the value.
            }
        });

        this.edit_genome.keyup(function() {
            geneSelect = [ "#edit_genome", "#edit_genome_busy"];
            coge.utils.wait_to_search(search_genomes, self.edit_genome.get(0));
        });
        this.edit_genome.click(function() {
        	$(this).autocomplete('search');
        });
    },

    is_valid: function() {
        if (document.getElementById('genomes').childElementCount < 4) {
            if (this.onError)
                this.onError('Please add two or more genomes before continuing.');
            return false;
        }
        return true;
    }
});

function GeneralOptionsView(opts) {
    this.data = {};
    this.sort = opts.sortby;
    this.min_syn = opts.min_syn;
    this.min_len = opts.min_len;
    this.c_eps = opts.c_eps;
    this.c_min = opts.c_min;
    this.ratio_type = opts.ratio;
    this.r_by = opts.r_by;
    this.r_min = opts.r_min;
    this.r_max = opts.r_max;
    this.initialize();
    this.onError = opts.onError;
}

$.extend(GeneralOptionsView.prototype, {
    initialize: function() {
        this.el = $($("#general-options-template").html());
    },

    render: function() {
    },

    is_valid: function() {
		return true;
    },

    get_options: function() {
        return {
            tdd: $('#tdd').val(),
            D: $('#D').val(),
            A: $('#A').val(),
            gm: $('#gm').val(),
            Dm: $('#Dm').val(),
            blast: $('#blast').val(),
            blast_option: $('#blast').val() == 6 ? null : $('#blast_option').children()[1].value,
            feat_type1: $('#feat_type1').val(),
            feat_type2: $('#feat_type2').val(),
            dsgid1: $('#dsgid1').val(),
            dsgid2: $('#dsgid2').val(),
            basename: pageObj.basename,
            dagchainer_type: $('#dagchainer_type:checked').val(),
            ks_type: $('#ks_type').val(),
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
            color_scheme: $('#color_scheme').val(),
            codeml_min: $('#codeml_min').val(),
            codeml_max: $('#codeml_max').val(),
            logks: $('#logks')[0].checked,
            csco: $('#csco').val()
        };
    }
});

function OptionsView(opts) {
    this.experiment = opts.experiment;
    this.admin = opts.admin;
    this.sortby = opts.sortby;
    this.min_syn = opts.min_syn;
    this.min_len = opts.min_len;
    this.c_eps = opts.c_eps;
    this.c_min = opts.c_min;
    this.ratio = opts.ratio;
    this.r_by = opts.r_by;
    this.r_min = opts.r_min;
    this.r_max = opts.r_max;
    this.vr = opts.vr;
    this.onError = opts.onError;
    this.title = "Options";
    this.initialize();
}

//AKB - Removed Admin Options (commented out in this extend)
$.extend(OptionsView.prototype, {
    initialize: function() {
        //this.admin_view = new AdminOptionsView();
        this.general_view = new GeneralOptionsView({
			onError: this.onError,
            sortby: this.sortby,
            min_syn: this.min_syn,
            min_len: this.min_len,
            c_eps: this.c_eps,
            c_min: this.c_min,
            ratio: this.ratio,
            r_by: this.r_by,
            r_min: this.r_min,
            r_max: this.r_max
		});
    	this.layout_view = new LayoutView({
        	template: "#options-layout-template",
			layout: {
        		"#general-options": this.general_view
        	}
    	});

		//if (this.admin)
		//    this.layout_view.updateLayout({"#admin-options": this.admin_view});

		this.el = this.layout_view.el;
    },

    // Validate and add all options to the experiment
    is_valid: function() {
        this.experiment.options = {};

        if (!this.general_view.is_valid())
            return false;

        $.extend(this.experiment.options, this.general_view.get_options());

        return true;
    },

    render: function() {
        // Render the views added to the layout view
        this.layout_view.renderLayout();
    }
});

function getBgColor(slideValue) {
    var colorCodes = ["#f5f5f5", "#dcdcdc", "#c4c4c4", "#ababab", "#939393",
                      "#7a7a7a", "#626262", "#494949", "#313131", "#181818"];
    return colorCodes[slideValue-1];
}

function showVisualizer(data) {
    $('#chartSvg').remove(); //removes histograms if a old version exists.
    $.getScript( "js/syn3d/" + synmapRenderer, function( data, textStatus, jqxhr ) {
        console.log( "Visualizer loaded." );
        //console.log( data ); // Data returned
        //console.log( textStatus ); // Success
        //console.log( jqxhr.status ); // 200
    }, false);
}

function launch(experiment) {
    var xgid = experiment.x_gid;
    var ygid = experiment.y_gid;
    var zgid = experiment.z_gid;
    final_experiment = experiment;
    //final_experiment.links = {};

    // Build Options-Based Name
    function buildOptionsName(exp) {
        var name = xgid + '_' + ygid + '_' + zgid + '_';
        var name_ext = [];
        // Add sort method.
        name_ext.push("sort=" + exp.options.sortby);
        // Add length parsing.
        name_ext.push('parse.len=' + exp.options.min_length);
        // Add synteny parsing.
        name_ext.push('parse.syn=' + exp.options.min_synteny);
        // Add clustering options.
        if (exp.options.cluster) {
            name_ext.push('cluster.eps=' + exp.options.c_eps);
            name_ext.push('cluster.min=' + exp.options.c_min);
        }
        // Add ratio options.
        if (exp.options.ratio) {
            name_ext.push('ratio.by=' + exp.options.r_by + '.' + exp.options.ratio);
            name_ext.push('ratio.min=' + exp.options.r_min);
            name_ext.push('ratio.max=' + exp.options.r_max);
        }
        // Sort & join name extensions, add to GID string.
        name_ext=name_ext.sort().join("_");
        name = name + name_ext;

        return name;
    }
    options_name = buildOptionsName(final_experiment);
    var graph_obj = options_name + '_graph.json';
    var log_obj = options_name + '_log.json';
    var download_obj = options_name + '_data.txt';
    final_experiment.download = download_obj;
    final_experiment.graph = graph_obj;
    final_experiment.log = log_obj;

    // Build URL for updating.
    function buildUrl(exp) {
        var url = "?x_gid=" + exp.x_gid + ";y_gid=" + exp.y_gid + ";z_gid=" + exp.z_gid;
        url += ";min_syn=" + exp.options.min_synteny + ";min_len=" + exp.options.min_length;
        url += ";sort=" + exp.options.sortby;
        if (exp.options.cluster) {
            url += ";cluster=" + exp.options.c_eps + "," + exp.options.c_min;
        }
        if (exp.options.ratio) {
            url += ";ratio=" + exp.options.ratio + "," + exp.options.r_by + "," + exp.options.r_min + "," + exp.options.r_max;
        }
        if (exp.options.vr) {
            url += ";vr=" + exp.options.vr;
        }

        return url
    }
    function getTiny(url) {
        var request_url = "https://genomevolution.org/r/yourls-api.php?signature=d57f67d3d9&action=shorturl&format=simple&url=" + url;
        var link;
        // $.get(request_url, function( data ) {
        //     link = data;
        //     return link;
        // });
        $.ajax({
            url: request_url,
            success: function(data) {
                link = data;
            },
            async: false
        });
        return link;
    }
    var urlUpdate = buildUrl(final_experiment);
    final_experiment.page_url = SERVER_URL + PAGE_NAME + urlUpdate;
    final_experiment.tiny_url = getTiny(SERVER_URL + PAGE_NAME + urlUpdate);

    // Build Link to SynMap Output (AKB Removed 5/25/16 - no longer needed)
    // function synmapOutputLink(id1, id2) {
    //     var fileDir = "/storage/coge/data/diags/";
    //     var fileTag = ".CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.Dm0.ma1.gcoords.ks";
    //     var link = '';
    //     id1 = id1.toString();
    //     id2 = id2.toString();
    //     if ( parseInt(id1.charAt(0)) <= parseInt(id2.charAt(0)) ) {
    //         link = fileDir + id1 + '/' + id2 + '/' + id1 + '_' + id2 + fileTag;
    //     }
    //     else {
    //         link = fileDir + id2 + '/' + id1 + '/' + id2 + '_' + id1 + fileTag;
    //     }
    //     return link;
    // }

    var synmap_request = {
        type: 'synmap3d',
        requester: {
		    page: PAGE_NAME
	    },
        parameters: {
            genome_id1: xgid,
            genome_id2: ygid,
            genome_id3: zgid,
            ks_type: "kn_ks",
            tinylink: final_experiment.tiny_url,
            //ksfile_xy: synmapOutputLink(xgid, ygid),
            //ksfile_xz: synmapOutputLink(xgid, zgid),
            //ksfile_yz: synmapOutputLink(ygid, zgid),
            sort: final_experiment.options.sortby,
            min_length: final_experiment.options.min_length,
            min_synteny: final_experiment.options.min_synteny,
            cluster: final_experiment.options.cluster.toString(),
            c_eps: final_experiment.options.c_eps,
            c_min: final_experiment.options.c_min,
            ratio: final_experiment.options.ratio.toString(),
            r_by: final_experiment.options.r_by,
            r_min: final_experiment.options.r_min,
            r_max: final_experiment.options.r_max,
            graph_out: graph_obj,
            log_out: log_obj,
            download: download_obj
        }
    };

    function makeSynmaps() {
        console.log(final_experiment);
        coge.progress.init({title: "Running SynMap3D Analysis",
                            onSuccess: function(results) {
                                showVisualizer(final_experiment);
                                console.log(final_experiment);
                            }
        });
        coge.progress.begin();

        coge.services.submit_job(synmap_request)
            .done(function(response) {
                if (!response) {
                    coge.progress.failed("Error: empty response from server");
                    return;
                }
                else if (!response.success || !response.id) {
                    coge.progress.failed("Error: failed to start workflow", response.error);
                    return;
                }

                //Start status update
                console.log("Launching SynMap3D Job (ID: " + response.id + ")");
                window.history.pushState({}, "Title", "SynMap3D.pl" + urlUpdate + ";wid=" + response.id); // Update URL with all options.
                //final_experiment.page_url = SERVER_URL + PAGE_NAME + urlUpdate; // Add that URL to final_experiment.
                coge.progress.update(response.id);
            })
            .fail(function(jqXHR, textStatus, errorThrown) {
                coge.progress.failed("Couldn't talk to the server: " + textStatus + ': ' + errorThrown);
            })
    }

    // Check if results files already exist, if so launch visualizer, if not start workflow.
    // Note: this is incredibly hacky to get around automatic redirects that result in 200 even if file exists. Might
    // need to update at some point if CoGe: Error page changes.
    var graphLoc = DATA_LOC + "/" + graph_obj;
    $.get(graphLoc)
        .done(function(d) {
            var testEl = document.createElement('html');
            testEl.innerHTML = d;
            var title = testEl.getElementsByTagName('TITLE');
            if (title.length > 0) {
                if (testEl.getElementsByTagName('TITLE')[0].text == 'CoGe: Error') makeSynmaps();
            } else {
                window.history.pushState({}, "Title", "SynMap3D.pl" + urlUpdate); // Update URL with all options.
                showVisualizer(final_experiment);
            }
        })
        .fail(function() { makeSynmaps() });
}

function reset_launch() { //AKB - Renamed from reset_load()
    window.history.pushState({}, "Title", PAGE_NAME);

    // Reset the wizard and set description step
    initialize_wizard({
        admin: IS_ADMIN,
        metadata: current_experiment.metadata,
        x_gid: current_experiment.x_gid,
        y_gid: current_experiment.y_gid,
        z_gid: current_experiment.z_gid,
        sortby: current_experiment.sortby,
        min_syn: MIN_SYN,
        min_len: MIN_LEN,
        c_eps: C_EPS,
        c_min: C_MIN,
        ratio: RATIO,
        r_by: R_BY,
        r_min: R_MIN,
        r_max: R_MAX,
	    vr: VR
    });

    $('#wizard-container').hide().fadeIn();
}

// Log Experiment for Confirmation
function logExperiment(exp) {
    console.log(exp);
    console.log("Sort By: " + exp.options.sortby);
    console.log("Minimum Synteny: " + exp.options.min_synteny);
    console.log("Minimum Length: " + exp.options.min_length);
    console.log("Clustering: " + exp.options.cluster);
    console.log("--> EPS: " + exp.options.c_eps);
    console.log("--> Minimum Cluster: " + exp.options.c_min);
    console.log("Ratio Parsing: " + exp.options.ratio);
    console.log("--> By: " + exp.options.r_by);
    console.log("--> Min: " + exp.options.r_min);
    console.log("--> Max: " + exp.options.r_max);
    console.log("VR Mode: " + exp.options.vr);
}

function initialize_wizard(opts) {
    current_experiment = {};
    var root = $("#wizard-container");
    var wizard = new Wizard({ 
    	onCompleted: launch,  //logExperiment,  //AKB, changed from 'load'
    	data: current_experiment, 
    	helpUrl: opts.helpUrl 
    });

    var expView = new GenomesView({
        experiment: current_experiment,
        metadata: opts.metadata,
        x_gid: opts.x_gid,
        y_gid: opts.y_gid,
        z_gid: opts.z_gid,
        onError: wizard.error_help.bind(wizard)
    });
    wizard.addStep(expView);

    wizard.addStep(new OptionsView({
		experiment: current_experiment,
        sortby: opts.sortby,
        min_syn: opts.min_syn,
        min_len: opts.min_len,
        c_eps: opts.c_eps,
        c_min: opts.c_min,
        ratio: opts.ratio,
        r_by: opts.r_by,
        r_min: opts.r_min,
        r_max: opts.r_max,
		vr: opts.vr,
		//admin: opts.admin,
		onError: wizard.error_help.bind(wizard) 
    }));

    wizard.addStep(new ConfirmationView(current_experiment));

    // Render all the wizard sub views
    wizard.render();

    // Add the wizard to the document
    root.html(wizard.el);

    // Validate pre-loaded genomes after elements are all loaded
    expView.validate_on_load();

    return wizard;
}
