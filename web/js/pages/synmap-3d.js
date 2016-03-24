/*global document,$,alert,window,JSON */

var concat = Array.prototype.concat;

// Global experiment data
var current_experiment = {};
var options_name;
var final_experiment;

function search_genomes (search_term) {
    var edit_genome = $(geneSelect[0]);
	edit_genome.autocomplete("close");
    var spinner = $(geneSelect[1]);
	spinner.show();

	coge.services.search_genomes(search_term, { fast: true })
		.done(function(response) { // success
			if (response && response.genomes) {
				var results = response.genomes.map(function(obj) {
					var label = obj.info.replace(/&reg;/g, "\u00ae"); // (R) symbol
					return { label: label, value: obj.id };
				});
				edit_genome
					.autocomplete({source: results})
					.autocomplete("search");
			}
			spinner.hide();
		})
		.fail(function() { // error
			//TODO
		});
}

function render_template(template, container) {
    container.empty()
    .hide()
    .append(template)
    .show();//.slideDown();
}

function ExperimentDescriptionView(opts) {
    this.experiment = opts.experiment;
    this.metadata = opts.metadata;
    this.x_gid = opts.x_gid;
    this.y_gid = opts.y_gid;
    this.z_gid = opts.z_gid;
    this.onError = opts.onError;
    this.title = "Select Genomes";
    this.initialize();
}

$.extend(ExperimentDescriptionView.prototype, {
    initialize: function() {
        this.el = $($("#genomes-template").html());
        this.edit_xgenome = this.el.find("#edit_xgenome");
		this.edit_ygenome = this.el.find("#edit_ygenome");
		this.edit_zgenome = this.el.find("#edit_zgenome");
        if (this.metadata) {
            //this.el.find('#edit_name').val(this.metadata.name);
            //this.el.find('#edit_description').val(this.metadata.description);
            //this.el.find('#edit_version').val(this.metadata.version);
            //this.edit_source.val(this.metadata.source_name);

            //if (!this.metadata.restricted)
            //    this.el.find('#restricted').removeAttr('checked');

            //this.el.find('#edit_genome').val(this.metadata.genome); AKB REMOVED 9/15
        }
    },

    render: function() {
        var self = this;
        
        var edit_xgenome = this.edit_xgenome;
		var edit_ygenome = this.edit_ygenome;
		var edit_zgenome = this.edit_zgenome;

        edit_xgenome.unbind().change(function() {
            // Reset gid when item has changed
            self.x_gid = undefined;
        });
        edit_ygenome.unbind().change(function() {
            // Reset gid when item has changed
            self.y_gid = undefined;
        });
        edit_zgenome.unbind().change(function() {
            // Reset gid when item has changed
            self.z_gid = undefined;
        });

        edit_xgenome.autocomplete({
            source:[],
            select: function(event, ui) {
                $(this).val(ui.item.label);
                self.x_gid = ui.item.value;
                return false; // Prevent the widget from inserting the value.
            },

            focus: function(event, ui) {
                //$("#edit_xgenome").val(ui.item.label);
                return false; // Prevent the widget from inserting the value.
            }
        });
        edit_ygenome.autocomplete({
            source:[],
            select: function(event, ui) {
                $(this).val(ui.item.label);
                self.y_gid = ui.item.value;
                return false; // Prevent the widget from inserting the value.
            },

            focus: function(event, ui) {
                //$("#edit_ygenome").val(ui.item.label);
                return false; // Prevent the widget from inserting the value.
            }
        });
        edit_zgenome.autocomplete({
            source:[],
            select: function(event, ui) {
                $(this).val(ui.item.label);
                self.z_gid = ui.item.value;
                return false; // Prevent the widget from inserting the value.
            },

            focus: function(event, ui) {
                //$("#edit_zgenome").val(ui.item.label);
                return false; // Prevent the widget from inserting the value.
            }
        });


        edit_xgenome.keyup(function() {
            geneSelect = [ "#edit_xgenome", "#edit_xgenome_busy"];
            coge.utils.wait_to_search(search_genomes, self.edit_xgenome.get(0));
        });
        edit_xgenome.click(function() {
        	$(this).autocomplete('search');
        });
        edit_ygenome.keyup(function() {
            geneSelect = ["#edit_ygenome", "#edit_ygenome_busy"];
            coge.utils.wait_to_search(search_genomes, self.edit_ygenome.get(0));
        });
        edit_ygenome.click(function() {
			$(this).autocomplete('search');
        });
        edit_zgenome.keyup(function() {
            geneSelect = [ "#edit_zgenome", "#edit_zgenome_busy"];
            coge.utils.wait_to_search(search_genomes, self.edit_zgenome.get(0));
        });
        edit_zgenome.click(function() {
			$(this).autocomplete('search');
        });
    },

    is_valid: function() {
        var xgenome = this.el.find('#edit_xgenome').val();
		var ygenome = this.el.find('#edit_ygenome').val();
		var zgenome = this.el.find('#edit_zgenome').val();
        
		if (!xgenome || xgenome === 'Search' || !this.x_gid) {
        	if (this.onError)
            	this.onError('Please specify an X-axis genome.');
            return false;
        }

		if (!ygenome || ygenome === 'Search' || !this.y_gid) {
        	if (this.onError)
            	this.onError('Please specify a Y-axis genome.');
            return false;
        }

		if (!zgenome || zgenome === 'Search' || !this.z_gid) {
        	if (this.onError)
            	this.onError('Please specify a Z-axis genome.');
            return false;
        }

       $.extend(this.experiment, {
            metadata: {
				xgenome: xgenome,
				ygenome: ygenome,
				zgenome: zgenome
            },

            x_gid: this.x_gid,
	    	y_gid: this.y_gid,
	    	z_gid: this.z_gid
        });
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

        this.sortby = this.el.find("#sortby");

        this.min_synteny = this.el.find("#min_synteny");
        this.min_length = this.el.find("#min_length");

        this.cluster = this.el.find("#enable_cluster");
        this.cluster_container = this.el.find("#cluster-container");

        this.ratio = this.el.find("#enable_ratio");
        this.ratio_container = this.el.find("#ratio-container");

        if (this.sort) {
            this.sortby.val(this.sort)
        }
        if (this.min_syn) {
            this.min_synteny.val(this.min_syn)
        }
        if (this.min_len) {
            this.min_length.val(this.min_len)
        }

        if (this.c_eps && this.c_min) {
            this.cluster.prop("checked", true);
            this.cluster_container.slideDown();
            this.el.find("#cluster_eps").val(this.c_eps);
            this.el.find("#cluster_min").val(this.c_min);
        }

        if (this.ratio_type && this.r_by && this.r_min && this.r_max) {
            this.ratio.prop("checked", true);
            this.ratio_container.slideDown();
            this.el.find("#ratio").val(this.ratio_type);
            this.el.find("#ratio_by").val(this.r_by);
            this.el.find("#ratio_min").val(this.r_min);
            this.el.find("#ratio_max").val(this.r_max);
        }

    },

    toggleCluster: function() {
        this.cluster_enabled = this.cluster.is(":checked");

        if (this.cluster_enabled)
            this.cluster_container.slideDown();
        else
            this.cluster_container.slideUp();
    },

    toggleRatio: function() {
        this.ratio_enabled = this.ratio.is(":checked");

        if (this.ratio_enabled)
            this.ratio_container.slideDown();
        else
            this.ratio_container.slideUp();
    },

    render: function() {
        //var self = this;

        // jQuery Events
        this.el.find("#enable_cluster").unbind().change(this.toggleCluster.bind(this));
        this.el.find("#enable_ratio").unbind().change(this.toggleRatio.bind(this));
    },

    is_valid: function() {
        this.data = {};

        // Sort By
        var sortby = this.sortby.val();
        this.data.sortby = sortby;

        // Minimum synteny per contig.
        var min_synteny = this.min_synteny.val();
        if (min_synteny < 0) {
	    	if (this.onError)
	    	this.onError('Minimum chromosome length must be greater than 0.');
	    	return false;
		}
        this.data.min_synteny = min_synteny;

        // Minimum contig length.
		var min_length = this.min_length.val();
		if (min_length < 0) {
	    	if (this.onError)
	    	this.onError('Minimum chromosome length must be greater than 0.');
	    	return false;
		}
        this.data.min_length = min_length;

        // Cluster parsing.
        if (this.cluster.is(":checked")) {
            this.data.cluster = true;
            this.data.c_eps = this.el.find("#cluster_eps").val();
            this.data.c_min = this.el.find("#cluster_min").val();
        } else {
            this.data.cluster = false;
        }

        // Mutation ratio parsing.
        if (this.ratio.is(":checked")) {
            this.data.ratio = this.el.find("#ratio").val();
            this.data.r_by = this.el.find("#ratio_by").val();
            this.data.r_min = this.el.find("#ratio_min").val();
            this.data.r_max = this.el.find("#ratio_max").val();
        } else {
            this.data.ratio = false;
        }

		return true;
    },

    get_options: function() {
        return this.data;
    }
});

function AdvancedOptionView(opts) {
    this.data = {};
    this.vr = opts.vr;
    this.initialize();
}

$.extend(AdvancedOptionView.prototype, {
    initialize: function() {
	    this.el = $($("#advanced-options-template").html());
        this.enable_vr = this.el.find("#vr");

        if (this.vr) {
            this.enable_vr.prop("checked", true);
        }
    },

    is_valid: function() {
		this.data.vr = this.enable_vr.is(":checked");
		return true;
    },

    get_options: function() {
		return this.data;
    }
});

function AdminOptionsView() {
    this.data = {};
    this.initialize();
}

$.extend(AdminOptionsView.prototype, {
    initialize: function() {
        this.el = $($("#admin-options-template").html());
    },

    is_valid: function() {
        return true;
    },

    get_options: function() {
        return this.data;
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
	
		this.advanced_view = new AdvancedOptionView({
            vr: this.vr
        });
    	this.layout_view = new LayoutView({
        	template: "#options-layout-template",
			layout: {
        		"#general-options": this.general_view,
				"#advanced-options": this.advanced_view
        	}
    	});

		//if (this.admin)
		//    this.layout_view.updateLayout({"#admin-options": this.admin_view});

		this.el = this.layout_view.el;
    },

    // Validate and add all options to the experiment
    is_valid: function() {
        this.experiment.options = {};

        if (!this.advanced_view.is_valid() || !this.general_view.is_valid())
            return false;

        var options = $.extend({}, this.general_view.get_options(), this.advanced_view.get_options());

		/*if (this.admin) {
            if (!this.admin_view.is_valid())
                return false;
            $.extend(options, this.admin_view.get_options());
        }*/

        $.extend(this.experiment.options, options);

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

var visVisible = false;
function showVisualizer(data) {
    console.log(data);

    if (visVisible) {
        //Refresh function can go here.
    } else {
        console.log(final_experiment);
        $.getScript( "js/syn3d/syn3d.js", function( data, textStatus, jqxhr ) {
            console.log( "Visualizer loaded." );
            //console.log( data ); // Data returned
            //console.log( textStatus ); // Success
            //console.log( jqxhr.status ); // 200
        });
        $('#analysis').css("display", "");
    }

    $(document).ready( function() {
        /* Halt scrolling when controlling visualization. */
        $('#canvas').hover( function() {
            $("body").css("overflow", "hidden");
        }, function() {
            $("body").css("overflow", "scroll");
        });

        /* Change background color with slider */
        var renderDiv = $("#rendering");
        var bg_slider = $("#bgslider");
        renderDiv.css("background-color", getBgColor(bg_slider.val()));
        bg_slider.change( function() {
            renderDiv.css("background-color", getBgColor(bg_slider.val()));
        });
    });
}

function launch(experiment) {
    var xgid = experiment.x_gid;
    var ygid = experiment.y_gid;
    var zgid = experiment.z_gid;
    final_experiment = experiment;
    final_experiment.links = {};

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

    // Build Link to SynMap Output
    function synmapOutputLink(id1, id2) {
        var fileDir = "/home/asherkhb/repos/coge/web/data/diags/";
        var fileTag = ".CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.Dm0.ma1.gcoords.ks";
        var link = '';
        id1 = id1.toString();
        id2 = id2.toString();
        if ( parseInt(id1.charAt(0)) <= parseInt(id2.charAt(0)) ) {
            link = fileDir + id1 + '/' + id2 + '/' + id1 + '_' + id2 + fileTag;
        }
        else {
            link = fileDir + id2 + '/' + id1 + '/' + id2 + '_' + id1 + fileTag;
        }
        return link;
    }
    
    var xy_request = {
        type: 'synmap',
        requester: {
		    page:      PAGE_NAME,
		    user_name: USER_NAME
	    },
        parameters: {
            genome_id1: xgid,
            genome_id2: ygid,
            ks_type: "kn_ks"
        }
    };
    var xz_request = {
    	type: 'synmap',
        requester: {
		    page:      PAGE_NAME,
		    user_name: USER_NAME
	    },
		parameters: {
			genome_id1: xgid,
			genome_id2: zgid,
			ks_type: 'kn_ks'
		}
    };
    var yz_request = {
    	type: 'synmap',
        requester: {
		    page:      PAGE_NAME,
		    user_name: USER_NAME
	    },
		parameters: {
			genome_id1: ygid,
			genome_id2: zgid,
			ks_type: 'kn_ks'
		}
    };

    // Compute XY Comparison!
    // onSuccess: Launch XZ Comparison
    function makeXY() {
        coge.progress.init({title: "Running SynMap Analysis 1/3",
                            onSuccess: function(results) {
                                coge.progress.end();
                                final_experiment.links.xy = synmapOutputLink(xgid, ygid);
                                makeXZ();
                            } 
                           });
        coge.progress.begin();

        coge.services.submit_job(xy_request)
            .done(function(response) {
                console.log(response.id);
                if (!response) {
                    coge.progress.failed("Error: empty response from server");
                    return;
                }
                else if (!response.success || !response.id) {
                    coge.progress.failed("Error: failed to start workflow", response.error);
                    return;
                }
                //Start status update
                window.history.pushState({}, "Title", "SynMap3D.pl" + "?wid=" + response.id); // Add workflow id to URL
                coge.progress.update(response.id);
                //coge.progress.update(response.id, response.site_url);
            })
            .fail(function(jqXHR, textStatus, errorThrown) {
                console.log("Fucked!");
                coge.progress.failed("Couldn't talk to the server: " + textStatus + ': ' + errorThrown);
            })
    }
    
    // Compute XZ Comparison
    // onSuccess: Launch YZ Comparison
    function makeXZ() { 
        coge.progress.init({title: "Running SynMap Analysis 2/3",
                            onSuccess: function(results) {
                                coge.progress.end();
                                final_experiment.links.xz = synmapOutputLink(xgid, zgid);
                                makeYZ();
                            } 
                           });
        coge.progress.begin();
        
        coge.services.submit_job(xz_request)
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
                window.history.pushState({}, "Title", "SynMap3D.pl" + "?wid=" + response.id); // Add workflow id to URL
                coge.progress.update(response.id);
                //coge.progress.update(response.id, response.site_url);
            })
            .fail(function(jqXHR, textStatus, errorThrown) {
                coge.progress.failed("Couldn't talk to the server: " + textStatus + ': ' + errorThrown);
            })
    }

    // Compute YZ Comparison
    // onSuccess: Analyses are complete 
    function makeYZ() {
        coge.progress.init({title: "Running SynMap Analysis 3/3",
                            onSuccess: function(results) {
                                coge.progress.end();
                                final_experiment.links.yz = synmapOutputLink(ygid, zgid);
                                runDotplotDots(final_experiment);
                            } 
                           });
        coge.progress.begin();
        
        coge.services.submit_job(yz_request)
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
                window.history.pushState({}, "Title", "SynMap3D.pl" + "?wid=" + response.id); // Add workflow id to URL
                coge.progress.update(response.id);
                //coge.progress.update(response.id, response.site_url);
            })
            .fail(function(jqXHR, textStatus, errorThrown) {
                coge.progress.failed("Couldn't talk to the server: " + textStatus + ': ' + errorThrown);
            })
    }

    function runDotplotDots(exp) {
        coge.progress.init({title: "Finding Syntenic Points",
                            onSuccess: function(results) {
                                showVisualizer(final_experiment);
                            }
                           });
        coge.progress.begin();

        $.ajax({
                dataType: "json",
                data: {
                    fname: "dotplot_dots",
                    genome_idX: exp.x_gid, // xgid.toString(),
                    genome_idY: exp.y_gid, // ygid.toString(),
                    genome_idZ: exp.z_gid, // zgid.toString(),
                    ksfile_xy: exp.links.xy,
                    ksfile_xz: exp.links.xz,
                    ksfile_yz: exp.links.yz,
                    sort: exp.options.sortby,
                    min_length: exp.options.min_length,
                    min_synteny: exp.options.min_synteny,
                    cluster: exp.options.cluster,
                    c_eps: exp.options.c_eps,
                    c_min: exp.options.c_min,
                    ratio: exp.options.ratio,
                    r_by: exp.options.r_by,
                    r_min: exp.options.r_min,
                    r_max: exp.options.r_max,
                    graph_out: graph_obj,
                    log_out: log_obj
                }
            })
            .done(function (response) {
                if (!response) {
                    coge.progress.failed("Error: empty response from server");
                    return;
                }
                else if (!response.success || !response.id) {
                    coge.progress.failed("Error: failed to start workflow", response.error);
                    return;
                }

                //Start status update
                window.history.pushState({}, "Title", "SynMap3D.pl" + "?wid=" + response.id); // Add workflow id to URL
                coge.progress.update(response.id);
                //coge.progress.update(response.id, response.site_url); TODO: Add tiny link (see SynMap.pm for example)

                final_experiment.links.xy_json = response.xy_json;
                final_experiment.links.xz_json = response.xz_json;
                final_experiment.links.yz_json = response.yz_json;
                final_experiment.links.graph = response.graph;
                final_experiment.links.log = response.log;
            })
            .fail(function (jqXHR, textStatus, errorThrown) {
                coge.progress.failed("Couldn't talk to the server: " + textStatus + ': ' + errorThrown);
            })

    }

    //function runDotplotDots() {
    //    coge.progress.init({title: "Finding Syntenic Points",
    //                        onSuccess: function(results) {
    //                            showVisualizer(final_experiment);
    //                        }
    //                       });
    //    coge.progress.begin();
    //
    //    $.ajax({
    //            dataType: "json",
    //            data: {
    //                fname: "dotplot_dots",
    //                genome_idX: xgid.toString(),
    //                genome_idY: ygid.toString(),
    //                genome_idZ: zgid.toString(),
    //                ksfile_xy: final_experiment.links.xy,
    //                ksfile_xz: final_experiment.links.xz,
    //                ksfile_yz: final_experiment.links.yz,
    //                option_name: options_name,
    //                hide: experiment.options.hide,
    //                min_len: experiment.options.min_len
    //            }
    //        })
    //        .done(function (response) {
    //            if (!response) {
    //                coge.progress.failed("Error: empty response from server");
    //                return;
    //            }
    //            else if (!response.success || !response.id) {
    //                coge.progress.failed("Error: failed to start workflow", response.error);
    //                return;
    //            }
    //
    //            //Start status update
    //            window.history.pushState({}, "Title", "SynMap3D.pl" + "?wid=" + response.id); // Add workflow id to URL
    //            coge.progress.update(response.id);
    //            //coge.progress.update(response.id, response.site_url); TODO: Add tiny link (see SynMap.pm for example)
    //
    //            final_experiment.links.xy_json = response.xy_json;
    //            final_experiment.links.xz_json = response.xz_json;
    //            final_experiment.links.yz_json = response.yz_json;
    //            final_experiment.links.merge = response.merge_json;
    //            final_experiment.links.histo = response.histo_json;
    //        })
    //        .fail(function (jqXHR, textStatus, errorThrown) {
    //            coge.progress.failed("Couldn't talk to the server: " + textStatus + ': ' + errorThrown);
    //        })
    //
    //}

    // Start Comparisons
    makeXY();
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

    wizard.addStep(new ExperimentDescriptionView({
        experiment: current_experiment,
        metadata: opts.metadata,
        x_gid: opts.x_gid,
        y_gid: opts.y_gid,
        z_gid: opts.z_gid,
        onError: wizard.error_help.bind(wizard)
    }));

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
    return wizard;
}
