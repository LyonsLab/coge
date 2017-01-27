function get_experiment_info() {
    $.ajax({
        data: {
            fname: 'get_experiment_info',
            eid: EXPERIMENT_ID,
        },
        success : function (data) {
            $('#experiment_info').html(data);
        }
    });
}

function edit_experiment_info() {
    $.ajax({
        data: {
            fname: 'edit_experiment_info',
            eid: EXPERIMENT_ID,
        },
        success : function(data) {
            var obj = jQuery.parseJSON(data);
            $("#experiment_info_edit_box").html(obj.output).dialog('open').dialog({width: '32em'});
        }
    });
}

function update_experiment_info () {
    var name = $('#edit_name').val();
    if (!name) {
        alert('Please specify a name.');
        return;
    }

    var version = $('#edit_version').val();
    if (!version) {
        alert('Please specify a version');
        return;
    }

    var source_id = parseInt( $('#edit_source_id').val() );
    if (!source_id) {
        alert('Please specify a source');
        return;
    }

    var desc = $('#edit_desc').val();
    var link = $('#edit_link').val();

    $.ajax({
        data: {
            fname: 'update_experiment_info',
            eid: EXPERIMENT_ID,
            name: name,
            desc: desc,
            source_id: source_id,
            version: version,
            link: link
        },
        success : function(val) {
            get_experiment_info();
            $("#experiment_info_edit_box").dialog('close');
        }
    });
}

function get_sources () {
    $.ajax({
        data: {
            fname: 'get_sources'
        },
        success : function(data) {
            var items = jQuery.parseJSON(data);
            $("#edit_source").autocomplete("option", "source", items).autocomplete("search");
        }
    });
}

function make_experiment_public () {
    $.ajax({
        data: {
            fname: 'make_experiment_public',
            eid: EXPERIMENT_ID
        },
        success : function(val) {
            get_experiment_info();
        }
    });
}

function make_experiment_private () {
    $.ajax({
        data: {
            fname: 'make_experiment_private',
            eid: EXPERIMENT_ID
        },
        success : function(val) {
            get_experiment_info();
        }
    });
}

function add_tag_to_experiment () {
    var name = $('#edit_tag_name').val();

    if (name) {
        $.ajax({
            data: {
                fname: 'add_tag_to_experiment',
                eid: EXPERIMENT_ID,
                name: name
            },
            success : function(val) {
                get_experiment_info();
                $("#experiment_tag_edit_box").dialog('close');
            },
        });
    }
    else { alert ('Error: Must have type name specified!');}
}

function reset_location() {
	window.history.pushState({}, "Title", PAGE_NAME + "?eid=" + EXPERIMENT_ID);
}

function download_files() {
    coge.progress.begin({
    	title: 'Preparing Download ...',
    	content: '<br>The download is being generated<br><br>',
    	width: '30em',
    	height: '10em'
    });

    $.ajax({
        data: {
            fname: "get_file_urls",
            eid: EXPERIMENT_ID
        },
        dataType: "json",
        success : function(json) {
            if (json.url) {  // finished successfully
                coge.progress.succeeded("File: " + json.filename);
                setTimeout(function() { coge.utils.open(json.url) }, 10);
            }
            else { // error occurred
                coge.progress.failed();
            }
        }
    });
}

function toggle_load_log() { // TODO: shared with GenomeInfo.pl, move to module
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
		            eid: EXPERIMENT_ID
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

function export_data() {
    coge.progress.begin({
    	title: 'Exporting ...',
    	content: 
    		"<br>Copying this experiment's data files to <br><br>" +
    		'<a class="bold" target="_blank" href="' + DISCOVERY_ENVIRONMENT.concat(IRODS_HOME) + '">'+IRODS_HOME+'</a>',
    	width: '30em',
    	height: '10em'
    });
    
    $.ajax({
        data: {
            fname: 'export_experiment_irods',
            eid: EXPERIMENT_ID
        },
        success : function(filename) {
            if (filename) {  // finished successfully
            	coge.progress.succeeded('<br><br>File: ' + filename);
            }
            else { // error occurred
                coge.progress.failed();
            }
        }
    });
}

function render_template(template, container) {
    container.empty()
        .hide()
        .append(template)
        .show();//.slideDown();
}

var snpMenu = {
    init: function() {
    	var self = this;
    	this.dialog = $('<div class="dialog_box small" style="padding:1em;"></div>');
    	var template = $($("#snp-template2").html());
    	
        var options_templates = {
            coge:     $($("#coge-snp-template").html()),
            samtools: $($("#samtools-snp-template").html()),
            platypus: $($("#platypus-snp-template").html()),
            gatk:     $($("#gatk-snp-template").html())
        };
        	
        var options_container = template.find('#snp-container');
        options_container.removeClass('hidden');
        
        this.method = template.find('#snp-method');
        this.method.removeAttr("disabled");
        
    	var render = function() {
            var selected = self.method.val();
            render_template(options_templates[selected], options_container);
            render_template(template, self.dialog);
        };
        
        $(document).on('change', '#snp-method', render);
        
        render();
        
    	this.dialog.dialog({ 
    		title: "Select SNP Analysis Options",
    		width: '32em',
    		autoOpen: false,
    		buttons: [
    	        {
    	            text: " Identify SNPs ",
    	            "class": "coge-button", //FIXME isn't working
    	            click: function() {
    	            	self.close();
    	            	setTimeout($.proxy(self.submit, self), 0);
    	            }
    	        }
    	    ]
    	});
    },
    
    open: function() {
    	this.dialog.dialog("open");
    },
    
    close: function() {
    	this.dialog.dialog("close");
    },

	is_valid: function () {
		var method = this.method.val();
		
		//TODO this can be automated
	    if (method === "coge") {
	        return { 
	            method: method,
	            'min-read-depth':   this.dialog.find("#min-read-depth").val(),
	            'min-base-quality': this.dialog.find("#min-base-quality").val(),
	            'min-allele-count': this.dialog.find("#min-allele-count").val(),
	            'min-allele-freq':  this.dialog.find("#min-allele-freq").val(),
	            scale: this.dialog.find("#scale").val()
	        };
	    } 
	    else if (method === "samtools") {
	    	return {
	            method: method,
	            'min-read-depth': this.dialog.find("#min-read-depth").val(),
	            'max-read-depth': this.dialog.find("#max-read-depth").val(),
	        };
	    } 
	    else if (method === "platypus") {
	    	return {
	            method: method
	        };
	    } 
	    else if (method === "gatk") {
	    	return {
	            method: method
	        };
	    }
	    return;
	},
	
	submit: function() { 
		var params = this.is_valid();
	
		coge.progress.begin({ 
			title: "Identifying SNPs ...",
			width: '60%',
	    	height: '50%'
		});
		newLoad = true;
	  
		// Build request
		var request = {
			type: 'analyze_snps',
			requester: {
				page: PAGE_NAME,
				url: PAGE_NAME + "?eid=" + EXPERIMENT_ID
			},
			parameters: {
				eid: EXPERIMENT_ID,
				snp_params: params
			}
		};
		
		// Submit request
		coge.services.submit_job(request)
			.done(function(response) {
		  		if (!response) {
		  			coge.progress.failed("Error: empty response from server");
		  			return;
		  		}
		  		else if (!response.success || !response.id) {
		  			coge.progress.failed("Error: failed to start workflow");
		  			return;
		  		}
		  		
		        // Start status update
		  		window.history.pushState({}, "Title", PAGE_NAME + "?eid=" + EXPERIMENT_ID + "&wid=" + response.id); // Add workflow id to browser URL
		  		coge.progress.update(response.id, response.site_url);
		    })
		    .fail(function(jqXHR, textStatus, errorThrown) {
		    	coge.progress.failed("Couldn't talk to the server: " + textStatus + ': ' + errorThrown);
		    });
	}
}

var expressionMenu = {
    init: function() {
    	var self = this;
    	this.dialog = $('<div class="dialog_box small" style="padding:1em;"></div>');
    	var template = $($("#expression-template").html()).find('#expression-container');
    	template.removeClass('hidden');
    	
    	var render = function() {
            render_template(template, self.dialog);
        };
        
        render();
        
    	this.dialog.dialog({ 
    		title: "Select Expression Analysis Options",
    		width: '30em',
    		autoOpen: false,
    		buttons: [
    	        {
    	            text: " Measure Expression ",
    	            "class": "coge-button", //FIXME isn't working
    	            click: function() {
    	            	self.close();
    	            	setTimeout($.proxy(self.submit, self), 0);
    	            }
    	        }
    	    ]
    	});
    },
    
    open: function() {
    	this.dialog.dialog("open");
    },
    
    close: function() {
    	this.dialog.dialog("close");
    },

	is_valid: function () {
	    return {
	    	'-Q': this.dialog.find("[id='-Q']").val()
	    };
	},
	
	submit: function() { 
		var params = this.is_valid();
	
		coge.progress.begin({ 
			title: "Measuring Expression ...",
			width: '60%',
	    	height: '50%'
		});
		newLoad = true;
	  
		// Build request
		var request = {
			type: 'analyze_expression',
			requester: {
				page: PAGE_NAME,
				url: PAGE_NAME + "?eid=" + EXPERIMENT_ID
			},
			parameters: {
				eid: EXPERIMENT_ID,
				expression_params: params
			}
		};
		
		// Submit request
		coge.services.submit_job(request)
			.done(function(response) {
		  		if (!response) {
		  			coge.progress.failed("Error: empty response from server");
		  			return;
		  		}
		  		else if (!response.success || !response.id) {
		  			coge.progress.failed("Error: failed to start workflow");
		  			return;
		  		}
		  		
		        // Start status update
		  		window.history.pushState({}, "Title", PAGE_NAME + "?eid=" + EXPERIMENT_ID + "&wid=" + response.id); // Add workflow id to browser URL
		  		coge.progress.update(response.id, response.site_url);
		    })
		    .fail(function(jqXHR, textStatus, errorThrown) {
		    	coge.progress.failed("Couldn't talk to the server: " + textStatus + ': ' + errorThrown);
		    });
	}
}

var metaplotMenu = {
    init: function() {
    	var self = this;
    	this.dialog = $('<div class="dialog_box small" style="padding:1em;"></div>');
    	var template = $($("#metaplot-template").html()).find('#metaplot-options');
    	template.removeClass('hidden');
    	
    	var render = function() {
            render_template(template, self.dialog);
        };
        
        render();
        
    	this.dialog.dialog({ 
    		title: "Select Metaplot Analysis Options",
    		width: '20em',
    		autoOpen: false,
    		buttons: [
    	        {
    	            text: " Create Metaplot ",
    	            "class": "coge-button", //FIXME isn't working
    	            click: function() {
    	            	self.close();
    	            	setTimeout($.proxy(self.submit, self), 0);
    	            }
    	        }
    	    ]
    	});
    },
    
    open: function() {
    	this.dialog.dialog("open");
    },
    
    close: function() {
    	this.dialog.dialog("close");
    },

	is_valid: function () {
        return { 
            'outer':  this.dialog.find("#metaplot-outer").val(),
            'inner':  this.dialog.find("#metaplot-inner").val(),
            'window': this.dialog.find("#metaplot-window").val()
        };
	},
	
	submit: function() { 
		var params = this.is_valid();
	
		coge.progress.begin({ 
			title: "Creating Metaplot ...",
			width: '60%',
	    	height: '50%'
		});
		newLoad = true;
	  
		// Build request
		var request = {
			type: 'analyze_metaplot',
			requester: {
				page: PAGE_NAME,
				url: PAGE_NAME + "?eid=" + EXPERIMENT_ID
			},
			parameters: {
				eid: EXPERIMENT_ID,
				methylation_params: {
					metaplot_params: params
				}
			}
		};
		
		// Submit request
		coge.services.submit_job(request)
			.done(function(response) {
		  		if (!response) {
		  			coge.progress.failed("Error: empty response from server");
		  			return;
		  		}
		  		else if (!response.success || !response.id) {
		  			coge.progress.failed("Error: failed to start workflow");
		  			return;
		  		}
		  		
		        // Start status update
		  		window.history.pushState({}, "Title", PAGE_NAME + "?eid=" + EXPERIMENT_ID + "&wid=" + response.id); // Add workflow id to browser URL
		  		coge.progress.update(response.id, response.site_url);
		    })
		    .fail(function(jqXHR, textStatus, errorThrown) {
		    	coge.progress.failed("Couldn't talk to the server: " + textStatus + ': ' + errorThrown);
		    });
	}
}

var diversityMenu = {
    init: function() {
    	var self = this;
    	this.dialog = $('<div class="dialog_box small" style="padding:1em;"></div>');
    	var template = $($("#diversity-template").html()).find('#diversity-container');
    	template.removeClass('hidden');
    	
    	var render = function() {
            render_template(template, self.dialog);
        };
        
        render();
        
    	this.dialog.dialog({ 
    		title: "Select Diversity Analysis Options",
    		width: '30em',
    		autoOpen: false,
    		buttons: [
    	        {
    	            text: " Measure Diversity ",
    	            "class": "coge-button", //FIXME isn't working
    	            click: function() {
    	            	self.close();
    	            	setTimeout($.proxy(self.submit, self), 0);
    	            }
    	        }
    	    ]
    	});
    },
    
    open: function() {
    	this.dialog.dialog("open");
    },
    
    close: function() {
    	this.dialog.dialog("close");
    },

	is_valid: function () {
	    return {
	    	//'-Q': this.dialog.find("[id='-Q']").val()
	    };
	},
	
	submit: function() { 
		var params = this.is_valid();
	
		coge.progress.begin({ 
			title: "Measuring Diversity ...",
			width: '60%',
	    	height: '50%'
		});
		newLoad = true;
	  
		// Build request
		var request = {
			type: 'analyze_diversity',
			requester: {
				page: PAGE_NAME,
				url: PAGE_NAME + "?eid=" + EXPERIMENT_ID
			},
			parameters: {
				eid: EXPERIMENT_ID,
				diversity_params: params
			}
		};
		
		// Submit request
		coge.services.submit_job(request)
			.done(function(response) {
		  		if (!response) {
		  			coge.progress.failed("Error: empty response from server");
		  			return;
		  		}
		  		else if (!response.success || !response.id) {
		  			coge.progress.failed("Error: failed to start workflow");
		  			return;
		  		}
		  		
		        // Start status update
		  		window.history.pushState({}, "Title", PAGE_NAME + "?eid=" + EXPERIMENT_ID + "&wid=" + response.id); // Add workflow id to browser URL
		  		coge.progress.update(response.id, response.site_url);
		    })
		    .fail(function(jqXHR, textStatus, errorThrown) {
		    	coge.progress.failed("Couldn't talk to the server: " + textStatus + ': ' + errorThrown);
		    });
	}
}

function analyze_diversity() {
	if (POPGEN_RESULT_URL) 
		coge.utils.open(POPGEN_RESULT_URL);
	else 
		diversityMenu.open();
}

function remove_experiment_tag (opts) {
    etid = opts.etid;
    $.ajax({
        data: {
            fname: 'remove_experiment_tag',
            eid: EXPERIMENT_ID,
            etid: etid,
        },
        success : function(val) {
            get_experiment_info();
        },
    });
}

function get_experiment_tags () {
    $.ajax({
        data: {
            fname: 'get_experiment_tags'
        },
        success : function(val) {
            var items = jQuery.parseJSON(val);
            $("#edit_tag_name").autocomplete("option", "source", items).autocomplete("search");
        },
    });
}

function toggle_favorite(img) {
	$.ajax({
		data: {
			fname: 'toggle_favorite',
			eid: EXPERIMENT_ID,
		},
		success :  function(val) {
			$(img).attr({ src: (val == '0' ? "picts/star-hollow.png" : "picts/star-full.png") });
		}
	});
}