/*global document,$,alert,window,JSON */

// Global experiment data
var current_experiment = {};

// Supported input file types
var POLY_FILES  = [ "vcf", "gcvf" ];
var ALIGN_FILES = [ "bam" ];
var SEQ_FILES   = [ "fastq", "fq", "sra" ];
var QUANT_FILES = [ "csv", "tsv", "bed", "wig", "bw", "gff", "gtf" ];
var SUPPORTED_FILE_TYPES = Array.prototype.concat.call(QUANT_FILES, ALIGN_FILES, SEQ_FILES, POLY_FILES);

function create_source() {
    var name = $('#edit_source_name').val();
    var desc = $('#edit_source_desc').val();
    var link = $('#edit_source_link').val();
    $.ajax({
        data: {
            fname: 'create_source',
            name: name,
            desc: desc,
            link: link,
        },
        success : function(name) {
            $('#create_new_source_dialog').dialog('close');
            if (name) {
                $('#edit_source').val(name);
            }
        }
    });
}

function search_genomes (search_term) {
	if (!search_term || search_term.length < 3)
		return;
	
	var edit_genome = $("#edit_genome");
	edit_genome.autocomplete("close");
	var spinner = $('#edit_genome_busy');
	spinner.show();
	
	coge.services.search_genomes(search_term, { fast: true, sort: true })
		.done(function(response) { // success
			if (response && response.genomes) {
				var results = response.genomes.map(function(obj) {
					var label = obj.info.replace(/&#x1f512;/g, "\uD83D\uDD12"); // Lock symbol // TODO make client-side like certified & favorite
					if (obj.certified)
					    label = '\u2705 ' + label;
                    if (obj.favorited)
					    label = '\u2B50 ' + label;
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

function search_users (search_term) {
	coge.services.search_users(search_term)
		.done(function(result) {
			if (result && result.users) {
				var transformed = result.users.map(function(obj) {
					return obj.user_name;
				});
				$("#edit_user")
					.autocomplete({source: transformed})
					.autocomplete("search");
			}
		})
		.fail(function() {
			//TODO
		});
}

function search_notebooks (search_term) {
	coge.services.search_notebooks(search_term)
		.done(function(data) { 
			if (data.notebooks) {
				var items = [];
				data.notebooks.forEach(function(notebook) {
					items.push({label:notebook.name,value:notebook.id});
			    });
			    $("#edit_notebook")
			    	.autocomplete({source: items})
			    	.autocomplete("search");
			}
		})
		.fail(function() {
			//TODO
		});
}

function activate_on_input(elements, target) {
	if (elements) {
		elements.forEach(function(e) {
			if ($('#'+e).val())
				$('#'+target).removeClass('ui-state-disabled');
		});
	}
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
    //this.gid = opts.gid;
    this.onError = opts.onError;
    this.sources = undefined;
    this.title = "Describe Experiment";
    this.initialize();
}

$.extend(ExperimentDescriptionView.prototype, {
    initialize: function() {
        this.el = $($("#description-template").html());
        this.edit_source = this.el.find("#edit_source");
        this.edit_genome = this.el.find("#edit_genome");
        this.get_sources();
    },

    get_sources: function() {
        var self = this;

        return $.ajax({
            dataType: "json",
            data: {
                fname: 'get_sources',
            },
            success : function(sources) {
                self.sources = sources;

                if (sources)
                    self.edit_source.autocomplete({source: sources});
            },
        });
    },

    render: function() {
        var self = this;

        this.isSRA = (this.experiment.data && this.experiment.data[0].file_type == 'sra');
        this.el.find('#edit_name,#edit_description,#edit_version,#edit_source').prop('disabled', (this.isSRA ? true : false));
        if (this.isSRA)
            this.el.find('#sra_message').show();
        else
            this.el.find('#sra_message').hide();

        // Set experiment metadata if from SRA
		if (this.isSRA) {
           this.metadata = {
                name: 'blah',
                description: '',
                version: '1',
                restricted: 0,
                source_name: 'NCBI-SRA'
            };
        }

        if (this.metadata) {
            this.el.find('#edit_name').val(this.metadata.name);
            this.el.find('#edit_description').val(this.metadata.description);
            this.el.find('#edit_version').val(this.metadata.version);
            this.edit_source.val(this.metadata.source_name);

            if (!this.metadata.restricted)
                this.el.find('#restricted').removeAttr('checked');

            /*this.el.find('#edit_genome').val(this.metadata.genome);*/
        }

        var edit_genome = this.edit_genome;

        edit_genome.unbind().change(function() {
            // Reset gid when item has changed
            self.gid = undefined;
        });

        edit_genome.autocomplete({
            source:[],
            select: function(event, ui) {
                $(this).val(ui.item.label);
                self.gid = ui.item.value;
                return false; // Prevent the widget from inserting the value.
            },
            focus: function(event, ui) {
                //$("#edit_genome").val(ui.item.label);
                return false; // Prevent the widget from inserting the value.
            }
        });
        
        edit_genome.keyup(function() {
        	coge.utils.wait_to_search(search_genomes, self.edit_genome.get(0));
        });
        edit_genome.click(function() {
        	$(this).autocomplete('search');
        });

        this.edit_source.autocomplete({source: this.sources});
    },

    is_valid: function() {
        var name = this.el.find('#edit_name').val();
        var description = this.el.find('#edit_description').val();
        var version = this.el.find('#edit_version').val();
        var restricted = this.el.find('#restricted').is(':checked');
        var genome = this.el.find('#edit_genome').val();

        if (!name) {
            if (this.onError)
            	this.onError('Please specify an experiment name.');
            return false;
        }
        
        if (!version) {
//        	if (this.onError)
//            	this.onError('Please specify an experiment version.');
//            return false;
            version = "1";
        }

        var source = $('#edit_source').val();
        if (!source) {
        	if (this.onError)
            	this.onError('Please specify a data source.');
            return false;
        }

        if (!genome || !this.gid) {
        	if (this.onError)
            	this.onError('Please specify a genome.');
            return false;
        }

       $.extend(this.experiment, {
            metadata: {
                name: coge.utils.removeSpecialChars(name),
                description: coge.utils.removeSpecialChars(description),
                version: coge.utils.removeSpecialChars(version),
                restricted: restricted,
                source_name: coge.utils.removeSpecialChars(source),
                genome: genome,
            },
            gid: this.gid
        });

        return true;
    },
});

function MethylationView(opts) {
	if (opts)
		this.format = opts.format;
    this.data = {};
    this.initialize();
};

$.extend(MethylationView.prototype, {
    initialize: function () {
        this.el = $($("#methyl-template").html());
        this.checkbox = this.el.find("#methyl");
        this.container = this.el.find("#methyl-container");
        this.templates = {
            bismark: $($("#bismark-methyl-template").html()),
            bwameth: $($("#bwameth-methyl-template").html())
        };
        this.metaplot_template = $($("#metaplot-template").html());
        this.metaplot_checkbox = this.metaplot_template.find("#metaplot-enable");
    },

    render: function () {
        this.checkbox.unbind().change(this.update.bind(this));
        this.metaplot_checkbox.unbind().change(this.update_metaplot.bind(this));
    },

    update: function() {
        var checkbox = this.checkbox;//$(ev.target);
        var enabled = checkbox.is(":checked"); 
        var selected = $("#alignment").val(); // FIXME pass alignment in as argument to constructor
        var template = this.templates[selected];
        var el = $(document.getElementById(selected));
        
        if (enabled && template) {
            this.data.methylation_params = $.extend({}, this.data.methylation_params, {method: selected});
        	el.show();
            this.container.slideDown();
        	render_template(template, this.container);
        	this.container.append(this.metaplot_template);
        }
        else {
            this.data.methylation_params = undefined;
            checkbox.attr('checked', false); // uncheck it
            this.container.hide();
            if (enabled && !template)
            	this.container.html('<span class="alert indent">Please select one of these two aligners above:  Bismark or BWAmeth</span>').show();
        }
    },
    
    update_metaplot: function() {
    	var metaplot_enabled = this.metaplot_checkbox.is(":checked");
    	
    	if (metaplot_enabled) {
    		this.metaplot_template.find('#metaplot-options').removeClass('hidden');
    	}
    	else {
    		this.metaplot_template.find('#metaplot-options').addClass('hidden');
    	}
    },

    is_valid: function() {
        var enabled = this.checkbox.is(":checked");
        var method = $("#alignment").val(); // FIXME pass alignment in as argument to constructor
        var paired = ( this.format && this.format.is_paired() );
        
        if (enabled) {
            if (method === "bismark") {
                this.data.methylation_params = {
                    method: method,
                    'bismark-deduplicate': this.el.find('#bismark-deduplicate').is(":checked"),
                    'bismark-min_converage': this.el.find('#bismark-min_coverage').val(),
                };
                if (paired) {
                    this.data.methylation_params['--ignore'] = this.el.find('#--ignore').val();
                    this.data.methylation_params['--ignore_3prime'] = this.el.find('#--ignore_3prime').val();
                    this.data.methylation_params['--ignore_r2'] = this.el.find('#--ignore_r2').val();
                    this.data.methylation_params['--ignore_3prime_r2'] = this.el.find('#--ignore_3prime_r2').val();
                }
                else { // single-ended
                    this.data.methylation_params['--ignore'] = this.el.find('#--ignore').val();
                    this.data.methylation_params['--ignore_3prime'] = this.el.find('#--ignore_3prime').val();
                }
            }
            else if (method === "bwameth") {
                this.data.methylation_params = {
                    method: method,
                    'picard-deduplicate': this.el.find('#picard-deduplicate').is(":checked"),
                    'pileometh-min_converage': this.el.find('#pileometh-min_coverage').val(),
                    '--OT': this.el.find('#--OT').val(),
                    '--OB': this.el.find('#--OB').val(),
                };
            }
            
            // Add metaplot params
            var metaplot_enabled = this.metaplot_checkbox.is(":checked");
            if (metaplot_enabled) {
	            this.data.methylation_params = $.extend(this.data.methylation_params, 
	            	{
	            		metaplot_params: {
	            			outer:  this.el.find('#metaplot-outer').val(),
	            			inner:  this.el.find('#metaplot-inner').val(),
	            	        window: this.el.find('#metaplot-window').val()
	            		}
	            	}
	            );
            }
        }
        return true;
    },

    get_options: function() {
    	return this.data;
    },
});

function FindSNPView() {
    this.data = {};
    this.initialize();
};

$.extend(FindSNPView.prototype, {
    initialize: function() {
        this.el = $($("#snp-template").html());
        this.snp_container = this.el.find("#snp-container");

        this.snp_templates = {
            coge:     $($("#coge-snp-template").html()),
            samtools: $($("#samtools-snp-template").html()),
            platypus: $($("#platypus-snp-template").html()),
            gatk:     $($("#gatk-snp-template").html())
        };
    },

    render: function() {
        var self = this;

        var method = this.el.find("#snp-method");
        this.el.find("#snps").unbind().change(this.update_snp.bind(this));

        // Events to rebind when the view is added to the dom
        method.unbind().change(function() {
            var selected = $(this).val();
            render_template(self.snp_templates[selected], self.snp_container);
        });

        if (this.data.snp_params) {
            method.val(this.data.snp_params.method);
        }
    },

    // Callback to display the selected snp pipeline
    update_snp: function (ev) {
        var enabled = $(ev.target).is(":checked"),
            method = this.el.find("#snp-method");

        var el = $(document.getElementById(method.val()));

        if (enabled) {
            this.data.snp_params = $.extend({}, this.data.snp_params, { method: method.val() });
            el.show();
            method.removeAttr("disabled");
            method.parent('div').slideDown();
            this.snp_container.slideDown();
            var selected = $("#snp-method").val();
            render_template(this.snp_templates[selected], this.snp_container);
        } 
        else {
            this.data.snp_params = undefined;
            method.attr("disabled", 1);
            method.parent('div').slideUp();
            this.snp_container.slideUp();
        }
    },

    is_valid: function() {
        // SNP pipeline
        var enabled = this.el.find("#snps").is(":checked");
        var method = this.el.find("#snp-method").val();

        if (enabled) {
        	//TODO this can be automated
            if (method === "coge") {
                this.data.snp_params = {
                    method: method,
                    'min-read-depth':   this.el.find("#min-read-depth").val(),
                    'min-base-quality': this.el.find("#min-base-quality").val(),
                    'min-allele-count': this.el.find("#min-allele-count").val(),
                    'min-allele-freq':  this.el.find("#min-allele-freq").val(),
                    scale: this.el.find("#scale").val()
                };
            } else if (method === "samtools") {
                this.data.snp_params = {
                    method: method,
                    'min-read-depth': this.el.find("#min-read-depth").val(),
                    'max-read-depth': this.el.find("#max-read-depth").val(),
                };
            } else if (method === "platypus") {
                this.data.snp_params = {
                    method: method
                };
            } else if (method === "gatk") {
                this.data.snp_params = {
                    method: method,
                    '-stand_call_conf': this.el.find("[id='-stand_call_conf']").val(),
                    '-stand_emit_conf': this.el.find("[id='-stand_emit_conf']").val()
                };
            }
        }
        return true;
    },

    get_options: function() {
        return this.data;
    },
});

function ReadFormatView() {
    this.initialize();
    this.data = {};
}

$.extend(ReadFormatView.prototype, {
    initialize: function() {
        this.el = $($("#format-template").html());
        this.container = this.el.find("#format-container");
    },

    render: function() {
    },

    is_valid: function() {
        return true;
    },

    get_options: function() {
        this.data.read_params = {
        	read_type: this.el.find("#read_type :checked").val(),
        	encoding: this.el.find("#encoding :checked").val(),
        };

        return this.data;
    },
    
    is_paired: function() {
    	return (this.get_options().read_type === 'paired');
    }
});

function TrimmingView() {
    this.data = {};
    this.initialize();
}

$.extend(TrimmingView.prototype, {
    initialize:function() {
        this.el = $($("#trim-template").html());
        this.container = this.el.find("#trim-container");
        this.templates = {
        	none:       $($("#none-template").html()),
            cutadapt:   $($("#cutadapt-template").html()),
            trimgalore: $($("#trimgalore-template").html()),
        };
    },

    render: function() {
        // jQuery events
        this.el.find("[name=trimmer]").unbind().click(this.update.bind(this));
        this.update();
    },

    // Callback to display the selected aligner
    update: function() {
        var selected = this.el.find("#trimming :checked").val();
        render_template(this.templates[selected], this.container);
    },

    is_valid: function() {
        var trimmer = this.el.find("#trimming :checked").val();
        // Pick the aligner and set the options
        if (trimmer === "cutadapt") {
            this.data = {
                trimming_params: {
                	'trimmer': 'cutadapt',
                    '-q': this.el.find("[id='-q']").val(),
                    '-m': this.el.find("[id='-m']").val()
                    //'--quality-base': this.el.find("[id='--quality-base']").val()
                }
            };
        } 
        else if (trimmer === "trimgalore") {
            this.data = {
                trimming_params: {
                	'trimmer': 'trimgalore',
                	'-q': this.el.find("[id='-q']").val(),
                }
            };
        } 
        else {
        	this.data = {}; // Skip trimming
        }

        return true;
    },

    get_options: function() {
        return this.data;
    }
});

function AlignmentView(opts) {
	if (opts) {
		this.onChange = opts.onChange;
	}
    this.data = {};
    this.initialize();
}

$.extend(AlignmentView.prototype, {
    initialize:function() {
        this.el = $($("#align-template").html());
        this.container = this.el.find("#align-container");
        this.templates = {
            gsnap:   $($("#gsnap-template").html()),
            bwa:     $($("#bwa-template").html()),
            bowtie2: $($("#bowtie2-template").html()),
            tophat:  $($("#tophat-template").html()),
            hisat2:  $($("#hisat2-template").html()),
            bismark: $($("#bismark-template").html()),
            bwameth: $($("#bwameth-template").html())
        };
    },

    render: function() {
    	var self = this;
    	
        var aligner = this.el.find("#alignment");
        aligner.unbind().change(this.update.bind(this));

        if (this.data.alignment_params) {
        	aligner.val(this.data.alignment_params.tool);
        }
        this.update();
    },

    update: function() {
        var selected = this.el.find("#alignment").val();
        var el = $(document.getElementById(selected));
        this.data.alignment_params = $.extend({}, this.data.alignment_params, { tool: selected });
        el.show();
        this.container.show();
        render_template(this.templates[selected], this.container);
        
        if (this.onChange)
        	this.onChange();
    },

    is_valid: function() {
        var aligner = this.el.find("#alignment :checked").val();
        // Pick the aligner and set the options
        if (aligner === "gsnap") {
            this.data = {
                alignment_params: { //TODO is there a way to automate this parameter passing?
                    tool: "gsnap",
                    '-N': this.el.find("[id='-N']").val(),
                    '-n': this.el.find("[id='-n']").val(),
                    '-Q': this.el.find("[id='-Q']").is(":checked"),
                    '--gap-mode': this.el.find("[id='--gap-mode']").val(),
                    '--nofails': this.el.find("[id='--nofails']").is(":checked"),
                }
            };
            
            if (this.el.find("[id='--max-mismatches-chk']").is(":checked")) {
            	this.data.alignment_params['--max-mismatches'] = this.el.find("[id='--max-mismatches']").val();
            }
        }
        else if (aligner === "bwa") {
            this.data = {
        		alignment_params: {
        			tool: "bwa",
        			'-M': this.el.find("[id='-M']").is(":checked"),
        			'-R': this.el.find("[id='-R']").val(),
        		}
        	}
        }
        else if (aligner === "bowtie2") {
        	this.data = {
        		alignment_params: {
        			tool: "bowtie2",
        			'presets': this.el.find("[id='presets']").val(),
        			'--rg-id':    this.el.find("[id='--rg-id']").val(),
        		}
        	}
        }
        else if (aligner === "tophat") {
            this.data = {
                alignment_params: {
                    tool: "tophat",
                    '-g': this.el.find("[id='-g']").val(),
                }
            }
        } 
        else if (aligner === "hisat2") {
        	this.data = {
        		alignment_params: {
        			tool: "hisat2"
        		}
        	}
        }
        else if (aligner === "bismark") {
        	this.data = {
        		alignment_params: {
        			tool: "bismark",
        			'-N': this.el.find("[id='-N']").val(),
        			'-L': this.el.find("[id='-L']").val(),
        		}
        	}
        }
        else if (aligner === "bwameth") {
        	this.data = {
        		alignment_params: {
        			tool: "bwameth"
        		}
        	}
        }
        else { // should never happen
        	console.error('Invalid aligner');
        }

        return true;
    },

    get_options: function() {
        return this.data;
    }
});

function PolymorphismView() {
    this.initialize();
}

$.extend(PolymorphismView.prototype, {
    initialize: function() {
        this.el = $($("#poly-template").html());
    },

    is_valid: function() {
        return true;
    },

    get_options: function() {
        return {};
    },
});

function AlignmentOptionView() {
    this.initialize();
}

$.extend(AlignmentOptionView.prototype, {
    initialize: function() {
        this.snp_view = new FindSNPView();
        this.expression_view = new ExpressionView();

        this.layout_view = new LayoutView({
            template: "#align-option-template",

            layout: {
                "#expression-view": this.expression_view,
                "#snp-view": this.snp_view
            }
        });

        this.el = this.layout_view.el;
    },

    render: function() {
        this.layout_view.renderLayout();
    },

    is_valid: function() {
        return this.snp_view.is_valid()
        	   && this.expression_view.is_valid();
    },

    get_options: function() {
        return $.extend(this.expression_view.get_options(),
        				this.snp_view.get_options());
    },
});

function QuantativeView(){
    this.initialize();
    this.data = {};
}

$.extend(QuantativeView.prototype, {
    initialize: function() {
        this.el = $($("#quant-template").html());
        this.container = this.el.find("#normalize_method");
    },

    render: function() {
        this.el.find("#normalize").unbind().change(this.toggleAnalysis.bind(this));
    },

    toggleAnalysis: function() {
        this.enabled = this.el.find("#normalize").is(":checked");

        if (this.enabled) {
            this.container.slideDown();
        } else {
            this.container.slideUp();
        }
    },

    is_valid: function() {
        this.data.normalize = this.el.find("#normalize").is(":checked");
        return true;
    },

    get_options: function() {
        if (this.enabled)
            this.data.normalize_method = this.el.find("#percentage").is(":checked") ? 'percentage' : this.el.find("#log10").is(":checked") ? 'log10' : this.el.find("#loge").is(":checked") ? 'loge' : null;
        return this.data;
    },
});

function ExpressionView() {
    this.initialize();
    this.data = {};
}

$.extend(ExpressionView.prototype, {
    initialize: function() {
        this.el = $($("#expression-template").html());
        this.enabled = false;
        this.container = this.el.find("#expression-container");
    },

    render: function() {
        this.el.find("#expression").unbind().change(this.toggleAnalysis.bind(this));
    },

    toggleAnalysis: function() {
        this.enabled = this.el.find("#expression").is(":checked");

        if (this.enabled) {
            this.container.slideDown();
        } else {
            this.container.slideUp();
        }
    },

    is_valid: function() {
        return true;
    },

    get_options: function() {
        if (this.enabled) {
            this.data.expression_params = {
                '-q': this.el.find("[id='-q']").val(),
                '-frag-bias-correct' : this.el.find("[id='-frag-bias-correct']").is(":checked"),
                '-multi-read-correct': this.el.find("[id='-multi-read-correct']").is(":checked")
            };
        }

        return this.data;
    },
});

function ChIPSeqView(opts) {
	this.experiment = opts.experiment;
    this.initialize();
    this.data = {};
}

$.extend(ChIPSeqView.prototype, {
    initialize: function() {
        this.el = $($("#chipseq-template").html());
        this.enabled = false;
        this.container = this.el.find("#chipseq-container");
        this.template = $(this.container.html());
    },

    render: function() {
        this.el.find("#chipseq").unbind().change(this.update.bind(this));
    },

    update: function() {
    	var checkbox = this.el.find("#chipseq");
    	
        var selected = $("#alignment").val(); // FIXME pass alignment in as argument to constructor
        if (selected != 'gsnap' && selected != 'bwa' && selected != 'bowtie2' && selected != 'tophat' && selected != 'hisat2') {
        	this.container.html('<span class="alert indent">Please select one of these aligners above: GSNAP, Bowtie2, TopHat2, or HISAT2</span>').show();
        	checkbox.attr('checked', false); // uncheck it
        	return;
        }
        
        var data = this.experiment.data;
        if (!data || data.length < 3) {
        	this.container.html('<span class="alert indent">This analysis requires 3 input files (the input and two replicates).</span>').show();
        	checkbox.attr('checked', false); // uncheck it
        	return;
        }

        this.enabled = checkbox.is(":checked");
        if (this.enabled) {
        	render_template(this.template, this.container);
        	
        	// Add input files to dropdown
            var select = this.el.find("#chipseq-input");
            var data = this.experiment.data;
            for(i = 0; i < data.length; i++) {
            	var file = data[i];
            	select.append("<option>" + file.name + "</option>");
            }
            
            this.container.slideDown();
        }
        else 
            this.container.slideUp();
    },

    is_valid: function() {
        return true;
    },

    get_options: function() {
        if (this.enabled) {
        	var notSelected = $('#chipseq-input option:not(:selected)').map(function(index, element) {
        		return $(element).val()
        	}).get().join(', ');
        	
            this.data.chipseq_params = {
            	'input':  this.el.find("#chipseq-input").val(),
            	'replicates': notSelected,
                '-size':  this.el.find("[id='-size']").val(),
                '-gsize': this.el.find("[id='-gsize']").val(),
                '-norm':  this.el.find("[id='-norm']").val(),
                '-fdr':   this.el.find("[id='-fdr']").val(),
                '-F':     this.el.find("[id='-F']").val(),
            };
        }

        return this.data;
    },
});

function FastqView(opts) {
	this.experiment = opts.experiment
    this.initialize();
}

$.extend(FastqView.prototype, {
    initialize: function() {
    	this.read_view = new ReadFormatView();
    	this.trim_view = new TrimmingView();
        this.expression_view = new ExpressionView();
        this.snp_view = new FindSNPView();
        this.methylation_view = new MethylationView({ format: this.read_view });
        this.chipseq_view = new ChIPSeqView({ experiment: this.experiment });
        this.align_view = new AlignmentView({ onChange: this.methylation_view.update.bind(this.methylation_view) });

        this.layout_view = new LayoutView({
            template: "#fastq-template",

            layout: {
                '#read-view': this.read_view,
                '#trim-view': this.trim_view,
                "#align-view": this.align_view,
                "#expression-view": this.expression_view,
                "#snp-view": this.snp_view,
                "#methylation-view": this.methylation_view,
                "#chipseq-view": this.chipseq_view
            }
        });

        // pass through to the layout
        this.el = this.layout_view.el;
    },

    render: function() {
        this.layout_view.renderLayout();
    },

    is_valid: function() {
    	return (   this.read_view.is_valid()
    			&& this.trim_view.is_valid()
    			&& this.align_view.is_valid()
    			&& this.snp_view.is_valid()
    			&& this.expression_view.is_valid()
    			&& this.methylation_view.is_valid()
    			&& this.chipseq_view.is_valid());
    },

    get_options: function() {
        return $.extend(this.read_view.get_options(),
		        		this.trim_view.get_options(),
		        		this.align_view.get_options(),
                        this.expression_view.get_options(),
                        this.snp_view.get_options(),
                        this.methylation_view.get_options(),
                        this.chipseq_view.get_options());
    },
});

function GeneralOptionsView(opts) {
    this.data = {};
    this.initialize();
    this.onError = opts.onError;
}

$.extend(GeneralOptionsView.prototype, {
    initialize: function() {
        this.el = $($("#general-options-template").html());
        this.edit_notebook = this.el.find("#edit_notebook");
        this.notebook_container = this.el.find("#notebook-container");
    },

    is_valid: function() {
        var notebook = this.edit_notebook.val();

        this.data.notebook = this.el.find("#notebook").is(":checked");
        this.data.notebook_type = this.el.find("[name=notebook] :checked").val();
        this.data.notebook_name = notebook;
        this.data.notebook_id = this.notebook_id;
        this.data.email = this.el.find("#email").is(":checked");

        if (this.data.notebook && this.data.notebook_type === "existing" && 
        		(!notebook || !this.notebook_id))
        {
        	if (this.onError)
            	this.onError('Please specify a notebook.');
            return false;
        }

        return true;
    },

    get_options: function() {
        return this.data;
    },
    
    render: function() {
        var self = this;

        // Events
        this.el.find("#notebook").unbind().change(this.toggleNotebook.bind(this));
        this.el.find("[name=notebook]").unbind().click(function() {
        	var option = $(this).val();
        	self.edit_notebook.prop("disabled", (option === 'new' ? true : false));
        });
        this.edit_notebook.unbind().change(function() {
            // Reset notebook_id when item has changed
            self.notebook_id = undefined;
        });

        // Setup "source" autocomplete
        this.edit_notebook.autocomplete({
            source:[],
            select: function(event, ui) {
                $(this).val(ui.item.label);
                self.notebook_id = ui.item.value;
                return false; // Prevent the widget from inserting the value.
            },

            focus: function(event, ui) {
                return false; // Prevent the widget from inserting the value.
            }
        });
    },
    
    toggleNotebook: function() {
        this.notebook_enabled = this.el.find("#notebook").is(":checked");

        if (this.notebook_enabled) 
            this.notebook_container.slideDown();
        else 
            this.notebook_container.slideUp();
    }
});

function AdminOptionsView() {
    this.data = {};
    this.initialize();
}

$.extend(AdminOptionsView.prototype, {
    initialize: function() {
        this.el = $($("#admin-options-template").html());
        this.edit_user = this.el.find("#edit_user");
    },

    render: function() {
        // jQuery UI
        this.edit_user.unbind().autocomplete({
            source:[],
            focus: function() { return false; },
        });
    },

    is_valid: function() {
        //var ignore_cb = this.el.find('#ignore_missing_chrs');
        //this.ignore_missing_chrs = ignore_cb.is(':checked');
        if (this.edit_user.val()) {
            this.data.user = this.edit_user.val();
        }

        return true;
    },

    get_options: function() {
        return this.data;
    },
});

function OptionsView(opts) {
    this.experiment = opts.experiment;
    this.admin = opts.admin;
    this.onError = opts.onError;
    this.title = "Select Options";
    this.initialize();
}

$.extend(OptionsView.prototype, {
    initialize: function() {
        this.admin_view = new AdminOptionsView();
        this.general_view = new GeneralOptionsView({
            experiment: this.experiment,
            onError: this.onError
         });
        this.layout_view = new LayoutView({
            template: "#options-layout-template",
            layout: {
                "#general-options": this.general_view
            },
            experiment: this.experiment
        });

        if (this.admin)
            this.layout_view.updateLayout({ "#admin-options": this.admin_view });

        this.el = this.layout_view.el;
        
        this.experiment.options = {};
    },

    // Validate and add all options to the experiment
    is_valid: function() {
        if (!this.analysis_view.is_valid() || !this.general_view.is_valid())
            return false;

        var options = $.extend({}, this.general_view.get_options(), this.analysis_view.get_options());
        if (this.admin) {
            if (!this.admin_view.is_valid())
                return false;
            $.extend(options, this.admin_view.get_options());
        }

        $.extend(this.experiment.options, options);
        return true;
    },

    render: function() {
        var file_type = this.experiment.data[0].file_type;
        if (!file_type) {
            if (this.onError)
            	this.onError("Please set the file type.");
            return;
        }

        //FIXME: An aggregate view should add analysis options for multiple file types
        if ($.inArray(file_type, POLY_FILES) > -1)
            this.analysis_view = new PolymorphismView();
        else if ($.inArray(file_type, SEQ_FILES) > -1)
            this.analysis_view = new FastqView({experiment: this.experiment});
        else if ($.inArray(file_type, QUANT_FILES) > -1)
            this.analysis_view = new QuantativeView();
        else if ($.inArray(file_type, ALIGN_FILES) > -1)
            this.analysis_view = new AlignmentOptionView();

        this.layout_view.updateLayout(
            {"#analysis-options": this.analysis_view}
        );

        // Render the views added to the layout view
        this.layout_view.renderLayout();
    },
});

function load(experiment) {
	coge.progress.begin();
    newLoad = true;

    // Set job type based on data
//    console.log(experiment);
//    var job_type = 'load_experiment';
//    if (experiment.data[0].type == 'sra') {
//        job_type = 'load_sra';
//    }

	// Convert request into format for job service
	var request = {
		type: job_type,
		requester: {
			page: PAGE_NAME
		},
		parameters: {
			genome_id:          experiment.gid,
			metadata:           experiment.metadata,
			read_params:    	experiment.options.read_params,
			trimming_params:    experiment.options.trimming_params,
			alignment_params:   experiment.options.alignment_params,
			expression_params:  experiment.options.expression_params,
			snp_params:         experiment.options.snp_params,
            methylation_params: experiment.options.methylation_params,
            chipseq_params:     experiment.options.chipseq_params,
			normalize:          experiment.options.normalize,
			normalize_method:   experiment.options.normalize_method,
			email:              experiment.options.email,
			notebook:           experiment.options.notebook,
			notebook_name:      experiment.options.notebook_name,
			notebook_id:        experiment.options.notebook_id,
			source_data:        experiment.data,
			load_id:            load_id
		}
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
            window.history.pushState({}, "Title", "LoadExperiment.pl" + "?wid=" + response.id); // Add workflow id to browser URL
            coge.progress.update(response.id, response.site_url);
	    })
	    .fail(function(jqXHR, textStatus, errorThrown) {
	    	coge.progress.failed("Couldn't talk to the server: " + textStatus + ': ' + errorThrown);
	    });
}

function reset_load() {
    window.history.pushState({}, "Title", PAGE_NAME);

    // Reset the wizard and set description step
    initialize_wizard({
        admin: IS_ADMIN,
        metadata: current_experiment.metadata,
        gid: current_experiment.gid
    });

    $('#wizard-container').hide().fadeIn();
}

function initialize_wizard(opts) {
    current_experiment = {
        gid: opts.gid
    };

    // Create wizard and add steps
    var wizard = new Wizard({ 
    	onCompleted: load, 
    	data: current_experiment,
    	helpUrl: opts.helpUrl
    });
    wizard.addStep(new DataView( current_experiment, {
        supportedFileTypes: SUPPORTED_FILE_TYPES,
        onError: wizard.error_help.bind(wizard)
    }));
    wizard.addStep(new OptionsView({
        experiment: current_experiment,
        admin: opts.admin,
        onError: wizard.error_help
    }));
    wizard.addStep(new ExperimentDescriptionView({
        experiment: current_experiment,
        metadata: opts.metadata,
        onError: wizard.error_help.bind(wizard)
    }));
    wizard.addStep(new ConfirmationView(current_experiment));

    // Render all the wizard sub views
    wizard.render();
    
    // Add the wizard to the document
    var root = $("#wizard-container");
    root.html(wizard.el);
    
    return wizard;
}
