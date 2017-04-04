/* 
 * Wizard Components, including DataView and ConfirmationView,
 * shared by LoadGenome and LoadExperiment.
 * 
 * Requires coge.utils, coge.fileSelect
 * 
 * TODO: this module isn't name-spaced as the other coge/ modules are.
 * 
 */

// Requires a done callback and a data object to pass to the callback
function Wizard(options) {
    this.onCompleted = options.onCompleted;
    this.data = options.data;
    this.steps = [];
    this.currentIndex = 0;
    this.initialize(options);
}

$.extend(Wizard.prototype, {
    initialize: function(options) {
        this.el = $($("#wizard-template").html());
        this.tabs = this.el.find(".sections");
        this.next = this.el.find(".next");
        this.prev = this.el.find(".prev");
        this.done = this.el.find(".done");
        this.viewer = this.el.find("#step-container");
        this.notifications = this.el.find("#error_help_text");
        this.help = this.el.find(".link");

        // jQuery events
        this.prev.unbind().click(this.movePrevious.bind(this));
        this.next.unbind().click(this.moveNext.bind(this));
        this.done.unbind().click(this.submit.bind(this));
        this.help.unbind().click(function() {
        	window.open(options.helpUrl);
        });
    },
    
    error_help: function(s) {
    	this.notifications
            .html(s)
            .show()
            .delay(10*1000)
            .fadeOut(1500);
        window.scrollTo(0, 0); // scroll to top where error message is displayed
    },

    at_first: function() {
        return this.currentIndex === 0;
    },

    at_last: function() {
        return this.currentIndex >= (this.steps.length - 1);
    },

    render: function(from) {
        var titles = this.steps.map(function(step) {
            return $("<div></div>", { text:  step.title });
        });

        this.tabs.html(titles);
        titles[this.currentIndex].addClass("active");
//        for (var i=0;i<=this.currentIndex;i++)
//        	titles[i].addClass("active");
//        for (var i=0;i<this.currentIndex;i++) {
//        	var self = this;
//        	var index = i;
//        	titles[i].click(function(){self.move(index)});
//        	titles[i].css('cursor','pointer');
//        }

        var step = this.steps[this.currentIndex];
        if (step.render)
            step.render(from);
        this.viewer.html(step.el);

        if (this.at_first())
            this.prev.hide();//this.prev.attr("disabled", 1);
        else
            this.prev.removeAttr("disabled").show();

        if (this.at_last()) {
            this.next.hide();
            this.done.show();
        } else {
            this.next.show();
            this.done.hide();
        }
    },

    move: function(index) {
        var from = this.currentIndex;
    	this.currentIndex = index;
    	this.render(from);
        this.notifications.stop(true, true).hide();
    },

    movePrevious: function() {
        if (!this.at_first())
        	this.move(this.currentIndex - 1);
    },

    moveNext: function() {
        var step = this.steps[this.currentIndex];

        if (!this.at_last() && step.is_valid())
        	this.move(this.currentIndex + 1);
    },

    message: function(message) {
        this.notifications.html(message)
            .show()
            .stop(true, true)
            .delay(10*1000)
            .fadeOut(1500);
    },
    
//    check_login: function() { // TODO move to web services
//    	var logged_in = false;
//
//    	$.ajax({
//    		async: false,
//    		data: {
//    			fname: 'check_login',
//    		},
//    		success : function(rc) {
//    			logged_in = rc;
//    		}
//    	});
//
//    	return logged_in;
//    },

    submit: function() {
//        if (!this.check_login()) {
//            this.message('Your session has expired, please log in again.');
//            return;
//        }

        if (this.at_last())
            this.onCompleted(this.data);
    },

    // Expects a view with render, is_valid methods and a element property el
    addStep: function(step, index) {
        if (index !== undefined && index < this.steps.length)
            this.steps.slice(index, 0, step);
        else
            this.steps.push(step);
    }
});

function DataView(experiment, opts) {
    this.experiment = experiment || {};
    this.onError = opts.onError;
    this.supportedFileTypes = opts.supportedFileTypes;
    this.maxSraItems = opts.maxSraItems || 10;
    this.disableMaxItemsCheck = opts.disableMaxItemsCheck;
    this.title = "Select Data";
    this.files = [];
    this.initialize();
}

$.extend(DataView.prototype, {
    initialize: function() {
    	var self = this;
    	
        this.el = $($("#data-template").html());
        this.file_selector = $($("#fileselect-template").html());
        this.selector_container = this.el.find("#selector_container");
        this.file_table = this.el.find('#file_table');
        
        if (this.supportedFileTypes)
        	this.FILE_TYPE_PATTERNS = new RegExp("(:?" + this.supportedFileTypes.join("|") + ")$");

        coge.fileSelect.init({
        	container: this.selector_container,
        	fileTable: this.file_table,
        	defaultTab: DEFAULT_TAB,
        	maxIrodsListFiles: MAX_IRODS_LIST_FILES,
        	maxIrodsTransferFiles: MAX_IRODS_TRANSFER_FILES,
        	maxFtpFiles: MAX_FTP_FILES,
        	fileSelectSingle: FILE_SELECT_SINGLE,
        	loadId: LOAD_ID,
        	fileSelectedCallback: function() {
        		self.el.find('#files').show();
        	},
        	fileFinishedCallback: function() {
        	    var files = coge.fileSelect.get_selected_files();
        	    if (!files || !files.length)
        	    	return;
        	    
        	    // Detect file type and show in view
        	    if (files[0].type === 'ncbi' || files[0].type === 'sra') {
        	    	self.el.find("#select_file_type").hide();
        	    	if (files[0].type === 'sra') // this is a kludge, why does it work for ncbi?
        	    		self.el.find("#file_type_selector").val('sra');
        	    }
        	    else {
        	        var ft_sel = self.el.find("#file_type_selector");
	        	    var file_type = self.autodetect_file_type(files[0].name);
	        	    if (file_type) {
	        	        // mdb changed 11/21/16 -- add support for file types with more than one possible ext (such as FASTA, FAA, FA and FASTQ, FQ)
	        	        ft_sel.find("option").each(function() {
                            var val = $(this).val();
                            val.split(',').forEach(function(type) {
                                if (file_type == type)
                                    ft_sel.val(val);
                            });
	        	        });
                    }
	        	    self.el.find("#select_file_type").show();
        	    }
        	},
        	fileCancelledCallback: function() {
        		var files = coge.fileSelect.get_selected_files();
        		if (files) {
        	    	if (files[0].type == 'ncbi' || files[0].type == 'sra')
        	    		self.el.find("#select_file_type").hide();
        	    	self.el.find("#select_file_type option:first").attr("selected", "selected");
        	    }
        		else {
        	    	self.el.find('#files').hide();
        	    	self.el.find("#select_file_type").hide();
        	    }
        	}
        });
    },
    
    // Returns the file extension detected or undefined
    autodetect_file_type: function(file) {
        var stripped_file = file.replace(/.gz|.bz2$/, '');
        if (this.FILE_TYPE_PATTERNS.test(stripped_file))
            return stripped_file.match(this.FILE_TYPE_PATTERNS)[0];
    },

    render: function() {
        this.analysis_view = null; // hack for load experiment
        //FIXME: This selector should be in another view
        var selector = this.file_selector.clone();
        this.selector_container.empty();
        selector.appendTo(this.selector_container);

        coge.fileSelect.render();
    },

    is_valid: function() {
        var items = coge.fileSelect.get_selected_files();
        console.log(items);
        if (!items || items.length === 0) {
            if (this.onError)
            	this.onError('Please select a valid data file of type: ' + this.supportedFileTypes.join(', ') + '.  Allow uploaded files to finish transferring.');
            return false;
        }
        
        items[0].file_type = this.el.find("#select_file_type option:selected").val();
        if (!items[0].file_type) {
        	if (this.onError)
            	this.onError("Please select the file type to continue");
            return false;
        }
        
    	// Prevent mix of NCBI and file data types
		var types = {};
		items.forEach(function(item) { types[item.type] = 1; });
		var isNCBI = 'ncbi' in types;
		var isSRA = 'sra' in types;
		if (Object.keys(types).length > 1 && isNCBI) {
			this.onError('Cannot mix NCBI data with other types of data.');
			return false;
		}
		if (Object.keys(types).length > 1 && isSRA) {
			this.onError('Cannot mix SRA data with other types of data.');
			return false;
		}

		// Prevent too many SRA items
		var totalItems = items.reduce(
		    function(a, b) {
		        if (b.name && b.name.startsWith('SRR'))
		            return a+1;
		        return a+b.size;
            },
            0
        );

		if (isSRA && !this.disableMaxItemsCheck && totalItems > this.maxSraItems) {
            this.onError('Too many SRA items specified (' + totalItems + '), the limit is ' + this.maxSraItems + '.  If you need to load more than that, please contact the CoGe Team at <a href="mailto:coge.genome@gmail.com">coge.genome@gmail.com</a>.');
			return false;
		}

        this.experiment.data = items;
        return true;
    }
});

function ConfirmationView(object) { // object is an experiment or genome
    this.object = object;
    this.initialize();
    this.title = "Review & Submit";
}

$.extend(ConfirmationView.prototype, {
    initialize: function() {
        this.el = $($("#confirm-template").html());
        this.description = this.el.find(".confirm-description");
        this.options = this.el.find(".confirm-options");
        this.data = this.el.find(".confirm-data");
        this.pair_template = $('<div><span class="name"></span>: <span class="data"></span></div>');
    },

    render: function() {
        this.renderDescription(this.object.metadata);
        if (this.object.data)
        	this.renderData(this.object.data);
        if (this.object.options) 
        	this.renderOptions(this.object.options);
    },

    // Render description summary
    renderDescription: function(description) {
        this.description.empty();
        
        var key, newpair;
        for(key in description) {
            if (description.hasOwnProperty(key)) {
            	var value = description[key];
            	if (typeof value === "boolean")
            		value = (value ? "yes" : "no");
            	var newpair = this.pair_template.clone();
                newpair.find(".name").html(coge.utils.ucfirst(key));
                newpair.find(".data").html(value);
                this.description.append(newpair);
            }
        }
    },

    // Render data files summary
    renderData: function(data) {
        this.data.empty();

        var index, newpair;
        for(index = 0; index < data.length; index++) {
            newpair = this.pair_template.clone();
            newpair.find(".name").html("File");
            var filename = data[index].name.replace(/^.*[\\\/]/, '')
            newpair.find(".data").html(filename);
            this.data.append(newpair);
        }
    },

    // Render options summary
    renderOptions: function(options) {
        this.options.empty();

        var key, newpair;
        for(key in options) {
            if (options.hasOwnProperty(key)) {
            	var val = options[key];
            	if (typeof val === 'undefined' || val === '')
            		continue;
            	else if (typeof val === 'object')
            		val = coge.utils.objToString(val);
            	else if (typeof val === 'boolean')
            		val = val ? 'yes' : 'no';
            	else
                	val = String(val);
                newpair = this.pair_template.clone();
                newpair.find(".name").html(coge.utils.ucfirst(key.replace('_', ' ')));
                newpair.find(".data").html(val);
                this.options.append(newpair);
            }
        }
    },

    // Validates the confirmation view (nothing to do here)
    is_valid: function() {
        return true;
    }
});

function LayoutView(options) {
    this.template = $(options.template);
    this.layout = options.layout;
    this.initialize();
}

$.extend(LayoutView.prototype, {
    initialize: function() {
        this.el = $(this.template.html());
    },

    renderLayout: function() {
        var elementId, section, view;

        for (elementId in this.layout) {
            if (this.layout.hasOwnProperty(elementId)) {
                section = this.el.find(elementId);
                view = this.layout[elementId];

                section.empty();
                section.html(view.el);

                if (view.render) {
                    view.render();
                }
            }
        }
    },

    updateLayout: function(layout) {
        this.layout = $.extend({}, this.layout, layout);
    }
});
