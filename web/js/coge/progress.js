/* 
 * CoGe Progress Dialog
 */

var coge = window.coge = (function(namespace) {
	// Methods
	namespace.progress = {
		init: function(opts) {
			var self = this;
			this.baseUrl = opts.baseUrl;
			this.userName = opts.userName;
			this.supportEmail = opts.supportEmail;
			this.success = opts.success;
			this.onError = opts.onError;
			this.onReset = opts.onReset;
			this.formatter = opts.formatter || this.default_formatter;
			
			var c = this.container = $('<div class="dialog_box progress"></div>');
			c.dialog({ 
				title: opts.title || 'Progress',
				autoOpen: false
			});
			
			c.dialog("widget").find('.ui-dialog-titlebar-close').hide(); // hide 'x' button
			
			var template = $($("#progress-template").html());
			this._render_template(template);
			
			// Setup resize handler
			var log = this.container.find(".log");
		    log.height( $( window ).height() * 0.5 );
		    c.dialog({
		    	modal: true,
		    	width: '60%',
		    	closeOnEscape: false,
		    	resizeStop: function(event, ui) {
		    		//console.log('resizeStop ' + ui.size.width + ' ' + ui.size.height + ' ' + ui.originalSize.width + ' ' + ui.originalSize.height);
		    		var widthChange = ui.size.width - ui.originalSize.width;
		    		var heightChange = ui.size.height - ui.originalSize.height;
		    		log.css({ width: log.width() + widthChange, height: log.height() + heightChange });
		    	}
		    });
		    
		    // Setup button handlers
		    c.find('.cancel,.ok').click( $.proxy(self.reset, self) );
		},
		
		reset: function() {
		    this.container.dialog('close');
		    if (this.onReset)
		    	this.onReset();
		},
		
		reset_log: function() {
			var c = this.container;
			c.find('.log,.progress-link').html('');
			c.find('.msg').show('');
		    c.find('.ok,.error,.finished,.done,.cancel,.logfile').hide();
		},
			
		begin: function(opts) {
			this.reset_log();
			
			if (opts && opts.title)
				this.container.dialog({title: opts.title});
			
			if (opts && opts.width)
				this.container.dialog({width: opts.width});
			
			if (opts && opts.height)
				this.container.find(".log").height(opts.height);
			
		    var log = this.container.find('.log');
		    if (opts && opts.content)
		    	log.html(opts.content);
		    else
		    	log.html('Initializing ...');
		    
		    this.container.dialog('open');
		    
		    this.startTime = new Date().getTime();
		},
		
		end: function() {
			this.container.dialog('close');
		},
		
		succeeded: function(string) {
			var c = this.container;
			
		    // Update dialog
		    c.find('.msg,.progress-link').hide();
		    c.find('.finished,.ok').fadeIn();
		    c.find('.log').append(string)
		    
		    // User callback
		    if (this.success)
		    	this.success();
		},

		failed: function(string) {
			var c = this.container;
			
			// Show error message
		    c.find('.log')
		    	.append('<div class="alert">' + string + '</div><br>')
		    	.append(
		    		'<div class="alert">' +
			        'The CoGe Support Team has been notified of this error but please ' +
			        'feel free to contact us at <a href="mailto:' + this.supportEmail + '">' +
			        this.supportEmail + '</a> ' +
			        'and we can help to determine the cause.' +
			        '</div>'
		    	);

		    // Show link to log file
		    var logfile = this.baseUrl + 'legacy/download?username=' + this.userName + '&' + 'wid=' + this.job_id;
		    $(".logfile a").attr("href", logfile);
		    $('.logfile').fadeIn();

		    // Update dialog
		    c.find('.msg,.progress-link').hide();
		    c.find('.error,.cancel').fadeIn();

// FIXME restore this email reporting
//		    if (newLoad) { // mdb added check to prevent redundant emails, 8/14/14 issue 458
//		        $.ajax({
//		            data: {
//		                fname: "send_error_report",
//		                load_id: load_id,
//		                job_id: WORKFLOW_ID
//		            }
//		        });
//		    }
		    
		    // User callback
		    if (this.onError)
		    	this.onError();
		},
		
		_refresh_interval: function() {
	        var interval = 2000;
	        
	        // Set refresh rate based on elapsed time
	        if (!this.startTime)
	        	this.startTime = new Date().getTime();
	        
	        var run_time = new Date().getTime() - this.startTime;
	        if (run_time > 10*60*1000)
	        	interval = 60*1000;
	        else if (run_time > 5*60*1000)
	        	interval = 30*1000;
	        else if (run_time > 60*1000)
	        	interval = 15*1000;
	        //console.log('Refresh run_time=' + run_time + ' refresh_interval=' + refresh_interval);
	        
	        return interval;
		},
		
		update: function(job_id, url) {
			var self = this;
			
			self._debug("update");
			
			if (job_id)
				self.job_id = job_id;
			if (url) {
				self.url = url;
				self.container.find('.progress-link').html('Link: <a href="'+url+'">'+url+'</a>').show();
			}
			
			var update_handler = function(json) {
				self._debug('update_handler');
				var c = self.container;
				
		        var workflow_status = $("<p></p>");
		        var log_content = $("<ul></ul>");
		        var results = [];

		        var refresh_interval = self._refresh_interval();
		        var retry_interval = 5*1000;
		        self._debug('refresh=' + refresh_interval + ', retry=' + retry_interval);
		        
		        // Sanity check -- progress dialog should be open
		        if (!c.dialog('isOpen')) {
		        	self._error('Error: progress dialog is closed');
		            return;
		        }

		        // Handle error
		        if (!json || json.error) {
		        	self.ajaxError++;
		            if ('Auth' in json.error) {
		            	c.find('.msg').html('Login required to continue');
		            	c.find('.log')
		            		.css({'font-size': '1em'})
		            		.html("<br>Your session has expired.<br><br>" + 
		            			"Please log in again by clicking " +
		            			"<a onclick='login_cas();' style='font-weight:bold'>here</a>.");
		            	return;
		            }
		            else {
			            self.alert('Server not responding ('+self.ajaxError+')');
			            setTimeout($.proxy(self.update, self), retry_interval);
			            return;
		            }
		        }
		        
		        self.ajaxError = 0;
		        self.alert();

		        // Retry on missing status (probably won't ever happen)
		        if (!json.status) {
		        	self._error('Error: missing status, retrying ...');
		        	setTimeout($.proxy(self.update, self), retry_interval);
		            return;
		        }
		        
		        // Render status
	            var current_status = json.status.toLowerCase();
	            self._debug('status=' + current_status);
	            workflow_status
	                .html("Workflow status: ")
	                .append( $('<span></span>').html(json.status) )
	                .addClass('bold');

	            // Retry on JEX error status -- mdb added 6/30/15
	            if (current_status == "error") {
	            	self.alert('JEX error status');
	            	setTimeout($.proxy(self.update, self), retry_interval);
	            	return;
	            }
	            
	            // Render tasks
		        if (json.tasks) {
		            var jobs = json.tasks;
		            for (var index = 0; index < jobs.length; index++) {
		                var item = self.formatter(jobs[index]);
		                if (item)
		                    results.push(item);
		            }
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
		            if (json.results && json.results.length) 
		            	self.succeeded(json.results);
		        }
		        else if (current_status == "failed"
		                //|| current_status == "error" // mdb removed 6/30/15 -- now handled first (see above code)
		                || current_status == "terminated"
		                || current_status == "cancelled")
		        {
		            workflow_status.find('span').addClass('alert');

		            if (json.results && json.results.length)
		                self.logfile = json.results[0].path;
		            self.failed();
		        }
		        else if (current_status == "notfound") {
		        	self._error('Error: status is "notfound"');
		        	setTimeout($.proxy(self.update, self), refresh_interval);
		            return;
		        }
		        else {
		            workflow_status.find('span').addClass('running');
		            setTimeout($.proxy(self.update, self), refresh_interval);
		        }

		        results.push(workflow_status);
		        log_content.append(results);
		        
		        if (json.results && json.results.length) {
		        	log_content.append("<div class='bold'>Here are the results (click to open):</div>");
		    	    json.results.forEach(function(result) {
		    	    	var html = self.format_result(result);
		    	    	log_content.append(html);
		    	    });
		        }
		        
		        self.container.find('.log').html(log_content);
		    };
			
		    setTimeout(
		    	function() { coge.services.fetch_job(self.job_id, update_handler, update_handler); },
		    	10
		    );
		},
		
		format_result: function(result) {
			if (result.type === 'experiment') {
				var url = 'ExperimentView.pl?eid=' + result.id;
				return "<div><a href='"+url+"'><img src='picts/testtube-icon.png' width='15' height='15'/> Experiment '"+result.name+"'</a></div>";
			}
			else if (result.type === 'notebook') {
				var url = 'NotebookView.pl?nid=' + result.id;
				return "<div><a href='"+url+"'><img src='picts/notebook-icon.png' width='15' height='15'/> Notebook '"+result.name+"'</a></div>";
			}
		},
		
		default_formatter: function(item) {
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

		    if (item.elapsed)
		        row.append(" in " + coge.utils.toPrettyDuration(item.elapsed));

		    if (item.log) {
		        var p = item.log.split("\n");

		        var pElements = p.map(function(item) {
		            var norm = item.replace(/\\t/g, " ").replace(/\\'/g, "'");
		            return $("<div></div>").html(norm);
		        });

		        var log = $("<div></div>").html(pElements).addClass("padded");
		        row.append(log);
		    }

		    return row;
		},
		
		alert: function(string) {
			var el = this.container.find(".alert");
			if (!string)
				el.html('').hide();
			else
				el.html(string).show();
		},
		
		_error: function(string) {
			console.error('progress: ' + string);
		},
		
		_debug: function(string) {
			console.log('progress: ' + string);
		},
		
		_render_template: function(template) {
		    this.container.empty()
		        .hide()
		        .append(template)
		        .show();//.slideDown();
		}
	
    };

    return namespace;
})(coge || {});