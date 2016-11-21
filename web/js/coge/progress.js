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
			this.onSuccess = opts.onSuccess;
			this.onError = opts.onError;
			this.onReset = opts.onReset;
			this.formatter = opts.formatter || this._default_formatter;
			this.buttonTemplate = opts.buttonTemplate;
			
			var c = this.container = $('<div class="dialog_box progress"></div>');
			c.dialog({ 
				title: opts.title || 'Progress',
				autoOpen: false
			});
			
			c.dialog("widget").find('.ui-dialog-titlebar-close').hide(); // hide 'x' button
			
			var template = $($("#progress-template").html());
			this._render_template(template);
			
			// Setup resize handler
			this.log = this.container.find(".log");
		    this.log.height( $( window ).height() * 0.5 );
		    c.dialog({
		    	modal: true,
		    	width: '60%',
		    	closeOnEscape: false,
		    	resizeStop: function(event, ui) {
		    		//console.log('resizeStop ' + ui.size.width + ' ' + ui.size.height + ' ' + ui.originalSize.width + ' ' + ui.originalSize.height);
		    		var widthChange = ui.size.width - ui.originalSize.width;
		    		var heightChange = ui.size.height - ui.originalSize.height;
		    		self.log.css({ width: self.log.width() + widthChange, height: self.log.height() + heightChange });
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
		
		_reset_log: function() {
			var c = this.container;
			c.find('.log,.progress-link').html('');
			c.find('.status').html('').show();
		    c.find('.ok,.cancel,.logfile,.buttons').hide();
		},
			
		begin: function(opts) {
			this._reset_log();

			if (!opts) opts = {};

			if (opts.title)
				this.container.dialog({title: opts.title});
			
			if (opts.width)
				this.container.dialog({width: opts.width});
			
			if (opts.height)
				this.log.height(opts.height);
			
		    if (opts.content)
		    	this.log.html(opts.content);
		    else
		    	this.log.html('Initializing ...');
		    
		    this.container.dialog('open');
		    
		    this.startTime = new Date().getTime();
		},
		
		end: function() {
			this.container.dialog('close');
		},
		
		succeeded: function(response) {
			var c = this.container;
			
		    // Update dialog
            c.find('.progress-link').hide();

		    // Show user-specified button template
		    if (this.buttonTemplate) {
		    	var template = $($("#"+this.buttonTemplate).html());
		    	c.find('.buttons').html(template).show();
		    }
		    else { // default buttons
		    	c.find('.ok').show();
		    }

            // Compute total duration of all tasks // FIXME replace with workflow elapsed time since this miscalculates concurrent tasks
            var totalDuration = response.tasks.reduce(function(a, b) {
                if (!b.elapsed) return a;
                return a + b.elapsed;
            }, 0);

		    this._set_status(response.status).append('<span class="info">' +
		        (totalDuration ? ' in ' + coge.utils.toPrettyDuration(totalDuration) : '') + '</span>');
		    
		    // User callback
		    if (this.onSuccess)
		    	this.onSuccess(response.results);
		},
		
		failed: function(string, error) {
			var c = this.container;
			
			var errorMsg = string + (error ? ': ' + this._errorToString(error) : '');
			
			// Show error message
		    this.log
		    	.append('<div class="alert">' + (errorMsg ? errorMsg : '') + '</div><br>')
		    	.append(
		    		'<div class="alert">' +
			        'The CoGe Support Team has been notified of this error but please ' +
			        'feel free to contact us at <a href="mailto:' + this.supportEmail + '">' +
			        this.supportEmail + '</a> ' +
			        'and we can help to determine the cause.' +
			        '</div>'
		    	);

		    // Show link to log file
		    var logfile = coge.services.download_url({ wid: this.job_id, attachment: 0 });
		    $(".logfile a").attr("href", logfile);
		    $('.logfile').show();

		    // Update dialog
		    c.find('.progress-link,.buttons').hide();
		    this._set_status('Failed');
		    c.find('.cancel').show();

// FIXME restore email reporting
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

		cancelled: function(response) {
		    // Update dialog
		    this.container.find('.progress-link,.buttons').hide();
		    this._set_status(response.status);
		    this.container.find('.cancel').show();
		},

        _errorToString: function(error) {
			var string = '';
			for (key in error) {
				string += error[key];
			}
			return string;
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
			
			var update_handler = function(response) {
				self._debug('update_handler');
				var c = self.container;
				
		        var refresh_interval = self._refresh_interval();
		        var retry_interval = 5*1000;
		        self._debug('refresh=' + refresh_interval + ', retry=' + retry_interval);
		        
		        // Sanity check -- progress dialog should be open
		        if (!c.dialog('isOpen')) {
		        	self._error('Error: progress dialog is closed');
		            return;
		        }

		        // Handle error
		        if (!response || response.error) {
		        	self.ajaxError++;
		            if ('Auth' in response.error) {
		            	c.find('.msg').html('Login required to continue');
		            	this.log
		            		.css({'font-size': '1em'})
		            		.html("<br>Your session has expired.<br><br>" + 
		            			"Please log in again by clicking " +
		            			"<a onclick='login_cas();' style='font-weight:bold'>here</a>.");
		            	return;
		            }
		            else {
			            self._alert('Server not responding ('+self.ajaxError+')');
			            setTimeout($.proxy(self.update, self), retry_interval);
			            return;
		            }
		        }
		        
		        self.ajaxError = 0;
		        self._alert();

		        // Retry on missing status (probably won't ever happen)
		        if (!response.status) {
		        	self._error('Error: missing status, retrying ...');
		        	setTimeout($.proxy(self.update, self), retry_interval);
		            return;
		        }
		        
		        // Render status
	            var current_status = response.status.toLowerCase();
	            self._debug('status=' + current_status);

	            // Retry on JEX error status -- mdb added 6/30/15
	            if (current_status == "error") {
	            	self._alert('JEX error status');
	            	setTimeout($.proxy(self.update, self), retry_interval);
	            	return;
	            }
	            
	            // Render tasks status
	            var log_content = $("<div/>");
		        if (response.tasks) {
                    var html = self.formatter(response.tasks);
                    if (html)
                        log_content = $(html);
		        }

                // Handle workflow status
                switch (current_status.toLowerCase()) {
                    case 'completed' :
                        self.succeeded(response);
                        break;
                    case 'failed' :
                    case 'terminated' :
                        if (response.results && response.results.length)
		                    self.logfile = response.results[0].path;
		                self.failed(response);
		                break;
                    case 'cancelled' :
                    case 'stopped' :
                        self.cancelled(response);
                        break;
                    case 'notfound' :
                        self._error('Error: status is "notfound"');
		        	    self.log.html('Error: status is "not found" ... retrying');
		        	    setTimeout($.proxy(self.update, self), refresh_interval);
		        	    return;
                    case 'running' :
                        setTimeout($.proxy(self.update, self), refresh_interval);
		                self._set_status( $('<span><img class="top" src="picts/ajax-loader.gif"/>&nbsp;&nbsp;</span>') ).append( self._format_status(response.status) );
                        break;
                    default :
                        self._error('Error: status is invalid');
		        	    self.log.html('Error: status is invalid ... retrying');
		        	    setTimeout($.proxy(self.update, self), refresh_interval);
		        	    return;
                }

		        // Render workflow results
		        if (response.results && response.results.length > 2) { // Ignore first two results (debug.log and workflow.log)
		        	log_content.append("<br><div class='header'>Results (as they are generated)</div>");
		    	    response.results.forEach(function(result) {
		    	    	log_content.append( self._format_result(result) );
		    	    });
		        }

                // Insert rendered log
                log_content.append('<br><br>');
		        self.log.html(log_content);

                // Scroll to bottom of log if finished to ensure results are shown
		        if (current_status == "completed")
		            self.log.scrollTop(self.log[0].scrollHeight);
		    };
			
		    setTimeout(
		    	function() { coge.services.fetch_job(self.job_id).always(update_handler); },
		    	10
		    );
		},

		_format_result: function(result) {
		    var formatted = $('<div></div>').addClass('progress result');

			if (result.type === 'experiment') {
				var url = 'ExperimentView.pl?eid=' + result.id;
				formatted.append("<a href='"+url+"'><img src='picts/testtube-icon.png' width='15' height='15'/> Experiment '"+result.name+"'</a>");
			}
			else if (result.type === 'notebook') {
				var url = 'NotebookView.pl?nid=' + result.id;
				formatted.append("<a href='"+url+"'><img src='picts/notebook-icon.png' width='15' height='15'/> Notebook '"+result.name+"'</a>");
			}
			else if (result.type === 'genome') {
				var url = 'GenomeInfo.pl?gid=' + result.id;
				formatted.append("<a href='"+url+"'><img src='picts/dna-icon.png' width='15' height='15'/> Genome '"+result.info+"'</a>");
			}
			else if (result.type === 'dataset') {
				var url = 'GenomeInfo.pl?gid=' + result.genome_id;
				formatted.append("<a href='"+url+"'> Dataset '"+result.info+"'</a>");
			}
			else if (result.type === 'popgen') {
				var url = 'PopGen.pl?eid=' + result.experiment_id;
				formatted.append("<a href='"+url+"'>"+result.name+"</a>");
			}
			else {
			    return;
			}

			return formatted;
		},

		_set_status: function(status) {
		    if (typeof status === 'string' || status instanceof String) {
                status = this._format_status(status);
		    }
		    else if (!status instanceof jQuery)
		        return;

		    return this.container.find('.status').html(status).show();
		},

		_format_status: function(status) {
            var el = $('<span></span>');

            switch (status.toLowerCase()) {
                case 'scheduled' : el.append(status).addClass('down bold');                 break;
                case 'completed' : el.append(status).addClass('completed bold');            break;
                case 'running'   : el.append(status).addClass('running bold');              break;
                case 'skipped'   : el.append("already generated").addClass('skipped bold'); break;
                case 'cancelled' :
                case 'stopped'   :
                case 'failed'    : el.append(status).addClass('alert bold');                break;
            }

            return el;
		},
		
		_default_formatter: function(tasks) {
		    var self = this;

            var max_col2_width = 170;
            var max_col1_width = self.container.width() - max_col2_width;
		    var table = $('<table></table>').css('width', '100%');

            tasks.forEach(function(task) {
                var row = $('<tr></tr>').append(
                    $('<td>' + task.description + '</td>')
                        .css({'width': max_col1_width, 'max-width': max_col1_width, 'white-space': 'nowrap', 'overflow': 'hidden', 'text-overflow': 'ellipsis'})
                );

                var status = self._format_status(task.status);

                var duration = $('<span></span>');
                if (task.elapsed)
                    duration.append(' in ' + coge.utils.toPrettyDuration(task.elapsed));

                row.append( $('<td></td>').append(status).css({'width': max_col2_width, 'max-width': max_col2_width, 'white-space': 'nowrap', 'text-align' : 'right'}).append(duration) );

                table.append(row);

                if (task.log) {
                    var cell = $('<td></td>').attr('colspan', 2);

                    var p = task.log.split("\n");
                    p.forEach(function(line) {
                        var norm = line.replace(/\\t/g, " ").replace(/\\'/g, "'");
                        cell.append('<div>'+norm+'</div>').addClass("indent");
                    });

                    table.append($('<tr></tr>').append(cell));
                }
            });

		    return table.html();
		},
		
		_alert: function(string) {
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
		        .show();
		}
	
    };

    return namespace;
})(coge || {});