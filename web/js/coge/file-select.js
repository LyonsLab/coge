/* 
 * CoGe File Selector
 * 
 * Requires coge.utils, coge.services
 * 
 */

var coge = window.coge = (function(namespace) {
	// Methods
	namespace.fileSelect = {
		init: function(opts) {
			var self = this;
			
			// Set options
			self.container             = opts.container;
			self.fileTable             = opts.fileTable;
			self.defaultTab 		   = opts.defaultTab || 0;
			self.maxIrodsListFiles     = opts.maxIrodsListFiles || 100;
			self.maxIrodsTransferFiles = opts.maxIrodsTransferFiles || 50;
			self.maxFtpFiles           = opts.maxFtpFiles || 30;
			self.fileSelectSingle      = opts.fileSelectSingle || false;
			self.loadId                = opts.loadId;
			self.fileSelectedCallback  = opts.fileSelectedCallback;
			self.fileFinishedCallback  = opts.fileFinishedCallback;
			self.fileCancelledCallback = opts.fileCancelledCallback;
			self.selectedFiles = new Array;
			
			// Verify options
			if (!self.container) {
				console.error('FileSelect widget error: container not defined!');
				return;
			}
			// if (!self.fileTable) {
			// 	console.error('FileSelect widget error: file table not defined!');
			// 	return;
			// }
			// if (!self.loadId) {
			// 	console.error('FileSelect widget error: loadId not defined!');
			// 	return;
			// }
		},
		
		resize: function(event, ui) {
			//console.log('resizeStop ' + ui.size.width + ' ' + ui.size.height + ' ' + ui.originalSize.width + ' ' + ui.originalSize.height);
    		var heightChange = ui.size.height - ui.originalSize.height;
    		var panel = $("#ids_panel");
    		panel.height(panel.height() + heightChange);
		},
		
		render: function() {
			var self = this;
			
			// Set default tab
			self.container.tabs({selected: self.defaultTab});
			
			self.container.resizable({
				handles: 's',
				ghost: true,
				stop: self.resize.bind(self)
			});
			
			// Initialize irods view
			self._irods_get_path();
			
			// Setup file list event handlers
			self.selectedFiles.forEach(function(file) {
				file.tr.find(".ui-icon-closethick")
					.unbind()
					.click(self._cancel_callback.bind(self, file));
			});
			
			// Setup menu event handlers
			self.container.find('.fileselect-home').click(function() {
				self._irods_get_path();
			});
			self.container.find('.fileselect-shared').click(function() {
				self._irods_get_path("/iplant/home/shared");
			});
			self.container.find('.fileselect-up').click(function() {
				if ($('.fileselect-up').css('cursor') === 'pointer')
					self._irods_get_path("..");
			});
			self.container.find('.fileselect-refresh').click(function() {
				self._irods_get_path(".");
			});
			self.container.find('.fileselect-getall').click(function() {
				self._irods_get_all_files(".");
			});
			self.container.find('.fileselect-mkdir').click(function() {
				self._irods_mkdir();
			});
			
			self.container.find('.fileselect-filter').unbind().bind('keyup', function() {
				var search_term = self.container.find('.fileselect-filter').val();
				self.container.find('#ids_table tr td:nth-child(1)').each(function() {
					var obj = $(this);
					if (obj.text().toLowerCase().indexOf(search_term.toLowerCase()) >= 0)
						obj.parent().show();
					else
						obj.parent().hide();
				});
			});
			self.container.find('.fileselect-filter').bind('search', function() {
				var search_term = self.container.find('.fileselect-filter').val();
				if (!search_term.length)
					self.container.find('#ids_table tr td:nth-child(1)').each(function() {
						$(this).parent().show();
					});
			});
			
			self.container.find('#ftp_get_button').bind('click', function() {
				self._load_from_ftp();
			});
			
			self.container.find('#input_url').bind('keyup focus click', function() {
	        	var button = self.container.find("#ftp_get_button"),
	        		disabled = !self.container.find("#input_url").val();

	        	button.toggleClass("ui-state-disabled", disabled);
	        });

			self.container.find('#input_url').bind('keyup focus click', function() {
				if ( self.container.find('#input_url').val() ) {
					self.container.find('#ftp_get_button').removeClass('ui-state-disabled');
				}
				else {
					self.container.find('#ftp_get_button').addClass('ui-state-disabled');
				}
			});

			self.container.find('#input_accn').bind('keyup focus click', function() {
				if ( self.container.find('#input_accn').val() && self.container.find('#input_accn').val().length >= 6 ) {
					self.container.find('#ncbi_get_button,#sra_get_button').removeClass('ui-state-disabled');
				}
				else {
					self.container.find('#ncbi_get_button,#sra_get_button').addClass('ui-state-disabled');
				}
			});
			
			self.container.find('#ncbi_get_button').bind('click', function() {
				self._load_from_ncbi();
			});
			
			self.container.find('#sra_get_button').bind('click', function() {
				self._load_from_sra();
			});

			if (self.container.find('#input_upload_file').length)
				self.container.find('#input_upload_file').fileupload({
			    	dataType: 'json',
			    	add:
			    		function(e, data) {
			    			var filename = data.files[0].name;
							if ( !self._add_file_to_list(filename, 'file://'+filename) ) {
								//alert('File already exists.');
							}
							else {
								self.container.find('#input_upload_file').fileupload('option', { formData: {
						    		fname: 'upload_file',
						    		load_id: self.loadId
						    	}});

								data.submit(); // what is this?
							}
			    		},
					done:
						function(e, data) {
							self._finish_file_in_list('file', 'file://'+data.result.filename, data.result.path, data.result.size);
						}
				});
		},
		
		get_selected_files: function() {
			var files = this.selectedFiles.filter(function(file) {
				return ( (file.path || file.url) && !file.isTransferring && !file.error );
			});
			if (!files || files.length == 0)
				return;
			return files;
		},

		has_file: function(filename) {
			return this._filenames.indexOf(filename) != -1;
		},
		
		_clear_filter: function() {
			this.container.find('.fileselect-filter').val('');
			return this;
		},

		_clear_list: function() {
			timestamps['ftp'] = new Date().getTime(); // Cancel ftp transfers
			this.fileTable.html('').hide();
			$('#ftp_get_button').removeClass('ui-state-disabled');
			return this;
		},
		
		_add_file_to_list: function(filename, url, username, password) { // username/password only for FTP
			var self = this;
			
			if (self.fileTable) {
				// Skip if file already exists in file table
				if (self._find_file_in_list(url)) {
					console.warn('_fine_file_in_list: file already exists');
					return 0;
				}
				
				// Add file to view
				var tr = $('<tr class="note middle" style="height:1em;"><td style="padding-right:15px;">' +
						   '<span class="text">Name:</span> ' + filename +
						   '</td>' + '</tr>');
				var td = $('<td style="float:right;"><img src="picts/ajax-loader.gif"/></td>');
				$(tr).append(td).fadeIn();

				if (self.fileSelectSingle)
					self.fileTable.empty(); // remove all rows

				self.fileTable.append(tr).show();
			}
			
			// Add file to internal list
			var file = { 
				name: filename, 
				url: url, 
				tr: tr,
				isTransferring: true
			};
			if (username) file.username = username;
			if (password) file.password = password;
			self.selectedFiles.push(file);

			// Call template user's generic hander if defined
			if (self.fileSelectedCallback)
				self.fileSelectedCallback(file);

			return 1;
		},

		_vilify_file_in_list: function(url, error) {
			// Find item in list and set error
			var file = self._find_file_in_list(url);
			if (!file) {
				console.error('_vilify_file_in_list: file not found');
				return;
			}
			file.error = true;
			
			// Remove spinner
			file.tr.children().last().remove();
			
			// Update view of file
			var td1 = $('<td><span style="margin-right:15px;">' + (error ? error : 'failed') + '</span></td>').fadeIn();
			var closeIcon = $('<span class="link ui-icon ui-icon-closethick" title="Remove file"></span>');
			closeIcon.click(self._cancel_callback.bind(self, file));
			var td2 = $('<td></td>').append(closeIcon).fadeIn();

			$(file.tr).append(td1, td2)
				.css({
					'font-style': 'normal',
					'color': 'red'
				});
		},

		_finish_file_in_list: function(type, url, path, size) {
			var self = this;
			
			// Find item in list and update fields
			var file = self._find_file_in_list(url);
			if (!file) {
				console.error('_finish_file_in_list: file not found');
				return;
			}
			file.path = path;
			file.type = type;
			file.size = size;
			file.isTransferring = false;
			
			// Create view of file
			var td1 = $('<td></td>');
			if (size)
				td1.html('<span style="margin-right:15px;"><span class="text">Size:</span> ' + self._units(size) + '</span>');
			td1.fadeIn();
			var closeIcon = $('<span class="link ui-icon ui-icon-closethick" title="Remove file"></span>');
			closeIcon.click(self._cancel_callback.bind(self, file));
			var td2 = $('<td></td>').append(closeIcon).fadeIn();
			
			// Animate completion in list
			var tr = $(file.tr);
			tr.children().last().remove(); // remove spinner
			if (file.type != 'ncbi')
				tr.append(td1); // add size
			tr.append(td2); // add close icon
			tr.css({
				'font-style': 'normal',
				'color': 'black'
			});

			// Call template user's generic hander if defined
			if (self.fileFinishedCallback)
				self.fileFinishedCallback(file);
		},

		_cancel_callback: function(file) {
			var self = this;
			
			// Remove file from list
			for (var i = self.selectedFiles.length - 1; i >= 0; i--) {
			    if (self.selectedFiles[i].url === file.url) {
			    	self.selectedFiles.splice(i, 1);
			    }
			}
			
			// Remove file from view
			$(file.tr).hide('fast',
				function() {
					$(this).remove();

					// Call template user's hander if defined
					if (self.fileCancelledCallback)
						self.fileCancelledCallback(file);
				}
			);
		},

		_activate_on_input: function(input_id, button_id) {
			var a;
			if (input_id instanceof Array) 
				a = input_id;
			else {
				a = new Array(1);
				a[0] = input_id;
			}

			var hide = a.length;
			for (var i in a) {
				if ( $('#'+a[i]).val() ) 
					hide--;
			}

			if (hide)
				$('#'+button_id).addClass('ui-state-disabled');
			else
				$('#'+button_id).removeClass('ui-state-disabled');
		},

		_resolve_path: function(path) {
			if (typeof path === 'undefined' || path === 'undefined')
				return '';
			else if (path == '.')
				return this.current_path;
			else if (path == '..')
				return this.parent_path;
			return path;
		},
		
		_irods_busy: function(busy) {
			var status = this.container.find('#ids_status');
			if (typeof busy === 'undefined' || busy) 
				status.html('<img src="picts/ajax-loader.gif"/>').show();
			else
				status.html('<img src="picts/ajax-loader.gif"/>').hide();
			return this;
		},
		
		_irods_error: function(message) {
			var status = this.container.find('#ids_status');
			status.html('<span class="alert">'+message+'</span>').show();
		},

		_irods_get_path: function(path) {
			var self = this;
			
			self._irods_busy();

			path = self._resolve_path(path);

			coge.services.irods_list(path)
				.done(function(result) { //TODO move out into named function
					self._irods_busy(false)
						._clear_filter();
					
					var table = self.container.find('#ids_table');
	
					if (result == null) {
						self._irods_error('No result');
						return;
					}
					else if (result.error) {
						if (result.error.IRODS && result.error.IRODS === 'Access denied') {
							self._irods_error('Access denied');
							return;
						}
						table
							.html('<tr><td><span class="alert">'
							+ 'The following error occurred while accessing the Data Store.<br>'
							+ coge.utils.objToString(result.error) + '<br>'
							+ 'We apologize for the inconvenience.  Our support staff have already been notified and will resolve the issue ASAP. '
							+ 'If you just logged into CoGe for the first time, give the system a few minutes to setup your Data Store connection and try again.  '
							+ 'Please contact <a href="mailto:<TMPL_VAR NAME=SUPPORT_EMAIL>"><TMPL_VAR NAME=SUPPORT_EMAIL></a> with any questions or comments.'
							+'</span></td></tr>');
						return;
					}
					else if (result.items.length > self.maxIrodsListFiles) {
						alert("Too many files (" + result.items.length + ") to display, the limit is " + self.maxIrodsListFiles + ".");
						return;
					}
	
					table.html('');
					var parent_path = result.path.replace(/\/$/, '').split('/').slice(0,-1).join('/') + '/';
	
					// Save for later resolve_path()
					self.parent_path = parent_path;
					self.current_path = result.path;
	
					$('#ids_current_path').html(result.path);
					if (path === '')
						self.home_path = result.path;
					var p = result.path;
					if (p.charAt(p.length - 1) == '/')
						p = p.substring(0, p.length - 1);
					var home = self.home_path.substring(0, self.home_path.lastIndexOf('/'));
					$('.fileselect-up').css('opacity', p === home || p === '/iplant/home/shared' ? 0.4 : 1.0);
					$('.fileselect-up').css('cursor', p === home || p === '/iplant/home/shared' ? 'default' : 'pointer');
	
					self._filenames = [];
					if (result.items.length == 0)
						table.append('<tr><td style="padding-left:20px;font-style:italic;color:gray;">(empty)</td></tr>');
	
					result.items.forEach(
						function(obj) {
							// Build row in to be displayed
							var icon;
							if (obj.type == 'directory')
								icon = '<span class="ui-icon ui-icon-folder-collapsed"></span>';
							else if (obj.type == 'link') 
								icon = '<span class="ui-icon ui-icon-link"></span>';
							else // assume file type
								icon = '<span class="ui-icon ui-icon-document"></span>';
							tr = $('<tr class="'+ obj.type +'" style="white-space:nowrap;user-select:none;-webkit-user-select:none;-moz-user-select:none;"><td>' 
									+ icon
									+ decodeURIComponent(obj.name) + '</td><td>'
									+ (obj.size ? decodeURIComponent(obj.size) : '') + '</td><td>' 
									+ (obj.timestamp ? decodeURIComponent(obj.timestamp) : '') + '</td></tr>'); // mdb added decodeURI 8/14/14 issue 441
							if (obj.type == 'directory' || obj.type == 'link') {
								$(tr).click(
									function() {
										self._irods_get_path(obj.path);
									}
								);
								//tr.contextmenu( function(e) { return coge.fileSelect._delete_menu(e, obj, coge.fileSelect._irods_rmdir.bind(self)); } );
							}
							else {
								$(tr).click(
									function() {
										if ( self._add_file_to_list(decodeURIComponent(obj.name), 'irods://'+obj.path) )
											//self._irods_get_file(obj.path);
											self._finish_file_in_list('irods', 'irods://'+obj.path, obj.path, obj.size);
									}
								);
								//tr.contextmenu( function(e) { return coge.fileSelect._delete_menu(e, obj, coge.fileSelect._irods_rm.bind(self)); } );
								self._filenames.push(obj.name);
							}
	
							$(tr).hover(
								function() { $(this).css({"cursor":"pointer", "background-color":"greenyellow"}); },
								function() { $(this).css("background-color", "white"); }
							);
	
							table.append(tr);
						}
					);
					$('#ids_panel').scrollTop(0);
				})
				.fail(function() {
					// TODO
				});
		},

		_delete_menu: function(e, obj, func) {
			var tr = $(e.target).closest('tr');
			var m = $('#fileselect_menu');
			m.children().children().html('Delete ' + obj.type);
			var position = $(e.target).position();
			m.css('left', position.left + e.offsetX - 3);
			m.css('top', position.top + e.offsetY - 3);
			m.show();
			m.one('click', function() { func(obj); });
			m.one('contextmenu', function() { func(obj); return false; });
			m.one("mouseleave", function() {
				m.off('mouseover');
				m.hide();
				tr.css('background-color', 'white');
			});
			m.on('mouseover', function() {tr.css('background-color', 'greenyellow')});
			e.preventDefault();
			e.stopPropagation();
			return false;
		},

		_irods_get_all_files: function(path) {
			var self = this;
			
			self._irods_busy();
			
			path = self._resolve_path(path);
			
			coge.services.irods_list(path) 
				.done(function(result) {
					self._irods_busy(false);
					
					if (!result || !result.items)
						return;
					if (result.items.length > self.maxIrodsTransferFiles) {
						alert("Too many files (" + result.items.length + ") to retrieve at one time, the limit is " + self.maxIrodsTransferFiles + ".");
						return;
					}
	
					var count = 0;
					result.items.forEach(
						function(obj) {
							if (obj.type == 'file') {
								setTimeout(
									function() {
										if ( self._add_file_to_list(obj.name, 'irods://'+obj.path) ) {
											//self._irods_get_file(obj.path);
											self._finish_file_in_list('irods', 'irods://'+obj.path, obj.path, obj.size);
										}
									},
									100 * count++
								);
							}
						}
					);				
				})
				.fail(function() {
					//TODO
				});
		},

		_irods_mkdir: function() {
			var self = this;
			get_dirname(function(path){
				path = $('#ids_current_path').html() + '/' + path;
				coge.services.irods_mkdir(path).done(function(result) {
					self._irods_busy(false);
					if (result.error)
						alert(result.error.Error);
					else
						self._irods_get_path(path);
				});
			});
		},

		// _irods_rm: function(obj) {
		// 	var self = this;
		// 	this._confirm('Delete File', 'Really delete file ' + obj.name + '?', function() {
		// 		coge.services.irods_rm(obj.path).done(function(result) {
		// 			self._irods_busy(false);
		// 			if (result.error)
		// 				alert(result.error.Error);
		// 			else
		// 				self._irods_get_path($('#ids_current_path').html());
		// 		});
		// 	});
		// },

		// _irods_rmdir: function(obj) {
		// 	var self = this;
		// 	this._confirm('Delete Directory', 'Really delete directory ' + obj.name + ' and everything in it?', function() {
		// 		coge.services.irods_rm(obj.path).done(function(result) {
		// 			self._irods_busy(false);
		// 			if (result.error)
		// 				alert(result.error.Error);
		// 			else
		// 				self._irods_get_path($('#ids_current_path').html());
		// 		});
		// 	});
		// },

		_units: function(val) {
			if (isNaN(val))
				return val;
			else if (val < 1024)
				return val;
			else if (val < 1024*1024)
				return Math.ceil(val/1024) + 'K';
			else if (val < 1024*1024*1024)
				return Math.ceil(val/(1024*1024)) + 'M';
			else
				return Math.ceil(val/(1024*1024*1024)) + 'G';
		},
		
		_find_file_in_list: function(url) {
			return this.selectedFiles.filter(
				function(file) {
					return (url === file.url);
				}	
			)[0];
		},

// mdb removed 8/24/15 COGE-644 -- HTTP/FTP files are now transferred in workflow
//		_ftp_get_file: function(url, username, password) {
//			var self = this;
//			$.ajax({
//				data: {
//					fname: 'ftp_get_file',
//					load_id: self.loadId,
//					url: url,
//					username: username,
//					password: password,
//				},
//				success : function(data) {
//					var obj = jQuery.parseJSON(data); //FIXME change ajax type to "json" and remove this
//					if (!obj || obj.error) {
//						console.error('ftp_get_file error: ' + (obj.error ? obj.error : 'null'));
//						self._vilify_file_in_list(url, obj.error);
//						return;
//					}
//					
//					self._finish_file_in_list('ftp', url, obj.path, obj.size);
//				},
//			});
//		},

		_load_from_ftp: function() {
			var self = this;
			
			var url = $('#input_url').val();
			var username = $('#input_username').val();
			var password = $('#input_password').val();

			$('#ftp_get_button').addClass('ui-state-disabled');
			$('#ftp_status').html('<img src="picts/ajax-loader.gif"/> Contacting host...');

			coge.services.ftp_list(url)
				.done(function(result) {
					var filelist = result.items;
					if (!filelist || filelist.length == 0) {
						alert("Location not found.");
						$('#ftp_get_button').removeClass('ui-state-disabled');
						$('#ftp_status').html('');
						return;
					}
					if (filelist.length > self.maxFtpFiles) {
						alert("Too many files (" + filelist.length + ") at specified location, limit is " + self.maxFtpFiles + ".");
						$('#ftp_get_button').removeClass('ui-state-disabled');
						$('#ftp_status').html('');
						return;
					}

					self.filecount = filelist.length;
					$('#ftp_status').html('<img src="picts/ajax-loader.gif"/> Retrieving '+self.filecount+' files');

					var count = 0;
					filelist.forEach(
						function(obj) {
							setTimeout(
								function() {
									if (self._add_file_to_list(obj.name, obj.url, username, password)) {
										//self._ftp_get_file(obj.url, username, password); // mdb removed 8/24/15 COGE-644
										self._finish_file_in_list('ftp', obj.url, obj.path, obj.size); // mdb added 8/24/15 COGE-644
									}
									if (--self.filecount == 0) { // FTP transfer complete
										$('#ftp_get_button').removeClass('ui-state-disabled');
										$('#ftp_status').html('');
									}
								},
								250 * count++
							);
						}
					);
				})
				.fail(function() {
					//TODO
				});
		},
		
		_load_from_ncbi: function() {
			var self = this;
			
			var accn = $('#input_accn').val();

			$('#ncbi_get_button').addClass('ui-state-disabled');
			$('#ncbi_status').html('<img src="picts/ajax-loader.gif"/> Contacting NCBI Nucleotide DB...');

			$.ajax({
				data: {
					fname: 'search_ncbi_nucleotide',
					accn: accn,
					load_id: self.loadId,
				},
				success : function(data) {
					var obj = jQuery.parseJSON(data); //FIXME change ajax type to "json" and remove this
					if (obj) {
						if (obj.error) {
							$('#ncbi_status').html(obj.error);
							return;
						}
						else if (typeof obj.id != 'undefined') {
							if (self._add_file_to_list(obj.name, 'ncbi://'+obj.id)) {
								self._finish_file_in_list('ncbi', 'ncbi://'+obj.id, obj.id, '');
							}
						}
					}

					$('#ncbi_get_button').removeClass('ui-state-disabled');
					$('#ncbi_status').html('');
				},
			});
		},
		
		_load_from_sra: function() {
			var self = this;
			
			var accn = $('#input_accn').val();
			if (!(accn && accn.toLowerCase().startsWith('srr'))) {
				$('#sra_status').html('Please enter an SRR accession (starts with "SRR")');
				return;
			}

			$('#sra_status').html('<img src="picts/ajax-loader.gif"/> Contacting NCBI SRA...');
			
		    var entrez = new Entrez({ database: 'sra' });
		    entrez.search(accn).then(function(id) {
		    	if (id) {
		    		entrez.fetch(id).then(function(result) {
		    			if (result) {
		    				if (self._add_file_to_list(result.accn, 'sra://'+result.accn)) {
								self._finish_file_in_list('sra', 'sra://'+result.accn, result.accn, result.size);
							}
							$('#sra_status').html('');
		    			}
		    			else {
		    				$('#sra_status').html('Item not found');
		    			}
		    		});
		    	}
		    	else {
		    		$('#sra_status').html('Item not found');
		    	}
		    });
		},

		_confirm: function(title, question, on_ok) {
			$('<div></div>').appendTo('body')
			    .html('<div><h6>' + question + '</h6></div>')
			    .dialog({
				    modal: true,
				    title: title,
				    resizable: false,
				    buttons: {
				        OK: function () {
				            on_ok();
				            $(this).dialog("close");
				        },
				        Cancel: function () {
				            $(this).dialog("close");
				        }
				    },
				    close: function (event, ui) {
				        $(this).remove();
				    }
				});
		}
    };

    return namespace;
})(coge || {});