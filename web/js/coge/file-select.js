/* 
 * CoGe File Selector
 * 
 * Requires coge.utils, coge.services
 * 
 */
const SPINNER_HTML = '<img src="picts/ajax-loader-medium.gif">';
const SPINNER_SMALL_HTML = '<img src="picts/ajax-loader.gif">';

class Table {
	constructor(table, options) {
		this.table = table;
		this.options = options || {};
	}
	busy() {
		this.table.empty();
		this.row([ SPINNER_HTML ]);
	}
	empty() {
		this.table.empty();
	}
	filter(text, col) {
		text = text.toLowerCase();
		this.rows().forEach(function(row){
			if ($(row.children[col]).text().toLowerCase().indexOf(text) >= 0)
				$(row).show();
			else
				$(row).hide();
		});
	}
	row(cells, attributes) {
		let row = $('<tr></tr>', attributes).appendTo(this.table);
		if (this.options.hover)
			row.hover(this.options.hover[0], this.options.hover[1]);
		cells.forEach(function(cell) {
			if (typeof cell === 'object')
				cell.appendTo($('<td style="height:1.5em;"></td>').appendTo($(row)));
			else
				$('<td style="height:1.5em;"></td>').html(cell).appendTo($(row));
		});
	}
	rows() {
		return Array.from(this.table[0].firstChild.children);
	}
	show_all_rows() {
		this.rows().forEach(function(row){
			$(row).show();
		});
	}
}
// TODO: generalize this for use with irods, add toolbar construction here, remove from template
class FtpSelect {
    constructor(table, path_div, toolbar, on_file_click, options) {
		this.table = new Table(table, {
			hover: [function() { $(this).css({"cursor":"pointer", "background-color":"greenyellow"}); }, function() { $(this).css("background-color", "white"); }]
		});
		this.path_div = path_div;
		this.toolbar = toolbar;
        this.on_file_click = on_file_click;
        this.options = options || {};
    }
	at_top() {
		let index = this.url.indexOf('/', 6);
		return index == -1 || index == this.url.length - 1;
	}
	click_all_files() {
		let self = this;
		let rows = this.table.rows();
		this.items.forEach(function(item,index) {
			if (item.type == 'f' && $(rows[index]).is(":visible"))
				self.on_file_click(item, self.options.username, self.options.password);
		});
	}
	go(url) {
        this.url = url;
		if (this.options.on_go_start)
			this.options.on_go_start();
		if (this.path_div)
			this.path_div.html(url);
        this.table.busy();
        let self = this;
		let full_url = this.options.username ? 'ftp://' + this.options.username + ':' + this.options.password + '@' + url.substring(6) : url;
        coge.services.ftp_list(full_url, true).done(function(result) {
			self.table.empty();
            self.items = result.items;
            if (!self.items || self.items.length == 0) {
				if (self.toolbar)
					self.toolbar.css('display', 'none');
				if (self.path_div)
					self.path_div.empty();
                alert("Location not found.");
			} else if (self.items.length == 1 && !self.items[0].type) {
				if (self.toolbar)
					self.toolbar.css('display', 'none');
				if (self.path_div)
					self.path_div.empty();
				self.on_file_click(self.items[0], self.options.username, self.options.password);
			} else {
				if (self.toolbar)
					self.toolbar.css('display', 'block');
				self.items.forEach(function(item) {
					if (item.type == 'f')
						self.table.row([$('<span><span class="ui-icon ui-icon-document"></span>' + item.name + '</span>'), item.size, item.time], { click: function(){self.on_file_click(item, self.options.username, self.options.password); return false;} });
					else
						self.table.row([$('<span><span class="ui-icon ui-icon-folder-collapsed"></span>' + item.name + '/</span>')], { click: function(){self.go(item.url); return false;} });
				});
			}
			if (self.options.on_go_finish)
				self.options.on_go_finish();
        });
    }
	go_home() {
		if (!this.at_top())
			this.go(this.url.substring(0, this.url.indexOf('/', 6) + 1));
	}
	go_up() {
		if (!this.at_top())
			this.go(this.url.substring(0, this.url.lastIndexOf('/', this.url.length - 2) + 1));
	}
	refresh() {
		this.go(this.url);
	}
	set(name, value) {
		this.options[name] = value;
		return this;
	}
}
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

			self.initialized = true;
		},
		
		resize: function(event, ui) {
			//console.log('resizeStop ' + ui.size.width + ' ' + ui.size.height + ' ' + ui.originalSize.width + ' ' + ui.originalSize.height);
    		var heightChange = ui.size.height - ui.originalSize.height;
    		var panel = $("#ids_panel");
    		panel.height(panel.height() + heightChange);
		},
		
		render: function() {
			let self = this;

			// Set default tab
			self.container.tabs({selected: self.defaultTab});

            // Make resizable by dragging bottom border
            self.container.resizable({
				handles: 's',
				ghost: true,
				stop: self.resize.bind(self)
			});

			// Initialize irods view
            self._irods_get_path(self.irods_path);

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

			let fetch_start = function() {
				self.container.find('#ftp_get_button').addClass('ui-state-disabled');
				self.container.find('.ftpselect-getall').addClass('ui-state-disabled');
				self.container.find('.ftpselect-home').addClass('ui-state-disabled');
				self.container.find('.ftpselect-refresh').addClass('ui-state-disabled');
				self.container.find('.ftpselect-up').addClass('ui-state-disabled');
				self.container.find('#input_url').prop('disabled', true);
			};
			let fetch_finish = function() {
				self.container.find('#ftp_get_button').removeClass('ui-state-disabled');
				self.container.find('.ftpselect-getall').removeClass('ui-state-disabled');
				self.container.find('.ftpselect-home').removeClass('ui-state-disabled');
				self.container.find('.ftpselect-refresh').removeClass('ui-state-disabled');
				self.container.find('.ftpselect-up').removeClass('ui-state-disabled');
				self.container.find('#input_url').prop('disabled', false);
			}
			self.container.find('#ftp_get_button').bind('click', function() {
				if (!$(this).hasClass('ui-state-disabled')) {
					fetch_start();
					self._load_from_ftp();
				}
			});
			self.container.find('#input_url').bind('keyup focus click', function() {
	        	self.container.find("#ftp_get_button").toggleClass("ui-state-disabled", !self.container.find("#input_url").val());
	        });
			self.container.find('#input_url').bind('keyup focus click', function() {
	        	self.container.find("#ftp_get_button").toggleClass("ui-state-disabled", !self.container.find("#input_url").val());
			});
			self._ftp_select = new FtpSelect(self.container.find('#ftp_table'), self.container.find('#ftp_current_path'), self.container.find('#ftp_list_head'), this._add_ftp_file.bind(this), {
				on_go_start: function() {
					fetch_start();
				},
				on_go_finish: function() {
					$('.ftpselect-filter').val('');
					fetch_finish();
					if (self._ftp_select.at_top()) {
						self.container.find('.ftpselect-home').addClass('ui-state-disabled');
						self.container.find('.ftpselect-up').addClass('ui-state-disabled');
					} else {
						self.container.find('.ftpselect-home').removeClass('ui-state-disabled');
						self.container.find('.ftpselect-up').removeClass('ui-state-disabled');
					}
				}
			});
			self.container.find('.ftpselect-getall').click(function() {
				if (!$(this).hasClass('ui-state-disabled'))
					self._ftp_select.click_all_files();
			});
			self.container.find('.ftpselect-home').click(function() {
				if (!$(this).hasClass('ui-state-disabled')) {
					fetch_start();
					self._ftp_select.go_home();
				}
			});
			self.container.find('.ftpselect-refresh').click(function() {
				if (!$(this).hasClass('ui-state-disabled')) {
					fetch_start();
					self._ftp_select.refresh();
				}
			});
			self.container.find('.ftpselect-up').click(function() {
				if (!$(this).hasClass('ui-state-disabled')) {
					fetch_start();
					self._ftp_select.go_up();
				}
			});
			self.container.find('.ftpselect-filter').unbind().bind('keyup', function() {
				let search_term = self.container.find('.ftpselect-filter').val();
				self._ftp_select.table.filter(search_term, 0);
			});
			self.container.find('.ftpselect-filter').bind('search', function() {
				let search_term = self.container.find('.ftpselect-filter').val();
				if (!search_term.length)
					self._ftp_select.table.show_all_rows();
			});

			self.container.find('#input_accn').bind('keyup focus click', function() {
				if ( self.container.find('#input_accn').val() && self.container.find('#input_accn').val().length > 0 ) {
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

			// Setup file list event handlers
			self.selectedFiles.forEach(function(file) {
				file.tr.find(".ui-icon-closethick")
					.unbind()
					.click(self._cancel_callback.bind(self, file));
			});

			// Setup file upload
			if (self.container.find('#input_upload_file').length) {
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
            }
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
				var td = $('<td style="float:right;">'+SPINNER_SMALL_HTML+'</td>');
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

		_finish_file_in_list: function(type, url, path, size, unit) {
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
				td1.html('<span style="margin-right:15px;"><span class="text">Size:</span> ' + (unit ? size + ' ' + unit : self._units(size)) + '</span>');
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
			if (typeof busy === 'undefined' || busy)
				this.container.find('#ids_table').html('<tr><td>'+SPINNER_HTML+'</td></tr>').show();
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
			self.irods_path = path;

			coge.services.irods_list(path)
				.done(function(result) { //TODO move out into named function
					self._clear_filter();
					
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
						if (result.error.message && result.error.message.endsWith('rc=4')) {
							table.html('<tr><td><span class="alert">Your CyVerse account does not appear to be configured for CoGe.<br>If you have not already done so, please go to <a href="https://user.cyverse.org">https://user.cyverse.org</a> and click "Request Access" in the box for CoGe.<br>Please contact us at <a href="mailto:coge.genome@gmail.com">coge.genome@gmail.com</a> if you continue to encounter problems.</span></td></tr>');
							return;
						}
						table
							.html('<tr><td><span class="alert">'
							+ 'The following error occurred while accessing the Data Store.<br>'
							+ coge.utils.objToString(result.error) + '<br>'
							+ 'We apologize for the inconvenience.  Our support staff have already been notified and will resolve the issue ASAP. '
							+ 'If you just logged into CoGe for the first time, give the system a few minutes to setup your Data Store connection and try again.  '
							+ 'Please contact <a href="mailto:coge.genome@gmail.com">coge.genome@gmail.com</a> with any questions or comments.'
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
							tr = $('<tr class="'+ obj.type +'" style="white-space:nowrap;user-select:none;-webkit-user-select:none;-moz-user-select:none;">'
							        + '<td style="height:1.5em;">'
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
			
			path = self._resolve_path(path);
			
			coge.services.irods_list(path) 
				.done(function(result) {
					if (!result || !result.items)
						return;
					if (result.items.length > self.maxIrodsTransferFiles) {
						alert("Too many files (" + result.items.length + ") to retrieve at one time, the limit is " + self.maxIrodsTransferFiles + ".");
						return;
					}
	
					var count = 0;
					let rows = Array.from($('#ids_table')[0].firstChild.children);
					result.items.forEach(
						function(obj, index) {
							if (obj.type == 'file' && $(rows[index]).is(":visible")) {
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

		_add_ftp_file: function(obj, username, password) {
			if (this._add_file_to_list(obj.name, obj.url, username, password))
				this._finish_file_in_list('ftp', obj.url, obj.path, obj.size);
		},

		_load_from_ftp: function() {
			let self = this;
			
			let url = $('#input_url').val();
			if (!url.startsWith('ftp://') && !url.startsWith('http://')) {
				$('#error_help_text')
					.html('URL must begin with ftp:// or http://')
					.show()
            		.delay(10*1000)
            		.fadeOut(1500);;
				return;
			}
			let username = $('#input_username').val();
			let password = $('#input_password').val();

			$('#ftp_get_button').addClass('ui-state-disabled');

			if (url.startsWith('ftp://')) {
				this._ftp_select.set('username', username).set('password', password).go(url);
				return;
			}

			$('#ftp_status').html(SPINNER_SMALL_HTML+' Contacting host...');

			coge.services.ftp_list(url)
				.done(function(result) {
					var filelist = result.items;
					self.filecount = filelist.length;
					$('#ftp_status').html(SPINNER_SMALL_HTML+' Retrieving '+self.filecount+' file');

					var count = 0;
					filelist.forEach(
						function(obj) {
							setTimeout(
								function() {
									if (self._add_file_to_list(obj.name, obj.url, username, password)) {
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
			$('#ncbi_status').html(SPINNER_SMALL_HTML+'< Contacting NCBI Nucleotide DB...');

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

			// A list of space-separated SRA accession numbers
			var accn_list = $('#input_accn').val();
			if (!accn_list) return;
            accn_list.replace(/[^0-9a-z]/gi, ''); // remove all non-alphanumeric chars (except for underscore)

            // Setup Entrez
            var sra = new SRA();
            self._sra_busy();

            accn_list.split(/\s+/).forEach(function (accn) {
                if (!accn) return;

                // Search SRA for accession and fetch corresponding records
                if (sra.type(accn) == sra.SRA_ACCN_TYPE_EXPERIMENT) {
                    sra.esearch(accn).then(function(ids) {
                        if (ids && ids.length) {
                            $(ids).each(function(index, id) {
                                sra.esummary(id).then(function(xml) {
                                    var item = sra.extract(xml);
                                    item.runs.each(function(index, item) {
                                        if (self._add_file_to_list(item.accn, 'sra://'+item.accn))
                                            self._finish_file_in_list('sra', 'sra://'+item.accn, item.accn, item.size);
                                        self._sra_status('');
                                    });
                                });
                            });
                        }
                        else {
                            self._sra_status('Item not found: ' + accn);
                        }
                    });
                }
                else if (sra.type(accn) == sra.SRA_ACCN_TYPE_PROJECT) {
                    sra.esearch(accn).then(function(ids) {
                        if (ids) {
                            sra.esummary(ids[0]).then(function(xml) {
                                var item = sra.extract(xml);
                                var title = (item.title ? item.title : accn);
                                if (self._add_file_to_list(title, 'sra://'+accn))
                                    self._finish_file_in_list('sra', 'sra://'+accn, accn, ids.length, 'items');
                                self._sra_status('');
                            });
                        }
                    });
                }
                else {
                    self._sra_error('Unsupported accession type: ' + accn);
                    return;
                }
            });
		},

		_sra_error: function(text) {
		    $('#sra_status').html(text).removeClass('note').addClass('alert');
		},

		_sra_status: function(text) {
		    $('#sra_status').html(text).removeClass('alert').addClass('note');
		},

		_sra_busy: function() {
		    $('#sra_status').html(SPINNER_SMALL_HTML+' Contacting NCBI-SRA...').removeClass('alert').addClass('note');
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
