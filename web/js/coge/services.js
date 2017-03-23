/* 
 * CoGe Web Services API
 */

var coge = window.coge = (function(namespace) {
	// Local constants
	var BASE_URL = "api/v1/";
	
	// Methods
	namespace.services = {
		init: function(opts) {
			if (opts.baseUrl)
				this.baseUrl = opts.baseUrl;
			else 
				this.baseUrl = BASE_URL;
			
			if (opts.userName)
				this.userName = opts.userName;
			
			this.debug = opts.debug; // display console debug messages
			this.debug = 1;
			
			this._debug('init completed');
		},
		
		download_url: function(opts) { // duplicated in Web.pm::download_url_for()
			var gid = opts.gid;
			var eid = opts.eid;
			var wid = opts.wid;
			var filename = opts.filename;
			var attachment = opts.attachment;

		    var params = '';
		    if (gid)
		        params += "gid=" + gid;
		    else if (eid)
		        params += "eid=" + eid;
		    else if (wid)
		        params += "wid=" + wid;

		    if (attachment != null)
		        params += "&attachment=" + attachment;
		    
		    if (filename)
		    	params += "&filename=" + filename;
		    
		    return this.baseUrl + "downloads/?" + params;			
		},
		
		search_global: function(search_term) {
			return this._ajax("GET", this.baseUrl + "global/search/" + encodeURIComponent(search_term) + "/");
		},
		
		search_organisms: function(search_term) {
			// TODO add param validation
			return this._ajax("GET", this.baseUrl + "organisms/search/" + search_term + "/");
		},

		add_organism: function(request) {
			return this._ajax("PUT", this.baseUrl + "organisms/", request);
		},
		
		search_genomes: function(search_term, opts) {
			// TODO add param validation
			return this._ajax("GET", this.baseUrl + "genomes/search/" + search_term + "/", null, opts);
		},
			
		search_notebooks: function(search_term) {
			// TODO add param validation
			return this._ajax("GET", this.baseUrl + "notebooks/search/" + search_term + "/");
		},

		search_users: function(search_term) {
			// TODO add param validation
			return this._ajax("GET", this.baseUrl + "users/search/" + search_term + "/");
		},
		
		submit_job: function(request) {
			return this._ajax("PUT", this.baseUrl + "jobs/", request);
		},

		fetch_job: function(id) {
			this._debug('fetch_job');
			if (!id) {
				this._error('fetch_job: missing id value');
				return;
			}
			return this._ajax("GET", this.baseUrl + "jobs/" + id);
		},
		
		fetch_logs: function(id, type) {
			this._debug('fetch_logs');
			if (!id || !type) {
				this._error('fetch_job: missing id/type value');
				return;
			}
			return this._ajax("GET", this.baseUrl + "logs/" + type + '/' + id);
		},
		
		fetch_notebook: function(id) {
		    return this._ajax("GET", this.baseUrl + "notebooks/" + id + "/");
		},

		ftp_list: function(url, dirs) {
			return this._ajax("GET", this.baseUrl + "ftp/list/", null, { url: url, dirs: dirs });
		},
		
		irods_list: function(path) {
			return this._ajax("GET", this.baseUrl + "irods/list/" + path);
		},
		
		irods_mkdir: function(path) {
			return this._ajax("POST", this.baseUrl + "irods/mkdir", null, { path: path });
		},
		
		// irods_rm: function(path) {
		// 	return this._ajax("DELETE", this.baseUrl + "irods/rm", null, { path: path });
		// },

		update: function(type, id, data) {
			return this._ajax("POST", this.baseUrl + type + 's/' + id, data);
		},
		
		_ajax: function(type, url, data, opts) { //, success, error) {
			var self = this;
			
			// Build param list
			var opts_str = '';
			if (this.userName)
				opts_str += "username=" + this.userName;
			for (key in opts) {
				if (opts_str)
					opts_str += "&";
				opts_str += key + "=" + opts[key];
			}
			url = url + (opts_str ? "?" + opts_str : '');
			
			// Return deferred
			if (this.debug) {
				this._debug('coge.services._ajax: ' + url);
				if (data)
					console.log(data);
			}
		    return $.ajax({
				    	type: type,
				    	url: url,
				    	dataType: "json",
				        contentType: "application/json",
				        // xhrFields: {
			         //        withCredentials: true
			         //    },
				        data: JSON.stringify(data),
				    })
				    .fail(function(jqXHR, textStatus, errorThrown) { 
				    	self._error('ajax error: status="' + jqXHR.status + '" statusText="' + jqXHR.statusText + '" url=' + url); 
				    });		
		},
		
		_error: function(string) {
			console.error('services: ' + string);
		},
		
		_debug: function(string) {
			if (this.debug)
				console.log('services: ' + string);
		}
	
    };

	// TODO add auth code
	
    return namespace;
})(coge || {});