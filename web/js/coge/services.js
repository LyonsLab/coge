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
				this.BASE_URL = opts.baseUrl;
			if (opts.userName)
				this.userName = opts.userName;
			this.debug = 1;
		},
			
		search_notebooks: function(search_term, user_name, success_callback, error_callback) {
			// TODO add param validation
			this._ajax("GET", BASE_URL + "notebooks/search/" + search_term + "/", null, success_callback);
		},
		
		submit_job: function(request, success_callback, error_callback) {
			this._ajax("PUT", BASE_URL + "jobs/", JSON.stringify(request), success_callback, error_callback);
		},

		fetch_job: function(id, success_callback, error_callback) {
			this._debug('fetch_job');
			if (!id) {
				this._error('fetch_job: missing id value');
				return;
			}
			this._ajax("GET", BASE_URL + "jobs/" + id, null, success_callback, error_callback);
		},
		
		irods_list: function(path, success_callback, error_callback) {
			this._ajax("GET", BASE_URL + "irods/list/" + path, null, success_callback, error_callback);
		},
		
		_ajax: function(type, url, data, success, error) {
			var self = this;
		    $.ajax({
		    	type: type,
		    	url: url + "?username=" + this.userName,
		    	dataType: "json",
		        contentType: "application/json",
		        xhrFields: {
	                withCredentials: true
	            },
		        data: data,
		        success: function(response) {
		            if (success)
		            	success(response);
		        },
		        error: function(jqXHR, textStatus, errorThrown) {
		        	self._error('ajax:error: ' + textStatus + (errorThrown ? errorThrown : ''));
		        	if (error)
		        		error(jqXHR, textStatus, errorThrown);
		        }
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