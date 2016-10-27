var coge = window.coge = (function (namespace) {
	var DEFAULT_INTERVAL = 5*1000; // ms
	
    namespace.logout = {

        init: function(params) {
        	if (!Cookies) {
        		console.error('logout.init: Cookies module is not loaded!');
        		return;
        	}
        	
        	if (!params) {
        		console.error('logout.init: Missing required parameter object');
        		return;
        	}
        	
        	if (params.loginCookieName) {
        		this.loginCookieName = params.loginCookieName;
        		console.log("logout.init: loginCookieName='" + this.loginCookieName + "'");
        	}
        	else {
        		console.warn('logout.init: Disabling logout detection due to missing parameter loginCookieName (this is not necessarily an error)');
        		return;
        	}
        	
        	if (params.loginCallback) {
        		this.loginCallback = params.loginCallback;
        	}
        	else {
        		console.error('logout.init: Missing required parameter: loginCallback');
        		return;
        	}
        	
        	this.interval = DEFAULT_INTERVAL;
        	if (params.interval)
        		this.interal = params.interval;
        	
        	// Start logout check timer if logged in
        	var isLoggedIn = this.isLoggedIn();
        	console.log('logout.init: logged in = ' + isLoggedIn)
        	if (isLoggedIn)
        		this.start();
        },
        
        start: function() {
        	this.timer = setInterval(this.checkLogin.bind(this), this.interval);
        },
        
        stop: function() {
        	clearInterval(this.timer);
        },
        
        isLoggedIn: function() {
        	var cookie = Cookies.get(this.loginCookieName);
    		if (!cookie)
    			return false;
    		return true;
        },
        
        checkLogin: function() {
        	var self = this;
        	
        	// Get login status
        	var isLoggedIn = this.isLoggedIn();
        	if (this.debug) console.log('logout.checkLogin: logged in = ', isLoggedIn);
        	if (isLoggedIn)
        		return;
        	
        	// Stop the check timer
    		this.stop();
    		
    		// Show alert dialog
    		var dialog = $("<div />") // TODO move rendering into separate function
    			.html("Your login session has expired. Please login or continue as a public user.")
    			.dialog({
    				title: "Your session has expired ...",
	    			dialogClass: "no-close",
	    			closeOnEscape: false,
	    			height: 'auto',
	    			modal: true,
	    	        resizable: false,
	    			buttons: [
		    			{
			    			text: "Login",
			    			click: function() {
			    				self.loginCallback();
			    			}
		    			},
		    			{
			    			text: "Continue",
			    			click: function() {
			    				window.location.reload(false);
			    			}
		    			}
	    			],
	    			open: function(event, ui) { // hide close 'X' button
	    		        $(".ui-dialog-titlebar-close", ui.dialog | ui).hide();
	    		    }
        		});
        },
        
        setDebug: function(enable) {
        	this.debug = enable;
        }

    };
    
    return namespace;
})(coge || {});
