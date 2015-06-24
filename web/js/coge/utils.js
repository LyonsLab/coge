/* global window, document, coge*/
var coge = window.coge = (function(ns) {
    ns.utils = {
        ascending: function(a, b) {
            return a < b ? -1 : a > b ? 1 : 0;
        },

        descending: function(a, b) {
            return a > b ? -1 : a < b ? 1 : 0;
        },

        open: function(link) {
            window.open(link, "_self");
        },
        
        shuffle: function(array) {
            var copy = array.slice(0), n = array.length, tmp, i;

            while (n) {
                i = Math.floor(Math.random() * n--);
                tmp = copy[n];
                copy[n] = copy[i];
                copy[i] = tmp;
            };

            return copy;
        },
        
        post: function(action, params) {
            var key,
                input,
                form = document.createElement("form");

            form.method = "post";
            form.action = action;
            form.setAttribute("target", "_blank");
            form.setAttribute("enctype", "multipart/form-data");

            for(key in params) {
                if (params.hasOwnProperty(key)) {
                    input = document.createElement("textarea");
                    input.name = key;
                    input.value = params[key];
                    form.appendChild(input);
                }
            }

            form.submit("action");
        },
        
        toPrettyDuration: function(seconds) {
            var fields = [
                [parseInt((seconds / 86400).toFixed(1), 10), " day"],
                [parseInt((seconds / 3600).toFixed(1), 10) % 24, " hour"],
                [parseInt(((seconds / 60) % 60).toFixed(1), 10), " minute"],
                [(seconds % 60).toFixed(2), " second"]
            ];

            return fields.filter(function(item) { return item[0] !== 0; })
                .map(function(item) {
                    var word = (item[0] > 1) ? item[1] + "s" : item[1];
                    return item[0] + word;
                })
                .join(", ");
        },

        log10: function(value) {
            return Math.log(value) / Math.log(10);
        },
        
        ucfirst: function(string) {
            return string.charAt(0).toUpperCase() + string.slice(1);
        },
        
        getURLParameters: function () { // returns object of query params
            var match,
                pl     = /\+/g,  // Regex for replacing addition symbol with a space
                search = /([^&=]+)=?([^&]*)/g,
                decode = function (s) { return decodeURIComponent(s.replace(pl, " ")); },
                query  = window.location.search.substring(1);

            var urlParams = {};
            while (match = search.exec(query))
               urlParams[decode(match[1])] = decode(match[2]);
            return urlParams;
        },
        
        toQueryString: function(obj, prefix) { // accepts an object of query params
		    var str = [];
		    for (var p in obj) {
		        if (obj.hasOwnProperty(p)) {
		        	var k = prefix ? prefix + "[" + p + "]" : p, v = obj[p];
		        	str.push(typeof v == "object" ?
		        			serialize(v, k) :
		        			encodeURIComponent(k) + "=" + encodeURIComponent(v));
		        }
		    }
		    return str.join("&");
    	}
    };

    return ns;
})(coge || {});
