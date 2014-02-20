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
            var form = document.createElement("form");
            form.method = "post";
            form.action = action;
            form.setAttribute("target", "_self");

            for(var key in params) {
                var input = document.createElement("input");
                input.name = key;
                input.value = params[key];

                form.appendChild(input);
            }

            form.submit("action");
        }
    };

    return ns;
})(coge || {});
