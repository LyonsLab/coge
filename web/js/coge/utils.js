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
            var copy = [], n = array.length, i = 0;

            while (n) {
                i = Math.floor(Math.random() * array.length);
                if (i in array) {
                    copy.push(array[i]);
                    delete array[i];
                    n--;
                }
            };

            return copy;
        }
    };

    return ns;
})(coge || {});
