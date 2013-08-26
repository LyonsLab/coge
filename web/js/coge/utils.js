var coge = (function(namespace) {
    namespace.ascending = function(a, b) {
        return a < b ? -1 : a > b ? 1 : 0;
    };

    namespace.descending = function(a, b) {
        return a > b ? -1 : a < b ? 1 : 0;
    };

    return namespace;
})(coge || {});
