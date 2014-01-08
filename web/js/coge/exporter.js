var coge = window.coge = (function (ns) {
    ns.exporter = function(selector, data, cb) {
            this.element = $(selector);
            this.data = data;

            this.export = function() {
                var startTime = new Date().getTime();
                var message = this.element.dialog("open");

                var close = function() {
                    message.dialog("close");
                };

                $.ajax({
                    data: this.data,
                    success: function() {
                        setTimeout(close, Math.max(0, 5000-(new Date().getTime()-startTime)));
                        if (cb !== undefined) return cb(arguments);
                    }
                });
            };

        };
    return ns;
})(coge || {});
