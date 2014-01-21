var coge = window.coge = (function(ns) {
    ns.progress = function() {
        var width = 500,
            title = "Running Analysis...",
            status_message = $("<div></div>"),
            data = $("<div></div>"),
            button = $("<span></span>"),
            link = $("<a></a>"),
            img = $("<img></img>");

        data.addClass("ui-corner-all")
            .addClass("ui-widget-content")
            .addClass("padded")
            .addClass("coge-progress-data");

        link.addClass("ui-button")
            .addClass("ui-corner-all")
            .addClass("ui-button-stop")
            .addClass("r")

        button.addClass("ui-button")
            .addClass("ui-corner-all")
            .addClass("ui-button-stop")
            .addClass("r")
            .html("Close")
            .on("click", function() { my.close(); });

        status_message.append()
            .addClass("coge-progress-status");

        var dialog = $("<div></div>").append(data)
            .append(status_message)
            .dialog({
                autoOpen: false,
                title: title,
                width: width,
                modal: true,
                dialogClass: "coge-progress-menu"
            });

        function my() {
        }

        my.element = function() {
            return dialog;
        }

        my.close = function() {
            dialog.dialog("close");
            return my;
        };

        my.status = function(response) {
            var message = $("<span></span>")
                .css("margin-right", "2%");

            status_message.hide()
                .html("")

            if (response.state === "complete") {
                img.attr("src", "./picts/thumbs_up.png");
                message.addClass("coge-progress-status-complete")
                    .html(response.message);

                status_message
                    .append(message)
                    .append(img);

            } else if (response.state === "error") {
                img.attr("src", "./picts/thumbs_down.png");
                message.addClass("coge-progress-status-error")
                    .html(response.message);
                status_message
                    .append(message)
                    .append(img);
            } else {
                img.attr("src", "./picts/ajax-loader.gif");
                message.html(response.message);
                status_message
                    .append(message)
                    .append(img);
            }

            if (response.button) {
                status_message.append(button);

                //NOTE: event needs to bound after being inserted into the dom
                button.unbind().on("click", my.close);
            }

            if (response.link !== undefined) {
                link.attr("href", response.link.href)
                    .attr("target", "_blank")
                    .css("margin-right", "2%")
                    .css("color", "black")
                    .html(response.link.message);

                status_message.append(link)
            }

            if (response.delay) {
                status_message.fadeIn(1000);
            } else {
                status_message.show();
            }

            return my;
        };

        my.button = function(_) {
            if (!arguments.length) return button;
            button.text(_);
            return my;
        }

        my.show = function() {
            dialog.dialog("open");
            return my;
        };

        my.title = function(_) {
            if (!arguments.length) return title;

            title = _;
            dialog.dialog("option", "title", title);
            return my;
        };

        my.data = function(_) {
            if (!arguments.length) return data;

            data.html(_);
            return my;
        };

        return my;
    }

    return ns;
})(coge || {});
