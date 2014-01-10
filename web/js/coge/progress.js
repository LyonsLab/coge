window.coge = window.coge || {};

coge.progress = (function() {
    function progress(selector) {
        var selection = selector,
            width = 500,
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

        button.addClass("ui-button")
            .addClass("ui-corner-all")
            .addClass("ui-button-go")
            .addClass("r")
            .html("Close")
            .on("click", my.close);

        status_message.append()
            .addClass("coge-progress-status");

        var dialog = $(selection).append(data)
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

        my.close = function() {
            console.log("TEST");
            dialog.dialog("close");
            return my;
        };

        my.status = function(response) {
            var message = $("<span></span>");
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
                status_message.append(message.html(response.message));
            }

            link.attr("href", response.link.href)
                .css("margin", "2%")
                .html(response.link.message);

            status_message
                .append(link)
                .append(button)
                .slideDown();

            return my;
        };

        my.show = function() {
            dialog.dialog("open");
            return my;
        };

        my.title = function(_) {
            if (!arguments.length) return title;

            title = _;
            return my;
        };

        my.data = function(_) {
            if (!arguments.length) return data;

            data.html(_);
            return my;
        };

        return my;
    }

    return progress;
})();
