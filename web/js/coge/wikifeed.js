window.wikifeed = function(url, element, size) {
    var posts = [],
        numberOfPosts = size || 5,
        slice = Array.prototype.slice;

    function getIndex(children) {
        return function(value, index, list) {
            return children.index(value);
        };
    }

    function getRange(value, index, list) {
        return [value, list[index + 1]];
    };

    function getData(children) {
        return function(value, index, list) {
            //offset by 1 to get anchor tag
            var offset = [value[0] - 1, value[1] - 1];
            return slice.apply(children, offset);
        };
    }

    function parseItem(value, index, list) {
        var elem = $(value);

        if(elem.is("div")) {
            var inner = $(value).find("a");
            var images = inner.find("img,embed,object");
            images.removeAttr("width");
            images.removeAttr("height");

            return inner;
        }

        return value;
    };

    function formatPost(value, index, list) {
        var id = $(value[0]).attr("id");

        var link = $("<a></a>");
        link.attr("href", url + "#" + id)
            .attr("target", "_blank")
            .addClass("ui-icon ui-icon-extlink")
            .css("float", "right");

        var content = $("<div></div>");

        var items = value.splice(3, value.length);
        items.unshift(link);

        var title = $(value[1]);
        var text = title.text().replace(/\[edit\]/gmi, "");
        var date = $("<span></span>").html($(value[2]).text())
            .addClass("date");

        var parsed = items.map(parseItem);

        return {
            title: title.text(text).append(date),
            content: content.html(parsed)
        };
    }

    function parsePosts(html) {
        var feed = $(html).find("#bodyContent");
        var children = feed.children();

        var headers = feed.find("h2").toArray();

        if (!headers.length) return console.warn("Could not fetch headers");

        var indexes = headers.map(getIndex(children));
        indexes.push(children.length);
        var ranges = indexes.map(getRange);
        var items = ranges.map(getData(children.toArray()));

        posts = items.map(formatPost).slice(0, -1);
        render();
    }

    function render() {
        var feed = $("<div></div>");

        for(var i = 0; i < numberOfPosts; i++) {
            feed.append(posts[i].title);
            feed.append(posts[i].content);
        }

        feed.accordion({
            autoHeight: false,
            icons: false
        });

        var footer = $("<a>...more...<a>").attr("href", url)
            .attr("target", "_blank")
            .css("text-align", "center")
            .css("display", "block")
            .css("margin", "3px");

        element.append(feed);
        $(".thumbcaption").remove();
        feed.append(footer);
    }

    $.ajax({
        url: url,
        dataType: "html",
        success: parsePosts
    });

    return this;
};
