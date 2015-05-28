window.wikifeed = function(url, element, size) {
	url = 'https://genomevolution.org/wiki/index.php/CoGepedia:Current_events';
    var num_posts = size || 5;
    var div = $("<div></div>");

    function add(header) {
    	var span = header.firstChild;
        var date = $("<span></span>").html($(header).next().html()).addClass("date");
        var link = $("<a></a>");
        link.attr("href", url + "#" + span.id).html(span.innerHTML);
        var title = $("<h4></h4>");
        div.append(title.html(link).append(date));
    }

    function parse(html) {
        var feed = $(html).find("#mw-content-text");
        var headers = feed.find("h2").toArray();
        if (!headers.length)
        	return console.warn("Could not fetch headers");
        
        for (var i=0; i<num_posts; i++)
        	add(headers[i]);
        var footer = $("<a>...more...<a>").attr("href", url)
	        .attr("target", "_blank")
	        .css("text-align", "center")
	        .css("display", "block")
	        .css("margin", "3px");
	    div.append(footer);
	    element.append(div);
    }

    $.ajax({
        url: url,
        dataType: "html",
        success: parse
    });

    return this;
};
