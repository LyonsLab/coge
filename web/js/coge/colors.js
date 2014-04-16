var coge = window.coge = (function(ns) {
    ns.color = {
        schemes: {
            "Rainbow 1": [
                "rgb(255, 255, 0)",
                "rgb(200, 200, 0)",
                "rgb(0,   200, 0)",
                "rgb(0,   100, 100)",
                "rgb(0,   200, 200)",
                "rgb(0,   0,   200)",
                "rgb(100, 0,   100)",
                "rgb(200, 0,   200)",
                "rgb(200, 0,   0)",
                "rgb(100, 0,   0)",
                "rgb(200, 100, 0)",
                "rgb(255, 126, 0)"
            ],
            "Rainbow 2": [
                "rgb(255, 255, 0)",
                "rgb(255, 0,   0)",
                "rgb(0,   255, 0)",
                "rgb(0,   255, 255)",
                "rgb(220, 0,   220)",
                "rgb(0,   0,   255)"
            ],
            "2.1xRainbow": [
                "rgb(255, 0,   0)",
                "rgb(255, 255, 0)",
                "rgb(0,   255, 0)",
                "rgb(0,   255, 255)",
                "rgb(220, 0,   220)",
                "rgb(0,   0,   255)",
                "rgb(255, 0,   0)",
                "rgb(255, 255, 0)",
                "rgb(0,   255, 0)",
                "rgb(0,   255, 255)",
                "rgb(220, 0,   220)",
                "rgb(0,   0,   255)"
            ],
            "2.2xRainbow": [
                "rgb(255, 0,   0)",
                "rgb(255, 255, 0)",
                "rgb(0,   255, 0)",
                "rgb(0,   255, 255)",
                "rgb(220, 0,   220)",
                "rgb(0,   0,   255)",
                "rgb(255, 0,   0)",
                "rgb(255, 255, 0)",
                "rgb(0,   255, 0)",
                "rgb(0,   255, 255)",
                "rgb(220, 0,   220)",
                "rgb(0,   0,   255)",
                "rgb(255, 0,   0)",
                "rgb(255, 255, 0)",
                "rgb(0,   255, 0)",
                "rgb(0,   255, 255)",
                "rgb(220, 0,   220)",
                "rgb(0,   0,   255)"
            ],
            "3x1Rainbow": [
                "rgb(255, 0,   0)",
                "rgb(255, 155, 0)",
                "rgb(0,   255, 0)",
                "rgb(0,   255, 255)",
                "rgb(220, 0,   220)",
                "rgb(0,   0,   255)",
                "rgb(255, 0,   0)",
                "rgb(255, 155, 0)",
                "rgb(0,   255, 0)",
                "rgb(0,   255, 255)",
                "rgb(220, 0,   220)",
                "rgb(0,   0,   255)"
            ],
            "3x2Rainbow": [
                "rgb(255, 0,   0)",
                "rgb(255, 155, 0)",
                "rgb(0,   255, 0)",
                "rgb(0,   255, 255)",
                "rgb(220, 0,   220)",
                "rgb(0,   0,   255)",
                "rgb(255, 0,   0)",
                "rgb(255, 155, 0)",
                "rgb(0,   255, 0)",
                "rgb(0,   255, 255)",
                "rgb(220, 0,   220)",
                "rgb(0,   0,   255)",
                "rgb(255, 0,   0)",
                "rgb(255, 155, 0)",
                "rgb(0,   255, 0)",
                "rgb(0,   255, 255)",
                "rgb(220, 0,   220)",
                "rgb(0,   0,   255)"
            ],
            "RYB": [
                "rgb(0, 0, 150)",
                "rgb(220, 220, 20)",
                "rgb(255, 0, 0)"
            ],
            "RYBG": [
                "rgb(0, 200, 0)",
                "rgb(0, 0, 200)",
                "rgb(220, 220, 20)",
                "rgb(255, 0, 0)"
            ],
            "Black-Red": [
                "rgb(255, 0, 0)",
                "rgb(0, 0, 0)"
            ],
            "3xRed-Blue": [
                "rgb(0, 255, 255)",
                "rgb(255, 0, 0)",
                "rgb(0, 255, 255)",
                "rgb(255, 0, 0)",
                "rgb(0, 255, 255)",
                "rgb(255, 0, 0)"
            ],
            "3xBlue-Orange": [
                "rgb(0, 0, 255)",
                "rgb(255, 99, 33)",
                "rgb(0, 0, 255)",
                "rgb(255, 99, 33)",
                "rgb(0, 0, 255)",
                "rgb(255, 99, 33)"
            ]
        },
        dropdown: function(element, colorSchemes) {
            var my = {},
                listeners = [],
                el = $(element),
                schemes = colorSchemes || coge.color.schemes,
                items;

            my.selected = function(callback) {
                listeners.push(callback);
            };

            el.on("change", function(e) {
                var colors = schemes[$(this).val()];
                my.handleSelected(colors);
            });

            my.handleSelected = function(colors) {
                listeners.forEach(function(callback) {
                    callback(colors);
                });
            }

            names = Object.keys(coge.color.schemes);
            items = names.map(function(name) {
                return $("<option>")
                    .attr("value", name)
                    .html(name);
            });

            // Add color schemes
            el.append(items);

            return my;
        }
    };

    return ns;
})(coge || {});
