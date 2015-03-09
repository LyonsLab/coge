(function(coge) {
    var color = coge.color = {};

    color.schemes = {
        "Black-Red": [
            "#ff0000",
            "#000000"
        ],
        "Rainbow 1": [
            "#ffff00",
            "#c8c800",
            "#00c800",
            "#006464",
            "#00c8c8",
            "#0000c8",
            "#640064",
            "#c800c8",
            "#c80000",
            "#640000",
            "#c86400",
            "#ff7e00"
        ],
        "Rainbow 2": [
            "#ffff00",
            "#ff0000",
            "#00ff00",
            "#00ffff",
            "#dc00dc",
            "#0000ff"
        ],
        "2.1xRainbow": [
            "#ff0000",
            "#ffff00",
            "#00ff00",
            "#00ffff",
            "#dc00dc",
            "#0000ff",
            "#ff0000",
            "#ffff00",
            "#00ff00",
            "#00ffff",
            "#dc00dc",
            "#0000ff"
        ],
        "2.2xRainbow": [
            "#ff0000",
            "#ffff00",
            "#00ff00",
            "#00ffff",
            "#dc00dc",
            "#0000ff",
            "#ff0000",
            "#ffff00",
            "#00ff00",
            "#00ffff",
            "#dc00dc",
            "#0000ff",
            "#ff0000",
            "#ffff00",
            "#00ff00",
            "#00ffff",
            "#dc00dc",
            "#0000ff"
        ],
        "3x1Rainbow": [
            "#ff0000",
            "#ff9b00",
            "#00ff00",
            "#00ffff",
            "#dc00dc",
            "#0000ff",
            "#ff0000",
            "#ff9b00",
            "#00ff00",
            "#00ffff",
            "#dc00dc",
            "#0000ff"
        ],
        "3x2Rainbow": [
            "#ff0000",
            "#ff9b00",
            "#00ff00",
            "#00ffff",
            "#dc00dc",
            "#0000ff",
            "#ff0000",
            "#ff9b00",
            "#00ff00",
            "#00ffff",
            "#dc00dc",
            "#0000ff",
            "#ff0000",
            "#ff9b00",
            "#00ff00",
            "#00ffff",
            "#dc00dc",
            "#0000ff"
        ],
        "RYB": [
            "#000063",
            "#dcdc14",
            "#ff0000"
        ],
        "RYBG": [
            "#00c800",
            "#0000c8",
            "#dcdc14",
            "#ff0000"
        ],
        "3xRed-Blue": [
            "#00ffff",
            "#ff0000",
            "#00ffff",
            "#ff0000",
            "#00ffff",
            "#ff0000"
        ],
        "3xBlue-Orange": [
            "#0000ff",
            "#ff6333",
            "#0000ff",
            "#ff6333",
            "#0000ff",
            "#ff6333"
        ]
    };

    color.dropdown = function(element, colorSchemes) {
        var my = {},
            listeners = [],
            el = $(element),
            schemes = colorSchemes || coge.color.schemes,
            names,
            items;

        my.selected = function(callback) {
            listeners.push(callback);
        };

        el.on("change", function() {
            var colors = schemes[$(this).val()];
            my.handleSelected(colors);
        });

        my.handleSelected = function(colors) {
            listeners.forEach(function(callback) {
                callback(colors);
            });
        };

        my.select = function(index) {
            var child = el.children()[index];
            if (child) {
                $(el.child).attr("selected", "selected");
                el.change();
            }
        };

        names = Object.keys(coge.color.schemes);
        items = names.map(function(name) {
            return $("<option>")
                .attr("value", name)
                .html(name);
        });

        // Add color schemes
        el.append(items);

        return my;
    };
}(window.coge || (window.coge = {})));
