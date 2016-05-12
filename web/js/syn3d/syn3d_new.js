/**
 * Created by senorrift on 4/29/16.
 */

/*-----------------------------------------------------------------------------------------------------------------
 ~~~~ GLOBAL DECLARATIONS ~~~~
 -----------------------------------------------------------------------------------------------------------------*/

// SynMap Global Variables [note: when adding, make sure to reset in renderSynMap().initialize()]
var camView = { "x": 0, "y": 0, "z": 120 };
var histData = {"kn": [], "ks": [], "knks": []};
var pointData = [];
var colors = {"kn": [], "ks": [], "knks": []};

// Updating Variables
var camUpdate = false;
var hQueue = ["hKnKs", "hKs", "hKn"];
var hCurrent = ["hKs", "ks"];
var gridState = true;
var colorScheme = "Jet";
var dataSubset;

// Color schemes variables
var schemes = {
    "Jet": [[0,'rgb(0,0,131)'], [0.125,'rgb(0,60,170)'], [0.375,'rgb(5,255,255)'],
        [0.625,'rgb(255,255,0)'], [0.875,'rgb(250,0,0)'], [1,'rgb(128,0,0)']],

    "Bluered": [[0,'rgb(0,0,255)'],[1,'rgb(255,0,0)']],

    "Portland": [[0,'rgb(12,51,131)'],[0.25,'rgb(10,136,186)'], [0.5,'rgb(242,211,56)'],
        [0.75,'rgb(242,143,56)'], [1,'rgb(217,30,30)']],

    "Viridis": [[0,'#440154'],[0.06274509803921569,'#48186a'], [0.12549019607843137,'#472d7b'],
        [0.18823529411764706,'#424086'], [0.25098039215686274,'#3b528b'],[0.3137254901960784,'#33638d'],
        [0.3764705882352941,'#2c728e'], [0.4392156862745098,'#26828e'], [0.5019607843137255,'#21918c'],
        [0.5647058823529412,'#1fa088'], [0.6274509803921569,'#28ae80'],[0.6901960784313725,'#3fbc73'],
        [0.7529411764705882,'#5ec962'], [0.8156862745098039,'#84d44b'], [0.8784313725490196,'#addc30'],
        [0.9411764705882353,'#d8e219'], [1,'#fde725']]
}; // Source: https://github.com/plotly/plotly.js/blob/master/src/components/colorscale/scales.js

/* CORE FUNCTION: Create coloring function */
function getPtColorFunc(minVal, maxVal, colorScheme) {
    var N = colorScheme.length,
        domain = new Array(N),
        range = new Array(N),
        si;

    for (var i = 0; i < N; i++) {
        si = colorScheme[i];
        domain[i] = minVal + si[0] * (maxVal - minVal);
        range[i] = si[1]
    }

    var calcColor = d3.scale.linear()
            .domain(domain)
            .interpolate(d3.interpolateRgb)
            .range(range);

    return calcColor
}

/* CORE FUNCTION: Load Data */
function loadData(graph_json) {
    var graphObj;

    function getData() {
        return $.getJSON(graph_json, function(json) { graphObj = json; });
    }

    return $.when( getData() ).then(function() {
        console.log(graphObj);
        return graphObj;
    });
}

/* CORE FUNCTION: Rotate Histograms */
function rotateHistogram(direction) {
    var n, c, type;
    if (direction == "R") {
        n = hQueue.pop();
        c = hQueue.pop();
        type = n.slice(1).toLowerCase();
        d3.select("#chartSvg").remove();
        document.getElementById(c).style.display = "none";
        document.getElementById(n).style.display = "";
        renderHistogram(n, histData[type]);
        hCurrent = [n, type];
        hQueue.unshift(n);
        hQueue.unshift(c);
    } else if (direction=="L") {
        n = hQueue.shift();
        c = hQueue.shift();
        type = n.slice(1).toLowerCase();
        d3.select("#chartSvg").remove();
        document.getElementById(c).style.display = "none";
        document.getElementById(n).style.display = "";
        renderHistogram(n, histData[type]);
        hCurrent = [n, type];
        hQueue.push(n);
        hQueue.push(c);
    } else {
        console.log("ERROR: Unrecognized Direction")
    }
}

/* CORE FUNCTION: Reset Camera */
function resetCamera() {
    camUpdate = true;
    camView = { "x": 0, "y": 0, "z": 80 };
}

/* CORE FUNCTION: Toggle Grid Display */
function toggleGrid() {
    gridState = !gridState;
}

/* CORE FUNCTION: Empty Renderings */
function emptyRenderings() {
    // Remove old SynMap
    var oldCanvas = document.getElementById("canvas");
    oldCanvas.parentNode.removeChild(oldCanvas);
    // Add new canvas
    var newCanvas = document.createElement('canvas');
    newCanvas.id = 'canvas';
    document.getElementById('rendering').appendChild(newCanvas);
    // Remove old histogram.
    d3.select("#chartSvg").remove();
}

/* CORE FUNCTION: Render Histograms */
function renderHistogram(element_id, values) {
    /* RESOURCES
        Simple Histogram: https://bl.ocks.org/mbostock/3048450
        Basic Mouseover: http://bl.ocks.org/phil-pedruco/9032348
        Tooltip Mouseover: http://bl.ocks.org/Caged/6476579
        Brushing (range selector): http://bl.ocks.org/cjfont/4420257
     */

    // Convert type of values.
    for (var i = 0; i < values.length; i++) {
        values[i] = +values[i];
    }

    // Format count function.
    var formatCount = d3.format(",.3f");
    var formatWide = d3.format(",.10f");
    // Bin Count.
    var binCount = 100.0;

    // Set margins & canvas size.
    var margin = {top: 10, right: 30, bottom: 30, left: 30},
        width = document.getElementById(element_id).clientWidth - margin.left - margin.right,
        height = (document.getElementById("rendering").clientHeight * .65) - margin.bottom - margin.top;
    
    // Get minimum & maximum values.
    var minVal = Math.min.apply(null, values);
    var maxVal = Math.max.apply(null, values);

    // Generate function to calculate point colors by colorscheme.
    var calcColor = getPtColorFunc(minVal, maxVal, schemes[colorScheme]);

    // Set X scale.
    var x = d3.scale.linear()
        .domain([minVal, maxVal])
        .range([0, width]);
    
    // Generate histogram bins.
    var data = d3.layout.histogram()
        .bins(binCount)
        (values);

    // Set Y scale.
    var y = d3.scale.linear()
        .domain([0, d3.max(data, function(d) { return d.y; })])
        .range([height, 0]);

    // Build X axis object.
    var xAxis = d3.svg.axis()
        .scale(x)
        .orient("bottom");

    // Build tip object.
    var tip = d3.tip()
        .attr('class', 'd3-tip')
        .offset([-10, 0])
        .html(function(d) {
            //console.log(d.x); // Returns bin minimum
            var range = formatCount(d3.min(d)) + " < log10(x) < " + formatCount(d3.max(d));
            var count = "Bin Count: <span style='color:red'>" + d.y + "</span>";
            return range + "<br>" + count;
        });

    // Draw SVG in specified element.
    var svg = d3.select('#' + element_id).append("svg")
        .attr("id", "chartSvg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
        .append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    // Add tip to SVG.
    svg.call(tip);

    // Draw histogram bars.
    // TODO: FIX overlapping bars (maybe has to do with SVG element size). Maybe use math.floor
    var bar = svg.selectAll(".bar")
        .data(data)
        .enter().append("g")
        .attr("class", "bar")
        .attr("transform", function(d) { return "translate(" + formatWide(x(formatWide(d.x))) + "," + y(d.y) + ")"; });

    bar.append("rect")
        .attr("x", 1)
        .attr("width", formatWide(width/(binCount)) )
        .attr("height", function(d) { return height - y(d.y); })
        .attr("fill", function(d) { return calcColor(d.x) })
        .on("mouseover", function(d) {
            tip.show(d);
            d3.select(this).attr("fill", "red");
        })
        .on("mouseout", function(d) {
            tip.hide(d);
            d3.select(this).attr("fill", function(d) { return calcColor(d.x) })
        });

    // Set up & animate selection brush.
    var brush, brushg;

    function brushstart() {}

    function brushmove() {}

    function brushend() {
        if (!brush.empty()) {
            dataSubset = brush.extent();
        } else {
            dataSubset = [minVal, maxVal]
        }
        console.log(dataSubset);
    }

    function resizePath(d) {
        var e = +(d == "e"),
            x = e ? 1 : -1,
            h = height * 2/3,
            y = ((height - h) / 2);
        return "M" + (.5 * x) + "," + y
            + "A6,6 0 0 " + e + " " + (6.5 * x) + "," + (y + 6)
            + "V" + (y + h - 6)
            + "A6,6 0 0 " + e + " " + (.5 * x) + "," + (y + h)
            + "Z"
            + "M" + (2.5 * x) + "," + (y + 8)
            + "V" + (y + h - 8)
            + "M" + (4.5 * x) + "," + (y + 8)
            + "V" + (y + h - 8);
    }

    brush = d3.svg.brush()
        .x(x)
        .on("brushstart", brushstart)
        .on("brush", brushmove)
        .on("brushend", brushend);

    brushg = svg.append("g")
        .attr("class", "brush")
        .call(brush);

    brushg.selectAll("rect")
        .attr("height", height);

    brushg.selectAll(".resize").append("path").attr("d", resizePath);

    if (brush) { brushmove(); }

    // Draw X axis.
    svg.append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + height + ")")
        .call(xAxis);
}

/* CORE FUNCTION: Render SynMap */
function renderSynMap(graph_object, element_id, color_by) {
    /*---------------------------------------------------------------------------------------------------------
     ~~~~VARIABLE DECLARATIONS~~~~
     --------------------------------------------------------------------------------------------------------*/

    // Rendering Variables
    var renderer, scene, camera, controls;
    var container = document.getElementById(element_id);
    var width, height, containerWidth, containerHeight;

    // Visualization & Data variables.
    var points;
    var xAxis, yAxis, zAxis;
    var grid;

    // Scaling variables.
    var axisWidth = 0.5;
    var pointSize = 8;

    // Raycasting variables.
    var raycaster, intersects, mouse, INTERSECTED;
    var ptHistory = ["NULL", "NULL"];
    var black = new THREE.Color("#000");


    /*---------------------------------------------------------------------------------------------------------
     ~~~~INITIALIZE & ANIMATE SYNMAP~~~~
     --------------------------------------------------------------------------------------------------------*/

    initialize();
    animate();

    /*---------------------------------------------------------------------------------------------------------
     ~~~~MINOR FUNCTIONS~~~~
     --------------------------------------------------------------------------------------------------------*/

    /* FUNCTION: Build Axis from Chromosome List */
    function buildAxisFromChr(label, chrNames, chrStarts, maxLen, focusAxis) {
        var axis = new THREE.Object3D();
        var colorToggle = true;

        var altAxes = ['x', 'y', 'z'];
        altAxes.splice(altAxes.indexOf(focusAxis), 1);

        var colorOne, colorTwo;
        switch (focusAxis) {
            case "x":
                colorOne = 'darkred';
                colorTwo = 'red';
                break;
            case "y":
                colorOne = 'darkblue';
                colorTwo = 'blue';
                break;
            case "z":
                colorOne = 'green';
                colorTwo = 'lightgreen';
                break;
            default:
                colorOne = '';
                colorTwo = '';
        }
        for (var i = 0; i<chrNames.length; i++) {
            // Set chromosome length.
            var chLen;
            if (i == chrNames.length - 1) {
                chLen = maxLen - chrStarts[i]
            } else {
                chLen = chrStarts[i+1] - chrStarts[i]
            }
            // Set chromosome geometry.
            var chGeo;
            if (focusAxis == 'x') {
                chGeo = new THREE.BoxGeometry( chLen, axisWidth, axisWidth );
            } else if (focusAxis == 'y') {
                chGeo = new THREE.BoxGeometry( axisWidth, chLen, axisWidth );
            } else if (focusAxis == 'z') {
                chGeo = new THREE.BoxGeometry( axisWidth, axisWidth, chLen );
            } else {
                console.log("Invalid Focus Axis, please select x y or z");
            }
            // Set chromosome material.
            var chMat;
            if (colorToggle) {
                chMat = new THREE.MeshBasicMaterial({color: colorOne});
                colorToggle = false;
            } else {
                chMat = new THREE.MeshBasicMaterial({color: colorTwo});
                colorToggle = true;
            }
            // Make chromosome mesh.
            var ch = new THREE.Mesh( chGeo, chMat );
            // Position chromosome.
            ch.position[focusAxis] = chrStarts[i] + (chLen / 2);
            ch.position[altAxes[0]] = -(axisWidth / 2);
            ch.position[altAxes[1]] = -(axisWidth / 2);
            // Add chromosome to axis.
            axis.add( ch );
        }

        // Draw axis labels.
        var labelSprite;
        if (focusAxis =='x') {
            labelSprite = makeTextSprite(label, { fontsize: 12, borderColor: {r:139, g:0, b:0, a:1.0}, backgroundColor: {r:255, g:100, b:100, a:0.65} });
            labelSprite.geometry.center();
            labelSprite.position.set(maxLen*1.5,-(axisWidth / 2),-(axisWidth / 2));
            //console.log(labelSprite);
        } else if (focusAxis == 'y') {
            labelSprite = makeTextSprite(label, { fontsize: 12, borderColor: {r:0, g:0, b:139, a:1.0}, backgroundColor: {r:0, g:100, b:255, a:0.65} });
            labelSprite.geometry.center();
            labelSprite.position.set(0,maxLen,0); //TODO (figure out which one is correct)
        } else if (focusAxis == 'z') {
            labelSprite = makeTextSprite(label, { fontsize: 12, borderColor: {r:0, g:128, b:0, a:1.0}, backgroundColor: {r:144, g:238, b:144, a:0.65} });
            labelSprite.geometry.center();
            labelSprite.position.set(0,0,maxLen);
        } else {console.log("Invalid focus axis, please select x y or z.")}

        axis.add(labelSprite);

        return axis;
    }

    /* FUNCTION: Generate GEvo Link */
    function genGevoLink(xDbId, yDbId, zDbId) {
        // https://genomevolution.org/wiki/index.php/Linking_to_GEvo
        //var linkbase = "http://genomevolution.org/coge/GEvo.pl?";
        var linkbase = "http://geco.iplantcollaborative.org/asherkhb/coge/GEvo.pl?";
        var seqX = 'fid1=' + xDbId + ';';
        var seqY = 'fid2=' + yDbId + ';';
        var seqZ = 'fid3=' + zDbId + ';';
        var range = 'apply_all=50000;';
        var numSeq = 'num_seqs=3';
        return linkbase + seqX + seqY + seqZ + range + numSeq;
    }

    /* FUNCTION: Histogram Change */
    function updateColors(select) {
        var geometry = points.geometry;
        geometry.colors = colors[select];
        geometry.colorsNeedUpdate = true;
    }

    /* FUNCTION: Make Text Sprites */
    function makeTextSprite( message, parameters ) {
        if ( parameters === undefined ) parameters = {};

        var fontface = parameters.hasOwnProperty("fontface") ?
            parameters["fontface"] : "Arial";

        var fontsize = parameters.hasOwnProperty("fontsize") ?
            parameters["fontsize"] : 18;

        var borderThickness = parameters.hasOwnProperty("borderThickness") ?
            parameters["borderThickness"] : 2;

        var borderColor = parameters.hasOwnProperty("borderColor") ?
            parameters["borderColor"] : { r:0, g:0, b:0, a:1.0 };

        var backgroundColor = parameters.hasOwnProperty("backgroundColor") ?
            parameters["backgroundColor"] : { r:255, g:255, b:255, a:1.0 };

        //var spriteAlignment = THREE.SpriteAlignment.topLeft;

        var canvas = document.createElement('canvas');
        var context = canvas.getContext('2d');
        context.font = "Bold " + fontsize + "px " + fontface;

        // get size data (height depends only on font size)
        var metrics = context.measureText( message );
        var textWidth = metrics.width;

        // background color
        context.fillStyle   = "rgba(" + backgroundColor.r + "," + backgroundColor.g + ","
                                      + backgroundColor.b + "," + backgroundColor.a + ")";
        // border color
        context.strokeStyle = "rgba(" + borderColor.r + "," + borderColor.g + ","
                                      + borderColor.b + "," + borderColor.a + ")";

        context.lineWidth = borderThickness;
        roundRect(context, borderThickness/2, borderThickness/2, textWidth + borderThickness, fontsize * 1.4 + borderThickness, 6);
        // 1.4 is extra height factor for text below baseline: g,j,p,q.

        // text color
        context.fillStyle = "rgba(0, 0, 0, 1.0)";

        context.fillText( message, borderThickness, fontsize + borderThickness);

        // canvas contents will be used for a texture
        var texture = new THREE.Texture(canvas);
        texture.needsUpdate = true;

        var spriteMaterial = new THREE.SpriteMaterial( { map: texture } );
        var sprite = new THREE.Sprite( spriteMaterial );

        sprite.scale.set(100,50,1.0);
        return sprite;
    } // Source: http://stemkoski.github.io/Three.js/Sprite-Text-Labels.html

    /* FUNCTION: Draw Rectangles (for text sprites) */
    function roundRect(ctx, x, y, w, h, r) {
        ctx.beginPath();
        ctx.moveTo(x+r, y);
        ctx.lineTo(x+w-r, y);
        ctx.quadraticCurveTo(x+w, y, x+w, y+r);
        ctx.lineTo(x+w, y+h-r);
        ctx.quadraticCurveTo(x+w, y+h, x+w-r, y+h);
        ctx.lineTo(x+r, y+h);
        ctx.quadraticCurveTo(x, y+h, x, y+h-r);
        ctx.lineTo(x, y+r);
        ctx.quadraticCurveTo(x, y, x+r, y);
        ctx.closePath();
        ctx.fill();
        ctx.stroke();
    } // Source: http://stemkoski.github.io/Three.js/Sprite-Text-Labels.html

    /*---------------------------------------------------------------------------------------------------------
     ~~~~MAJOR FUNCTIONS: CREATE 3D OBJECTS~~~~
     --------------------------------------------------------------------------------------------------------*/

    /* FUNCTION: Draw Points */
    function drawPoints(points) {
        // Load point image.
        var loader = new THREE.TextureLoader();
        var sprite = loader.load("picts/ball.png");
        // Define plot geometry & material & color array.
        var plotGeo = new THREE.Geometry();
        //plotGeo.addAttribute(size, new THREE.BufferAttribute( size, 1 ));
        var plotMat = new THREE.PointsMaterial( {size: pointSize, sizeAttenuation: false, map: sprite,
            vertexColors: THREE.VertexColors, alphaTest: 0.5, transparent: true });
        //var colors = [];

        //Set up point counter.
        var ptCount = points.length;
        var pointMutData = {"kn": [], "ks": [], "knks": []};

        // Build points.
        for (var i=0; i<ptCount; i++ ) {
            // Assign position
            var vertex = new THREE.Vector3();
            vertex.x = points[i][0];
            vertex.y = points[i][1];
            vertex.z = points[i][2];
            plotGeo.vertices.push(vertex);

            // Add reference to pointData for retrieving feature info.
            pointData.push(i);

            // Get Substitution Data
            // NOTE: CURRENTLY IMPLEMENTED FOR XY!!!!!!!!!!!!!!!!!!!!!!!!!!
            var substitution_indexes = {"xy": [6, 7], "xz": [8, 9], "yz": [10, 11]};
            var kn, ks;
            if (color_by != "xyz") {
                kn = points[i][substitution_indexes[color_by][0]];
                ks = points[i][substitution_indexes[color_by][1]];
            } else if (color_by == "xyz") {
                kn = ( points[i][6] + points[i][8] + points[i][10] ) / 3;
                ks = ( points[i][7] + points[i][9] + points[i][11] ) / 3;
            } else {
                console.log("Error: Unrecognized comparison")
            }

            // var kn = points[i][6];
            // var ks = points[i][7];
            var knks = Math.log10(Math.pow(10, kn) / Math.pow(10, ks));

            if (!isNaN(kn)) {
                histData.kn.push(kn);
                pointMutData.kn.push(kn);
            } else {
                pointMutData.kn.push("NULL");
            }
            if (!isNaN(ks)) {
                histData.ks.push(ks);
                pointMutData.ks.push(ks);
            } else {
                pointMutData.ks.push("NULL");
            }
            if (!isNaN(knks)) {
                histData.knks.push(knks);
                pointMutData.knks.push(knks);
            } else {
                pointMutData.knks.push("NULL");
            }
            
        }

        // Build colors lists.
        var cList = ["ks", "kn", "knks"];
        for (var s=0; s<cList.length; s++) {
            var c = cList[s];
            var minVal = Math.min.apply(null, histData[c]);
            var maxVal = Math.max.apply(null, histData[c]);
            var calcColor = getPtColorFunc(minVal, maxVal, schemes[colorScheme]);

            for (var j=0; j<ptCount; j++ ) {
                if (pointMutData[c][j] != 'NULL') {
                    colors[c][j] = new THREE.Color(calcColor(pointMutData[c][j]))
                } else {
                    colors[c][j] = new THREE.Color("#fff");
                }
            }
        }

        // Apply colors to geometry.
        plotGeo.colors = colors[hCurrent[1]];

        // Create points object.
        var pts = new THREE.Points(plotGeo, plotMat);

        // Shift points object to rotation offset.
        pts.translateX(-graph_object.x[2] / 2);
        pts.translateY(-graph_object.y[2] / 2);
        pts.translateZ(-graph_object.z[2] / 2);

        // Report point count.
        console.log("POINTS RENDERED: " + ptCount);
        return pts;
    }

    /* FUNCTION: Draw Chromosomes */
    function drawChromosomes() {
        /* 3d Object: Axes */
        xAxis = buildAxisFromChr(graph_object.x[1], graph_object.axes[0], graph_object.axes[1], graph_object.x[2], 'x');
        yAxis = buildAxisFromChr(graph_object.y[1], graph_object.axes[2], graph_object.axes[3], graph_object.y[2], 'y');
        zAxis = buildAxisFromChr(graph_object.z[1], graph_object.axes[4], graph_object.axes[5], graph_object.z[2], 'z');

        // Shift axes to place graph center in middle.
        xAxis.translateX(-graph_object.x[2] / 2);
        xAxis.translateY(-graph_object.y[2] / 2);
        xAxis.translateZ(-graph_object.z[2] / 2);
        yAxis.translateX(-graph_object.x[2] / 2);
        yAxis.translateY(-graph_object.y[2] / 2);
        yAxis.translateZ(-graph_object.z[2] / 2);
        zAxis.translateX(-graph_object.x[2] / 2);
        zAxis.translateY(-graph_object.y[2] / 2);
        zAxis.translateZ(-graph_object.z[2] / 2);

        return [xAxis, yAxis, zAxis];
    }

    /* FUNCTION: Draw Grid */
    function drawGrid(graph_object) {
        var g = new THREE.Object3D();
        var gridWidth = 0.25;
        var gridMat = new THREE.MeshBasicMaterial({color: 'grey'});

        var xA = graph_object.axes[1];
        var xL = graph_object.x[2];
        xA.push(xL);
        var yA = graph_object.axes[3];
        var yL = graph_object.y[2];
        yA.push(yL);
        var zA = graph_object.axes[5];
        var zL = graph_object.z[2];
        zA.push(zL);

        for (var x=1; x<xA.length; x++) {
            var xyGeo = new THREE.BoxGeometry(gridWidth, yL, gridWidth);
            var xzGeo = new THREE.BoxGeometry(gridWidth, gridWidth, zL);
            var xyGrid = new THREE.Mesh( xyGeo, gridMat );
            xyGrid.position.set(xA[x], yL/2, 0);
            g.add(xyGrid);
            var xzGrid = new THREE.Mesh( xzGeo, gridMat );
            xzGrid.position.set(xA[x], 0, zL/2);
            g.add(xzGrid);
        }

        for (var y=1; y<yA.length; y++) {
            var yxGeo = new THREE.BoxGeometry(xL, gridWidth, gridWidth);
            var yzGeo = new THREE.BoxGeometry(gridWidth, gridWidth, zL);
            var yxGrid = new THREE.Mesh( yxGeo, gridMat );
            yxGrid.position.set(xL/2, yA[y], 0);
            g.add(yxGrid);
            var yzGrid = new THREE.Mesh( yzGeo, gridMat );
            yzGrid.position.set(0, yA[y], zL/2);
            g.add(yzGrid);
        }

        for (var z=1; z<zA.length; z++) {
            var zxGeo = new THREE.BoxGeometry(xL, gridWidth, gridWidth);
            var zyGeo = new THREE.BoxGeometry(gridWidth, yL, gridWidth);

            var zxGrid = new THREE.Mesh( zxGeo, gridMat );
            zxGrid.position.set(xL/2, 0, zA[z]);
            g.add(zxGrid);

            var zyGrid = new THREE.Mesh( zyGeo, gridMat );
            zyGrid.position.set(0, yL/2, zA[z]);
            g.add(zyGrid);
        }

        g.translateX(-xL / 2);
        g.translateY(-yL / 2);
        g.translateZ(-zL / 2);

        return g
    }

    /*---------------------------------------------------------------------------------------------------------
     ~~~~INITIALIZE~~~~
     --------------------------------------------------------------------------------------------------------*/

    function initialize() {
        /* Reset Globals */
        //camView = { "x": 0, "y": 0, "z": 80 };
        histData = {"kn": [], "ks": [], "knks": []};
        pointData = [];
        colors = {"kn": [], "ks": [], "knks": []};

        /* Set sizing variables */
        width = document.getElementById("rendering").clientWidth;
        height = window.innerHeight * 0.9;
        containerWidth = container.clientWidth;
        containerHeight = container.clientHeight;

        /* Create a three.js scene */
        scene = new THREE.Scene();

        /* Setup three.js WebGL renderer */
        renderer = new THREE.WebGLRenderer( { antialias: true, alpha: true, canvas: container} );
        renderer.setSize( width, height );

        /* Create a three.js camera */
        camera = new THREE.PerspectiveCamera( 75, width / height, 1, 1000 );
        camera.position.x = camView.x;
        camera.position.y = camView.y;
        camera.position.z = camView.z;

        /* Apply FlatOrbitControls to camera. */
        controls = new THREE.FlatOrbitControls( camera );

        /* Add raycaster */
        raycaster = new THREE.Raycaster();
        mouse = new THREE.Vector2();

        /* Create chromosome axes */
        var axes = drawChromosomes();
        for (var i = 0; i < axes.length; i++) {
            scene.add(axes[i])
        }

        /* Create grid */
        grid = drawGrid(graph_object);
        scene.add(grid);

        /* Create points geometry */
        points = drawPoints(graph_object.points);
        scene.add(points);

    }

    /*---------------------------------------------------------------------------------------------------------
     ~~~~RENDER~~~~
     --------------------------------------------------------------------------------------------------------*/

    function render() {
        // Render scene.
        renderer.render( scene, camera );
    }

    /*---------------------------------------------------------------------------------------------------------
     ~~~~ANIMATE~~~~
     --------------------------------------------------------------------------------------------------------*/
    var current = hCurrent[1];
    var gridS = gridState;

    function animate() {
        /* Update controls. */
        controls.update();

        /* Camera control */
        if (camUpdate) {
            camera.position.x = camView.x;
            camera.position.y = camView.y;
            camera.position.z = camView.z;
            camUpdate = false;
        } else {
            /* Update camera position record. */
            camView = { "x": camera.position.x, "y": camera.position.y, "z": camera.position.z };
        }


        /* Check for recoloring */
        if (hCurrent[1] != current) {
            updateColors(hCurrent[1]);
            current = hCurrent[1];
        }

        /* Check for grid toggling */
        if (gridState != gridS) {
            if (gridState) {
                scene.add(grid);
            } else {
                scene.remove(grid);
            }
            gridS = gridState;
        }

        /* Render the scene. */
        render();

        /* Create continuous animation. */
        requestAnimationFrame( animate );
    }

    /*---------------------------------------------------------------------------------------------------------
     ~~~~HANDLERS~~~~
     --------------------------------------------------------------------------------------------------------*/

    /* Window Resize */
    function onWindowResize() {
        width = document.getElementById("rendering").clientWidth;
        height = window.innerHeight * 0.9;
        renderer.setSize( width, height);
        camera.aspect = width / height;
        camera.updateProjectionMatrix();

        containerWidth = container.clientWidth;
        containerHeight = container.clientHeight;

        // Re-render Histograms on page resize
        d3.select("#chartSvg").remove();
        renderHistogram(hCurrent[0], histData[hCurrent[1]]);
    }
    window.addEventListener( 'resize', onWindowResize, false );

    /* Mouse Click */
    function onDocumentMouseDown( event ) {
        // Specify object of interest (points).
        var canvas = $("canvas");
        var geometry = points.geometry;

        // Calculate mouse location on canvas.
        var offset_l = canvas.offset().left - $(window).scrollLeft();
        var offset_t = canvas.offset().top - $(window).scrollTop();
        mouse.x = ( (event.clientX - offset_l) / container.clientWidth ) * 2 - 1;
        mouse.y = -( (event.clientY - offset_t) / container.clientHeight ) * 2 + 1;

        // Cast ray & find intersects.
        raycaster.setFromCamera( mouse, camera );
        intersects = raycaster.intersectObject( points );

        // Actions when intersects are found.
        if (intersects.length > 0) {
            // Set INTERSECTED to closest intersected point.
            INTERSECTED = intersects[ 0 ].index;
            // Assign INTERSECTED point data to pt_select.
            var pt_select = pointData[ INTERSECTED ];
            // Assign FIDs
            var xFid = graph_object.points[pt_select][3];
            var yFid = graph_object.points[pt_select][4];
            var zFid = graph_object.points[pt_select][5];
            // Log a GEvo link.
            console.log(genGevoLink(xFid, yFid, zFid));

            /*
            // Show point info panels & GEvo button.
            $(".gene-info-display").css("display", "");
            $("#gevo_div").css("display", "");
            // Select info panels.
            var xdisplay = $("#x-gene-display");
            var ydisplay = $("#y-gene-display");
            var zdisplay = $("#z-gene-display");
            // Empty old content from info panels.
            xdisplay.empty();
            ydisplay.empty();
            zdisplay.empty();
            // Print intersected point data to respective info panels.
            xdisplay.append(displayGenePair(pt_select, final_experiment.x_gid));
            ydisplay.append(displayGenePair(pt_select, final_experiment.y_gid));
            zdisplay.append(displayGenePair(pt_select, final_experiment.z_gid));

            // Build GEvo link.
            var gevoLink = genGevoLink(
                pt_select[final_experiment.x_gid].db_feature_id,
                pt_select[final_experiment.y_gid].db_feature_id,
                pt_select[final_experiment.z_gid].db_feature_id);
            // Assign link to click event.
            document.getElementById("gevo_button").onclick = function() { window.open(gevoLink) };
            */
            // Revert previously selected point color.
            if (ptHistory[0] != "NULL") { geometry.colors[ ptHistory[0] ] = ptHistory[1]; }
            // Record new point index & color.
            ptHistory = [INTERSECTED, geometry.colors[INTERSECTED]];
            // Change intersected point color to black.
            geometry.colors[ INTERSECTED ] = black;
            // Trigger color update.
            geometry.colorsNeedUpdate = true;

        } else if ( INTERSECTED !== null ) {
            // Clear INTERSECTED variable if no intersections on click.
            INTERSECTED = null;
        }
    }
    window.addEventListener('mousedown', onDocumentMouseDown, false);

}

/*-----------------------------------------------------------------------------------------------------------------
 ~~~~ POPULATE PAGE ~~~~
 -----------------------------------------------------------------------------------------------------------------*/

$(document).ready( function() {
    var d;
    // TODO: Start Spinny Wheely(s) & end after loading
    // Load data & launch initial visualizations
    $.when(loadData(final_experiment.links.graph)).done(function(data) {
        d = data;
        // Render initial SynMap & Histogram
        renderSynMap(data, "canvas", "xy");
        renderHistogram(hCurrent[0], histData[hCurrent[1]]);
    });

    /* Monitor mutation ratio coloring option & update visualizations on change. */
    var colorBySelect = $("#color_by");
    colorBySelect.change( function () {
        emptyRenderings();
        // Draw new SynMap & histogram.
        renderSynMap(d, "canvas", colorBySelect.val());
        renderHistogram(hCurrent[0], histData[hCurrent[1]]);
    });

    /* Monitor mutation ratio coloring option & update visualizations on change. */
    var colorSchemeSelect = $("#color_scheme");
    colorSchemeSelect.change( function () {
        emptyRenderings();
        colorScheme = colorSchemeSelect.val();
        // Draw new SynMap & histogram.
        renderSynMap(d, "canvas", colorBySelect.val());
        renderHistogram(hCurrent[0], histData[hCurrent[1]]);
    });

});
