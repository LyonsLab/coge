/**
 * Created by senorrift on 4/29/16.
 */

console.log("Syn3D_New Loaded!");

var camView = { "x": 0, "y": 0, "z": 80 };
var histData = {"kn": [], "ks": [], "knks": []};
var colorScheme = "Jet";
var pointData = [];
var colors = {"kn": [], "ks": [], "knks": []};
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
var hQueue = ["hKnKs", "hKs", "hKn"];
var hCurrent = ["hKs", "ks"];
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

/* CORE FUNCTION: Render Histograms */
function renderHistogram(element_id, values) {
    /* RESOURCES
        Simple Histogram: https://bl.ocks.org/mbostock/3048450
        Basic Mouseover: http://bl.ocks.org/phil-pedruco/9032348
        Tooltip Mouseover: http://bl.ocks.org/Caged/6476579
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
    // TODO: FIX overlapping bars (maybe has to do with SVG element size).
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

    // Append count labels to bars.
    // bar.append("text")
    //     .attr("dy", ".75em")
    //     .attr("y", 6)
    //     .attr("x", width/binCount/2)
    //     .attr("text-anchor", "middle")
    //     .text(function(d) { return formatCount(d.y); });

    // Draw X axis.
    svg.append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + height + ")")
        .call(xAxis);
}

/* CORE FUNCTION: Render SynMap */
function renderSynMap(graph_object, element_id) {
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
    function buildAxisFromChr(chrNames, chrStarts, maxLen, focusAxis) {
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
            var kn = points[i][6];
            var ks = points[i][7];
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

        // Assign colors.
        // var compar = "ks";
        // var minVal = Math.min.apply(null, histData[compar]);
        // var maxVal = Math.max.apply(null, histData[compar]);
        // var calcColor = getPtColorFunc(minVal, maxVal, schemes[colorScheme]);
        //
        // for (var j=0; j<ptCount; j++ ) {
        //     if (pointMutData[compar][j] != 'NULL') {
        //         colors[j] = new THREE.Color(calcColor(pointMutData[compar][j]))
        //     } else {
        //         colors[j] = new THREE.Color("#fff");
        //     }
        // }

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
        plotGeo.colors = colors.ks;

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
        xAxis = buildAxisFromChr(graph_object.axes[0], graph_object.axes[1], graph_object.x[2], 'x');
        yAxis = buildAxisFromChr(graph_object.axes[2], graph_object.axes[3], graph_object.y[2], 'y');
        zAxis = buildAxisFromChr(graph_object.axes[4], graph_object.axes[5], graph_object.z[2], 'z');

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

    /* TODO: FUNCTION: Draw Grid */
    // function drawGrid(xLengths, xMax, yLengths, yMax, zLengths, zMax) {
    //     var g = new THREE.Object3D();
    //     var gridWidth = 0.25;
    //     var gridMat = new THREE.MeshBasicMaterial({color: 'grey'});
    //
    //     for (var i=1; i<xLengths.length; i++) {
    //         var start;
    //         if (i == xLengths.length - 1) {
    //             start =
    //         } else {
    //             start = xLengths[i]
    //         }
    //         var yxGridGeo = new THREE.BoxGeometry( startPositions.x, gridWidth, gridWidth );
    //         var yxGrid = new THREE.Mesh( yxGridGeo, gridMat );
    //         yxGrid.position.x = (startPositions.x / 2);
    //         yxGrid.position.y = e.end;
    //         yxGrid.position.z = 0;
    //         g.add(yxGrid);
    //
    //         var yzGridGeo = new THREE.BoxGeometry (gridWidth, gridWidth, startPositions.z);
    //         var yzGrid = new THREE.Mesh( yzGridGeo, gridMat);
    //         yzGrid.position.x = 0;
    //         yzGrid.position.y = e.end; //(yStartPos / 2);
    //         yzGrid.position.z = (startPositions.z / 2);
    //         g.add(yzGrid)
    //     }
    //
    //     yChr.forEach( function(e) {
    //         var yxGridGeo = new THREE.BoxGeometry( startPositions.x, gridWidth, gridWidth );
    //         var yxGrid = new THREE.Mesh( yxGridGeo, gridMat );
    //         yxGrid.position.x = (startPositions.x / 2);
    //         yxGrid.position.y = e.end;
    //         yxGrid.position.z = 0;
    //         g.add(yxGrid);
    //
    //         var yzGridGeo = new THREE.BoxGeometry (gridWidth, gridWidth, startPositions.z);
    //         var yzGrid = new THREE.Mesh( yzGridGeo, gridMat);
    //         yzGrid.position.x = 0;
    //         yzGrid.position.y = e.end; //(yStartPos / 2);
    //         yzGrid.position.z = (startPositions.z / 2);
    //         g.add(yzGrid)
    //     });
    //
    //     xChr.forEach( function(e) {
    //         var hGridGeo = new THREE.BoxGeometry( gridWidth, gridWidth, startPositions.z );
    //         var hGrid = new THREE.Mesh( hGridGeo, gridMat );
    //         hGrid.position.x = e.end;
    //         hGrid.position.y = 0;
    //         hGrid.position.z = (startPositions.z / 2);
    //         g.add(hGrid);
    //
    //         var vGridGeo = new THREE.BoxGeometry (gridWidth, startPositions.y, gridWidth);
    //         var vGrid = new THREE.Mesh( vGridGeo, gridMat);
    //         vGrid.position.x = e.end;
    //         vGrid.position.y = (startPositions.y /2);
    //         vGrid.position.z = 0;
    //         g.add(vGrid)
    //     });
    //
    //     zChr.forEach( function(e) {
    //         var hGridGeo = new THREE.BoxGeometry( startPositions.x, gridWidth, gridWidth );
    //         var hGrid = new THREE.Mesh( hGridGeo, gridMat );
    //         hGrid.position.x = (startPositions.x / 2);
    //         hGrid.position.y = 0;
    //         hGrid.position.z = e.end;
    //         g.add(hGrid);
    //
    //         var vGridGeo = new THREE.BoxGeometry (gridWidth, startPositions.y, gridWidth);
    //         var vGrid = new THREE.Mesh( vGridGeo, gridMat);
    //         vGrid.position.x = 0;
    //         vGrid.position.y = (startPositions.y / 2);
    //         vGrid.position.z = e.end;
    //         g.add(vGrid)
    //     });
    //
    //     g.translateX(-startPositions.x / 2);
    //     g.translateY(-startPositions.y / 2);
    //     g.translateZ(-startPositions.z / 2);
    //     //scene.add(grid);
    //     return g;
    // }

    /*---------------------------------------------------------------------------------------------------------
     ~~~~INITIALIZE~~~~
     --------------------------------------------------------------------------------------------------------*/

    function initialize() {
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

        // Clear point data.
        //pointData = [];

        /* Create chromosome axes */
        var axes = drawChromosomes();
        for (var i = 0; i < axes.length; i++) {
            scene.add(axes[i])
        }

        /* TODO: Create grid */
        // grid = drawGrid();
        // scene.add(grid);

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
    var current = "ks";
    function animate() {
        /* Update controls. */
        controls.update();

        /* Update camera position record. */
        camView = { "x": camera.position.x, "y": camera.position.y, "z": camera.position.z };

        /* Check for recoloring */
        if (hCurrent[1] != current) {
            updateColors(hCurrent[1]);
            current = hCurrent[1];
        }

        /* Render the scene. */
        render();

        /* Create continuous animation. */
        requestAnimationFrame( animate );
    }

    /*---------------------------------------------------------------------------------------------------------
     ~~~~HANDLERS~~~~
     --------------------------------------------------------------------------------------------------------*/

    /* Histogram Change */
    function updateColors(select) {
        var geometry = points.geometry;
        geometry.colors = colors[select];
        geometry.colorsNeedUpdate = true;
    }
    //document.getElementById("").onclick = updateColors(?);
    //window.addEventListener()

    /* Window Resize */
    function onWindowResize() {
        width = document.getElementById("rendering").clientWidth;
        height = window.innerHeight * 0.9;
        renderer.setSize( width, height);
        camera.aspect = width / height;
        camera.updateProjectionMatrix();

        containerWidth = container.clientWidth;
        containerHeight = container.clientHeight;

        // TODO: Rerender Histograms on page resize
        d3.select("#chartSvg").remove();
        renderHistogram(hCurrent[0], histData[hCurrent[1]]);
    }
    window.addEventListener( 'resize', onWindowResize, false );

    /* Mouse Click */
    function onDocumentMouseDown( event ) {
        // Specify object of interest (points).
        var geometry = points.geometry;

        // Calculate mouse location on canvas.
        var offset_l = $("canvas").offset().left - $(window).scrollLeft();
        var offset_t = $("canvas").offset().top - $(window).scrollTop();
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

    /* TODO: Reinstate Grid Toggle */
    // var grid_state = 1;
    // function onGridToggle() {
    //     if (grid_state == 1) {
    //         grid_state = 0;
    //         scene.remove(grid);
    //     } else {
    //         grid_state = 1;
    //         scene.add(grid);
    //     }
    // }
    // var GridToggle = document.getElementById("grid_toggle");
    // GridToggle.addEventListener("click", onGridToggle);
}


$(document).ready( function() {
    // TODO: Start Spinny Wheely(s)
    $.when(loadData(final_experiment.links.graph)).done(function(data) {
        // Log data (development)
        //console.log(data);
        // TODO: End spinny wheel(s)

        // Render SynMap
        renderSynMap(data, "canvas");
        // Render Histograms
        renderHistogram("hKs", histData.ks);
    });
});
