/*-----------------------------------------------------------------------------------------------------------------
 ~~~~SELECTABLE OPTIONS~~~~
 -----------------------------------------------------------------------------------------------------------------*/

var visSelect;// = "chm"; // Can be "chm", "hco"
var knKsSelect; // = "Ks"; // Can be "Ks" or "Kn"
var comparisonSelect; // = "mean"; // Can be "mean", "xy", "xz", or "yz"

/*-----------------------------------------------------------------------------------------------------------------
 ~~~~FUNCTION DEFINITIONS~~~~
 -----------------------------------------------------------------------------------------------------------------*/

/* FUNCTION: Render Histogram*/
function renderHistogram(histogramObj) {
    function subBin(val, min, max) {
        var diff = Math.abs(max - min) / 5;
        var bins = [min + diff, min+(diff*2), min+(diff*3), min+(diff*4), max];
        if (val <= bins[0]) {
            return 0;
        } else if (val <= bins[1]) {
            return 1;
        } else if (val <= bins[2]) {
            return 2;
        } else if (val <= bins[3]) {
            return 3;
        } else if (val <= bins[4]) {
            return 4;
        } else {
            console.log("value: " + val);
            console.log("bins: " + min + '-' + max);
            console.log("Subbinning Error");
        }
    }

    function buildXaxisLabels(histogramInfo) {
        var d = Math.abs(histogramInfo['max'] - histogramInfo['min']) / 25;
        var ax = [];
        for (i = 0; i < 25; i++) {
            ax[i] = histogramInfo['min'] + 0.5*d + i*d
        }
        return ax
    }

    var plot = document.getElementById("histogram");
    var histogramInfo = histogramObj.logten.info[knKsSelect][comparisonSelect];
    var histogramData = histogramObj.logten.data[knKsSelect][comparisonSelect];
    var binCounts = [[0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0]];
    histogramData.forEach( function(value) {
        if (value <= histogramInfo.bin1[1]) {
            binCounts[0][subBin(value, histogramInfo.bin1[0], histogramInfo.bin1[1])] += 1;
        } else if (value <= histogramInfo.bin2[1]) {
            binCounts[1][subBin(value, histogramInfo.bin2[0], histogramInfo.bin2[1])] += 1;
        } else if (value <= histogramInfo.bin3[1]) {
            binCounts[2][subBin(value, histogramInfo.bin3[0], histogramInfo.bin3[1])] += 1;
        } else if (value <= histogramInfo.bin4[1]) {
            binCounts[3][subBin(value, histogramInfo.bin4[0], histogramInfo.bin4[1])] += 1;
        } else if (value <= histogramInfo.bin5[1]) {
            binCounts[4][subBin(value, histogramInfo.bin5[0], histogramInfo.bin5[1])] += 1;
        }
    });

    var xax = buildXaxisLabels(histogramInfo);
    var yax = binCounts.reduce( function(a, b) {
        return a.concat(b);
    });

    var histPlotData = [{
        x: xax,
        y: yax,
        // TODO: Make this a list comprehension so its not so repetitive.
        marker: {color: ['#2ECC71', '#2ECC71', '#2ECC71', '#2ECC71', '#2ECC71',
            '#3498DB', '#3498DB', '#3498DB', '#3498DB', '#3498DB',
            '#9B59B6', '#9B59B6', '#9B59B6', '#9B59B6', '#9B59B6',
            '#E67E22', '#E67E22', '#E67E22', '#E67E22', '#E67E22',
            '#E74C3C', '#E74C3C', '#E74C3C', '#E74C3C', '#E74C3C']},
        type: "bar"
    }];
    var histPlotLayout = {
        title: 'log10(Kn/Ks) Distribution',
        bargap: 0,
        xaxis: {title: 'log10 Ks/Kn', fixedrange: true},
        yaxis: {title: 'Counts', fixedrange: true}
    };

    Plotly.Plots.purge(plot);
    Plotly.plot(plot, histPlotData, histPlotLayout)
}

/* FUNCTION: Render Points */
//TODO: Split renderPoints into this here badboy
function renderPoints(matches, scene) {
    var dotGeo = new THREE.SphereGeometry(0.05);
            for (x_ch in matches) {
                for (y_ch in matches[x_ch]) {
                    for (z_ch in matches[x_ch][y_ch]) {
                        for (index in matches[x_ch][y_ch][z_ch]) {
                            var coordinate = matches[x_ch][y_ch][z_ch][index];
                            var xPos = coordinate[0];
                            var yPos = coordinate[1];
                            var zPos = coordinate[2];
                            var knks = coordinate[3][comparisonSelect][knKsSelect];
                            var logKnKs = knks != 0 ? Math.log10(knks) : 'NA';
                            var dotMat = new THREE.MeshBasicMaterial(
                                    {color: getColorFromHist(logKnKs, logKnKsBins)}
                            );
                            var dot = new THREE.Mesh (dotGeo, dotMat);
                            dot.position.x = xChrIndexed[x_ch]['start'] + (xPos * scaleFactor);
                            dot.position.y = yChrIndexed[y_ch]['start'] + (yPos * scaleFactor);
                            dot.position.z = zChrIndexed[z_ch]['start'] + (zPos * scaleFactor);
                            matchDots.add(dot)
                        }
                    }
                }
            }

            matchDots.translateX(-startPositions.x / 2);
            matchDots.translateY(-startPositions.y / 2);
            matchDots.translateZ(-startPositions.z / 2);
            scene.add(matchDots);
}

/* FUNCTION: Render Comparison */
function renderComparison(xid, yid, zid, xyjson, xzjson, yzjson, mergejson, histo) {
    /*---------------------------------------------------------------------------------------------------------
     ~~~~SETUP CANVAS & VR ENVIRONMENT~~~~
     --------------------------------------------------------------------------------------------------------*/
    var width = window.innerWidth * 0.689;
    var height = window.innerHeight * 0.9;

    /* Setup three.js WebGL renderer */
    var renderer = new THREE.WebGLRenderer( { antialias: true, alpha: true,
                                              canvas: document.getElementById('canvas')} );
    renderer.setSize( width, height);

    /* Create a three.js scene */
    var scene = new THREE.Scene();

    /* Create a three.js camera */
    var camera = new THREE.PerspectiveCamera( 75, width / height, 1, 100 );
    camera.position.x = 0;
    camera.position.y = 0;
    camera.position.z = 30;

    /* Apply FlatOrbitControls. */
    var controls = new THREE.FlatOrbitControls( camera );

    /*---------------------------------------------------------------------------------------------------------
     ~~~~ASSIGN VARIABLES FROM OPTIONS~~~~
     --------------------------------------------------------------------------------------------------------*/

    /* Scaling Options */
    //TODO: Scale dynamically based on size
    var scaleFactor = 0.00000001;
    var axisWidth = 0.05;

    /* Assign Files by Selected Visualization */
    var xID = xid;
    var yID = yid;
    var zID = zid;
    var xy_json = xyjson;
    var xz_json = xzjson;
    var yz_json = yzjson;
    var threeway_json = mergejson;
    var histogram = histo;

    // AKB Removed Switch from Demo, Added Above assignments
    /*var xID, yID, zID, xy_json, xz_json, yz_json, threeway_json, histogram;
    switch(visSelect) {
        case "chm":
            xID = '22736';
            yID = '25747';
            zID = '7073';
            xy_json = './demo_data/chm/parsed_22736_25747_synteny.json';
            xz_json = './demo_data/chm/parsed_22736_7073_synteny.json';
            yz_json = './demo_data/chm/parsed_25747_7073_synteny.json';
            threeway_json = './demo_data/chm/parsed_chicken_human_mouse.json';
            histogram = './demo_data/chm/parsed_chm_histogram.json';
            break;
        case "trn":
            xID = '25869';
            yID = '24777';
            zID = '20192';
            xy_json = './demo_data/trn/parsed_24777_25869_synteny.json';
            xz_json = './demo_data/trn/parsed_20192_25869_synteny.json';
            yz_json = './demo_data/trn/parsed_20192_24777_synteny.json';
            threeway_json = './demo_data/trn/parsed_thaliana_rapa_napus.json';
            histogram = './demo_data/trn/parsed_trn_histogram.json';
            break;
        case "ctz":
            xID = '22736';
            yID = '11338';
            zID = '6825';
            xy_json = './demo_data/ctz/parsed_11338_22736_synteny.json';
            xz_json = './demo_data/ctz/parsed_22736_6825_synteny.json';
            yz_json = './demo_data/ctz/parsed_11338_6825_synteny.json';
            threeway_json = './demo_data/ctz/parsed_chicken_turkey_zebrafinch.json';
            histogram = './demo_data/ctz/parsed_ctz_histogram.json';
            break;
        case "hco":
            xID = '25747';
            yID = '11691';
            zID = '9642';
            xy_json = './demo_data/hco/parsed_11691_25747_synteny.json';
            xz_json = './demo_data/hco/parsed_25747_9642_synteny.json';
            yz_json = './demo_data/hco/parsed_11691_9642_synteny.json';
            threeway_json = './demo_data/hco/parsed_human_chimp_orang.json';
            histogram = './demo_data/hco/parsed_hco_histogram.json';
            break;
        case "roc":
            xID = '24668';
            yID = '24777';
            zID = '25847';
            xy_json = './demo_data/roc/parsed_24668_24777_synteny.json';
            xz_json = './demo_data/roc/parsed_24668_25847_synteny.json';
            yz_json = './demo_data/roc/parsed_24777_25847_synteny.json';
            threeway_json = './demo_data/roc/parsed_rapa_oleraceae_capitata.json';
            histogram = './demo_data/roc/parsed_roc_histogram.json';
            break;
        default:
            console.log("No Dataset Selected");
    }*/

    /*---------------------------------------------------------------------------------------------------------
     ~~~~GLOBAL 3D OBJECT VARIABLES~~~~
     --------------------------------------------------------------------------------------------------------*/

    /* Global Object Variables */
    var xChr = new Array();
    var xAxis = new THREE.Object3D();

    var yChr = new Array();
    var yAxis = new THREE.Object3D();

    var zChr = new Array();
    var zAxis = new THREE.Object3D();

    var grid = new THREE.Object3D();
    //var gridToggle = true; TODO: Reinstate gridToggle

    var matches = new Array();
    var matchDots = new THREE.Object3D();

    var histogram_data = new Array();

    /*---------------------------------------------------------------------------------------------------------
     ~~~~FUNCTIONS~~~~
     --------------------------------------------------------------------------------------------------------*/

    /* FUNCTION: Index Chromosome Info by Name */
    function indexChr(chr) {
        var chrIndexed = {};
        chr.forEach( function(e) {
            chrIndexed[e.name] = e;
        });
        return chrIndexed
    }

    /* FUNCTION: Sort Chromosome List */
    function sortChr(a,b){
        first = parseInt(a.name);
        second = parseInt(b.name);

        if (isNaN(first) && isNaN(second)) {
            if (a.name < b.name) {
                return -1;
            } else if (a.name > b.name) {
                return 1;
            }
        } else if (isNaN(first)) {
            return 1;
        } else if (isNaN(second)) {
            return -1;
        } else {

            if(first == second) {
                if (a.name < b.name) {
                    return -1;
                } else if (a.name > b.name) {
                    return 1;
                }
            }

            else if(first < second) {
                return -1;
            }

            else if(first > second) {
                return 1;
            }}
    }

    /* FUNCTION: Build Axis from Chromosome List */
    var startPositions = {"x": 0, "y": 0, "z": 0};
    function buildAxisFromChr(chr, focusAxis) {
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

        chr.forEach(function (e) {
            var chLen = e.length * scaleFactor;

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

            var chMat;
            if (colorToggle) {
                chMat = new THREE.MeshBasicMaterial({color: colorOne});
                colorToggle = false;
            } else {
                chMat = new THREE.MeshBasicMaterial({color: colorTwo});
                colorToggle = true;
            }

            var ch = new THREE.Mesh( chGeo, chMat );

            ch.position[focusAxis] = startPositions[focusAxis] + (chLen / 2);
            ch.position[altAxes[0]] = -(axisWidth / 2);
            ch.position[altAxes[1]] = -(axisWidth / 2);

            axis.add( ch );

            e["start"] = startPositions[focusAxis];
            startPositions[focusAxis] += chLen;
            e["end"] = startPositions[focusAxis];
        });

        return axis;
    }

    /* FUNCTION: Generate Cutoff List */
    function getHistogramBins(Ks_or_Kn, logBool, comparison, histogram_object) {
        var logFactor;
        if (logBool) {
            logFactor = "logten"
        } else {
            logFactor = "values"
        }

        histInfo = histogram_object[logFactor]["info"][Ks_or_Kn][comparison];
        var binCutoffs = [histInfo.bin1[1], histInfo.bin2[1],
                          histInfo.bin3[1], histInfo.bin4[1], histInfo.bin5[1]];
        return binCutoffs;
    }

    /* FUNCTION: Get Color From Histogram */
    function getColorFromHist(value, cutoff_list) {
        var binOneColor = '#2ECC71';
        var binTwoColor = '#3498DB';
        var binThreeColor = '#9B59B6';
        var binFourColor = '#E67E22';
        var binFiveColor = '#E74C3C';
        var defaultColor = 'black';
        if (isNaN(value)) {
            return defaultColor;
        } else if (value <= cutoff_list[0]) {
            return binOneColor;
        } else if (value <= cutoff_list[1]) {
            return binTwoColor;
        } else if (value <= cutoff_list[2]) {
            return binThreeColor;
        } else if (value <= cutoff_list[3]) {
            return binFourColor;
        } else if (value <= cutoff_list[4]) {
            return binFiveColor;
        } else {
            return defaultColor;
        }
    }

    /*---------------------------------------------------------------------------------------------------------
     ~~~~GET DATA & CREATE 3D OBJECTS~~~~
     --------------------------------------------------------------------------------------------------------*/
    function loadXY() {
        return $.ajax({
                       type: 'GET',
                       //isLocal: true,
                       url: xy_json,
                       dataType: 'json',
                       success: function(data) {
                           xChr = data["genomes"][xID]["chromosomes"];
                           yChr = data["genomes"][yID]["chromosomes"];
                       },
                       error: function(xhr, status, error) {
                           console.warn("XY JSON Error: " + error);
                       }
                     })
    }

    function loadXZ() {
        return $.ajax({
                       type: 'GET',
                       url: xz_json,
                       dataType: 'json',
                       success: function(data) {
                           zChr = data["genomes"][zID]["chromosomes"];
                       },
                       error: function(xhr, status, error) {
                           console.warn("XZ JSON Error: " + error);
                       }
                     })
    }
    //TODO: Combine XY and XZ genome/chromosome info into one file (faster loading) - generate w/ three_dots_merge.py

    function loadThree() {
        return $.ajax({
                       type: 'GET',
                       url: threeway_json,
                       dataType: 'json',
                       success: function(data) {
                           matches = data;
                       },
                       error: function(xhr, status, error) {
                           console.warn("Threeway JSON Error: " + error);
                       }
                     })
    }

    function loadHistogram() {
        return $.ajax({
                       type: 'GET',
                       url: histogram,
                       dataType: 'json',
                       success: function(data) {
                           histogram_data = data;
                       },
                       error: function(xhr, status, error) {
                           console.log("Histogram JSON Error: " + error);
                       }
                     })
    }

    $.when( loadXZ(), loadXY(), loadThree(), loadHistogram()).done( function() {
        /* Sort and Index Chromosomes */
        xChr.sort(sortChr);
        yChr.sort(sortChr);
        zChr.sort(sortChr);
        var xChrIndexed = indexChr(xChr);
        var yChrIndexed = indexChr(yChr);
        var zChrIndexed = indexChr(zChr);

        /* Establish Color Bin Cutoffs */
        var logKnKsBins = getHistogramBins(knKsSelect, true, comparisonSelect, histogram_data);

        /* 3d Object: Axes */
        xAxis = buildAxisFromChr(xChr, 'x');
        yAxis = buildAxisFromChr(yChr, 'y');
        zAxis = buildAxisFromChr(zChr, 'z');

        xAxis.translateX(-startPositions.x / 2);
        xAxis.translateY(-startPositions.y / 2);
        xAxis.translateZ(-startPositions.z / 2);

        yAxis.translateX(-startPositions.x / 2);
        yAxis.translateY(-startPositions.y / 2);
        yAxis.translateZ(-startPositions.z / 2);

        zAxis.translateX(-startPositions.x / 2);
        zAxis.translateY(-startPositions.y / 2);
        zAxis.translateZ(-startPositions.z / 2);

        scene.add(xAxis);
        scene.add(yAxis);
        scene.add(zAxis);

        /* 3d Object: Grid */
        var gridWidth = 0.025;
        var gridMat = new THREE.MeshBasicMaterial({color: 'grey'});

        yChr.forEach( function(e) {
            var yxGridGeo = new THREE.BoxGeometry( startPositions.x, gridWidth, gridWidth );
            var yxGrid = new THREE.Mesh( yxGridGeo, gridMat );
            yxGrid.position.x = (startPositions.x / 2);
            yxGrid.position.y = e.end;
            yxGrid.position.z = 0;
            grid.add(yxGrid);

            var yzGridGeo = new THREE.BoxGeometry (gridWidth, gridWidth, startPositions.z);
            var yzGrid = new THREE.Mesh( yzGridGeo, gridMat);
            yzGrid.position.x = 0;
            yzGrid.position.y = e.end; //(yStartPos / 2);
            yzGrid.position.z = (startPositions.z / 2);
            grid.add(yzGrid)
        });

        xChr.forEach( function(e) {
            var hGridGeo = new THREE.BoxGeometry( gridWidth, gridWidth, startPositions.z );
            var hGrid = new THREE.Mesh( hGridGeo, gridMat );
            hGrid.position.x = e.end;
            hGrid.position.y = 0;
            hGrid.position.z = (startPositions.z / 2);
            grid.add(hGrid);

            var vGridGeo = new THREE.BoxGeometry (gridWidth, startPositions.y, gridWidth);
            var vGrid = new THREE.Mesh( vGridGeo, gridMat);
            vGrid.position.x = e.end;
            vGrid.position.y = (startPositions.y /2);
            vGrid.position.z = 0;
            grid.add(vGrid)
        });

        zChr.forEach( function(e) {
            var hGridGeo = new THREE.BoxGeometry( startPositions.x, gridWidth, gridWidth );
            var hGrid = new THREE.Mesh( hGridGeo, gridMat );
            hGrid.position.x = (startPositions.x / 2);
            hGrid.position.y = 0;
            hGrid.position.z = e.end;
            grid.add(hGrid);

            var vGridGeo = new THREE.BoxGeometry (gridWidth, startPositions.y, gridWidth);
            var vGrid = new THREE.Mesh( vGridGeo, gridMat);
            vGrid.position.x = 0;
            vGrid.position.y = (startPositions.y / 2);
            vGrid.position.z = e.end;
            grid.add(vGrid)
        });

        grid.translateX(-startPositions.x / 2);
        grid.translateY(-startPositions.y / 2);
        grid.translateZ(-startPositions.z / 2);
        scene.add(grid);

        /* 3d Object: 3-Way Comparison */
        function drawPoints() {
            var dotGeo = new THREE.SphereGeometry(0.05);
            for (x_ch in matches) {
                for (y_ch in matches[x_ch]) {
                    for (z_ch in matches[x_ch][y_ch]) {
                        for (index in matches[x_ch][y_ch][z_ch]) {
                            var coordinate = matches[x_ch][y_ch][z_ch][index];
                            var xPos = coordinate[0];
                            var yPos = coordinate[1];
                            var zPos = coordinate[2];
                            var knks = coordinate[3][comparisonSelect][knKsSelect];
                            var logKnKs = knks != 0 ? Math.log10(knks) : 'NA';
                            var dotMat = new THREE.MeshBasicMaterial(
                                    {color: getColorFromHist(logKnKs, logKnKsBins)}
                            );
                            var dot = new THREE.Mesh (dotGeo, dotMat);
                            dot.position.x = xChrIndexed[x_ch]['start'] + (xPos * scaleFactor);
                            dot.position.y = yChrIndexed[y_ch]['start'] + (yPos * scaleFactor);
                            dot.position.z = zChrIndexed[z_ch]['start'] + (zPos * scaleFactor);
                            matchDots.add(dot)
                        }
                    }
                }
            }

            matchDots.translateX(-startPositions.x / 2);
            matchDots.translateY(-startPositions.y / 2);
            matchDots.translateZ(-startPositions.z / 2);
            scene.add(matchDots);
        }
        drawPoints();

        /* Histogram */
        renderHistogram(histogram_data);

    });

    /*---------------------------------------------------------------------------------------------------------
     ~~~~ANIMATE LOOP~~~~
     --------------------------------------------------------------------------------------------------------*/

    function animate() {
        /* Update controls. */
        controls.update();

        /* Render the scene. */
        renderer.render( scene, camera );

        /* Create continuous animation. */
        requestAnimationFrame( animate );
    }
    animate();

    /*---------------------------------------------------------------------------------------------------------
     ~~~~WINDOW RESIZE HANDLER~~~~
     --------------------------------------------------------------------------------------------------------*/

    function onWindowResize() {
        width = window.innerWidth * 0.689;
        height = window.innerHeight * 0.9;
        renderer.setSize( width, height);
        camera.aspect = width / height;
        camera.updateProjectionMatrix();
        //renderHistogram(histogram_data); TODO: Re-render histogram on page resize, but without blanking out
    }
    window.addEventListener( 'resize', onWindowResize, false );
} //END renderIt()

/*-----------------------------------------------------------------------------------------------------------------
 ~~~~ POPULATE PAGE ~~~~
 -----------------------------------------------------------------------------------------------------------------*/
$(document).ready( function() {
    /* Monitor Visualization Selector */
    var VisSelector = $("#dataset_select");
    visSelect = VisSelector.val();
    VisSelector.change( function() {
        visSelect = VisSelector.val();
        renderComparison();
    });

    /* Monitor KnKsSelect */
    var KnKsSelector = $("#knks_select");
    knKsSelect = KnKsSelector.val();
    KnKsSelector.change( function() {
        knKsSelect = KnKsSelector.val();
        renderComparison();
    });

    /* Monitor Comparison Select */
    var ComparisonSelector = $("#comparison_select");
    comparisonSelect = ComparisonSelector.val();
    ComparisonSelector.change( function() {
        comparisonSelect = ComparisonSelector.val();
        renderComparison();
    });

    /* Render Visualization on Page Load */
    renderComparison(
        final_experiment.x_gid,
        final_experiment.y_gid,
        final_experiment.z_gid,
        final_experiment.links.xy_json,
        final_experiment.links.xz_json,
        final_experiment.links.yz_json,
        final_experiment.links.merge,
        final_experiment.links.histo
    );

}); //END $(document).ready()
