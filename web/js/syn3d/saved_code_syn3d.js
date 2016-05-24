/**
 * Created by senorrift on 1/22/16.
 */

/* FUNCTION: Render Histogram DEPRICATED (see new version below)*/
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
    Plotly.plot(plot, histPlotData, histPlotLayout, {showLink: false, displayModeBar: false});
}



/* 3d Object: 3-Way Comparison */
function drawPoints() {
    var calcPtColor = getPtColorFunc(histogram_data.logten.data["Ks"][comparisonSelect], JetColorScheme);
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
                    //OLD dotMat
                    //var dotMat = new THREE.MeshBasicMaterial(
                    //        {color: getColorFromHist(logKnKs, logKnKsBins)}
                    //);
                    // NEW dotMat
                    var dotMat = new THREE.MeshBasicMaterial(
                        {color: calcPtColor(logKnKs)}
                    );
                    var dot = new THREE.Mesh(dotGeo, dotMat);
                    dot.position.x = xChrIndexed[x_ch]['start'] + (xPos * scaleFactor);
                    dot.position.y = yChrIndexed[y_ch]['start'] + (yPos * scaleFactor);
                    dot.position.z = zChrIndexed[z_ch]['start'] + (zPos * scaleFactor);
                    //dot.on('click', function() { console.log(coordinate[4]); } );
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