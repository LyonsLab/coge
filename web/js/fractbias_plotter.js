/**
 * Created by senorrift on 7/8/16.
 */

// NOTE: Make sure to include these three sources (Plot.ly, D3 (v.4.1.0), and jQuery 3)
// <script type="text/javascript" src="https://cdn.plot.ly/plotly-latest.min.js"></script>
// <script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/d3/4.1.0/d3.min.js"></script>
// <script type="text/javascript" src="https://code.jquery.com/jquery-3.0.0.min.js" integrity="sha256-JmvOoLtYsmqlsWxa7mDSLMwa6dZ9rrIdtrrVYRnDRH0=" crossorigin="anonymous"></script>

// NOTE: Make sure to have the following three classes defined.
// <style type="text/css">
//     .leftPlot {
//         float: left;
//     }
//     .rightPlot {
//         float: right;
//     }
//     .plot {
//         width: 48.5%;
//         text-align: center;
//         /*height: 300px;*/
//     }
// </style>

// data_json originally "parsed_output_dict.json"

function plotFractBias(parentdivid, data_json, target_genome) {
    var fbDiv = $("#" + parentdivid);
    var iterList = function(cap) {
        var l = [];
        for (var z=1; z<=cap; z++) { l.push(z) };
        return l;
    };
    var add = function(a, b) {
        return a + b;
    };

    // Load data
    d3.json(data_json, function(error, data) {
        // List target chromosomes
        var targetChromosomes = Object.keys(data);

        // Create color scale for chromosomes
        //var c20 = d3.scaleOrdinal(d3.schemeCategory20); // https://github.com/d3/d3/blob/master/CHANGES.md#scales-d3-scale
        var c20 = d3.scale.category20();

        // Add grid of divs to hold each plot
        var rowCount = Math.ceil(targetChromosomes.length / 2);
        for (var r=0; r < (rowCount*2); r+=2) {
            fbDiv.append("<div style='clear:both'><div id='plot" + r + "'></div><div id='plot" + (r+1) + "'></div></div>");
            $("#plot" + r).addClass("leftPlot plot");
            $("#plot" + (r+1)).addClass("rightPlot plot");
        }

        for (var t=0; t<targetChromosomes.length; t++) {
            var plotDiv = document.getElementById("plot" + t);
            var tChr = targetChromosomes[t];
            var qChrs = Object.keys(data[tChr]);

            var layout = {
                title: target_genome + '<br>Target Chromosome: ' + tChr,
                xaxis: { title: "Window Iteration (Gene Number)" },
                yaxis: { title: "% Retention<br>(# Syntenic Genes/Window Size" },
                annotations: [{
                    text: "Query Chr.",
                    xref: 'paper', x: 1.02, xanchor: 'left',
                    yref: 'paper', y: 1, yanchor: 'bottom',
                    showarrow: false
                }]
            };

            var traces = [];
            for (var q=0; q<qChrs.length; q++) {
                var qChr = qChrs[q];
                var trackData = data[tChr][qChr];
                var sum = trackData.reduce(add, 0);
                if (sum > 0) {
                    var trace = {
                        x: iterList(trackData.length),
                        y: trackData,
                        name: "Chr. " + qChr,
                        marker: { color: c20(qChr) }
                    };
                    traces.push(trace);
                }

            }

            // Plot
            Plotly.plot(plotDiv, traces, layout);
        }

    });
}