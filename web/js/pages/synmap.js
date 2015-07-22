//TODO: Create a proper helper function
function get_gc (dsgid, divid)
{
    $('#'+divid).removeClass('link').html('loading...');
    get_dsg_gc(['args__dsgid','args__'+dsgid,'args__text','args__1'],[divid]);
}

//TODO: Replace with proper promise chain
function get_organism_chain(type,val,i)
{
	$('#feattype_menu'+i).hide();
	$('#dsgid'+i).html('');
	
    $('#org_list').html('<input type=hidden id = "org_id"+i><font class="loading"></font>');
  //  if (type == 'name')
    	get_orgs(['args__search','args__'+val,'args__i','args__'+i], ['org_list'+i]);
//    else if (type == 'desc') 
//   	get_orgs(['args__desc','args__'+val,'args__i','args__'+i], ['org_list'+i]);
    //$('#dsg_info'+i).html('<div class="loading dna_small small">loading. . .</div>');
    $('#dsg_info'+i).html('<div class="small note indent">loading... <img src="picts/ajax-loader.gif"/></div>');
    ajax_wait("gen_dsg_menu(['args__oid','org_id"+i+"', 'args__num','args__"+i+"'],['dsg_menu"+i+"', 'genome_message"+i+"']);");
    ajax_wait("get_genome_info(['args__dsgid','dsgid"+i+"','args__org_num','args__"+i+"'],[handle_dsg_info]);");
}

//TODO: Replace with proper promise chain
function get_genome_info_chain(i) {
	//$('#dsg_info'+i).html('<div class=dna_small class=loading class=small>loading. . .</div>');
	$('#dsg_info'+i).html('<div class="small note indent">loading... <img src="picts/ajax-loader.gif"/></div>');
    // ajax_wait("gen_dsg_menu(['args__oid','org_id"+i+"', 'args__num','args__"+i+"'],['dsg_menu"+i+"','genome_message"+i+"']);");
    gen_dsg_menu(['args__oid','org_id'+i, 'args__num','args__'+i],['dsg_menu'+i, 'genome_message'+i]);
    $('#depth_org_1').html($('#org_id1 option:selected').html());
    $('#depth_org_2').html($('#org_id2 option:selected').html());

    ajax_wait("get_genome_info(['args__dsgid','dsgid"+i+"','args__org_num','args__"+i+"'],[handle_dsg_info]);");
}

function search_bar(text, divId) {
    console.log(divId);
    var el = $(divId);
    el.val(el.val() + " " + text);
}

function close_dialog(dialog) {
    dialog.dialog('close');
    dialog.find('#text').empty();
    dialog.find('#progress').show();
    dialog.find('#dialog_error').hide();
    dialog.find('#dialog_success').hide();
}

function load_results() {
    $('#intro').hide();
    $('#log_text').hide();
    $('#results').fadeIn();
}

function update_params(val) {
    var cmd;
    var params;
    var type;

    if (val) {
        params = val.split('_');
    } else {
        params = $('#prev_params').val()[0].split('_');;
    }

    if ($('#org_id1').val() == params[3]) {
        cmd = "$('#feat_type1').attr('value', '"+params[5]+"');$('#feat_type2').attr('value', '"+params[8]+"');";
        $('#dsgid1').attr('value',params[4]);
        $('#dsgid2').attr('value',params[7]);
    } else {
        cmd = "$('#feat_type2').attr('value', '"+params[5]+"');$('#feat_type1').attr('value', '"+params[8]+"');";
        $('#dsgid2').attr('value',params[4]);
        $('#dsgid1').attr('value',params[7]);
    }

    get_genome_info(['args__dsgid','dsgid1','args__org_num','args__1'],[handle_dsg_info]);
    get_genome_info(['args__dsgid','dsgid2','args__org_num','args__2'],[handle_dsg_info]);
    ajax_wait(cmd);

    $('#blast').attr('value',params[9]);

    if (params[10] == 'Distance') {
        $("input[name='dagchainer_type']:nth(1)").attr("checked","checked");
        type=" bp";
    } else {
        $("input[name='dagchainer_type']:nth(0)").attr("checked","checked");
        type= " genes";
    }

    display_dagchainer_settings([params[1],params[2]],type);
    $('#c').val(params[11]);
    merge_select_check();
    depth_algo_check();
}    

function handle_dsg_info(dsg_html, feat_menu, genome_message, length, org_num, org_name, seq_id) {
    $('#dsg_info'+org_num).html(dsg_html);
    
    $('#feattype_menu'+org_num).html(feat_menu);
    if (dsg_html)
    	$('#feattype_menu'+org_num).show();
    
    $('#genome_message'+org_num).html(genome_message);

    if (org_num == 1) {
        pageObj.org_length1 = length;
        pageObj.org_name1 = org_name;
        pageObj.seq_type1 = seq_id;
    } else {
        pageObj.org_length2 = length;
        pageObj.org_name2 = org_name;
        pageObj.seq_type2 = seq_id;
    }
}

function set_dagchainer_defaults(params, type) {
    var settings = $('#dagchainer_default').val();

    if (!(params && type)) {
        if ($('#dagchainer_type')[0].checked) {
            params = [20,5,0,0];
            type = " genes";
        } else {
            if (settings == 1) { // for plant genomes
                params = [120000, 5, 96000,480000];
            } else if (settings == 2) { // for microbe genomes
                params = [2000, 5, 4000, 8000];
            }

            type=" bp";
        }
    }

    if (!params) {
        return;
    }

    $('#D').val(params[0]);
    $('#A').val(params[1]);

    if (typeof(params[2]) == 'undefined') {
        params[2] = 4*params[0];
    }

    if (typeof(params[3]) == 'undefined') {
        params[3] = 4*params[1];
    }

    $('#gm').val(params[2]);
    $('#Dm').val(params[3]);
    $('.distance_type').html(type);
}

function ajax_wait (val){
    if (ajax.length) {
        setTimeout("ajax_wait("+'"'+val+'"'+")",100);
        return;
    }

    eval(val);
}

function timing(val, val2){
    var searchterm;
    namere = /name/;
    descre = /desc/;

    if (namere.exec(val)) {
        searchterm = $('#'+val).val();
    } else if (descre.exec(val)) {
        searchterm = $('#'+val).val();
    }

    if (!searchterm) {
        val=0;
    }

    if(searchterm == "Search") {
        searchterm = "";
    }

    pageobjsearch = "search"+val;
    pageobjtime = "time"+val;

    if (pageObj.pageobjsearch && pageObj.pageobjsearch == searchterm+val) {
    //    return;
    }

    pageObj.pageobjsearch=searchterm+val;

    if (pageObj.pageobjtime){
        clearTimeout(pageObj.pageobjtime);
    }

    re = /(\d+)/;
    i = re.exec(val);

    if (namere.exec(val)) {
        if (val2) {
            get_organism_chain('search',$('#'+val).val(),i[0])
        } else {
            pageObj.pageobjtime = setTimeout("get_organism_chain('search',$('#"+val+"').val(),i[0])",500);
        }
    } else if (descre.exec(val)) {
        if (val2) {
            get_organism_chain('search',$('#'+val).val(),i[0])
        } else {
            pageObj.pageobjtime = setTimeout("get_organism_chain('search',$('#"+val+"').val(),i[0])",200);
        }
    }
}

function display_dagchainer_settings(params,type) {

    if ($('#dagchainer_type')[0].checked) {
        $('#dagchainer_distance').hide(0);
    } else {
        $('#dagchainer_distance').show(0);
    }

    set_dagchainer_defaults(params, type);
}

function address_validity_check(validity) {
    if (validity) {
        if(validity == 'invalid') {
            $('#email_error').show(0);
        } else {
            $('#email_error').hide(0);
        }
    } else {
        check_address_validity(['email'],[address_validity_check]);
    }
}

function fill_jobtitle(){
    var title;
    var org1 = $('#org_id1 option:selected').html() || 0;
    var org2 = $('#org_id2 option:selected').html() || 0;

    if (org1 != 0) {
        org1 = org1.replace(/\s+\(id\d+\)$/,"");
    }

    if (org2 != 0) {
        org2 = org2.replace(/\s+\(id\d+\)$/,"");
    }

    if (org1 != 0 && org2 != 0) {
        title = org1 + " v. " + org2;
    } else if (org1 != 0) {
        title = org1;
    } else if (org2 != 0) {
        title = org2;
    } else {
        return;
    }

    $('#jobtitle').val(title);
}

function update_basename(basename){
    pageObj.basename=basename;
}

function reset_basename(){
    if(pageObj.basename) pageObj.basename=0;
}

function synteny_zoom(dsgid1, dsgid2, basename, chr1, chr2, ksdb) {
    var url = 'dsg1='+dsgid1+';dsg2='+dsgid2+';chr1='+chr1+';chr2='+chr2+';base='+basename;
    var loc = $('#map_loc').val();
    var width = $('#zoom_width').val();
    var min = $('#zoom_min').val();
    var max = $('#zoom_max').val();
    var am = $('#axis_metric').val();
    var fid1=0;
    if (pageObj.fid1) {fid1 = pageObj.fid1;}
    var fid2=0;
    if (pageObj.fid2) {fid2 = pageObj.fid2;}
    var ct = $('#color_type').val();
    var loc = pageObj.loc;

    if (!loc) {loc=1;}

    loc++;
    pageObj.loc=loc;

    win = window.open ('DisplayMessage.pl', 'win'+loc,'width=400,height=200,scrollbars=1');
    win.focus();

    get_dotplot(
        ['args__url','args__'+url, 'args__loc','args__'+loc, 'args__flip','args__'+$('#flip')[0].checked,'args__regen_images','args__'+$('#regen_images')[0].checked, 'args__width', 'args__'+width, 'args__ksdb','args__'+ksdb,'args__kstype','ks_type','args__min', 'args__'+min,'args__max', 'args__'+max, 'args__am', 'args__'+am, 'args__ct','args__'+ct, 'args__bd', 'args__'+$('#box_diags')[0].checked, 'args__color_scheme','color_scheme', 'args__am','axis_metric', 'args__ar','axis_relationship', 'args__fid1','args__'+fid1, 'args__fid2', 'args__'+fid2],[open_window]);
}

function open_window (url, loc, width, height) {
    if (!loc) {
        loc = pageObj.loc;
    }

    if (!loc) {
        loc=1;
    }
    my_window = window.open(url,"win"+loc,'"width='+width+',height='+height+', scrollbars=1"');
    my_window.resizeTo(width,height);
}

function merge_select_check () {
    var merge_algo = $('#merge_algo').val();

    if (merge_algo == 0) {
        $('#merge_algo_options').hide();
    } else if (merge_algo == 1) {
        $('#merge_algo_options').show();
        $('#max_dist_merge').hide();
    } else {
        $('#merge_algo_options').show();
        $('#max_dist_merge').show();
    }
}

function depth_algo_check() {
   var depth_algo = $('#depth_algo').val();

    if (depth_algo == 0) {
        $('#depth_options').hide();
    } else if (depth_algo == 1) {
        $('#depth_options').show();
    }
}

function post_to_grimm(seq1, seq2) {
    var url = "http://nbcr.sdsc.edu/GRIMM/grimm.cgi#report";
    var query_form = document.createElement("form");
    var input1 = document.createElement("textarea");
    var input2 = document.createElement("textarea");

    seq1 = seq1.replace(/\|\|/g,"\n");
    seq2 = seq2.replace(/\|\|/g,"\n");

    input1.name="genome1";
    input1.value=seq1;

    input2.name="genome2";
    input2.value=seq2;

    query_form.method="post" ;
    query_form.action=url;
    query_form.setAttribute("target", "_blank");
    query_form.setAttribute("name", "genomeForm");
    query_form.appendChild(input1);
    query_form.appendChild(input2);
    query_form.submit("action");
}

function add(a, b) {
    return a + b;
};

function sum(iterable) {
    return iterable.reduce(add);
}

function pick(attribute) {
    return function(object, index, list) {
        return object[attribute];
    };
};

function scan(func, initial, iterable) {
    var values = [initial],
        index;

    for(index = 0; index < iterable.length; index++) {
        values.push(func(values[index], iterable[index]));
    }

    return values;
}

function toObject(pairs) {
    return pairs.reduce(function(object, value, index) {
        object[value[0]] = value[1];
        return object
    }, {});
}

function compareAlphaNumeric(a, b) {
    var regAlpha = /[^a-zA-Z_]/g,
        regNumeric = /[^0-9]/g,
        aString = a.toString(),
        bString = b.toString(),
        aAlpha,
        aNum,
        bAlpha,
        bNum;

    aAlpha = aString.replace(regAlpha, "");
    bAlpha = bString.replace(regAlpha, "");

    if(aAlpha === bAlpha) {
        aNum = parseInt(aString.replace(regNumeric, ""), 10);
        bNum = parseInt(bString.replace(regNumeric, ""), 10);

        return aNum === bNum ? 0 : aNum > bNum ? 1 : -1;
    }

    return aAlpha.localeCompare(bAlpha);
}

function sortBy(attribute, cmp) {
    return function (a, b) {
        if (cmp) {
            return cmp(a[attribute], b[attribute]);
        } else {
            return a[attribute] - b[attribute];
        }
    };
}

function identity(x) {
    return x;
}

function inverse(func) {
    return function(a, b) {
        return func(b, a);
    }
}

function getValues(object) {
    var key,
        coll = [];

    for(key in object) {
        coll.push(object[key]);
    }

    return coll;
}

function checkRequestSize(url) {
    return $.ajax({
        type: "HEAD",
        url: url
    });
}

(function(root, $, _) {
    var synmap = root.synmap = {};

    var ArrayProto = Array.prototype.slice;
    var slice = ArrayProto.slice;

    var lastElement = function(x) { return _.last(x); };

    var Resizable = synmap.Resizable = function(view, classes) {
        var el = $(view.el()),
            resize = _.debounce(view.resize, 300);

        classes || (classes = "ui-widget-content ui-corner-all");

        el.resizable({
            minWidth: 300,
            minHeight: 300,
            resize: function(e, ui) {
                resize(ui.size);
            }
        }).addClass(classes);

        return view;
    };

    var Draggable = synmap.Draggable = function(view) {
        var el = $(view.el());

        el.draggable({
            handle: "h3"
        });

        return view;
    };

    var PlotViewer = synmap.PlotViewer = function(element, metric) {
        var genomes,
            pairs = [],
            plots = [],
            sort,
            filter = Filter("syntenic_pairs"),
            builder = PlotBuilder(),
            changedEvent = Event();
            plotElement = $("<div></div>", { id: _.uniqueId("plot") });
            my = {};

        plotElement.appendTo(element);

        my.onChanged = function(callback) {
            changedEvent.attach(callback);
        };

        my.el = function(_) {
            return element;
        };

        my.loadPlots = function(json) {
            var keys,
                type,
                layer,
                source,
                reference;

            for(layer in json.layers) {
                type = _.keys(json.layers[layer].data)[0];
                for (reference in json.layers[layer].data[type]) {
                    for (source in json.layers[layer].data[type][reference]) {
                        pairs.push([reference, source]);
                    }
                }
            }

            keys = pairs.pop();

            // Create a new dotplot
            builder.loadJSON(json);
            builder.setAxisMetric(metric);
            builder.setChromosomeSort(sort);
            builder.setReference(keys[0]);
            builder.setSource(keys[1]);
            data = builder.get();

            var handler =  function() {
                var active = [],
                    index,
                    layer,
                    new_layer,
                    lines,
                    rects,
                    points,
                    ids = _.keys(data.layers).reverse();

                for (index = 0; index < ids.length; index++) {
                    if (data.layers[ids[index]].enabled) {
                        layer = data.layers[ids[index]];
                        new_layer = _.clone(layer);

                        lines = layer.lines || {};
                        rects = layer.rects || {};
                        points = layer.points || {};

                        if (ids[index] === filter.name && !filter.empty()) {
                            new_layer.lines = _.values(_.pick(lines, filter.indices));
                            new_layer.rects = _.values(_.pick(rects, filter.indices));
                            new_layer.points = _.values(_.pick(points, filter.indices));
                        } else {
                            new_layer.lines = _.values(lines);
                            new_layer.rects = _.values(rects);
                            new_layer.points = _.values(points);
                        }

                        active.push(new_layer);
                    }
                }

                return active;
            };

            plotElement.empty();
            var plot = new DotPlot(plotElement.attr("id"), {
                size: { width: 1000, height: 800 },
                genomes: [{
                    name: data.xtitle,
                    length: data.xtotal,
                    chromosomes: data.xlabels,
                    fetchDataHandler: handler
                },
                {
                    name: data.ytitle,
                    length: data.ytotal,
                    chromosomes: data.ylabels,
                    fetchDataHandler: handler
                }],
                fetchDataHandler: handler,
                style: {
                    position: "relative",
                    minHeight: "800px"
                }
            });

            plot.redraw();

            plots.length = 0;
            plots.push({
                "data": data,
                "plot": plot
            });

            changedEvent.notify();
        };

        my.setSort = function(func) {
            sort = func;
        };

        my.toggleLayer = function(layerId, enabled) {
            var index;

            for(index = 0; index < plots.length; index++) {
                if (plots[index].data.layers[layerId]) {
                    plots[index].data.layers[layerId].enabled = enabled;
                    plots[index].plot.redraw();
                }
            }

            changedEvent.notify();
        }

        my.filter = function(indices) {
            filter.clear();
            filter.add(indices);
            my.update();
        };

        my.update = function() {
            for(index = 0; index < plots.length; index++) {
                plots[index].plot.redraw();
            }
        };

        my.getLayer = function(layerId) {
            return plots[0].data.layers[layerId];
        };

        my.reset = function() {
            var index;

            for(index = 0; index < plots.length; index++) {
                plots[index].plot.reset();
            }
        };

        return my;
    };

    var PlotBuilder = synmap.PlotBuilder = function() {
        var my = {},
            genomes = {},
            layers = {},
            by = undefined,
            metric = "nucleotides",
            source,
            reference;

        my.setChromosomeSort = function(func) {
            by = func;
        };

        my.setAxisMetric = function(new_metric) {
            metric = new_metric;
        };

        my.setSource = function(src) {
            source = src;
        }

        my.setReference = function(ref) {
            reference = ref;
        };

        my.addGenome = function(id, data) {
            genomes[id] = data;
        };

        my.addLayer = function(id, data) {
            layers[id] = data;
        };

        my.loadJSON = function(json) {
            // Add genomes
            _.each(json.genomes, function(data, id) {
                my.addGenome(id, data);
            });

            // Add layers
            _.each(json.layers, function(layer, id) {
                my.addLayer(id, layer);
            });
        };

        my.get = function() {
            var xlabels,
                ylabels,
                xtotal,
                ytotal,
                new_layers,
                genome1 = genomes[reference],
                genome2 = genomes[source];

            // Check if we have any genomes
            if (genome1 === undefined || genome2 === undefined) {
                return;
            };

            // Construct labels with property axis and offsets
            // x is flipped to sort by maximum value
            // y is not flipped because increasing value is down in canvas
            xlabels = generateLabels(genome1.chromosomes, metric, inverse(by));
            ylabels = generateLabels(genome2.chromosomes, metric, inverse(by));

            // Generate the layers
            new_layers = createAllLayers(xlabels, ylabels);

            return {
                "xtitle": genome1.name,
                "ytitle" : genome2.name,
                "xlabels": xlabels,
                "ylabels": ylabels,
                "xtotal": sum(xlabels.map(pick("length"))),
                "ytotal": sum(ylabels.map(pick("length"))),
                "layers": new_layers
            };
        };

        function generateLabels(chromosomes, metric, by) {
            return orderBy(getLengths(chromosomes, metric), by);
        }

        function getLengths(chromosomes, prop) {
            return chromosomes.map(function(x) {
                return { name : x.name, length: x[prop] };
            });
        }

        function orderBy(chromosomes, by) {
            return chromosomes.sort(by);
        };

        function offsets(chromosomes) {
            return scan(add, 0, chromosomes.map(pick("length")));
        }

        function generateIndexBy(collection, by) {
            var index = {};

            for(i = 0; i < collection.length; i++) {
                index[by(collection[i])] = i;
            };

            return index;
        }

        function createAllLayers(xlabels, ylabels) {
            var xindex = generateIndexBy(xlabels, pick("name")),
                yindex = generateIndexBy(ylabels, pick("name")),
                xoffsets = offsets(xlabels),
                // For Canvas coordinates the chromosome should be offset
                // starting from the end of its length
                yoffsets = offsets(ylabels),//.slice(1),
                output = {},
                create,
                layerId,
                rawLayers,
                transformedLayers;

            create = createLayer.bind(null, xoffsets, yoffsets, xindex, yindex);

            layerIds = _.keys(layers);
            rawLayers = getValues(layers);
            transformedLayers = rawLayers.map(create);

            return _.object(_.zip(layerIds, transformedLayers));
        }

        function createLayer(xoffset, yoffset, xindex, yindex, layer) {
            var style = layer.style || (style = {}),
                enabled = layer.enabled || false,
                points,
                rects,
                lines,
                raw;

            if (layer.data.points) {
                raw = layer.data.points[reference][source];
                points = transformBy(raw, xoffset, yoffset, xindex, yindex, transformPoint);
            }

            if (layer.data.lines) {
                raw = layer.data.lines[reference][source]
                lines = transformBy(raw, xoffset, yoffset, xindex, yindex, transformLine);
            }

            if (layer.data.rects) {
                raw = layer.data.rects[reference][source];
                rects = transformBy(raw, xoffset, yoffset, xindex, yindex, transformRect);
            }

            return {
                "enabled": enabled,
                "style": style,
                "points": points,
                "lines": lines,
                "rects": rects
            };
        }

        function transformBy(gridLayer, rowOffset, colOffset, rowIndex, colIndex, transform) {
            var row,
                col,
                r,
                c,
                data,
                transformed,
                collection = [];

            for(row in gridLayer) {
                for (col in gridLayer[row]) {
                    r = rowOffset[rowIndex[row]];
                    c = colOffset[colIndex[col]];
                    data = gridLayer[row][col];

                    if (data !== undefined)  {
                        indexes = _.keys(data);
                        values = _.values(data);
                        transformed = values.map(transform.bind(undefined, r, c));
                        collection = [].concat.call(collection, _.zip(indexes, transformed));
                    }
                }
            }

            return _.object(collection);
        }

        function transformPoint(rowoffset, colOffset, point) {
            return {
                x: point[0] + rowOffset,
                y: point[1] + colOffset //toCanvasRow(colOffset, point[1]), // mdb
            };
        }

        function transformLine(rowOffset, colOffset, line) {
            return {
                x1: line[0] + rowOffset,
                x2: line[1] + rowOffset,
                y1: line[2] + colOffset, //toCanvasRow(colOffset, line[2]), // mdb
                y2: line[3] + colOffset  //toCanvasRow(colOffset, line[3])  // mdb
            };
        }

        function transformRect(rowOffset, colOffset, rect) {
            var width = Math.abs(rect[0] - rect[1]),
                height = Math.abs(rect[2] - rect[3]),
                x = Math.min(rect[0], rect[1]);
                y = Math.min(rect[2], rect[3]);

            return [x + rowOffset, y + colOffset, width, height];
        }

        return my;
    };

    var Event = function() {
        var listeners = [];
        return {
            attach: function(observer) {
                listeners.push(observer);
            },
            notify: function(args) {
                for(var i = 0; i < listeners.length; i++) {
                    listeners[i](args);
                }
            }
        };
    };

    var Histogram = synmap.Histogram = function(element, model, config) {
        var histogram = bioplot.histogram(),
            selectionChange = Event(),
            transform,
            scheme,
            my = {};

        my.onSelection = function(observer) {
            selectionChange.attach(observer);
        }

        my.el = function() {
            return $(element)[0];
        };

        my.resize = function(size) {
            var width = size.width - 25,
                height = size.height - 120;

            histogram.setSize(width, height);
            var pairs = model.get({pairs: true});
            var bins = histogram.bin(pairs, lastElement);
            histogram(element, bins);
        };

        my.setModel = function(newModel) {
            model = newModel;
            my.render();
        };

        my.setTransform = function(func, index) {
            model.setTransform(func, index);
            my.render();
        };

        my.getSelectedTransform = function() {
            return model.getTransform();
        }

        my.colors = function() {
            var pairs = model.get({pairs: true});
            var colors = _.map(_.map(pairs, _.last), histogram.color());

            return _.object(_.map(pairs, _.first), colors);
        }

        my.setColorScheme = function(colors) {
            histogram.setColors(colors);
            my.render();
        };

        my.select = function(selection) {
            model.get({selected: selection});
        }

        my.render = function() {
            var pairs = model.get({pairs: true});
            var bins = histogram.bin(pairs, lastElement);
            histogram(element, bins);
        }

        histogram.on("selected", function(selection) {
            // Publish the list of selected indices
            selectionChange.notify(_.map(selection, _.first));
        });

        return my;
    };

    var Dataset = synmap.Dataset = function(dataset) {
        var my = {},
            transform = undefined,
            transformIndex,
            keys = _.keys(dataset);

        my.layer = function(dataset) {
            return dataset.layer;
        };

        my.setTransform = function(func, index) {
            if (_.isFunction(func)) {
                transform = func;
                transformIndex = index;
            }
        };

        my.getTransform = function() {
            return transformIndex;
        };

        my.get = function(attr) {
            var values;

            if (transform === undefined) {
                values = _.values(dataset);
            } else {
                values = transform(_.values(dataset));
            }

            if (_.isObject(attr) && attr.pairs) {
                return _.zip(keys, values);
            }

            if (_.isNumber(attr)){
                return values[attr];
            }

            return values;
        };

        return my;
    };

    var Dropdown = synmap.Dropdown = function(element, datasets) {
        var my = {},
            listeners = [],
            selectedIndex = 0,
            el = $(element);

        el.on("change", function(e) {
            selectedIndex = el.find("option:selected").val();

            my.handleSelected(datasets[$(this).val()].data);
        });

        my.selected = function(callback) {
            listeners.push(callback);
        };

        my.handleSelected = function(selected) {
            listeners.forEach(function(callback) {
                callback(selected, selectedIndex);
            });
        };

        my.datasets = function() {
            return datasets;
        };

        my.select = function(index) {
            var child = el.children()[index];
            if (child) {
                $(child).attr("selected", "selected");
                el.change();
            }
        };

        var options = _.keys(datasets).map(function(dataset, index) {
            return $("<option>")
                .attr("value", index)
                .html(datasets[dataset].title);
        });

        el.append(options);

        return my;
    };

    var Filter = synmap.Filter = function(name, type) {
        var id = _.uniqueId();

        return {
            id: id,
            name:  name || ("Unnamed filter " + id),
            indices: [],
            enabled: true,
            type: type || ("inclusive"),
            toggle: function() {
                this.enabled = !this.enabled;
            },
            empty: function() {
                return !this.indices.length;
            },
            add: function(indice) {
                if (_.isArray(indice)) {
                    this.indices = _.union(indice);
                } else {
                    this.indices = _.union(arguments);
                }
            },
            remove: function() {
                this.indices = _.without(indices, arguments);
            },
            clear: function() {
                this.indices.length = 0;
            }
        };
    };

    var FilterChain = synmap.FilterChain = function() {
        var filters = [];

        return {
            indices: function() {
                var enabled,
                    grouped,
                    inclusive,
                    exclusive;

                // Only use filters that are enabled
                enabled = _.where(filters, {enabled: true});

                grouped = _.groupBy(enabled, "type");
                inclusive = _.flatten(_.pluck(grouped.inclusive, "indices"));
                exclusive = _.flatten(_.pluck(grouped.exclusive, "indices"));

                return _.difference(_.union(inclusive), _.union(exclusive));
            },
            addFilter: function(filter) {
                var added = _.where(filters, {id: filter.id});

                if (added.length === 0) {
                    filters.push(filter);
                }
            },
            removeFilter: function(id) {
                filters = _.without(filters, this.get(id));
            },
            getFilter: function(id) {
                return _.findWhere(filters, {"id": String(id)});
            },
            toggleFilter: function(id) {
                var filter = this.get(id);
                if (filter !== undefined) {
                    filter = filter.toggle();
                }
            },
            all: function() {
                return slice.call(filters, 0);
            }
        };
    };

    synmap.setup = function(options) {
        var root = options.rootElement;

        $.ajaxSetup({
            type: "GET",
            url: options.page,
            dataType: "html",
            cache: false
        });

        $(".dialog_box").dialog({
            autoOpen: false,
            width: 500,
        });

        $("#synmap_dialog").dialog({modal: true});

        if($('#org_name1').val() != "Search") {
            $('#org_name1').css({fontStyle: "normal"});
            timing('org_name1',1);
        }
        if($('#org_desc1').val() != "Search") {
            $('#org_desc1').css({fontStyle: "normal"});
            timing('org_desc1',1);
        }
        if($('#org_name2').val() != "Search") {
            $('#org_name2').css({fontStyle: "normal"});
            timing('org_name2',1);
        }
        if($('#org_desc2').val() != "Search") {
            $('#org_desc2').css({fontStyle: "normal"});
            timing('org_desc2',1);
        }

        if ($('#assemble')[0].checked) {
            $('#assemble_info').toggle();
        }

        merge_select_check();
        depth_algo_check();

        $("#pair_info").draggable();
        $("#tabs").tabs({selected:0});
        $(".resizable").resizable();
        $('#depth_org_1').html($('#org_id1 option:selected').html());
        $('#depth_org_2').html($('#org_id2 option:selected').html());

        if (options.autostart) {
            run_synmap(true, $('#regen_images')[0].checked);
        }

        $("#tabs").removeClass("invisible");

        // track analysis
        $("#synmap_go").on("click", function() {
            run_synmap(false, $('#regen_images')[0].checked);
            ga('send', 'event', 'synmap', 'run');
        });
    };

    //TODO: Remove basename generation from client
    function rand () {
        return ( Math.floor ( Math.random ( ) * 99999999 + 1 ) );
    }

    //TODO: Create a proper state object
    function populate_page_obj(basefile) {
        if (!basefile) {
            basefile = "SynMap_"+rand();
        }
        pageObj.basename = basefile;
        pageObj.nolog = 0;
        pageObj.waittime = 1000;
        pageObj.runtime = 0;
        pageObj.fetch_error = 0;
        pageObj.error = 0;
        pageObj.engine = $("<span></span>", {
            "class": "alert",
            text: "The job engine has failed."
        });
        pageObj.failed = "The workflow could not be submitted.";
    }

    //TODO: Replace helper with state checking
    function has_organisms () {
        return ($('#org_id1').val() != "") &&
            ($('#org_id2').val() != "");
    }

    function run_synmap(scheduled, regenerate){
        populate_page_obj();

        var org_name1 = pageObj.org_name1;
        var org_name2 = pageObj.org_name2;
        //feat_type 1 == CDS, 2 == genomic
        var feat_type1 = $('#feat_type1').val();
        var feat_type2 = $('#feat_type2').val();
        var org_length1 = pageObj.org_length1;
        var org_length2 = pageObj.org_length2;
        //seq_type == 1 is unmasked
        var seq_type1 = pageObj.seq_type1;
        var seq_type2 = pageObj.seq_type2;
        //check to see if we will allow this run
        var max_size = 50 * 1000 * 1000;
        // console.log (org_name1, org_length1, seq_type1);
        // console.log (org_name2, org_length2, seq_type2);
        if (( org_length1 > max_size && feat_type1 == 2 && seq_type1 == 1) &&
            ( org_length2 > max_size && feat_type2 == 2 && seq_type2 == 1) ) {
            var message = "You are trying to compare unmasked genomic sequences that are large!  This is a bad idea.  Chances are there will be many repeat sequences that will cause the entire pipeline to take a long time to complete.  This usually means that the analyses will use a lot of RAM and other resources.  As such, these jobs are usually killed before they can complete.  Please contact coge.genome@gmail.com for assistance with your analysis.";
                alert(message);
                $('#log_text').hide(0);
                    $('#results').show(0).html("<span class=alert>Analysis Blocked:</span>"+message);
                    return;
        }

        if (!has_organisms())
            return;

        if ($('#blast').val() == 5 && (feat_type1 != 1 || feat_type2 != 1) ) {
            alert('BlastP only works if both genomes have protein coding sequences (CDS) AND CDS is selected for both!');
            return;
        }

        var overlay = $("#overlay").show();

        //TODO Scale polling time linearly with long running jobs
        var duration = pageObj.waittime;

        $('#results').hide();

        if (regenerate) {
            return schedule(get_params("go", regenerate));
        }

        return $.ajax({
            type: 'GET',
            data: get_params("get_results"),
            dataType: "json",
            success: function(data) {
                if (!data.error) {
                    $('.box').css("float", "left");
                    overlay.hide();
                    $("#synmap_zoom_box").draggable();
                    $('#results').html(data.html).slideDown();
                } else {
                    schedule(get_params("go"))
                }
            },
            error: schedule.bind(this, get_params("go"))
        });
    }

    synmap.submit_assembly = function(e, input, gid1, gid2,flip) {
        e.preventDefault();

        var promise = $.ajax({
            dataType: "json",
            data: {
                fname: "generate_assembly",
                jquery_ajax: 1,
                input: input,
                gid1: gid1,
                gid2: gid2,
		flip: flip
            }
        });

        $("#dialog").dialog({
            autoOpen: true,
            position: {
                my: "top",
                at: "top+150",
            }
        });

        promise.then(wait_for_assembly)
               .then(download_file, report_error);
    };

    check_status = function(id) {
        return $.getJSON("jex/status/" + id);
    };

    wait_for_assembly = function(response) {
        var deferred = $.Deferred();

        if (response.success) {
            wait_for_job(response.id, deferred, response.output);
        } else {
            deferred.reject(undefined);
        }

        return deferred;
    };

    wait_for_job = function(id, promise, args) {
        check_status(id).then(function(response) {
            switch(response.status) {
                case "Completed":
                    promise.resolve(args);
                    break;
                case "Failed":
                    promise.reject("The workflow has failed");
                    break;
                default:
                    setTimeout(function() {
                        wait_for_job(id, promise, args);
                    }, 3000);
                    break;
            }
        });
    };


    report_error = function() {
        $("#dialog").html("The pseudo assembly could not be generated");
    };

    download_file = function(url) {
        $("#dialog").dialog("close");
        window.open(url, "_self");
    };

    function schedule(params) {
        pageObj.nolog=1;
        var status_dialog = $('#synmap_dialog');
        close_dialog(status_dialog);

        return $.ajax({
            type: "post",
            dataType: 'json',
            data: params,
            success: function(data) {
                $('.box').css("float", "left");
                $("#overlay").hide();
                status_dialog.dialog('open');

                if (data.error) {
                    status_dialog.find('#text').html(data.error).addClass("alert");
                    status_dialog.find('#progress').hide();
                    status_dialog.find('#dialog_error').slideDown();
                } else if (data.success) {
                    var link = "Return to this analysis: <a href="
                    + data.link + " onclick=window.open('tiny')"
                    + "target = _new>" + data.link + "</a>";

                    var logfile = '<a href="tmp/SynMap/'
                    + pageObj.basename + '.log">Logfile</a>';

                    $('#dialog_log').html(logfile);
                    $('#synmap_link').html(link);

                    update_dialog(data.request, "#synmap_dialog", synmap_formatter,
                            get_params("get_results"));
                } else {
                    status_dialog.find('#text').html(pageObj.failed).addClass("alert");
                    status_dialog.find('#progress').hide();
                    status_dialog.find('#dialog_error').slideDown();
                }
            },
            error: function(err) {
                $('.box').css("float", "left");
                $("#overlay").hide();
                status_dialog.dialog('open');
                status_dialog.find('#text').html(pageObj.failed).addClass("alert");
                status_dialog.find('#progress').hide();
                status_dialog.find('#dialog_error').slideDown();
            }
        })
    }

    function get_params(name, regenerate) {
        return {
            fname: name,
            tdd: $('#tdd').val(),
            D: $('#D').val(),
            A: $('#A').val(),
            beta: pageObj.beta,
            gm: $('#gm').val(),
            Dm: $('#Dm').val(),
            blast: $('#blast').val(),
            feat_type1: $('#feat_type1').val(),
            feat_type2: $('#feat_type2').val(),
            dsgid1: $('#dsgid1').val(),
            dsgid2: $('#dsgid2').val(),
            jobtitle: $('#jobtitle').val(),
            basename: pageObj.basename,
            email: $('#email').val(),
            regen_images: regenerate,
            width: $('#master_width').val(),
            dagchainer_type: $('#dagchainer_type').filter(':checked').val(),
            ks_type: $('#ks_type').val(),
            assemble: $('#assemble')[0].checked,
            axis_metric: $('#axis_metric').val(),
            axis_relationship: $('#axis_relationship').val(),
            min_chr_size: $('#min_chr_size').val(),
            spa_ref_genome: $('#spa_ref_genome').val(),
            show_non_syn: $('#show_non_syn')[0].checked,
            color_type: $('#color_type').val(),
            box_diags: $('#box_diags')[0].checked,
            merge_algo: $('#merge_algo').val(),
            depth_algo: $('#depth_algo').val(),
            depth_org_1_ratio: $('#depth_org_1_ratio').val(),
            depth_org_2_ratio: $('#depth_org_2_ratio').val(),
            depth_overlap: $('#depth_overlap').val(),
            fid1: pageObj.fid1,
            fid2: pageObj.fid2,
            show_non_syn_dots: $('#show_non_syn_dots')[0].checked,
            flip: $('#flip')[0].checked,
            clabel: $('#clabel')[0].checked,
            skip_rand: $('#skiprand')[0].checked,
            color_scheme: $('#color_scheme').val(),
            chr_sort_order: $('#chr_sort_order').val(),
            codeml_min: $('#codeml_min').val(),
            codeml_max: $('#codeml_max').val(),
            logks: $('#logks')[0].checked,
            csco: $('#csco').val(),
            jquery_ajax: 1,
        };
    }

    function synmap_formatter(item) {
        var msg;
        var row = $('<li>'+ item.description + ' </li>');
        row.addClass('small');

        var job_status = $('<span></span>');

        if (item.status == 'scheduled') {
            job_status.append(item.status);
            job_status.addClass('down');
            job_status.addClass('bold');
        } else if (item.status == 'completed') {
            job_status.append(item.status);
            job_status.addClass('completed');
            job_status.addClass('bold');
        } else if (item.status == 'running') {
            job_status.append(item.status);
            job_status.addClass('running');
            job_status.addClass('bold');
        } else if (item.status == 'skipped') {
            job_status.append("already generated");
            job_status.addClass('skipped');
            job_status.addClass('bold');
        } else if (item.status == 'cancelled') {
            job_status.append(item.status);
            job_status.addClass('alert');
            job_status.addClass('bold');
        } else if (item.status == 'failed') {
            job_status.append(item.status);
            job_status.addClass('alert');
            job_status.addClass('bold');
        } else {
            return;
        }

        row.append(job_status);

        if (item.elapsed)  {
            row.append(" in " + coge.utils.toPrettyDuration(item.elapsed));
        }

        /*
        if (item.status == "skipped") {
            row.append("<p>The analyses previously was generated</p>");
        }
        */

        return row;
    }

    function update_dialog(request, identifier, formatter, args) {
        var get_status = function () {
            $.ajax({
                type: 'GET',
                url: request,
                dataType: 'json',
                success: update_callback,
                error: update_callback,
            });
        };

        var get_poll_rate = function() {
            pageObj.runtime += 1;

            if (pageObj.runtime <= 5) {
                return 1000;
            } else if (pageObj.runtime <= 60) {
                return 2000;
            } else if (pageObj.runtime <= 300) {
                return 5000;
            } else if (pageObj.runtime <= 1800) {
                return 30000;
            } else if (pageObj.runtime <= 10800) {
                return 60000;
            } else {
                return 300000;
            }
        };

        var fetch_results = function(completed) {
            var request = window.location.href.split('?')[0];
            args.fname = 'get_results';
            dialog = $(identifier);

            $.ajax({
                type: 'GET',
                url: request,
                data: args,
                dataType: "json",
                success: function(data) {
                    if (completed && !data.error) {
                        $("#synmap_zoom_box").draggable();
                        $('#results').html(data.html);
                        dialog.find('#progress').hide();
                        dialog.find('#dialog_success').slideDown();
                    } else {
                        $('#results').html(data.error);
                        dialog.find('#progress').hide();
                        dialog.find('#dialog_error').slideDown();
                    }
                },
                error: function(data) {
                    if (pageObj.fetch_error >= 3) {
                        dialog.find('#progress').hide();
                        dialog.find('#dialog_error').slideDown();
                    } else {
                        pageObj.fetch_error += 1;
                        console.log("error");
                        var callback = function() {fetch_results(completed)};
                        setTimeout(callback, 100);
                    }
                }
            });
        }

        var update_callback = function(json) {
            var dialog = $(identifier);
            var workflow_status = $("<p></p>");
            var data = $("<ul></ul>");
            var results = [];
            var current_status;
            var timeout = get_poll_rate();

            var callback = function() {
                update_dialog(request, identifier, formatter, args);
            }

            if (json.error) {
                pageObj.error++;
                if (pageObj.error > 3) {
                    workflow_status.html(pageObj.engine);
                    dialog.find('#text').html(workflow_status);
                    dialog.find('#progress').hide();
                    dialog.find('#dialog_error').slideDown();
                    return;
                }
            } else {
                pageObj.error = 0;
            }

            if (json.status) {
                current_status = json.status.toLowerCase();
                workflow_status.html("Workflow status: ");
                workflow_status.append($('<span></span>').html(json.status));
                workflow_status.addClass('bold');
            } else {
                setTimeout(callback, timeout);
                return;
            }

            if (json.jobs) {
                var jobs = json.jobs;
                for (var index = 0; index < jobs.length; index++) {
                    var item = formatter(jobs[index]);
                    if (item) {
                        results.push(item);
                    }
                }
            }

            if (!dialog.dialog('isOpen')) {
                return;
            }

            //FIXME Update when a workflow supports elapsed time
            if (current_status == "completed") {
                var total = json.jobs.reduce(function(a, b) {
                    if (!b.elapsed) return a;

                    return a + b.elapsed;
                }, 0);

                var duration = coge.utils.toPrettyDuration(total);

                workflow_status.append("<br>Finished in " + duration);
                workflow_status.find('span').addClass('completed');
                fetch_results(true);
            } else if (current_status == "failed" || current_status == "error"
                    || current_status == "terminated"
                    || current_status == "cancelled") {
                workflow_status.find('span').addClass('alert');
                fetch_results(false);
            } else if (current_status == "notfound") {
                setTimeout(callback, timeout);
                return;
            } else {
                workflow_status.find('span') .addClass('running');
                setTimeout(callback, timeout);
            }

            results.push(workflow_status);
            data.append(results);
            dialog.find('#text').html(data);
        };

        get_status();
    }

}(this.coge || (this.coge = {}), jQuery, _));
