var spinner = '<img src="picts/ajax-loader.gif"/>';

function query_results(options) {
    options.fname = "get_results";

    return $.ajax({
        type: "GET",
        data: options,
        dataType: "json"
    });
}

function launch(dialog, results, options) {
    var _results = $(results);
    var status_dialog = $(dialog).dialog({
        title: 'Running SynFind ....',
        modal: true,
        autoOpen: true
    });

    options.fname = 'go_synfind';
    // Sends to JEX and blocks until completion
    var submit = $.ajax({
        type: 'POST',
        dataType: 'json',
        data: options,
        success : function(response) {
            status_dialog.unbind().on("dialogclose", function() {
                _results.removeClass('hidden').slideDown('slow');

                // reset dialog
                status_dialog.find(".dialog-error,.dialog-complete").hide();
                status_dialog.find(".dialog-running").show();
                status_dialog.find(".dialog-text").html("");
                status_dialog.find(".dialog-link").html("");
                status_dialog.find(".dialog-log").html("");
            });

            if(response.success) {
                pageObj.runtime = 0;
                pageObj.fetch_error = 0;
                pageObj.error = 0;
                pageObj.engine = "<span class=\"alert\">The job engine has failed.</span><br>Please use the link below to use the previous version of SynMap.";

                var link = $("<a></a>")
                    .attr("href", response.link)
                    .html(response.link);

                var link_message = $("<span></span>")
                    .html("Return to this analysis: ")
                    .append(link);

                var logfile = $("<a></a>")
                    .attr("href", response.logfile)
                    .html("Logfile");

                status_dialog.find(".dialog-link").html(link_message);
                status_dialog.find(".dialog-log").html(logfile);

                options.fname = "get_results";
                options.logfile = response.logfile;

                update_dialog(response.request, dialog, results, formatter, options);
            } else {

                var error = $("<div></div>")
                    .addClass("alert")
                    .html(response.message);

                _results.append(results, error);
                status_dialog.find(".dialog-error").slideDown();
                status_dialog.find(".dialog-running").hide();
            }
        },
    });
}

function formatter(item) {
    var msg;
    var row = $('<li>'+ item.description + ' </li>');

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

    return row;
}

function update_dialog(request, identifier, result, formatter, args) {
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

    var fetch_results = function(completed, attempts) {
        dialog = $(identifier);

        $.ajax({
            type: 'GET',
            url: 'SynFind.pl',
            data: args,
            dataType: "json",
            success: function(data) {
                if (completed && data.success) {
                    handle_results(result, data.html);
                    dialog.find('.dialog-running').hide();
                    dialog.find('.dialog-complete').slideDown();
                } else {
                    handle_results(result, data.message);
                    dialog.find('.dialog-running').hide();
                    dialog.find('.dialog-error').slideDown();
                }
            },
            error: function(data) {
                if (attempts >= 3) {
                    dialog.find('.dialog-running').hide();
                    dialog.find('.dialog-error').slideDown();
                } else {
                    var callback = function() {fetch_results(completed, attempts + 1)};
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
            update_dialog(request, identifier, result, formatter, args);
        }

        if (json.error) {
            pageObj.error++;
            if (pageObj.error > 3) {
                workflow_status.html(pageObj.engine);
                dialog.find('.dialog-text').html(workflow_status);
                dialog.find('.dialog-running').hide();
                dialog.find('.dialog-error').slideDown();
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
            fetch_results(true, 1);
        } else if (current_status == "failed" || current_status == "error"
                || current_status == "terminated"
                || current_status == "cancelled") {
            workflow_status.find('span').addClass('alert');
            fetch_results(false, 1);
        } else if (current_status == "notfound") {
            setTimeout(callback, timeout);
            return;
        } else {
            workflow_status.find('span') .addClass('running');
            setTimeout(callback, timeout);
        }

        results.push(workflow_status);
        data.append(results);
        dialog.find('.dialog-text').html(data);
    };

    get_status();
}

$.fn.getLength = function(val){
  var opt_length;
  var searchable;
  var blanked=0; //otherwise get math problems later...boo javascript
  this.each(
    function()
    {
        var opts = this.options;
        opt_length = opts.length;
        if (opt_length == 0) {return opt_length;}
        searchable = opts[0].id;
        if (searchable == 'blank') {blanked++;} //Counts the number of instances of blank
        if (val){
          for(var i=1;i<opts.length;i++)
          {
            searchable += ","+opts[i].id;
            //need to chopoff last comma
          }
        }
    }
  );
  if(val) return searchable;
  if (blanked) {opt_length-=blanked;} //subtract elements that are classified as blank above
  return opt_length;
};

$.fn.sortSelect = function(){
  this.each(
      function()
      {
        if(this.nodeName.toLowerCase() != "select"){
          return;}
        var opts = this.options;
        var sortArray = [];
        for(var i=0;i<opts.length;i++)
        {
            sortArray[i] = {
                    v: opts[i].value,
                    t: opts[i].text,
                    d: opts[i].id,
                    }
        }
        sortArray.sort(
                function(obj1,obj2)
                {
                    obj1t = obj1.t.toLowerCase(),
                    obj2t = obj2.t.toLowerCase();
                    if(obj1t == obj2t){
                      return 0;}
                    return obj1t < obj2t ? -1 : 1;
                }
        );
        for(var i=0;i<opts.length;i++)
        {
            opts[i].id = sortArray[i].d;
            opts[i].text = sortArray[i].t;
            opts[i].value = sortArray[i].v;
        }
      }
    );
    return this;
};

function handle_results(selector, html) {
    $(selector).html(html);
    init_table_sorter();
    setup_button_states();
}

/*------------------------------------------------------------------------------
    Feature Search
------------------------------------------------------------------------------*/

function search_chain(val) {
    $('#data_table').show(0);
    $('#anno,#accn_list,#Source,#FeatType').hide();

    var accnminlen = 3;
    var annominlen = 8;

    if ( val == 1
        || $('#accn').val().length > accnminlen
        || $('#annosearch').val().length > annominlen )
    {
        $('#anno').html(spinner).show();
        go_cogefeatsearch();
    }
    else {
        $('anno').html('Search not run.').show();
    }
}

function go_cogefeatsearch() {
    $('#fid').html('');
    var accn = $('#accn').val();
    var annosearch = $('#annosearch').val();
    var org_id = $('#org_id_feat').val();
    var org_name = $('#org_name_feat').val();
    var org_desc = $('#org_desc_feat').val();

    //cogefeatsearch(['args__accn','accn', 'args__anno','annosearch',
    //'args__org_id','org_id_feat', 'args__org_name','org_name_feat',
    //'args__org_desc','org_desc_feat'],[source_search_chain]);
    $.ajax({
        data: {
            fname:      'cogefeatsearch',
            accn:       accn,
            anno:       annosearch,
            org_id:     org_id,
            org_name:   org_name,
            org_desc:   org_desc
        },
        success : function(data) {
            source_search_chain(data);
        },
    });
}

function source_search_chain(val) {
    if (val) {
        $('#accn_list').html(val).show();
    }
    else {
        $('#fid').html('');
    }
    $('#Source').html(spinner);
    $('#FeatType').html(spinner);
    $('#anno').html('');

    var accn = $('#accn_select').val()[0];
    var org_id = $('#org_id_feat').val();
    var org_name = $('#org_name_feat').val();
    var org_desc = $('#org_desc_feat').val();

    //source_search(['args__accn','accn_select', 'args__org_id','org_id_feat',
    //'args__org_name','org_name_feat','args__org_desc','org_desc_feat'], [get_types_chain]);
    $.ajax({
        data: {
            fname:      'source_search',
            accn:       accn,
            org_id:     org_id,
            org_name:   org_name,
            org_desc:   org_desc
        },
        success : function(data) {
            get_types_chain(data);
        },
    });
}

function get_types_chain(val) {
    if (val) {
        $('#Source').html(val).show();
    }

    var accn = $('#accn_select').val()[0];
    var dsgid = $('#feat_dsgid').val()[0];

    //get_types(['args__accn','accn_select', 'args__dsgid','feat_dsgid'],[get_anno_chain]);
    $.ajax({
        data: {
            fname:      'get_types',
            accn:       accn,
            dsgid:      dsgid,
        },
        success : function(data) {
            var obj = jQuery.parseJSON(data);
            if (obj) {
                get_anno_chain(obj.html, obj.dsgid);
            }
        },
    });
}

function get_anno_chain(val, dsgid, fid) {
    if (val) {
        $('#FeatType').html(val).show();
        if ($('#add_all').is(":hidden")) {$('#add_all').show(0);}
        if ($('#remove').is(":hidden")) {$('#remove').show(0);}
        if ($('#clear').is(":hidden")) {$('#clear').show(0);}
        if ($('#send').is(":hidden")) {$('#send').show(0);}
    }

    var accn;
    if ( $('#accn_select').length )
        accn = $('#accn_select').val()[0];
    var type;
    if ( $('#type_name').length ) {
        var type_name = $('#type_name').val();
        if (type_name)
            type = type_name[0];
        else {
            alert('The selected genome doesn\'t have any CDS features. Please select a different genome.');
            return;
		}
    }

    $('#anno').html(spinner).show();

    //get_anno(['args__accn','accn_select', 'args__type','type_name', 'args__dsgid','args__'+dsgid],[show_anno]);
    $.ajax({
        data: {
            fname:  'get_anno',
            fid:    fid,
            accn:   accn,
            type:   type,
            dsgid:  dsgid
        },
        success : function(data) {
			if (data) {
            	var obj = jQuery.parseJSON(data);
            	if (obj)
                	show_anno(obj.anno, obj.fid);
            } else
				$('#anno').html('');
        },
    });
}

function show_anno(anno, fid) {
    $('#anno').html(anno);
    $('#fid').html(fid);
    setup_button_states();
}

function search_org_feat(val){
    if (val == 'name') {
        $('#org_desc_feat').val("");
        var searchterm = $('#org_name_feat').val();
        get_orgs_feat('name', searchterm); //pageObj.time = setTimeout("get_orgs_feat(['args__type','args__name','args__search','org_name_feat'],['org_list_feat'])",500);
    }
    else if (val == 'desc') {
        $('#org_name_feat').val("");
        var searchterm = $('#org_desc_feat').val();
        get_orgs_feat('desc', searchterm); //pageObj.time = setTimeout("get_orgs_feat(['args__type','args__desc','args__search','org_desc_feat'],['org_list_feat'])",500);
    }
}

function get_orgs_feat(type, searchterm) {
    $.ajax({
        data: {
            fname:  'get_orgs_feat',
            type:   type,
            search: searchterm,
        },
        success : function(html) {
            $('#org_list_feat').html(html);
        },
    });
}

function onEnter(e){
    if (!e) {
        var e = window.event;
    }
    if (e.which == 13){
        search_chain(1);
    }
}

function update_info_box(featid) {
    generate_feat_info(featid); //generate_feat_info(['args__'+featid],['feature_info_popup']);
    $('#feature_info_popup').dialog('open');
}

function generate_feat_info(featid) {
    $.ajax({
        data: {
            fname:  'generate_feat_info',
            featid: featid
        },
        success : function(html) {
            $('#feature_info_popup').html(html);
        },
    });
}

function init_table_sorter() {
    $(function() {
        $("#syntelog_table").tablesorter({
            sortClassAsc: 'headerSortUp',       // Class name for ascending sorting action to header
            sortClassDesc: 'headerSortDown',    // Class name for descending sorting action to header
            headerClass: 'header',              // Class name for headers (th's)
            widgets: ['zebra'],
            headers: {
                5: { sorter: 'digit' }
            },
            textExtraction: 'complex',
        });
    });
}

function get_master(link) {
    pad = $('#pad').val();
    pad = pad*1;
    if (pad) {
        link += ";pad="+pad;
    }
    window.open(link);
}
