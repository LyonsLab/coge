var spinner = '<img src="picts/ajax-loader.gif"/>';

//function loading(id, msg) {
//  $('#'+id).html('<font class="loading">Loading '+msg+' . . .</font>');
//}

//function counting(){
//  var count = $('#genome_choice').getLength();
//  if (count == 0) {
//      $('#genome_choice').html('<option id=blank value=null>No Organism Selected</option>');
//  }
//  $('#selected_genome_count').html(count);
//}

//function update_basename(basename){
//    pageObj.basename = basename;
//}
//
//function reset_basename(){
//    if (pageObj.basename) pageObj.basename = 0;
//}
//
//function generate_basefile() {
//  $.ajax({
//      data: {
//          fname: 'generate_basefile',
//      },
//      success : function(filename) {
//          update_basename(filename);
//      }
//  });
//}

function launch(dialog, results) {
    var fid = $('#fid').html();
    if (!fid) {
        alert('Please search and select a feature for this analysis.')
        return;
    }

//  if (!pageObj.basename) {
//      setTimeout("launch()", 100);
//      return;
//  }

    var check = $('#genome_choice').getLength();
    if( $('#blank').val() || check == 0 ){
        alert('Please select at least one genome to search.');
        return;
    }

    var _results = $(results).hide().html('');
    var status_dialog = $(dialog).dialog({
        title: 'Running SynFind ....',
        modal: true,
        autoOpen: true
    });

//  pageObj.nolog = 1;
    pageObj.waittime = 1000;
//  monitor_log();
//

    var selected_genomes = $('#genome_choice').getLength(1);
    var qdsgid = $('#feat_dsgid').val();
    var ws = $('#ws').val();
    var co = $('#co').val();
    var sf = $('#sf').val();
    var algo = $('#algo').val();
    var sd = $('#sd').val();

    //go(['args__dsgids','args__'+selected_genomes,'args__fid','fid',
    //'args__qdsgid','feat_dsgid', 'args__basename','args__'+pageObj.basename,
    //'args__window_size','ws', 'args__cutoff','co', 'args__scoring_function','sf',
    //'args__algo','algo', 'args__depth','sd'],[handle_results],'POST');

    // Sends to JEX and blocks until completion
    $.ajax({
        type: 'POST',
        dataType: 'html',
        data: {
            fname:              'go_synfind',
            dsgids:             selected_genomes,
            fid:                fid,
            qdsgid:             qdsgid,
            window_size:        ws,
            cutoff:             co,
            scoring_function:   sf,
            algo:               algo,
            depth:              sd
        },
        success : function(html) {
            if(!html) {
                status_dialog.find(".error").slideDown();
                status_dialog.find(".running").hide();
                return;
            }

            status_dialog.find(".running").hide();
            status_dialog.find(".complete").slideDown();

            _results.html(html);

            status_dialog.unbind().on("dialogclose", function() {
                _results.removeClass('hidden').slideDown();

                // reset dialog
                status_dialog.find(".error,.complete").hide();
                status_dialog.find(".running").show();
            });

            init_table_sorter();
            setup_button_states();
        },
    });
}

//function show_add() {
//  if($('#add').is(":hidden")) {
//    $('#remove').hide(0);
//    $('#add').show(0);
//  }
//}

//function org_search(desc_search){
//  if (pageObj.time){
//   clearTimeout(pageObj.time);
//  }
//  name = $('#org_name').val();
//  desc = $('#org_desc').val()
//  if (desc_search) {
//    if ( pageObj.prev_desc == desc ) {return;}
//    pageObj.time = setTimeout("get_orgs(['args__desc','args__'+desc], ['org_list']);",500);
//    setTimeout("pageObj.prev_desc=desc",500);
//    pageObj.prev_name="";
//   }
//  else
//   {
//    if ( pageObj.prev_name == name) {return;}
//    pageObj.time = setTimeout("get_orgs(['args__name','args__'+name], ['org_list']);",500);
//    setTimeout("pageObj.prev_name=name",500);
//    pageObj.prev_desc="";
//
//   }
//  setTimeout("seq_type_search()",600);
//}

//function seq_type_search() {
// if (ajax.length)
//    {
//     setTimeout("seq_type_search()",500);
//    }
//  else if ($('#org_id').val()) {
//    gen_dsg_menu(['args__oid', 'org_id'],['org_seq_types']);
//    show_add();
//   }
//  else {
//   $('#org_seq_types').html('');
//   $('#add').hide(0);
//  }
//}

//function monitor_log(log) {
//  var fasta = "";
//  var blast = 0;
//  var results = 0;
//  var done = 0;
//  var match;
//  pageObj.finished = 0;
//
//  if (log) {
//      if (log.match(/\*\s+fasta\sfile/i)) {
//          fasta="Generating blastable databases . . . <br/>";
//          match=log.match(/\*(.+)\*\sblastdb/g);
//          if (match) {
//              var crappy_workaround;
//              for (i=0;i<match.length;i++) {
//                  crappy_workaround = match[i].match(/\*(.+)\*/);
//                  fasta += "&nbsp;&nbsp;&nbsp;"+crappy_workaround[1]+" database built!<br/>";
//              }
//          }
//      }
//
//      if (log.match(/running/i)) {
//          fasta += "Blastable databases generation done!<br/>";
//
//          var blast_program = log.match(/-p\s(\w+)\s/);
//          if (!blast_program)
//              blast_program = log.match(/(blastz)/);
//          blast="Running blast algorithm . . . <br/>";
//
//          match=log.match(/\*(.+)\*\sblast\sanalysis\scomplete/g);
//          if (match) {
//              var crappy_workaround;
//              for (i=0;i<match.length;i++) {
//                  crappy_workaround = match[i].match(/\*(.+)\*/);
//                  blast += "&nbsp;&nbsp;&nbsp;"+crappy_workaround[1]+" complete!<br/>";
//              }
//          }
//      }
//
//      if (log.match(/Results/i)) {
//          blast += "Analysis complete!<br/>";
//          results = "Collating results. . . ";
//      }
//      if (log.match(/Finished/i)) {
//          results += "done!<br/>";
//          done = "Displaying results.  This may take your browser some time to process the HSP table.  Please be patient."
//          pageObj.finished = 1;
//      }
//  }
//  else {
//      pageObj.nolog += 1;
//  }
//
//  var message = "Initializing search . . . ";
//  if (fasta) message += "done!<br/>"+fasta;
//  if (blast) message += blast;
//  if (results) message = results;
//  if (done) message += done;
//
//  if ( !pageObj.finished && pageObj.nolog < 20 ) {
//      pageObj.waittime = pageObj.waittime*2;
//      if (pageObj.waittime > 60*1000)
//          pageObj.waittime = 60*1000;
//      message += "<br/>Next progress check in "+pageObj.waittime/1000+" seconds.";
//      setTimeout("monitor_log()", pageObj.waittime);
//  }
//  if (message)
//      $('#log_text').html("<div class='small'>"+message+'</div>');
//}

//function clear_org_list()
// {
// $('#genome_choice').empty();
// counting();
// }

//function remove_selected_orgs() {
// $('#genome_choice option:selected').each(function(){
//    $('#'+$(this).val()).remove();
//    });
// counting();
// }

//function add_all_orgs() {
// var ids;
// var count =0;
// $('#org_id option').each(function(){
//     ids = ids+","+$(this).val();
//     if (count == 100)
//      {
//         get_dsg_for_search_menu(['args__orgid','args__'+ids],[add_to_list])
//         ids = "";
//         count =0;
//      }
//
//     count++;
//   });
// get_dsg_for_search_menu(['args__orgid','args__'+ids],[add_to_list])
// sort_genome_choice();
//}

//function add_selected_orgs() {
// var ids;
// var count =0;
// $('#org_id option:selected').each(function(){
//     ids = ids+","+$(this).val();
//     if (count == 100)
//      {
//         get_dsg_for_search_menu(['args__orgid','args__'+ids],[add_to_list]);
//         ids = "";
//         count =0;
//      }
//
//     count++;
//    });
//  if (count == 1)
//   {
//        if ($('#dsgid').length > 0) { // mdb tempfix for issue #21
//            get_dsg_for_search_menu(['args__dsgid','dsgid'],[add_to_list]);
//        }
//   }
//  else
//   {
//    get_dsg_for_search_menu(['args__orgid','args__'+ids],[add_to_list]);
//   }
//  sort_genome_choice();
//}

//function add_to_list(stuff){
//  var orgs = stuff.split(':::');
//  for (var i=0; i < orgs.length; i++)
//   {
//     var item = orgs[i].split('::');
//     id = item[0];
//     org = item[1];
//     if (!id && !org) {continue;}
//     var check = $('#'+id).val();
//     if (check){ continue; }
//     var html = '<option id='+id+' value='+id+' ondblclick="remove_selected_orgs();">'+org+'</option>';
//     $('#blank').remove();
//     $('#genome_choice').append(html);
//   }
//  counting();
//}

//function sort_genome_choice() {
//    if (ajax.length)
//     {
//       setTimeout("sort_genome_choice()",100);
//       return;
//     }
//   $('#genome_choice').append().sortSelect();
//}

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

function handle_results(html) {
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
    if ( $('#type_name').length )
        type = $('#type_name').val()[0];

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
            var obj = jQuery.parseJSON(data);
            if (obj) {
                show_anno(obj.anno, obj.fid);
            }
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
            //textExtraction: 'complex',
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
