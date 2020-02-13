var Gobe = {
    WEDGE: 'wedge',
    LINE: 'line',
    DIV: 'flashcontent',
/*    'get_movie': function () {
          if (swfobject.ua.ie) {
              return window[Gobe.DIV];
          } else {
              return document[Gobe.DIV];
          }
    },*/
    'handle_html': function(html){
        document.getElementById('querybox').innerHTML = html;
    },
    'set_genespace': function(updown, idx, value){
        var eid = 'dr' + updown + idx;
        document.getElementById(eid).value = value;
    },
    'clear': function(){
        // dont change this!
//        Gobe.swf.clear_wedges();
        $('svg').find('rect').remove();
        $('svg').find('line').remove();
        $('svg').find('polygon').remove();
        $('#svg')[0].connectors = {};
    },
    'set_linewidth': function(linewidth){
        // dont change this. could check for 0 <= lw <= ?
//        Gobe.swf.set_linewidth(linewidth);
    },
    'set_connector': function(contype){
        // dont change this!
        contype = contype.toLowerCase();
        if (! (contype == Gobe.WEDGE || contype == Gobe.LINE)){
            alert('must be one of: ' + Gobe.WEDGE + ' or ' + Gobe.LINE);
        }
//        Gobe.swf.set_connector(contype);
    }

};
