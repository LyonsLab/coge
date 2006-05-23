/* written by: Brent Pedersen */
/*        BSD license         */
/*  ----------------------    */


// ahah (only simpler and target can be a function)
// jfetch('script.pl?par=val',myFunc/*,obj*/);
// jfetch('test.html','divId');
function jfetch(url,target,context) {
  var req = (XMLHttpRequest)?new XMLHttpRequest():new ActiveXObject("Microsoft.XMLHTTP");
  req.open("GET",url,true);
  req.onreadystatechange = function() {
    if(req.readyState == 4){
      var rsp = req.responseText;
      if(target.constructor == Function){
        target.call(context,rsp);
      }else{
        var t = document.getElementById(target);
        if(t.value==undefined){t.innerHTML=rsp;}else{t.value=rsp;}
      }
    }
  };
  req.send(null);
}

// like toggle, but simpler yo.
function yotoggle(id) {
    var es = document.getElementById(id).style;
    es.display = (es.display != 'none')?'none':'';
}

// simplified version of many out there
function getElementsByClassName(c,t) {
    var els = document.getElementsByTagName(t || '*');
    var l = els.length;
    c = new RegExp('\\b' + c + '\\b');
    var ans = [];
    while(l--) {
        if ( c.test(els[l].className) ) {
          ans.push(els[l]);
        }
    }
    return ans;
}

/* addEvent, handleEvent, fixEvent, removeEvent        */
/* written by D. Edwards and J. Resig                  */
/* http://dean.edwards.name/weblog/2005/10/add-event/  */
/* modified to include scope override, custom events   */
/* USEAGE: addEvent(elm,'click',myFunc)                */
/* --------------------------------------------------- */
function addEvent(el, type, handler) {

    // bind all events on a particular html el to an object
    // by setting el.scope = someObj.
    // to bind only 1 function, send in a 4th arg to this function
    var s = el.scope || arguments[3] || false;
    if(s){ handler = bindEvent(s,handler); }

    // TODO: test custom events
    el = (el!=null) ? el : customEvent;

    if (!handler.$$guid) handler.$$guid = addEvent.guid++;
    if (!el.events) el.events = {};

    var handlers = el.events[type];
    if (!handlers) {
        handlers = el.events[type] = {};
        if (el["on" + type]) { handlers[0] = el["on" + type]; }
     }
     handlers[handler.$$guid] = handler;
     el["on" + type] = handleEvent;
}
addEvent.guid = 1;

var customEvent = { type : 'custom' };

function removeEvent(element, type, handler) {
    if (element.events && element.events[type]) {
        delete element.events[type][handler.$$guid];
    }
};

// modified from jquery.js by John Resig 
// new 3rd arg allows custom events (!!)
// to receive the last event;
function fireEvent(el,type,evt) {
    if ( !el["on" + type] ){ return ; }
    evt = (evt) ? evt : {type:type}; 
    evt.customtype = (type) ? type : evt.type;
    el["on" + type](evt);
}

// modified
function handleEvent(event) {
  var returnValue = true;
  event = extendEvent(event || fixEvent(window.event));

  var handlers = this.events[event.type];
  for (var i in handlers) {
    this.$$handleEvent = handlers[i];
    // make sure we have the element even when overriding scope
    event.element = this;
    if (this.$$handleEvent(event) === false) {
      returnValue = false;
    }
  }
  return returnValue;
};

function fixEvent(event) {
    event.preventDefault =  $prevent;
    event.stopPropagation = $cancel; 
    return event;
};
$prevent = function() { this.returnValue = false; };
$cancel = function() { this.cancelBubble = true; };
$stop = function() { this.preventDefault(); this.stopPropagation(); }
/* End Edwards and Resig Stuff */


/* Extended Event Stuff   */ 
// Brent Pedersen BSD license
/*------------------------*/

// add cross-browserness for where the event occured
// and the location of the element in which the event
// occured. 
function extendEvent(event){
    // cross browser element;

    // not needed:
   // fixed above: event.element = event.target || event.srcElement || 'custom';

    event.stop = $stop;
    // give the location of the element that was clicked
    event.getTargetXY = getElementXY;
    // where in page did evt occur;
    event.getPageXY = _getPageXY;
    // where in the element did evt occur?
    event.getItemXY = _getItemXY;
    return event;
} 

// added as method to event to get location of
// a click within the html element fo the click
function _getItemXY(){
   var elXY = this.getTargetXY();
   var evXY = this.getPageXY();
   return [evXY[0] - elXY[0], evXY[1] - elXY[1]];
}

// bind a function so that the object o can be 
// access as 'this' whenever that function is fired;
function bindEvent(o,f){
    return function(){
        f.apply(o,arguments);
    }
    // f=o=null;
}

// modified from prototype.js (MIT license) sam stephenson
//  called by the extendEvent object
function _getPageXY() {
    return [ 
    this.pageX || (this.clientX + (document.documentElement.scrollLeft || document.body.scrollLeft))
    , 
    this.pageY || (this.clientY + (document.documentElement.scrollTop || document.body.scrollTop)) 
    ];
}
/* end prototype  */

// find the position of an object either from 
// an event or as an obj returns an array
// should add a cache...
function getElementXY(o) {
    o = (o != undefined )?o:this.element;
    var loc = [0,0];
    if (o.offsetParent) {
        while (o.offsetParent){
            loc[0] += o.offsetLeft;
            loc[1] += o.offsetTop;
            o = o.offsetParent;
        }
    }
    else if (o.y || o.x){
        loc[0] += o.x;
        loc[1] += o.y;
    }
    return loc;
}

function falseFunction(){
  return false;
}
