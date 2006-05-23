/* adapted by: Brent Pedersen       */
/*        BSD license               */
/*  ----------------------          */
/* to take advantage of the event distribution   */
/* in addEvent by dean edwards and john resig    */
/* and the event normalization and leakage       */
/* prevention in prototype.js by sam stephenson  */

// addEvent:
/* http://dean.edwards.name/weblog/2005/10/add-event/  */
// jquery:
/* http://jquery.com/                                  */
// prototype:
/* http://prototype.conio.net/                         */

/******************************************************************
   
   addEvent(element,'click',myFunc);
   
then the event object is 'normalized' before being sent to myFunc:
   
   function myFunc(e){
      // no more e || window.event;
      // where in the page did the event occur.
      e.getPageXY()

      // cross-browser stoppage
      e.stop()

      // the name of the html element that is the source/target of
      // the event
      e.element.name
   }

the 4th argument in subscribe is the object to bind the fn to.
in that case, the element that triggered the event is available
as this.element:
   
   myObj = { a: 'hello', b:'goodbye' };
   Event.subscribe(el2,'blur', testFn, myObj);

   function testFn(e){
       assert(this.a == 'hello')
       assert(this.element==el2);
    }

finally, custom events are simple:
      
     addEvent(element,'drag',onDrag);

     // fire an event by name...
     fireEvent(element,'drag')

     function onDrag(e){
       // do some stuff
     }
   
******************************************************************/
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
    el["on" + type]({type:type});
}

// modified
function handleEvent(event) {
  var returnValue = true;
  event = extendEvent(event || window.event);

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

/* End Edwards and Resig Stuff */



function extendEvent(event){
    for (var k in eventExtenders){
        event[k] = eventExtenders[k];
    }
    return event;
} 

// bind a function so that the object o can be 
// access as 'this' whenever that function is fired;
function bindEvent(o,f){
    return function(){
        f.apply(o,arguments);
    }
}

 
// modified from prototype.js (MIT license) sam stephenson
//  called by the extendEvent object
var eventExtenders = {

  isLeftClick: function() {
    return (((event.which) && (event.which == 1)) ||
            ((event.button) && (event.button == 1)));
  },
  // get where in the element it occurred;
  getPageXY: function() {
     return [this.pageX || (this.clientX + (document.documentElement.scrollLeft || document.body.scrollLeft)) , this.pageY    || (this.clientY + (document.documentElement.scrollTop || document.body.scrollTop))];
  },
  getItemXY: function(){
     var elXY = this.getElementXY();
     var evXY = this.getPageXY();
     return [evXY[0] - elXY[0], evXY[1] - elXY[1]];
  },
  stop: function() {
    if (this.preventDefault) { 
      this.preventDefault(); 
      this.stopPropagation(); 
    } else {
      this.returnValue = false;
      this.cancelBubble = true;
    }
  },
  getElementXY : getElementXY
};


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
