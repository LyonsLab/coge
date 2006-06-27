falsefunc = function(){return 0};
falseevent = function(e){e.stop()};

function jfetch(url,target,context) {
  var req = _req_();
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
var _req_ = 
    (typeof XMLHttpRequest != 'undefined')
   ? function(){ return new XMLHttpRequest()}
   : function(){ return new ActiveXObject("Microsoft.XMLHTTP"); };


function parseURL(){
    var params = location.search.substr(1).split('&');
    var qvpairs = {};
    var pair;var i=0;
    for (;pair=params[i]; i++) {
        var tmp = pair.split('=');
        tmp[1] = sp(unescape(tmp[1]));
        qvpairs[tmp[0]] = qvpairs[tmp[0]] ? [qvpairs[tmp[0]],tmp[1]] : tmp[1];
    }
    return qvpairs;

    function sp(str,sep){
        return str.split(sep || ",")
    }
}
