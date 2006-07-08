var kaj = {

    $ : function(id){
      if(id.constructor == Array){ return kaj.$(id[0]); /* TODO: implement */ }
      if(id.constructor != String){ return id;}
      return document.getElementById(id) || document.forms[0].elements[id] || false;
    },
    
    // Set and Get
    GS : function(id /* , value, append */){
      var el = kaj.$(id);
      var args = arguments;
      var len = arguments.length;
      if(el.type == 'select-multiple'){
         return kaj.getSelect(id); 
      }
      if(el.type == 'radio'){
        return kaj.getRadio(id);
      }
      if(el.type =='text' || el.type=='textarea' || el.type=='hidden'){
        if(len==2){
          el.value = (args[2])?el.value + args[1]:args[1];
          return 1;
        }
        return el.value;
      }
      if(len==2){
        el.innerHTML = (args[2])?el.innerHTML + args[1]:args[1];
          return 1;
      }
      return el.innerHTML;
    },
    
    G : function(id){ return kaj.GS(id) },
    S : function(id,val){ return kaj.GS(id,val) },
    
    getRadio : function(id){
      var el = $(id);
      if(el.type != 'radio'){return;}
      var ans =[];
      var els = document.getElementsByTagName('input');
      var endk = els.length;
      for(var k=0;k<endk;k++){
        if(els[k].type=='radio' && els[k].checked && (els[k].id==id || els[k].name==id)){
            ans.push(elm[k].value);
         }
      }
      return ans;
    },
    
    getSelect : function(id){
      var sel = kaj.$(id);
      if(sel.type != 'select-multiple'){
        return sel.options[sel.selectedIndex].value;
      }
      var l = sel.length;
      var ans = [];
      for (var i=0;i<l;i++) {
        if (sel[i].selected) {
          ans.push(sel[i].value);
        }
      }
      return ans;
    },
    
    map : function(arr,f){
      var l = arr.length;
      var aret = new Array(l);
      for(var i=0;i<l;i++){
        aret[i]=f(arr[i]);
      }
      return aret;
    },
    
    grep : function(ary,pattern){
      var arr = kaj.flatten(ary);
      var l = arr.length;
      var k = 0;
      var aret = [];
      for(var i=0;i<l;i++){
        var aVal = arr[i].toString();
        if(aVal.match(pattern)){
          aret[k++] = aVal;
        }
      }
      return aret;
    },

    Dumper : function(item){
      if(item.constructor == Array){
        return kaj._arrDumper(item);
      }
    },
    // internal
    _arrDumper : function(arr){
      var str = '[';
      for(var i=0;i<arr.length;i++){
        str +=  (arr[i].constructor==Array)?kaj.Dumper(arr[i]) + ',':arr[i] + ',';
      }
      return str.substr(0,str.length-1) + ']';
    },

    flatten : function(arr){
      return arr.toString().split(',');
    },
    ajax : [],
    _check_url_local : function(url){
      if(url.indexOf('http:')==-1){return url;}
    	return 'wget.pl?url=' + encodeURIComponent(url);
    	url =url.replace(/\?/g,'__kQj__');
    	url =url.replace(/&/g,'__kAj__');
    	url =url.replace(/\//g,'__kSj__');
    	return 'wget.pl?url=' + url;
    },
    
    wget_max_wait : 3000, // seconds;
    _cont : [],
    wget : function(url /*,target, handle_ajax_return */) {
        var l = kaj.ajax.length;
        var async = (arguments[1])?true:false;
    		url = kaj._check_url_local(url);
        kaj.ajax[l] = {};
        kaj.ajax[l].target = arguments[1] || null;
        kaj.ajax[l].req = (XMLHttpRequest)?new XMLHttpRequest():new ActiveXObject("Microsoft.XMLHTTP");
        kaj.ajax[l].req.open("GET", url, true);
        kaj.ajax[l].req.send(null);
        if(async){
          kaj.ajax[l].req.onreadystatechange = arguments[2] || kaj._handle_return;
        }else{
          var t0 = new Date();
          while(kaj.ajax[l].req.readyState !=4){
            if((new Date() - t0) > kaj.wget_max_wait){
              kaj.ajax.splice(l--,1);
              return 'Timed Out: Set kaj.wget_max_wait to a larger value'; 
            }
          }
          var res = kaj.ajax[l].req.responseText;
          kaj.ajax.splice(l--,1);
          return res;
        }
    },    
    
    // internal
    _handle_return : function(){
      for(var l=0;l<kaj.ajax.length;l++){
         var req = kaj.ajax[l].req;
         if(req.readyState == 4){
           var t = kaj.ajax[l].target;
           if(t == Function){
              t(req.responseText);
           }else{
              kaj.GS(t,req.responseText);
           }
           kaj.ajax.splice(l--,1);
         }
      }
    },
    
    // Pollute the namespace.
    // Export functions so kaj.$() can be
    // Accessed via: $();
    pollute : function(){
      for (var i in kaj){
        if(i == 'pollute'){continue;}
        if(i == 'ajax'){continue;}
        if(i.substr(0,1) == '_'){continue;}
        eval(i + '= kaj.' + i);
      }
   }
};
    
String.prototype.escapeHTML = function(){
  var d = document.createElement('DIV');
  var t = document.createTextNode(this);
  d.appendChild(t);
  return d.innerHTML;
};
