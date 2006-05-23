// it's ajax yo.

/* 
USAGE: 
to return the result to a div:

jfetch('myscript.pl?var=something','divId');

to return the result to a fuction:

jfetch('myscript.php?n=1234',myFunc);


where:
function myFunc(){
  alert(arguments[0]);
}


where myscript.php something like:
<?php
  $answer = $_REQUEST['n'];
  print 'i got ' + $answer;
?>

report bugs/questions to:
brent
bpederse@nature.berkeley.edu

*/

function jfetch(url,target,context) {
   var req = (XMLHttpRequest)?new XMLHttpRequest():new ActiveXObject("Microsoft.XMLHTTP");
   req.open("GET",url,true);
   req.onreadystatechange = function() {
      if(req.readyState==4){
			  if(target.constructor == Function){
				   target.call(context,req.responseText);
				}else{
					document.getElementById(target).innerHTML = req.responseText;
        }
      }
   };
   req.send(null);
}
