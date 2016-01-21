/**
* Copyright (c) 2010, CyVerse Collaborative, University of Arizona, Cold Spring Harbor Laboratories, University of Texas at Austin
* This software is licensed under the CC-GNU GPL version 2.0 or later.
* License: http://creativecommons.org/licenses/GPL/2.0/
*/
var SSWAP = function() {
   var pipelineUI = "http://sswap.iplantcollaborative.org/ipc/create-pipeline-with-rrg";

   function discover(jsonRRG, containerId) {
      $(document).ready(function() {
         var jsonString = null;

         if (typeof(jsonRRG) == "object") {
           jsonString = JSON.stringify(jsonRRG);
         }
         else {
           jsonString = jsonRRG;
         }

         var container = $(containerId);

         var form = $("<form>");
         form.attr("method","post").attr("action",pipelineUI);

         var rrgField = $("<input>");
rrgField.attr("type","hidden").attr("name","rrg").attr("value",jsonString);

         form.append(rrgField);

         var submit = $("<input>");
submit.attr("type","image").attr("src","http://sswap.iplantcollaborative.org/discover.png").attr("width","32").attr("height","32").attr("value","Start Pipeline with this data");

         form.append(submit);

         container.append(form);
      });
   }

   return {
       discover : discover
   };
}();
