define([
           'dojo/_base/declare',
           'JBrowse/Store/SeqFeature',
           'JBrowse/Model/SimpleFeature'
       ],
       function( declare, SeqFeatureStore, SimpleFeature ) {

return declare( SeqFeatureStore, {

    constructor: function( args ) {
        // perform any steps to initialize your new store.  
    },

    getGlobalStats: function( statsCallback, errorCallback ) {
    	statsCallback({
    		"scoreMin" : -1,
    		"scoreMax" : 1
    	});
    },

    getFeatures: function( query, featureCallback, finishCallback, errorCallback ) {
    	var eid = this.config.eid;
    	var ref = '"' + query.ref + '"';
    	this.config.data.forEach(function(hit) {
    	   	if (hit.startsWith(ref)) {
    	   		var tokens = hit.split(',');
    	   		var start = parseInt(tokens[1]);
    	   		var end = parseInt(tokens[2]);
    	   		if (query.start <= end && query.end >= start) {
    	   			if (start == end)
    	   				end++;
    	   			var strand = parseInt(tokens[3]);
    	   			if (strand == 0)
    	   				strand = -1;
    	   			var score = strand * parseFloat(tokens[4]);
    	   			featureCallback(new SimpleFeature({ data: { id: eid, start: start, end: end, score: score } }));
    	   		}
    	   	}
    	});
    	finishCallback();
    }
});
});