define([
           'dojo/_base/declare',
           'JBrowse/Store/SeqFeature'
       ],
       function( declare, SeqFeatureStore ) {

return declare(SeqFeatureStore, {

    constructor: function(args) {
        // perform any steps to initialize your new store.  
    },

    getGlobalStats: function(statsCallback, errorCallback) {
    	statsCallback({
    		"scoreMin" : -1,
    		"scoreMax" : 1
    	});
    },

    getFeatures: function(query, featureCallback, finishCallback, errorCallback) {
    	var features = this.config.results.get_features(query.ref, query.start, query.end);
    	features.forEach(function(feature){ featureCallback(feature); });
    	finishCallback();
    }
});
});