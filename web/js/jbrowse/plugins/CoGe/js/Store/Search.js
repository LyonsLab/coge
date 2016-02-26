define([
           'dojo/_base/declare',
           'dojo/_base/array',
           'dojo/request/xhr',
           'JBrowse/Store/SeqFeature',
           'JBrowse/Model/SimpleFeature'
       ],
       function( declare, array, xhr, SeqFeatureStore, SimpleFeature ) {

return declare( SeqFeatureStore, {

    constructor: function( args ) {
        // perform any steps to initialize your new store.  
    },

    getGlobalStats: function( statsCallback, errorCallback ) {
        var thisB = this;
        xhr.get( this.config.baseUrl+'my/webservice/url',
                 { handleAs: 'json' }
               ).then(
                   function( data ) {
                       data = thisB._transformStatsDataSomehow( data );
                       statsCallback( data );
                   },
                   errorCallback
               );
    },

    getRegionStats: function( query, statsCallback, errorCallback ) {
        var thisB = this;
        xhr.get( this.config.baseUrl+'my/other/webservice/url',
                 { handleAs: 'json', query: query }
               ).then(
                   function( data ) {
                       data = thisB._transformStatsDataSomehow( data );
                       statsCallback( data );
                   },
                   errorCallback
               );
    },

    getFeatures: function( query, featureCallback, finishCallback, errorCallback ) {
        var thisB = this;
        xhr.get( this.config.baseUrl+'my/features/webservice/url',
                 { handleAs: 'json', query: query }
               ).then(

                   function( featuredata ) {

                       // transform the feature data into feature
                       // objects and call featureCallback for each
                       // one. for example, the default REST
                       // store does something like:
                       array.forEach( featuredata || [],
                           function( featureKeyValue ) {
                               var feature = new SimpleFeature({
                                       data: featureKeyValue
                                   });
                               featureCallback( feature );
                           });

                       // call the endCallback when all the features
                       // have been processed
                       finishCallback();
                   },

                   errorCallback
               );

    }
});
});