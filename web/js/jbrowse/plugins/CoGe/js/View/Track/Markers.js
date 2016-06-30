
var coge_markers;
define( [
             'dojo/_base/declare',
             'dijit/Dialog',
            'JBrowse/View/Track/HTMLFeatures'
         ],

         function(
             declare,
             Dialog,
             HTMLFeatures
         ) {
return declare( [ HTMLFeatures ], {
    constructor: function() {
        this.inherited(arguments); // call superclass constructor
        coge_markers = this;
        this.browser = arguments[0].browser;
    },

    // ----------------------------------------------------------------

    _trackMenuOptions: function() {
        var options = this.inherited(arguments);
        var track = this;

        if (track.config.coge.type == 'notebook')
            return options;

        if (!track.config.coge.search_track)  {
            options.push({
                label: 'Find Markers in Features',
                onClick: function(){coge_plugin.features_overlap_search_dialog(track, 'Markers', 'markers');}
            });
        }
        options.push({
            label: 'Export Track Data',
            onClick: function(){coge_plugin.export_dialog(track);}
        });
        if (track.config.coge.search)
            options.push({
                label: 'Save Results as New Experiment',
                onClick: function(){coge_plugin.save_as_experiment_dialog(track)}
            });
        return options;
    },

    // ----------------------------------------------------------------

    updateStaticElements: function( coords ) {
        this.inherited( arguments );
        coge_plugin.adjust_nav(this.config.coge.id)
    }
});
});

