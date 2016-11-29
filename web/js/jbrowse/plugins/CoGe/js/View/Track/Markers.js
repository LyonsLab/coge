
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

        options.push({ type: 'dijit/MenuSeparator' });

		if (track.config.coge.type == 'experiment')
			options.push({
				label: 'ExperimentView',
				onClick: function(){window.open( 'ExperimentView.pl?eid=' + track.config.coge.id );}
			});

        if (track.config.coge.type != 'search')  {
            options.push({
                label: 'Find Data that Overlaps Features',
                onClick: function(){coge_plugin.features_overlap_search_dialog(track, 'Markers', 'markers');}
            });
        }
        options.push({
            label: 'Export Track Data',
            onClick: function(){coge_plugin.export_dialog(track);}
        });
        if (track.config.coge.search)
            options.push({
                label: 'Merge Markers',
                onClick: function(){coge_plugin.markers_merge_dialog(track)}
            });
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

