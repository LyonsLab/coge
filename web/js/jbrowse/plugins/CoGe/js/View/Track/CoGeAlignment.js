define( [
            'dojo/_base/declare',
            'dojo/_base/array',
            'dojo/_base/lang',
            'dojo/promise/all',
            'JBrowse/Util',
            'JBrowse/View/Track/Alignments2',
            'CoGe/View/ColorDialog'
        ],
        function(
            declare,
            array,
            lang,
            all,
            Util,
            Alignment2Track,
            ColorDialog
        ) {

    return declare( [ Alignment2Track ], {
        //FIXME The track color should be fetched from the browser
        constructor: function (args) {
            this.inherited(arguments);
            var id = this.config.coge.id;
            var cookie = this.browser.cookie('track-' + this.name);

            if (!this.config.style.featureColor) {
                this.config.style.featureColor = {};
            }

            if (cookie) {
                this.config.style = lang.mixin(this.config.style, dojo.fromJson(cookie));
                this.config.histograms.color = this.config.style.featureColor[id];
            }

            if(!this.config.style.featureColor[id]) {
                this.config.style.featureColor[id] = this._getFeatureColor(id);
                this.config.histograms.color = this.config.style.featureColor[id];
            }
        },

        _getFeatureColor: function(id) {
            if (this.config.style.featureColor && this.config.style.featureColor[id])
                return this.config.style.featureColor[id];
            return coge_plugin.calc_color(id);
        },

        _trackMenuOptions: function() {
            var opts = this.inherited(arguments);
            var track = this;
            var config = this.config;

            return opts.then(function(options) {
                options.push({ type: 'dijit/MenuSeparator' });
                if (config.coge.type == 'experiment')
                    options.push({
                        label: 'ExperimentView',
                        onClick: function(){window.open( 'ExperimentView.pl?eid=' + config.coge.id );}
                    });
                options.push.apply(options, [
                    {
                    label: 'Change colors',
                    onClick: function(event) {
                        if (!track.colorDialog) {
                            track.colorDialog = new ColorDialog({
                                title: "Change colors",
                                style: { width: '230px' },
                                items: [track.config.coge],
                                featureColor: track.config.style.featureColor,
                                callback: function(id, color) {
                                    var curColor = track.config.style.featureColor[id];
                                    if (!curColor || curColor != color) {
                                        track.config.style.featureColor[id] = color;

                                        // Save color choice
                                        track.config.histograms.color = color;

                                        // Use cookie to persist color choice - mdb added 1/13/14, issue 279
                                        var cookieName = 'track-' + track.name;
                                        track.browser.cookie(cookieName, config.style);

                                        // Repaint track
                                        track.changed();

                                        //FIXME TrackList should update itself
                                        track.browser.publish('/jbrowse/v1/c/tracks/show', [track.config]);
                                    }
                                }
                            });
                        }
                        track.colorDialog.show();
                    }
                }]);

                if (track.config.coge.type == 'notebook')
                    return options;

                if (track.config.coge.type != 'search')  {
                    options.push({
                        label: 'Find Data that Overlaps Features',
                        onClick: function(){coge_plugin.features_overlap_search_dialog(track, 'Alignments');}
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
            });
        },

        // ----------------------------------------------------------------

        updateStaticElements: function( coords ) {
            this.inherited( arguments );
            coge_plugin.adjust_nav(this.config.coge.id)
        }
    });
});
