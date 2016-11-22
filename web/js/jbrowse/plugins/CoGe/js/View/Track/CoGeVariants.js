/**
 * Sub-class of HTMLFeatures that adds CoGe-specific features.
 * Added by mdb 1/15/2014
 */

var coge_variants;
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
        coge_variants = this;
        this.browser = arguments[0].browser;
    },

    // ----------------------------------------------------------------

    _exportFormats: function() {
        return [
            {name: 'GFF3', label: 'GFF3', fileExt: 'gff3'},
            {name: 'BED', label: 'BED', fileExt: 'bed'},
            { name: 'SequinTable', label: 'Sequin Table', fileExt: 'sqn' },
        	{ name: 'VCF4.1', label: 'VCF', fileExt: 'vcf' }
        ];
    },

    // ----------------------------------------------------------------

    _search_types: function() {
        var types = coge_plugin.get_checked_values('coge_search_types', 'SNP types');
        if (!types)
            return;
        var ref_seq = dojo.byId('coge_ref_seq');
        var chr = ref_seq.options[ref_seq.selectedIndex].innerHTML;
        var div = dojo.byId('coge-track-search-dialog');
        dojo.empty(div);
        div.innerHTML = '<img src="picts/ajax-loader.gif">';
        var search = {type: 'SNPs', chr: chr, snp_type: types};
        this._track.config.coge.search = search;
        var eid = this._track.config.coge.id;
        var url = api_base_url + '/experiment/' + eid + '/snps/' + chr + '?snp_type=' + types.replace(new RegExp('>', 'g'), '-');
        dojo.xhrGet({
            url: url,
            handleAs: 'json',
            load: dojo.hitch(this, function(data) {
                if (this._search_dialog)
	  				this._search_dialog.hide();
                if (data.error) {
                    coge_plugin.error('Search', data);
                    return;
                }
                if (data.length == 0) {
                    coge_plugin.error('Search', 'no SNPs found');
                    return;
                }
                coge_plugin.new_search_track(this._track, data);
            }),
            error: dojo.hitch(this, function(data) {
                if (this._search_dialog)
	  				this._search_dialog.hide();
                coge_plugin.error('Search', data);
            })
        })
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
                onClick: function(){coge_plugin.features_overlap_search_dialog(track, 'Alignments');}
            });
	        options.push({
                label: 'Find SNPs in Features',
                onClick: function(){coge_plugin.features_overlap_search_dialog(track, 'SNPs', 'snps');}
            });
            options.push({
                label: 'Find types of SNPs',
                onClick: function(){coge_variants._types_search_dialog(track);}
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

    _types_search_dialog: function(track) {
        this._track = track;
        var content = '<div id="coge-track-search-dialog"><table><tr><tr><td>Chromosome:</td><td>';
        content += coge_plugin.build_chromosome_select('Any');
        content += '</td></tr><tr><td style="vertical-align:top;">SNPs:</td><td id="coge_search_types">';
        ['A>C','A>G','A>T','C>A','C>G','C>T','G>A','G>C','G>T','T>A','T>C','T>G','deletion','insertion'].forEach(function(t) {
            content += '<div><input name="type" type="checkbox"> <label>' + t + '</label></div>';
        });
        content += '</td></tr></table>';
        content += coge_plugin.build_buttons('coge_variants._search_types()', 'coge_variants._search_dialog.hide()');
        content += '</div>';
        this._search_dialog = new Dialog({
            title: "Find types of SNPs",
            content: content,
            onHide: function() {
                this.destroyRecursive();
                coge_variants._search_dialog = null;
            }
        });
        this._search_dialog.show();
    },

    // ----------------------------------------------------------------

    updateStaticElements: function( coords ) {
        this.inherited( arguments );
        coge_plugin.adjust_nav(this.config.coge.id)
    }
});
});
