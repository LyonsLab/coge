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

    _create_search_dialog: function() {
    	var content = '<div><table><tr><tr><td>Chromosome:</td><td><select id="coge_ref_seq"><option>Any</option>';
    	this.browser.refSeqOrder.forEach(function(rs){
    		content += '<option>' + rs + '</option>';
    	})
    	content += '</select></td></tr><tr><td style="vertical-align:top;">Features:</td><td id="coge_search_features">';
    	var features = coge_track_list._track_configs.filter(function(f) {
    		return (f.coge.type && f.coge.type == 'features');
    	});
    	features.forEach(function(t) {
    		content += '<div><input type="checkbox"> <label>' + t.coge.id + '</label></div>';
    	});
    	content += '</td></tr><tr><td></td><td><button onClick="coge_track_list._check_all(this.parentNode.parentNode.parentNode,true)">check all</button> <button onClick="coge_track_list._check_all(this.parentNode.parentNode.parentNode,false)">uncheck all</button></td></tr></table>';
    	content += '<div class="dijitDialogPaneActionBar"><button data-dojo-type="dijit/form/Button" type="button" onClick="coge_variants._search_features()">OK</button><button data-dojo-type="dijit/form/Button" type="button" onClick="coge_variants._search_dialog.hide()">Cancel</button></div></div>';
        this._search_dialog = new Dialog({
            title: "Find SNPs in Features",
            content: content,
            onHide: function(){
            	this.destroyRecursive();
            	coge_variants._search_dialog = null;
            },
            style: "width: 300px"
        });
    	this._search_dialog.show();
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

	_search_features: function() {
    	var features = document.getElementById('coge_search_features').getElementsByTagName('INPUT');
    	var types = [];
    	for (var i=0; i<features.length; i++)
    		if (features[i].checked)
    			types.push("'" + features[i].nextElementSibling.innerText + "'");
    	if (!types.length) {
    		coge.error('Search', 'Please select one or more feature types to search.');
    		return;
    	}
    	var url = api_base_url + '/genome/' + gid + '/snps?features=' + (types.length == features.length ? 'all' : types.join());
    	var ref_seq = dojo.byId('coge_ref_seq');
    	if (ref_seq.selectedIndex > 0)
    		url += '&chr=' + ref_seq.options[ref_seq.selectedIndex].innerHTML;
    	dojo.xhrGet({
    		url: url,
    		handleAs: 'json',
	  		load: dojo.hitch(this, function(data) {
	  			if (data.error) {
	  				coge.error('Search', data);
	  				return;
	  			}
	  			if (data.length == 0) {
	  				coge.error('Search', 'no SNPs found');
	  				return;
	  			}
    		}),
    		error: function(data) {
    			coge.error('Search', data);
    		}
    	})
    },

    // ----------------------------------------------------------------

    _trackMenuOptions: function() {
        var options = this.inherited(arguments);
        var track = this;

        options.push({
	        label: 'Search',
	        onClick: function(){coge_variants._create_search_dialog(track);}
        });

        return options;
    }
});
});
