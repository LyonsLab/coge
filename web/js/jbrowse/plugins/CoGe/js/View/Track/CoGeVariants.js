/**
 * Sub-class of HTMLVariants that adds CoGe-specific features.
 * Added by mdb 1/15/2014
 */

define( [
             'dojo/_base/declare',
             'JBrowse/View/Track/HTMLVariants'
         ],

         function(
             declare,
             HTMLVariants
         ) {
return declare( [ HTMLVariants ], {
    _exportFormats: function() {
        return [ {name: 'GFF3', label: 'GFF3', fileExt: 'gff3'}, {name: 'BED', label: 'BED', fileExt: 'bed'}, { name: 'SequinTable', label: 'Sequin Table', fileExt: 'sqn' },
        	{ name: 'VCF4.1', label: 'VCF', fileExt: 'vcf' } ];
    },
});
});
