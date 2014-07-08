define([
           'dojo/_base/declare',
           'dijit/focus',
           'JBrowse/View/Dialog/WithActionBar',
           'dojo/on',
           'dojo/store/Memory',
           'dijit/form/ComboBox',
           'dijit/form/Button',
           'dijit/ColorPalette',
           'dijit/form/CheckBox',
           'dijit/registry'
       ],
       function( declare, focus, ActionBarDialog, on, Memory, ComboBox, dijitButton, ColorPalette, CheckBox, registry ) {

return declare( ActionBarDialog,

    /**
     * Dijit Dialog subclass that pops up a yes/no confirmation
     * more pleasant for use as an information popup.
     * @lends JBrowse.View.ConfirmDialog
     */
{
    autofocus: false,

    constructor: function( args ) {
        //this.message = args.message || 'Select a color ...';
        this.items = args.items;
        this.featureColor = args.featureColor;
        this.callback = args.callback || function() {};
    },

    _fillActionBar: function( actionBar ) {
        var thisB = this;
        var firstItem = thisB.items[0];

        this.itemBox = new ComboBox({
            //id: "itemSelect",
            //value: firstItem.name,
            item: firstItem,
            store: new Memory({
                data: thisB.items
            }),
            style: {
            	'margin-bottom': '10px'
            },
            onChange: function() {
            	var itemId = thisB.itemBox.item.id;
            	//console.log('comboxbox change: ', itemId);
            	if (thisB.featureColor[itemId]) {
            		thisB.palette.set('value', thisB.featureColor[itemId]);
            	}
            	else {
            		registry.byNode(thisB.defaultCB.domNode).set('value', true);
            	}
            }
        }).placeAt(actionBar);

        this.palette = new ColorPalette({
            palette: "7x10",
            value: firstItem.featureColor,
            onChange: function(val) {
            	//console.log('palette change: '+val);
            	if (val) {
            		registry.byNode(thisB.defaultCB.domNode).set('value', false);
            		thisB.callback(thisB.itemBox.item.id, val);
                	//thisB.hide();
            	}
            }
        }).placeAt(actionBar);

		this.defaultCB = new CheckBox({
			//id: "defaultCB",
			checked: !firstItem.featureColor,
			style: {
				'margin-top': '10px',
			},
			onChange: function(b) {
				//console.log('checkbox change: '+b);
				if (b) {
					thisB.palette.set('value', null);
					thisB.callback(thisB.itemBox.item.id, null);
				}
			}
		}).placeAt(actionBar);

		actionBar.appendChild(
			dojo.create("label",
					{ 'for': 'defaultCB',
					  innerHTML: ' Default',
					}
			)
		);
    },

    show: function( ) {
        this.inherited( arguments );
        focus.focus( this.closeButtonNode );
    },

    _getFeatureColor: function(id) {
    	return '#' + ((((id * 1234321) % 0x1000000) | 0x444444) & 0xe7e7e7 ).toString(16); //FIXME: dup'ed in CoGe.js
    }

});
});
