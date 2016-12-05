define([
           'dojo/_base/declare',
           'dojo/dom-construct',
           'dijit/focus',
           'JBrowse/View/Dialog/WithActionBar',
           'dojo/on',
           'dojo/store/Memory',
           'dijit/form/Select',
           'dijit/form/Button',
           'dijit/ColorPalette',
           'dijit/form/CheckBox',
           'dijit/registry'
       ],
       function( declare, domConstruct, focus, ActionBarDialog, on, Memory, Select, dijitButton, ColorPalette, CheckBox, registry ) {

return declare( ActionBarDialog,

    /**
     * Dijit Dialog subclass that pops up a yes/no confirmation
     * more pleasant for use as an information popup.
     * @lends JBrowse.View.ConfirmDialog
     */
{
    autofocus: false,

    constructor: function( args ) {
        this.track = args.track;
        this.items = args.items;
        this.featureColor = args.featureColor;
        this.callback = args.callback || function() {};
    },

    _fillActionBar: function( actionBar ) {
        var self = this;
 
        if (this.items.length == 1)
            domConstruct.place(domConstruct.toDom('<div style="margin-bottom:10px">' + this.items[0].label + '</div>'), actionBar);
        else
            this._select = new Select({
                options: this.items,
                style: { 'margin-bottom': '10px' },
                onChange: function() {
                    var item_id = self._select.value;
                    if (self.featureColor[item_id])
                        self.palette.set('value', self.featureColor[item_id]);
                    else
                        registry.byNode(self.defaultCB.domNode).set('value', true);
                }
            }).placeAt(actionBar);

        this.palette = new ColorPalette({
            palette: "7x10",
            value: self.featureColor[self.items[0].value],
            style: { 'margin-bottom': '10px' },
			onChange: function(val) {
            	if (val) {
                    if (self.defaultCB)
                		registry.byNode(self.defaultCB.domNode).set('value', false);
            		self.callback(self._select ? self._select.value : self.items[0].value, val);
            	}
            }
        }).placeAt(actionBar);

		this.defaultCB = new CheckBox({
			checked: !self.featureColor[self.items[0].value],
			onChange: function(b) {
				if (b) {
					self.palette.set('value', null);
					self.callback(self._select ? self._select.value : self.items[0].value, null);
				}
			}
		}).placeAt(actionBar);

		actionBar.appendChild(
			dojo.create("label", { 'for': 'defaultCB', innerHTML: ' Default', })
		);
    },

    onHide: function() {
        this.destroyRecursive();
        this.track._color_dialog = null;
    },

    show: function( ) {
        this.inherited( arguments );
        focus.focus( this.closeButtonNode );
    },
});
});
