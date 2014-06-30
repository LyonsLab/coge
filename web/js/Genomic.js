/* Copyright (c) 2006 MetaCarta, Inc., published under the BSD license.
 * See http://svn.openlayers.org/trunk/openlayers/license.txt for the full
 * text of the license. */

/**
 * @class
 *
 * @requires OpenLayers/Layer/Grid1D.js
 */

OpenLayers.Layer.Grid.prototype.initGriddedTiles = function(bounds) {
    if(this.grid.length == 0){ this.grid.push([]); }

    var viewSize = this.map.getSize();
    var minCols = Math.ceil(viewSize.w/this.tileSize.w) +
            Math.max(1, 2 * this.buffer);

    var extent = this.map.getMaxExtent();
    var resolution = this.map.getResolution();

    var tileLayout = this.calculateGridLayout(bounds, extent, resolution);

    var tileoffsetx   = Math.round(tileLayout.tileoffsetx);
    var tileoffsetlon = tileLayout.tileoffsetlon;
    var tilelon       = tileLayout.tilelon;

    this.origin = new OpenLayers.Pixel(tileoffsetx, 0);

    var startX   = tileoffsetx;
    var startLon = tileoffsetlon;

    var colidx = 0;

    var layerContainerDivLeft = parseInt(this.map.layerContainerDiv.style.left);
    var layerContainerDivTop = parseInt(this.map.layerContainerDiv.style.top);

    // follow the inner loop in the original OpenLayers version.
    do {
        var tileBounds = new OpenLayers.Bounds(tileoffsetlon,
                                                0,
                                                tileoffsetlon + tilelon,
                                                1);

        //var x = tileoffsetx;
        //x -= layerContainerDivLeft;
        //var y = 0;

        var px = new OpenLayers.Pixel(tileoffsetx - layerContainerDivLeft, 0);

        var tile = this.grid[0][colidx++];
        if (!tile) {
            tile = this.addTile(tileBounds, px);
            this.addTileMonitoringHooks(tile);
            this.grid[0].push(tile);
        } else {
            tile.moveTo(tileBounds, px, false);
        }

        tileoffsetlon += tilelon;
        tileoffsetx += this.tileSize.w;
    } while ((tileoffsetlon <= bounds.right + tilelon * this.buffer)
        || colidx < minCols)

    this.removeExcessTiles(2, colidx);
    //now actually draw the tiles
    //this.linearTileLoad();
    this.spiralTileLoad();
};
OpenLayers.Layer.Grid.prototype.buffer = 1;

OpenLayers.Layer.Genomic = OpenLayers.Class.create();
OpenLayers.Layer.Genomic.prototype =
  OpenLayers.Class.inherit( OpenLayers.Layer.Grid, {

    /** Hashtable of default parameter key/value pairs
     * @final @type Object */
    DEFAULT_PARAMS: {

    },

    /**
    * @constructor
    *
    * @param {String} name
    * @param {String} url
    * @param {Object} params
    * @param {Object} options Hashtable of extra options to tag onto the layer
    */
    initialize: function(name, url, params, options) {
        var newArguments = [];
        //uppercase params
        newArguments.push(name, url, params, options);
        OpenLayers.Layer.Grid.prototype.initialize.apply(this, newArguments);
        OpenLayers.Util.applyDefaults(
                       this.params,
                       this.DEFAULT_PARAMS
                       );

        // unless explicitly set in options, if the layer is transparent,
        // it will be an overlay
        if (options == null || options.isBaseLayer == null) {
            this.isBaseLayer = ((this.params.TRANSPARENT != "true") &&
                                (this.params.TRANSPARENT != true));
        }
    },

    /**
     *
     */
    destroy: function() {
        // for now, nothing special to do here.
        OpenLayers.Layer.Grid.prototype.destroy.apply(this, arguments);
    },

    /**
     * @param {Object} obj
     *
     * @returns An exact clone of this OpenLayers.Layer.Genomic
     * @type OpenLayers.Layer.Genomic
     */
    clone: function (obj) {

        if (obj == null) {
            obj = new OpenLayers.Layer.Genomic(this.name,
                                           this.url,
                                           this.params,
                                           this.options);
        }

        //get all additions from superclasses
        obj = OpenLayers.Layer.Grid.prototype.clone.apply(this, [obj]);

        // copy/set any non-init, non-simple values here

        return obj;
    },

    /**
     * @param {OpenLayers.Bounds} bounds
     *
     * @returns A string with the layer's url and parameters and also the
     *           passed-in bounds and appropriate tile size specified as
     *           parameters
     * @type String
     */
    getZoom: function(){
        return this.map.resolutions.length - this.map.getZoom() -1;
    },
    getURL: function (bounds) {
        var zoom = this.getZoom();
        var res = this.getResolution();
        var xmin = Math.floor(Math.floor(bounds.left/res) * res);
        var xmax = Math.floor(Math.floor(bounds.right/res) * res);
        return this.getFullRequestString(
                     {
                     xmin : xmin,
                     xmax : xmax,
                     width: this.map.tileSize.w
                      });
    },

    getResolution: function(){
        return this.map.resolutions[this.map.getZoom()];
    },

    /**
    * addTile creates a tile, initializes it, and
    * adds it to the layer div.
    *
    * @param {OpenLayers.Bounds} bounds
    *
    * @returns The added OpenLayers.Tile.Image
    * @type OpenLayers.Tile.Image
    */
    addTile:function(bounds,position) {
        url = this.getURL(bounds);
        return new OpenLayers.Tile.Image(this, position, bounds,
                                             url, this.tileSize);
    },

    getExtent: function(resolution) {
        var extent = null;
        var center = this.map.getCenter();
        if (center != null) {

            if (resolution == null) {
                resolution = this.getResolution();
            }
            var size = this.map.getSize();
            var w_deg = size.w * resolution;

            extent = new OpenLayers.Bounds(center.lon - w_deg / 2,
                                           0,
                                           center.lon + w_deg / 2,
                                           1);

        }

        return extent;
    },

    /**
     * Catch changeParams and uppercase the new params to be merged in
     *  before calling changeParams on the super class.
     *
     * Once params have been changed, we will need to re-init our tiles
     *
     * @param {Object} newParams Hashtable of new params to use
     */
    mergeNewParams:function(newParams) {
        var upperParams = OpenLayers.Util.upperCaseObject(newParams);
        var newArguments = [upperParams];
        OpenLayers.Layer.Grid.prototype.mergeNewParams.apply(this,
                                                             newArguments);

        if (this.map != null) {
            this._initTiles();
        }
    },

    /** combine the layer's url with its params and these newParams.
    *
    *    Add the SRS parameter from projection -- this is probably
    *     more eloquently done via a setProjection() method, but this
    *     works for now and always.
    *
    * @param {Object} newParams
    *
    * @type String
    */
    getFullRequestString:function(newParams) {

        return OpenLayers.Layer.Grid.prototype.getFullRequestString.apply(
                                                    this, arguments);
    },

    /** @final @type String */
    CLASS_NAME: "OpenLayers.Layer.Genomic"
});
