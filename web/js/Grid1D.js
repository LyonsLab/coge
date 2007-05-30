/* Copyright (c) 2006 MetaCarta, Inc., published under the BSD license.
 * See http://svn.openlayers.org/trunk/openlayers/license.txt for the full
 * text of the license. */

/**
 * @class
 * 
 * @requires OpenLayers/Layer/HTTPRequest.js
 */
OpenLayers.Layer.Grid1D = OpenLayers.Class.create();
OpenLayers.Layer.Grid1D.prototype = 
  OpenLayers.Class.inherit( OpenLayers.Layer.Grid, {
    
    /** This function is called whenever the map is moved. All the moving
     * of actual 'tiles' is done by the map, but moveTo's role is to accept
     * a bounds and make sure the data that that bounds requires is pre-loaded.
     * 
     * @param {OpenLayers.Bounds} bounds
     * @param {Boolean} zoomChanged
     * @param {Boolean} dragging
     */
    moveTo:function(bounds, zoomChanged, dragging) {
        OpenLayers.Layer.HTTPRequest.prototype.moveTo.apply(this, arguments);
        
        if (bounds == null) {
            bounds = this.map.getExtent();
        }
        if (bounds != null) {
            if (!this.grid.length || zoomChanged 
                || !this.getGridBounds().containsBounds(bounds, true)) { 
                this._initTiles();
            } else {
                var buffer = (this.buffer) ? this.buffer*1.5 : 1;
                while (true) {
                    var tlLayer = this.grid[0].position;
                    var tlViewPort = this.map.getViewPortPxFromLayerPx(tlLayer);
                    if (tlViewPort.x > -this.tileSize.w * (buffer - 1)) {
                        this.shiftColumn(true);
                    } else if (tlViewPort.x < -this.tileSize.w * buffer) {
                        this.shiftColumn(false);
                    } else {
                        break;
                    }
                }

                if (this.buffer == 0) {
                        for (var c=0, cl=row.length; c<cl; c++) {
                            var tile = this.grid[c];
                            if (!tile.drawn && tile.bounds.intersectsBounds(bounds, false)) {
                                tile.draw();
                            }
                        }
                }
            }
        }
    },
    
    /**
     * @private
     * 
     * @returns A Bounds object representing the bounds of all the currently 
     *           loaded tiles (including those partially or not at all seen 
     *           onscreen)
     * @type OpenLayers.Bounds
     */
    getGridBounds:function() {
        
        var LeftTile = this.grid[0];

        var right = this.grid.length - 1; 
        var RightTile = this.grid[right];

        return new OpenLayers.Bounds(LeftTile.bounds.left, 
                                     0,
                                     RightTile.bounds.right, 
                                     1);
    },

    /**
     * @private
     */
    _initTiles:function() {
        var viewSize = this.map.getSize();
        var minCols = Math.ceil(viewSize.w/this.tileSize.w) + 1;
        var bounds = this.map.getExtent();
        var extent = this.map.getMaxExtent();
        var resolution = this.map.getResolution();
        var tilelon = resolution * this.tileSize.w;
        
        var offsetlon = bounds.left - extent.left;
        var tilecol = Math.floor(offsetlon/tilelon) - this.buffer;
        var tilecolremain = offsetlon/tilelon - tilecol;
        var tileoffsetx = -tilecolremain * this.tileSize.w;
        var tileoffsetlon = extent.left + tilecol * tilelon;

       /* changed HERE */ 
        
        tileoffsetx = Math.round(tileoffsetx); // heaven help us

        var startX = tileoffsetx; 
        var startLon = tileoffsetlon;


        var colidx = 0;

        do {
            var tileBounds = new OpenLayers.Bounds(tileoffsetlon, 
                                                    0, 
                                                    tileoffsetlon + tilelon,
                                                    1);

            var x = tileoffsetx;
            x -= parseInt(this.map.layerContainerDiv.style.left);

            var y = 0;

            var px = new OpenLayers.Pixel(x, y);
            var tile = this.grid[colidx++];
            if (!tile) {
                tile = this.addTile(tileBounds, px);
                this.grid.push(tile);
            } else {
                tile.moveTo(tileBounds, px, false);
            }
    
            tileoffsetlon += tilelon;       
            tileoffsetx += this.tileSize.w;
        } while ((tileoffsetlon <= bounds.right + tilelon * this.buffer)  
            || colidx < minCols)
        
        while (this.grid.length > colidx){
            var tile = this.grid.pop()
            tile.destroy();;
        }
        //now actually draw the tiles
        //this.linearTileLoad();
        this.outwardTileLoad();
    },
    
    /** 
     * @private 
     * 
     *   Starts at the top right corner of the grid and proceeds in a spiral 
     *    towards the center, adding tiles one at a time to the beginning of a 
     *    queue. 
     * 
     *   Once all the grid's tiles have been added to the queue, we go back 
     *    and iterate through the queue (thus reversing the spiral order from 
     *    outside-in to inside-out), calling draw() on each tile. 
     */
    linearTileLoad: function() {
        for(var i=0; i < this.grid.length; i++) {
            var tile = this.grid[i]
            tile.draw();
        }
    },

    outwardTileLoad: function(){
        var ntiles = this.grid.length;
        var center = Math.floor(ntiles/2);
        this.grid[center].draw();
        var i = -1, j=1
        if( center != ntiles/2) {
            this.grid[center+1].draw();
            j = 2;
        }

        
        for(;j+center < this.grid.length ; i--,j++){
            this.grid[center+i].draw()
            this.grid[center+j].draw()
        }
    },

    /** go through and remove all tiles from the grid, calling
     *    destroy() on each of them to kill circular references
     * 
     * @private
     */
    clearGrid:function() {
        if (this.grid) {
            for(var i=0; i < this.grid.length; i++) {
                this.grid[i].destroy();
            }
        }        
    },


    /**
     * @private
     * 
     * @param {Boolean} prepend if true, prepend to beginning.
     *                          if false, then append to end
     */
    shiftColumn: function(prepend) {
        var deltaX = (prepend) ? -this.tileSize.w : this.tileSize.w;
        var resolution = this.map.getResolution();
        var deltaLon = resolution * deltaX;

        var modelTileIndex = (prepend) ? 0 : (this.grid.length - 1);
        var modelTile = this.grid[modelTileIndex];
        
        var bounds = modelTile.bounds.clone();
        var position = modelTile.position.clone();
        bounds.left = bounds.left + deltaLon;
        bounds.right = bounds.right + deltaLon;
        position.x = position.x + deltaX;

        var tile = prepend ? this.grid.pop() : this.grid.shift()
        tile.moveTo(bounds, position);
        if (prepend) {
            this.grid.unshift(tile);
        } else {
            this.grid.push(tile);
        }
    },
    
    /** @final @type String */
    CLASS_NAME: "OpenLayers.Layer.Grid1D"
});
