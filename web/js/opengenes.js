OpenLayers.GenomeBrowser = OpenLayers.Map;

OpenLayers.GenomeBrowser.prototype.pan = function(dx, dy) {
    // getCenter
    var centerPx = this.getViewPortPxFromLonLat(this.getCenter());

    // adjust
    var newCenterPx = centerPx.add(dx, 0);

    // only call setCenter if there has been a change
    if (newCenterPx.x != centerPx.x) {
        var newCenterLonLat = this.getLonLatFromViewPortPx(newCenterPx);
        // this is required or the tools wont cause map to pan
        newCenterLonLat.lat = 0
        this.setCenter(newCenterLonLat);
    }
};

OpenLayers.GenomeBrowser.prototype.centerLayerContainer = function(lonlat){
            var originPx = this.getViewPortPxFromLonLat(this.layerContainerOrigin);
            var newPx = this.getViewPortPxFromLonLat(lonlat);

            if ((originPx != null) && (newPx != null)) {
                this.layerContainerDiv.style.left = (originPx.x - newPx.x) + "px";
            }

};

OpenLayers.Control.PanZoom.prototype.draw = function(px) {
    // initialize our internal div
    OpenLayers.Control.prototype.draw.apply(this, arguments);
    px = this.position.clone();

    // place the controls
    this.buttons = new Array();

    var sz = new OpenLayers.Size(18,18);
    var centered = new OpenLayers.Pixel(px.x+sz.w/2, px.y);
    px.y = centered.y+sz.h;
    this._addButton("panleft", "west-mini.png", px, sz);
    this._addButton("panright", "east-mini.png", px.add(sz.w, 0), sz);
    this._addButton("zoomin", "zoom-plus-mini.png",
                    centered.add(0, sz.h*3+5), sz);
    this._addButton("zoomworld", "zoom-world-mini.png",
                    centered.add(0, sz.h*4+5), sz);
    this._addButton("zoomout", "zoom-minus-mini.png",
                    centered.add(0, sz.h*5+5), sz);
    return this.div;
};
OpenLayers.Control.PanZoomBar.prototype.draw = function(px) {
    // initialize our internal div
    OpenLayers.Control.prototype.draw.apply(this, arguments);
    px = this.position.clone();

    // place the controls
    this.buttons = new Array();
    var sz = new OpenLayers.Size(18,18);
    var centered = new OpenLayers.Pixel(px.x+sz.w/2, px.y);

    px.y = centered.y+sz.h;
    this._addButton("panleft", "west-mini.png", px, sz);
    this._addButton("panright", "east-mini.png", px.add(sz.w, 0), sz);
    this._addButton("zoomin", "zoom-plus-mini.png", centered.add(0, sz.h*3+5), sz);
    centered = this._addZoomBar(centered.add(0, sz.h*4 + 5));
    this._addButton("zoomout", "zoom-minus-mini.png", centered, sz);
    return this.div;
};

OpenLayers.Control.MousePosition.prototype.redraw = function(evt) {

    var lonLat;
    if (evt == null) {
        lonLat = new OpenLayers.LonLat(0, 0);
    } else {
        if (this.lastXy == null ||
            Math.abs(evt.xy.x - this.lastXy.x) > this.granularity)
        {
            this.lastXy = evt.xy;
            return;
        }

        lonLat = this.map.getLonLatFromPixel(evt.xy);
        this.lastXy = evt.xy;
    }

    var digits = parseInt(this.numdigits);
    var newHtml = parseInt(lonLat.lon) + this.map.units;

    if (newHtml != this.element.innerHTML) {
        this.element.innerHTML = newHtml;
    }
};

/**************************************************************/
/*  WHEN GETTING A PINK TILE, RELOAD FROM A DIFFERENT SERVER  */
/**************************************************************/
OpenLayers.IMAGE_RELOAD_ATTEMPTS = 3;
OpenLayers.Util.onImageLoadError = function() {
    this._attempts = (this._attempts) ? (this._attempts + 1) : 1;
    if(this._attempts <= OpenLayers.IMAGE_RELOAD_ATTEMPTS) {
        var urls = this.layer.url;
        if(urls instanceof Array && urls.length > 1){
            var src = this.src.toString();
            var current_url, k;
            for (k = 0; current_url = urls[k]; k++){
                if(src.indexOf(current_url) != -1){
                    break;
                }
            }
            var guess = Math.floor(urls.length * Math.random())
            var new_url = urls[guess];
            k = 0;
            while(new_url == current_url && k++ < 4){
                guess = Math.floor(urls.length * Math.random())
                new_url = urls[guess];
            }
            this.src = src.replace(current_url, new_url);
        } else {
            this.src = this.src;
        }
    } else {
        this.style.backgroundColor = OpenLayers.Util.onImageLoadErrorColor;
    }
    this.style.display = "";
};
OpenLayers.Tile.Image.prototype.destroy = function() {
        if (this.imgDiv != null)  {
            OpenLayers.Event.stopObservingElement(this.imgDiv.id);
            if (this.imgDiv.parentNode == this.frame) {
                this.frame.removeChild(this.imgDiv);
                this.imgDiv.map = null;
            }
            this.imgDiv.layer = null;
        }
        this.imgDiv = null;
        if ((this.frame != null) && (this.frame.parentNode == this.layer.div)) {
            this.layer.div.removeChild(this.frame);
        }
        this.frame = null;
        OpenLayers.Tile.prototype.destroy.apply(this, arguments);
    };

OpenLayers.Tile.Image.prototype.draw = function() {
        if (this.layer != this.layer.map.baseLayer && this.layer.reproject) {
            this.bounds = this.getBoundsFromBaseLayer(this.position);
        }
        if (!OpenLayers.Tile.prototype.draw.apply(this, arguments)) {
            return false;
        }

        if (this.isLoading) {
            //if we're already loading, send 'reload' instead of 'loadstart'.
            this.events.triggerEvent("reload");
        } else {
            this.isLoading = true;
            this.events.triggerEvent("loadstart");
        }

        if (this.imgDiv == null) {
            this.initImgDiv();
        }
        // ADDED FOR CHANGING URLS:
        this.imgDiv.layer = this.layer;

        this.imgDiv.viewRequestID = this.layer.map.viewRequestID;

        this.url = this.layer.getURL(this.bounds);
        // position the frame
        OpenLayers.Util.modifyDOMElement(this.frame,
                                         null, this.position, this.size);

        var imageSize = this.layer.getImageSize();
        if (this.layerAlphaHack) {
            OpenLayers.Util.modifyAlphaImageDiv(this.imgDiv,
                    null, null, imageSize, this.url);
        } else {
            this.imgDiv.src = this.url;
            OpenLayers.Util.modifyDOMElement(this.imgDiv,
                    null, null, imageSize) ;
        }
        return true;
    };

/********************************************************/

OpenLayers.BasePair = function(bp){
    var ll = new OpenLayers.LonLat(bp,0);
    ll.basepair = ll.x;
    return ll;
}

function jfetch(url,t,o) {
  var req = jfetch.xhr();
  req.open("GET",url,true);
  req.onreadystatechange = function() {
    if(req.readyState == 4){
      var rsp = req.responseText;
      if(t.constructor == Function) return t.apply(o,[rsp]);
      t = document.getElementById(t);
      t[t.value ==undefined ? 'innerHTML': 'value'] = rsp;
      req = null;
    }
  };
  req.send(null);
}
jfetch.xhr =
    (window.ActiveXObject)
   ? function(){ return new ActiveXObject("Microsoft.XMLHTTP"); }
   : function(){ return new XMLHttpRequest()};
