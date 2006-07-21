/* STYLE rules i hope to follow:
   method_names
   _private_method_names
   propertyNames
   CONSTANTS 
   EVENT_EVTNAME_CLICK
   
*/

/* CONVENTIONS:
   rw is used to mean "Real World" coordinates as
       opposed to [p]ixel coordinates

   TIC: [T]he [I]nner [C]ontainer holds images,
       gets dragged.

   $() is a short cut for document.getElementById()

*/
function Tiler (containerId,config_arr) {
    this.div = $(containerId);
    var tmp = getElementDimensions(this.div);
    this.w = tmp.w;
    this.h = tmp.h;
    this.ZOOM_POWER = 2; 
    this.tiler = "http://biocon.berkeley.edu/CoGe/tiler.pl?";
    this._config(config_arr);
    this._initialize_TIC();
    this._initialize_tiles();

    var nc = this.rw2pix(this.INITIAL_CENTER[0],this.INITIAL_CENTER[1]);
    this.slide_by( - this.w/2 + nc[0],this.h/2 - nc[1])
    
    // TODO: bind the HTML elements with class = Tiler 
    // so that the tiler object is always 'this'
}

Tiler.prototype = {
  /**
   * Loads parameters from config the obects
   * in the configarr
   * @return {boolean} always true for now
   */
  _config: function(configarr){
      // TODO: make it so you acutally can have multiple
      // layers... 
      var queryparams = parseURL();
      var config;
      for(var c = 0;config=configarr[c]; c++){
          for(var k in config){
            var j = k;
            // number vars start with 'n'
            if(k.substr(0,1) == 'n'){
              j = k.substr(1); 
              config[j] = parseFloat(config[k]);
            }
            this[j] = queryparams[j] || config[k];
          }
          this.zoomLevel = this.MAX_ZOOM;
          try { this.zoomLevel = $(this.ZOOM_PAR_NAME).value } catch(e){}
          this.zoomLevel = queryparams[this.ZOOM_PAR_NAME] || this.zoomLevel;
          // the config file overrides the html value;
          try {$(this.ZOOM_PAR_NAME).value = this.zoomLevel;}catch(e){}
           
          /* This determines which are x,y coordinates
           * because it could be a bbox or just 2 coords.
           */
          var spars = this.SPATIAL_PARS_NAME;
          var xpars =  ypars = [];
          if(spars.length > 2){
            xpars = [spars[0],spars[2]]
            ypars = [spars[1],spars[3]]
          }else{
            xpars = spars[0];
            ypars = spars[1] || this.SPATIAL_PARS_SEP ? 'AGGREGATE' : 'ONE_DIM';
          }
          this.SPATIAL_PARS_NAME[0] = xpars;
          this.SPATIAL_PARS_NAME[1] = ypars;
          //alert(this.SPATIAL_PARS_NAME);
          var il = this.INITIAL_CENTER;
          var tmp = this.normalize_xy(il[0],il[1]);
          this._left_rw = tmp[0];
          this._top_rw = tmp[1];
        return true;
      }  
  },
  /**
   *
   * Initializes the layer to hold images
   * TIC = The Image Container
   * @return {boolean} for success 
   */
  _initialize_TIC : function(){
      var st = this.div.style;
      st.overflow = 'hidden';
       //  st.position = 'relative';
      // bigger buffer. 
      var SUPER_BUFFER = 0;
      this.nTilesWide = 1 + Math.ceil((this.w + 2 * this._div_buffer)/ this.TILE_WIDTH) + SUPER_BUFFER;
      this.nTilesHigh = (this.ONE_DIM) ? 1 
        : 1 + Math.ceil((this.h + 2 * this._div_buffer) 
          / this.TILE_HEIGHT) + SUPER_BUFFER * !this.ONE_DIM; 

      // create TIC [T]he [I]nner [C]ontainer)
      var tic = document.createElement("div");
      tic.id = 'TIC';
      tic.oncontextmenu = falsefunc;
      tic.ondragstart = falsefunc;
      connect(tic,'onmouseover',falseevent);
      connect(tic,'onhover',falseevent);
      connect(tic,'onselectstart',falseevent);

      // style stuff
      var t = tic.style;
      // TIC is aligned with the UL corner of this.div
      t.position = 'absolute';
      t.left = '0px';
      t.top  = '0px';
      // add to the div container
      this.div.appendChild(tic);

    //  tic.w = this.nTilesWide*this.TILE_WIDTH-this.w
    //  tic.h = this.nTilesHigh*this.TILE_HEIGHT-this.h
      this.TIC = tic;

      this.dragcontrol = new Panner('TIC',this.div,this.ONE_DIM);
      
      connect(this.TIC,'EVENT_PAN_START',this,'panStart');
      connect(this.TIC,'EVENT_PAN',this,'pan'); 
      connect(this.TIC,'EVENT_PAN_END',this,'panEnd'); 
  },

  get_pixel_offset : function(){
    return [ parseInt(this.TIC.style.left)
          ,  parseInt(this.TIC.style.top)];
       
  },    

  get_real_world_offset: function(){
      var pxy = this.get_pixel_offset();
      return this.pix2rw(pxy[0],pxy[1]);
  },
  slide_by: function(px,py){
     var n = 1;
     if(Math.abs(px) > 30   || Math.abs(py) > 30) n = 6;
     if(Math.abs(px) > 1000 || Math.abs(py) > 1000) n = 20;
     if(Math.abs(px) > 4000 || Math.abs(py) > 4000) n = 40;
     var sumx = px; var sumy = py;
     var dx = Math.round(px/n); var dy = Math.round(py/n);
     while(n--){
         sumx -= dx; sumy -= dy; 
         this.dragcontrol._slideBy(dx,dy); 
     }
     log('sumx,sumy: ' + sumx + ',' + sumy);
     if(Math.abs(sumx) > 3 || Math.abs(sumy) > 3){
          this.slide_by(sumx,sumy);
     }
     log('slide_by: ' + px);

  },
  setCenter: function(x,y){
      var old_center_rw = this.pix2rw(this.w/2,this.h/2);
      var diff_rw = [x - old_center_rw[0], y - old_center_rw[1]];
      var upt = this.units_per_tile();
      var diff_tiles = [Math.round(diff_rw[0]/upt[0]),
                        Math.round(diff_rw[1]/upt[1]) ];

      this._shuffle_tile_N(parseInt(diff_tiles[0]),parseInt(diff_tiles[1]));
      
      old_center_rw = this.pix2rw(this.w/2,this.h/2);
      diff_rw = [x - old_center_rw[0], y - old_center_rw[1]];
      var slide_px = this.rw_distance_to_px(diff_rw[0]);
      this.slide_by(slide_px,0);
  },

  // zoom in is tiler.zoom(1), zoomout is tiler.zoom(-1);
 zoom : function(dir){
      if(this.zoomLevel-dir < 0){return;}
      if(this.zoomLevel-dir > this.MAX_ZOOM){return;}
      var old_zoom = this.zoomLevel;
      var old_center_rw = this.pix2rw(this.w/2,this.h/2);
      this.zoomLevel-= dir;

      var tmp = this.normalize_xy(this._left_rw,this._top_rw);
      this._left_rw = tmp[0];
      this._top_rw = tmp[1]; 


      try { $(this.ZOOM_PAR_NAME).value = this.zoomLevel } catch(e){}

      var upt = this.units_per_tile();
      var new_center_rw = this.pix2rw(this.w/2,this.h/2);
      
      var diff_rw = [ new_center_rw[0]-old_center_rw[0], 
                      new_center_rw[1]-old_center_rw[1] ];

      var zre = new RegExp(this.ZOOM_PAR_NAME+'='+old_zoom,'ig' )
      var znew = this.ZOOM_PAR_NAME + '=' + this.zoomLevel;

      var rwname = this.SPATIAL_PARS_NAME[0];
      var jstr = this.magic_rw.left.join("|");
      var RE_left = new RegExp("[?&]" + rwname + '=' + "(" + jstr + ")[&$]",'g');
      var from_src;
      var sortedNodes = sorted(this.TIC.childNodes,function(a,b){
          return parseInt(a.style.left)- parseInt(b.style.left);
      });

      var seenrw = {left:{},top:{}};
      var i =0;
      forEach(sortedNodes,function(node){
           var ns = node.src;
           var l = parseInt(node.style.left );
           //var lnew = l + (tile_shift[0] * this.TILE_WIDTH);
           var t = parseInt(node.style.top );
           //var tnew = t + (tile_shift[1] * this.TILE_HEIGHT);
           //var newrw = this.pix2rw (lnew,tnew,1);
           var newrw = this.pix2rw(l,t,0);
           newrw = this.normalize_xy(newrw[0],newrw[1]);
           seenrw.left[newrw[0]]++;
           seenrw.top[ newrw[1]]++;
           ns = ns.replace (zre,znew);
           ns.match(RE_left);

           var RE_old = new RegExp(rwname + '=' + RegExp.$1, 'g');
           ns = ns.replace(RE_old,rwname + '=' + newrw[0]);
           node.src=ns;
      },this); 

      this.magic_rw = {left:[],top:[]};
      for(var k in seenrw.left){ this.magic_rw.left.push(k); }
      for(var k in seenrw.top) { this.magic_rw.top.push(k); }
      this.magic_rw.left = this.magic_rw.left.sort();
      this.magic_rw.top  = this.magic_rw.top.sort();
      log('magic_px.left: ' + this.magic_px.left);
      this.setCenter(old_center_rw[0],0);
      new_center_rw = this.pix2rw(this.w/2,this.h/2);
      var diff_rw = [ new_center_rw[0]-old_center_rw[0], 
                      new_center_rw[1]-old_center_rw[1] ];
      var diff_px = this.rw_distance_to_px(diff_rw[0]);  
      log('diff_px:' + diff_px);
      this.slide_by(-diff_px,0);

      new_center_rw = this.pix2rw(this.w/2,this.h/2);
      log("oldcenter: " + old_center_rw[0])
      log("newcenter: " + new_center_rw[0]);
      var i=3; while(this.pan() && i){i--;}
      return false;

  },
  rw_distance_to_px : function(rwdist){
     return this.rw2pix(rwdist,0,1)[0] - this.rw2pix(0,0,1)[0];
  },

  /**
   * converts pixel coords to real world 
   * coords accounting for dragging offset
   * this is used __often__ to create img.src
   * etc,
   * 
   * @return {x,y} in real world coords
   */
  rw2pix: function(gx,gy,ignoreOffset){
      var uppx = this.BASE_UNITS_PER_PIXEL_X;
      uppx *= Math.pow(this.ZOOM_POWER,this.zoomLevel);  

      gx -= parseInt(this._left_rw);
      gy -= parseInt(this._top_rw);
      px = gx/ uppx  ;
      py = gy/ uppx  ;  // TODO: for uppy  

      if(!ignoreOffset){
          var oxy = this.get_pixel_offset();
          px += oxy[0];
          py += oxy[1];
      }
      return [Math.round(px),Math.round(py)];

  },
  pix2rw: function(px,py,ignoreOffset){
      var uppx = this.BASE_UNITS_PER_PIXEL_X;
      //TODO: do for uppy
      uppx *= Math.pow(this.ZOOM_POWER,this.zoomLevel);
   //   uppt = this.units_per_tile();
      if(!ignoreOffset){
          var oxy = this.get_pixel_offset();
          px -= oxy[0];
          py -= oxy[1];
       }
      return [

               parseInt(uppx*px)+this._left_rw
             , parseInt(uppx*py)+this._top_rw
//Math.rount replacement?               parseInt(uppx*px)+this._left_rw
//             , parseInt(uppx*py)+this._top_rw
           ];
  },
  _initialize_tiles: function(){
      var nWide = this.nTilesWide;
      var nHigh = this.nTilesHigh;
      var ymin = (nHigh==1)?0:-1;

      var seenleft = {};
      var seentop = {};

      this.magic_px = {left:[],top:[]};
      this.magic_rw = {left:[],top:[]};

      this.TIC.innerHTML = '';
      for(var tx = -1; tx < nWide; tx++){
          for(var ty = ymin; ty < nHigh; ty++){
             var lt = this._initialize_and_append_tile(tx,ty); 
             seenleft[lt[0]]=1;
             seentop[lt[1]]=1;
          }
      }
      for(var k in seenleft){ 
          this.magic_px.left.push(k); 
          this.magic_rw.left.push(this.pix2rw(k,0)[0]);
      }
      for(var k in seentop){ 
          this.magic_px.top.push(k); 
          this.magic_rw.top.push(this.pix2rw(0,k)[1]);
      }

      connect(this.TIC,'ondragstart',falseevent);
      connect(this.div,'ondragstart',falseevent);
      connect(this.div,'onselectstart',falseevent);
      connect(this.TIC,'onselectstart',falseevent);
  },

  _initialize_and_append_tile: function(tilex,tiley){

      var left = tilex * this.TILE_WIDTH;
      var top = tiley * this.TILE_HEIGHT;
       
      var img = IMG({
         'height': this.TILE_HEIGHT,
         'width' : this.TILE_WIDTH,
         'class' : 'tile',
         'style' :'position:absolute;left:'+left+'px;top:'+top+'px',
         'alt'   : 'loading...',
         'src'   : this.src_from_pxy(left,top),
         'galleryimg' : 'no'
      });
      this.TIC.appendChild(img);
      return [left,top];
  },


  /**
   * if the edges of drawn tiles are within 
   * _div_buffer pixels the div edge, draw the
   * next tile.
   */
  _div_buffer: 30,

  /**
   * given a left and top in pixels,
   * return the source of the image.
   * @return {string} of img.src
   */
  src_from_tilexy : function(tl,tt){
      return this.src_from_pxy(tl*this.TILE_WIDTH,
                               tt*this.TILE_HEIGHT,1);
  },
  src_from_pxy : function(left,top){
     var url = this.tiler;
  //   var tmp = this.normalize_xy(left,top);
//     left = tmp[0]; top = tmp[1];
     var rwxy = this.pix2rw(left,top);
     rwxy  = this.normalize_xy(rwxy[0],rwxy[1]);


     // other parameters are gather from HTML elements with
     // NON_SPATIAL_PARS_NAME
     var urlpars = this.NON_SPATIAL_PARS_NAME; 
     url += this.serialize(urlpars || []);

     url += this.ZOOM_PAR_NAME ? this.serialize([this.ZOOM_PAR_NAME]) : '';

     // now add the spatial pars to the url
     var xpars = this.SPATIAL_PARS_NAME[0];
     xpars = (xpars.constructor != Array)?[xpars]:xpars;
     //TODO test for array for Y
     var ypars = this.SPATIAL_PARS_NAME[1];

     url += '&'+xpars[0]+'='+ rwxy[0];
     if(xpars.length > 1 ){     
        var locxy = this.pix2rw(left+this.TILE_WIDTH,top);
        url += '&'+xpars[1]+'='+ locxy[0];
     }
     // the Y dimension might not be needed...
     if(ypars!='ONE_DIM'){
        ypars = (ypars.constructor != Array)?[ypars]:ypars;
        url += '&'+ypars[0]+'='+ rwxy[1];
        if(ypars.constructor == Array ){     
            var rwxy =
              this.pix2rw(left,top+this.TILE_HEIGHT);
            url += '&'+xpars[1]+'='+ rwxy[1];
        }
     } 
     return url + '&';
  },
  cgi2url: function(cgi){
      var url = cgi.replace(/&/g,'/').replace(/=/g,'__').replace(/tiler\.pl\?/,'_cache_') + '.png';
      //TODO: fix for splitting large numbers;
      return url;
  },
  url2cgi: function(cgi){
     var cgi = url.replace(/\//g,'&').replace(/=/g,'__').replace(/_cache_/,'tiler.pl?');
     return cgi;
  },
     
  serialize: function(arr){
      var str = '';
      for(var i=0;i<arr.length;i++){
          var val;
          try { val = $(arr[i]).value;     }
          catch(e){ val = config[arr[i]]   }
          str+='&'+arr[i]+'='+ val;
      }
      return str;
  },    

  // cache the units per tile 
  // but what's cost of hash look up vs calcs?
  upt_cache: [],
  units_per_tile: function(){
     if(! this.upt_cache[this.zoomLevel]){
        this.upt_cache[this.zoomLevel] = this.BASE_UNITS_PER_PIXEL_X * Math.pow(this.ZOOM_POWER,this.zoomLevel) * this.TILE_WIDTH; 
     }
     return [this.upt_cache[this.zoomLevel],this.upt_cache[this.zoomLevel]];
  }, 
  tile_from_xyz: function(x,y,z){
     z = (z) ? z : this.zoomLevel;

     // units_per_tile
     var upt = this.units_per_tile(); 

     var tile_x = Math.floor( x / upt[0] ); 
     var tile_y = Math.floor( y / upt[1] ); 
     /* return the upt as well for use in
      * normalize_xy : ugly   
      * use object? {tx : tilex, ty:tiley,upt:upt} 
      */
     return [  tile_x,  tile_y, upt ];

  },
  xy_from_tile : function(tilex,tiley,z){
      z = (z) ? z : this.zoomLevel;

     // units_per_tile
      var upt = this.units_per_tile();
      return [upt[0]*tilex,upt[1]*tiley];
  },
  normalize_xy : function(x,y,z){
     z = z ? z : this.zoomLevel;
     var txy = this.tile_from_xyz(x,y,z);
     // txy[2] is units per tile;
     var upt = txy[2];
     var x = txy[0]*upt[0];
     var y = txy[1]*upt[1];
     return [x, y];
  },

  panStart : function(e){
  },

  /* dir: +1 takes the left-most tile
   *      and moves it to the right.
   *    : -1 takes the right-most tile
   *      and moves it to the left.
   *    dont look if you dont like the ternary operator
   *    and/or innerHTML. 
   */
 //   pn * this.TILE_WIDTH 
  _shuffle_tile_N: function(ntX,ntY){
     ntY = 0;

     var upt = this.units_per_tile();
     this._left_rw -= -ntX * upt[0];
     this._top_rw -= -ntY * upt[1];

      sortedNodes = sorted(this.TIC.childNodes,function(a,b){
             return parseInt(a.style.left)- parseInt(b.style.left);
      });
      var seenrw = {left:{},top:{}};
      forEach(sortedNodes,function(node){
          var ns = node.src;    

          var oldpxl = parseInt(node.style.left);
          var oldpxt  = parseInt(node.style.top)

          newrw = this.pix2rw(oldpxl,oldpxt,1);
          var rwname = this.SPATIAL_PARS_NAME[0];
          //TODO: this aint 2D compatible;
          var srcmatch = new RegExp('&' + rwname + '=([^&$]+)&', 'ig');
          ns.match(srcmatch);
          var oldrw = RegExp.$1;
          log('oldrw: ' + oldrw);
          log('newrw:'  + newrw);
          seenrw.left[newrw[0]]++;
          node.src = node.src.replace(new RegExp('&'+ rwname +'=' + oldrw + '&','g'),'&' +rwname +'=' + newrw[0] +'&');
          log(ns + '\n' + node.src);
      },this);

      this.magic_rw.left = []; 
      for(var k in seenrw.left) { this.magic_rw.left.push(k); }
      this.magic_rw.left = this.magic_rw.left.sort(function(a,b){ return a - b });

  },
  _shuffle_tile: function(xdir,ydir,imgstr){ 
     imgstr = (imgstr)? imgstr :this.TIC.innerHTML
     if(xdir > 1 || ydir > 1){ return this._shuffle_tile_N(xdir,ydir) }
     var dir = parseInt((xdir)?xdir/Math.abs(xdir):ydir/Math.abs(ydir));
     var pxdir = (xdir)?'left':'top';
     var WorH = parseInt((xdir) ? this.TILE_WIDTH: this.TILE_HEIGHT);
     // to index the [x,y] array returned from pix2rw()
     var len = this.magic_px[pxdir].length;

     // old pixel value
     var from_px = (dir==-1)?this.magic_px[pxdir].pop():this.magic_px[pxdir].shift();
     var from_src = (dir==-1)?this.magic_rw[pxdir].pop():this.magic_rw[pxdir].shift();
     
     // faster way to do this. how??
     var to_px =
     (dir==-1)?parseInt(this.magic_px[pxdir][0])-WorH:parseInt(this.magic_px[pxdir][len-2])+WorH;
     var idx = (xdir)?0:1;
     var to_src = this.pix2rw(to_px,0,true)[idx];

     if(dir==-1){
        this.magic_px[pxdir].unshift(to_px);
        this.magic_rw[pxdir].unshift(to_src);
     }else{
        this.magic_px[pxdir].push(to_px);
        this.magic_rw[pxdir].push(to_src);
     }

     var nrw = this.magic_rw.left;

     imgstr = imgstr.replace(/&amp;/g,'&');

     var rwname = (xdir) ? this.SPATIAL_PARS_NAME[0] : this.SPATIAL_PARS_NAME[1];
     // the regexp to move images around
     var reg_px  = new RegExp(pxdir + ':\\s*'+from_px+'px','ig');
     // the regexp to change the src 
     var reg_src = new RegExp('&' + rwname + '='+from_src + '[&$]', 'ig');
     imgstr = imgstr.replace(reg_src,'&' + rwname + '='+to_src+'&' )
     imgstr = imgstr.replace(reg_px, pxdir + ':' + to_px+'px')

     return imgstr;
  },
  brXY : function(){ return this.pix2rw(this.w-this.TILE_WIDTH,this.h-this.TILE_HEIGHT);},


  pan: function(){
      // the amout TIC has moved since last shuffle;
      var nW = this.nTilesWide -1;
      var p2rw = this.pix2rw(-30,-30);
      var brxy = this.brXY();
      var xchange = 0;
      if(p2rw[0]<this.magic_rw.left[0]){
          // rotate right-most tile to left;
          --xchange;
          log('*** right ***');
      }
      // account for the fact that we wand the right edge of tile.
      
      if(brxy[0] > this.magic_rw.left[nW]){
          // rotate left-most tile to right; 
          log('*** left ***'); 
          ++xchange;
      } 
      if(xchange) return this.TIC.innerHTML = this._shuffle_tile(xchange,0);

      // dont need to check up down if its
      // not a 2D thingy
      if(this.ONE_DIM){ return false }


      var ychange = 0;
      var nH = this.nTilesHigh -1;
      if(p2rw[1]<this.magic_rw.top[0]){
          // rotate right-most tile to left;
          ++ychange;
          log('*** up ***');
      }
      // account for the fact that we want the top edge of tile.
      
      if(brxy[1] > this.magic_rw.top[nH]){
          // rotate left-most tile to right; 
          log('*** down ***'); 
          --ychange;
      } 


      if(ychange) return this.TIC.innerHTML = this._shuffle_tile(0,ychange);
      return false;

  }, 
  panEnd : function(){

  }
  

};
