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

   G() is a short cut for document.getElementById()

*/

function Tiler (containerId,config_arr) {
    this.div = G(containerId);
    this.w = parseInt(Dom.getStyle(this.div,'width'));
    this.h = parseInt(Dom.getStyle(this.div,'height'));
    this.ZOOM_POWER = 2; 
    this.tiler = "tiler.pl?";

    this._config(config_arr);
    this._initialize_TIC();
    this._initialize_tiles();
    
    // bind the HTML elements so that the tiler object is always
    // available as 'this' in anything with class='tileClass';
    var els = getElementsByClassName('tileClass');
    for(var i=0;i<els.length;i++){
        els[i].scope = this;
     }
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
      for(var c = 0; c<configarr.length; c++){
          var config = configarr[c];
          for(var k in config){
            var j = k;
            // number vars start with 'n'
            if(k.substr(0,1) == 'n'){
              j = k.substr(1); 
              config[j] = parseFloat(config[k]);
            }
            this[j] = config[k];
          }
          this.zoomLevel = G(this.ZOOM_PAR_NAME).value;
          this.zoomLevel = this.zoomLevel ? this.zoomLevel:  this.INITIAL_ZOOM ;
          // the config file overrides the html value;
          G(this.ZOOM_PAR_NAME).value = this.zoomLevel;
          
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
            ypars = spars[1] || 'ONE_DIM';
          }

          this.SPATIAL_PARS_NAME_X = xpars;
          this.SPATIAL_PARS_NAME_Y = ypars;
          var il = (G('x') != undefined )? G('x').value : this.INITIAL_LEFT;
          var it = this.INITIAL_TOP || 0;
          var tmp = this.normalize_xy(il,it,this.zoomLevel);
          this.INITIAL_LEFT = tmp[0];
          this.INITIAL_TOP = tmp[1];
          this._2D = (this.ONE_DIM)?0:1;
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
      st.position = 'relative';
      // bigger buffer. 
      var SUPER_BUFFER = 1;
      this.nTilesWide = 1 + Math.ceil((this.w + 2 * this._div_buffer)/ this.TILE_WIDTH) + SUPER_BUFFER;
      this.nTilesHigh = !this.ONE_DIM + Math.ceil((this.h + 2 * this._div_buffer) / this.TILE_HEIGHT) + SUPER_BUFFER * !this.ONE_DIM; 
      // create TIC [T]he [I]nner [C]ontainer)
      // TODO: see about leakage with mixing append and
      // innerHTML
      var tic = document.createElement("div");
      tic.id = 'TIC';

      // style stuff
      var t = tic.style;
      // TIC is aligned with the UL corner of this.div
      t.position = 'absolute';
      t.left = '0px';
      t.top  = '0px';
      /* could apply this change to each img style.left and .top ...
      t.left = (this.w - this.nTilesWide * this.TILE_WIDTH)/2;
      t.top  = (this.h - this.nTilesHigh * this.TILE_Height)/2 * !this.ONE_DIM;
      */

      // add to the div container
      this.div.appendChild(tic);
    //  tic.w = this.nTilesWide*this.TILE_WIDTH-this.w
    //  tic.h = this.nTilesHigh*this.TILE_HEIGHT-this.h
      this.TIC = tic;

      this.dragcontrol = new Panner('TIC',this.div,1);
      
      //TODO: subscribe to events
        addEvent(this.TIC,'EVENT_PAN_START',this.panStart,this);
        addEvent(this.TIC,'EVENT_PAN',this.pan,this); 
        addEvent(this.TIC,'EVENT_PAN_END',this.panEnd,this); 
  },
  get_pixel_offset : function(){
    return [ parseInt(this.TIC.style.left) - -this.left_shift
          ,  parseInt(this.TIC.style.top) - -this.top_shift];
       
  },    
  left_shift:0,top_shift:0,

  get_real_world_offset: function(){
      var pxy = this.get_pixel_offset();
      return this.pix_to_real_world(pxy[0],pxy[1]);
  },

  // zoom in is tiler.zoom(1), zoomout is tiler.zoom(-1);
  zoom : function(dir){
      //TODO: see TODO below re xshift et al.
      if(this.zoomLevel-dir < 0){return;}
      var old_zoom = this.zoomLevel;
      var octr = this.pix_to_real_world(this.w/2,this.h/2)[0];


      var nTIC = this.TIC.innerHTML;
      DEBUG(this.magic_rw.left.join(",") + "<br/>");

      var old_center_rw = this.pix_to_real_world(this.w/2,this.h/2);
      this.zoomLevel-=dir;

      var new_center_rw = this.pix_to_real_world(this.w/2,this.h/2);

      var new_center_rw_normalized 
          = this.normalize_xy(new_center_rw[0]
          ,                   new_center_rw[1]);
      var diff_rw = [
          new_center_rw[0]-old_center_rw[0], 
          new_center_rw[1]-old_center_rw[1] 
      ];


      var modx = diff_rw[0];
      var abspx = this.rw_distance_to_px(modx);
      var modpx = abspx % this.TILE_WIDTH;
      var xshift = abspx - modpx;


      this.left_shift += modpx;

      //TODO: repeat above for Y

     var zre = new RegExp( this.ZOOM_PAR_NAME + '=' + old_zoom,'g' )
     nTIC = nTIC.replace(zre , this.ZOOM_PAR_NAME + '=' + this.zoomLevel);



      var rwname = this.SPATIAL_PARS_NAME[0];

      var from_src;
      DEBUG("<br/>old: " + this.magic_rw.left.join(",") + '<br/>');
      for(var m = 0; from_rw=this.magic_rw.left[m];m++){
          var from_px = this.magic_px.left[m];
          var to_px = from_px - -xshift;
          var to_rw = this.pix_to_real_world(to_px,0,1)[0];
         
          var reg_rw = new RegExp(rwname + '='+ from_rw,'g');
          nTIC = nTIC.replace(reg_rw,rwname + '=' + to_rw);
 
          reg_px = new RegExp('left:\\s*' + from_px,'g'); 
          nTIC = nTIC.replace(reg_px,'left:' + to_px);
          this.magic_px.left[m]=to_px;
          this.magic_rw.left[m]=to_rw;
      } 
      //TODO: start here, it's the combination of xshift/abssfhit,
      // this.TIC.style.left,this.left_shift;
      this.TIC.style.left = parseInt(this.TIC.style.left) 
         -  -parseInt(xshift) + 'px';
      
      this.TIC.innerHTML = nTIC; 
     while( this.pan() ){}
      DEBUG("<br/>new: " + this.magic_rw.left.join(",") + '<br/>');
      DEBUG(
      'modx:' + modx +  '<br/>' +
      'abspx:' + abspx +'<br/>' +
      'modpx:' + modpx +'<br/>' +
      'xshift:' + xshift 
      );


      
      var nctr = this.pix_to_real_world(this.w/2,this.h/2)[0];
      DEBUG("<br/>octr: " + octr + '\nnctr: ' + nctr + '<br/>');
      //TODO: use these to readjust in the
      // while (this.pan() loop )... hacky, but 
      // easy...

      return false;

  },

  rw_distance_to_px : function(rwdist){
     var BX = this.BASE_UNITS_PER_PIXEL_X; 
     return Math.ceil(rwdist/(BX * Math.pow(this.ZOOM_POWER,this.zoomLevel)));
  },


  /**
   * converts pixel coords to real world 
   * coords accounting for dragging offset
   * this is used __often__ to create img.src
   * etc,
   * 
   * @return {x,y} in real world coords
   */
  pix_to_real_world:  function(px,py,ignoreOffset){
      var uppx = this.BASE_UNITS_PER_PIXEL_X;
      //TODO: do for uppy
      uppx *= Math.pow(this.ZOOM_POWER,this.zoomLevel);
   //   uppt = this.units_per_tile();
      if(!ignoreOffset){
          var oxy = this.get_pixel_offset();
          px -= oxy[0];
        //  upt[0] -= parseInt(this.TIC.style.left);
          py -= oxy[1];
         // upt[1] -= parseInt(this.TIC.style.top);
       }
      return [ parseInt(uppx*px)+this.INITIAL_LEFT
             , parseInt(uppx*py)+this.INITIAL_TOP
           ];
  },
  _initialize_tiles: function(){
      var nWide = this.nTilesWide;
      var nHigh = this.nTilesHigh;
      var ymin = (nHigh==1)?0:-1;
      this.magic_px = {left:[],top:[]};
      this.magic_rw = {left:[],top:[]};

      for(var tx = -1; tx < nWide; tx++){
          for(var ty = ymin; ty < nHigh; ty++){
             this._initialize_and_append_tile(tx,ty); 

          }
      }
      var txt = this.TIC.innerHTML;
      var re;
      var i = 0;
      // make sure we start at the beginning
      RegExp.lastIndex = 0;
      // to check for dups.
      var seenhash = {};
      // initialize the magic_px
      while(/(left|top):\s*(-*\d+)px?/g.exec(txt) != null) {
        var key = RegExp.$1 + RegExp.$2;
        if(! seenhash[key] ){
            seenhash[key] = 1;
            this.magic_px[RegExp.$1].push(RegExp.$2); 
            if(RegExp.$1 == 'left'){
              this.magic_rw.left.push(this.pix_to_real_world(RegExp.$2,0)[0]);
            }else{
              this.magic_rw.top.push(this.pix_to_real_world(0,RegExp.$2)[1]);
            }
        }
      }
  },
  draw_from_magic_rw: function(new_magic){
      var TIC = this.TIC.innerHTML;
      var rwname = this.SPATIAL_PARS_NAME[0];
      var rwleft;
      for(var i=0;rwleft=this.magic_rw.left[i];i++){
          var reg = new RegExp(rwname + '=' + rwleft);
          TIC = TIC.replace(reg,rwname + '=' + new_magic.left[i]);
          this.magic_rw.left[i] = new_magic.left[i];
      }
      // TODO: for y;
      this.TIC.innerHTML = TIC;          
  },

  magic_rw_replace : function(from,to,TIC,rwname){
      var rwname = this.SPATIAL_PARS_NAME[0];
      var reg =  new RegExp(rwname + '=' + from ,'g');
      return TIC.replace(reg,rwname +'=' + to);
  },

  _initialize_and_append_tile: function(tilex,tiley){

      var left = tilex * this.TILE_WIDTH;
      var top = tiley * this.TILE_HEIGHT;


      var img = "<img height="+this.TILE_HEIGHT;
      img+= ' class="tile"';
//    img+= " id="+left+'_'+top;
      img+= ' style=position:absolute;left:'+ left +'px;top:'+top +'px;';
      img+= ' src=' + this.src_from_pxy(left,top) ;

      var rw = this.debug_passer;
      img+= " width="+this.TILE_WIDTH + ' alt="loading... "';
      img += '/>';
      this.TIC.innerHTML+=img;
      k=100000
      while(k--){}
  },
  // horrible hack to pass the rw coords
  debug_passer:'',


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
  src_from_pxy : function(left,top){
     var url = this.tiler;
   //     url = 'http://biocon.berkeley.edu/CoGe/GenomePNG.pl?&simple=1';
     var rwxy = this.pix_to_real_world(left,top);

     // other parameters are gather from HTML elements with
     // NON_SPATIAL_PARS_NAME
     var urlpars = this.NON_SPATIAL_PARS_NAME; 
     url += this.serialize(urlpars);
     url += this.serialize([this.ZOOM_PAR_NAME]);

     // now add the spatial pars to the url
     var xpars = this.SPATIAL_PARS_NAME_X;
     var ypars = this.SPATIAL_PARS_NAME_Y;

     this.debug_passer = rwxy[0];
     url += '&'+xpars[0]+'='+ rwxy[0];
     if(xpars.constructor == Array ){     
        var locxy = this.pix_to_real_world(left+this.TILE_WIDTH,top);
        url += '&'+xpars[1]+'='+ locxy[0];
     }

     // the Y dimension might not be needed...
     if(ypars!='ONE_DIM'){
        url += '&'+ypars[0]+'='+ rwxy[1];
        if(ypars.constructor == Array ){     
            var rwxy =
              this.pix_to_real_world(left,top+this.TILE_HEIGHT);
             this.debug_passer = rwxy.join(",");
            url += '&'+xpars[1]+'='+ rwxy[1];
        }
     } 
    // alert(url);
     return url;
  },
  serialize: function(arr){
      var str = '';
      for(var i=0;i<arr.length;i++){
          str+='&'+arr[i]+'='+G(arr[i]).value;
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

     // units_per_pixel
     //var upp = this.BASE_UNITS_PER_PIXEL_X * Math.pow(this.ZOOM_POWER,z); 

     // units_per_tile
     var upt = this.units_per_tile(); // upp * this.TILE_WIDTH;

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

     // units_per_pixel
     //var upp = this.BASE_UNITS_PER_PIXEL_X * Math.pow(this.ZOOM_POWER,z); 
      
     // units_per_tile
     var upt = this.units_per_tile(); //upp * this.TILE_WIDTH;
     return [upt[0]*tilex,upt[1]*tiley];
  },
  normalize_xy : function(x,y,z){
     z = z ? z : this.zoomLevel;
     var txy = this.tile_from_xyz(x,y,z);
     // txy[2] is units per tile;
     var upt = txy[2];
     var x = txy[0]*upt[0];
     var y = txy[1]*upt[1];
     return [ x, y ];
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
//TODO: change this to be called as: ('left',-1);
// then can use                      (pxdir,pn) as:
 //   pn * this.TILE_WIDTH 
 
  _shuffle_tile: function(xdir,ydir){ 
     var dir = (xdir)?xdir:ydir;
     var pxdir = (xdir)?'left':'top';
     var WorH = (xdir) ? this.TILE_WIDTH: this.TILE_HEIGHT;
     // to index the [x,y] array returned from pix_to_real_world()
     var len = this.magic_px[pxdir].length;

     var orw = this.magic_rw.left;
     orw = orw.toString();
     // old pixel value
     var from_px = (dir==-1)?this.magic_px[pxdir].pop():this.magic_px[pxdir].shift();
     var from_src = (dir==-1)?this.magic_rw[pxdir].pop():this.magic_rw[pxdir].shift();
     
     // faster way to do this. how??
     var to_px =
     (dir==-1)?parseInt(this.magic_px[pxdir][0])-WorH:parseInt(this.magic_px[pxdir][len-2])+WorH;
     var idx = (xdir)?0:1;
     var to_src = this.pix_to_real_world(to_px,0,true)[idx];

     if(dir==-1){
        this.magic_px[pxdir].unshift(to_px);
        this.magic_rw[pxdir].unshift(to_src);
     }else{
        this.magic_px[pxdir].push(to_px);
        this.magic_rw[pxdir].push(to_src);
     }

     var nrw = this.magic_rw.left;
     // DEBUG('<br/>' + orw + '<br/>' + nrw.toString());
     imgstr = this.TIC.innerHTML;

//  G('loc').innerHTML += 'from src: ' + from_src +'\nto src: '+ to_src + '<br/>';
     
     var rwname = (xdir) ? this.SPATIAL_PARS_NAME[0] : this.SPATIAL_PARS_NAME[1];
     // the regexp to move images around
     var reg_px  = new RegExp(pxdir + ':\\s*'+from_px+'px','g');
     // the regexp to change the src 
     var reg_src = new RegExp(rwname + '='+from_src, 'g');

     imgstr = imgstr.replace(reg_src,rwname + '='+to_src+'')
     imgstr = imgstr.replace(reg_px, pxdir + ':' + to_px+'px')

     this.TIC.innerHTML = imgstr;
//    G('loc').innerHTML += imgstr.escapeHTML().replace(/&amp;/g,'&');
  },

  pan : function(){
      // the amout TIC has moved since last shuffle;
     
      var tlXY = this.pix_to_real_world(-30,-30);
      var nW = this.nTilesWide -1;
      var nH = this.nTilesHigh -1;
      if(tlXY[0] < this.magic_rw.left[0]){
          // rotate right-most tile to left;
          this._shuffle_tile(-1,0);
          return 'right';
      }
      // account for the fact that we wand the right edge of tile.
      var brXY = this.pix_to_real_world(this.w-this.TILE_WIDTH,this.h-this.TILE_HEIGHT);
      if(brXY[0] > this.magic_rw.left[nW]){
          // rotate left-most tile to right; 
          this._shuffle_tile(1,0);
          return 'left';
      } 

      // dont need to check up down if its
      // not a 2D thingy
      if(this.ONE_DIM){ return false }
      var yo = 0;

      if(yo < this.tbBuffer - this.TIC.h + dbuff ){
          // rotate top-most tile to bottom;
          alert( 'top');
          return 'bottom';
      }

      if(yo  > this.tbBuffer - dbuff ){
          // rotate bottom-most tile to top;
          alert( 'top');
          return 'top';
      }
      return false;

  }, 
  panEnd : function(){

  }
  

};

function G(id){return document.getElementById(id);}

String.prototype.escapeHTML = function(){
    var div = document.createElement('div');
    var text = document.createTextNode(this);
    div.appendChild(text);
    return div.innerHTML;
};
var DEBUGLOC = 'loc';
var NDEBUG=0;

function DEBUG(str){
//  NDEBUG++;
//  if(NDEBUG>10){ NDEBUG=0;G(DEBUGLOC).innerHTML='';}
//  G(DEBUGLOC).innerHTML += str;
}
