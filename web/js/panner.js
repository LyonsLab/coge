function Panner(div_id,container,constrainY){
  this.div = $(div_id);
  container = (container)?container:div_id;
  this.container = container;
  this.elXY = [];
  this.moveY = (constrainY)? 0 : 1;
  this.qMouseDown = 0;
  connect(container,'onmousedown',this,'_startPan',this);
  connect(container,'onmousemove',this,'_onPan',this);
  connect(container,'onmouseup',this,'_endPan',this);
  connect(container,'onmouseout',this,'_endPan',this);
}

Panner.prototype = {
    _startPan: function(e){
        this.qMouseDown = 1; 
        e.stop();
        this.elXY = getElementPosition(this.div,this.container); 
        this.startXY = e.mouse().client;
        signal(this.div,'EVENT_PAN_START');
    },
    _onPan : function(e){       
        if(!this.qMouseDown){ return 0; }
        var pageXY = e.mouse().client;
        e.stop();
        var dX = this.startXY.x -  pageXY.x;
        var dY = this.startXY.y -  pageXY.y;
        setElementPosition(this.div, {x: this.elXY.x - dX, y: (this.elXY.y - dY) * this.moveY});
        signal(this.div,'EVENT_PAN');
      
    },
    _endPan: function(e){
        if(!this.qMouseDown){ return; }
        e.stop();
        this.qMouseDown = 0;
        signal(this.div,'EVENT_PAN_END');
    },
    _slideTo: function(px,py){
        setElementPosition(this.div,{x:px,y:py});
        signal(this.div,'EVENT_PAN');
    },
    _slideBy: function(px,py){
        var pos = getElementPosition(this.div,this.container);
        pos.x -= px;
        pos.y -= py;
        pos.y *= 0;
        setElementPosition(this.div,{x:pos.x,y:pos.y});
        signal(this.div,'EVENT_PAN');
    }

};
