function Panner(div_id,container,constrainY){
  this.div = G(div_id);
  container = (container)?container:div_id;
  this.elXY = [];
  this.moveY = (constrainY)? 0 : 1;

  this.qMouseDown = this.startX = this.startY = this.currentX = this.currentY = this.lastX = this.lastY = 0;

  addEvent(container,'mousedown',this._startPan,this);
  addEvent(container,'mousemove',this._onPan,this);
  addEvent(container,'mouseup',this._endPan,this);
  addEvent(container,'mouseout',this._endPan,this);

}

Panner.prototype = {
    _startPan: function(e){
        this.qMouseDown = 1; 
        this.elXY = Dom.getXY(this.div); 
        this.startXY = e.getPageXY();
        fireEvent(this.div,'EVENT_PAN_START');
        e.stop();
    },
    _onPan : function(e){       
        if(!this.qMouseDown){ return 0; }
        var pageXY = e.getPageXY();
        var dX = this.startXY[0] -  pageXY[0];
        var dY = this.startXY[1] -  pageXY[1];
        this.x_is_plus = (dX > 0)? 1 : 0;
        this.y_is_plus = (dY > 0)? 1 : 0;
        dY *= this.moveY;
        Dom.setXY(this.div, [this.elXY[0] - dX, this.elXY[1] - dY]);
        fireEvent(this.div,'EVENT_PAN');
        e.stop();
      
    },
    _endPan: function(e){
        if(!this.qMouseDown){ return; }
        this.qMouseDown = 0;
         fireEvent(this.div,'EVENT_PAN_END');
         e.stop();
    }

};

