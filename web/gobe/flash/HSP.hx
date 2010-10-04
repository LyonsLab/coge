import flash.display.Sprite;
import flash.display.Shape;
import flash.events.MouseEvent;
import flash.geom.Point;
import flash.utils.Timer;
import flash.events.TimerEvent;
import Gobe;

/*
TODO:
1. handle clicks on highlighted HSPs
2. link to gui for toggle between line and wedge.
3. highlight currently queried HSP
*/

class HSP extends Sprite {
    public var gobe:Gobe;
    public var pair:Array<SimRect>;
    public var panel:Sprite;

    public var coords1:Array<Int>;
    public var coords2:Array<Int>;
    public var line_color:Int;
    public var as_wedge:Bool;    
    public var wedge_alpha:Float;    

    public var wedge:Wedge;    
    
    public var db_ids:Array<Int>;
 

    
    public function new(panel:Sprite, coords1:Array<Int>, coords2:Array<Int>, 
                                img1:GImage, img2:GImage,
                                line_color:Int, line_width:Int, wedge_alpha:Float, as_wedge:Bool){
        super();
        this.db_ids = [0, 0];


        this.panel = panel;
        this.coords1 = coords1;
        this.coords2 = coords2;
        this.line_color = line_color;
        this.as_wedge = as_wedge;
        this.wedge_alpha = wedge_alpha;
        this.panel.addChild(this);

        var rect1 = this.make_rect(coords1, img1);
        var rect2 = this.make_rect(coords2, img2);
        this.addChild(rect1);
        this.addChild(rect2);
        this.pair = [rect1, rect2];

        this.db_ids[0] = coords1[4];
        this.db_ids[1] = coords2[4];
        
        this.wedge = new Wedge(coords1, coords2, img1, img2, line_color, line_width, wedge_alpha, as_wedge);
        wedge.hsp = this;
        this.addChild(wedge); 
        
    }


    public function make_rect(coords:Array<Int>, img:GImage):SimRect {
        var xy = img.localToGlobal(new flash.geom.Point(coords[0],coords[1]));
        
        var db_id = coords[4]; // this links to the id in the image_data table

        var w = coords[2] - coords[0];
        var h = coords[3] - coords[1];
    
        var r = new SimRect(xy.x, xy.y, w, h, img);
        r.hsp = this;
        return r;
    }
    public function redraw(){
        this.pair[0].draw();
        this.pair[1].draw();
        this.wedge.wedge_alpha = this.wedge_alpha;
        this.wedge.as_wedge = this.as_wedge;
        this.wedge.draw();
    }
}

class MouseOverableSprite extends Sprite {
    public var hsp:HSP;
    public var mouse_over:Bool;
    public function new(){
        super();
        //this.addEventListener(MouseEvent.ROLL_OVER, onMouseOver);
        //this.addEventListener(MouseEvent.ROLL_OUT, onMouseOut);
        this.addEventListener(MouseEvent.CLICK, onClick);
    }
    public function onDblClick(e:MouseEvent){
        trace('dbl click');
    }
    public function onMouseOver(e:MouseEvent){
        //if(this.mouse_over){ return;}
        //this.mouse_over = true;
        this.hsp.pair[0].draw(true);
        this.hsp.pair[1].draw(true);
        //e.updateAfterEvent();
    }
    public function onMouseOut(e:MouseEvent){
        //if(!this.mouse_over){ return;}
        //this.mouse_over = false;
        this.hsp.pair[0].draw(false);
        this.hsp.pair[1].draw(false);
        //this.removeEventListener(MouseEvent.ROLL_OUT, onMouseOut);
        //e.updateAfterEvent();
    }

    public function onClick(e:MouseEvent){
        if(! e.shiftKey){
            this.hsp.panel.removeChild(this.hsp);
        }
        else {
           this.hsp.gobe.follow(e, this.getBounds(flash.Lib.current)); 
        }
    }
}

class Wedge extends MouseOverableSprite {
    public var line_color:Int;
    public var line_width:Int;
    public var strand:Int;
    public var as_wedge:Bool;
    public var wedge_alpha:Float;
    public var xy1a:Point;
    public var xy1b:Point;
    public var xy2a:Point;
    public var xy2b:Point;

    public function new(coords1:Array<Int>, coords2:Array<Int>, 
                                img1:GImage, img2:GImage,
                                line_color:Int, line_width:Int, wedge_alpha:Float, as_wedge:Bool){
        super();
        this.xy1a = img1.localToGlobal(new flash.geom.Point(coords1[0] - 0.75, coords1[1]));
        this.xy1b = img1.localToGlobal(new flash.geom.Point(coords1[2] - 0.75, coords1[3]));
        
        this.xy2a = img2.localToGlobal(new flash.geom.Point(coords2[0] - 0.75, coords2[1]));
        this.xy2b = img2.localToGlobal(new flash.geom.Point(coords2[2] - 0.75, coords2[3]));

        this.line_width = line_width;
        this.line_color = line_color;
        this.strand = coords1[5] < 0 ? -1 : 1;
        this.as_wedge = as_wedge;
        this.wedge_alpha = wedge_alpha;


        this.draw();
    }

    public function draw(highlight:Bool=false){
        this.graphics.clear();
        if(!this.as_wedge){
            this.draw_line(xy1a, xy1b, xy2a, xy2b);
        }
        else {
            this.draw_wedge(xy1a, xy1b, xy2a, xy2b, strand);
        }
    }

    public function draw_wedge(xy1a:Point, xy1b:Point, xy2a:Point, xy2b:Point, strand:Int){
            var g = this.graphics;
            g.beginFill(this.line_color, this.wedge_alpha);
            g.lineStyle(0, this.line_color, this.wedge_alpha > 0.5 ? 1.0: 0.6);
            g.moveTo(xy1a.x, xy1b.y);
            if (strand == 1){ 
                // go from bl1->tl2->tr2->br1->bl1 then fillRect
                g.lineTo(xy2a.x, xy2a.y);
                g.lineTo(xy2b.x, xy2a.y);
                g.lineTo(xy1b.x, xy1b.y);
            }    
            else {
                g.lineTo(xy2b.x, xy2a.y);
                g.lineTo(xy2a.x, xy2a.y);
                g.lineTo(xy1b.x, xy1b.y);
            }    
            g.endFill();
    }

    public function draw_line(xy1a:Point, xy1b:Point, xy2a:Point, xy2b:Point){
            var x1mid = (xy1a.x + xy1b.x) / 2;
            var y1mid = (xy1a.y + xy1b.y) / 2;

            var x2mid = (xy2a.x + xy2b.x) / 2;
            var y2mid = (xy2a.y + xy2b.y) / 2;

            this.graphics.lineStyle(this.line_width, this.line_color, 0.65);
            this.graphics.moveTo(x1mid, y1mid);
            this.graphics.lineTo(x2mid, y2mid);
    }
}

class SimRect extends MouseOverableSprite {
    public var color:Int;
    public var w:Float;
    public var h:Float;
    public var img:GImage;

    public function new(x:Float, y:Float, w:Float, h:Float, img:GImage) {
        super();
        this.x = x -1.15;
        this.y = y;

        this.w = w + 1.0;
        this.h = h;
        this.img = img;
        this.draw();
        
    }
    public function draw(highlight:Bool=false){
        var g = this.graphics;
        g.clear();
        g.beginFill(0x000000, 0.0);
        if(!highlight){
            g.lineStyle(0, 0x000000);
            g.drawRect(0, 0, this.w, this.h);
        }
        else {
            g.lineStyle(0, 0xaaaaaa);
            g.drawRect(0, 0, this.w, this.h);
            g.lineStyle(1, 0xcccccc);
            g.drawRect(1, 1, this.w - 2, this.h - 2);
             
        }
        g.endFill();
    }
}
