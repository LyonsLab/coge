import java.awt.*;
import javax.swing.*;
import java.util.*;
import java.awt.event.*;
import RangeTree.*;

public class XYplot extends JInternalFrame {
    
    private boolean SEE = false;
    private boolean showDiagonals = true;
    private EvalueSlider maxLogscore_adjust = null;
    private float maxLogScore = -5; // init value
    private int maxNumHits = 20; // init value.
    
    private Vector matches = null;
    
    private Hashtable coordToMatch = null;
    private Hashtable selectedMatches; // keys match ID to boolean.
    private Hashtable matchCounter = null;

    private RangeTree rangeTree = null;

    private int max_x_bp_coord = 0;
    private int max_y_bp_coord = 0;
    
    // current image width and height, influenced by zoom factor.
    private int image_height;
    private int image_width;
    private int maxInitDimension; // max initial dimension for image.
    
    // original image width and height, needed to reset zoom to original view.
    private int original_image_width;
    private int original_image_height;



    private XYcanvas plotArea;
    private JScrollPane plotScroller;
    private SequenceTicker y_axis_ticker;
    private SequenceTicker x_axis_ticker;
    
    private float zoomFactor = 1;
    
    public static boolean TEST = false;

    // viewport parameters
    private int viewPort_x;
    private int viewPort_y;
    
    //////////////////////////////////////////////////////////////////////////////
    // Constructor.
    public XYplot (int maxD, String title) {
	
	super (title, true, true, true, true);

	this.maxInitDimension = maxD;
	
	// add a menu bar:
	createMenuBar();
	
	// add a canvas to paint to.
	plotArea = new XYcanvas();
		
	plotScroller = new JScrollPane(plotArea);
	
	y_axis_ticker = new SequenceTicker(this, SequenceTicker.VERTICAL);
	x_axis_ticker = new SequenceTicker(this, SequenceTicker.HORIZONTAL); 
	
	
	plotScroller.setColumnHeaderView(x_axis_ticker);
	plotScroller.setRowHeaderView(y_axis_ticker);
	plotScroller.setBackground(new Color(174,114,228));
	
	this.getContentPane().add(plotScroller);
	
	this.coordToMatch = new Hashtable();
	
	this.selectedMatches = new Hashtable();
	this.matchCounter = new Hashtable();
	
	// for unit testing purposes only.
	if (TEST) {
	    image_height = maxD;
	    image_width = image_height;
	    original_image_width = image_width;
	    original_image_height = image_height;
	    this.setSize(image_width,image_height);
	    this.setVisible(true);

	    JFrame testFrame = new JFrame("testing XYplot");

	    testFrame.addWindowListener(new WindowAdapter() {
		    public void windowClosing(WindowEvent e) {
			System.exit(0);
		    }	
		});
	    
	    JDesktopPane desktop = new JDesktopPane();
	    desktop.add(this);

	    testFrame.getContentPane().add(desktop);
	    
	    testFrame.setSize(image_width+200, image_height + 300);
	    testFrame.setVisible(true);
	    
	}
	
    }
    
    //////////////////////////////////////////////////////////////////////////
    public Dimension getPlotDimension () {
	return (new Dimension(image_width, image_height));
    }
    

    ///////////////////////////////////////////////////////////////////////////
    public void setMatches (Vector matches) {
	
	this.matches = matches;
	
	int numMatches = matches.size();
	RTpoint [] pointList = new RTpoint [numMatches];
	for (int i = 0; i < numMatches; i++) {
	    
	    DatabaseMatch match = (DatabaseMatch) matches.elementAt(i);
	    
	    SequenceFeature feat1 = match.feat1();
	    SequenceFeature feat2 = match.feat2();
	    
	    String acc1 = feat1.get_accession();
	    String acc2 = feat2.get_accession();
	    if (matchCounter.containsKey(acc1)) {
		Integer currentCount = (Integer) matchCounter.get(acc1);
		int val = currentCount.intValue();
		matchCounter.put (acc1, new Integer( ++val) );
	    } else {
		matchCounter.put (acc1, new Integer(1));
	    }
	    
	    if (matchCounter.containsKey(acc2)) {
		Integer currentCount = (Integer) matchCounter.get(acc2);
		int val = currentCount.intValue();
		matchCounter.put (acc2, new Integer( ++val ) );
	    } else {
		matchCounter.put (acc2, new Integer(1));
	    }
	    	    
	    
	    int end5_1 = (int) feat1.end5();
	    int end3_1 = (int) feat1.end3();
	    if (end5_1 > max_y_bp_coord || end3_1 > max_y_bp_coord) {
		if (end5_1 > max_y_bp_coord) { max_y_bp_coord = end5_1; }
		if (end3_1 > max_y_bp_coord) { max_y_bp_coord = end3_1; }
	    }

	    int end5_2 = (int) feat2.end5();
	    int end3_2 = (int) feat2.end3();
	    if (end5_2 > max_x_bp_coord || end3_2 > max_x_bp_coord) {
		if (end5_2 > max_x_bp_coord) { max_x_bp_coord = end5_2; }
		if (end3_2 > max_x_bp_coord) { max_x_bp_coord = end3_2; }
	    }

	    // track the (x,y) coordinates
	    int mid_bp_y = (int) (end5_1 + end3_1)/2;
	    int mid_bp_x = (int) (end5_2 + end3_2)/2;
	    
	    String pointKey = String.valueOf(mid_bp_x) + "," + String.valueOf(mid_bp_y);
	    coordToMatch.put(pointKey, match);

	    pointList[i] = new RTpoint(mid_bp_x, mid_bp_y);
	    
	}
	
	// add to range tree.
	RangeTree.mergeSort(pointList, 0);
	rangeTree = new RangeTree(pointList);
	
	
	// calculate image dimensions:
	float pixels_per_bp;
	if (max_x_bp_coord > max_y_bp_coord) {
	    pixels_per_bp = (float) maxInitDimension/max_x_bp_coord;
	} else {
	    pixels_per_bp = (float) maxInitDimension/max_y_bp_coord;
	}
	
	image_width = (int) (max_x_bp_coord * pixels_per_bp);
	original_image_width = image_width;
	
	image_height = (int) (max_y_bp_coord * pixels_per_bp);
	original_image_height = image_height;
	
	if (SEE) {
	    System.out.println ("Original image size: w:" + String.valueOf(image_width) + ", h:" + 
				String.valueOf(image_height));
	}
	
	// ticker management.
	x_axis_ticker.set_max_seq_length(max_x_bp_coord);;
	y_axis_ticker.set_max_seq_length(max_y_bp_coord);
	
	y_axis_ticker.setPreferredHeight(image_height);
	x_axis_ticker.setPreferredWidth(image_width);
	

	this.setSize(getPreferredSize());
    
	boolean visible = this.plotArea.isVisible();
	
	if (SEE) {
	    System.out.println ("visible canvas: " +  String.valueOf(visible));
	}
	
	//this.plotArea.setSize(image_width, image_height);
	//this.plotArea.repaint();
	repaint();
	
    }
    


    ///////////////////////////////////////////////////////////////////////////
    public void setXYaxisLabels (String xLabel, String yLabel) {
	
	if (x_axis_ticker != null) {
	    x_axis_ticker.setAxisLabel(xLabel);
	}

	if (y_axis_ticker != null) {
	    y_axis_ticker.setAxisLabel(yLabel);
	}
    }
    

    ////////////////////////////////////////////////////////////////////////////
    private class XYcanvas extends JPanel  
	implements MouseListener, MouseMotionListener {
	
	boolean startBand = false;
	int startX, startY;
	int dragX, dragY;
	
	int mouseButton;

	Rectangle bandedRegion;

	public XYcanvas () {
	    
	    bandedRegion = new Rectangle();
	    addMouseListener(this);
	    addMouseMotionListener(this);
	    
	}
	

	/////////////////////////////////////////////////////////////////////
	public void paintComponent (Graphics g) {
	    
	    super.paintComponent(g);

	    // get viewport dimensions
	    Dimension viewDimension = XYplot.this.plotScroller.getViewport().getExtentSize();
	    float viewable_height_ratio = (float) viewDimension.height / viewDimension.width;
	     

	    if (SEE) {
		System.out.println ("Viewable region: (width:" + String.valueOf(viewDimension.width) +
				    ", height: " + String.valueOf(viewDimension.height));
	    }
	    

	    // color the entire panel area in purple.  
	    // (size diff from image size only when image is smaller than panel size)
	    int viewable_width = this.getWidth();
	    int viewable_height = this.getHeight();
	    g.setColor(new Color(174,114,228));
	    g.fillRect(0,0,viewable_width,viewable_height);
	    

	    if (SEE) {
		System.out.println ("Panel region: (width:" + String.valueOf(viewable_width) +
				    ", height: " + String.valueOf(viewable_height));
	    }

	    int width = image_width; 
	    int height = image_height; 
	    
	    // fill xy-plot area in white
	    g.setColor(Color.white);
	    g.fillRect(0,0,width,height);
	    

	    if (SEE) {
		System.out.println ("Plot image (width: " + String.valueOf(width) + ", height: "
				    + String.valueOf(height) + ")\n\n");
	    }
	    
	    if (TEST) {
		// draw an X
		g.setColor(Color.black);
		g.drawLine(0,0,width,height);
		g.drawLine(0,height, width,0);
	    }

	    if (matches != null) {
		g.setColor(Color.black);
		//draw each match:
		Boolean selected;
		int numMatches = matches.size();
		for (int i = 0; i < numMatches; i++) {
		    
		    DatabaseMatch match = (DatabaseMatch) matches.elementAt(i);
		    
		    if (! scoreWithinRange(match)) {
			continue;
		    }
		    
		    Point pixels = matchToPixels(match, width, height);

		    g.drawOval(pixels.x, pixels.y, 1, 1);
		    
		    if (showDiagonals && match.inDiagonal()) {
			g.setColor(match.getColor());
			g.fillOval(pixels.x-2, pixels.y-2, 5, 5);
			g.drawOval(pixels.x-6, pixels.y-6, 12, 12);
			g.setColor(Color.black); // reset to normal draw mode.
		    }
		    
		    selected = (Boolean) selectedMatches.get(new Integer(match.id()));
		    if (selected != null) {
			g.setColor(Color.blue);
			g.drawOval(pixels.x-2, pixels.y-2, 6, 6);
			g.setColor(Color.black); // reset to normal draw mode.
		    }

		    
		}

	    }
	    
	    
	    // repaint the tickers so they correspond to new canvas/coordinate dimensions
	    x_axis_ticker.repaint();
	    y_axis_ticker.repaint();
	    
	    

	    if (startBand) {
		g.setColor(Color.red);
		// determine coordinates for rectangle drawing.
		int x,y,w,h;
		if (startX > dragX) {
		    x = dragX;
		    w = startX - dragX;
		} else {
		    x = startX;
		    w = dragX - startX;
		}

		if (startY > dragY) {
		    y = dragY;
		    h = startY - dragY;
		} else {
		    y = startY;
		    h = dragY - startY;
		}
		
		// restrict banded region to same ratio as viewable area:
		h = (int) (w * viewable_height_ratio);
		

		// draw band selection
		bandedRegion.setRect(x,y,w,h);
		g.drawRect(x,y,w,h);
		if (SEE) {
		    System.out.println ("\n\ndrawing rect: (" + String.valueOf(x) +
					", " + String.valueOf(y)
					+ ", " + String.valueOf(w)
					+ ", " + String.valueOf(h) + ", " 
					+ String.valueOf(viewable_height_ratio) + ")\n\n");
		}
	    }
	}
	

	////////////////////////////////////////////
	private Point matchToPixels (DatabaseMatch match, int width, int height) {
	    
	    SequenceFeature feat1 = match.feat1();
	    SequenceFeature feat2 = match.feat2();
	    
	    long end5_1 = feat1.end5();
	    long end3_1 = feat1.end3();
	    
	    long end5_2 = feat2.end5();
	    long end3_2 = feat2.end3();
	    
	    int mid_bp_y = (int) (end5_1 + end3_1)/2;
	    int mid_bp_x = (int) (end5_2 + end3_2)/2;
	    
	    int pixel_x = transform_x_coord(mid_bp_x, width);
	    int pixel_y = transform_y_coord(mid_bp_y, height);


	    return (new Point (pixel_x, pixel_y));
	}
	


	//////////////////////////////////////////////////////////////////////////////////
	public Dimension getPreferredSize() {
	    return (new Dimension(image_width,image_height));
	}
	
	
	//////////////////////////////////////////////////////////////////////////////////
	private int transform_x_coord (int coord, int width) {
	    
	    return ((int) ((coord/(float)max_x_bp_coord) * width));
	}
	

	/////////////////////////////////////////////////////////////////////////////////
	private int transform_y_coord (int coord, int height) {
	    
	    // coord = max_y_bp_coord - coord + 1; // no longer revcomp coordinates.
	    
	    return ((int) (coord/(float)max_y_bp_coord * height));
	    
	}
	



	// Mouse event handling:
	//////////////////////////////////////////////////////////////////////////////////
	public void mousePressed (MouseEvent evt) {
	    
	    // initiating band selection.
	    startX = dragX = evt.getX();
	    startY = dragY = evt.getY();
	    
	    System.out.println (evt);
	    

	    if (SEE) {
		System.out.println("Mouse pressed (" + String.valueOf(startX) 
				   + ", " + String.valueOf(startY));
	    }


	    if (startBand) {
		return;
	    }
	    	   
	    
	    startBand = true;



	}

	//////////////////////////////////////////////////////////////////////////////////
	public void mouseReleased (MouseEvent evt) {
	    
	    String modifiers = evt.getMouseModifiersText(evt.getModifiers());
	    System.out.println ("modifiers: " + modifiers);

	    int button;
	    mouseButton = evt.getButton();
	    if (modifiers.indexOf("Button1") > -1 ) { //mouseButton == evt.BUTTON1) {
		System.out.println ("mouse button 1 pressed.");
		button = 1;
	    } else if (modifiers.indexOf("Button2") > -1 ) { // mouseButton == evt.BUTTON2) {
		System.out.println ("mouse button 2 pressed.");
		button = 2;
	    } else if (modifiers.indexOf("Button3") > -1) { //mouseButton == evt.BUTTON3) {
		System.out.println ("mouse button 3 pressed.");
		button = 3;
	    } else {
		System.out.println ("mouse button ? pressed.");
		button = 0;
	    }
	    
	    startBand = false;
	    if (button == 3) {
		XYplot.this.zoomOnCoordinates(bandedRegion);
	    } else if (button == 1) {
		XYplot.this.selectedRegion(bandedRegion);
		
	    }
	}
	
	
	////////////////////////////////////////////////////////////////////////////////
	public void mouseDragged (MouseEvent evt) {
	    
	    dragX = evt.getX();
	    dragY = evt.getY();
	    repaint();
	
	} 

	// empty methods for mouseEvent, mouseListener
	public void mouseClicked(MouseEvent evt) { }  
	public void mouseEntered(MouseEvent evt) { }  
	public void mouseExited(MouseEvent evt) { }   
	public void mouseMoved(MouseEvent evt) { }
	
	
 	
    }

    
    ///////////////////////////////////////////////////////////////////////////////////
    private void createMenuBar () {
	
	//instantiate menubar
    	JMenuBar menuBar = new JMenuBar();
    
    	// Create a JMenu
    	JMenu menuZoom = new JMenu("Zoom...");
	
    	JMenuItem menuZoomItem = new JMenuItem("zoom (+)");
    	
    	ActionListener zoomListener = new ActionListener() { 
    		public void actionPerformed(ActionEvent e) {  
     			System.out.println ("Zoom (+) option selected.");
			XYplot.this.zoom(2);
     		} 
	    };
    	
    	menuZoomItem.addActionListener(zoomListener);
    	menuZoom.add(menuZoomItem);
	
	
	JMenuItem menuZoomItem2 = new JMenuItem("zoom (-)");
    	
    	ActionListener zoomListener2 = new ActionListener() { 
    		public void actionPerformed(ActionEvent e) {  
		    System.out.println ("Zoom (-) option selected.");
		    XYplot.this.zoom(0.5f);
     		} 
    	};
    	
    	menuZoomItem2.addActionListener(zoomListener2);
    	menuZoom.add(menuZoomItem2);
	
	JMenuItem menuZoomItem3 = new JMenuItem("reset");
    	
    	ActionListener zoomListener3 = new ActionListener() { 
    		public void actionPerformed(ActionEvent e) {  
		    System.out.println ("Zoom Reset option selected.");
		    XYplot.this.resetImage();
     		} 
    	};
    	
    	menuZoomItem3.addActionListener(zoomListener3);
    	menuZoom.add(menuZoomItem3);
	
	menuBar.add(menuZoom);

	// Create separate entry for showing diagonals turn on/off.
	JMenu diagMenu = new JMenu("Diagonals");
	JMenuItem toggleDiagsItem = new JMenuItem("Toggle show diagonals.");
	ActionListener toggleDiagsItemListener = new ActionListener() {
		public void actionPerformed ( ActionEvent e) {
		    boolean currentSetting = XYplot.this.showDiagonals;
		    if (currentSetting == true) {
			XYplot.this.showDiagonals = false;
		    } else {
			XYplot.this.showDiagonals = true;
		    }
		    XYplot.this.repaint();
		}
	    };
	toggleDiagsItem.addActionListener(toggleDiagsItemListener);
	diagMenu.add(toggleDiagsItem);
	menuBar.add(diagMenu);
	
	
	// Create separate entry for showing diagonals turn on/off.
	JMenu matchMenu = new JMenu("MatchDisplay");
	JMenuItem matchMenuItem = new JMenuItem("Show matches by E-value.");
	ActionListener matchMenuItemListener = new ActionListener() {
		public void actionPerformed ( ActionEvent e) {
		    if (maxLogscore_adjust == null) {
			XYplot.this.launch_maxLogscoreAdjuster();
		    }
		}
	    };
	matchMenuItem.addActionListener(matchMenuItemListener);
	matchMenu.add(matchMenuItem);
	menuBar.add(matchMenu);
	this.setJMenuBar(menuBar);

    }
    
    
    ////////////////////////////////////////////////////////////////////////////////////////
    public void zoom (float zoomFactor) {

	this.zoomFactor = zoomFactor;
	image_width = (int)(image_width*zoomFactor);
	image_height = (int)(image_height*zoomFactor);
	
	this.plotArea.setSize(image_width, image_height);
	this.x_axis_ticker.setPreferredWidth(image_width);
	this.y_axis_ticker.setPreferredHeight(image_height);
	
	System.out.println("Zooming by factor of " + String.valueOf(zoomFactor));
	repaint();
    }
    
    
    ////////////////////////////////////////////////////////////////////////////////////////
    public void resetImage () {
	image_width = original_image_width;
	image_height = original_image_height;
	zoom(1f);
    }

    
    ///////////////////////////////////////////////////////////////////////////////////////
    public void selectedRegion (Rectangle r) {
	System.out.println ("Selected range of hits.");
	
	selectedMatches.clear();

	int x1_pix = r.x;
	int y1_pix = r.y;
	int x2_pix = x1_pix + r.width;
	int y2_pix = y1_pix + r.height;

	// determine pixels per bp:
	double pixelsPerBP_width = (double) image_width / max_x_bp_coord;
	double pixelsPerBP_height = (double) image_height / max_y_bp_coord;
	
	// convert pixel coordinates to bp coordinates:
	
	int x1 = (int) (x1_pix / pixelsPerBP_width);
	int x2 = (int) (x2_pix / pixelsPerBP_width);
	int y1 = (int) (y1_pix / pixelsPerBP_height);
	int y2 = (int) (y2_pix / pixelsPerBP_height);
	
	System.out.println ("Selected region is: (" + String.valueOf(x1) + "," 
			    + String.valueOf(y1) + ") , ("
			    + String.valueOf(x2) + "," +
			    String.valueOf(y2)  + ")");
	
	// find matches within selected region:
	
	RTpoint[] hits = rangeTree.query2D(x1,x2,y1,y2); 

	if (hits != null) {
	    System.out.println ("number of hits: " + hits.length);
	    
	    for (int i = 0; i < hits.length; i++) {
		//System.out.println (hits[i]);
		int x = hits[i].get(0);
		int y = hits[i].get(1);
		String matchKey = String.valueOf(x) + "," + String.valueOf(y);
		DatabaseMatch match = (DatabaseMatch) coordToMatch.get(matchKey);
		if (scoreWithinRange(match)) {
		    System.out.println(match);
		    selectedMatches.put(new Integer(match.id()), new Boolean(true));
		}
				
	    }
	} else {
	    System.out.println ("No hits selected.");
	}
	
	repaint();
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////
    public void zoomOnCoordinates (Rectangle r) {

	// determine zoom value:
	// zoom ratio based on viewport size to select band ratio
	Dimension viewDimension = XYplot.this.plotScroller.getViewport().getExtentSize();


	float zoomRatio = (float)viewDimension.width / r.width;
	
	// get upmost coordinates
	int viewX = (int) (zoomRatio * r.x);
	int viewY = (int) (zoomRatio * r.y);
	
	zoom(zoomRatio);
	
	JViewport viewer = XYplot.this.plotScroller.getViewport();
	viewer.setViewPosition(new Point(viewX, viewY));
	
	repaint(); 
    }

    /////////////////////////////////////////////////////////////////////////////////

    public Dimension getPreferredSize() {
	return (new Dimension(image_width+200, image_height+300));
    }


    public void launch_maxLogscoreAdjuster() {
	
	maxLogscore_adjust = new EvalueSlider(this);
	repaint();

    }

    public void disposedEvalueSlider () {
	this.maxLogscore_adjust = null;
	repaint();
    }
    

    private boolean scoreWithinRange (DatabaseMatch m) {
	if (maxLogscore_adjust == null) {
	    return (true);
	} else {
	    double score = m.getScore();
	    
	    double logValue = Math.log(score) / Math.log(10);
	    
	    String acc1 = m.feat1().get_accession();
	    String acc2 = m.feat2().get_accession();
	    
	    int numHits1 = ( (Integer) matchCounter.get(acc1)) . intValue();
	    int numHits2 = ( (Integer) matchCounter.get(acc2)) . intValue();
	    


	    if (logValue <= maxLogScore && (numHits1 <= maxNumHits && numHits2 <= maxNumHits) ) {
		return (true);
	    }
	}
	return (false);
    }
    
    public void setMaxLogScore (float score) {
	this.maxLogScore = score;
	repaint();
    }

    public void setMaxNumHits (int numHits) {
	this.maxNumHits = numHits;
	repaint();
    }
    
    
    //////////////////////////////////////////////////////////////////////////////////////////
    public static void main (String args[]) {
	
	XYplot.TEST = true;
	
	new XYplot(500, "Testing");
	
    }

}
