import java.awt.*;
import javax.swing.*;
import java.lang.Math;

/* Adapted from: Rule.java is used by ScrollDemo.java. */

public class SequenceTicker extends JComponent {

    private boolean SEE = false;
    
    public static final int HORIZONTAL = 0;
    public static final int VERTICAL = 1;
    public static final int SIZE = 100;
    
    private XYplot adjacentPanel;
    private int max_seq_length = 0;
    public int orientation;
    
    private int increment = 10; //pixels for small tick line.
    private int units = 100; // pixels for big tick line
        
    private String axisLabel = "label";

    //constructor
    public SequenceTicker (XYplot panel, int orient) {
        orientation = orient;
	adjacentPanel = panel;
    }
    
    
    public void setAxisLabel (String label) {
	axisLabel = label;
    }


    public void set_max_seq_length (int l) {
	max_seq_length = l;
    }
    
    public void setPreferredHeight(int height) {
        setPreferredSize(new Dimension(SIZE, height));
    }

    public void setPreferredWidth(int width) {
        setPreferredSize(new Dimension(width, SIZE));
    }
    
    
    public void paintComponent(Graphics g) {
	
	float pixels_per_bp;

	if (max_seq_length == 0) { 
	    return;
	}
	
	float bp_per_pixel = 0;
	Dimension plotDimension = adjacentPanel.getPlotDimension();

	if (orientation == HORIZONTAL) {
	    int width = plotDimension.width;
	    if (width > 0) {
		bp_per_pixel = (float) max_seq_length/width;
	    }
	} else {
	    int height = plotDimension.height;
	    if (height > 0) {
		bp_per_pixel = (float) max_seq_length / height;
	    }
	}

	if (SEE) {
	    System.out.println ("bp per pixel: " + String.valueOf(bp_per_pixel));
	}

	if (bp_per_pixel == 0) {
	    return;
	}
	
	Rectangle drawHere = g.getClipBounds();
	
        // Fill clipping area with dirty brown/orange.
        g.setColor(new Color(174,114,228));
        g.fillRect(drawHere.x, drawHere.y, drawHere.width, drawHere.height);
	
	if (SEE) {
	    System.out.println ("drawHere: " + String.valueOf(drawHere.x) + ", " + String.valueOf(drawHere.y)
			    + ", " + String.valueOf(drawHere.width) + ", " + String.valueOf(drawHere.height) 
				+ " type: " + String.valueOf(this.orientation));
	}

	
        // Do the ruler labels in a small font that's black.
        g.setFont(new Font("SansSerif", Font.PLAIN, 10));
        g.setColor(Color.white);

        // Some vars we need.
        int end = 0;
        int start = 0;
        int tickLength = 0;
        String text = null;
	
        // Use clipping bounds to calculate first and last tick locations.
        if (orientation == HORIZONTAL) {
            
	    start = drawHere.x;
	    end = drawHere.x + drawHere.width;
	    if (end > plotDimension.width) {
		end = plotDimension.width;
	    }

	    
	} else {
            start = drawHere.y;
            end = drawHere.y + drawHere.height;
	    if (end > plotDimension.height) {
		end = plotDimension.height;
	    }
	}
	
	int midPt = (start+end)/2;
	// Make a special case of 0 to display the number
        // within the rule and draw a units label.
	tickLength = 10;
	if (orientation == HORIZONTAL) {
	    g.drawLine(0, SIZE-1, 0, SIZE-tickLength-1); 
	    g.drawString(axisLabel, midPt, 21); //2
	} else {
	    g.drawLine(SIZE-1, 0, SIZE-tickLength-1, 0);
	    g.drawString(axisLabel, 9, midPt); //10
	}
	text = null;
		
        // ticks and labels
	// walk up to increment distance.
	int startAdj = start % increment;
	if (startAdj != 0) {
	    start += increment - startAdj;
	}
	
        for (int i = start; i < end; i += increment) {
            if (i % units == 0)  {
                tickLength = 10;
		// calculate bp position based on pixel coordinate.
		int coordPosition = (int) (i * bp_per_pixel);
		text = String.valueOf(coordPosition);

		if (SEE) {
		    System.out.println ("Label: (bp_per_pixel: " + String.valueOf(bp_per_pixel) + "), (i = " + String.valueOf(i) 
					+ "), " + text);
		}
		
		
	    } else {
                tickLength = 7;
                text = null;
            }
	    
            if (tickLength != 0) {
                if (orientation == HORIZONTAL) {
                    g.drawLine(i, SIZE-1, i, SIZE-tickLength-1); //y is constant
                    if (text != null)
                        g.drawString(text, i-3, 75);
                } else {
                    g.drawLine(SIZE-1, i, SIZE-tickLength-1, i); // x is constant
                    if (text != null)
                        g.drawString(text, 30, i+3);
                }
            }
        }
    }
}


