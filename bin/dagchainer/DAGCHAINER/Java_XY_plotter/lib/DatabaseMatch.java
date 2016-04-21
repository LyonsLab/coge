import java.util.*;
import java.text.*;
import java.awt.*;

public class DatabaseMatch {
    
    private static int matchCounter = 0;

    private SequenceFeature feat1;
    private SequenceFeature feat2;
    
    private int diagonal; // set to diagonal number; 0 if not in diagonal.
    private int id;

    private double score = -1; //default.
    DecimalFormat df;
    
    private Color diagColor;

    // include other parameters such as percent ID, E-value, etc.


    public DatabaseMatch (SequenceFeature feat1, SequenceFeature feat2) {
	
	this.feat1 = feat1;
	this.feat2 = feat2;
	diagonal = 0;
	matchCounter++;
	this.id = matchCounter;
	df = new DecimalFormat("0.0#E0");
    }
    
    public void setDiagonal (int d) {
	diagonal = d;
    }

    
    public boolean inDiagonal () {
	if (diagonal != 0) {
	    return (true);
	} else {
	    return (false);
	}
    }
    
    public int diagonal () {
	return (this.diagonal);
    }

    
    public SequenceFeature feat1 () {
	return (feat1);
    }

    public SequenceFeature feat2 () {
	return (feat2);
    }
    
    public int id () {
	return (this.id);
    }
    

    public String toString() {
	String retText = feat1.toString() + "\t" + feat2.toString();
	if (score > -1) {
	    String scoreString = df.format(this.score); //String.valueOf(position);

	    retText += "\t" + scoreString;
	}
	return (retText);
    }
    
    
    public void setScore (double s) {
	this.score = s;
    }
    
    public double getScore () {
	return (this.score);
    }
    
    
    public void setColor (Color c) {
	this.diagColor = c;
    }
    
    public Color getColor () {
	return (diagColor);
    }
    
}


