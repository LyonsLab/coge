import java.awt.*;
import javax.swing.*;
import java.util.*;
import java.awt.event.*;
import javax.swing.event.*;

public class EvalueSlider extends JFrame 
    implements ChangeListener  {
    
    private XYplot parentPlot = null;
    private JSlider slider = null;
    private JSlider hitNumslider = null;

    public EvalueSlider (XYplot xyPlot) {
	
	this.parentPlot = xyPlot;
	// init settings.
	xyPlot.setMaxLogScore (-5);
	xyPlot.setMaxNumHits(20);
	
	
	addWindowListener(new WindowAdapter() {
		public void windowClosing(WindowEvent e) {
		    EvalueSlider.this.parentPlot.disposedEvalueSlider(); // callback.
		    dispose();
		}	
	    });
	
	//Create the slider and its label
        JLabel sliderLabel = new JLabel("Maximum Log E-value", JLabel.CENTER);
        sliderLabel.setAlignmentX(Component.CENTER_ALIGNMENT);

        slider = new JSlider(JSlider.HORIZONTAL,-100, 1, -5);
        slider.addChangeListener(this);
	
        //Turn on labels at major tick marks.
        slider.setMajorTickSpacing(20);
        slider.setMinorTickSpacing(5);
        slider.setPaintTicks(true);
        slider.setPaintLabels(true);
        slider.setBorder(
	BorderFactory.createEmptyBorder(0,0,10,0));


	// add slider for max number of hits adjustment
	JLabel hitNumsliderLabel = new JLabel("Max # of Matches", JLabel.CENTER);
        hitNumsliderLabel.setAlignmentX(Component.CENTER_ALIGNMENT);

        hitNumslider = new JSlider(JSlider.HORIZONTAL,1, 20, 20);
        hitNumslider.addChangeListener(this);
	
        //Turn on labels at major tick marks.
        hitNumslider.setMajorTickSpacing(5);
        hitNumslider.setMinorTickSpacing(1);
        hitNumslider.setPaintTicks(true);
        hitNumslider.setPaintLabels(true);
        hitNumslider.setBorder(
	BorderFactory.createEmptyBorder(0,0,10,0));


	//Put everything in the content pane.
        Container contentPane = getContentPane();
        contentPane.setLayout(new BoxLayout(contentPane, BoxLayout.Y_AXIS));
        contentPane.add(sliderLabel);
        contentPane.add(slider);
	contentPane.add(hitNumsliderLabel);
	contentPane.add(hitNumslider);
		
	pack();
	setVisible(true);
    }
    
    public void stateChanged(ChangeEvent e) {
	JSlider source = (JSlider)e.getSource();
	if (source == slider) {
	    if (!source.getValueIsAdjusting()) {
		int logVal = (int)source.getValue();
		System.out.println ("Log Evalue set to: " + String.valueOf(logVal));
		this.parentPlot.setMaxLogScore(logVal);
	    }
	} else if (source == hitNumslider) {
	    if (!source.getValueIsAdjusting()) {
		int maxNumHits = (int)source.getValue();
		System.out.println ("Max # hits set to: " + String.valueOf(maxNumHits));
		this.parentPlot.setMaxNumHits(maxNumHits);
	    }
	}
    }
    
    

}    


