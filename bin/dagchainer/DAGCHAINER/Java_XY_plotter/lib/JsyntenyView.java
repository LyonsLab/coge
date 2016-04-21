import java.io.*;
import java.util.*;
import java.awt.*;
import javax.swing.*;
import java.awt.event.*;

public class JsyntenyView 
    implements Comparator {

    private Hashtable matchMap; // holds match pairs
    
    private Hashtable allMatches;

    private Hashtable chromoKeyMap;
    
    // constructor
    public JsyntenyView () {
	
	allMatches = new Hashtable(); /* key each match pair using a key constructed from:
                                     chromoA + accA + chromoB + accB
                                  */
	
	matchMap = new Hashtable(); // keys a vector of match pairs to a chromosome pair.
    
	chromoKeyMap = new Hashtable(); /* provides chromsome pair identity based on key.
                                       ...could just parse the hash key to get it, but 
                                       the paired info is cleaner. */
    
    }
    
    
    public static void main (String args[]) {
        
        int numArgs = args.length;
        if (numArgs < 1) {
            System.err.println("usage: matchFile [diagonalFile]");
            System.exit(1);
        }
        
        JsyntenyView viewer = new JsyntenyView();
        
        String matchFile = args[0];
        
        viewer.parseMatchFile(matchFile);
        
        
        
        // see if diagonalFile is available
        if (numArgs > 1) {
            String diagonalFile = args[1];
            viewer.parseDiagonalFile(diagonalFile);
            
            // viewer.describeDiagonalList();
        }
        
        
        final JFrame appWindow = new JFrame("testing XYplot");
        
        appWindow.addWindowListener(new WindowAdapter() {
                public void windowClosing(WindowEvent e) {
                    appWindow.dispose();
                    System.exit(0);
                }	
            });
        
        final JDesktopPane desktop = new JDesktopPane();
        
        
        
        appWindow.getContentPane().add(desktop);
        
        
        // get screen dimensions:
        Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
        int appWidth = (int) (0.6f * screenSize.width);
        int appHeight = (int) (0.8f * screenSize.height);
        
        int maxPlotDimension = (int) (0.75f * ((appWidth+appHeight)/2.0f));
        
        appWindow.setSize(appWidth, appHeight);
        
        
        Enumeration chromoPairs = viewer.matchMap.keys();
        
        while (chromoPairs.hasMoreElements()) {
            
            String chromoPairKey = (String) chromoPairs.nextElement();
            
            //System.out.println ("Sorting matches between " + chromoPairKey);
            System.out.println ("Launching " + chromoPairKey);
            
            
            StringPair chromos = (StringPair) viewer.chromoKeyMap.get(chromoPairKey);
            
            
            
            Vector matchList = (Vector) viewer.matchMap.get(chromoPairKey);
            
            // Collections.sort (matchList, viewer);
            
            //viewer.describeMatchList();
            
            
            XYplot genomePlot = new XYplot(maxPlotDimension, chromoPairKey);
            
            genomePlot.setXYaxisLabels("chr" + chromos.b, "chr" + chromos.a);
            
            genomePlot.setMatches(matchList);
            
            genomePlot.setVisible(true);
            
            desktop.add(genomePlot);
            
            
        }
        
        appWindow.setVisible(true);
    }
    
    
    
    private void parseMatchFile (String filename) {
        
        try {
            BufferedReader bufferedFileReader = new BufferedReader (new FileReader(filename));
            
            
            try {
                String inLine;
                while ((inLine = bufferedFileReader.readLine()) != null) {
                    // System.out.println(inLine);
                    
                    StringTokenizer stringParser = new StringTokenizer(inLine, "\t", false);
                    
                    int numTokens = stringParser.countTokens();
                    if (numTokens < 8) {
                        System.err.println("Line: " + inLine + " lacks required data\n");
                        continue;
                    }
                    String chromo1 = stringParser.nextToken();
                    String acc1 = stringParser.nextToken();
                    long end5_1 = Long.parseLong(stringParser.nextToken());
                    long end3_1 = Long.parseLong(stringParser.nextToken());
                    
                    SequenceFeature feature1 = new SequenceFeature(acc1, end5_1, end3_1);
                    feature1.set_moleculeName(chromo1);
                    
                    
                    String chromo2 = stringParser.nextToken();
                    String acc2 = stringParser.nextToken();
                    long end5_2 = Long.parseLong(stringParser.nextToken());
                    long end3_2 = Long.parseLong(stringParser.nextToken());
                    
                    SequenceFeature feature2 = new SequenceFeature (acc2, end5_2, end3_2);
                    feature2.set_moleculeName(chromo2);
                    
                    String chromoA = chromo1;
                    String chromoB = chromo2;
                    String accA = acc1;
                    String accB = acc2;
                    
                    
                    double score = -1;
                    if (numTokens > 8) {
                        // parse score
                        score = Double.parseDouble(stringParser.nextToken());
                    }
                    
                    if (chromo1.compareTo(chromo2) > 0 || (chromo1.compareTo(chromo2) == 0 && end5_1 > end5_2) ) {
                        // swap the features:
                        SequenceFeature tempFeat = feature1;
                        feature1 = feature2;
                        feature2 = tempFeat;
                        
                        chromoA = chromo2;
                        chromoB = chromo1;
                        accA = acc2;
                        accB = acc1;
                    }
                    
                    
                    String matchKey = chromoA + "$" +  accA + "$"
                        + chromoB + "$" + accB;
                    
                    DatabaseMatch matchPair = (DatabaseMatch) allMatches.get(matchKey);
                    
                    
                    if (matchPair == null) {
                        
                        matchPair = new DatabaseMatch (feature1, feature2);
                        
                        Vector matchList = get_match_list(chromo1, chromo2);
                        
                        matchList.add(matchPair);
                        allMatches.put(matchKey, matchPair);
                        
                        // duplicate entry for chromosme self search:
                        if (chromoA.equals(chromoB)) {
                            String matchKey2 = chromoB + "$" +  accB + "$"
                                + chromoA + "$" + accA;
                            allMatches.put(matchKey2, matchPair);
                        }			
                    }
                    if (score > -1) {
                        matchPair.setScore(score);
                    }
                    
                }
                bufferedFileReader.close();
            }
            catch (IOException ioe) {
                System.err.println(ioe.toString());
            }
            
        }
        catch (FileNotFoundException fnfe) {
            
            System.err.println("Can't find file " + filename);
            System.exit(1);
        }
    }
    
    
    private void parseDiagonalFile (String filename) {
        
        try {
            BufferedReader bufferedFileReader = new BufferedReader (new FileReader(filename));
            
            
            try {
                
                String inLine;
                
                int currDiagonal = 0;
                
                Color diagColor = Color.red; // init 
                
                while ((inLine = bufferedFileReader.readLine()) != null) {
                    // System.out.println(inLine);
                    
                    StringTokenizer stringParser = new StringTokenizer(inLine, "\t", false);
                    
                    int numTokens = stringParser.countTokens();
                    if (numTokens < 8) {
                        // diagonal separator
                        currDiagonal++;
                        diagColor = getRandomColor();
                        System.out.println("Reading diagonal: " + String.valueOf(currDiagonal));
                        continue;
                    }
                    
                    String chromo1 = stringParser.nextToken();
                    String acc1 = stringParser.nextToken();
                    long end5_1 = Long.parseLong(stringParser.nextToken());
                    long end3_1 = Long.parseLong(stringParser.nextToken());
                    
                    
                    String chromo2 = stringParser.nextToken();
                    String acc2 = stringParser.nextToken();
                    long end5_2 = Long.parseLong(stringParser.nextToken());
                    long end3_2 = Long.parseLong(stringParser.nextToken());
                    
                    
                    
                    
                    String chromoA = chromo1;
                    String chromoB = chromo2;
                    String accA = acc1;
                    String accB = acc2;
                    
                    if (chromo1.compareTo(chromo2) > 0 || (chromo1.compareTo(chromo2) == 0 && end5_1 > end5_2) ) {
                                                
                        chromoA = chromo2;
                        chromoB = chromo1;
                        accA = acc2;
                        accB = acc1;
                        
                    }
                    
                    String matchKey = chromoA + "$" +  accA + "$"
                        + chromoB + "$" + accB;
                    DatabaseMatch matchPair = (DatabaseMatch) allMatches.get(matchKey);
                    
                    if (matchPair != null) {
                        if (matchPair.inDiagonal()) {
                            System.out.println ("Match already in diagonal: " + matchPair);
                        } else {
                            matchPair.setDiagonal(currDiagonal);
                            matchPair.setColor(diagColor);
                        }
                    } else {
                        System.err.println ("Error, can't find match entry for: " + matchKey);
                    }
                    
                }
                bufferedFileReader.close();
            }
            catch (IOException ioe) {
                System.err.println(ioe.toString());
            }
            
        }
        catch (FileNotFoundException fnfe) {
            
            System.err.println("Can't find file " + filename);
            System.exit(1);
        }
    }
    
    
    private void describeMatchList () {
        
        Enumeration chromoPairs = this.matchMap.keys();
        
        while (chromoPairs.hasMoreElements()) {
            
            String chromoPairKey = (String) chromoPairs.nextElement();
            
            Vector matchList = (Vector) this.matchMap.get(chromoPairKey);
            
            for (int i=0; i < matchList.size(); i++) {
                
                System.out.println (chromoPairKey + "  " + matchList.elementAt(i).toString());
            }
        }
        
    }
    
    
    
    public Vector get_match_list (String c1, String c2) {
        
        String key = get_chromo_pair_key (c1,c2);
        
        Vector matchList = (Vector) matchMap.get(key);
        if (matchList == null) {
            matchList = new Vector();
            matchMap.put(key, matchList);
        }
        return (matchList);
    }
    
    
    
    public String get_chromo_pair_key (String c1, String c2) {
        
        String key, a, b;
        if (c1.compareTo(c2) <= 0) {
            a = c1; b = c2;
        } else {
            a = c2; b = c1;
        }
        
        key = a + "_||_" + b;
        
        chromoKeyMap.put(key, new StringPair(a,b));
        
        return (key);
    }
    
    
    
    public int compare (Object o1, Object o2) {
        
        DatabaseMatch m1 = (DatabaseMatch) o1;
        DatabaseMatch m2 = (DatabaseMatch) o2;
        
        int midY1 = (int) (m1.feat1().midPt());
        int midX1 = (int) (m1.feat2().midPt());
        
        int midY2 = (int) (m2.feat1().midPt());
        int midX2 = (int) (m2.feat2().midPt());
        
        
        // Sorting by rows, then columns.
        if (midY1 < midY2) {
            return (-1);
        } else if (midY1 > midY2) {
            return (1);
        } else if (midY1 == midY2) {
            // must examine X values:
            if (midX1 < midX2) {
                return (-1);
            } else if (midX1 > midX2) {
                return (1);
            }
        }
        
        // default:
        // same values for both.
        // midY1 == midY2 && midX1 == midX2
        
        return (0);
    }
    
    
    
    private class StringPair {
        
        public String a,b;
        
        public StringPair (String a, String b) {
            this.a = a;
            this.b = b;
        }
    }
    
    
    
    public static Color getRandomColor () {
        // choose a random color, cache it.
        
        int[] rgb = {0,0,0};
        
        rgb[0] = (int) (Math.random() * 256);
        rgb[1] = (int) (Math.random() * 256);
        rgb[2] = (int) (Math.random() * 256);
        
        
        
        // chose two randomly, set to ff then 00
        int i = (int) (Math.random() * 3);
        int j = i;
        while ( j != i) {
            j = (int) (Math.random() * 3);
        }
        rgb[i] = 255;
        rgb[j] = 0;
        
        
        
        Color randColor = new Color (rgb[0], rgb[1], rgb[2]);
        return (randColor);
    }
    
}

