import java.util.*;

public class MatchDiagonal {

    private Vector matchList;

    public MatchDiagonal () {
	
	matchList = new Vector();
    }


    public void addDatabaseMatch (DatabaseMatch match) {
	
	this.matchList.add (match);
	
    }


    public final Vector getMatchList () {
	
	return (matchList);

    }

    public String toString() {
	
	String out = "";

	for (int i = 0; i < matchList.size(); i++) {
	    out += matchList.elementAt(i).toString();
	}

	return (out);
    }
}

