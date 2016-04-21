
public class SequenceFeature {
    
    private String accession;
    private long end5;
    private long end3;
    
    private String organismName = null;
    private String moleculeName = null;


    // constructor
    public SequenceFeature (String accession, long end5, long end3) {
	
	this.accession = accession;
	this.end5 = end5;
	this.end3 = end3;
	
    }
    
    
    public void set_organismName (String name) {
	this.organismName = name;
    }

    public String get_organismName () {
	return (this.organismName);
    }

    public void set_moleculeName (String name) {
	this.moleculeName = name;
    }

    public String get_moleculeName () {
	return (this.moleculeName);
    }


    public String get_accession () {
	return (this.accession);
    }

    public long end5() {
	return (this.end5);
    }

    public long end3() {
	return (this.end3);
    }
    
    public long midPt() {
	return ( (this.end5 + this.end3)/2);
    }



    public String toString() {
	
	String text = moleculeName + "\t" + accession + "\t" + String.valueOf(end5) + "\t" + String.valueOf(end3); 
	
	return (text);
    }
}



	
