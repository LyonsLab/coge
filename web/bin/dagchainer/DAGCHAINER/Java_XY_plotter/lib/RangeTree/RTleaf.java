package RangeTree;

// A part of the 2d Orth Search.



class RTleaf 
{
	private Object Data=null;  //this will hold the data (Tree object or Point)
	private RTleaf lc=null,rc=null; //left and right children
	private int el=0;   //this will help in findSplitRTleaf and Query1D/2D
						// Will hold the X or Y coordinate of the object in the tree.
	
	//constructors
	
	public RTleaf() 	{;} 
	
	public RTleaf(Object d) 
	{
		this.setData(d);
	}
	
	// functions
	
	public void setEl(int e) //set apropriate coordinate dipending on the tree type.
	{
		this.el=e;
	}//setElemnt
	
	public int getEl() 
	{
		return this.el;
	}//getEl
	
	public boolean childless() //checks if left and right children doesnt not exist. 
	{
		return ((this.lc==null)&&(this.rc==null));
	}
	
	public void setData(Object d) 
	{
		this.Data=d;
	}
	
	public Object getData()
	{
		return this.Data;
	}// getData
	
	public RTleaf getLc()
	{
		return this.lc;
	}//Left Child
	
	public RTleaf getRc()
	{
		return this.rc;
	}//right child
	
	public void setLeft(RTleaf O)
	{
		this.lc=O;
	}
	
	public void setRight(RTleaf O)
	{
		this.rc=O;
	}
}//class RTleaf
	
