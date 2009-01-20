package RangeTree;

/**
* this class represent a 2d point in the plane. <br>
*/
public class RTpoint
{
// ***** private data *****
	private int[] p=new int[2];
	
// ****** constructors ******
   public RTpoint (int x1, int y1)
   {
      p[0] = x1;
      p[1] = y1;
   }
 /** copy constructor */
   public RTpoint (RTpoint pp)
   {
      this.p[0] = pp.get(0);
      this.p[1] = pp.get(1);
   }
 
   // ***** public methodes *****
   public int get(int i)
   {return this.p[i];}

   public void set(int i,int x) 
   { p[i]= x;}

   public boolean equals (RTpoint pp)
   {
   	boolean eq=true;
   	int i=0;
   	while ((i<pp.length()) && eq) eq=(this.p[i]==pp.get(i));

   	return eq;
   }
   public String toString()
   {
		return this.p[0]+" "+this.p[1];
   }
   
   public int length()
   {
   		return this.p.length;
   }
 
 	public Object clone()
 	{
		return new RTpoint(this);
 	}  
}// class RTpoint
