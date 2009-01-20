package RangeTree;

// A part of the 2D Orth. Search


public class RangeTree 
{
	RTleaf root=null;
	int el=0;
	
	//constructors
	
	public RangeTree() {;}
	
	public RangeTree(RTpoint[] A,int dim)
	{
		if (dim==0) root=build(A,0,A.length-1);  else root=build1D(A,0,A.length-1);
		root.setEl(A[A.length/2].get(dim));
	}
	
	public RangeTree(RTpoint[] A)
	{
		this(A,0);
	}
	
	//functions
	
	public boolean isEmpty() 
	{
		return (root==null);
	}// isEmpry
	
	public RTleaf getRoot() 
	{
		return this.root;
	}
		
	public int getEl() 
	{
		return root.getEl();
	}
		
	public RTleaf build1D(RTpoint[] A,int start,int end)    //this will build a 1D tree.
	{													// Assuming the A is sorted by Y.
		int mid=(start+end)/2;
		RTleaf tmp=new RTleaf(A[mid]);
		tmp.setEl(A[mid].get(1));// this will set the elemnt of the leaf(Y is the apropriate)
	
		if (start!=end) 
		{
			tmp.setLeft(build1D(A,start,mid));
			tmp.setRight(build1D(A,mid+1,end));
		}//if start!=end
		return tmp;
	}//build1D
	
	public RTleaf build(RTpoint[] A,int start,int end)	//this will build a 2D tree, assuming A 
	{												// is sorted by X.
		int mid=(start+end)/2;
		RTpoint[] B=new RTpoint[end-start+1];// build an array for the assoc. tree.
		
		for (int i=start;i<=end;i++) B[i-start]=(RTpoint)A[i].clone();
		mergeSort(B,1);//sorting by Y.
			
		RTleaf tmp=new RTleaf(new RangeTree(B,1));//creating an associated tree by build1D.
		tmp.setEl(A[mid].get(0));//set elemnt (apropriate is X coordinate.)
		
		if (start!=end) 
		{
			tmp.setLeft(build(A,start,mid));//seting children
			tmp.setRight(build(A,mid+1,end));
		} 
		
		return tmp;
	}// build 2D
	
	private RTleaf findSplitRTleaf(int low,int high){ //finding split leaf, in this tree.
		RTleaf v=this.root;
		if(v!=null)
		{
			while (!v.childless()&&(v.getEl()<low||v.getEl()>high))
			{
				if(v.getEl()<low) v=v.getRc(); else v=v.getLc();
			}//while
		}//if
		
		return v;
	}//findSplitLeaF

	public RTpoint[] query2D(int xlow, int xhigh, int ylow, int yhigh)
	{  //this will query in 2D range
		RTpoint[] ans=null; // and array for points found.

		RTleaf v=this.findSplitRTleaf(xlow,xhigh); //to start from?
		if (v.childless()){
			RTleaf l=(RTleaf)(v.getData());
			RangeTree t=(RangeTree)(l.getData());
			if ((v.getEl()>=xlow)&&(v.getEl()<=xhigh)) ans=concat(ans,t.query1D(ylow,yhigh));
										//^^ found a childless point? > check in Y dim.
		} else
		{
			RangeTree t;
			RTleaf l=v.getLc();// check left side of split leaf
			while(!l.childless()){
				if(l.getEl()>=xlow)  //if left is good, then all the right side is good too.
				{
					t=(RangeTree)(l.getRc().getData());
					ans=concat(ans,t.query1D(ylow,yhigh));//now check Y dim.
					l=l.getLc();
				} else l=l.getRc();//if left is not good, go alittle right.
			}//while
			t=(RangeTree)(l.getData());
			if ((l.getEl()>=xlow)&&(l.getEl()<=xhigh)) ans=concat(ans,t.query1D(ylow,yhigh));
								//check last leaf.
			
			l=v.getRc(); //now check in the same method, the right side.
			while(!l.childless()){
				if(l.getEl()<=xhigh)
				{
					t=(RangeTree)(l.getLc().getData());
					ans=concat(ans,t.query1D(ylow,yhigh));
					l=l.getRc();
				} else l=l.getLc();
			}//while
			t=(RangeTree)(l.getData());
			if ((l.getEl()>=xlow)&&(l.getEl()<=xhigh)) ans=concat(ans,t.query1D(ylow,yhigh));
		}//else
		
		if (ans!=null) return checkRTpoints(ans,xlow,xhigh,ylow,yhigh);//will return only good points
		else return null;
	}//query2D
	
	private RTpoint[] checkRTpoints(RTpoint[] A,int xlow, int xhigh, int ylow, int yhigh)
	{  //this will check each point in the array and return only good points array.
		int count=0;
		for (int i=0;i<A.length;i++) if ((A[i].get(0)>=xlow)&&(A[i].get(0)<=xhigh)&&
										(A[i].get(1)>=ylow)&&(A[i].get(1)<=yhigh)) count++;
		RTpoint[] tmp=new RTpoint[count];
		count=0;
		
		for (int i=0;i<A.length;i++) if ((A[i].get(0)>=xlow)&&(A[i].get(0)<=xhigh)&&
										(A[i].get(1)>=ylow)&&(A[i].get(1)<=yhigh)) 
										{
											tmp[count]=(RTpoint)A[i].clone();
											count++;
										}
		return tmp;
	}//checkRTpoints
	
	private RTpoint[] query1D(int ylow,int yhigh) //this will query in Y dim.
	{
		RTpoint[] tmp=null;
		RTleaf v=this.findSplitRTleaf(ylow,yhigh);
		
		if (v.childless()&&(v.getEl()<=yhigh)&&(v.getEl()>=ylow))
		{				//^^ a childless leaf that is good?
			tmp=new RTpoint[1];
			RTpoint p=(RTpoint)(v.getData());
			tmp[0]=(RTpoint)p.clone();
			return tmp;
		}
		
		if (!v.childless()) 
		{  //same algorithm as explained in query2D function.
			RTleaf l=v.getLc();
			
			while(!l.childless())
			{
				if (l.getEl()>=ylow)
				{
					tmp=concat(tmp,reportSubRangeTree(l.getRc()));					
	  				l=l.getLc();
				}//if
				else l=l.getRc();
			}//while
			if (l.getEl()>=ylow) tmp=concat(tmp,reportSubRangeTree(l));

			l=v.getRc();	
			while(!l.childless())
			{
				if(l.getEl()<=yhigh)
				{
					tmp=concat(tmp,reportSubRangeTree(l.getLc()));					
					l=l.getRc();
				}//if
				else l=l.getLc();
			}//while
 	  		if (l.getEl()<=yhigh) tmp=concat(tmp,reportSubRangeTree(l));
		} // if !childless	
		return tmp;
	}//query1D
	
	private static RTpoint[] reportSubRangeTree(RTleaf l)
	{   		//will report all the leaf below the RTleaf 'l'.
		RTpoint[] tmp=null;

		if (l.childless())
		{
			tmp=new RTpoint[1];
			tmp[0]=(RTpoint)((RTpoint)(l.getData())).clone();
		} else 
		{
			if (l.getLc()!=null) tmp=concat(tmp,reportSubRangeTree(l.getLc()));
			if (l.getRc()!=null) tmp=concat(tmp,reportSubRangeTree(l.getRc()));
		}
		
		return tmp;
	}// reportSubRangeTree

	public static void mergeSort(RTpoint[] A,int dim) 
	{			//regular mergeSort algorithm
				//added as a static method, so no additional class is needed.
		ms_divide(A,0,A.length/2,dim);
		ms_divide(A,(A.length/2)+1,A.length-1,dim);
		ms_conq(A,0,A.length/2,A.length-1,dim);
	} //mergeSort
	
	private static void ms_divide(RTpoint[] A, int start,int end,int dim)
	{
		if (start<end) 
		{
			ms_divide(A,start,(start+end)/2,dim);
			ms_divide(A,((start+end)/2)+1,end,dim);
			ms_conq(A,start,(start+end)/2,end,dim);
		}
	} // ms_divide;
	
	private static void ms_conq(RTpoint[] A,int start,int mid,int end,int dim)
	{
		RTpoint[] tmp=new RTpoint[end-start+1];
		int a=start;
		int b=mid+1;
		int c=0;
		for (int i=0;i<tmp.length;i++) 
		{
			if ((b>end)||((a<=mid)&&(A[a].get(dim)<A[b].get(dim)))) 
			{
				tmp[i]=(RTpoint)A[a].clone();
				a++;
			} else 
			{
				tmp[i]=(RTpoint)A[b].clone();
				b++;
			}
		}//for
		
		for (int i=0;i<tmp.length;i++) A[i+start]=(RTpoint)tmp[i].clone();
	}//ms_conq	
				
	private static RTpoint[] concat(RTpoint[] A,RTpoint[] B) 
	{ 			// this will merge two RTpoint arrays to one.
		if ((A==null)&&(B==null)) return null;
		if (A==null) return B;
		if (B==null) return A;
		
		RTpoint[] tmp=new RTpoint[A.length+B.length];
		for (int i=0;i<A.length;i++) tmp[i]=(RTpoint)A[i].clone();
		for (int i=0;i<B.length;i++) tmp[i+A.length]=(RTpoint)B[i].clone();
		
		return tmp;
	}//concat
	
}//class RangeTree
