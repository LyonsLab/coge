
                 /*************************\
                 *                         *
                 *     DIALIGN 2           *
                 *                         *
                 *     functions.c         *
                 *                         *
                 \*************************/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "define.h"
#include "dialign.h"
#include "alig_graph_closure.h"


extern int iter_cond_prob , col_score  ;
extern short **cont_it_p; 
extern char input_name[NAME_LEN];  
extern int anchors ; 
extern int frg_mult_file ;
extern int frg_mult_file_v ;
extern short crick_strand ;
/* extern char dia_pa_name[NAME_LEN]; */
extern int pr_av_max_nd , wgt_type ;
extern int mask ;
extern char prn[ NAME_LEN ];
extern char input_line[ NAME_LEN ];
extern int print_status;
extern char pst_name[NAME_LEN];
extern char clust_sim[ NAME_LEN ] ;
extern FILE *fp_dia, *fp_dpa;
extern float tot_weight, av_len ;
extern int dia_num, max_dia_num, msf_file, cw_file; 
extern int istep , anc_num;
extern int fasta_file;
extern char *upg_str;
extern int plot_num;
extern int argnum ;
extern char *seq[MAX_SEQNUM];
extern int  *seqlen;
extern char *seq_name[MAX_SEQNUM];
extern char printname[NAME_LEN];
extern char amino_acid[22];
extern int par_count;
extern int *num_dia_bf;  
extern int *num_dia_af;  

extern float pairalignsum;    
extern int pairalignlen ;  

extern struct multi_frag *this_it_dia;
extern struct multi_frag *all_it_dia;

extern CLOSURE *clos;
extern float **glob_sim;

extern int ***open_pos;

extern float **wgt_prot, **wgt_dna, **wgt_trans ;

extern int sim_score[21][21];

extern int min_dia, max_dia;
extern int long_output ; 
extern int seqnum;
extern short dots;
extern float  threshold; 
extern int num_all_it_dia;

extern int num_dia_p, overlap_weights ;

extern int **amino;
extern int **amino_c;


extern int **shift; 
extern double **tp400_prot; 
extern double **tp400_dna; 
extern double **tp400_trans; 


int num_test( char *cp )
  { 
    int result = 1;
    int i;
    char *strng;
    
    strng = cp;
    
    for(i = 0 ; i < strlen(strng) ; i++ )
      if( ! isdigit(strng[i]) && ( strng[i] != '.' ) )
        {
          result = 0;
     /*     printf("\n %c is no digit !!!\n", strng[i]);   */
        } 
         
    return result ;    
  }



float minf2(float a, float b)
{
 if (a<b)
    return a;
 else return b;
}

float maxf2(float a, float b)
 {
  if (a>b)
     return a;
  else return b;
 }


int mini2(int a, int b)
  {
    if(a<b)
      return a;
    else 
      return b;
  }




int maxi2(int a, int b)
 {
  if (a>b)
     return a;
  else return b;
 }

int mini3(int a, int b, int c)
 {
  return mini2(a, mini2(b,c));
 }

void minf(float *a, float b)
 {
  if (*a > b)
     *a = b;
 }

void mini(int *a, int b)
 {
  if (*a > b)
     *a = b;
 }

void maxf(float *a, float b)
 {
  if (*a < b)
     *a = b;
 }

void maxi(int *a, int b)
 {
  if (*a < b)
     *a = b;
 }


char invert ( char c1 )
{
  char c2 = c1 ;

  if(c1 == 'T')
    c2 = 'A' ; 
  if(c1 == 'A')
    c2 = 'T' ;
  if(c1 == 'C')
    c2 = 'G' ;
  if(c1 == 'G')
    c2 = 'C' ;
  
  return( c2 );

}

int translate(char c1, char c2 ,char c3, int seqno, int pos)
{
 /* translation of triplets into amino acids */
  

 int  amac;          /* resulting amino acid */

 if(c1 == 'T')
   {
    if(c2 == 'T')
      {
       if(c3 == 'T') amac = 18;
       if(c3 == 'C') amac = 18;
       if(c3 == 'A') amac = 16;
       if(c3 == 'G') amac = 16;
      }
    if(c2 == 'C') amac = 2;
    if(c2 == 'A')
      {
       if(c3 == 'T') amac = 19;
       if(c3 == 'C') amac = 19;
       if(c3 == 'A') amac = 0;     /* stop codon */
       if(c3 == 'G') amac = 0;
      }
    if(c2 == 'G')
      {
       if(c3 == 'T') amac = 1;
       if(c3 == 'C') amac = 1;
         if(c3 == 'A') amac = 20;
      
       if(c3 == 'G') amac = 20;
      }
   }

 if(c1 == 'C')
   {
    if(c2 == 'T') amac = 16;
    if(c2 == 'C') amac = 4;
    if(c2 == 'A')
      {
       if(c3 == 'T') amac = 11;
       if(c3 == 'C') amac = 11;
       if(c3 == 'A') amac = 10;
       if(c3 == 'G') amac = 10;
      }
   if(c2 == 'G') amac = 12;
   }

 if(c1 == 'A')
   {
    if(c2 == 'T')
      {
       if(c3 == 'T') amac = 15;
       if(c3 == 'C') amac = 15;
       if(c3 == 'A') amac = 15;
       if(c3 == 'G') amac = 14;
      }
    if(c2 == 'C') amac = 3;
    if(c2 == 'A')
      {
       if(c3 == 'T') amac = 7;
       if(c3 == 'C') amac = 7;
       if(c3 == 'A') amac = 13;
       if(c3 == 'G') amac = 13;
      }
    if(c2 == 'G')
      {
       if(c3 == 'T') amac = 2;
       if(c3 == 'C') amac = 2;
       if(c3 == 'A') amac = 12;
       if(c3 == 'G') amac = 12;
      }
   }

 if(c1 == 'G')
   {
    if(c2 == 'T') amac = 17;
    if(c2 == 'C') amac = 5;
    if(c2 == 'A')
      {
       if(c3 == 'T') amac = 8;
       if(c3 == 'C') amac = 8;
       if(c3 == 'A') amac = 9;
       if(c3 == 'G') amac = 9;
      }
    if(c2 == 'G') amac = 6;
   }


if( 
    ( c1 != 'A'  &&  c1 != 'T'  &&  c1 != 'G'  &&  c1 != 'C' )  ||
    ( c2 != 'A'  &&  c2 != 'T'  &&  c2 != 'G'  &&  c2 != 'C' )  ||
    ( c3 != 'A'  &&  c3 != 'T'  &&  c3 != 'G'  &&  c3 != 'C' )  
  ) 
 return( 0 );
else
 return( amac );

 }  /*  translate */

int int_test(float f)
{
 int i = f;

 if(i == f) 
     return (1);
 else  
     return (0);
}



/*==========================================================
 *         OLD SORT FUNCTION (BUBBLE-SORT)
 *=========================================================*/


void change(struct multi_frag *a, struct multi_frag *b)
{
  struct multi_frag c, *an, *bn;
  
  c = *a;
  an = a->next;
  bn = b->next;
  
  *a = *b;
  *b = c;
  
    a->next = an;
    b->next = bn;
}


void pair_change(struct seq_pair *a, struct seq_pair *b)
{
  struct seq_pair c;
  
  c = *a;
  *a = *b;
  *b = c;
}


void ow_bubble_sort( int number , struct multi_frag *dp )
{ 
  /* sorting diagonals in multiple alignment according to their
     overlap weights */
  
  struct multi_frag *hp;
  int hv1, hv2;
  
  FILE *fp_st;

  for( hv1 = 1 ; hv1 < number ; hv1++ )
    { 
      hp = dp;
      
      if( print_status )
        if( ( hv1 % 100 ) == 0 )
          {         
            fp_st = fopen( pst_name ,"w");
	    
            fprintf(fp_st,"\n\n\n    Status of the program run:\n");  
            fprintf(fp_st,"    ==========================\n\n");  
            fprintf(fp_st,"      %s \n\n", input_line);
            fprintf(fp_st,"      iteration step %d in multiple alignment\n", istep );
            fprintf(fp_st,"      overlap weight sorting of diagonals\n");
            fprintf(fp_st,"      current diagonal = %d\n\n", hv1 );
            fprintf(fp_st,"      total number of"); 
            fprintf(fp_st," diagonals: %d\n\n\n\n", number);
            fclose(fp_st);
          }
      
      
      for( hv2 = hv1 ; hv2 < number ; hv2++ )
	{
	  if( hp->ow < (hp->next)->ow )
	    change( hp , hp->next );
            hp = hp->next;
	} 
    }
} /*  ow_bubble_sort */





void bubble_sort( int number , struct multi_frag *dp )
{ 
  /* sorting diagonals in multiple alignment according to their
       individual weights */
  
  struct multi_frag *hp;
  int hv1, hv2;
  FILE *fp_st;

  for( hv1 = 1 ; hv1 < number ; hv1++ )
    { 
      hp = dp;
      
      if( print_status )
        if( ( hv1 % 100 ) == 0 )
          {   
            fp_st = fopen( pst_name ,"w");
	    
            fprintf(fp_st,"\n\n\n    Status of the program run:\n");  
            fprintf(fp_st,"    ==========================\n\n");  
            fprintf(fp_st,"      %s \n\n", input_line);
            fprintf(fp_st,"      iteration step %d\n", istep );
            fprintf(fp_st,"      ind. weight sorting of diagonals\n");
            fprintf(fp_st,"      current diagonal = %d\n\n", hv1 );
            fprintf(fp_st,"      total number of"); 
            fprintf(fp_st," diagonals: %d\n\n\n\n", number);
            fclose(fp_st);
          }
      
      
      for( hv2 = hv1 ; hv2 < number ; hv2++ )
	{ 
	  if( hp->weight < (hp->next)->weight )
	    change( hp , hp->next );
	  hp = hp->next;
	}  
    }
  
} /*  bubble_sort */



/*==========================================================
 *         NEW SORT FUNCTION (QUICK-SORT)
 *=========================================================*/

/***********************************************************
 *                   change()
 ***********************************************************/
void change_struct_el(struct multi_frag **a, int l, int r)
{
  struct multi_frag *dummy;
  dummy = a[l];
  a[l]  = a[r];
  a[r]  = dummy;  
}
/***********************************************************
 *                    change_first()
 ***********************************************************/
void change_first(struct multi_frag *a, struct multi_frag *b)
  {
    struct multi_frag c, *an, *bn;

    if(a==b)
      {  /* Change the first list-element with the second one (old first-el.). */
	c  = *a;
	an = a->next;
	bn = a->next->next;

	*a = *(a->next);
	a->next = bn;
	
	*an = c;
	an->next = a;
      }
    else /* Change the new first list-el. with the old first-el. */
      {
	c  = *a;              /* Make a copy of the new first-listelement a. */
	an = a->next;         /* Make a copy of the pointer at the second-el. */ 
	bn = b->next->next;   /* Make a copy of the pointer old first-el. shows at. */ 
	
	*a      = *(b->next); /* Whrite the value of the old first-el. on the place of the new first-el.*/ 
	a->next = bn;         /* Bend his "next" pointer at the next el. of the old first-el. */ 
	
	*(b->next)    = c;    /* Whrite the value of the new fist-el. on the place of the old first-el. */ 
	b->next->next = an;   /* Bend his "next" pointer at the next el. of the new first-el. */ 

	b->next       = a;    
      }
  }


/***********************************************************
 *                    quicksort_ow()
 ***********************************************************/
void quicksort_ow(struct multi_frag **array,int left, int right)
{
  int l = left, r = right;
  struct multi_frag *element;
  element = array[(left+right)/2];

  do
    {
      while(array[l]->ow > element->ow)  l++;
      while(element->ow  > array[r]->ow) r--;
	  
      if(l < r)  change_struct_el(array,l,r);
      if(l <= r) {l++; r--;}
    }while(l<=r);
  
  if(left < r)  quicksort_ow(array, left, r);
  if(l < right) quicksort_ow(array, l, right);

}/*quicksort_ow*/



/***********************************************************
 *                    quicksort_weight()
 ***********************************************************/
void quicksort_weight(struct multi_frag **array,int left, int right)
{
  int l = left, r = right;
  struct multi_frag *element;
  element = array[(left+right)/2];
  
  do
    {
      while(array[l]->weight > element->weight)  l++;
      while(element->weight  > array[r]->weight) r--;
      
      if(l < r)  change_struct_el(array,l,r);
      if(l <= r) {l++; r--;}
    }while(l<=r);
  
  if(left < r)  quicksort_weight(array, left, r);
  if(l < right) quicksort_weight(array, l, right);

}/*quicksort_weight*/


/***********************************************************
 *                    assemble_list()
 ***********************************************************/
void assemble_list(struct multi_frag **array, struct multi_frag *dp,int number)
{
  int i,index=0;
  for (i = 0; i< number-1; i++)
    {
      if(dp==array[i])
	index = i;
      array[i]->next = array[i+1];
    }
  
  array[number-1]->next = 0;
  if(dp==array[number-1])
    index = number-1;
  
  if(index!=0)
    change_first(array[0],array[index-1]);
} /* assemble_list */


/***********************************************************
*                      frag_sort() 
************************************************************/
void frag_sort(int number , struct multi_frag *dp , int olw ) 
{
  int i=1;

  struct multi_frag **array;
  if((array = (struct multi_frag**)calloc(number+1,sizeof(struct multi_frag*)))==0)
    {
      printf(" problems with memory allocation for `all_clades'\n \n");
      exit(1);
    }
  
  array[0] = dp;
  while(array[i-1]->next)
      {array[i] = array[i-1]->next; i++;}

  if( olw )
    quicksort_ow(array,0,number);
  else 
    quicksort_weight(array,0,number);
  
  
  assemble_list(array, dp, number+1);

}/* frag_sort */



void ow_add( struct multi_frag *sm1 , struct multi_frag *sm2 )
{
 /* increasing the overlap weights of two diagonals, if they
    have any overlap */

 int trans, i, j, k, s1, s2, b1, b2, conslen, dif, match;
 float add_wgt;

 trans = sm1->trans;

 for(i=0;i<2;i++)
 for(j=0;j<2;j++)
 if( sm1->s[i] == sm2->s[j] )
 if( sm1->s[j] != sm2->s[i] )
 if( sm1->b[i] < sm2->b[j] + sm2->ext && 
     sm2->b[j] < sm1->b[i] + sm1->ext )
      {
        conslen = mini2( sm1->b[i] + sm1->ext, sm2->b[j] + sm2->ext)
                  - maxi2( sm1->b[i] , sm2->b[j] );        
        if( 
            ( trans == 0 ) ||
            ( ( conslen % 3 ) == 0 )
          )    { 
 
          s1 = sm1->s[(i+1)%2];
          s2 = sm2->s[(j+1)%2];
     
          b1 = sm1->b[(i+1)%2];
          dif = sm2->b[j] - sm1->b[i];
          if (dif > 0)
            b1 = b1 + dif;

          b2 = sm2->b[(j+1)%2];
          dif = sm1->b[i] - sm2->b[j];
          if (dif > 0)
            b2 = b2 + dif;

          match = 0;

          for( k = 0 ; k < conslen ; k++ ) {
            if( 
                ( wgt_type == 0 ) ||
                ( trans && ( ( k % 3 ) == 0 ) ) 
              ) 
              match = match
                + sim_score[ amino[ s1 ][ b1 + k ] ][ amino[s2][ b2 + k ] ];
            else
              match = match + ( seq[ s1 ][ b1 + k ] == seq[ s2 ][ b2 + k ] );
          }


          if( wgt_type == 0 )
            add_wgt = wgt_prot[ conslen ][ match ];
          else
            if( trans )
              add_wgt = wgt_trans[ conslen / 3 ][ match ] ;
            else
              add_wgt = wgt_dna[ conslen ][ match ] ;   

          sm1->ow = sm1->ow + add_wgt ; 
          sm2->ow = sm2->ow + add_wgt ; 

        }
      }     
}    /*  ow_add  */

void seq_shift()
 {
  int i, hv;
  
  for(i = 0 ; i < seqnum ; i++)
  for(hv = seqlen[i]+1 ; hv > 0 ; hv--)
     seq[i][hv] = seq[i][hv-1];
 }

     


void filter(int *number, struct multi_frag *diagonal)
{
   /* checks diagonals one by one, if they are consistent with the 
      diagonals already included into the alignment. If a new diagonal 
      is consistent, it is included into the alignment and the frontiers 
      in clos (when GABIOS is used) are changed accordingly */

 int i, j, k, l, sv, hv, ab[2], as[2], ae[2], aext, nv;
 float awgt ; 

 int test;        /* = 1 if current diagonal consistent; = 0 otherwise */
 int number_bf;   /* number of diagonals before filter */ 

 FILE *fp_st, *fp_cap ;

 float lb, ub;
 
 struct multi_frag *dia; 
 char cap_file_name[ NAME_LEN ] ;
 
 if( ( istep > 0 ) && ( iter_cond_prob == 0 ) )  
 for( i = 0 ; i < seqnum ; i++ )
 for( j = 0 ; j < seqnum ; j++ )
   cont_it_p[i][j] = 0 ;

 dia = diagonal;
 number_bf = *number;

 if( ( istep == 0 ) && anchors && ( seqnum > 2 ) ) {
   strcpy( cap_file_name , input_name );
   strcat( cap_file_name , ".cap" );
   fp_cap = fopen( cap_file_name ,"w");
 }
   

 for(nv = 0 ; nv < number_bf ; nv++ )
   {
     ab[0] = dia->b[0];         /* begin of n-th diagonal in 1. sequence */
     ab[1] = dia->b[1];         /* begin of n-th diagonal in 2. sequence */
     as[0] = dia->s[0];         /* 1. sequence of n-th diagonal */
     as[1] = dia->s[1];         /* 2. sequence of n-th diagonal */
     aext  = dia->ext;          /* length of n-th diagonal */
     awgt  = dia->weight;          /* length of n-th diagonal */
     ae[0] = ab[0] + aext - 1;  /* end of n-th diagonal in 1. sequence */
     ae[1] = ab[1]+aext-1;      /* end of n-th diagonal in 2. sequence */


     if( print_status )
     if( ( ( nv + 1 ) % 10 ) == 0 )
       {  
         fp_st = fopen( pst_name ,"w");

         fprintf(fp_st,"\n\n\n    Status of the program run:\n");  
         fprintf(fp_st,"    ==========================\n\n");  
         fprintf(fp_st,"      %s \n\n", input_line);
         fprintf(fp_st,"      iteration step %d \n", istep );
         fprintf(fp_st,"      checking diagonal %d for ", nv + 1);
         fprintf(fp_st,"consistency\n\n      total number of"); 
         fprintf(fp_st," diagonals = %d \n\n\n\n", number_bf);
         fclose(fp_st);
       }

     test = alignableSegments(clos, as[0], ab[0], as[1], ab[1], aext);
       
     if(test)      /* i.e current diagonal consistent with the diagonals
                      already included into the alignment */
       {

       addAlignedSegments(clos, as[0], ab[0], as[1], ab[1], aext);
  
         if( istep )  
         for(hv=0;hv<aext;hv++)
         for(i=0;i<2;i++) 
           {
             j = (i+1)%2;
             open_pos[ as[i] ][ as[j] ][ ab[i]+hv ] = 0;
           }
   
         dia->sel = 1;   
         glob_sim[ as[0] ][ as[1] ] =  
           glob_sim[ as[0] ][ as[1] ] + dia->weight;

         if( istep )
           tot_weight = tot_weight + dia->weight;
   

       }   /* if test, i.e. current diagonal consistent */   
     else  /* no consistency */
       {
         (*number)--;
         dia->sel = 0;   
         cont_it_p[ as[0] ][ as[1] ] = 1 ; 
       }

     if( ( istep == 0 ) && anchors && ( seqnum > 2 ) ) {
       fprintf( fp_cap, " anchor %d %d %d %d %d %f " , as[0] + 1, as[1] + 1 , ab[0], ab[1], aext , awgt);
       if( dia->sel == 0 )
          fprintf( fp_cap , " inconsistent ");
       fprintf( fp_cap , "\n");
     }


     dia = dia->next;

   }     /*  for(hv = 0 ; hv < number_bf ; hv++ )  */
                    
  if( ( istep == 0 ) && anchors && ( seqnum > 2 ) ) 
    fclose(  fp_cap  ) ;

}          /*   filter(  )  */




void sel_test()
  {
    int hv;
    struct multi_frag *hp;   

    hp = this_it_dia;   

    for( hv = 0 ; hv < num_dia_af[ istep ] ; hv++ )
      {
        if( hp->sel == 0 )   
          {  
            printf("\n \n \n   WARNING:   \n \n \n");
            printf(" sel[%d] = %d \n", hv, hp->sel);
            exit(2);
          }
        hp = hp->next;
      }
  }

    








void throw_out( float *weight_sum )
  {
    int nc;
    short consist_found = 0; 
  
    struct multi_frag *cp;  /* current diagonal */
    struct multi_frag *hp;  /* predecedor of cp */ 
    
    hp = ( struct multi_frag *) calloc( 1 , sizeof( struct multi_frag ) ); 
    cp = this_it_dia;
    hp = NULL;
    *weight_sum = 0;


    for( nc = 0 ; nc <  num_dia_bf[ istep ] ; nc++ )
      {
        if( cp->sel )
          {
            *weight_sum = *weight_sum + cp->weight;
            consist_found = 1;

            hp = cp;
            cp = cp->next;  
          } 
        else
          {
            cp = cp->next;
            if( consist_found ) 
              {
                free(hp->next);
                hp->next = cp;              
              } 
            else
              {
                free( this_it_dia);
                this_it_dia = cp;
              }
          }
      }
  }  /* throw_out */










void new_shift(int s, int p, int dif)
  /* shifts the elements of sequence s starting with  position p 
     for dif elements to the right */
  {
    int hv;
    int shift_dif; /* length of a gap (if existing) between position hv
                      and position hv+1. In case of gaps, the function
                      `new_shift' diminishs the lengths of the gaps instead
                      of shifting further sequence elements to the right  */

    for(hv=p ; ( hv<seqlen[s]+1 ) && (dif>0) ; hv++)
      {
        shift_dif = shift[s][hv+1] - shift[s][hv] - 1;
        shift[s][hv] = shift[s][hv] + dif;
        dif = dif -  shift_dif;
      }
  }

wgt_type_count( int num , int e_len, int *plus_cnt, int *minus_cnt,
           int *nuc_cnt , int *frg_inv, struct multi_frag *dia )      {

int i, dc, pc, s1, pos;

for( dc = 0 ; dc < num ; dc++ )
      {

        for( pc = 0 ; pc < dia->ext ; pc++ )
          {
            i  = dia->b[0] + pc;
            s1 = dia->s[0];
            pos = shift[s1][i];
            if ( dia->trans )
              if ( dia->cs )
                minus_cnt[ pos ] = minus_cnt[ pos ] + 1 ;
              else
                plus_cnt[ pos ] = plus_cnt[ pos ] + 1 ;
            else {
              nuc_cnt[ pos ] = nuc_cnt[ pos ] + 1 ;
            }
            frg_inv[ pos ] = frg_inv[ pos ] + 1 ;
          }
        dia = dia->next;
      }
}



plot_calc( int num , int e_len, float *w_count, float *pl, 
           struct multi_frag *dia , FILE *fp_csc )
  {
    int i, dc, pc, s1, pos;
    float max_weight = 0;     /* maximum value of `weight_count' */
    float shrink, shrink_csc, hsc ;
   
    for( dc = 0 ; dc < num ; dc++ )
      {
  
        for( pc = 0 ; pc < dia->ext ; pc++ )
          {
            i  = dia->b[0] + pc;
            s1 = dia->s[0];
            pos = shift[s1][i]; 
            w_count[ pos ] = w_count[ pos ] + dia->weight;  
          }   
        dia = dia->next;
      }
  
  
    for( i = 0 ; i <= e_len ; i++ )
      if( max_weight < w_count[i] )
        max_weight = w_count[i];
  
  
    if( max_weight )
      {
        shrink = plot_num / max_weight;
        shrink_csc = MAX_CSC / max_weight;
  
        for( i = 0 ; i <= e_len ; i++ )
          pl[i] = w_count[i] * shrink;

        if( col_score ) {
          printf(" e_len = %d \n\n", e_len) ; 
          for( i = 0 ; i <= e_len ; i++ ) {
            hsc = w_count[i] * shrink_csc ;
            fprintf( fp_csc , "%5.1f\t0\n", hsc ) ;
          }
        }
      }
    else { 
      for( i = 0 ; i <= e_len ; i++ )
        pl[i] = 0 ; 

      printf(" e_len = %d \n\n", e_len) ; 
      printf(" no max weight\n\n"); 
    }

  }  



void av_tree_print()
  {
    int i, j, k, connect, max_pair[2], cv, m1, m2;
    struct subtree *all_clades;
    double **clade_similarity, new_similarity; 
    double max_sim; 
    char *string,  l_name[2][20];   
    char tree_name[NAME_LEN];
    float max_seq_sim, branch_len[2], depth; 

    FILE *t_file;



    if( (all_clades = (struct subtree *) 
        calloc( seqnum , sizeof( struct subtree ) )) == NULL)
      {
        printf(" problems with memory allocation for `all_clades'\n \n");
        exit(1);
      }


    if( (clade_similarity = (double **) 
         calloc( seqnum , sizeof( double* ) )) == NULL)
      exit(1);

    for(i = 0 ; i < seqnum ; i++ )
    if( (clade_similarity[i] = (double *) 
         calloc( seqnum , sizeof( double ) )) == NULL)
      exit(1);

    if( (string = (char *) 
      calloc( seqnum * 100 , sizeof(char) )) == NULL)
      {
        printf(" problems with memory allocation for `string'\n \n");
        exit(1);
      }




    for(i = 0 ; i < seqnum ; i++ ) 
       { 
         if( (all_clades[i].member = (int *) 
             calloc( seqnum , sizeof( int ) )) == NULL)
           {    
             printf(" problems with memory allocation for `all_clades'\n \n");
             exit(1);
           }


         if( (all_clades[i].name = (char *) 
             calloc( seqnum * 100 , sizeof( char ) )) == NULL)
           {    
             printf(" problems with memory allocation for `all_clades'\n \n");
             exit(1);
           }
         
         strcpy( all_clades[i].name , seq_name[i] );
         all_clades[i].member_num = 1;
         all_clades[i].member[0] = i;
         all_clades[i].valid = 1;
         all_clades[i].depth = 0;
       } 



    for(i = 0 ; i < seqnum ; i++ ) 
    for(j = i + 1 ; j < seqnum ; j++ ) 
      {
        clade_similarity[i][j] =  glob_sim[i][j];
        clade_similarity[j][i] =  glob_sim[i][j];
      } 


    for(connect = 1 ; connect < seqnum ; connect++)
      {
        max_sim = - 1;

   

        for(i = 0 ; i < seqnum ; i++ ) 
        for(j = 0 ; j < seqnum ; j++ ) 
        if( i != j )
        if( all_clades[i].valid && all_clades[j].valid )    
        if( clade_similarity[i][j] > max_sim )
          {
            max_sim =  clade_similarity[i][j];
            max_pair[0] = i;  
            max_pair[1] = j;  
          }  

  
        depth = 1 / ( max_sim + 1 ) ; 

          {
            m1 = max_pair[0];  
            m2 = max_pair[1];  

            for( i = 0 ; i < seqnum ; i++ )
            if( all_clades[i].valid )        
            if( i != m1 )
            if( i != m2 )      
              {
                if( ! strcmp(clust_sim , "av") )
                  new_similarity = 
                    ( 
                      clade_similarity[i][m1] * all_clades[m1].member_num + 
                      clade_similarity[i][m2] * all_clades[m2].member_num  
                    ) /
                     ( all_clades[m1].member_num + all_clades[m2].member_num );
 
                if( ! strcmp(clust_sim , "max") )
                  new_similarity = 
                    maxf2(  clade_similarity[i][m1] ,  clade_similarity[i][m2] );
		    
                if( ! strcmp(clust_sim , "min") )
                  new_similarity =
                    minf2(  clade_similarity[i][m1] ,  clade_similarity[i][m2] );


                clade_similarity[i][m1] = new_similarity;
                clade_similarity[m1][i] = new_similarity;
              }


            all_clades[m2].valid = 0;

            for(k = 0 ; k <  all_clades[m2].member_num  ; k++)
              all_clades[m1].member[ all_clades[m1].member_num + k ] = 
              all_clades[m2].member[ k ] ;

            all_clades[m1].member_num = 
               all_clades[m1].member_num + all_clades[m2].member_num;


            for(k = 0 ; k < 2 ; k++)
              {
                branch_len[k] = depth - all_clades[ max_pair[k] ].depth;
                sprintf( l_name[k],":%f", branch_len[k]);
              } 


            all_clades[m1].depth = depth;

       
            strcpy(string,"(");
            strcat(string, all_clades[m1].name); 
            strcat(string,l_name[0]); 
/*            strcat(string,",\n");   */ 
            strcat(string, all_clades[m2].name); 
            strcat(string,l_name[1]); 
            strcat(string,")");

            strcpy( all_clades[m1].name , string ); 
          }
      }


    strcat(string, ";"); 

    i = strlen( string ) + 2;

    if( (upg_str = (char *) calloc( i , sizeof(char) )) == NULL)
      {
        printf(" problems with memory allocation for `upg_str'\n \n");
        exit(1);
      }

    for(i = 0 ; i <= strlen( string ) ; i++ )
      upg_str[i] = string[i] ;


  }   /*  av_tree_print  */



void print_log( struct multi_frag *d , FILE *fp_l , FILE *fp_fs )
  {
    int i, j, mn, pv, percent, this_frag_trans , frg_count = 0 ;
    struct multi_frag *diagonal;
    char hc ;

    if(long_output)
      {
        fprintf(fp_l," \n \n  Iteration %d:\n", istep );

        if( istep < 10 ) 
          fprintf(fp_l,"  ------------");
        else
          fprintf(fp_l,"  -------------");
      } 


    for(i= 0 ; i<seqnum ; i++)
    for(j= i+1 ; j<seqnum; j++)
      {
        if(long_output) {
          if( seqnum > 2 ) {
            fprintf(fp_l, "\n \n \n \n  Pairwise alignment ");
            fprintf(fp_l, "%d/%d", i + 1, j + 1); 
            fprintf(fp_l, " (%s / %s) \n" ,seq_name[i] ,seq_name[j] );           
            fprintf(fp_l, "  =========================");
            fprintf(fp_l, "===================== ");
          }
          fprintf(fp_l, " \n \n \n");
        }

        pairalignsum = 0;
        pairalignlen = 0;

        diagonal = d;
        while(diagonal != NULL)
          {
            frg_count++ ;
            if( diagonal->s[0] == i && diagonal->s[1] == j)
              {
                if(diagonal->sel)
                  {
                    if(long_output)
                      {
                        fprintf(fp_l,"   *");
                        fprintf(fp_l," (%3d,", diagonal->b[0]);
                      }

                    pairalignsum = pairalignsum + diagonal->weight;
                    pairalignlen = pairalignlen + diagonal->ext;
                  }
                else
                if(long_output)
                  fprintf(fp_l,"     (%3d,", diagonal->b[0]);

                if(long_output)
                  {
                    fprintf(fp_l,"%3d)  ", diagonal->b[1]);
                    fprintf(fp_l," wgt:%7.3f ", diagonal->weight);
                    if(seqnum > 2) 
                    if(overlap_weights)
                      fprintf(fp_l," olw:%7.3f ", diagonal->ow);
                    fprintf(fp_l,"len: %2d", diagonal->ext);
                    if( ( wgt_type == 3 ) || crick_strand ) {
                      if( diagonal->trans )
                        fprintf(fp_l,"  P-frg" );
                      else
                        fprintf(fp_l,"  N-frg" );
                    }

                    if( diagonal->trans )
                    if( crick_strand ) {
                      if( diagonal->cs )
                        fprintf(fp_l,", CRICK strand " );
                      else
                        fprintf(fp_l,", WATSON strand " );
                    }

                  }

                if( frg_mult_file_v ) {
                  fprintf(fp_fs,"FRG %d ", frg_count ); 
                  fprintf(fp_fs,"name: %s %s ", seq_name[i] , seq_name[j] ) ;  
 
                  fprintf(fp_fs,"seq: %d %d ", i + 1 , j + 1 ) ;  
                  fprintf(fp_fs,"beg: %d %d ", diagonal->b[0], diagonal->b[1]); 
                  fprintf(fp_fs,"len: %d ", diagonal->ext);

                  fprintf(fp_fs,"wgt:%7.3f ", diagonal->weight);
                  if(diagonal->sel) 
                    fprintf(fp_fs," CONS  "); 
                  else   
                  fprintf(fp_fs," NON-CONS ");
                  fprintf(fp_fs,"\n") ; 
                  fprintf(fp_fs,"SEG1   ");
                  for(pv = 0 ; pv < diagonal->ext ; pv ++)
                    fprintf(fp_fs,"%c", seq[i][ diagonal->b[0] + pv ] );
                  fprintf(fp_fs,"\n"); 

                  fprintf(fp_fs,"SEG2   ");
                  for(pv = 0 ; pv < diagonal->ext ; pv ++)
                    fprintf(fp_fs,"%c", seq[j][ diagonal->b[1] + pv ] );
                  fprintf(fp_fs,"\n"); 
                  fprintf(fp_fs,"\n"); 
                }
                if( frg_mult_file & ! frg_mult_file_v ) { 
                  if( diagonal->sel ) {
                    fprintf(fp_fs," %d %d ", i + 1 , j + 1 ) ;  
                    fprintf(fp_fs," %d %d ", diagonal->b[0], diagonal->b[1]); 
                    fprintf(fp_fs," %d \n", diagonal->ext);
                  } 
                }                

                if(long_output)
                  {
                    fprintf(fp_l,"\n");
 
                    if( 
                        wgt_type == 2 || 
                        ( ( wgt_type == 3 ) && diagonal->trans ) 
                      )
                      this_frag_trans = 1;
                    else
                      this_frag_trans = 0;
                           
                    if( this_frag_trans ) 
                      {
                        fprintf(fp_l,"\n           ");
                        for(pv = 0 ; pv < diagonal->ext ; pv ++)
                          { 
                            hc = amino_acid[ amino[i][ diagonal->b[0] + pv - 1 ] ] ;
                            if( crick_strand )
                            if( diagonal->cs )
                              hc = amino_acid[ amino_c[i][ diagonal->b[0] + pv - 1 ] ] ;
                          
                            if( ( pv % 3 ) == 0 )
                              fprintf(fp_l,"/");
                            if( ( pv % 3 ) == 1 )
                              fprintf(fp_l,"%c", hc ) ; 
                            if( ( pv % 3 ) == 2 )
                              fprintf(fp_l,"\\");

                          } 
                      }

                    fprintf(fp_l,"\n           "); 
                    for(pv = 0 ; pv < diagonal->ext ; pv ++)
                      fprintf(fp_l,"%c", seq[i][ diagonal->b[0] + pv ] );
                    fprintf(fp_l,"\n"); 


                    fprintf(fp_l,"           ");
                    for(pv = 0 ; pv < diagonal->ext ; pv ++)
                      fprintf(fp_l,"%c", seq[j][ diagonal->b[1] + pv ] );


                    if( this_frag_trans )
                      {
                        fprintf(fp_l,"\n           ");
                        for(pv = 0 ; pv < diagonal->ext ; pv ++)
                          {
                            hc = amino_acid[ amino[j][ diagonal->b[1] + pv - 1 ]
 ] ;
                            if( crick_strand )
                            if( diagonal->cs )
                              hc = amino_acid[ amino_c[j][ diagonal->b[1] + pv - 1 ] ] ;

                            if( ( pv % 3 ) == 0 )
                              fprintf(fp_l,"\\");
                            if( ( pv % 3 ) == 1 )
                              fprintf(fp_l,"%c", hc ) ;
                            if( ( pv % 3 ) == 2 )
                              fprintf(fp_l,"/");

                          }
                      }

                    fprintf(fp_l,"\n \n");
                  }   
              }   /*  if( diagonal->s[0] == i && diagonal->s[1] == j)  */

            diagonal = diagonal->next;

          }        /*  while(diagonal != NULL) */

        percent = pairalignlen*100/mini2(seqlen[i],seqlen[j]);

        if(long_output)
          {
            fprintf(fp_l,"\n      Sum of diagonal scores: %f\n", pairalignsum);
            fprintf(fp_l,"      Aligned residues: %d\n", pairalignlen);
            fprintf(fp_l,"      (%d percent of the shorter", percent);
            fprintf(fp_l," sequence aligned)\n");
          }
      }               /* for(i = 0     ; i < seqnum ; i++)
                         for(j = i + 1 ; j < seqnum ; j++)  */

  }    /*  print_log  */


   
