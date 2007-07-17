

                 /*******************\
                 *                   *
                 *     DIALIGN 2     *
                 *                   *
                 *     output.c      *
                 *                   *
                 \*******************/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "define.h"
#include "dialign.h"
#include "alig_graph_closure.h"

extern int cd_gobics , wgt_type_plot , col_score ; 
extern int ref_seq, anchors, speed_optimized, online ;
extern short crick_strand ;
extern float sf_mat_thr ; 
extern int **amino; 
extern char amino_acid[ 22 ] ; 
extern int quali_num , wgt_plot , mask, lgs_option ;
extern char input_line[ NAME_LEN ];
extern char clust_sim[ NAME_LEN ] ;
extern int msf_file, cw_file; 
extern int lmax;
extern float  threshold;
extern float av_len ;
extern int pr_av_max_nd , wgt_type ;
extern int num_dia_p, overlap_weights ;
extern int fasta_file;
extern char *upg_str;
extern int plot_num;
extern char *seq[MAX_SEQNUM];
extern int  *seqlen;
extern int  maxlen;
extern char *seq_name[MAX_SEQNUM];
extern char *full_name[MAX_SEQNUM];
extern int **shift; 

extern struct multi_frag *pair_dia;
extern struct multi_frag *this_it_dia;
extern struct multi_frag *all_it_dia;

extern CLOSURE *clos;

extern int max_sim_score ;
extern int max_dia;
extern int seqnum ;
extern int num_all_it_dia;
extern int frg_count ;

 extern int mini2(int a, int b) ;
 extern void new_shift(int s, int p, int dif) ;
 extern void mini(int *a, int b);
 extern void maxi(int *a, int b);
 extern int int_test(float f);
 extern plot_calc( int num , int e_len, float *w_count, float *pl,
           struct multi_frag *dia , FILE *fp_csc ) ;
 extern wgt_type_count( int num , int e_len, int *plus_cnt, int *minus_cnt,
           int *nuc_cnt , int *frg_inv, struct multi_frag *dia ) ;  



void subst_mat( char *file_name, int fragno , struct multi_frag *frg ) {

  int s0 , s1 , i , j , frg_count ; 
  short a0 , a1 ; 
  int ****sbsmt ;
  struct multi_frag *frag ; 
  char mat_file_name[ NAME_LEN ] ;
  FILE *fp_mat; 


  if( ( sbsmt = (int **** ) calloc( seqnum , sizeof(int ***))) == NULL) { 
    printf("Problems with memory allocation for sbsmt\n");
    exit(1);
    }

  for( i = 0 ; i < seqnum ; i++ ) 
    if( ( sbsmt[i] = (int *** ) calloc( seqnum , sizeof(int **))) == NULL) { 
      printf("Problems with memory allocation for sbsmt\n");
      exit(1);
    }

  for( i = 0 ; i < seqnum ; i++ )
  for( j = 0 ; j < seqnum ; j++ ) 
    if( ( sbsmt[i][j] = (int ** ) calloc( 21 , sizeof( int* ) ) ) == NULL) { 
      printf("Problems with memory allocation for sbsmt\n");
      exit(1);
    }

  for( i = 0 ; i < seqnum ; i++ )
  for( j = 0 ; j < seqnum ; j++ ) 
  for( a0 = 0 ; a0 < 21 ; a0++ ) 
    if( ( sbsmt[i][j][a0] = (int * ) calloc( 21 , sizeof( int ) ) ) == NULL) { 
      printf("Problems with memory allocation for sbsmt\n");
      exit(1);
    }

  for( i = 0 ; i <seqnum ; i++ )
  for( j = 0 ; j <seqnum ; j++ )
  for( a0 = 0 ; a0 <= 20 ; a0++ )
  for( a1 = 0 ; a1 <= 20 ; a1++ ) 
     sbsmt[ i ][ j ][ a0 ][ a1 ] = 0 ;


  strcpy( mat_file_name , file_name );
  strcat( mat_file_name , ".mat" );

  fp_mat = fopen( mat_file_name, "w") ;



  frag = frg ;  

  for( frg_count = 0 ; frg_count < fragno ; frg_count++ ) {
    if( frag->weight > sf_mat_thr )
    for( i = 0 ; i < frag->ext ; i++ ) {
      a0 = amino[ frag->s[0] ][ frag->b[0] + i ] ; 
      a1 = amino[ frag->s[1] ][ frag->b[1] + i ] ; 
      s0 = frag->s[0] ; 
      s1 = frag->s[1] ;
      sbsmt[ s0 ][ s1 ][ a0 ][ a1 ]++ ;
      sbsmt[ s1 ][ s0 ][ a1 ][ a0 ]++ ;

    }
    frag = frag->next ;  
  }


  fprintf( fp_mat, "taxanumber: %d ;\n", seqnum) ;
  fprintf( fp_mat, "description: DIALIGN alignment ;\n" ) ;
  fprintf( fp_mat, "description: %s;\n", input_line ) ;


  for( i = 0 ; i < seqnum ; i++ ) 
     fprintf( fp_mat, "taxon: %.3d  name: %s  ;\n", i + 1 , full_name[i] ) ;


  for( s0 = 0 ; s0 < seqnum ; s0++ )
  for( s1 = s0 + 1  ; s1 < seqnum ; s1++ ) 
  for ( a0 = 1 ; a0 <= 20 ; a0++ ) 
  for( a1 = 1 ; a1 < 21 ; a1++ )  {    
    fprintf( fp_mat, "pair: %.3d %.3d ", s0 + 1, s1 + 1 );  
    fprintf( fp_mat, " acids: %c%c ", amino_acid[a0] , amino_acid[a1] );  
    fprintf( fp_mat, " number: %d ;\n", sbsmt[ s0 ][ s1 ][ a0 ][ a1 ] ); 
  }
} /* subst_mat */


void print_fragments( struct multi_frag *d , FILE *fp_ff2 ) {

  struct multi_frag *fragment ;

  fragment = d;
  while( fragment != NULL ) {
    if( fragment->it ){
      frg_count++ ;
      fprintf( fp_ff2, "%6d) ", frg_count );
      fprintf( fp_ff2, "seq: %3d %3d  ", fragment->s[0] + 1 , fragment->s[1] + 1 );
      fprintf( fp_ff2, "beg: %7d %7d ", fragment->b[0] , fragment->b[1] );
      fprintf( fp_ff2, "len: %3d ", fragment->ext  );
      fprintf( fp_ff2, "wgt: %6.2f ", fragment->weight  );
      fprintf( fp_ff2, "olw: %6.2f ", fragment->ow );

      fprintf( fp_ff2, "it: %d ", fragment->it  );
      if( fragment->sel )
        fprintf( fp_ff2, "cons   " );
      else 
        fprintf( fp_ff2, "incons " );

      if( ( wgt_type == 3 ) || crick_strand ) { 
        if( fragment->trans )
          fprintf( fp_ff2, " P-frg" );
        else
          fprintf( fp_ff2, " N-frg" );
        if( fragment->trans )
        if( crick_strand ) 
          if( fragment->cs )
            fprintf( fp_ff2, " -" );
          else
            fprintf( fp_ff2, " +" );
      } 


      fprintf( fp_ff2, "\n" );
    }
    fragment = fragment->next ;
  }
}


void weight_print( float **wgt )
{
 int i, j , l, s ; 
 FILE *fp;

 fp = fopen("weight_table","w");
 

 fprintf(fp,"  len1 = %d, len2 = %d\n\n",seqlen[0], seqlen[1] );
 fprintf(fp,"  \n   %s \n\n", input_line );
 for( l = 1 ; l <= max_dia ; l++ )
   for( s = 0 ; s <= l * max_sim_score ; s++ )
     fprintf(fp," %d %d %7.8f \n", l, s, wgt[l][s] );

 fclose(fp);

}  /* weight_print */





void ali_arrange( int fragno , struct multi_frag *d, FILE *fp, FILE *fp2, FILE *fp3 , FILE *fp4 , FILE *fp_col_score  )
{
 int block_no, char_no ;
 int shift_cond, endlen;
 int  p, pn, i, j, k, l, hv,  bc, lc,  max_p;
 int b1, b2, s1, s2, e, dif, sv, lv, add, msf_lines;

 char sim_char;
 float weak_wgt_type_thr = WEAK_WGT_TYPE_THR ; 
 float strong_wgt_type_thr = STRONG_WGT_TYPE_THR  ; 
 float frac_plus, frac_minus, frac_nuc, f_inv ;
 
 char **endseq;
 char **hseq;
 char *clear_seq;
 float *weight_count;
 int *plus_count;
 int *minus_count;
 int *nuc_count;
 int *frg_involved;
 float *plot;        /* plot[i] = sum of weights of fragments involved at
                        position i normalizet such that the maximum value */
 
 char gap_char = '-';
 char ambi_char = ' ';
 int *begin, *end, *b_len, *first_pos, pl_int ;
 int b_size;         /* size of fragments */
 struct multi_frag *fragments, *dia; 
 int **inv_shift;     
 int char_per_line;  /* number of residues per line in output file        */
 char aligned; 
 char_per_line = ( ( PAPER_WIDTH - 18 ) / 11) * 10;

 dia = d;

 if( (endseq = (char **) calloc( seqnum , sizeof(char *) )) == NULL)
   {
     printf(" problems with memory allocation for `endseq' !  \n \n");
     exit(1);
   }

 if( (hseq = (char **) calloc( seqnum , sizeof(char *) )) == NULL)
   {
     printf(" problems with memory allocation for `hseq' !  \n \n");
     exit(1);
   }

 if( (begin = (int *) calloc( seqnum , sizeof(int) )) == NULL)
   {
     printf(" problems with memory allocation for `begin' !  \n \n");
     exit(1);
   }

 if( (end = (int *) calloc( seqnum , sizeof(int) )) == NULL)
   {
     printf(" problems with memory allocation for `end' !  \n \n");
     exit(1);
   }

 if( (b_len = (int *) calloc( seqnum , sizeof(int) )) == NULL)
   {
     printf(" problems with memory allocation for `b_len' !  \n \n");
     exit(1);
   }

 if( ( first_pos = (int *) calloc( seqnum , sizeof(int) )) == NULL)
   {
     printf(" problems with memory allocation for `first_pos' !  \n \n");
     exit(1);
   }

 if( (shift = (int **) calloc( seqnum , sizeof(int *) )) == NULL ) 
    {
       printf("not enough memory available for `shift' !!!!\n");   
       fprintf(fp,"not enough memory available for `shift' !\n");   
       exit(1); 
    }

 for(hv=0 ; hv<seqnum ; hv++)
 if( (shift[hv] = (int *) calloc( (seqlen[hv]+2) , sizeof(int) )) == NULL ) 
    {
       printf("not enough memory available for `shift' !!!!\n");   
       fprintf(fp,"not enough memory available for `shift' !\n");   
       exit(1); 
    }


   if( fragno >= 0 )
     {

       for(hv=0;hv<seqnum;hv++)
         {
           begin[hv] = seqlen[hv];
           end[hv] = 1;
         }


       if( fragno > 0 )
       if( ( fragments = calloc( fragno , sizeof(struct multi_frag) )) == NULL )
         {
           printf("not enough memory available for fragments!\n");   
           fprintf(fp,"not enough memory available for fragments!\n");   
           exit(1);
         } 

       for( hv = 1 ; hv <= fragno ; hv++)
         {
           fragments[hv-1] = *dia;
           dia = dia->next;
         }

       for( hv = 0 ; hv < fragno ; hv++ )
       for( j = 0 ; j < 2 ; j++ )
         {
           mini( &begin[ fragments[hv].s[j] ] , fragments[hv].b[j] );
           maxi( &end[ fragments[hv].s[j] ] , fragments[hv].b[j] + fragments[hv].ext );
         }

       for(hv=0;hv<seqnum;hv++)
         {
           begin[hv] = 1;
           end[hv] = seqlen[hv]+1; 
         }  
      
       b_size = 0;

       for(i=0;i<seqnum;i++)
         {
           b_len[i] = end[i] - begin[i];
           maxi(&b_size,b_len[i]);
         }  

       for(i=0;i<seqnum;i++)
       for(hv=0;hv<b_len[i];hv++)
         shift[i][ begin[i]+hv ] = hv;

       shift_cond = 1;
    
       while(shift_cond)
             {
               shift_cond = 0;
 
               for( hv = 0 ; hv < fragno ; hv++ )
                 for(j=0;j<2;j++)
                   {
                     k = (j+1)%2;
                     s1 = fragments[hv].s[j]; 
                     s2 = fragments[hv].s[k]; 
                     b1 = fragments[hv].b[j]; 
                     b2 = fragments[hv].b[k]; 
                     e = fragments[hv].ext; 

                     for(l = e-1;l>=0;l--)
                       {
                         dif =  shift[s2][b2+l] - shift[s1][b1+l]; 
                         if (dif > 0 )
                           {
                             new_shift(s1,b1+l,dif);
                             shift_cond = 1;
                           }   
                       }
                   }
             }       /*  while (shift_cond)  */







       endlen = 0;

       for(hv=0;hv<seqnum;hv++)
         maxi(&endlen,shift[hv][ end[hv]-1 ] + 1);  

       for(hv=0;hv<seqnum;hv++)
       if( (endseq[hv] = calloc(endlen, sizeof(char) )) == NULL )
         {
           printf(" not enough memory available for printing results!\n");   
           fprintf(fp," not enough memory available");   
           fprintf(fp," for printing results!\n");   
           exit(1);
         }
 

       if( (inv_shift = (int **) calloc( seqnum , sizeof(int *) )) == NULL ) 
         {
           printf("not enough memory available for `inv_shift' !!!!\n");   
           fprintf(fp,"not enough memory available for `inv_shift' !\n");   
           exit(1); 
         }

       for(hv=0 ; hv<seqnum ; hv++)
       if( (inv_shift[hv] = (int *) calloc( (endlen+2) , sizeof(int) )) 
            == NULL ) 
         {
           printf("not enough memory available for `inv_shift' !!!!\n");   
           fprintf(fp,"not enough memory available for `inv_shift' !\n");   
           exit(1); 
         }

       if( (clear_seq = (char *) calloc( (endlen+1) , sizeof(char) )) == NULL)
         {
           printf(" problems with memory allocation for `clear_seq' !  \n \n");
           exit(1);
         }

       if( (weight_count = 
             (float *) calloc( ( endlen + 2 ) , sizeof(float) )) == NULL)
         {
           printf(" problems with memory allocation for `weight_count' !\n \n");
           exit(1);
         }

       if( (plot = (float *) calloc( ( endlen + 2 ) , sizeof(float) )) == NULL)
         {
           printf(" problems with memory allocation for `plot' ! \n \n");
           exit(1);
         }

       if( (plus_count = 
             (int *) calloc( ( endlen + 2 ) , sizeof( int ) )) == NULL)
         {
           printf(" problems with memory allocation for `plus_count' !\n \n");
           exit(1);
         }

       if( (minus_count = 
             (int *) calloc( ( endlen + 2 ) , sizeof( int ) )) == NULL)
         {
           printf(" problems with memory allocation for `minus_count' !\n \n");
           exit(1);
         }

       if( (nuc_count = 
             (int *) calloc( ( endlen + 2 ) , sizeof( int ) )) == NULL)
         {
           printf(" problems with memory allocation for `nuc_count' !\n \n");
           exit(1);
         }

       if( (frg_involved = 
             (int *) calloc( ( endlen + 2 ) , sizeof( int ) )) == NULL)
         {
           printf(" problems with memory allocation for `frg_involved ' !\n \n");
           exit(1);
         }

   

       for(hv=0 ; hv<seqnum ; hv++)
       for(p=1 ; p <= seqlen[hv] ; p++)
         inv_shift[hv][ shift[hv][p] ] = p; 

       for(hv=0;hv<seqnum;hv++)
       if( (hseq[hv] = calloc( (maxlen+1), sizeof(char) )) == NULL )
         {
           printf("not enough memory available for printing results! \n");   
           fprintf(fp,"not enough memory available");   
           fprintf(fp," for printing results! \n");   
           exit(1);
         }
/*
printf("endlen = %d \n\n", endlen); 
*/ 

       for(hv=0;hv<seqnum;hv++)
       for(i=0;i<endlen;i++)
         endseq[hv][i] = gap_char;

       for(hv=0;hv<seqnum;hv++)
       for(i=begin[hv];i<end[hv];i++)
         hseq[hv][i] = tolower(seq[hv][i]);

       for( hv = 0 ; hv < fragno ; hv++ )
       for(k=0;k<2;k++)
       for(i = fragments[hv].b[k] ; i < fragments[hv].b[k] + fragments[hv].ext ; i++)
         hseq[ fragments[hv].s[k] ][i] = seq[ fragments[hv].s[k] ][i];

       for(hv=0;hv<seqnum;hv++)
       for(i = begin[hv] ; i < end[hv] ; i++)
         endseq[hv][ shift[hv][i] ] = hseq[hv][i];

       for(i=0;i<endlen;i++)
         clear_seq[i] = ' ';






       for(p=0;p<endlen;p++)
         {
           s1 = 0;
           while( 
                  ( endseq[s1][p] == tolower( endseq[s1][p] ) )
                   && (s1 < (seqnum - 1) )    /* no capital letter */
                )
           s1++;

           if(s1 < (seqnum - 1) )
             {
               for(s2 = s1+1 ; s2 < seqnum ; s2++)   
                 {
                   if( endseq[s2][p] != tolower( endseq[s2][p] ) )
                      /* endseq[s2][p] capital letter */ 
                     {
                         aligned = alignedPositions(clos,s1,inv_shift[s1][p],s2,
                            succFrontier(clos,s1,inv_shift[s1][p],s2));

                         if (!aligned)
                          /* i.e.endseq[s1][p] not aligned with end seq[s2][p]*/ 
                               clear_seq[p] =ambi_char;
   	             }
                 }
             }
         }


       if( mask )
       for(sv = 0 ; sv < seqnum ; sv++)
       for(hv = 0 ; hv < endlen ; hv++ )
       if( endseq[sv][hv] != gap_char )
       if( endseq[sv][hv] == tolower( endseq[sv][hv] ) )
         endseq[sv][hv] = '*' ;


       if( col_score ){
         fprintf(fp_col_score , "# 1 %d \n" , endlen  );
         fprintf(fp_col_score,"# %s \n", upg_str);
       }

       plot_calc( num_all_it_dia , endlen , weight_count , plot , all_it_dia , fp_col_score  );

       wgt_type_count( num_all_it_dia , endlen , plus_count, minus_count, nuc_count , frg_involved, all_it_dia );
      

       lc = (endlen-1)/char_per_line;
       for(hv=0;hv<seqnum;hv++)
         first_pos[hv] = begin[hv] ;


       for( k = 0 ; k <= lc ; k++ )
         {
           for( hv = 0 ; hv < seqnum ; hv++ )
             { 
               fprintf(fp, "%s", seq_name[hv] );

               fprintf(fp,"%8d  ", first_pos[hv]);


               for(i=0;i<mini2(char_per_line,endlen-k*char_per_line);i++)
                 {
                   if(!(i%10))fprintf(fp, " ");
                   fprintf(fp, "%c",endseq[hv][k*char_per_line+i]);
                   if(endseq[hv][k*char_per_line+i] != gap_char)
                   first_pos[hv]++;
                 }
               fprintf(fp, " \n");
             }

           fprintf(fp,"         ");
           for( i = 0 ; i < mini2( char_per_line , endlen-k*char_per_line )
                      ; i++ )
             {
               if(!(i%10))fprintf(fp, " ");
               fprintf(fp, "%c",clear_seq[k*char_per_line+i]);
             }

           if( plot_num )
             fprintf(fp, " \n");


      
           if( quali_num == 0 )  
           for( pn = 0 ; pn < plot_num ; pn ++ ) 
             { 
               fprintf(fp,"                      ");
               for(i=0;i<mini2(char_per_line,endlen-k*char_per_line);i++)
                 {
                   if( !(i%10) )fprintf(fp, " ");
                     if( plot[ k*char_per_line + i ]  >  pn )
                       fprintf(fp, "*");
                     else
                       fprintf(fp, " ");
                 }
               fprintf(fp, " \n");

               if( plot_num == 1 )  
                 fprintf(fp, " \n");
             }


           if( quali_num ) {
             for( i = 0 ; i < SEQ_NAME_LEN ; i++ ) { 
               fprintf(fp," ");
             }
            
             fprintf(fp,"          ");
             for( i = 0 ; i < mini2( char_per_line , endlen-k*char_per_line ) ; i++ ) { 
               if( !(i%10) )fprintf(fp, " ");
               pl_int = 9 * plot[ k * char_per_line + i ] / plot_num ;
               fprintf(fp, "%d", pl_int );
             } 
             fprintf(fp, " \n");
           } 

           /***********************************************************************

           fprintf(fp, " \n");
           if( wgt_type > 1 ) {
             for( i = 0 ; i < SEQ_NAME_LEN ; i++ ) { 
               fprintf(fp," ");
             }
            
             fprintf(fp,"  plus    ");
             for( i = 0 ; i < mini2( char_per_line , endlen-k*char_per_line ) ; i++ ) { 
               if( !(i%10) )fprintf(fp, " ");
               fprintf(fp, "%d", plus_count[ k * char_per_line + i ] );
             } 
             fprintf(fp, " \n");
           } 

           if( wgt_type > 1 ) {
             for( i = 0 ; i < SEQ_NAME_LEN ; i++ ) { 
               fprintf(fp," ");
             }
            
             fprintf(fp,"  minus   ");
             for( i = 0 ; i < mini2( char_per_line , endlen-k*char_per_line ) ; i++ ) { 
               if( !(i%10) )fprintf(fp, " ");
               fprintf(fp, "%d", minus_count[ k * char_per_line + i ] );
             } 
             fprintf(fp, " \n");
           } 

           if( wgt_type > 1 ) {
             for( i = 0 ; i < SEQ_NAME_LEN ; i++ ) { 
               fprintf(fp," ");
             }
            
             fprintf(fp,"  nuc     ");
             for( i = 0 ; i < mini2( char_per_line , endlen-k*char_per_line ) ; i++ ) { 
               if( !(i%10) )fprintf(fp, " ");
               fprintf(fp, "%d", nuc_count[ k * char_per_line + i ] );
             } 
             fprintf(fp, " \n");
             fprintf(fp, " \n");
           } 
      
           ************************************************************************/

           if( wgt_type_plot ) 
           if( wgt_type == 3 ) {

             fprintf(fp,"sim. level");

             for( i = 0 ; i < SEQ_NAME_LEN ; i++ ) { 
               fprintf(fp," ");
             }
            
             for( i = 0 ; i < mini2( char_per_line , endlen-k*char_per_line ) ; i++ ) { 
               if( !(i%10) )fprintf(fp, " ");
               sim_char = '.' ; 

               if( frg_involved[ k * char_per_line + i ] ) {   

                  f_inv = frg_involved[ k * char_per_line + i] ; 
                  frac_plus =  plus_count[ k * char_per_line + i ] /  f_inv ;
                  frac_minus =  minus_count[ k * char_per_line + i ] /  f_inv ;
                  frac_nuc =  nuc_count[ k * char_per_line + i ] /  f_inv ;

               if ( frac_plus > weak_wgt_type_thr )
                 if( crick_strand ) 
                   sim_char = 'f' ;
                 else  
                   sim_char = 'p' ;
               if ( frac_plus > strong_wgt_type_thr )
                 if( crick_strand ) 
                   sim_char = 'F' ; 
                 else 
                   sim_char = 'P' ; 
               if ( frac_minus > weak_wgt_type_thr )
                 sim_char = 'r' ; 
               if ( frac_minus > strong_wgt_type_thr )
                 sim_char = 'R' ; 

               if ( frac_nuc > weak_wgt_type_thr )
                 sim_char = 'n' ; 
               if ( frac_nuc > strong_wgt_type_thr )
                 sim_char = 'N' ; 

               }  
               fprintf(fp, "%c", sim_char );
             } 
             fprintf(fp, " \n");
             fprintf(fp, " \n");
           } 

       
       

	   
           fprintf(fp, " \n");
   
         }    /*   for(k=0;k<=lc;k++)  */
    

       if( fasta_file )
         { 
           for(sv = 0 ; sv < seqnum ; sv++ )
             {
               fprintf(fp2,">%s", full_name[sv]);
               for(i = 0 ; i < endlen ; i++)
                 {
                   if( ! ( i % 50 ) )   
                     fprintf(fp2,"\n"); 
                   fprintf(fp2,"%c", endseq[sv][i]);  
                 }
                    
               fprintf(fp2,"\n ");         
               if( sv < ( seqnum - 1 ) )
                 fprintf(fp2,"\n");             
             }
         }     
   
       
       if( cw_file )
	 {   
           block_no = 0;
 
           fprintf(fp4,"DIALIGN 2.1 multiple sequence alignment \n\n");
           fprintf(fp4,"// \n\n\n");
          
           while( block_no * 60 < endlen )
             {
               char_no = mini2( 60 ,  ( endlen - block_no * 60 ) ) ;
               for( sv = 0 ; sv < seqnum ; sv++ )
                 {
                   fprintf(fp4,"%s ", seq_name[sv] );
                   for( i = 0 ; i < char_no ; i++)
                     fprintf(fp4,"%c", endseq[sv][ block_no * 60 + i ] );
                   fprintf(fp4,"\n");
                 }
               fprintf(fp4,"\n\n");
               block_no++; 
             } 


	 }
	 

       if( msf_file )
         { 
           msf_lines = endlen / 50;
           if(endlen % 50)
             msf_lines = msf_lines + 1;
          

           fprintf(fp3,"DIALIGN 2\n\n\n");
           fprintf(fp3,"   MSF: %d \n\n", endlen);

           for( sv = 0 ; sv < seqnum ; sv++ )
             fprintf(fp3," Name: %s    Len: %d \n", seq_name[sv], seqlen[sv] );
           fprintf(fp3,"\n// \n\n");

           for(lv = 0 ; lv < msf_lines ; lv++ )
             {
               add = lv * 50;
               max_p = mini2( endlen - add , 50 );
           
               for( sv = 0 ; sv < seqnum ; sv++ )
                 {
                   fprintf(fp3, "%s", seq_name[sv] );
                   for(i=0 ; i < 4 ; i++ )
                     fprintf(fp3, " "); 
               
                   for(i = 0 ; i < max_p ; i++)
                     {
                       if( !(i%10) )fprintf(fp3, " ");
                       if(  endseq[sv][add + i] == '-' )
                         fprintf(fp3,".");
                       else
                         fprintf(fp3,"%c", endseq[sv][add + i]);
                     } 
                   fprintf(fp3,"\n"); 
		 }
               fprintf(fp3,"\n\n");
             } 

         }


       if( ( seqnum > 2 ) && ( ref_seq == 0 ) )  
         { 
           fprintf(fp,"\n \n \n   Sequence tree:\n");
           fprintf(fp,"   ==============\n\n");

           if( ! strcmp( clust_sim , "av" ) )
             fprintf(fp,"Tree constructed using UPGMA");
             fprintf(fp,"based on DIALIGN fragment weight scores");

           if( ! strcmp( clust_sim , "max" ) )
             fprintf(fp,"Tree constructed using maximum linkage clustering");


           if( ! strcmp( clust_sim , "min" ) )
             fprintf(fp,"Tree constructed using minimum linkage clustering");


           fprintf(fp,"\n \n%s", upg_str);
         }

       fprintf(fp,"\n \n \n");
           
       for(hv=0;hv<seqnum;hv++)
         free(hseq[hv]);  

       for(hv=0;hv<seqnum;hv++)
         free(endseq[hv]);

       if( fragno > 0 )
         free( fragments );

       free(plot);

       free(weight_count);


     }    /* for(bc=0;bc<1;bc++) */


  for(hv=0;hv<seqnum;hv++)
    free(shift[hv]);  






}     /*  ali_arrange  */                     


void para_print( char *s_f, FILE *fpi )
  {
      int p_count = 1;
      int hv, i ; 
   

        {
 


          if( cd_gobics ) {
            fprintf(fpi," \n                          CHAOS / DIALIGN  \n");
            fprintf(fpi,"                          ***************\n \n");
 
            if( BETA )
              fprintf(fpi,"                           beta version\n\n"); 

            fprintf(fpi,"          Program code written by \n");
            fprintf(fpi,"          Burkhard Morgenstern, Said Abdeddaim and Michael Brudno \n\n");
            fprintf(fpi,"             e-mail contact: ");
            fprintf(fpi,"dialign (at) gobics (dot) de \n \n");
            fprintf(fpi,"          Published research assisted");
            fprintf(fpi," by CHAOS / DIALIGN should cite:  \n \n");
            fprintf(fpi,"             Michael Brudno et al.");
            fprintf(fpi," (2003)\n");
            fprintf(fpi,"             \"Fast and sensitive multiple alignment");
            fprintf(fpi," of large genomic sequences\" \n"); 
            fprintf(fpi,"             BMC Bioinformatics 4:66 \n");
            fprintf(fpi,"             http://www.biomedcentral.com/1471-2105/4/66 \n\n");
          }
          else {
            fprintf(fpi," \n                           DIALIGN 2.2.1 \n");
            fprintf(fpi,"                           *************\n \n");
 
            if( BETA )
              fprintf(fpi,"                           beta version\n\n"); 

            fprintf(fpi,"          Program code written by Burkhard");
            fprintf(fpi," Morgenstern and Said Abdeddaim \n");
            fprintf(fpi,"             e-mail contact: ");
            fprintf(fpi,"dialign (at) gobics (dot) de \n \n");
            fprintf(fpi,"          Published research assisted");
            fprintf(fpi," by DIALIGN 2 should cite:  \n \n");
            fprintf(fpi,"             Burkhard Morgenstern");
            fprintf(fpi," (1999).\n");
          
            fprintf(fpi,"             DIALIGN 2: improvement of the");
            fprintf(fpi," segment-to-segment\n             approach");
            fprintf(fpi," to multiple sequence alignment.\n");
            fprintf(fpi,"             Bioinformatics 15,");
            fprintf(fpi," 211 - 218. \n\n");
          }

          fprintf(fpi,"          For more information, please visit");
          fprintf(fpi," the DIALIGN home page at \n\n             ");
          fprintf(fpi,"http://bibiserv.techfak.uni-bielefeld.de/dialign/");
          fprintf(fpi," \n \n");

            fprintf(fpi,"         ************************************************************\n \n");
        } 

/*

      fprintf(fpi,"   Options:\n");
      fprintf(fpi,"   ========\n \n");


      if( wgt_type )
        fprintf(fpi,"    %2d) nucleic acid sequences aligned \n", p_count++);
      else
        fprintf(fpi,"    %2d) protein sequences aligned \n", p_count++);

      if( wgt_type == 2 )
        {
          fprintf(fpi,"    %2d) translation",p_count++);
          fprintf(fpi," of nucleotide fragments");
          fprintf(fpi," into peptide fragments\n");
        }

      if( wgt_type == 1 )
        {
          fprintf(fpi,"    %2d) no translation of",p_count++);
          fprintf(fpi," of nucleotide fragments");
          fprintf(fpi," into peptide fragments\n");
        }

      if( wgt_type == 3 )
        {
          fprintf(fpi,"    %2d) mixed alignment consisting", p_count++);
          fprintf(fpi," of P-fragments and N-fragments \n");
        }


      if( seqnum > 2 )
      if( overlap_weights )
        fprintf(fpi,"    %2d) overlap weights used \n", p_count++);
      else
        fprintf(fpi,"    %2d) overlap weights NOT used \n", p_count++);

      if( threshold )
        {
          fprintf(fpi,"    %2d) threshold T =", p_count++);
          fprintf(fpi," %2.2f\n",threshold);
        }

      if( mask )
        {
          fprintf(fpi,"    %2d) non-aligned residues masked", p_count++);
          fprintf(fpi," by `*' \n", p_count++);
        }


      if( lgs_option)
        {
          fprintf(fpi,"    %2d) option for long genomic ", p_count++);
          fprintf(fpi,"sequences used \n");
        }

      if( crick_strand )
        {
          fprintf(fpi,"    %2d) translation of Watson and Crick strand \n", p_count++);
        }

 
 
      if( lmax != MAX_DIA )
        {
          fprintf(fpi,"    %2d) maximum length of fragments = %d", p_count++ , lmax );   
          if( wgt_type == 0)
            fprintf(fpi," residues ");
          if( wgt_type == 1)
            fprintf(fpi," residues ");
          if( wgt_type == 0)
            fprintf(fpi," codons ");
          if( wgt_type == 0)
            fprintf(fpi," codons / residues ");
          fprintf(fpi," \n");
        }



      if( fasta_file )
        fprintf(fpi,"    %2d) separate file in FASTA format \n", p_count++);

      if( msf_file )
        fprintf(fpi,"    %2d) separate file in msf format \n", p_count++);

      if( cw_file )
        fprintf(fpi,"    %2d) separate file in clustal format \n", p_count++);

     if( plot_num )
       {
         fprintf(fpi,"    %2d) %d \"*\" characters", p_count++, plot_num);
         fprintf(fpi," for regions of maximum similarity\n");
       }
*/

 if( online ) { 
   fprintf(fpi,"\n\n    The following options have been used: \n\n") ;
   fprintf(fpi,"     - sequences are");
   if( wgt_type == 0 ) 
     fprintf(fpi," protein sequences \n");     
   if( wgt_type == 1 ) 
     fprintf(fpi," nucleic acid sequences without translation option\n");     
   if( wgt_type == 2 ) 
     fprintf(fpi," nucleic acid sequences with translation option\n");     
   if( speed_optimized ) 
     fprintf(fpi,"     - speed optimized,"); 
     fprintf(fpi," see user guide for details \n"); 
   if( anchors )
     fprintf(fpi,"     - anchor points used\n" ); 
   fprintf(fpi,"\n"); 
 }
 else
   fprintf(fpi,"\n\n   %s \n\n", input_line );

      fprintf(fpi," \n");

      fprintf(fpi,"   Aligned sequences:          length:\n");
      fprintf(fpi,"   ==================          =======\n \n");
       
      for(hv=0;hv<seqnum;hv++)
        {
          fprintf(fpi, " %3d) ", hv + 1 );
          fprintf(fpi, "%s", seq_name[hv] );
          fprintf(fpi, "         %9d\n",seqlen[hv]);
        }
   
   

      fprintf(fpi, "\n   Average seq. length:" );
      fprintf(fpi, "      %9.1f \n", av_len );

      fprintf(fpi,"\n\n   Please note that only upper-case letters are");
      fprintf(fpi," considered to be aligned. \n");

      fprintf(fpi,"\n\n   Alignment (DIALIGN format):\n");
      fprintf(fpi,"   ===========================\n \n");

  }  /* para_print */

