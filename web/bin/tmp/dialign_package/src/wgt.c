

                 /*******************\
                 *                   *
                 *     DIALIGN 2     *
                 *                   *
                 *     wgt.c         *
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

extern int max_sim_score, istep, *seqlen ;
extern float **min_weight , threshold , av_sim_score_pep , av_sim_score_nuc ;
extern int  max_dia , wgt_type ;
extern double **tp400_prot,  **tp400_dna, **tp400_trans; 
extern float **wgt_prot, **wgt_dna, **wgt_trans;
extern char par_dir[NAME_LEN];

void rel_wgt_calc( int l1 , int l2 , float **rel_wgt )
  {
    int  l , m, mss;
    float ent, factor, l1f, l2f, **wgt , av_sim_score ;
    double t_pr, pr400, **tpr ;

/*
    printf(" it %d, rel_wgt_calc: len = %d , %d \n", istep , l1 , l2 );
*/

    if( rel_wgt == wgt_prot ) {
      tpr = tp400_prot ;
      mss = max_sim_score ;
      av_sim_score = av_sim_score_pep ; 
    } 
    
    if( rel_wgt == wgt_dna ) {
      tpr = tp400_dna ;
      mss = 1 ;
      av_sim_score = av_sim_score_nuc ; 
    }
      
    if( rel_wgt == wgt_trans ) {
      tpr = tp400_trans ;
      mss = max_sim_score ;
      av_sim_score = av_sim_score_pep ; 
    }
      
  

    l1f = l1;
    l2f = l2;

    factor = ( l1f * l2f ) / 400.00;


    for( l = 1 ; l <= max_dia           ; l++ )
    for( m = 0 ; m <= l * mss ; m++ )
      {
        rel_wgt[l][m] = 0;


        if( tpr[l][m] )
        if( m > av_sim_score * l )

          {
            pr400 = tpr[l][m];

            if( pr400 > 0.0000000001 )
              t_pr = 1 - pow( 1 - pr400 , factor );
            else
              t_pr = pr400 * factor;

            ent = 0;

            if(t_pr)
              ent = -log( t_pr );

            if( ent > threshold )
              rel_wgt[l][m] = ent;
          }
      } 
  }  /*  rel_wgt_calc */




void wgt_prnt_prot( ) {
  int  i, j ; 
  printf(" \n\n  weight scores for PROTEIN fragments\n\n" ); 
  printf("  sequence lengths = %d , %d \n\n", seqlen[0] , seqlen[1] ) ;  
  for( i = 1 ; i <= max_dia ; i++ ) {
    for( j = 0 ; j <= ( i * 15 ) ; j++ ) 
      printf(" %3d %3d %f  \n", i , j , wgt_prot[ i ][ j ] );
  }
}

void wgt_prnt_dna( ) {
  int  i, j ; 
  printf(" \n\n  weight scores for NON-TRANSLATED DNA fragments\n\n" ); 
  printf("  sequence lengths = %d , %d \n\n", seqlen[0] , seqlen[1] ) ;  
  for( i = 1 ; i <= max_dia ; i++ ) {
    for( j = 0 ; j <= i ; j++ ) 
      printf(" %3d %3d %f  \n", i , j , wgt_dna[ i ][ j ] );
  }
}


void wgt_prnt_trans( ) {
  int  i, j ; 
  printf(" \n\n  weight scores for TRANSLATED DNA fragments\n\n" ); 
  printf("  sequence lengths = %d , %d \n\n", seqlen[0] , seqlen[1] ) ;  
  for( i = 1 ; i <= max_dia ; i++ ) {
    for( j = 0 ; j <= ( i * 15 ) ; j++ ) 
      printf(" %3d %3d %f  \n", i , j , wgt_trans[ i ][ j ] );
  }
}



void wgt_prnt( ) {
  if (wgt_type == 0 )  
    wgt_prnt_prot( );

  if (wgt_type % 2 )  
    wgt_prnt_dna( );

  if (wgt_type > 1 )  
    wgt_prnt_trans( );
}



void mem_alloc( ) {
      /* allocates memory for `tp400_xxx', `wgt_xxx' */ 

  int i;
 
  if( wgt_type == 0 ) {
    if( (tp400_prot = (double **) calloc( ( max_dia + 1 ) , sizeof(double*) )) 
      == NULL) { 
      printf(" problems with memory allocation for `tp400_prot' !  \n \n");
      exit(1);
    }

    if( ( wgt_prot = (float **) calloc( (max_dia+1) , sizeof(float*) ))
      == NULL) {
      printf(" problems with memory allocation for `weights' !  \n \n");
      exit(1);
    }
  }

  if( wgt_type % 2 ) { 
    if( (tp400_dna = (double **) calloc( ( max_dia + 1 ) , sizeof(double*) ))
      == NULL) {
      printf(" problems with memory allocation for `tp400_dna' !  \n \n");
      exit(1);
    }

    if( ( wgt_dna = (float **) calloc( (max_dia+1) , sizeof(float*) ))
      == NULL) {
      printf(" problems with memory allocation for `weights' !  \n \n");
      exit(1);
    }
  }

 
  if( wgt_type > 1 ) {
    if( (tp400_trans = (double **) calloc( ( max_dia + 1 ) , sizeof(double*) ))
      == NULL) {
      printf(" problems with memory allocation for `tp400_trans' !  \n \n");
      exit(1);
    }
 
    if( ( wgt_trans = (float **) calloc( (max_dia+1) , sizeof(float*) ))
      == NULL) {
      printf(" problems with memory allocation for `weights' !  \n \n");
      exit(1);
    }
  }
  
  for( i = 1 ; i <= max_dia ; i++ ){
   
 
    if( wgt_type == 0 ) {
      if( (tp400_prot[i] = 
        (double *) calloc( ((i + 1) * max_sim_score ) , sizeof(double) )) 
        == NULL) { 
          printf(" problems with memory allocation for `tp400_prot' !  \n \n");
          exit(1);
      }

      if( (wgt_prot[i] =
        (float *) calloc( ((i+1) * max_sim_score ) , sizeof(float) ))
        == NULL) {
           printf(" problems with memory allocation for `weights'!\n\n");
           exit(1);
      }
    }
  

    if( wgt_type % 2 ) {
      if( (tp400_dna[i] =
        (double *) calloc( ((i + 1) ) , sizeof(double) ))
        == NULL) {
          printf(" problems with memory allocation for `tp400_dna' !  \n \n");
          exit(1);
      }

      if( (wgt_dna[i] =
        (float *) calloc( ((i+1) ) , sizeof(float) ))
        == NULL) {
           printf(" problems with memory allocation for `weights'!\n\n");
           exit(1);
      }
    }


    if( wgt_type > 1 ) {
      if( (tp400_trans[i] =
        (double *) calloc( ((i + 1) * max_sim_score ) , sizeof(double) ))
        == NULL) {
          printf(" problems with memory allocation for `tp400_trans' %d !  \n \n", i);
          exit(1);
      }

      if( (wgt_trans[i] =
        (float *) calloc( ((i+1) * max_sim_score ) , sizeof(float) ))
        == NULL) {
           printf(" problems with memory allocation for `weights'!\n\n");
           exit(1);
      }
    }
  }     
}      /*  void memory_allocation  */







