
                 /*******************\
                 *                   *
                 *     DIALIGN 2     *
                 *                   *
                 *     anchor.c      *
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


extern int  anc_num, *seqlen ;
extern int seqnum;
extern char *seq[MAX_SEQNUM];
extern struct multi_frag *anchor_frg ;

void anchor_check( int s1, int s2, int b1, int b2, int l ,  float scr ) {  

      if(  
        ( s1 < 1 ) || 
        ( s1 > seqnum ) 
      ) { 
        printf(" \n\n  wrong sequence # %d in anchoring file\n\n", s1 );
        printf("  data set consists only of %d sequences \n\n", seqnum );
        printf("  PROGRAM TERMINATED \n\n" ) ;
        exit( 1 ) ; 
      }  
      if(  
        ( s2 < 1 ) || 
        ( s2 > seqnum )  
      ) { 
        printf(" \n\n  wrong sequence # %d in anchoring file\n\n", s2 );
        printf("  data set consists only of %d sequences \n\n", seqnum );
        printf("  PROGRAM TERMINATED \n\n" ) ;
        exit( 1 ) ; 
      }


      if( s1 == s2 ) { 
        printf("\n strange data in anchoring file:\n");
        printf(" sequence # %d anchored with itself.\n\n", s1 );
        printf("  PROGRAM TERMINATED \n\n" ) ;
        exit(1) ; 
      }
 


/*
      if(  
        ( b1 < 1 ) || 
        ( b1 + l - 1 > seqlen[ s1 - 1 ] )  
      ) { 
        printf(" \n\n anchor # %d starts", anc_num + 1 ) ;
        printf(" at position %d in sequence %d and has a length of %d.\n", b1, s1, l ) ;
        printf(" This does not fit into sequence # %d " , s1 );
        printf(" (sequence length = %d) \n\n", seqlen[ s1 - 1 ] ) ; 
        printf("  PROGRAM TERMINATED \n\n" ) ;
        exit( 1 ) ; 
      } 
*/

      if( 
        ( b1 < 1 ) ||
        ( b1 + l - 1 > seqlen[ s1 - 1 ] )
      ) {
        printf(" \n\n  WARNING:"); 
        printf(" \n\n  anchor # %d starts", anc_num + 1 ) ;
        printf(" at position %d in sequence %d\n ", b1, s1 ) ;
        printf(" and is %d residues in length.\n", l ) ;
        printf("  However, sequence %d" , s1 );
        printf(" is only %d residues in length \n\n", seqlen[ s1 - 1 ] ) ;
        printf("  PROGRAM TERMINATED \n\n" ) ;
        exit( 1 ) ;
      }

      if( 
        ( b2 < 1 ) ||
        ( b2 + l - 1 > seqlen[ s2 - 1 ] )
      ) {
        printf(" \n\n  WARNING:"); 
        printf(" \n\n  anchor # %d starts", anc_num + 1 ) ;
        printf(" at position %d in sequence %d\n ", b2, s2 ) ;
        printf(" and is %d residues in length.\n", l ) ;
        printf("  However, sequence %d" , s2 );
        printf(" is only %d residues in length \n\n", seqlen[ s2 - 1 ] ) ;
        printf("  PROGRAM TERMINATED \n\n" ) ;
        exit( 1 ) ;
      }


}


int multi_anc_read( char *file_name ) {

  char anc_file_name[ NAME_LEN ] ;
  FILE *fp;
  struct multi_frag *current_frg ;
  char line[ 10000 ] ;
  int i, len, beg1, beg2, sv = 0, wrdl, hv, word_num  ;
  int seq1, seq2 ;
  float wgt; 

  strcpy( anc_file_name , file_name );
  strcat( anc_file_name , ".anc" );

  if( (fp = fopen( anc_file_name, "r")) == NULL)
    erreur("\n\n cannot find file with anchor points \n\n\n");

    if( ( anchor_frg = ( struct multi_frag * ) calloc( 1 , sizeof( struct multi_frag ) ))
      == NULL) {
      printf(" problems with memory allocation for `anchor fragments' !  \n \n");
      exit(1);
    }

  current_frg = anchor_frg ; 


  while( fgets( line , MLINE , fp ) != NULL ) {

    if(  word_count( line ) == 6  ) {   
      sscanf(line,"%d %d %d %d %d %f ", &seq1 , &seq2 , &beg1, &beg2 , &len , &wgt );

      anchor_check( seq1 , seq2 , beg1, beg2 , len , wgt ) ; 
 
      seq1 = seq1 - 1 ; 
      seq2 = seq2 - 1 ; 

      current_frg->s[0] = seq1 ;
      current_frg->s[1] = seq2 ;
      current_frg->b[0] = beg1 ;
      current_frg->b[1] = beg2 ;
      current_frg->ext  = len ;
      current_frg->weight  = wgt;


      current_frg->next = (struct multi_frag *)
                       calloc( 1 , sizeof(struct multi_frag) );

      current_frg = current_frg->next;
      anc_num++; 
    }
    else { 
      if( word_count( line ) != 0 ){
        printf("\n\n  Anchor file has wrong format. ");
        printf("\n  Each line must contain 6 numbers! \n");
        printf("\n  Anchor file contains line \n\n");
        printf("         %s \n", line); 
        printf("  PROGRAM TERMINATED \n\n" ) ;
        exit(1) ;
      } 
    }
  }
} /* multi_anc_read  */ 



