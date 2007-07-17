
                 /*******************\
                 *                   *
                 *     DIALIGN 2     *
                 *                   *
                 *     regex.c       *
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

extern float max_mot_offset , mot_offset_factor , mot_factor ;
extern int char_num[ MAX_REGEX ] ;
extern int regex_len , mot_len ; 
extern char *mot_char[ MAX_REGEX ] ;
extern int seqnum, *seqlen ;
extern char *seq[MAX_SEQNUM];
extern short **mot_pos ;
extern FILE *fp_mot ; 

void regex_complain( char *regex ) {
  printf("\n   bracket structure in regular expression makes no sense \n");
  printf("\n          %s  \n\n", regex) ; 
  printf("   program terminated\n\n"); 
  exit(1);  
}

void struc_check( char *regex ) {
  int p, bracket_count = 0 ;

  for( p = 0 ; p < MAX_REGEX; p++ ) { 
    char_num[ p ] = 0 ;
  }

  for( p = 0 ; p < regex_len ; p++ ) {

    if( regex[ p ] == '[' ) 
      bracket_count++ ; 
       
    if( ( regex[ p ] != '[' ) && ( regex[ p ] != ']' ) ) { 
      char_num[ mot_len ]++ ;  
      regex[ p ] = toupper( regex[ p ] ) ; 
    }

    if( regex[ p ] == ']' ) 
      bracket_count-- ; 

    if( ( regex[ p ] == ']' ) || ( bracket_count == 0 ) )  
      mot_len++ ; 


    if( ( bracket_count < 0 ) || ( bracket_count > 1 ) )   
      regex_complain( regex ) ; 


  }

  if( bracket_count != 0 ) 
    regex_complain( regex ) ; 

}
 

void regex_parse( char *mot_regex ) {

  int i, p , mp = 0 ; 
  int in_bracket = 0;  
  int char_c = 0 ; 


  if( ( mot_pos = ( short  ** ) calloc( seqnum , sizeof( short *) ) ) == NULL) {
    printf(" problems with memory allocation");
    printf(" for `mot_pos' !  \n \n");
    exit(1);
  }

  for( i = 0 ; i < seqnum ; i++ )
    if( ( mot_pos[i] = ( short *) calloc( ( seqlen[i] + 2 ) , sizeof( short ) ) ) == NULL) {
    printf(" problems with memory allocation");
    printf(" for `mot_pos[%d]' !  \n \n", i);
    exit(1);
  }




  struc_check( mot_regex ) ; 

/*
  printf("  \n  regex_len = %d\n", regex_len) ; 
  printf("  mot_len = %d\n", mot_len) ; 
  printf("\n"); 

    for( p = 0 ; p < mot_len ; p++ ) {
      printf("  %d ", char_num[ p ] );
    }
    printf("\n\n");  
*/
 
  for( p = 0 ; p < mot_len ; p++ ) {
    mot_char[ p ] = (char *) calloc( char_num[ p ] , sizeof(char) );
  }


  /* PROBLEM */ 


  for( p = 0 ; p < regex_len ; p++ ) {

    if( mot_regex[ p ] == '[' ) { 
      in_bracket = 1 ; 
    }

    if( mot_regex[ p ] == ']' ) {  
      in_bracket = 0 ; 
      char_c = 0 ; 
      mp++ ;
    }

    if( ( mot_regex[ p ] != '[' ) && ( mot_regex[ p ] != ']' ) ) {  /* char */ 
      if( in_bracket ) {
        mot_char[ mp ][ char_c ] = mot_regex[ p ] ; 
        char_c++;
      }
      else {     /* not in bracket */                   
        char_c = 0 ; 
        mot_char[ mp ][ 0 ] = mot_regex[ p ] ; 
        mp++ ; 
      }
    }
  }

/*
  for( mp = 0 ; mp < mot_len ; mp++ ) {
    printf("  position %d   ", mp + 1 ); 
    for( p = 0 ; p < char_num[ mp ] ; p++ ) {
      printf(" %c ", mot_char[ mp ][ p ]  ) ;
    }
    printf("\n"); 
  }
*/ 

}
 
seq_parse( char *mot_regex ) { 
  int sn, ok , i ; 
  int sp, ap, rp, hv, match;
  max_mot_offset = sqrt ( - log ( 0.1 ) *  10 / mot_factor ) * mot_offset_factor; 


  for( sn = 0 ; sn < seqnum ; sn++ ) 
  for( sp = 0 ; sp < ( seqlen[ sn ] - mot_len + 1 ) ; sp++ ) { 
    ok = 1 ;
    rp = 0 ;  
    while( ok && ( rp < mot_len ) ) {
      if( mot_char[ rp ][ 0 ]  != 'X' ) { 
        match = 0 ;  
        for( hv = 0 ; hv < char_num[ rp ] ; hv++ ) { 
          if( mot_char[ rp ][ hv ] == seq[ sn ][ sp + rp ] ) { 
            match = 1 ; 
          }
        }
      }
      ok = match ; 
      rp++;    
    }
    if( ok ) { 
      printf( " motif in seq %d at pos %d  \n", sn + 1 , sp + 1 ) ;   
      mot_pos[ sn ][ sp + 1 ] = 1 ;  
    }
    else 
      mot_pos[ sn ][ sp + 1 ] = 0 ;  
  }


  printf("\n") ; 

/*
  for( sn = 0 ; sn < seqnum ; sn++ ) { 
        printf("     %s \n", seq[ sn ] ) ;
        printf("     "); 
    for( i = 1 ; i <= seqlen[ sn ] ; i++ ) { 
     

      if( mot_pos[ sn ][ i ]  ) 
        printf("*");
      else 
        printf(" ");
    }
    printf("\n\n" ) ; 
  }
  printf("\n" ) ; 
*/ 

}


void regex_format_complain() { 
  printf("\n \n   Arguments in command line don't make sense! \n");
  printf("   (Motifs not properly specified) \n \n");
  printf("   With the motif-search option, the program call is:\n\n");
  printf("      ./dialign2-2 [para] -mot <regex> <fct1> <fct2> ");
  printf("[para] <seq> \n\n");
  printf("   where \n      <regex>  is a regular expression,");
  printf(" e.g. \"AT[CG]XT\",\n");
  printf("      <fct1>    is a weighting factor \n");
  printf("      <fct2>    is a weighting factor \n");
  printf("      <seq>    is the input sequence file and \n");
  printf("      [para]   are (optional)");
  printf(" additional program parameters\n\n" );
  exit(1);
}
 
float mot_dist_factor( int offset , float parameter ) {
  float mdf , parameter2, factor1 ;
  int offset2 ; 

  offset2 = offset * offset ;
  parameter2 = parameter * parameter ;

  factor1 = (float) offset2 / ( parameter2 * 10 ) ;
  mdf = exp( - ( offset2 ) / ( parameter2 * 10 ) ) ; 

  return mdf ;

}
 
