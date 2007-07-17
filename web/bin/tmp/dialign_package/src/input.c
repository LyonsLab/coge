
                 /*******************\
                 *                   *
                 *     DIALIGN 2     *
                 *                   *
                 *     input.c       *
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

extern int max_dia , self_comparison ; 
extern int sim_score[21][21]; 
extern int max_sim_score ;
extern float av_sim_score_pep ;
extern float av_sim_score_nuc ;
extern char par_dir[ NAME_LEN ] ;
extern double **tp400_prot, **tp400_dna, **tp400_trans ; 
extern int *seqlen ;
extern int seqnum ; 


int word_count( char *str ) {

  short word = 0 ; 
  int i ; 
  int word_len = 0 ; 

  for( i = 0 ; i < strlen( str ) - 1  ; i++ ) { 
    if( ( str[i] != ' ' ) && ( str[i] != '\t' ) ) {  
      if( ! word ) { 
        word_len++ ; 
        word = 1 ; 
      }
    }
    else 
      word = 0 ; 

  }
  return( word_len ) ; 
} 


void exclude_frg_read( char *file_name , int ***exclude_list) {

  char exclude_file_name[ NAME_LEN ] ;
  FILE *fp;
  char line[ 10000 ] ;
  int i, len, beg1, beg2, seq1, seq2; 
  int sv = 0, hv, word_num  ;

  strcpy( exclude_file_name , file_name );
  strcat( exclude_file_name , ".xfr" );

  if( (fp = fopen( exclude_file_name, "r")) == NULL)
    erreur("\n\n cannot find file with excluded fragments \n\n\n");

  

  while( fgets( line , MLINE , fp ) != NULL ) {
    if( strlen( line ) > 4 ) {   
      sscanf(line,"%d %d %d %d %d", &seq1, &seq2, &beg1, &beg2 , &len  );

      if( seq1 > seqnum ){
        printf ("\n\n exclueded fragment makes no sense!\n\n");
        printf (" wrong sequence no %d in fragment\n\n", seq1 );
        printf ("%d %d %d %d %d \n\n ", seq1, seq2, beg1, beg2 , len  );
        exit(1) ;
      }
 
      if( seq2 > seqnum ){
        printf ("\n\n    excluded fragment makes no sense!\n\n");
        printf ("    wrong sequence no %d in fragment\n\n", seq2 );
        printf ("    %d %d %d %d %d \n\n", seq1, seq2, beg1, beg2 , len );
        exit(1) ;
      }

/*
      seq1 = seq1 - 1; 
      seq2 = seq2 - 1;
*/

      if( beg1 + len > seqlen[ seq1 - 1 ] + 1 ){
        printf ("\n\n    excluded fragment makes no sense!\n");
        printf ("    fragment");
        printf ("     \" %d %d %d %d %d \"\n", seq1, seq2, beg1, beg2 , len );
        printf ("    doesn't fit into sequence %d:\n", seq1 );
        printf ("    sequence %d has length =  %d\n\n", seq1 , seqlen[ seq1 - 1 ] );
        exit(1) ;
      }


 
      for( i = 0 ; i < len ; i++ ) {  
        exclude_list[ seq1 - 1 ][ seq2 - 1 ][ beg1 + i ] = beg2 + i ;
      }
    }
  }
  
} /* excluded_frg_read  */ 







void ws_remove( char *str ) {
  int pv = 0 ;

  while( ( str[ pv ] == ' ' ) || ( str[ pv ] == '\t' ) ) {
    pv++ ;
  }

  strcpy( str , str + pv );
}

void n_clean( char *str ) {
  int pv = 0 ;
  char *char_ptr ;

  while( ( str[ pv ] == ' ' ) || 
         ( str[ pv ] == '\t' ) || 
         ( str[ pv ] == '>' ) ) {
    pv++ ;
  }
  strcpy( str , str + pv ) ;

  if( ( char_ptr = strchr( str ,' ') ) != NULL)
    *char_ptr = '\0';
  if( ( char_ptr = strchr( str ,'\t') ) != NULL)
    *char_ptr = '\0';
  if( ( char_ptr = strchr( str ,'\n') ) != NULL)
    *char_ptr = '\0';




}


void fasta_test( char *seq_file ) {

  int test = 1;
  int pv = 0;
  FILE *fp;

  char line[ MAX_INPUT_LINE ] ;

  if( (fp = fopen( seq_file , "r")) == NULL) { 
    printf("\n\n Cannot find sequence file %s \n\n\n", seq_file );
    exit(1) ;
  }

  while( test ) {
    fgets( line , MAX_INPUT_LINE , fp );

    ws_remove( line );

    if( line[0] != '\n' )
    if( line[0] == '>' ) 
      test = 0;
    else
      erreur("\n\n  file not in FASTA format  \n\n");
  }
  
  fclose( fp );
}


int seq_read( char *seq_file , char *sq[MAX_SEQNUM] , char **sqn , char **fsqn ) {
 
  char line[ MAX_INPUT_LINE ] ;

  char *nom_seq;
  char *char_ptr;
  int  sn, i, j, k , pv , crc ;
  FILE *fp;
  int max_char[ MAX_SEQNUM ] ;

  if( (fp = fopen( seq_file , "r")) == NULL) { 
    printf("\n\n Cannot find sequence file %s \n\n\n", seq_file );
    exit(1) ;
  }
  fasta_test( seq_file );

  sn = -1 ;
  while( fgets( line , MAX_INPUT_LINE , fp ) != NULL ) {


  ws_remove( line );

    if( line[0] == '>' ) {
      sn++;

   
      n_clean( line );


      fsqn[ sn ] = ( char * ) calloc( strlen( line ) + 3 , sizeof ( char ) ); 

      strcpy( fsqn[ sn ] , line ) ; 


      max_char[ sn ] = 0;
      sqn[ sn ]  = ( char * ) calloc( SEQ_NAME_LEN + 3 , sizeof ( char ) );

      for( crc = 0 ; crc < SEQ_NAME_LEN ; crc++ )  
        if( crc < strlen(line) ) 
          sqn[ sn ][ crc ] =  line[ crc ] ;
        else 
          sqn[ sn ][ crc ] =  ' ';

      sqn[ sn ][ SEQ_NAME_LEN ] = '\0';  

       

    }

 
    else  
      max_char[ sn ] = max_char[ sn ] + strlen( line ) - 1 ; 
  }

  for( i = 0 ; i <= sn ; i++ ) {
    sq[ i ]  = ( char * ) calloc( max_char[ i ] + 1 , sizeof ( char ) );
  }

  if( (seqlen = (int *) calloc( ( sn + 1 ) , sizeof(int) )) == NULL)
    erreur("\n\n problems with memory allocation for `seqlen' \n\n");


  fclose( fp );


  /******************************************/

  if( self_comparison == 1 ) {
    if( sn != 0 ) {
      printf("\n\n With option \"self comparison\" input file must contain one single sequence \n\n" ); 
      exit(1) ;
    }

    sq[ 1 ]  = ( char * ) calloc( max_char[ 0 ] + 1 , sizeof ( char ) );

    sqn[ 1 ]  = ( char * ) calloc( strlen( line ) + 3 , sizeof ( char ) );
    strcpy( sqn[ 1 ] , sqn[ 0 ] ) ;
  }

  /******************************************/


  if( (fp = fopen( seq_file , "r")) == NULL) 
    erreur("\n\n no seq file \n\n");
 
  sn = -1 ;
  while( fgets( line , MAX_INPUT_LINE , fp ) != NULL ) {
    ws_remove( line ); 
    if( line[0] == '>' ) {
      sn++;
      j = 0;
    }
    else
      for( k = 0 ; k < strlen( line )  ; k++ )
      if( 
          ( line[ k ] >= 65 ) && ( line[ k ] <= 90 ) || 
          ( line[ k ] >= 97 ) && ( line[ k ] <= 122 ) 
        ) 
      sq[ sn ][ j++ ] = toupper( line[ k ] ) ; 
  }
  
  sn++; 

  for( i = 0 ; i <  sn ; i++ ) {
    seqlen[ i ] = strlen ( sq[ i ] ) ;
  }

  if( self_comparison ) {
    seqlen[ 1 ] = seqlen[ 0 ] ;  
    for( i = 0 ; i <=  seqlen[ 0 ] ; i++ ) 
      sq[ 1 ][ i ] = sq[ 0 ][ i ] ;    
    sn++; 
  }

  fclose( fp );

  return( sn );
}



void matrix_read( FILE *fp_mat ) {
  int i, j;
  char line[MLINE], dummy[MLINE];
 
  fgets( line , MLINE , fp_mat );
  fgets( line , MLINE , fp_mat );


  for( i = 1 ; i <= 20 ; i++ ) {
    for(j=i;j<=20;j++) {
      fscanf( fp_mat , "%d" , &sim_score[i][j]);
      sim_score[j][i] = sim_score[i][j];  
      if ( sim_score[i][j] > max_sim_score )
        max_sim_score = sim_score[i][j] ;
    }

    fscanf( fp_mat, "%s\n", dummy);
  }

  fclose(fp_mat);

  for( i = 0 ; i <= 20 ; i++ ) {
    sim_score[i][0] = 0 ;
    sim_score[0][i] = 0 ;
  }

/*
 sim_score[0][0] = max_sim_score ;
*/

}
 


void tp400_read( int w_type , double **pr_ptr ) {  
 
  /* reads probabilities from file */
   /* w_type = 0 (protein), 1 (dna w/o transl.), 2 (dna with transl.) */  

  char line[MLINE], file_name[MLINE], suffix[10], str[MLINE] ;
  int sum, len, max_sim, i ;
  double pr;

  FILE *fp;
 
  if ( w_type == 0 ) {
    strcpy( suffix , "prot" );
  }

  if ( w_type == 1 ) {
    strcpy( suffix , "dna" );
  }  

  if ( w_type == 2 ) {
    strcpy( suffix , "trans" );
  }
 
  strcpy( file_name , par_dir ); 
  strcat( file_name , "/tp400_" );
  strcat( file_name , suffix );


 if ( ( fp = fopen( file_name , "r" ) ) == NULL ) { 
   printf("\n\n Cannot find the file %s \n\n", file_name );    
   printf(" Make sure the environment variable DIALIGN2_DIR points\n");
   printf(" to a directory containing the files \n\n");
   printf("   BLOSUM \n   tp400_dna\n   tp400_prot \n   tp400_trans \n\n" );
   printf(" These files should be contained in the DIALIGN package \n\n\n" ) ;
   exit(1) ;
 }


  if ( fgets( line , MLINE , fp ) == NULL ) 
    erreur("\n\n problem with file %s  \n\n", file_name );
  else
    if( w_type % 2 )  
      av_sim_score_nuc = atof( line );
    else
      av_sim_score_pep = atof( line );
     

  while( fgets( line , MLINE , fp ) != NULL )
   {
      sscanf(line,"%d %d %s", &len, &sum, str  );

      pr = atof(str);
      pr_ptr[len][sum] = pr;

    }


}    /*  tp400_read  */




