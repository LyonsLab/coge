
                 /************************\
                 *                        *
                 *     DIALIGN 2.2.1      *
                 *                        *
                 *       dialign.c        * 
                 *                        *
                 *       written by       *
                 *                        *
                 *    B. Morgenstern      *
                 *                        *
                 \************************/




#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "dialign.h"
#include "define.h"
#include "alig_graph_closure.h"



FILE *fp_dia, *fp_dpa, *fp_frg , *fp_mot ; 
struct multi_frag *anchor_frg ;

int col_score = 0; 
int char_num[ MAX_REGEX ] ;
char *mot_char[ MAX_REGEX ] ; 
int regex_len , mot_len = 0 ; 


clock_t beg_pa , end_pa , beg_ali , end_ali , beg_ts , end_ts ;
float time_diff_pa , time_diff_ali , perc_pa_time , time_diff_srt ; 
float total_pa_time = 0 ; 


float mot_factor , mot_offset_factor , max_mot_offset ; 

int wgt_type_plot = 0 , motifs = 0 ; 
int bubblesort = 0 , cd_gobics = 0 ; 
int nas = 0 , ref_seq = 0 , i_max ; 
int speed_optimized = 0 ; 
int online = 0 ; 
int time_stamps = 0 ; 
int break1 = 0 ; 
int break2 = 0 ; 
int wgt_print = 0 ; 
int wgt_print_x = 0 ; 
short max_itnum = MAX_ITNUM ; 
int quali_num = 1 ; 
int wgt_plot = 0 ; 
int self_comparison = 0;
short exclude_frg = 0; 
int ***exclude_list ;
int max_sim_score = -2000 ; 
int sf_mat = 0 ; 
char nuc1, nuc2, nuc3 ;
short crick_strand = 0;
int frg_count = 0; 
int dna_speed = 0;
char pst_name[NAME_LEN];
int cont_it = 1 , wgt_type = 0  ;
int mask = 0, strict = 0 , textual_alignment = 1;
char prn[ NAME_LEN ] ;
int redundant, print_max_nd = 1;
int lmax = MAX_DIA;
char **arguments;
int  pr_av_nd = 0, pr_av_max_nd ;
char input_line[ NAME_LEN ];
char input_parameters[ NAME_LEN ];
int print_status = 0 ;
char clust_sim[NAME_LEN] ;
float tot_weight = 0, av_len;
int anchors = 0;
int pa_only = 0;
int dia_num = 0;
int max_dia_num = 0;
float av_dia_num = 0;
float av_max_dia_num = 0;
int afc_file = 0;
int afc_filex = 0;
int dia_pa_file = 0;
int frag_file = 0;
int argnum;
int standard_out = 0;
int plot_num = 4 ;
int default_name = 1;
int fasta_file = 0;
int cw_file = 0; 
int msf_file = 0;
char *upg_str;
int dcount = 0;


int **shift; 
int   thr_sim_score = 4 ;
char *seq[MAX_SEQNUM];   /* sequences */
char *newseq[MAX_SEQNUM];   /* sequences */
int sim_score[21][21];  /* similarity matrix */
float av_sim_score_pep ;
float av_sim_score_nuc ;
float **glob_sim;        /* overall similarity between any two sequences */
float **wgt_prot  ;      /* `weight' of diagonals */
float **wgt_dna   ;      /* `weight' of diagonals */
float **wgt_trans ;      /* `weight' of diagonals */
float **min_weight;      /* `weight' of diagonals */
int min_dia = MIN_DIA ;             /* minimum length of diagonals */
int max_dia = MAX_DIA ;  /* maximum length of diagonals */
int iter_cond_prob = 0;
int *seqlen;                /* lengths of sequences */
char *full_name[MAX_SEQNUM] ;
float **pair_score;
short **cont_it_p; 
float score;
int maxlen;              /* maximum length of sequences */
int seqnum;              /* number of sequences */
int *num_dia_bf;         /* num_dia_bf[ istep ] = number of diagonals from
                            all pairwise alignments BEFORE FILTER
                            PROCEDURE in iteration step `istep' */     
int *num_dia_af;         /* num_dia_af[istep] = number of diagonals from
                            all pairwise alignments AFTER FILTER 
                            PROCEDURE in iteration step `it' */     
int num_dia_anc;         /* number of diagonals definde by anchored 
                            regions */
int num_all_it_dia = 0;  /* total number of diagonals in multiple alignment 
                            in all iteration steps */
float weight_sum_bf;     /* sum of weights  of diagonals in multiple 
                            alignment before filter procedure */  
float weight_sum_af;     /* sum of weights  of diagonals in multiple 
                            alignment after fliter procedure*/
float threshold = 0.0 ;  /* threshold T */
int num_dia_p;           /* number of diagonals in pairwise alignment */ 
int long_output = 0;     /* if long_output = 1, a log-file is produced.  */   
int frg_mult_file = 0 ; 
int frg_mult_file_v = 0 ; 
int overlap_weights = 1 ;  
int ow_force = 0 ;
int anc_num = 0;          /* number of anchored regions 
                            (specified in file *.anc) */
int par_count;           /* number of parameters       */
float pairalignsum;      /* sum of weights in pairwise alignment */ 
int pairalignlen;        /* sum of aligned residues in pairwise alignment */
char amino_acid[22];
int istep;  
struct multi_frag        /* pointer to first diagonal in multiple alignment */
      *this_it_dia;      /* in current iteration step */  
struct multi_frag        /* pointer to first diagonal in multiple alignment */
      *all_it_dia;       /* in all iteration step */
struct multi_frag *end_dia;  
                         /* pointer to last diagonal in multiple alignment */

char par_dir[NAME_LEN];
char *seq_name[MAX_SEQNUM];
char mat_name[NAME_LEN];         /* name of file containing similarity matrix */
char mat_name_p[NAME_LEN];
char anc_name[NAME_LEN];  /* anchored regions */
char seq_file[NAME_LEN];
char input_name[NAME_LEN];
char tmp_str[NAME_LEN];
char output_name[NAME_LEN];
char printname[NAME_LEN];
char mot_regex[MAX_REGEX] ; 

char *par_file;

short **mot_pos ;       /* positions of pre-defined motifs */ 

int **amino;           /* amino acid residues in protein sequences or 
                          translated DNA sequences, respective */

int **amino_c;         /* amino acid residues on crick strand */ 
 
CLOSURE *clos;         /* closure data structure for GABIOS-LIB */

int ***open_pos;           /* open_pos[i][j][p] = 1, if the p-th residue of 
                          sequence i is not yet directly (by one diagonal) 
                          aligned with any residue of sequence j and 
                          open_pos[i][j][r] = 0 otherwise. So, at the
                          beginning of the first iteration step, all values 
                          are 1. In the subsequent iteration steps,
                          only those parts of the sequence are considered,  
                          that are not yet aligned. */     

  
struct multi_frag *pair_dia;   /* diagonals in pairwise alignemnt */


double **tp400_prot ;    /* propability distribution for sums of similarity
                       socores in diagonals occurring in comparison matrix
                       (by random experiments and approximation  */

double **tp400_dna ;    /* propability distribution for sums of similarity
                       socores in diagonals occurring in comparison matrix
                       (by random experiments and approximation  */

double **tp400_trans ;    /* propability distribution for sums of similarity
                       socores in diagonals occurring in comparison matrix
                       (by random experiments and approximation  */


char dia_pa_name[NAME_LEN];
char frag_file_name[NAME_LEN];
char mot_file_name[NAME_LEN];


/********************************/
/* prototypes                   */
/********************************/

 extern float mot_dist_factor ( int offset , float parameter ) ;
 extern int word_count( char *seq ) ; 
 extern void subst_mat(char *file_name, int fragno , struct  multi_frag *smp );
 extern int seq_read( char *in_file , char *sq[MAX_SEQNUM] , char **sqn , char **fsqn) ;
 extern int anc_read( char *file_name ) ;
 extern int multi_anc_read( char *file_name ) ;
 extern void randomize( int r_numb , FILE *fp1 );
 extern int mini2(int a, int b);
 extern int maxi2(int a, int b);
 extern int mini3(int a, int b, int c);
 extern int num_test( char *cp );
 extern void mini(int *a, int b);
 extern void maxi(int *a, int b);
 extern void filter( int *num, struct multi_frag *vector );
 extern void throw_out( float *weight_sum );
 extern void sel_test();
 extern float frag_chain( int n1 , int n2 , FILE *fp , FILE *fp2, int *num );
 extern void para_read( int num , char **arg ); 
 extern void frag_sort( int number , struct multi_frag *dp , int olw );
 extern void ow_frag_sort( int number , struct multi_frag *dp , int olw );
 extern void bubble_sort( int number , struct multi_frag *dp );
 extern void ow_bubble_sort( int number , struct multi_frag *dp );
 extern void seq_shift();
 extern int translate(char c1, char c2, char c3, int s , int i);
 extern char invert( char c1 ) ;
 extern int int_test(float f);
 extern int match_test( struct multi_frag *dia, int mn);
 
 extern void para_print(char *s_f, FILE *f);
 extern void ali_arrange(int fragno , struct  multi_frag *smp, FILE *fp, FILE *fp2, FILE *fp3 , FILE *fp4 , FILE *fp_csc );
 extern void print_log( struct multi_frag *d , FILE *fp_l , FILE *fp_fs);
 extern void print_fragments( struct multi_frag *d , FILE *fp_frg );
 extern void tp400_read( int wgt_type , double **pr_ptr );
 extern void ow_add(struct  multi_frag *sm1 , struct  multi_frag *sm2);
 extern void av_tree_print();
 extern void matrix_read( FILE *fp_mat ) ;
 extern void mem_alloc( ) ;
 

                    /******************************/
                    /*           main             */
                    /******************************/



main(int argc, char **argv)
{
 int k,  anc1, dia_counter, tmpi1, tmpi2 ;

 struct multi_frag *current_dia, *diagonal1, *diagonal2, *anc_dia;  
                        /* pointers to diagonals in multiple alignment */ 

 char str[NAME_LEN], dist_name[NAME_LEN]; 
 char par_str[NAME_LEN];  
 char *char_ptr;
 char prn2[NAME_LEN];
 char logname[NAME_LEN];
 char fsm_name[NAME_LEN];
 char dia_name[NAME_LEN];
 char csc_name[NAME_LEN];
 char itname[NAME_LEN], itname2[NAME_LEN], itname3[NAME_LEN];
 char itname4[NAME_LEN];
 char dialign_dir[NAME_LEN];

 int i, j, hv, sv, fv; 


 FILE *fp_ali, *fp2, *fp3, *fp4, *fp_log, *fp_fsm, *fp_st , *fp_csc ; 
 FILE *fp_matrix ;               /* file containing similarity matrix */

 strcpy(mat_name,MATNAME);
 strcpy( clust_sim , "av" );
 
 par_file = (char *) calloc((size_t) NAME_LEN , sizeof(char) );


 if( time_stamps ) 
   beg_ali = clock() ; 


  par_file = "/opt/apache/CoGe/bin/dialign2_dir/";

 /* This stuff was commented out and hardcoded above by Eric and Josh
  strcpy ( dialign_dir , "DIALIGN2_DIR" );
  if ((par_file = getenv(dialign_dir)) == NULL)
    {
      printf("\n \n \n    Please set the environmentvariable DIALIGN2_DIR \n");
      printf("    as described in the README file \n");
      exit(1);
    }
 */

 argnum = argc;

 strcpy( par_dir , par_file );

if(argc == 1)
  {
    printf("\n    usage: %s [ options ] <seq_file> \n\n", argv[0] );
    printf("    <seq_file> contains input sequences in FASTA format.\n"); 
    printf("    Per default, sequences are assumed to be protein sequences.\n" ) ;
    printf("    For DNA alignment, please use one of these options: \n\n");
    printf("     -n    DNA sequences; similarity calculated at the nucleotide level \n\n"); 
    printf("     -nt   DNA sequences; similarity calculated at the peptide level\n");
    printf("           (by translation using the genetic code) \n\n");
    printf("     -lgs  long genomic sequences: Both nucleotide and peptide\n");
    printf("           similarities calculated \n\n");  
    printf("    Many more options are available, please consult the \n");
    printf("    DIALIGN USER_GUIDE that should come with the DIALIGN package.\n");
    printf("    For more information on DIALIGN, please visit the DIALIGN\n"); 
    printf("    home page at BiBiServ (Bielefeld Bioinformatic Server): \n\n") ;
    printf("        http://bibiserv.techfak.uni-bielefeld.de/dialign/ \n\n");  
    printf (par_dir);  
    exit(1) ;
  }

 arguments = ( char ** ) calloc( argnum , sizeof ( char * ) );

 for( i = 0 ; i < argnum ; i++ )
   {
     arguments[i] = ( char *)  calloc( NAME_LEN , sizeof (char) );
     strcpy( arguments[i] , argv[i] );
   }
 


 strcpy( input_name , argv[ argc - 1 ] );
  
 threshold = 0.0 ;


 para_read( argnum , arguments );

 if( ( textual_alignment == 0 ) && ( col_score == 1 ) ) { 
   printf("\n\n   Option -csc makes sense only if \"textual alignment\"");
   printf(" is produced. \n");
   printf("   This can be enforced with option -ta \n\n");
   printf("   program terminated \n\n\n");
   exit(1) ;
 } 


 if( cd_gobics ) {
 strcpy( input_line , "program parameters:  " ) ; 
 for( i = 1 ; i < ( argnum -1 ) ; i++ ) {
     strcat( input_line , argv[i] );
     strcat( input_line , " " );
   }
 }
 else {
 strcpy( input_line , "program call:  " ) ; 
 for( i = 0 ; i < argnum ; i++ ) {
     strcat( input_line , argv[i] );
     strcat( input_line , " " );
   }
 }


 if ( wgt_type > 0 )  
   strict = 1 ; 

 strcpy( seq_file , input_name );

 if(
        ( ! strcmp( input_name + strlen( input_name ) - 4 , ".seq" ) )
     || ( ! strcmp( input_name + strlen( input_name ) - 3 , ".fa" ) )
     || ( ! strcmp( input_name + strlen( input_name ) - 6 , ".fasta" ) )
   )
 if( ( char_ptr = strrchr(input_name,'.') ) != NULL)
   *char_ptr = '\0';


 strcpy( anc_name , input_name );
 strcat( anc_name , ".anc" );

 seqnum = seq_read( seq_file , seq , seq_name , full_name ) ;

 if ( motifs )
   regex_parse( mot_regex ) ; 


 if( ( seqnum == 2 ) && ( iter_cond_prob == 0 ) ) 
   max_itnum = 1 ; 

 
     if(  ( ow_force == 0 ) && ( seqnum > OVERLAP_THRESHOLD )  )
       overlap_weights = 0;
     if( seqnum == 2 )
       overlap_weights = 0;

  if( seqnum < 2 ) { 

    if( cd_gobics ) {
      printf("\n\n         Something is wrong with your sequence file. Maybe you entered a\n");
      printf("         MS WORD or RFT file or your file contains only one single sequence.\n");
      printf("         Please note that our server only accepts plain text files. \n\n");  
      printf("         For more information, please consult our online manual \n");
      printf("         at the CHAOS/DIALIGN home page:\n\n");  
      printf("             http://dialign.gobics.de/chaos-dialign-manual");
    }

    else { 
      printf("\n\n         Your sequence file containes only a single sequence.\n");
      printf("         Please make sure your input file contains at least two sequences.\n\n");
      printf("         For more information, please consult the online manual \n");
      printf("         at the DIALIGN home page: \n\n");
      printf("             http://bibiserv.techfak.uni-bielefeld.de/dialign/manual.html ");
    }



    printf("\n       \n \n \n \n");
    exit(1);
  }

  maxlen = 0;

  

  if( (pair_score = (float **) calloc( seqnum , sizeof(float *) )) == NULL)
       {       
           printf(" problems with memory allocation for `pair_score' !  \n \n");
           exit(1);
       }

  for(i=0;i<seqnum;i++)
  if( (pair_score[i] = (float *) calloc( seqnum , sizeof(float) )) == NULL)
       {       
           printf(" problems with memory allocation for `pair_score' !  \n \n");
           exit(1);
       }


  if(( cont_it_p = (short **) calloc( seqnum , sizeof( short *))) == NULL ){
    printf(" problems with memory allocation for `cont_it_p ' !  \n \n");
    exit(1);
  }

  for( i = 0 ; i < seqnum ; i++ )  
  if( (cont_it_p[i] = (short *) calloc( seqnum , sizeof(short) )) == NULL) {   
    printf(" problems with memory allocation for `cont_it_p' !  \n \n");
    exit(1);
  }

  for( i = 0 ; i < seqnum ; i++ )
  for( j = 0 ; j < seqnum ; j++ )
    cont_it_p[i][j] = 1 ; 





  for( i = 0 ; i < seqnum ; i++ )
   {
    av_len = av_len + seqlen[i];

    if( seqlen[i] == 0 )
      {
        printf("\n \n \n                       WARNING: \n \n");
        printf("          Sequence %d contains no residues.\n",i+1);
        printf("          Please inspect the sequence file.\n \n ");
        printf("\n \n          Program terminated \n \n \n " );     

        exit(1);
      }
 
    if(maxlen < seqlen[i])
       maxlen = seqlen[i];
   }

  av_len = av_len / seqnum;

  if ( motifs )
    seq_parse( mot_regex ) ; 
  
  seq_shift();


       if( (glob_sim = 
           (float **) calloc( seqnum , sizeof(float*))) == NULL) 
            {
                printf("Problems with memory allocation for glob_sim\n"); 
                exit(1); 
            } 

   for(i=0;i<seqnum;i++)
     {

       if( (glob_sim[i] = 
       (float *) calloc( seqnum , sizeof(float))) == NULL) 
         { 
           printf("Problems with memory allocation for glob_sim \n"); 
           exit(1); 
         } 
 
      }

   strcpy(par_str,"sdfsdf");

   if( argc > 1 )
   {
   strcpy(str,par_dir);
   strcat(str,"/");
   strcat(str,mat_name);
   strcpy(mat_name_p,str);
   
   if( (fp_matrix = fopen(mat_name_p, "r")) == NULL)
   {


   printf("\n\n Cannot find the file %s \n\n", mat_name );
   printf(" Make sure the environment variable DIALIGN2_DIR points\n");
   printf(" to a directory containing the files \n\n");
   printf("   BLOSUM \n   tp400_dna\n   tp400_prot \n   tp400_trans \n\n" );
   printf(" These files should be contained in the DIALIGN package \n\n\n" ) ;
   exit(1) ;




     printf("\n \n \n \n              ATTENTION ! \n \n");
     printf("\n   There is no similarity matrix `%s'. \n", mat_name);
     printf("   in the directory \n \n");
     printf("           %s\n \n", par_dir);
     exit(1);
   }
   }


    if( wgt_type != 1 )
      matrix_read( fp_matrix );

    mem_alloc(  );


    if( wgt_type != 1 )  
    if( (amino = (int **) calloc( seqnum , sizeof(int *) ) ) == NULL)
      {
         printf(" problems with memory allocation");
         printf(" for `amino' !  \n \n");
         exit(1);
      }

    if( wgt_type != 1 )
    for( i = 0 ; i < seqnum ; i++ ) 
    if( (amino[i] = (int *) calloc( ( seqlen[i]+5 ) , sizeof(int) ) ) == NULL)
      {
         printf(" problems with memory allocation");
         printf(" for `amino[%d]' !  \n \n", i);
         exit(1);
      }




 
    if( crick_strand ) { 
      if( (amino_c = (int **) calloc( seqnum , sizeof(int *) ) ) == NULL) {
        printf(" problems with memory allocation");
        printf(" for `amino_c' !  \n \n");
        exit(1);
      }

      for( i = 0 ; i < seqnum ; i++ )
        if( (amino_c[i] = (int *) calloc( ( seqlen[i]+5 ) , sizeof(int) ) ) == NULL) {
          printf(" problems with memory allocation");
          printf(" for `amino_c[%d]' !  \n \n", i);
          exit(1);
        }
    }
 

             /******************************************************  
             *                                                     *      
             *  read file, that contains data of anchored regions  *
             *                                                     *      
             ******************************************************/  



if( anchors ) {
  multi_anc_read( input_name );
}

if( exclude_frg ) { 

  if( ( exclude_list = (int ***) calloc( seqnum , sizeof(int **) )) == NULL) {
    printf(" problems with memory allocation for 'exclude_list' \n \n");
    exit(1);
  } 

  for(i = 0 ; i < seqnum ; i++ ) 
  if( ( exclude_list[ i ] = (int **) calloc( seqnum , sizeof(int *) )) == NULL) {
    printf(" problems with memory allocation for 'exclude_list' \n \n");
    exit(1);
  } 
  
  for(i = 0 ; i < seqnum ; i++ ) 
  for(j = 0 ; j < seqnum ; j++ ) 
  if( ( exclude_list[ i ][ j ]  = (int *) calloc( seqlen[ i ] + 1 , sizeof(int) )) == NULL) {
    printf(" problems with memory allocation for 'exclude_list' \n \n");
    exit(1);
  } 

  exclude_frg_read ( input_name , exclude_list ) ;
}



   if( wgt_type == 0 ) 
     tp400_read( 0 , tp400_prot);
   if( wgt_type % 2 )
     tp400_read( 1 , tp400_dna );
   if( wgt_type > 1 )
     tp400_read( 2 , tp400_trans );



           /****************************\
           *                            * 
           *    Name of output files    *  
           *                            * 
           \****************************/

   if( default_name )
     {
       strcpy( printname , input_name);
       strcpy( prn , printname);
     } 
   else
     {  
       strcpy( printname , output_name );
       strcpy( prn , printname);
     }
    

   strcpy(prn2 , prn); 
  
   if( default_name )
     strcat(prn,".ali");

   strcat(prn2,".fa");  
    


   strcpy(logname,printname);
   strcat(logname,".log");

   strcpy(fsm_name , printname);
   strcat(fsm_name,".fsm");

   if( print_status ) {
     strcpy( pst_name , printname );
     strcat( pst_name,".sta");
   }    

   if( afc_file )
     {
       strcpy( dia_name , printname );  
       strcat( dia_name , ".afc" );
       fp_dia = fopen( dia_name , "w" );
       fprintf(fp_dia,"\n #  %s \n\n  seq_len: " , input_line );
       for( i = 0 ; i < seqnum ; i++ )
         fprintf(fp_dia,"  %d ", seqlen[i] );
       fprintf(fp_dia,"\n\n");

     }

   if( col_score ) { 
     strcpy( csc_name , printname );  
     strcat( csc_name , ".csc" );
     fp_csc = fopen( csc_name , "w" );
   }

   if( dia_pa_file )
     {
       strcpy( dia_pa_name , printname );  
       strcat( dia_pa_name , ".fop" );

       fp_dpa = fopen( dia_pa_name , "w" );


       fprintf(fp_dpa,"\n #  %s \n\n  seq_len: " , input_line );
       for( i = 0 ; i < seqnum ; i++ ) 
         fprintf(fp_dpa,"  %d ", seqlen[i] ); 
       fprintf(fp_dpa,"\n\n");
       fclose( fp_dpa ) ;
     }


   if( motifs ) {
     strcpy( mot_file_name , printname );  
     strcat( mot_file_name , ".mot" );
     fp_mot = fopen( mot_file_name , "w" );
      
     fprintf(fp_mot,"\n #  %s \n\n   " , input_line );
     fprintf(fp_mot," motif: %s \n\n", mot_regex ); 
     fprintf(fp_mot," max offset for motifs = %d \n\n", (int) max_mot_offset ); 
     fprintf(fp_mot," the following fragments contain the motif: \n\n" ); 
     fprintf(fp_mot,"   seq1 seq2    beg1 beg1 len    wgt" ); 
     fprintf(fp_mot,"   # mot    mot_wgt  \n\n" ); 
   }


   if( frag_file ) {
     strcpy( frag_file_name , printname );  
     strcat( frag_file_name , ".frg" );
     fp_frg = fopen( frag_file_name , "w" );
      
     fprintf(fp_frg,"\n #  %s \n\n  seq_len: " , input_line );
     for( i = 0 ; i < seqnum ; i++ )
       fprintf(fp_frg,"  %d ", seqlen[i] );
     fprintf(fp_frg,"\n  sequences: " );
     for( i = 0 ; i < seqnum ; i++ )
       fprintf(fp_frg,"  %s ", seq_name[i] );

     fprintf(fp_frg ,"\n\n");
   }



  clos = newAligGraphClosure(seqnum, seqlen, 0, NULL);

  if( (open_pos = (int *** ) calloc( seqnum , sizeof(int **))) == NULL)
     {
       printf("Problems with memory allocation for open_pos\n"); 
       exit(1);
     }

  for(i=0;i<seqnum;i++)
      {
        if( (open_pos[i] = 
        (int ** ) calloc( seqnum , sizeof(int *))) == NULL)
          { 
             printf("Problems with memory allocation for open_pos\n"); 
             exit(1);
          }
      }


  for(i=0;i<seqnum;i++)
  for(j=0;j<seqnum;j++)
      {
       if( (open_pos[i][j] = 
       (int * ) calloc( ( seqlen[i]+2) , sizeof(int) ) ) == NULL)
          { 
             printf("Problems with memory allocation for open_pos\n"); 
             exit(1);
          }
      }

  for( i = 0 ; i <seqnum ; i++)
  for( j = 0 ; j <seqnum ; j++)
  for( hv = 1 ; hv <= seqlen[i] ; hv++)
     open_pos[i][j][hv] = 1;


   	  /**************************************
          *                                     *
          *      definition of  `amino'         *       
    	  *                                     *
          **************************************/




  if( wgt_type > 1 ) 
    for(hv=0;hv<seqnum;hv++)
    for(i=1;i<=seqlen[hv]-2;i++)
      {


        if( translate( seq[hv][i],seq[hv][i+1],seq[hv][i+2],hv,i ) == -1)
          exit(1);


        amino[hv][i] = translate( seq[hv][i],seq[hv][i+1],seq[hv][i+2],hv,i);
   
        if( crick_strand ) { 
          nuc1 = invert( seq[hv][i+2] );
          nuc2 = invert( seq[hv][i+1] );
          nuc3 = invert( seq[hv][i] );
 
          amino_c[hv][i] = translate( nuc1 , nuc2 , nuc3 , hv , i);
        }
      }


   if( wgt_type == 0 ) 
   for(hv=0;hv<seqnum;hv++)
   for(i=1;i<=seqlen[hv];i++)
    {
     if( seq[hv][i] == 'C' ) amino[hv][i] = 1;           
     if( seq[hv][i] == 'S' ) amino[hv][i] = 2;           
     if( seq[hv][i] == 'T' ) amino[hv][i] = 3;           
     if( seq[hv][i] == 'P' ) amino[hv][i] = 4;           
     if( seq[hv][i] == 'A' ) amino[hv][i] = 5;           
     if( seq[hv][i] == 'G' ) amino[hv][i] = 6;           
     if( seq[hv][i] == 'N' ) amino[hv][i] = 7;           
     if( seq[hv][i] == 'D' ) amino[hv][i] = 8;           
     if( seq[hv][i] == 'E' ) amino[hv][i] = 9;           
     if( seq[hv][i] == 'Q' ) amino[hv][i] = 10;           
     if( seq[hv][i] == 'H' ) amino[hv][i] = 11;           
     if( seq[hv][i] == 'R' ) amino[hv][i] = 12;           
     if( seq[hv][i] == 'K' ) amino[hv][i] = 13;           
     if( seq[hv][i] == 'M' ) amino[hv][i] = 14;           
     if( seq[hv][i] == 'I' ) amino[hv][i] = 15;           
     if( seq[hv][i] == 'L' ) amino[hv][i] = 16;           
     if( seq[hv][i] == 'V' ) amino[hv][i] = 17;           
     if( seq[hv][i] == 'F' ) amino[hv][i] = 18;           
     if( seq[hv][i] == 'Y' ) amino[hv][i] = 19;           
     if( seq[hv][i] == 'W' ) amino[hv][i] = 20;           
    }


     
     amino_acid[0] = 'X';           
     amino_acid[1] = 'C';           
     amino_acid[2] = 'S';           
     amino_acid[3] = 'T';           
     amino_acid[4] = 'P';           
     amino_acid[5] = 'A';           
     amino_acid[6] = 'G';           
     amino_acid[7] = 'N';           
     amino_acid[8] = 'D';           
     amino_acid[9] = 'E';           
     amino_acid[10] = 'Q';           
     amino_acid[11] = 'H';           
     amino_acid[12] = 'R';           
     amino_acid[13] = 'K';           
     amino_acid[14] = 'M';           
     amino_acid[15] = 'I';           
     amino_acid[16] = 'L';           
     amino_acid[17] = 'V';           
     amino_acid[18] = 'F';           
     amino_acid[19] = 'Y';           
     amino_acid[20] = 'W';



num_dia_anc = anc_num * (seqnum-1);




if ( anchors ) {

  if( time_stamps )
    beg_ts = clock() ;

  if ( nas == 0 ) 
    if( bubblesort ) 
      bubble_sort ( anc_num , anchor_frg ) ;  
    else
      frag_sort ( anc_num , anchor_frg , 0 ) ;


  if( time_stamps) {
    end_ts = clock() ;
    time_diff_srt = (float) ( end_ts - beg_ts ) / CLOCKS_PER_SEC ;
    if( time_stamps )
      printf (" for anc: time_diff_srt = %f \n", time_diff_srt );
  }


  filter( &anc_num , anchor_frg);
/*  exit(1) ; 
*/
}


if(long_output)
  {
    fp_log = fopen(logname,"w");
    fprintf(fp_log,"\n #  %s \n\n   " , input_line );
  }

if(frg_mult_file) {
  fp_fsm = fopen(fsm_name,"w");
  fprintf(fp_fsm,"\n #  %s \n\n" , input_line );
}





if( 
    (  num_dia_bf = (int *) calloc( ( max_itnum + 1 )  ,  sizeof( int ) )  )
      == NULL
  )
     {
          printf(" problems with memory allocation for `num_dia_bf' !  \n \n");
          exit(1);
     }


if( 
    (  num_dia_af = (int *) calloc( ( max_itnum + 1 )  ,  sizeof( int ) )  )
      == NULL
  )
     {
          printf(" problems with memory allocation for `num_dia_af' !  \n \n");
          exit(1);
     }


all_it_dia = (struct multi_frag *) calloc( 1 , sizeof(struct multi_frag) ); 
current_dia = all_it_dia;


  strcpy(itname,printname);
  strcpy(itname2,printname);
  strcpy(itname3,printname);
  strcpy(itname4,printname);
  sprintf(str,".ali");
  
  if( default_name )
  strcat(itname,str);
  
  sprintf(str,".fa");
  strcat(itname2,str);


  if( msf_file )
  strcat(itname3,".ms");
  
  
  if( cw_file )
  strcat(itname4,".cw");
  





        if( textual_alignment )
          fp_ali = fopen(itname,"w");
 
        if( standard_out )
          fp_ali = stdout;

   
        if( textual_alignment )
        if(fasta_file)
          fp2 = fopen(itname2,"w");

        if(msf_file)
          fp3 = fopen(itname3,"w");

        if(cw_file)
          fp4 = fopen(itname4,"w");

        if( textual_alignment )
          para_print(seq_file , fp_ali);


                    /***************************\
                    *                           * 
                    *      ITERATION START      *   
                    *                           * 
                    \***************************/









istep = 0 ; 
 while( ( cont_it == 1 ) && ( istep < max_itnum ) ) 
  {

    cont_it = 0 ;
    istep++ ; 

/* printf("\n  istep = %d \n", istep ); */


    this_it_dia = current_dia;
 
    strcpy(itname,printname);
    strcpy(itname2,printname);
    strcpy(itname3,printname);
    strcpy(itname4,printname);
    sprintf(str,".ali");

        
    if( default_name )
      strcat(itname,str); 

    sprintf(str,".fa");
    strcat(itname2,str); 

    if( msf_file )
      strcat(itname3,".ms"); 
   
    
    if( cw_file )
      strcat(itname4,".cw"); 
  
    weight_sum_af = 0;
    num_dia_bf[ istep ] = 0;
 
    if( time_stamps ) 
      beg_pa = clock(); 


    if( ref_seq == 0 ) 
      i_max = seqnum ; 
    else
      i_max = 1 ; 

    for(i = 0 ; i < i_max ; i++)
      {

        for(j = i + 1 ; j < seqnum ; j++)
          {


              /****************************************\
              *                                        * 
              *          PAIRWISE  ALIGNMENT           *
              *                                        * 
              \****************************************/

                if( cont_it_p[ i ][ j] ) {

/*
                  printf("\n out of frc it %d : wgt 20 = %f \n", istep ,  wgt_dna[ 20 ][ 20 ] ) ;
*/
                  score = frag_chain( i , j , fp_ali, fp_mot, &num_dia_p );
                }
                else {
                  score = 0 ; 
                  num_dia_p = 0 ; 
                }

                if( istep == 1 )
                  {
                    pair_score[j][i] = score;  
                    pair_score[i][j] = score;  
                  }


                for(k=0;k<num_dia_p;k++)
                  {
                    *current_dia = pair_dia[k];

                    current_dia->next 
                    = (struct multi_frag *) calloc( 1 , sizeof(struct multi_frag) );
                    end_dia = current_dia;   
                    current_dia = current_dia->next;
                    current_dia->pred = end_dia; 
                  }

                num_dia_bf[ istep ] = num_dia_bf[ istep ] + num_dia_p;

                for(hv=0; hv<num_dia_p;hv++)
                  weight_sum_af = weight_sum_af + (pair_dia[hv]).weight;


            if(num_dia_p)
              free(pair_dia);

          }    /*    for(j = i+1 ; j<seqnum ; j++) */

      }        /*   for(i = 0 ; i<seqnum ; i++) */


    if( time_stamps ) { 
      end_pa = clock(); 

      time_diff_pa = (float) ( end_pa - beg_pa ) / CLOCKS_PER_SEC ; 
      if( time_stamps )
        printf (" time_diff_pa = %f \n", time_diff_pa ); 
      total_pa_time = total_pa_time + time_diff_pa;
    }




    if( break1 ) {
      printf("\n  break1\n");
      exit(1) ;
    }


/*
    if( pa_only ) {
      printf("\n\n istep = %d, pa finished - exit \n\n", istep );       
      exit(1);
    }
*/

    if(overlap_weights)
      {
        diagonal1 = this_it_dia;
        dia_counter = 0;

        if( diagonal1 != NULL )   
        while( diagonal1->next != NULL )   
          {
            dia_counter++;
            if( print_status )
            if( ( dia_counter % 100 ) == 0 )
              {                
                fp_st = fopen( pst_name ,"w");

                fprintf(fp_st," dsd  %s \n", input_line);
                fprintf(fp_st,"\n\n\n    Status of the program run:\n");  
                fprintf(fp_st,"    ==========================\n\n");  
                if( seqnum > 2 ) { 
                  fprintf(fp_st,"      iteration step %d in ", istep); 
                  fprintf(fp_st,"multiple alignment\n" );
                }
                fprintf(fp_st,"      calculating overlap weight for diagonals\n");
                fprintf(fp_st,"      current diagonal = %d\n\n", dia_counter );
                fprintf(fp_st,"      total number of"); 
                fprintf(fp_st," diagonals: %d\n\n\n\n", num_dia_bf[ istep ]);
                fclose(fp_st);
              }

            diagonal2 = diagonal1->next;

            while(diagonal2->next != NULL) 
              {
                if( diagonal1->trans == diagonal2->trans ) 
                  ow_add(diagonal1 , diagonal2); 
                diagonal2 = diagonal2->next;        
              }
            diagonal1 = diagonal1->next;        
          }
        if( bubblesort )   
          ow_bubble_sort( num_dia_bf[ istep ] , this_it_dia ); 
        else 
          frag_sort( num_dia_bf[ istep ] , this_it_dia , overlap_weights ); 
      }
    else /* no overlap_weights */ {
      beg_ts = clock() ; 

      if( bubblesort ) 
        bubble_sort( num_dia_bf[ istep ] , this_it_dia );
      else 
        frag_sort( num_dia_bf[ istep ] , this_it_dia , overlap_weights );

      end_ts = clock() ; 
      time_diff_srt = (float) ( end_ts - beg_ts ) / CLOCKS_PER_SEC ;
      if( time_stamps )
        printf (" time_diff_srt = %f \n", time_diff_srt );
  }


    num_dia_af[ istep ] = num_dia_bf[ istep ];
    weight_sum_bf = weight_sum_af;

    pairalignsum = 0;
    pairalignlen = 0;


    filter( num_dia_af + istep , this_it_dia); 
    num_all_it_dia = num_all_it_dia + num_dia_af[ istep ];


/*
    if( pa_only == 0 ) {
      printf("\n\n istep = %d, filter finished - exit \n\n", istep );       
      exit(1);
    }
*/


    weight_sum_af = 0;
       
    print_log( this_it_dia , fp_log , fp_fsm );

    if( frag_file )
      print_fragments( this_it_dia , fp_frg );

    throw_out( &weight_sum_af );

    sel_test( );

     

    threshold = threshold ;

    if( break2 ) {
      printf("\n  break2\n");
      exit(1) ;
    }


  } /* while ( cond_it == 1 ) */  


                    /***************************\
                    *                           * 
                    *       ITERATION END       *   
                    *                           * 
                    \***************************/

strcpy( dist_name , printname);
strcat(dist_name , ".dst");





if ( ref_seq == 0 ) 
  av_tree_print();



      if( standard_out )
        fp_ali = stdout;
      
if(sf_mat){
  subst_mat( input_name , num_all_it_dia , all_it_dia ) ;
}


if( textual_alignment )
  ali_arrange( num_all_it_dia , all_it_dia , fp_ali , fp2, fp3, fp4, fp_csc );


if(long_output) 
  {
/*     fprintf(fp_log "\n\n thr = %f , lmax = %d , speed = %f  */  
    fprintf(fp_log, "\n\n    total sum of weights: %f \n\n\n", tot_weight);
    fclose(fp_log);
  }




if( argnum == 1 )
  {
    printf("\n     Program terminated normally\n");
    printf("     Results are contained in file `%s' \n \n \n", itname);
  }


av_dia_num = 2 * dia_num ;
av_dia_num = av_dia_num / ( seqnum * ( seqnum - 1) ) ;

av_max_dia_num = 2 * max_dia_num ;
av_max_dia_num = av_max_dia_num / ( seqnum * ( seqnum - 1) ) ;



tmpi1 = av_dia_num ;
tmpi2 = av_max_dia_num ;

if(pr_av_nd)
  printf(" %d ", tmpi1 );

if(pr_av_max_nd)
  printf(" %d ", tmpi2 );



if(pr_av_nd)
  fprintf(fp_ali, "    %d fragments considered for alignment \n", tmpi1 );

if(pr_av_max_nd)
  fprintf(fp_ali, "    %d fragments simultaneously stored \n\n", tmpi2 );

if( textual_alignment )
  fclose(fp_ali);


  if( time_stamps ){ 
    end_ali = clock() ; 
    time_diff_ali = (float) ( end_ali - beg_ali ) / CLOCKS_PER_SEC ;

    perc_pa_time = total_pa_time / time_diff_ali * 100 ; 
    printf (" time_diff_ali = %f \n", time_diff_ali ); 
    printf (" total_pa_time = %f \n", total_pa_time );
    printf (" corresponds to %f percent \n\n", perc_pa_time );
  } 
}    /* main */



