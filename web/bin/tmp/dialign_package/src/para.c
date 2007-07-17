      
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "define.h"
#include "dialign.h"

int    lgs_option = 0 ; 
float  sf_mat_thr = 0 ;

extern int col_score , regex_len , wgt_type_plot, cd_gobics, motifs ; 
extern int bubblesort , nas , ref_seq, speed_optimized , online, time_stamps ; 
extern int thr_sim_score, lmax, redundant, seqnum, strict, sf_mat;
extern int quali_num , anchors , mask, textual_alignment ;
extern int pa_only , dna_speed , frg_mult_file ,  frg_mult_file_v ; 
extern int pr_av_nd, pr_av_max_nd, wgt_print , wgt_print_x ;
extern int standard_out, self_comparison;
extern int default_name;
extern char mot_regex[ MAX_REGEX ] ,   output_name[ NAME_LEN ], clust_sim[ NAME_LEN ]; 
extern int afc_file, afc_filex , frag_file; 
extern int dia_pa_file, iter_cond_prob , fasta_file, cw_file;
extern int overlap_weights, ow_force, long_output ;
extern int msf_file, wgt_type,  print_status ;
extern int  plot_num , break1, break2 ; 
extern float threshold , mot_factor , mot_offset_factor ; 
extern int num_test( char *cp );
extern short crick_strand, exclude_frg , max_itnum ;

extern void regex_format_complain() ; 


void para_read( int num , char ** arg )      
  {

    int an = 1;

    
    while( an < num - 1 )
      { 

        if( 
            strcmp( arg[an] , "-afc")  &&      /* create file containing 
                                                  ALL fragments considered 
                                                  for alignment */
            strcmp( arg[an] , "-afc_v") &&      /* like -afc with fragments
                                                  explicitly printed */
            strcmp( arg[an] , "-b1")   &&      /* break */ 
            strcmp( arg[an] , "-b2")   &&      /* break */ 
            strcmp( arg[an] , "-bs")   &&      /* bubble sort */ 
            strcmp( arg[an] , "-csc") && /* column score output */ 
            strcmp( arg[an] , "-cs")   &&      /* crick strand */ 
            strcmp( arg[an] , "-cw")   &&      /* additional output file
                                                  in clustalw format */
            strcmp( arg[an] , "-d1w")  &&      /* old weight fkt */

            strcmp( arg[an] , "-ds")   &&    
            strcmp( arg[an] , "-fa")   &&      /* separate file with 
                                                   alignment in fasta format */
            strcmp( arg[an] , "-ff")   &&      /* fragment file  */
            strcmp( arg[an] , "-fn")   &&      /* name of output file  */
            strcmp( arg[an] , "-fop")  &&      /* create file containing 
                                                  fragments selected for 
                                                  optimal pairwise alignment */
            strcmp( arg[an] , "-fsm")  &&      /* create file containing 
                                                  consistent fragments in 
                                                  multiple alignment (in 
                                                  format needed for -xfr ) */
            strcmp( arg[an] , "-fsmv") &&      /* same as -fsm but verbose */ 
            strcmp( arg[an] , "-cd_gobics") && /* chaos + dialign @ gobics */
            strcmp( arg[an] , "-lgs_t")   &&      /* genomic sequences, transl. */
            strcmp( arg[an] , "-istep")   &&   /* max iteration steps */ 
            strcmp( arg[an] , "-it")   &&      /* iteration */ 
            strcmp( arg[an] , "-iw")   &&      /* ind. weights */ 
            strcmp( arg[an] , "-lgs")  &&      /* genomic sequences  */
            strcmp( arg[an] , "-lgsx")  &&      /* genomic sequences, accurate + textual alignment   */
            strcmp( arg[an] , "-lmax") &&      /* max. length of diag. */  
            strcmp( arg[an] , "-lo")   &&      /* long output */ 
            strcmp( arg[an] , "-ma")   &&      /* mixed weights */
            strcmp( arg[an] , "-anc")  &&      /* anchor regions */ 
            strcmp( arg[an] , "-mask") && 
            strcmp( arg[an] , "-mat")  &&      /* calc. subst. freq. matrix */
            strcmp( arg[an] , "-mat_thr")  &&      /* thr for sbst. fr. mat. */
            strcmp( arg[an] , "-max_link") &&  /* max. linkage clustering */ 
            strcmp( arg[an] , "-min_link") &&  /* min. linkage clustering */ 
            strcmp( arg[an] , "-mot")  &&      /* motifs considered */ 
            strcmp( arg[an] , "-msf")  &&      /* separate file with 
                                                   alignment in msf format */
            strcmp( arg[an] , "-n")    &&      /* DNA/RNA sequences */
            strcmp( arg[an] , "-nas")    &&      /* no anchor sorting */
            strcmp( arg[an] , "-nt")   &&      /* DNA/RNA sequences with 
                                                   translation option */
            strcmp( arg[an] , "-nta")  &&      /* no textual alignment */ 
            strcmp( arg[an] , "-o")    &&      /* optimized  */
            strcmp( arg[an] , "-online")    && /* online */
            strcmp( arg[an] , "-ow")   &&      /* overlap weights */ 
            strcmp( arg[an] , "-pamnd") &&  /* print av. max. number of frg. */
            strcmp( arg[an] , "-pand") &&      /* print av. number of diag. */
            strcmp( arg[an] , "-pao")  &&      /* pairw. alignments only */
            strcmp( arg[an] , "-ref_seq")  &&  /* seq_2, ... , seq_n 
                                                 aligned to seq_1 */
            strcmp( arg[an] , "-stars")&&      /* maximum number of stars under 
                                                alignment indicating relative similarity*/
            strcmp( arg[an] , "-pst")  &&      /* print status */
            strcmp( arg[an] , "-sc")  &&      /* self comparison */
            strcmp( arg[an] , "-smin")  &&  
            strcmp( arg[an] , "-stdo")  &&      /* standard output */
	    strcmp( arg[an] , "-ta")  &&       /* textual alignment*/
	    strcmp( arg[an] , "-thr") &&        /* threshold */
	    strcmp( arg[an] , "-ts") &&        /* time stamps */
	    strcmp( arg[an] , "-wgtpr") &&        /* weight print */
	    strcmp( arg[an] , "-wgtprx") &&        /* weight print */
	    strcmp( arg[an] , "-wtp") &&        /* weight type plot */
	    strcmp( arg[an] , "-xfr")         /* excluded fragments */

          )
          {
            printf("\n \n   Arguments in command line make no sense! \n \n");
            printf("\n   Unknown option %s \n \n \n \n",  arg[an] );
            exit(1);
          }  

        if( !strcmp( arg[an] , "-afc") )
          afc_file = 1;

        if( !strcmp( arg[an] , "-afc_v") ) { 
          afc_file = 1;
          afc_filex = 1 ;
        }

        if( !strcmp( arg[an] , "-b1") )
          break1 = 1;

        if( !strcmp( arg[an] , "-b2") )
          break2 = 1;

        if( !strcmp( arg[an] , "-bs") )
          bubblesort = 1;
        
        if( !strcmp( arg[an] , "-csc") )
          col_score = 1;

        if( !strcmp( arg[an] , "-cd_gobics") )
          cd_gobics = 1;

        if( !strcmp( arg[an] , "-cs") )
          crick_strand = 1;

        if( !strcmp( arg[an] , "-cw") )
          cw_file = 1;
	
        if( !strcmp( arg[an] , "-ds") )
          dna_speed = 1 ;

        if( !strcmp( arg[an] , "-fa") )
          fasta_file = 1;

        if( !strcmp( arg[an] , "-ff") )
          frag_file = 1;

        if( !strcmp( arg[an] , "-fop") )
          dia_pa_file = 1;

        if( !strcmp( arg[an] , "-fsm") )
          frg_mult_file = 1;

        if( !strcmp( arg[an] , "-fsmv") ) { 
          frg_mult_file = 1;
          frg_mult_file_v = 1;
        }

        if( !strcmp( arg[an] , "-it") )
          iter_cond_prob = 1;

        if( !strcmp( arg[an] , "-iw") )
          overlap_weights = 0;

        if( !strcmp( arg[an] , "-lgs") ) {
          wgt_type = 3 ;
/*          iter_cond_prob = 1 ; 
*/ 
          threshold = 2.0 ;
          lmax = 30 ;
          thr_sim_score = 8 ;
          strict = 1 ;
          textual_alignment = 0 ;
          /* dia_pa_file = 1; */ 
          frag_file = 1 ;
          dna_speed = 1 ;
          crick_strand = 1 ;
          lgs_option = 1 ;
          print_status = 1 ; 
        }

        if( !strcmp( arg[an] , "-lgs_t") ) {
          wgt_type = 2 ;
          iter_cond_prob = 1 ;
          threshold = 0.0 ;
          lmax = 30 ;
          thr_sim_score = 8 ;
          strict = 1 ;
          textual_alignment = 0 ;
          dia_pa_file = 1;
          frag_file = 1 ;
          dna_speed = 1 ;
          print_status = 1 ; 
        }

        if( !strcmp( arg[an] , "-lgsx") ) {
          wgt_type = 3 ;
          iter_cond_prob = 1 ;
          strict = 1 ;
          frag_file = 1 ;
          crick_strand = 1 ;
          lgs_option = 1 ;
          print_status = 1 ; 
        }

        if( !strcmp( arg[an] , "-lo") )
          long_output = 1;

        if( !strcmp( arg[an] , "-ma") ) {
          wgt_type = 3;
        }

        if( !strcmp( arg[an] , "-anc") )
          anchors = 1;

        if( !strcmp( arg[an] , "-mask") )
          mask = 1;

        if( !strcmp( arg[an] , "-max_link") )
          strcpy (clust_sim , "max" );

        if( !strcmp( arg[an] , "-min_link") )
          strcpy (clust_sim , "min" );

        if( !strcmp( arg[an] , "-msf") )
          msf_file = 1;

        if( !strcmp( arg[an] , "-n") ) {
          wgt_type = 1;
        }

        if( !strcmp( arg[an] , "-nas") ) {
          nas = 1;
        }

        if( !strcmp( arg[an] , "-nt") )
          wgt_type = 2;

        if( !strcmp( arg[an] , "-nta") )
          textual_alignment = 0;

        if( !strcmp( arg[an] , "-o") )
          {
            speed_optimized = 1 ; 
            threshold = 0.5 ;
            lmax = 30 ;
            thr_sim_score = 8 ;
          }

        if( !strcmp( arg[an] , "-ow") )
          ow_force = 1;

        if( !strcmp( arg[an] , "-pao") )
          pa_only = 1;

        if( !strcmp( arg[an] , "-pamnd") )
          pr_av_max_nd = 1;

        if( !strcmp( arg[an] , "-pand") )
          pr_av_nd = 1;

        if( !strcmp( arg[an] , "-pst") )
          print_status = 1;

        if( !strcmp( arg[an] , "-red") )
          redundant = 1;

        if( !strcmp( arg[an] , "-mat") )
          sf_mat = 1;

        if( !strcmp( arg[an] , "-online") )
          online = 1;
 
        if( !strcmp( arg[an] , "-ref_seq") )
          ref_seq = 1;
 
        if( !strcmp( arg[an] , "-sc") )
          self_comparison = 1;

        if( !strcmp( arg[an] , "-stdo") )
          standard_out = 1;

        if( !strcmp( arg[an] , "-strict") )
          strict = 1;
 
        if( !strcmp( arg[an] , "-ta") )
          textual_alignment = 1 ;

        if( !strcmp( arg[an] , "-ts") )
          time_stamps = 1 ;

        if( !strcmp( arg[an] , "-wgtpr") )
          wgt_print = 1 ;

        if( !strcmp( arg[an] , "-wgtprx") )  
          wgt_print_x = 1 ;

        if( !strcmp( arg[an] , "-wtp") )  
          wgt_type_plot = 1 ;

        if( !strcmp( arg[an] , "-xfr") )
          exclude_frg = 1 ;


 
 
	/********************************************************************/


        if( !strcmp( arg[an] , "-fn") )
        if( an + 2 < num )  
          { 
            strcpy( output_name , arg[++an] );
            default_name = 0;
          } 
        else
          {
            printf("\n \n   Arguments in command line don't make sense! \n");
            printf("   (Name of output file not properly specified) \n \n");
            exit(1);
          }  



	/********************************************************************/


        if( !strcmp( arg[an] , "-istep") )
        if( ( an + 2 < num ) && num_test( arg[an + 1] ) )  
          max_itnum = atoi( arg[++an] ); 
        else
          {
            printf("\n \n   Arguments in command line don't make sense! \n");
            printf("   (max_itnum not properly specified) \n \n");
            exit(1);
          }  


	/********************************************************************/


        if( !strcmp( arg[an] , "-lmax") )
        if( ( an + 2 < num ) && num_test( arg[an + 1] ) )  
          lmax = atoi( arg[++an] ); 
        else
          {
            printf("\n \n   Arguments in command line don't make sense! \n");
            printf("   (lmax not properly specified) \n \n");
            exit(1);
          }  


        /********************************************************************/


        if( !strcmp( arg[an] , "-stars") )
        if( ( an + 2 < num ) && num_test( arg[an + 1] ) ) { 
          plot_num = atoi( arg[++an] ); 
          quali_num = 0 ; 
        }
        else
          {
            printf("\n \n   Arguments in command line don't make sense! \n");
            printf("   (Number of \"*\" characters not properly specified) \n \n");
            exit(1);
          }  


	/********************************************************************/


        if( !strcmp( arg[an] , "-smin") )
        if( (an + 2 < num) && num_test( arg[an + 1] ) )  
          thr_sim_score = atoi( arg[++an] );
        else
          {
            printf("\n \n   Arguments in command line don't make sense! \n");
            printf("   (Speed not properly specified) \n \n");
            exit(1);
          }  


        /********************************************************************/


        if( !strcmp( arg[an] , "-thr") )
        if( (an + 2 < num) && num_test( arg[an + 1] ) )  
          {
            threshold = atof( arg[++an] );
          }
        else
          {
            printf("\n \n   Arguments in command line don't make sense! \n");
            printf("   (Threshod not properly specified) \n \n");
            exit(1);
          }  


	/********************************************************************/


        if( !strcmp( arg[an] , "-mat_thr") )
        if( (an + 2 < num) && num_test( arg[an + 1] ) )  
          {
            sf_mat_thr = atof( arg[++an] );
          }
        else
          {
            printf("\n \n   Arguments in command line don't make sense! \n");
            printf("   (subst. mat. threshod not properly specified) \n \n");
            exit(1);
          }  


	/********************************************************************/


        if( !strcmp( arg[an] , "-mot") )
        if(  ( an + 4 < num )            && 
             num_test( arg[ an + 2 ] )   && 
             num_test( arg[ an + 3 ] ) 
          )  
          {
            motifs = 1 ;  
            strcpy( mot_regex , arg[++an] );
            mot_factor = atof( arg[++an] ) ;
            mot_offset_factor = atof( arg[++an] ) ;
            regex_len = strlen( mot_regex ) ;
          }
        else
          regex_format_complain();

        /********************************************************************/

        an++;
       } 
  }



  
