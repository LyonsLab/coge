
#define PAPER_WIDTH      80 
#define MLINE 1000 
#define MAX_REGEX 1000 
#define NAME_LEN 1000   
#define SEQ_NAME_LEN 12   
#define MAX_SEQNUM 10000
#define MAX_ITNUM 3 
#define MAX_INPUT_LINE 10000 
#define MIN_MOT_WGT 0.1 
#define MAX_CSC 10 


         /**************************\
         *                          * 
         *    default parameters    *
         *                          * 
         \**************************/


#define BETA              0
#define WEB               0
#define OVERLAP_THRESHOLD 35
#define MIN_DIA           1
#define MAX_DIA          40
#define MATNAME          "BLOSUM"     
#define WEAK_WGT_TYPE_THR      0.5 
#define STRONG_WGT_TYPE_THR    0.75 


struct pair_frag {int b1, b2, ext; float weight; short trans, cs; 
                   struct pair_frag *prec, *last; float sum; };
     /* 
           fragments in function `pairalign' 

           b1, b2:    begin of the diagonal
           ext:       length of the diagonal
           weight:    weight of the diagonal
           prec:      preceding diagonal in dot matrix
           last:      last diagonal ending in the same column 
           sum:       sum of weights accumulated  
           cs:        crick strand 
           trans:     translation
     */ 

struct multi_frag {int b[2], s[2], ext, it; float weight, ow; short sel, trans;
                   short cs; struct multi_frag *next, *pred;};
     /*
           fragments outside function `pairalign' 

           b[0], b[1]:  begin of the diagonal
           s[0], s[1]:  sequences, to which diagonal belongs
           ext:         length of the diagonal
           weight:      individual weight of the diagonal
           ow:          overlap weight of the diagonal
           sel:         1, if accepted in filter proces, 0 else
           trans:       translation
           cs:          crick strand 
           it:          iteration step 
           *next:       next diagonal 
     */

struct leaf {int s1, s2, clade;};
struct seq_pair {int s1, s2; float weight;};      

struct subtree { int member_num, valid ; int *member; char *name ;
                 float depth; };         


