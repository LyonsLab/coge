
/* version 3.6. (c) Copyright 1993-2004 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#include <float.h>

#include "phylip.h"
#include "seq.h"

#define initialv        0.1     /* starting value of branch length          */

#define over            60   /* maximum width of a tree on screen */

/* Define this to print messages when free_trans() is misused.
 * See FIXME note in free_trans() */
/* #define TRANS_DEBUG */

#ifndef OLDC
/* function prototypes */
void   restml_inputnumbers(void);
void   getoptions(void);
void   allocrest(void);
void   setuppie(void);
void   doinit(void);
void   inputoptions(void);
void   restml_inputdata(void);
void   restml_sitesort(void);
void   restml_sitecombine(void);
void   makeweights(void);

void   restml_makevalues(void);
void   getinput(void);
void   copymatrix(transmatrix, transmatrix);
void   maketrans(double, boolean);
void   branchtrans(long, double);
double evaluate(tree *, node *);
boolean nuview(node *);
void   makenewv(node *);
void   update(node *);
void   smooth(node *);

void   insert_(node *p, node *);
void   restml_re_move(node **, node **);
void   restml_copynode(node *, node *);
void   restml_copy_(tree *, tree *);
void   buildnewtip(long , tree *);
void   buildsimpletree(tree *);
void   addtraverse(node *, node *, boolean);
void   rearrange(node *, node *);
void   restml_coordinates(node *, double, long *,double *, double *);
void   restml_fprintree(FILE *fp);
void   restml_printree(void);

double sigma(node *, double *);
void   fdescribe(FILE *, node *);
void   summarize(void);
void   restml_treeout(node *);
static phenotype2 restml_pheno_new(long endsite, long sitelength);
static void restml_pheno_delete(phenotype2 x2);
void initrestmlnode(node **p, node **grbg, node *q, long len, long nodei,
    long *ntips, long *parens, initops whichinit, pointarray treenode,
    pointarray nodep, Char *str, Char *ch, FILE *intree);
static void restml_unroot(node* root, node** nodep, long nonodes);
void   inittravtree(tree* t,node *);
static void adjust_lengths_r(node *p);
void   treevaluate(void);
void   maketree(void);
void   globrearrange(void); 
void   adjust_lengths(tree *);
double adjusted_v(double v);
sitelike2 init_sitelike(long sitelength);
void   free_sitelike(sitelike2 sl);
void   copy_sitelike(sitelike2 dest, sitelike2 src,long sitelength);
void   reallocsites(void);
static void set_branchnum(node *p, long branchnum);
void   alloctrans(tree *t, long nonodes, long sitelength);
long   get_trans(tree* t);
void   free_trans(tree* t, long trans);
void   free_all_trans(tree* t);
void   alloclrsaves(void);
void   freelrsaves(void);
void   resetlrsaves(void);
void   cleanup(void);
/* function prototypes */
#endif

Char   infilename[FNMLNGTH];
Char   outfilename[FNMLNGTH];
Char   intreename[FNMLNGTH];
Char   outtreename[FNMLNGTH];

long nonodes2, sites, enzymes, weightsum, sitelength, datasets, 
        ith, njumble, jumb=0;
long inseed, inseed0;

/* User options */

boolean  global;        /* Perform global rearrangements? */
boolean  jumble;        /* Randomize input order? */
boolean  lengths;       /* Use lengths from user tree? */
boolean  weights;
boolean  trout;         /* Write tree to outtree? */
boolean  trunc8;
boolean  usertree;      /* Input user tree? Or search. */
boolean  progress;      /* Display progress */
boolean  mulsets;       /* Use multiple data sets */

/* Runtime state */

boolean  firstset;
boolean  improve;
boolean  smoothit;
boolean  inserting = false;

double bestyet;
tree curtree, priortree, bestree, bestree2;
longer seed;
long *enterorder;
steptr aliasweight;
char *progname;
node *qwhere,*addwhere;
/* local rearrangements need to save views.  
   created globally so that reallocation of the same variable is unnecessary */
node **lrsaves;

/* Local variables for maketree, propagated globally for C version: */
long       nextsp, numtrees, maxwhich, col, shimotrees;
double      maxlogl;
boolean     succeeded, smoothed;
#define NTEMPMATS 7
transmatrix *tempmatrix, tempslope, tempcurve;
sitelike2    pie;
double      *l0gl;
double     **l0gf;
Char ch;

/* variables added to keep treeread2() happy */
boolean goteof;
double trweight;
node *grbg = NULL;

#ifdef DEBUG
void checknode(node *p) {
  assert(p != NULL);
  assert(p->back != NULL);
  assert(p->back->back == p);
  //assert(p->v == p->back->v);
  //assert(p->branchnum > 0);
  //assert(p->branchnum == p->back->branchnum);
}

void checktree_r(node *p) {
  node *q;

  checknode(p);
  if (!p->tip) {
    for (q = p->next; q != p; q = q->next) {
      checktree_r(q->back);
    }
  }
}

void checktree() {
  checktree_r(curtree.start);
  checktree_r(curtree.start->back);
}
#endif /* DEBUG */


static void
set_branchnum(node *p, long branchnum)
{
  assert(p != NULL);
  assert(branchnum > 0);
  p->branchnum = branchnum;
}
  

void alloctrans(tree *t, long nonodes, long sitelength)
{
  /* used by restml */
  long i, j;

  t->trans = (transptr)Malloc(nonodes*sizeof(transmatrix));
  for (i = 0; i < nonodes; ++i){
    t->trans[i] = (transmatrix)Malloc((sitelength + 1) * sizeof(double *));
    for (j = 0;j < sitelength + 1; ++j)
      t->trans[i][j] = (double *)Malloc((sitelength + 1) * sizeof(double));
  }
 
  t->freetrans = Malloc(nonodes* sizeof(long));
  for ( i = 0; i < nonodes; i++ )
    t->freetrans[i] = i+1;
  t->transindex = nonodes - 1;
}  /* alloctrans */


long get_trans(tree* t)
{
  long ret;
  assert(t->transindex >= 0);
  ret = t->freetrans[t->transindex];
  t->transindex--;
  return ret;
}


void free_trans(tree* t, long trans)
{
  long i;

  /* FIXME This is a temporary workaround and probably slows things down a bit.
   * During rearrangements, this function is sometimes called more than once on
   * already freed nodes, causing the freetrans array to overrun other data. */
  for ( i = 0 ; i < t->transindex; i++ ) {
    if ( t->freetrans[i] == trans ) {
#ifdef TRANS_DEBUG
      printf("ERROR: trans %ld has already been freed!!\n", trans);
#endif
      return;
    }
  }
  /* end of temporary fix */
  
  t->transindex++;
  t->freetrans[t->transindex] = trans;
}


void free_all_trans(tree* t) 
{
  long i;

  for ( i = 0; i < nonodes2; i++ )
    t->freetrans[i] = i;
  t->transindex = nonodes2 - 1;
}


sitelike2 init_sitelike(long sitelength) 
{
  return Malloc((sitelength+1) * sizeof(double));
}


void free_sitelike(sitelike2 sl) 
{
  free(sl);
}


void copy_sitelike(sitelike2 dest, sitelike2 src,long sitelength)
{
  memcpy(dest,src,(sitelength+1)*sizeof(double));
}


void restml_inputnumbers()
{
  /* read and print out numbers of species and sites */
  fscanf(infile, "%ld%ld%ld", &spp, &sites, &enzymes);
  nonodes2 = spp * 2 - 1;
}  /* restml_inputnumbers */


void getoptions()
{
  /* interactively set options */
  long loopcount, loopcount2;
  Char ch;

  fprintf(outfile, "\nRestriction site Maximum Likelihood");
  fprintf(outfile, " method, version %s\n\n",VERSION);
  putchar('\n');
  sitelength = 6;
  trunc8 = true;
  global = false;
  improve = false;
  jumble = false;
  njumble = 1;
  lengths = false;
  outgrno = 1;
  outgropt = false;
  trout = true;
  usertree = false;
  weights = false;
  printdata = false;
  progress = true;
  treeprint = true;
  interleaved = true;
  loopcount = 0;
  for (;;) {
    cleerhome();
    printf("\nRestriction site Maximum Likelihood");
    printf(" method, version %s\n\n",VERSION);
    printf("Settings for this run:\n");
    printf("  U                 Search for best tree?  %s\n",
           (usertree ? "No, use user trees in input file" : "Yes"));
    if (usertree) {
      printf("  N          Use lengths from user trees?  %s\n",
              (lengths ? "Yes" : "No"));
    }
    printf("  A               Are all sites detected?  %s\n",
           (trunc8 ? "No" : "Yes"));
    if (!usertree) {
      printf("  S        Speedier but rougher analysis?  %s\n",
             (improve ? "No, not rough" : "Yes"));
      printf("  G                Global rearrangements?  %s\n",
             (global ? "Yes" : "No"));
      printf("  J   Randomize input order of sequences?  ");
      if (jumble)
        printf("Yes (seed =%8ld,%3ld times)\n", inseed0, njumble);
      else
        printf("No. Use input order\n");
    }
    printf("  L                          Site length?%3ld\n",sitelength);
    printf("  O                        Outgroup root?  %s%3ld\n",
           (outgropt ? "Yes, at sequence number" :
                       "No, use as outgroup species"),outgrno);

    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld sets\n", datasets);
    else
      printf("  No\n");
    printf("  I          Input sequences interleaved?  %s\n",
           (interleaved ? "Yes" : "No, sequential"));
    printf("  0   Terminal type (IBM PC, ANSI, none)?  %s\n",
           ibmpc ? "IBM PC" : ansi  ? "ANSI" : "(none)");
    printf("  1    Print out the data at start of run  %s\n",
           (printdata ? "Yes" : "No"));
    printf("  2  Print indications of progress of run  %s\n",
           (progress ? "Yes" : "No"));
    printf("  3                        Print out tree  %s\n",
           (treeprint ? "Yes" : "No"));
    printf("  4       Write out trees onto tree file?  %s\n",
           (trout ? "Yes" : "No"));
    printf("\n  Y to accept these or type the letter for one to change\n");
    fflush(stdout);
    scanf("%c%*[^\n]", &ch);
    getchar();
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
    if (ch == 'Y')
      break;
    if (((!usertree) && (strchr("UNASGJOTMI01234", ch) != NULL))
          || (usertree && ((strchr("UNASLOTMI01234", ch) != NULL)))){
      switch (ch) {

      case 'A':
        trunc8 = !trunc8;
        break;
        
      case 'S':
        improve = !improve;
        break;

      case 'G':
        global = !global;
        break;
        
      case 'J':
        jumble = !jumble;
        if (jumble)
          initjumble(&inseed, &inseed0, seed, &njumble);
        else njumble = 1;
        break;
        
      case 'L':
        loopcount2 = 0;
        do {
          printf("New Sitelength?\n");
          fflush(stdout);
          scanf("%ld%*[^\n]", &sitelength);
          getchar();
          if (sitelength < 1)
            printf("BAD RESTRICTION SITE LENGTH: %ld\n", sitelength);
          countup(&loopcount2, 10);
        } while (sitelength < 1);
        break;
        
      case 'N':
        lengths = !lengths;
        break;

      case 'O':
        outgropt = !outgropt;
        if (outgropt)
          initoutgroup(&outgrno, spp);
        else outgrno = 1;
        break;
        
      case 'U':
        usertree = !usertree;
        break;
        
      case 'M':
        mulsets = !mulsets;
        if (mulsets)
          initdatasets(&datasets);
        break;
        
      case 'I':
        interleaved = !interleaved;
        break;
        
      case '0':
        initterminal(&ibmpc, &ansi);
        break;
        
      case '1':
        printdata = !printdata;
        break;
        
      case '2':
        progress = !progress;
        break;
        
      case '3':
        treeprint = !treeprint;
        break;
        
      case '4':
        trout = !trout;
        break;
      }
    } else
      printf("Not a possible option!\n");
    countup(&loopcount, 100);
  }
}  /* getoptions */

void reallocsites() 
{
  long i;
  for (i = 0; i < spp; i++) {
    free(y[i]);
    y[i] = (Char *)Malloc(sites*sizeof(Char));
  }
  free(weight);
  free(alias);
  free(aliasweight);
  weight = (steptr)Malloc((sites+1)*sizeof(long));
  alias = (steptr)Malloc((sites+1)*sizeof(long));
  aliasweight = (steptr)Malloc((sites+1)*sizeof(long));
}


void allocrest()
{
  long i;

  y = (Char **)Malloc(spp*sizeof(Char *));
  for (i = 0; i < spp; i++)
    y[i] = (Char *)Malloc(sites*sizeof(Char));
  nayme = (naym *)Malloc(spp*sizeof(naym));
  enterorder = (long *)Malloc(spp*sizeof(long));
  weight = (steptr)Malloc((sites+1)*sizeof(long));
  alias = (steptr)Malloc((sites+1)*sizeof(long));
  aliasweight = (steptr)Malloc((sites+1)*sizeof(long));
}  /* allocrest */


void freelrsaves()
{
  long i,j;
  for ( i = 0 ; i < NLRSAVES ; i++ ) {
    for (j = 0; j < endsite; j++)
      free(lrsaves[i]->x2[j]);
    free(lrsaves[i]->x2);
    free(lrsaves[i]->underflows);
    free(lrsaves[i]);
  }
  free(lrsaves);
}


void resetlrsaves() 
{
  freelrsaves();
  alloclrsaves();
}


void alloclrsaves()
{
  long i,j;
  lrsaves = Malloc(NLRSAVES * sizeof(node*));
  for ( i = 0 ; i < NLRSAVES ; i++ ) {
    lrsaves[i] = Malloc(sizeof(node));
    lrsaves[i]->x2 = Malloc((endsite + 1)*sizeof(sitelike2));
    for ( j = 0 ; j <= endsite ; j++ ) {
      lrsaves[i]->x2[j] = Malloc((sitelength + 1) * sizeof(double));
    }
  }
} /* alloclrsaves */


void setuppie()
{
  /* set up equilibrium probabilities of being a given
     number of bases away from a restriction site */
  long i;
  double sum;

  pie = init_sitelike(sitelength);
  pie[0] = 1.0;
  sum = pie[0];
  for (i = 1; i <= sitelength; i++) {
    pie[i] = 3 * pie[i - 1] * (sitelength - i + 1) / i;
    sum += pie[i];
  }
  for (i = 0; i <= sitelength; i++)
    pie[i] /= sum;
}  /* setuppie */


void doinit()
{
  /* initializes variables */
  long i,j;

  restml_inputnumbers();
  getoptions();
  if (!usertree)
    nonodes2--;
  if (printdata)
    fprintf(outfile, "%4ld Species, %4ld Sites,%4ld Enzymes\n",
            spp, sites, enzymes);
  tempmatrix = Malloc(NTEMPMATS * sizeof(transmatrix));
  for ( i = 0 ; i < NTEMPMATS ; i++ ) {
    tempmatrix[i] = Malloc((sitelength+1) * sizeof(double *));
    for ( j = 0 ; j <= sitelength ; j++)
      tempmatrix[i][j] = (double *)Malloc((sitelength+1) * sizeof(double));
  }
  tempslope = (transmatrix)Malloc((sitelength+1) * sizeof(double *));
  for (i=0; i<=sitelength; i++)
    tempslope[i] = (double *)Malloc((sitelength+1) * sizeof(double));
  tempcurve = (transmatrix)Malloc((sitelength+1) * sizeof(double *));
  for (i=0; i<=sitelength; i++)
    tempcurve[i] = (double *)Malloc((sitelength+1) * sizeof(double));
  setuppie();
  alloctrans(&curtree, nonodes2, sitelength);
  alloctree(&curtree.nodep, nonodes2, usertree);
  allocrest();
  if (usertree)
    return;
  alloctrans(&bestree, nonodes2, sitelength);
  alloctree(&bestree.nodep, nonodes2, 0);
  alloctrans(&priortree, nonodes2, sitelength);
  alloctree(&priortree.nodep, nonodes2, 0);
  if (njumble == 1) return;
  alloctrans(&bestree2, nonodes2, sitelength);
  alloctree(&bestree2.nodep, nonodes2, 0);
}  /* doinit */


void cleanup() {
  long i, j;

  for (i = 0; i < NTEMPMATS; i++) {
    for (j = 0; j <= sitelength; j++)
      free(tempmatrix[i][j]);
    free(tempmatrix[i]);
  }
  free(tempmatrix);
  tempmatrix = NULL;
  for (i = 0; i <= sitelength; i++) {
    free(tempslope[i]);
    free(tempcurve[i]);
  }
  free(tempslope);
  tempslope = NULL;
  free(tempcurve);
  tempcurve = NULL;
}


void inputoptions()
{
  /* read the options information */
  Char ch;
  long i, extranum, cursp, curst, curenz;

  if (!firstset) {
    if (eoln(infile))
      scan_eoln(infile);
    fscanf(infile, "%ld%ld%ld", &cursp, &curst, &curenz);
    if (cursp != spp) {
      printf("\nERROR: INCONSISTENT NUMBER OF SPECIES IN DATA SET %4ld\n",
             ith);
      exxit(-1);
    }
    if (curenz != enzymes) {
      printf("\nERROR: INCONSISTENT NUMBER OF ENZYMES IN DATA SET %4ld\n",
             ith);
      exxit(-1);
    }
    sites = curst;
  }
  if ( !firstset )
    reallocsites();

  for (i = 1; i <= sites; i++)
    weight[i] = 1;
  weightsum = sites;
  extranum = 0;
  readoptions(&extranum, "W");
  for (i = 1; i <= extranum; i++) {
    matchoptions(&ch, "W");
    if (ch == 'W')
      inputweights2(1, sites+1, &weightsum, weight, &weights, "RESTML");
  }
  fprintf(outfile, "\n  Recognition sequences all%2ld bases long\n",
          sitelength);
  if (trunc8)
    fprintf(outfile,
      "\nSites absent from all species are assumed to have been omitted\n\n");
  if (weights)
    printweights(outfile, 1, sites, weight, "Sites");
}  /* inputoptions */


void restml_inputdata()
{
  /* read the species and sites data */
  long i, j, k, l, sitesread, sitesnew=0;
  Char ch;
  boolean allread, done;

  if (printdata)
    putc('\n', outfile);
  j = nmlngth + (sites + (sites - 1) / 10) / 2 - 5;
  if (j < nmlngth - 1)
    j = nmlngth - 1;
  if (j > 39)
    j = 39;
  if (printdata) {
    fprintf(outfile, "Name");
    for (i = 1; i <= j; i++)
      putc(' ', outfile);
    fprintf(outfile, "Sites\n");
    fprintf(outfile, "----");
    for (i = 1; i <= j; i++)
      putc(' ', outfile);
    fprintf(outfile, "-----\n\n");
  }
  sitesread = 0;
  allread = false;
  while (!(allread)) {
    /* eat white space -- if the separator line has spaces on it*/
    do {
      ch = gettc(infile);
    } while (ch == ' ' || ch == '\t');
    ungetc(ch, infile);
    if (eoln(infile))
      scan_eoln(infile);
    i = 1;
    while (i <= spp ) {
      if ((interleaved && sitesread == 0) || !interleaved)
        initname(i - 1);
      if (interleaved)
        j = sitesread;
      else
        j = 0;
      done = false;
      while (!done && !eoff(infile)) {
        if (interleaved)
          done = true;
        while (j < sites && !(eoln(infile) || eoff(infile))) {
          ch = gettc(infile);
          if (ch == '\n' || ch == '\t')
            ch = ' ';
          if (ch == ' ')
            continue;
          uppercase(&ch);
          if (ch != '1' && ch != '0' && ch != '+' && ch != '-' && ch != '?') {
            printf(" ERROR: Bad symbol %c", ch);
            printf(" at position %ld of species %ld\n", j+1, i);
            exxit(-1);
          }
          if (ch == '1')
            ch = '+';
          if (ch == '0')
            ch = '-';
          j++;
          y[i - 1][j - 1] = ch;
        }
        if (interleaved)
          continue;
        if (j < sites) 
          scan_eoln(infile);
        else if (j == sites)
          done = true;
      }
      if (interleaved && i == 1)
        sitesnew = j;
      scan_eoln(infile);
      if ((interleaved && j != sitesnew ) || ((!interleaved) && j != sites)){
        printf("ERROR: SEQUENCES OUT OF ALIGNMENT\n");
        exxit(-1);}
      i++;
    }
    if (interleaved) {
      sitesread = sitesnew;
      allread = (sitesread == sites);
    } else
      allread = (i > spp);
  }
  if (printdata) {
    for (i = 1; i <= ((sites - 1) / 60 + 1); i++) {
      for (j = 0; j < spp; j++) {
        for (k = 0; k < nmlngth; k++)
          putc(nayme[j][k], outfile);
        fprintf(outfile, "   ");
        l = i * 60;
        if (l > sites)
          l = sites;
        for (k = (i - 1) * 60 + 1; k <= l; k++) {
          putc(y[j][k - 1], outfile);
          if (k % 10 == 0 && k % 60 != 0)
            putc(' ', outfile);
        }
        putc('\n', outfile);
      }
      putc('\n', outfile);
    }
    putc('\n', outfile);
  }
  putc('\n', outfile);
}  /* restml_inputdata */


void restml_sitesort()
{
  /* Shell sort keeping alias, aliasweight in same order */
  long gap, i, j, jj, jg, k, itemp;
  boolean flip, tied;

  gap = sites / 2;
  while (gap > 0) {
    for (i = gap + 1; i <= sites; i++) {
      j = i - gap;
      flip = true;
      while (j > 0 && flip) {
        jj = alias[j];
        jg = alias[j + gap];
        flip = false;
        tied = true;
        k = 1;
        while (k <= spp && tied) {
          flip = (y[k - 1][jj - 1] > y[k - 1][jg - 1]);
          tied = (tied && y[k - 1][jj - 1] == y[k - 1][jg - 1]);
          k++;
        }
        if (tied) {
          aliasweight[j] += aliasweight[j + gap];
          aliasweight[j + gap] = 0;
        }
        if (!flip)
          break;
        itemp = alias[j];
        alias[j] = alias[j + gap];
        alias[j + gap] = itemp;
        itemp = aliasweight[j];
        aliasweight[j] = aliasweight[j + gap];
        aliasweight[j + gap] = itemp;
        j -= gap;
      }
    }
    gap /= 2;
  }
}  /* restml_sitesort */


void restml_sitecombine()
{
  /* combine sites that have identical patterns */
  long i, j, k;
  boolean tied;

  i = 1;
  while (i < sites) {
    j = i + 1;
    tied = true;
    while (j <= sites && tied) {
      k = 1;
      while (k <= spp && tied) {
        tied = (tied && y[k - 1][alias[i] - 1] == y[k - 1][alias[j] - 1]);
        k++;
      }
      if (tied && aliasweight[j] > 0) {
        aliasweight[i] += aliasweight[j];
        aliasweight[j] = 0;
        alias[j] = alias[i];
      }
      j++;
    }
    i = j - 1;
  }
}  /* restml_sitecombine */


void makeweights()
{
  /* make up weights vector to avoid duplicate computations */
  long i;

  for (i = 1; i <= sites; i++) {
    alias[i] = i;
    aliasweight[i] = weight[i];
  }
  restml_sitesort();
  restml_sitecombine();
  sitescrunch2(sites + 1, 2, 3, aliasweight);
  for (i = 1; i <= sites; i++) {
    weight[i] = aliasweight[i];
    if (weight[i] > 0)
      endsite = i;
  }
  weight[0] = 1;
}  /* makeweights */


void restml_makevalues()
{
  /* set up fractional likelihoods at tips */
  long i, j, k, l, m;
  boolean found;

  for (k = 1; k <= endsite; k++) {
    j = alias[k];
    for (i = 0; i < spp; i++) {
      for (l = 0; l <= sitelength; l++)
        curtree.nodep[i]->x2[k][l] = 1.0;
      switch (y[i][j - 1]) {

      case '+':
        for (m = 1; m <= sitelength; m++)
          curtree.nodep[i]->x2[k][m] = 0.0;
        break;

      case '-':
        curtree.nodep[i]->x2[k][0] = 0.0;
        break;

      case '?':
        /* blank case */
        break;
      }
    }
  }
  for (i = 0; i < spp; i++) {
    for (k = 1; k <= sitelength; k++)
      curtree.nodep[i]->x2[0][k] = 1.0;
    curtree.nodep[i]->x2[0][0] = 0.0;
  }
  if (trunc8)
    return;
  found = false;
  i = 1;
  while (!found && i <= endsite) {
    found = true;
    for (k = 0; k < spp; k++)
      found = (found && y[k][alias[i] - 1] == '-');
    if (!found)
      i++;
  }
  if (found) {
    weightsum += (enzymes - 1) * weight[i];
    weight[i] *= enzymes;
  }
}  /* restml_makevalues */


void getinput()
{
  /* reads the input data */
  inputoptions();
  restml_inputdata();
  if ( !firstset ) freelrsaves();
  makeweights();
  alloclrsaves();
  if (!usertree) {
    setuptree2(&curtree);
    setuptree2(&priortree);
    setuptree2(&bestree);
    if (njumble > 1) 
      setuptree2(&bestree2);
  }
  allocx2(nonodes2, endsite+1, sitelength, curtree.nodep, usertree);
  if (!usertree) {
    allocx2(nonodes2, endsite+1, sitelength, priortree.nodep, 0);
    allocx2(nonodes2, endsite+1, sitelength, bestree.nodep, 0);
    if (njumble > 1)
      allocx2(nonodes2, endsite+1, sitelength, bestree2.nodep, 0);
  }
  restml_makevalues();
}  /* getinput */


void copymatrix(transmatrix tomat, transmatrix frommat)
{
  /* copy a matrix the size of transition matrix */
  int i,j;

  for (i=0;i<=sitelength;++i){
    for (j=0;j<=sitelength;++j)
      tomat[i][j] = frommat[i][j];
  }
} /* copymatrix */


void maketrans(double p, boolean nr)
{
  /* make transition matrix, product matrix with change
     probability p.  Put the results in tempmatrix, tempslope, tempcurve */
  long i, j, k, m1, m2;
  double sump, sums=0, sumc=0, pover3, pijk, term;
  sitelike2 binom1, binom2;

  binom1 = init_sitelike(sitelength);
  binom2 = init_sitelike(sitelength);
  pover3 = p / 3;
  for (i = 0; i <= sitelength; i++) {
    if (p > 1.0 - epsilon)
      p = 1.0 - epsilon;
    if (p < epsilon)
      p = epsilon;
    binom1[0] = exp((sitelength - i) * log(1 - p));
    for (k = 1; k <= sitelength - i; k++)
      binom1[k] = binom1[k - 1] * (p / (1 - p)) * (sitelength - i - k + 1) / k;
    binom2[0] = exp(i * log(1 - pover3));
    for (k = 1; k <= i; k++)
      binom2[k] = binom2[k - 1] * (pover3 / (1 - pover3)) * (i - k + 1) / k;
    for (j = 0; j <= sitelength; ++j) {
      sump = 0.0;
      if (nr) {
        sums = 0.0;
        sumc = 0.0;
      }
      if (i - j > 0)
        m1 = i - j;
      else
        m1 = 0;
      if (sitelength - j < i)
        m2 = sitelength - j;
      else
        m2 = i;
      for (k = m1; k <= m2; k++) {
        pijk = binom1[j - i + k] * binom2[k];
        sump += pijk;
        if (nr) {
          term = (j-i+2*k)/p - (sitelength-j-k)/(1.0-p) - (i-k)/(3.0-p);
          sums += pijk * term;
          sumc += pijk * (term * term
                            - (j-i+2*k)/(p*p)
                            - (sitelength-j-k)/((1.0-p)*(1.0-p))
                            - (i-k)/((3.0-p)*(3.0-p)) );
        }
      }
      tempmatrix[0][i][j] = sump;
      if (nr) {
        tempslope[i][j] = sums;
        tempcurve[i][j] = sumc;
      }
    }
  }
  free_sitelike(binom1);
  free_sitelike(binom2);
}  /* maketrans */


void branchtrans(long i, double p)
{
  /* make branch transition matrix for branch i with probability of change p */
  boolean nr;

  nr = false;
  maketrans(p, nr);
  copymatrix(curtree.trans[i - 1], tempmatrix[0]);
}  /* branchtrans */


double evaluate(tree *tr, node *p)
{
  /* evaluates the likelihood, using info. at one branch */
  double sum, sum2, y, liketerm, like0, lnlike0=0, term;
  long i, j, k,branchnum;
  node *q;
  sitelike2 x1, x2;

  x1 = init_sitelike(sitelength);
  x2 = init_sitelike(sitelength);
  sum = 0.0;
  q = p->back;

  nuview(p);
  nuview(q);

  y = p->v;
  branchnum = p->branchnum;
  copy_sitelike(x1,p->x2[0],sitelength);
  copy_sitelike(x2,q->x2[0],sitelength);
  if (trunc8) {
    like0 = 0.0;
    for (j = 0; j <= sitelength; j++) {
      liketerm = pie[j] * x1[j];
      for (k = 0; k <= sitelength; k++)
        like0 += liketerm * tr->trans[branchnum-1][j][k] * x2[k];
    }
    lnlike0 = log(enzymes * (1.0 - like0));
  }
  for (i = 1; i <= endsite; i++) {
    copy_sitelike(x1,p->x2[i],sitelength);
    copy_sitelike(x2,q->x2[i],sitelength);
    sum2 = 0.0;
    for (j = 0; j <= sitelength; j++) {
      liketerm = pie[j] * x1[j];
      for (k = 0; k <= sitelength; k++)
        sum2 += liketerm * tr->trans[branchnum-1][j][k] * x2[k];
    }
    term = log(sum2);
    if (trunc8)
      term -= lnlike0;
    if (usertree && (which <= shimotrees))
      l0gf[which - 1][i - 1] = term;
    sum += weight[i] * term;
  }
/* *** debug  put a variable "saveit" in evaluate as third argument as to
   whether to save the KHT suff */
  if (usertree) {
    if(which <= shimotrees)
      l0gl[which - 1] = sum;
    if (which == 1) {
      maxwhich = 1;
      maxlogl = sum;
    } else if (sum > maxlogl) {
      maxwhich = which;
      maxlogl = sum;
    }
  }
  tr->likelihood = sum;
  free_sitelike(x1);
  free_sitelike(x2);
  return sum;
}  /* evaluate */


boolean nuview(node *p)
{
  /* recompute fractional likelihoods for one part of tree */
  long i, j, k, lowlim;
  double sumq;
  node *q, *s;
  sitelike2 xq, xr, xp;
  double **tempq = NULL;
  double *tq     = NULL;

  if (p->tip)
    return false;
  xq = init_sitelike(sitelength);
  xr = init_sitelike(sitelength);
  xp = init_sitelike(sitelength);
  for (s = p->next; s != p; s = s->next) {
    if ( nuview(s->back) )
      p->initialized = false;
  }

  if (p->initialized)
    return false;

  lowlim = trunc8 ? 0 : 1;

  for (i = lowlim; i <= endsite; i++) {
    xp = p->x2[i];
    for (j = 0; j <= sitelength; j++)
      xp[j] = 1.0;
    for (s = p->next; s != p; s = s->next) {
      q = s->back;
      tempq = curtree.trans[q->branchnum - 1];
      xq = q->x2[i];
      for (j = 0; j <= sitelength; j++) {
        sumq = 0.0;
        tq = tempq[j];
        for (k = 0; k <= sitelength; k++)
          sumq += tq[k] * xq[k];
        xp[j] *= sumq;
      }
    }
  }

  return true;
}  /* nuview */


void makenewv(node *p)
{
  /* Newton-Raphson algorithm improvement of a branch length */
  long i, j, k, lowlim, it, ite;
  double sum, sums, sumc, like, slope, curve, liketerm, liket,
         y, yold=0, yorig, like0=0, slope0=0, curve0=0, oldlike=0, temp;
  boolean done, nr, firsttime, better;
  node *q;
  sitelike2 xx1, xx2;
  double *tm, *ts, *tc;

  q = p->back;
  y = p->v;
  yorig = y;
  if (trunc8)
    lowlim = 0;
  else
    lowlim = 1;
  done = false;
  nr = true;
  firsttime = true;
  it = 1;
  ite = 0;
  while ((it < iterations) && (ite < 20) && (!done)) {
    like = 0.0;
    slope = 0.0;
    curve = 0.0;
    maketrans(y, nr);
    for (i = lowlim; i <= endsite; i++) {
      xx1 = p->x2[i];
      xx2 = q->x2[i];
      sum = 0.0;
      sums = 0.0;
      sumc = 0.0;
      for (j = 0; j <= sitelength; j++) {
        liket = xx1[j] * pie[j];
        tm = tempmatrix[0][j];
        ts = tempslope[j];
        tc = tempcurve[j];
        for (k = 0; k <= sitelength; k++) {
          liketerm = liket * xx2[k];
          sum += tm[k] * liketerm;
          sums += ts[k] * liketerm;
          sumc += tc[k] * liketerm;
        }
      }
      if (i == 0) {
        like0 = sum;
        slope0 = sums;
        curve0 = sumc;
      } else {
        like += weight[i] * log(sum);
        slope += weight[i] * sums/sum;
        temp = sums/sum;
        curve += weight[i] * (sumc/sum-temp*temp);
      }
    }
    if (trunc8 && fabs(like0 - 1.0) > 1.0e-10) {
      like -= weightsum * log(enzymes * (1.0 - like0));
      slope += weightsum * slope0 /(1.0 - like0);
      curve += weightsum * (curve0 /(1.0 - like0)
                            + slope0*slope0/((1.0 - like0)*(1.0 - like0)));
    }
    better = false;
    if (firsttime) {
      yold = y;
      oldlike = like;
      firsttime = false;
      better = true;
    } else {
      if (like > oldlike) {
        yold = y;
        oldlike = like;
        better = true;
        it++;
      }
    }
    if (better) {
      y = y + slope/fabs(curve);
      if (y < epsilon)
        y = 10.0 * epsilon;
      if (y > 0.75)
        y = 0.75;
    } else {
        if (fabs(y - yold) < epsilon)
          ite = 20;
        y = (y + yold) / 2.0;
    }
    ite++;
    done = fabs(y-yold) < epsilon;
  }
  smoothed = (fabs(yold-yorig) < epsilon) && (yorig > 1000.0*epsilon);
  p->v = yold;
  q->v = yold;
  branchtrans(p->branchnum, yold);
  curtree.likelihood = oldlike;
}  /* makenewv */


void update(node *p)
{
  /* improve branch length and views for one branch */
  nuview(p);
  nuview(p->back);

  if ( !(usertree && lengths) ) {
    makenewv(p);
    if (smoothit ) {
      inittrav(p);
      inittrav(p->back);
    }
    else {
     if (inserting && !p->tip) {
        p->next->initialized = false;
        p->next->next->initialized = false;
      }
    }
  }
}  /* update */


void smooth(node *p)
{
  /* update nodes throughout the tree, recursively */

  smoothed = false;
  update(p);
  if (!p->tip) {
    if (smoothit && !smoothed) {
      smooth(p->next->back);
    }
    if (smoothit && !smoothed) {
      smooth(p->next->next->back);
    }
  }
}  /* smooth */


void insert_(node *p, node *q)
{
  /* insert a subtree into a branch, improve lengths in tree */
  long i;
  node *r;

  r = p->next->next;
  hookup(r, q->back);
  hookup(p->next, q);
  if (q->v >= 0.75)
    q->v = 0.75;
  else
    q->v = 0.75 * (1 - sqrt(1 - 1.333333 * q->v));
  if ( q->v < epsilon)
    q->v = epsilon;

  q->back->v = q->v;
  r->v = q->v;
  r->back->v = r->v;
 
  set_branchnum(q->back, q->branchnum);
  set_branchnum(r, get_trans(&curtree));
  set_branchnum(r->back, r->branchnum);
  branchtrans(q->branchnum, q->v);
  branchtrans(r->branchnum, r->v);

  if ( smoothit ) {
    inittrav(p);
    inittrav(p->back);
  } 
  p->initialized = false;

  i = 1;
  inserting = true;
  while (i <= smoothings) {
    smooth(p);
    if (!p->tip) {
      smooth (p->next->back);
      smooth (p->next->next->back);
    }
    i++;
  }
  inserting = false;
}  /* insert_ */


void restml_re_move(node **p, node **q)
{
  /* remove p and record in q where it was */
  long i;

  *q = (*p)->next->back;
  hookup(*q, (*p)->next->next->back);
  free_trans(&curtree,(*q)->back->branchnum);
  set_branchnum((*q)->back, (*q)->branchnum);
  (*q)->v = 0.75*(1 - (1 - 1.333333*(*q)->v) * (1 - 1.333333*(*p)->next->v));

  if ( (*q)->v > 1 - epsilon)
    (*q)->v = 1 - epsilon;
  else if ( (*q)->v < epsilon)
    (*q)->v = epsilon;
                  
  (*q)->back->v = (*q)->v;
  branchtrans((*q)->branchnum, (*q)->v);

  (*p)->next->back = NULL;
  (*p)->next->next->back = NULL;
 
  if ( smoothit ) {
    inittrav((*q)->back);
    inittrav(*q);
  }
  
  if ( smoothit ) {
    for ( i = 0 ; i < smoothings ; i++ ) {
      smooth(*q);
      smooth((*q)->back);
    }
  }
  else ( smooth(*q));
}  /* restml_re_move */


void restml_copynode(node *c, node *d)
{
  /* copy a node */
  long i;

  set_branchnum(d, c->branchnum);
  
  for ( i = 0 ; i <= endsite ; i++)
    copy_sitelike(d->x2[i],c->x2[i],sitelength);
  d->v = c->v;
  d->iter = c->iter;
  d->xcoord = c->xcoord;
  d->ycoord = c->ycoord;
  d->ymin = c->ymin;
  d->ymax = c->ymax;
  d->initialized = c->initialized;
}  /* restml_copynode */


void restml_copy_(tree *a, tree *b)
{
  /* copy tree a to tree b */
  long i,j;
  node *p, *q;

  for (i = 0; i < spp; i++) {
    restml_copynode(a->nodep[i], b->nodep[i]);
    if (a->nodep[i]->back) {
      if (a->nodep[i]->back == a->nodep[a->nodep[i]->back->index - 1])
        b->nodep[i]->back = b->nodep[a->nodep[i]->back->index - 1];
      else if (a->nodep[i]->back
                == a->nodep[a->nodep[i]->back->index - 1]->next)
          b->nodep[i]->back = b->nodep[a->nodep[i]->back->index - 1]->next;
        else
          b->nodep[i]->back
                         = b->nodep[a->nodep[i]->back->index - 1]->next->next;
    }
    else b->nodep[i]->back = NULL;
  }
  for (i = spp; i < nonodes2; i++) {
    p = a->nodep[i];
    q = b->nodep[i];
    for (j = 1; j <= 3; j++) {
      restml_copynode(p, q);
      if (p->back) {
        if (p->back == a->nodep[p->back->index - 1])
          q->back = b->nodep[p->back->index - 1];
        else if (p->back == a->nodep[p->back->index - 1]->next)
          q->back = b->nodep[p->back->index - 1]->next;
        else
          q->back = b->nodep[p->back->index - 1]->next->next;
      }
      else
        q->back = NULL;
      p = p->next;
      q = q->next;
    }
  }
  b->likelihood = a->likelihood;
  for (i=0;i<nonodes2;++i)
    copymatrix(b->trans[i],a->trans[i]);
  b->transindex = a->transindex;
  memcpy(b->freetrans,a->freetrans,nonodes*sizeof(long));

  b->start = a->start;
}  /* restml_copy */


void buildnewtip(long m,tree *tr)
{
  /* set up a new tip and interior node it is connected to */
  node *p;
  long i, j;

  p = tr->nodep[nextsp + spp - 3];
  for (i = 0; i <= endsite; i++) {
    for (j = 0; j < sitelength; j++) { /* trunc8 */
      p->x2[i][j] = 1.0; 
      p->next->x2[i][j] = 1.0; 
      p->next->next->x2[i][j] = 1.0; 
    }
  }

  hookup(tr->nodep[m - 1], p);
  p->v = initialv;
  p->back->v = initialv;
  set_branchnum(p, get_trans(tr));
  set_branchnum(p->back, p->branchnum);
  branchtrans(p->branchnum, initialv);
}  /* buildnewtip */


void buildsimpletree(tree *tr)
{
  /* set up and adjust branch lengths of a three-species tree */
  long branch;

  hookup(tr->nodep[enterorder[0] - 1], tr->nodep[enterorder[1] - 1]);
  tr->nodep[enterorder[0] - 1]->v = initialv;
  tr->nodep[enterorder[1] - 1]->v = initialv;
  branchtrans(enterorder[1], initialv);
  branch = get_trans(tr);
  set_branchnum(tr->nodep[enterorder[0] - 1], branch);
  set_branchnum(tr->nodep[enterorder[1] - 1], branch);
  buildnewtip(enterorder[2], tr);
  insert_(tr->nodep[enterorder[2] - 1]->back, tr->nodep[enterorder[1] - 1]);
  tr->start = tr->nodep[enterorder[2]-1]->back;
}  /* buildsimpletree */


void addtraverse(node *p, node *q, boolean contin)
{
  /* try adding p at q, proceed recursively through tree */
  double like, vsave = 0;
  node *qback =NULL;

  if (!smoothit) {
    copymatrix (tempmatrix[1], curtree.trans[q->branchnum - 1]);
    vsave = q->v;
    qback = q->back;
  }
  insert_(p, q);
  like = evaluate(&curtree, p);
  if (like > bestyet) {
    bestyet = like;
    if (smoothit) {
      restml_copy_(&curtree, &bestree);
      addwhere = q;
    }
    else
      qwhere = q;
    succeeded = true;
  }
  if (smoothit)
    restml_copy_(&priortree, &curtree);
  else {
    hookup (q, qback);
    q->v = vsave;
    q->back->v = vsave;
    free_trans(&curtree,q->back->branchnum);
    set_branchnum(q->back, q->branchnum);
    copymatrix (curtree.trans[q->branchnum - 1], tempmatrix[1]);
#ifdef DEBUG
    evaluate(&curtree, curtree.start);
    assert( close(bestyet, curtree.likelihood) );
#endif
    /* curtree.likelihood = bestyet; */
    evaluate(&curtree, curtree.start);
  }
  if (!q->tip && contin) {
    /* assumes bifurcation (OK) */
    addtraverse(p, q->next->back, contin);
    addtraverse(p, q->next->next->back, contin);
  }
  if ( contin && q == curtree.root ) {
    /* FIXME!! curtree.root->back == NULL?  curtree.root == NULL? */
    addtraverse(p,q->back,contin);
  }
}  /* addtraverse */


void globrearrange() 
{
  /* does global rearrangements */
  tree globtree;
  tree oldtree;
  int i,j,k,l,num_sibs,num_sibs2;
  node *where,*sib_ptr,*sib_ptr2;
  double oldbestyet = curtree.likelihood;
  int success = false;
  printf("\n   ");
  alloctree(&globtree.nodep,nonodes2,0);
  alloctree(&oldtree.nodep,nonodes2,0);
  alloctrans(&globtree, nonodes2, sitelength);
  alloctrans(&oldtree, nonodes2, sitelength);
  setuptree2(&globtree);
  setuptree2(&oldtree);
  allocx2(nonodes2, endsite + 1, sitelength,globtree.nodep, 0);
  allocx2(nonodes2, endsite + 1, sitelength,oldtree.nodep, 0);
  restml_copy_(&curtree,&globtree);
  restml_copy_(&curtree,&oldtree);
  bestyet = curtree.likelihood;
  for ( i = spp ; i < nonodes2 ; i++ ) {
    num_sibs = count_sibs(curtree.nodep[i]);
    sib_ptr  = curtree.nodep[i];
    if ( (i - spp) % (( nonodes2 / 72 ) + 1 ) == 0 )
      putchar('.');
    fflush(stdout);
    for ( j = 0 ; j <= num_sibs ; j++ ) {
      restml_re_move(&sib_ptr,&where);
      restml_copy_(&curtree,&priortree);
      qwhere = where;
      
      if (where->tip) {
        restml_copy_(&oldtree,&curtree);
        restml_copy_(&oldtree,&bestree);
        sib_ptr=sib_ptr->next;
        continue;
      }
      else num_sibs2 = count_sibs(where);
      sib_ptr2 = where;
      for ( k = 0 ; k < num_sibs2 ; k++ ) {
        addwhere = NULL;
        addtraverse(sib_ptr,sib_ptr2->back,true);
        if ( !smoothit ) {
          if (succeeded && qwhere != where && qwhere != where->back) {
            insert_(sib_ptr,qwhere);
            smoothit = true;
            for (l = 1; l<=smoothings; l++) {
              smooth (where);
              smooth (where->back);
            }
            smoothit = false;
            success = true;
            restml_copy_(&curtree,&globtree);
          }
          restml_copy_(&priortree,&curtree);
        }
        else if ( addwhere && where != addwhere && where->back != addwhere
              && bestyet > globtree.likelihood) {
            restml_copy_(&bestree,&globtree);
            success = true;
        }
        sib_ptr2 = sib_ptr2->next;
      } 
      restml_copy_(&oldtree,&curtree);
      restml_copy_(&oldtree,&bestree);
      sib_ptr = sib_ptr->next;
    }
  }
  restml_copy_(&globtree,&curtree);
  restml_copy_(&globtree,&bestree);
  if (success && globtree.likelihood > oldbestyet)  {
    succeeded = true;
  }
  else  {
    succeeded = false;
  }
  bestyet = globtree.likelihood;
  freetrans(&globtree.trans, nonodes2, sitelength);
  freetrans(&oldtree.trans, nonodes2, sitelength);
  freetree2(globtree.nodep,nonodes2);
  freetree2(oldtree.nodep,nonodes2);
}

void printnode(node* p);
void printnode(node* p)
{
  if (p->back)
  printf("p->index = %3ld, p->back->index = %3ld, p->branchnum = %3ld,evaluates"
      " to %f\n",p->index,p->back->index,p->branchnum,evaluate(&curtree,p));
  
  else 
  printf("p->index = %3ld, p->back->index =none, p->branchnum = %3ld,evaluates"
      " to nothing\n",p->index,p->branchnum);
}

void printvals(void);
void printvals(void) 
{
  int i;
  node* p;

  for ( i = 0 ; i < nextsp ; i++ )  {
    p = curtree.nodep[i];
    printnode(p);
  }
  for ( i = spp ; i <= spp + nextsp - 3 ; i++ ) {
    p = curtree.nodep[i];
    printnode(p);
    printnode(p->next);
    printnode(p->next->next);
  }
}

void rearrange(node *p, node *pp)
{
  /* rearranges the tree locally */
  long i;
  node *q;
  node *r;
  node *rnb;
  node *rnnb;

  if (p->tip) return;
  if (p->back->tip) {
    rearrange(p->next->back, p);
    rearrange(p->next->next->back, p);

    return;

  } else /* if !p->tip && !p->back->tip */ {
    /*
    evaluate(&curtree, curtree.start);
    bestyet = curtree.likelihood;
    */

    if (p->back->next != pp)
      r = p->back->next;
    else
      r = p->back->next->next;

    if (smoothit) {
      /* Copy the whole tree, because we may change all lengths */
      restml_copy_(&curtree, &bestree);
      restml_re_move(&r, &q);
      nuview(p->next);
      nuview(p->next->next);
      restml_copy_(&curtree, &priortree);
      addtraverse(r, p->next->back, false);
      addtraverse(r, p->next->next->back, false);
      restml_copy_(&bestree, &curtree);
    } else {
      /* Save node data and matrices, so we can undo */
      rnb = r->next->back;
      rnnb = r->next->next->back;

      restml_copynode(r,lrsaves[0]);
      restml_copynode(r->next,lrsaves[1]);
      restml_copynode(r->next->next,lrsaves[2]);
      restml_copynode(p->next,lrsaves[3]);
      restml_copynode(p->next->next,lrsaves[4]);

      copymatrix (tempmatrix[2], curtree.trans[r->branchnum - 1]);
      copymatrix (tempmatrix[3], curtree.trans[r->next->branchnum - 1]);
      copymatrix (tempmatrix[4], curtree.trans[r->next->next->branchnum-1]);
      copymatrix (tempmatrix[5], curtree.trans[p->next->branchnum-1]);
      copymatrix (tempmatrix[6], curtree.trans[p->next->next->branchnum-1]);

      restml_re_move(&r, &q);
      nuview(p->next);
      nuview(p->next->next);
      qwhere = q;
      addtraverse(r, p->next->back, false);
      addtraverse(r, p->next->next->back, false);
      if (qwhere == q) {
        hookup(rnb,r->next);
        hookup(rnnb,r->next->next);
        restml_copynode(lrsaves[0],r);
        restml_copynode(lrsaves[1],r->next);
        restml_copynode(lrsaves[2],r->next->next);
        restml_copynode(lrsaves[3],p->next);
        restml_copynode(lrsaves[4],p->next->next);
        r->back->v = r->v;
        r->next->back->v = r->next->v;
        r->next->next->back->v = r->next->next->v;
        p->next->back->v = p->next->v;
        p->next->next->back->v = p->next->next->v;
        set_branchnum(r->back, r->branchnum);
        set_branchnum(r->next->back, r->next->branchnum);
        set_branchnum(p->next->back, p->next->branchnum);
        set_branchnum(p->next->next->back, p->next->next->branchnum);
        copymatrix (curtree.trans[r->branchnum-1], tempmatrix[2]);
        copymatrix (curtree.trans[r->next->branchnum-1], tempmatrix[3]);
        copymatrix (curtree.trans[r->next->next->branchnum-1], tempmatrix[4]);
        copymatrix (curtree.trans[p->next->branchnum-1], tempmatrix[5]);
        copymatrix (curtree.trans[p->next->next->branchnum-1], tempmatrix[6]);
        curtree.likelihood = bestyet;
      } else {
        smoothit = true;
        insert_(r, qwhere);
        
        for (i = 1; i<=smoothings; i++) {
          smooth (r);
          smooth (r->back);
        }
        smoothit = false;
      }
    }
  }
}  /* rearrange */


void restml_coordinates(node *p, double lengthsum, long *tipy,
                double *tipmax, double *x)
{
  /* establishes coordinates of nodes */
  node *q, *first, *last;

  if (p->tip) {
    p->xcoord = (long)(over * lengthsum + 0.5);
    p->ycoord = (*tipy);
    p->ymin = (*tipy);
    p->ymax = (*tipy);
    (*tipy) += down;
    if (lengthsum > (*tipmax))
      (*tipmax) = lengthsum;
    return;
  }
  q = p->next;
  do {
    (*x) = -0.75 * log(1.0 - 1.333333 * q->v);
    restml_coordinates(q->back, lengthsum + (*x),tipy,tipmax,x);
    q = q->next;
  } while ((p == curtree.start || p != q) &&
           (p != curtree.start || p->next != q));
  first = p->next->back;
  q = p;
  while (q->next != p)
    q = q->next;
  last = q->back;
  p->xcoord = (long)(over * lengthsum + 0.5);
  if (p == curtree.start)
    p->ycoord = p->next->next->back->ycoord;
  else
    p->ycoord = (first->ycoord + last->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* restml_coordinates */


void restml_fprintree(FILE *fp)
{
  /* prints out diagram of the tree */
  long tipy,i;
  double scale, tipmax, x;

  putc('\n', fp);
  if (!treeprint)
    return;
  putc('\n', fp);
  tipy = 1;
  tipmax = 0.0;
  restml_coordinates(curtree.start, 0.0, &tipy,&tipmax,&x);
  scale = 1.0 / (tipmax + 1.000);
  for (i = 1; i <= tipy - down; i++)
    fdrawline2(fp, i, scale, &curtree);
  putc('\n', fp);
}  /* restml_fprintree */


void restml_printree()
{
  restml_fprintree(outfile);
}
  


double sigma(node *q, double *sumlr)
{
  /* get 1.95996 * approximate standard error of branch length */
  double sump, sumr, sums, sumc, p, pover3, pijk, Qjk, liketerm, f;
  double  slopef,curvef;
  long i, j, k, m1, m2;
  sitelike2 binom1, binom2;
  transmatrix Prob, slopeP, curveP;
  node *r;
  sitelike2 x1, x2;
  double term, TEMP;

  x1 = init_sitelike(sitelength);
  x2 = init_sitelike(sitelength);
  binom1 = init_sitelike(sitelength);
  binom2 = init_sitelike(sitelength);
  Prob   = (transmatrix)Malloc((sitelength+1) * sizeof(double *));
  slopeP = (transmatrix)Malloc((sitelength+1) * sizeof(double *));
  curveP = (transmatrix)Malloc((sitelength+1) * sizeof(double *));
  for (i=0; i<=sitelength; ++i) {
    Prob[i]   = (double *)Malloc((sitelength+1) * sizeof(double));
    slopeP[i] = (double *)Malloc((sitelength+1) * sizeof(double));
    curveP[i] = (double *)Malloc((sitelength+1) * sizeof(double));
  }
  p = q->v;
  pover3 = p / 3;
  for (i = 0; i <= sitelength; i++) {
    binom1[0] = exp((sitelength - i) * log(1 - p));
    for (k = 1; k <= (sitelength - i); k++)
      binom1[k] = binom1[k - 1] * (p / (1 - p)) * (sitelength - i - k + 1) / k;
    binom2[0] = exp(i * log(1 - pover3));
    for (k = 1; k <= i; k++)
      binom2[k] = binom2[k - 1] * (pover3 / (1 - pover3)) * (i - k + 1) / k;
    for (j = 0; j <= sitelength; j++) {
      sump = 0.0;
      sums = 0.0;
      sumc = 0.0;
      if (i - j > 0)
        m1 = i - j;
      else
        m1 = 0;
      if (sitelength - j < i)
        m2 = sitelength - j;
      else
        m2 = i;
      for (k = m1; k <= m2; k++) {
        pijk = binom1[j - i + k] * binom2[k];
        sump += pijk;
        term = (j-i+2*k)/p - (sitelength-j-k)/(1.0-p) - (i-k)/(3.0-p);
        sums += pijk * term;
        sumc += pijk * (term * term
                          - (j-i+2*k)/(p*p)
                          - (sitelength-j-k)/((1.0-p)*(1.0-p))
                          - (i-k)/((3.0-p)*(3.0-p)) );
      }
      Prob[i][j] = sump;
      slopeP[i][j] = sums;
      curveP[i][j] = sumc;
    }
  }
  (*sumlr) = 0.0;
  sumc = 0.0;
  sums = 0.0;
  r = q->back;
  for (i = 1; i <= endsite; i++) {
    f = 0.0;
    slopef = 0.0;
    curvef = 0.0;
    sumr = 0.0;
    copy_sitelike(x1,q->x2[i],sitelength);
    copy_sitelike(x2,r->x2[i],sitelength);
    for (j = 0; j <= sitelength; j++) {
      liketerm = pie[j] * x1[j];
      sumr += liketerm * x2[j];
      for (k = 0; k <= sitelength; k++) {
        Qjk = liketerm * x2[k];
        f += Qjk * Prob[j][k];
        slopef += Qjk * slopeP[j][k];
        curvef += Qjk * curveP[j][k];
      }
    }
    (*sumlr) += weight[i] * log(f / sumr);
    sums += weight[i] * slopef / f;
    TEMP = slopef / f;
    sumc += weight[i] * (curvef / f - TEMP * TEMP);
  }
  if (trunc8) {
    f = 0.0;
    slopef = 0.0;
    curvef = 0.0;
    sumr = 0.0;
    copy_sitelike(x1,q->x2[0],sitelength);
    copy_sitelike(x2,r->x2[0],sitelength);
    for (j = 0; j <= sitelength; j++) {
      liketerm = pie[j] * x1[j];
      sumr += liketerm * x2[j];
      for (k = 0; k <= sitelength; k++) {
        Qjk = liketerm * x2[k];
        f += Qjk * Prob[j][k];
        slopef += Qjk * slopeP[j][k];
        curvef += Qjk * curveP[j][k];
      }
    }
    (*sumlr) += weightsum * log((1.0 - sumr) / (1.0 - f));
    sums += weightsum * slopef / (1.0 - f);
    TEMP = slopef / (1.0 - f);
    sumc += weightsum * (curvef / (1.0 - f) + TEMP * TEMP);
  }
  for (i=0;i<=sitelength;++i){
    free(Prob[i]);
    free(slopeP[i]);
    free(curveP[i]);
  }
  free(Prob);
  free(slopeP);
  free(curveP);
  free_sitelike(x1);
  free_sitelike(x2);
  free_sitelike(binom1);
  free_sitelike(binom2);
  if (sumc < -1.0e-6)
    return ((-sums - sqrt(sums * sums - 3.841 * sumc)) / sumc);
  else
    return -1.0;
}  /* sigma */


void fdescribe(FILE *fp, node *p)
{
  /* print out information on one branch */
  double sumlr;
  long i;
  node *q;
  double s;
  double realv;

  q = p->back;
  fprintf(fp, "%4ld      ", q->index - spp);
  fprintf(fp, "    ");
  if (p->tip) {
    for (i = 0; i < nmlngth; i++)
      putc(nayme[p->index - 1][i], fp);
  } else
    fprintf(fp, "%4ld      ", p->index - spp);
  if (q->v >= 0.75)
    fprintf(fp, "     infinity");
  else {
    realv = -0.75 * log(1 - 4.0/3.0 * q->v);
    fprintf(fp, "%13.5f", realv);
  }
  if (p->iter) {
    s = sigma(q, &sumlr);
    if (s < 0.0)
      fprintf(fp, "     (     zero,    infinity)");
    else {
      fprintf(fp, "     (");
      if (q->v - s <= 0.0)
        fprintf(fp, "     zero");
      else
        fprintf(fp, "%9.5f", -0.75 * log(1 - 1.333333 * (q->v - s)));
      putc(',', fp);
      if (q->v + s >= 0.75)
        fprintf(fp, "    infinity");
      else
        fprintf(fp, "%12.5f", -0.75 * log(1 - 1.333333 * (q->v + s)));
      putc(')', fp);
      }
    if (sumlr > 1.9205)
      fprintf(fp, " *");
    if (sumlr > 2.995)
      putc('*', fp);
    }
  else
    fprintf(fp, "            (not varied)");
  putc('\n', fp);
  if (!p->tip) {
    for (q = p->next; q != p; q = q->next)
      fdescribe(fp, q->back);
  }
}  /* fdescribe */


void summarize()
{
  /* print out information on branches of tree */
  node *q;

  fprintf(outfile, "\nremember: ");
  if (outgropt)
    fprintf(outfile, "(although rooted by outgroup) ");
  fprintf(outfile, "this is an unrooted tree!\n\n");
  fprintf(outfile, "Ln Likelihood = %11.5f\n\n", curtree.likelihood);
  fprintf(outfile, " \n");
  fprintf(outfile, " Between        And            Length");
  fprintf(outfile, "      Approx. Confidence Limits\n");
  fprintf(outfile, " -------        ---            ------");
  fprintf(outfile, "      ------- ---------- ------\n");
  for (q = curtree.start->next; q != curtree.start; q = q->next)
    fdescribe(outfile, q->back);
  fdescribe(outfile, curtree.start->back);
  fprintf(outfile, "\n     *  = significantly positive, P < 0.05\n");
  fprintf(outfile, "     ** = significantly positive, P < 0.01\n\n\n");
}  /* summarize */


void restml_treeout(node *p)
{
  /* write out file with representation of final tree */
  long i, n, w;
  Char c;
  double x;

  if (p->tip) {
    n = 0;
    for (i = 1; i <= nmlngth; i++) {
      if (nayme[p->index - 1][i - 1] != ' ')
        n = i;
    }
    for (i = 0; i < n; i++) {
      c = nayme[p->index - 1][i];
      if (c == ' ')
        c = '_';
      putc(c, outtree);
    }
    col += n;
  } else {
    putc('(', outtree);
    col++;
    restml_treeout(p->next->back);
    putc(',', outtree);
    col++;
    if (col > 45) {
      putc('\n', outtree);
      col = 0;
    }
    restml_treeout(p->next->next->back);
    if (p == curtree.start) {
      putc(',', outtree);
      col++;
      if (col > 45) {
        putc('\n', outtree);
        col = 0;
      }
      restml_treeout(p->back);
    }
    putc(')', outtree);
    col++;
  }
  if (p->v >= 0.75)
    x = -1.0;
  else
    x = -0.75 * log(1 - 1.333333 * p->v);
  if (x > 0.0)
    w = (long)(0.43429448222 * log(x));
  else if (x == 0.0)
    w = 0;
  else
    w = (long)(0.43429448222 * log(-x)) + 1;
  if (w < 0)
    w = 0;
  if (p == curtree.start)
    fprintf(outtree, ";\n");
  else {
    fprintf(outtree, ":%*.5f", (int)(w + 7), x);
    col += w + 8;
  }
}  /* restml_treeout */


static phenotype2
restml_pheno_new(long endsite, long sitelength)
{
  phenotype2 ret;
  long k, l;

  endsite++;
  ret = (phenotype2)Malloc(endsite*sizeof(sitelike2));
  for (k = 0; k < endsite; k++) {
    ret[k] = Malloc((sitelength + 1) * sizeof(double));
    for (l = 0; l < sitelength; l++)
      ret[k][l] = 1.0;
  }
  
  return ret;
}


static void
restml_pheno_delete(phenotype2 x2)
{
  long k;

  for (k = 0; k < endsite+1; k++)
    free(x2[k]);
  free(x2);
}


void initrestmlnode(node **p, node **grbg, node *q, long len, long nodei,
                        long *ntips, long *parens, initops whichinit,
                        pointarray treenode, pointarray nodep, Char *str, 
                        Char *ch, FILE *intree)
{
  /* initializes a node */
  boolean minusread;
  double valyew, divisor;

  switch (whichinit) {
  case bottom:
    gnu(grbg, p);
    (*p)->index = nodei;
    (*p)->tip = false;
    (*p)->branchnum = 0;
    (*p)->x2 = restml_pheno_new(endsite, sitelength);
    nodep[(*p)->index - 1] = (*p);
    break;
  case nonbottom:
    gnu(grbg, p);
    (*p)->x2 = restml_pheno_new(endsite, sitelength);
    (*p)->index = nodei;
    break;
  case tip:
    match_names_to_data (str, nodep, p, spp);
    break;
  case iter:
    (*p)->initialized = false;
    (*p)->v = initialv;
    (*p)->iter = true;
    if ((*p)->back != NULL){
      (*p)->back->iter = true;
      (*p)->back->v = initialv;  
      (*p)->back->initialized = false;
    }
    break;
  case length:
    processlength(&valyew, &divisor, ch, &minusread, intree, parens);
    (*p)->v = valyew / divisor;
    (*p)->iter = false;
    if ((*p)->back != NULL) {
      (*p)->back->v = (*p)->v;
      (*p)->back->iter = false;
    }
    break;
  case hsnolength:
    break;
  default:        /* cases hslength, treewt, unittrwt */
    break;
  }
} /* initrestmlnode */


static void
restml_unroot(node* root, node** nodep, long nonodes) 
{
  node *p,*r,*q;
  double newl;
  long i;
  long numsibs;
  
  numsibs = count_sibs(root);

  if ( numsibs > 2 ) {
    q = root;
    r = root;
    while (!(q->next == root))
      q = q->next;
    q->next = root->next;

    /* FIXME?
    for(i=0 ; i < endsite ; i++){
      free(r->x[i]);
      r->x[i] = NULL;
    }
    free(r->x);
    r->x = NULL;
    */

    chuck(&grbg, r);
    curtree.nodep[spp] = q;
  } else { /* Bifurcating root - remove entire root fork */
    /* Join v on each side of root */
    newl = root->next->v + root->next->next->v;
    root->next->back->v = newl;
    root->next->next->back->v = newl;

    /* Connect root's children */
    hookup(root->next->back, root->next->next->back);

    /* Move nodep entries down one and set indices */
    for ( i = spp; i < nonodes-1; i++ ) {
      p = nodep[i+1];
      nodep[i] = p;
      nodep[i+1] = NULL;
      if ( nodep[i] == NULL ) /* This may happen in a
                                 multifurcating intree */
        break;
      do {
        p->index = i+1;
        p = p->next;
      } while (p != nodep[i]);
    }

    /* Free protx arrays from old root */
    /*
    for(i=0 ; i < endsite ; i++){
      free(root->x[i]);
      free(root->next->x[i]);
      free(root->next->next->x[i]);
      root->x[i] = NULL;
      root->next->x[i] = NULL;
      root->next->next->x[i] = NULL;
    }
    free(root->x);
    free(root->next->x);
    free(root->next->next->x);
    */

    chuck(&grbg,root->next->next);
    chuck(&grbg,root->next);
    chuck(&grbg,root);
  }
} /* dnaml_unroot */

void inittravtree(tree* t,node *p)
{ /* traverse tree to set initialized and v to initial values */
  node* q; 

  if ( p->branchnum == 0) {
    set_branchnum(p, get_trans(t));
    set_branchnum(p->back, p->branchnum);
  }
  p->initialized = false;
  p->back->initialized = false;
  if ( usertree && (!lengths || p->iter)) {
    branchtrans(p->branchnum, initialv);
    p->v = initialv;
    p->back->v = initialv;
  }
  else branchtrans(p->branchnum, p->v);

  if ( !p->tip ) {
    q = p->next;
    while ( q != p ) {
      inittravtree(t,q->back);
      q = q->next;
    }
  }
} /* inittravtree */


double adjusted_v(double v) {
  return 3.0/4.0 * (1.0-exp(-4.0/3.0 * v));
}


static void
adjust_lengths_r(node *p)
{
  node *q; 

  p->v = adjusted_v(p->v);
  p->back->v = p->v;
  if (!p->tip) {
    for (q = p->next; q != p; q = q->next)
      adjust_lengths_r(q->back);
  }
}


void adjust_lengths(tree *t)
{
  assert(t->start->back->tip);
  adjust_lengths_r(t->start);
}


void treevaluate()
{
  /* find maximum likelihood branch lengths of user tree */
  long i;

  if ( lengths)  
    adjust_lengths(&curtree);
  nonodes2--;

  inittravtree(&curtree,curtree.start);
  inittravtree(&curtree,curtree.start->back);
  smoothit = true;
  for (i = 1; i <= smoothings * 4; i++)  {
    smooth (curtree.start);
    smooth (curtree.start->back);
  }
  evaluate(&curtree, curtree.start);
  
  nonodes2++;
}  /* treevaluate */


void maketree()
{
  /* construct and rearrange tree */
  long i,j;
  long nextnode;

  if (usertree) {
    /* Open in binary: ftell() is broken for UNIX line-endings under WIN32 */
    openfile(&intree,INTREE,"input tree file","rb",progname,intreename);
    numtrees = countsemic(&intree);
    if(numtrees > MAXSHIMOTREES)
      shimotrees = MAXSHIMOTREES;
    else
      shimotrees = numtrees;
    if (numtrees > 2)
      initseed(&inseed, &inseed0, seed);
    l0gl = (double *) Malloc(shimotrees * sizeof(double));
    l0gf = (double **) Malloc(shimotrees * sizeof(double *));
    for (i=0; i < shimotrees; ++i)
      l0gf[i] = (double *)Malloc(endsite * sizeof(double));
    if (treeprint) {
      fprintf(outfile, "User-defined tree");
      if (numtrees > 1)
        putc('s', outfile);
      fprintf(outfile, ":\n\n");
    }
    which = 1;
    while (which <= numtrees) {
      /*
      treeread2 (intree, &curtree.start, curtree.nodep, lengths, &trweight, 
                      &goteof, &haslengths, &spp,false,nonodes2);
      */

      /* These initializations required each time through the loop
         since multiple trees require re-initialization */
      nextnode         = 0;
      goteof           = false;

      treeread(intree, &curtree.start, NULL, &goteof, NULL,
               curtree.nodep, &nextnode, NULL, &grbg, 
               initrestmlnode, false, nonodes2);

      restml_unroot(curtree.start, curtree.nodep, nonodes2);

      if ( outgropt )
        curtree.start = curtree.nodep[outgrno - 1]->back;
      else
        curtree.start = curtree.nodep[0]->back;

      treevaluate();
      restml_fprintree(outfile);
      summarize();
      if (trout) {
        col = 0;
        restml_treeout(curtree.start);
      }
      clear_connections(&curtree,nonodes2); 

      which++;
    }
    FClose(intree);
    if (numtrees > 1 && weightsum > 1 )
    standev2(numtrees, maxwhich, 0, endsite-1, maxlogl, l0gl, l0gf,
              aliasweight, seed);
  } else {
    free_all_trans(&curtree);
    for (i = 1; i <= spp; i++)
      enterorder[i - 1] = i;
    if (jumble)
      randumize(seed, enterorder);
    if (progress) {
      printf("\nAdding species:\n");
      writename(0, 3, enterorder);
    }
    nextsp = 3;
    smoothit = improve;
    buildsimpletree(&curtree);
    curtree.start = curtree.nodep[enterorder[0] - 1]->back;
    nextsp = 4;
    while (nextsp <= spp) {
      buildnewtip(enterorder[nextsp - 1], &curtree);
      /* bestyet = - nextsp*sites*sitelength*log(4.0); */
      bestyet = -DBL_MAX;
      if (smoothit)
        restml_copy_(&curtree, &priortree);
      addtraverse(curtree.nodep[enterorder[nextsp - 1] - 1]->back,
                  curtree.nodep[enterorder[0]-1]->back, true);
      if (smoothit)
        restml_copy_(&bestree, &curtree);
      else {
        smoothit = true;
        insert_(curtree.nodep[enterorder[nextsp - 1] - 1]->back, qwhere);
        for (i = 1; i<=smoothings; i++) {
          smooth (curtree.start);
          smooth (curtree.start->back);
        }
        smoothit = false;
        /* bestyet = curtree.likelihood; */
      }
      if (progress)
        writename(nextsp - 1, 1, enterorder);
      if (global && nextsp == spp) {
        if (progress) {
          printf("Doing global rearrangements\n");
          printf("  !");
          for (j = spp ; j < nonodes2 ; j++)
            if ( (j - spp) % (( nonodes2 / 72 ) + 1 ) == 0 )
              putchar('-');
          putchar('!');
        }
      }
      succeeded = true;
      while (succeeded) {
        succeeded = false;
        if (global && nextsp == spp)
          globrearrange();
        else 
          rearrange(curtree.start, curtree.start->back);
      }
      nextsp++;
    }
    if (global && progress) {
      putchar('\n');
      fflush(stdout);
    }
    restml_copy_(&curtree, &bestree);
    if (njumble > 1) {
      if (jumb == 1)
        restml_copy_(&bestree, &bestree2);
      else
        if (bestree2.likelihood < bestree.likelihood)
          restml_copy_(&bestree, &bestree2);
    }
    if (jumb == njumble) {
      if (njumble > 1)
        restml_copy_(&bestree2, &curtree);
      curtree.start = curtree.nodep[outgrno - 1]->back;
      restml_fprintree(outfile);
      summarize();
      if (trout) {
        col = 0;
        restml_treeout(curtree.start);
      }
    }
  }
  if ( jumb < njumble )
    return;
  freex2(nonodes2, curtree.nodep);
   if (!usertree) {
     freex2(nonodes2, priortree.nodep);
     freex2(nonodes2, bestree.nodep);
     if (njumble > 1)
       freex2(nonodes2, bestree2.nodep);
   } else {
     free(l0gl);
     for (i=0;i<shimotrees;++i)
       free(l0gf[i]);
     free(l0gf);
   }
  if (jumb == njumble) {
    if (progress) {
      printf("\nOutput written to file \"%s\"\n\n", outfilename);
      if (trout)
        printf("Tree also written onto file \"%s\"\n", outtreename);
      putchar('\n');
    }
  }
}  /* maketree */


int main(int argc, Char *argv[])
{  /* maximum likelihood phylogenies from restriction sites */
#ifdef MAC
  argc = 1;                /* macsetup("Restml","");        */
  argv[0] = "RestML";
#endif
  init(argc, argv);
  progname = argv[0];
  openfile(&infile,INFILE,"input file","r",argv[0],infilename);
  openfile(&outfile,OUTFILE,"output file","w",argv[0],outfilename);
  ibmpc = IBMCRT;
  ansi = ANSICRT;
  mulsets = false;
  datasets = 1;
  firstset = true;
  doinit();
  if (trout)
    openfile(&outtree,OUTTREE,"output tree file","w",argv[0],outtreename);
  for (ith = 1; ith <= datasets; ith++) {
    if (datasets > 1) {
      fprintf(outfile, "Data set # %ld:\n",ith);
      if (progress)
        printf("\nData set # %ld:\n",ith);
    }
    getinput();
    if (ith == 1)
      firstset = false;
    for (jumb = 1; jumb <= njumble; jumb++)
      maketree();
  }
  cleanup();
  FClose(infile);
  FClose(outfile);
  FClose(outtree);
#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(outtreename);
#endif
  printf("Done.\n\n");
#ifdef WIN32
  phyRestoreConsoleAttributes();
#endif
  return 0;
}  /* maximum likelihood phylogenies from restriction sites */
