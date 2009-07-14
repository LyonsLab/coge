/* PHYLIP version 3.6. (c) Copyright 1986-2007 by the University of Washington
 * and by Joseph Felsenstein.  Written by Joseph Felsenstein.  Permission is
 * granted to copy and use this program provided no fee is charged for it and
 * provided that this copyright notice is not removed. */

#include <float.h>
#include <stdlib.h>

#include "phylip.h"
#include "seq.h"
#include "mlclock.h"
#include "printree.h"

#define over            60       /* Maximum xcoord of tip nodes */

/* These are redefined from phylip.h */
/* Fractional accuracy to which node tymes are optimized */
#undef epsilon
double epsilon = 1e-3;

/* Number of (failed) passes over the tree before giving up */
#undef smoothings
#define smoothings      4

#undef initialv
#define initialv        0.3


typedef struct valrec {
  double rat, ratxi, ratxv, orig_zz, z1, y1, z1zz, z1yy, xiz1,
         xiy1xv;
  double *ww, *zz, *wwzz, *vvzz; 
} valrec;

struct options {
  boolean auto_;
  boolean ctgry;
  long    categs;
  long    rcategs;
  boolean freqsfrom;
  boolean gama;
  boolean invar;
  boolean global;
  boolean hypstate;
  boolean jumble;
  long    njumble;
  double  lambda;
  double  lambda1;
  boolean lengthsopt;
  boolean trout;
  double  ttratio;
  boolean ttr;
  boolean usertree;
  boolean weights;
  boolean printdata;
  boolean dotdiff;
  boolean progress;
  boolean treeprint;
  boolean interleaved;
};


typedef double contribarr[maxcategs];

valrec ***tbl;

#ifndef OLDC
/* function prototypes */
static void    menuconf(void);
static void    allocrest(void);
static void    doinit(void);
static void    inputoptions(void);
static void    makeweights(void);
static void    getinput(void);
static void    inittable_for_usertree (FILE *);
static void    inittable(void);
static void    alloc_nvd(long, nuview_data *);
static void    free_nvd(nuview_data *);
static node   *invalid_descendant_view(node *);
static boolean nuview(node *);
static double  dnamlk_evaluate(node *);
static boolean update(node *);
static boolean smooth(node *);
static void    restoradd(node *, node *, node *, double);
static void    dnamlk_add(node *, node *, node *);
static void    dnamlk_re_move(node **, node **, boolean);
static void    tryadd(node *, node **, node **);
static void    addpreorder(node *, node *, node *, boolean, boolean);
static boolean tryrearr(node *);
static boolean repreorder(node *);
static void    rearrange(node **);
static void    initdnamlnode(node **, node **, node *, long, long, long *, long *,
                initops, pointarray, pointarray, Char *, Char *, FILE *);
static double    tymetrav(node *p);
static void    reconstr(node *, long);
static void    rectrav(node *, long, long);
static void    summarize(FILE *fp);
static void    dnamlk_treeout(node *);
static void    init_tymes(node *p, double minlength);
static void    treevaluate(void);
static void    maketree(void);
static void    reallocsites(void);
static void    freetable(void);
#endif /* OLDC */


Char infilename[FNMLNGTH], outfilename[FNMLNGTH], intreename[FNMLNGTH], outtreename[FNMLNGTH],
     catfilename[FNMLNGTH], weightfilename[FNMLNGTH];
double *rrate;
long sites, weightsum, categs, datasets, ith, njumble, jumb, numtrees, shimotrees;
/*  sites = number of sites in actual sequences
  numtrees = number of user-defined trees */
long inseed, inseed0, mx, mx0, mx1;
boolean freqsfrom, global, global2=0, jumble, trout, usertree, weights, 
          rctgry, ctgry, ttr, auto_, progress, mulsets, firstset, hypstate, 
          smoothit, polishing, justwts, gama, invar;
boolean lengthsopt = false;             /* Use lengths in user tree option */
boolean lngths     = false;             /* Actually use lengths (depends on
                                           each input tree) */
tree curtree, bestree, bestree2;
node *qwhere, *grbg;
double *tymes;
double xi, xv, ttratio, ttratio0, freqa, freqc, freqg, freqt, freqr,
              freqy, freqar, freqcy, freqgr, freqty, fracchange, sumrates,
              cv, alpha, lambda, lambda1, invarfrac;
long *enterorder;
steptr aliasweight;
double *rate;
longer seed;
double *probcat;
contribarr *contribution;
char *progname;
long rcategs;
long **mp;
char basechar[16] = "acmgrsvtwyhkdbn";


/* Local variables for maketree, propagated globally for C version: */
long    k, maxwhich, col;
double  like, bestyet, maxlogl;
boolean lastsp;

boolean smoothed; /* set true before each smoothing run, and set false
                     each time a branch cannot be completely optimized. */
double  *l0gl;
double  expon1i[maxcategs], expon1v[maxcategs],
        expon2i[maxcategs], expon2v[maxcategs];
node   *there;
double  **l0gf;
Char ch, ch2;

static void menuconf()
{
  /* interactively set options */
  long i, loopcount, loopcount2;
  Char ch;
  boolean done;
  boolean didchangecat, didchangercat;
  double probsum;

  fprintf(outfile, "\nNucleic acid sequence\n");
  fprintf(outfile, "   Maximum Likelihood method with molecular ");
  fprintf(outfile, "clock, version %s\n\n", VERSION);
  putchar('\n');
  auto_ = false;
  ctgry = false;
  didchangecat = false;
  rctgry = false;
  didchangercat = false;
  categs = 1;
  rcategs = 1;
  freqsfrom = true;
  gama = false;
  invar = false;
  global = false;
  hypstate = false;
  jumble = false;
  njumble = 1;
  lambda = 1.0;
  lambda1 = 0.0;
  lengthsopt = false;
  trout = true;
  ttratio = 2.0;
  ttr = false;
  usertree = false;
  weights = false;
  printdata = false;
  dotdiff = true;
  progress = true;
  treeprint = true;
  interleaved = true;
  loopcount = 0;
  do {
    cleerhome();
    printf("\nNucleic acid sequence\n");
    printf("   Maximum Likelihood method with molecular clock, version %s\n\n",
           VERSION);
    printf("Settings for this run:\n");
    printf("  U                 Search for best tree?");
    if (usertree)
     printf("  No, use user trees in input file\n");
    else
      printf("  Yes\n");
    if (usertree) {
      printf("  L           Use lengths from user tree?");
      if (lengthsopt)
        printf("  Yes\n");
      else
        printf("  No\n");
    }
    printf("  T        Transition/transversion ratio:");
    if (!ttr)
      printf("  2.0\n");
    else
      printf("  %8.4f\n", ttratio);
    printf("  F       Use empirical base frequencies?");
    if (freqsfrom)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  C   One category of substitution rates?");
    if (!ctgry)
      printf("  Yes\n");
    else
      printf("  %ld categories\n", categs);
    printf("  R           Rate variation among sites?");
    if (!rctgry)
      printf("  constant rate\n");
    else {
      if (gama)
        printf("  Gamma distributed rates\n");
      else {
        if (invar)
          printf("  Gamma+Invariant sites\n");
        else
          printf("  user-defined HMM of rates\n");
      }
      printf("  A   Rates at adjacent sites correlated?");
      if (!auto_)
        printf("  No, they are independent\n");
      else
        printf("  Yes, mean block length =%6.1f\n", 1.0 / lambda);
    }
    if (!usertree) {
      printf("  G                Global rearrangements?");
      if (global)
        printf("  Yes\n");
      else
        printf("  No\n");
    }
    printf("  W                       Sites weighted?  %s\n",
           (weights ? "Yes" : "No"));
    if (!usertree) {
      printf("  J   Randomize input order of sequences?");
      if (jumble)
        printf("  Yes (seed = %8ld, %3ld times)\n", inseed0, njumble);
      else
        printf("  No. Use input order\n");
    }
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld %s\n", datasets,
               (justwts ? "sets of weights" : "data sets"));
    else
      printf("  No\n");
    printf("  I          Input sequences interleaved?");
    if (interleaved)
      printf("  Yes\n");
    else
      printf("  No, sequential\n");
    printf("  0   Terminal type (IBM PC, ANSI, none)?");
    if (ibmpc)
      printf("  IBM PC\n");
    if (ansi)
      printf("  ANSI\n");
    if (!(ibmpc || ansi))
      printf("  (none)\n");
    printf("  1    Print out the data at start of run");
    if (printdata)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  2  Print indications of progress of run");
    if (progress)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  3                        Print out tree");
    if (treeprint)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  4       Write out trees onto tree file?");
    if (trout)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  5   Reconstruct hypothetical sequences?  %s\n",
           (hypstate ? "Yes" : "No"));
    printf("\nAre these settings correct? "
        "(type Y or the letter for one to change)\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    fflush(stdout);
    scanf("%c%*[^\n]", &ch);
    getchar();
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
    done = (ch == 'Y');
    if (!done) {
      uppercase(&ch);
      if (((!usertree) && (strchr("JUCRAFWGTMI012345", ch) != NULL))
          || (usertree && ((strchr("UCRAFWLTMI012345", ch) != NULL)))){
        switch (ch) {

        case 'C':
          ctgry = !ctgry;
          if (ctgry) {
            printf("\nSitewise user-assigned categories:\n\n");
            initcatn(&categs);
            if (rate){
              free(rate);
            }
            rate    = (double *) Malloc( categs * sizeof(double));
            didchangecat = true;
            initcategs(categs, rate);
          }
          break;

        case 'R':
          if (!rctgry) {
            rctgry = true;
            gama = true;
          } else {
            if (gama) {
              gama = false;
              invar = true;
            } else {
              if (invar)
                invar = false;
              else
                rctgry = false;
            }
          }
          break;

        case 'A':
          auto_ = !auto_;
          if (auto_) {
            initlambda(&lambda);
            lambda1 = 1.0 - lambda;
          }
          break;

        case 'F':
          freqsfrom = !freqsfrom;
          if (!freqsfrom)
            initfreqs(&freqa, &freqc, &freqg, &freqt);
          break;

        case 'G':
          global = !global;
          break;

        case 'W':
          weights = !weights;
          break;

        case 'J':
          jumble = !jumble;
          if (jumble)
            initjumble(&inseed, &inseed0, seed, &njumble);
          else njumble = 1;
          break;

        case 'L':
          lengthsopt = !lengthsopt;
          break;

        case 'T':
          ttr = !ttr;
          if (ttr)
            initratio(&ttratio);
          break;

        case 'U':
          usertree = !usertree;
          break;

        case 'M':
          mulsets = !mulsets;
          if (mulsets) {
            printf("Multiple data sets or multiple weights?");
            loopcount2 = 0;
            do {
              printf(" (type D or W)\n");
#ifdef WIN32
              phyFillScreenColor();
#endif
              fflush(stdout);
              scanf("%c%*[^\n]", &ch2);
              getchar();
              if (ch2 == '\n')
                ch2 = ' ';
              uppercase(&ch2);
              countup(&loopcount2, 10);
            } while ((ch2 != 'W') && (ch2 != 'D'));
            justwts = (ch2 == 'W');
            if (justwts)
              justweights(&datasets);
            else
              initdatasets(&datasets);
            if (!jumble) {
              jumble = true;
              initjumble(&inseed, &inseed0, seed, &njumble);
            }
          }
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

        case '5':
          hypstate = !hypstate;
          break;
        }
      } else
        printf("Not a possible option!\n");
    }
    countup(&loopcount, 100);
  } while (!done);
  if (gama || invar) {
    loopcount = 0;
    do {
      printf(
"\nCoefficient of variation of substitution rate among sites (must be positive)\n");
      printf(
  " In gamma distribution parameters, this is 1/(square root of alpha)\n");
#ifdef WIN32
      phyFillScreenColor();
#endif
      fflush(stdout);
      scanf("%lf%*[^\n]", &cv);
      getchar();
      countup(&loopcount, 10);
    } while (cv <= 0.0);
    alpha = 1.0 / (cv * cv);
  }
  if (!rctgry)
    auto_ = false;
  if (rctgry) {
    printf("\nRates in HMM");
    if (invar)
      printf(" (including one for invariant sites)");
    printf(":\n");
    initcatn(&rcategs);
    if (probcat){
      free(probcat);
      free(rrate);
    }
    probcat = (double *) Malloc(rcategs * sizeof(double));
    rrate   = (double *) Malloc(rcategs * sizeof(double));
    didchangercat = true;
    if (gama)
      initgammacat(rcategs, alpha, rrate, probcat); 
    else {
      if (invar) {
        loopcount = 0;
        do {
          printf("Fraction of invariant sites?\n");
          fflush(stdout);
          scanf("%lf%*[^\n]", &invarfrac);
          getchar();
          countup(&loopcount, 10);
        } while ((invarfrac <= 0.0) || (invarfrac >= 1.0));
        initgammacat(rcategs-1, alpha, rrate, probcat); 
        for (i = 0; i < rcategs-1; i++)
          probcat[i] = probcat[i]*(1.0-invarfrac);
        probcat[rcategs-1] = invarfrac;
        rrate[rcategs-1] = 0.0;
      } else {
        initcategs(rcategs, rrate);
        initprobcat(rcategs, &probsum, probcat);
      }
    }
  }
  if (!didchangercat){
    rrate      = (double *)Malloc( rcategs*sizeof(double));
    probcat    = (double *)Malloc( rcategs*sizeof(double));
    rrate[0]   = 1.0;
    probcat[0] = 1.0;
  }
  if (!didchangecat){
    rate       = (double *)Malloc( categs*sizeof(double));
    rate[0]    = 1.0;
  }
}  /* menuconf */

static void reallocsites(void) 
{
  long i;

  for (i = 0; i < spp; i++) {
    free(y[i]);
    y[i] = (char *)Malloc(sites * sizeof(char));
  }
  free(weight);
  free(category);
  free(alias);
  free(aliasweight);
  free(ally);
  free(location);
  
  weight      = (long *)Malloc(sites*sizeof(long));
  category    = (long *)Malloc(sites*sizeof(long));
  alias       = (long *)Malloc(sites*sizeof(long));
  aliasweight = (long *)Malloc(sites*sizeof(long));
  ally        = (long *)Malloc(sites*sizeof(long));
  location    = (long *)Malloc(sites*sizeof(long));
} /* reallocsites */


static void allocrest()
{
  long i;

  y     = (Char **)Malloc(spp*sizeof(Char *));
  nayme  = (naym *)Malloc(spp*sizeof(naym));
  for (i = 0; i < spp; i++)
    y[i] = (char *)Malloc(sites * sizeof(char));
  enterorder  = (long *)Malloc(spp*sizeof(long));
  weight      = (long *)Malloc(sites*sizeof(long));
  category    = (long *)Malloc(sites*sizeof(long));
  alias       = (long *)Malloc(sites*sizeof(long));
  aliasweight = (long *)Malloc(sites*sizeof(long));
  ally        = (long *)Malloc(sites*sizeof(long));
  location    = (long *)Malloc(sites*sizeof(long));
  tymes       = (double *)Malloc((nonodes - spp) * sizeof(double));
}  /* allocrest */


static void doinit()
{
  /* initializes variables */

  inputnumbers(&spp, &sites, &nonodes, 1);
  menuconf();
  if (printdata)
    fprintf(outfile, "%2ld species, %3ld  sites\n", spp, sites);
  alloctree(&curtree.nodep, nonodes, usertree);
  allocrest();
  if (usertree)
    return;
  alloctree(&bestree.nodep, nonodes, 0);
  if (njumble <= 1)
    return;
  alloctree(&bestree2.nodep, nonodes, 0);
}  /* doinit */


static void inputoptions()
{
  long i;

  if (!firstset && !justwts) {
    samenumsp(&sites, ith);
    reallocsites();
  }

  for (i = 0; i < sites; i++)
    category[i] = 1;
  for (i = 0; i < sites; i++)
    weight[i] = 1;
  
  if (justwts || weights)
    inputweights(sites, weight, &weights);
  weightsum = 0;
  for (i = 0; i < sites; i++)
    weightsum += weight[i];
  if (ctgry && categs > 1) {
    inputcategs(0, sites, category, categs, "DnaMLK");
    if (printdata)
      printcategs(outfile, sites, category, "Site categories");
  }
  if (weights && printdata)
    printweights(outfile, 0, sites, weight, "Sites");
}  /* inputoptions */


static void makeweights()
{
  /* make up weights vector to avoid duplicate computations */
  long i;

   for (i = 1; i <= sites; i++) {
    alias[i - 1] = i;
    ally[i - 1] = 0;
    aliasweight[i - 1] = weight[i - 1];
    location[i - 1] = 0;
  }
  sitesort2(sites, aliasweight);
  sitecombine2(sites, aliasweight);
  sitescrunch2(sites, 1, 2, aliasweight);
  for (i = 1; i <= sites; i++) {
    if (aliasweight[i - 1] > 0)
      endsite = i;
  }
  for (i = 1; i <= endsite; i++) {
    ally[alias[i - 1] - 1] = alias[i - 1];
    location[alias[i - 1] - 1] = i;
  }
  contribution = (contribarr *) Malloc( endsite*sizeof(contribarr));
}  /* makeweights */


static void getinput()
{
  /* reads the input data */
  inputoptions();
  if (!freqsfrom)
    getbasefreqs(freqa, freqc, freqg, freqt, &freqr, &freqy, &freqar, &freqcy,
                   &freqgr, &freqty, &ttratio, &xi, &xv, &fracchange,
                   freqsfrom, true);
  if (!justwts || firstset)
    inputdata(sites);
  makeweights();
  setuptree2(&curtree);
  if (!usertree) {
    setuptree2(&bestree);
    if (njumble > 1)
      setuptree2(&bestree2);
  }
  allocx(nonodes, rcategs, curtree.nodep, usertree);
  if (!usertree) {
    allocx(nonodes, rcategs, bestree.nodep, 0);
    if (njumble > 1)
      allocx(nonodes, rcategs, bestree2.nodep, 0);
  }
  makevalues2(rcategs, curtree.nodep, endsite, spp, y, alias);
  if (freqsfrom) {
    empiricalfreqs(&freqa, &freqc, &freqg, &freqt, aliasweight, curtree.nodep);
    getbasefreqs(freqa, freqc, freqg, freqt, &freqr, &freqy, &freqar, &freqcy,
                   &freqgr, &freqty, &ttratio, &xi, &xv, &fracchange,
                   freqsfrom, true);
  }
  if (!justwts || firstset)
    fprintf(outfile, "\nTransition/transversion ratio = %10.6f\n\n", ttratio);
}  /* getinput */


static void inittable_for_usertree (FILE *intree)
{
  /* If there's a user tree, then the ww/zz/wwzz/vvzz elements need
     to be allocated appropriately. */
  long num_comma;
  long i, j;

  /* First, figure out the largest possible furcation, i.e. the number
     of commas plus one */
  countcomma (&intree, &num_comma);
  num_comma++;
  
  for (i = 0; i < rcategs; i++) {
    for (j = 0; j < categs; j++) {
      /* Free the stuff allocated assuming bifurcations */
      free (tbl[i][j]->ww);
      free (tbl[i][j]->zz);
      free (tbl[i][j]->wwzz);
      free (tbl[i][j]->vvzz);

      /* Then allocate for worst-case multifurcations */
      tbl[i][j]->ww   = (double *) Malloc( num_comma * sizeof (double));
      tbl[i][j]->zz   = (double *) Malloc( num_comma * sizeof (double));
      tbl[i][j]->wwzz = (double *) Malloc( num_comma * sizeof (double));
      tbl[i][j]->vvzz = (double *) Malloc( num_comma * sizeof (double));
    }
  }
}  /* inittable_for_usertree */


static void freetable()
{
  long i, j;
 
  for (i = 0; i < rcategs; i++) {
    for (j = 0; j < categs; j++) {
      free(tbl[i][j]->ww);
      free(tbl[i][j]->zz);
      free(tbl[i][j]->wwzz);
      free(tbl[i][j]->vvzz);
    }
  }
  for (i = 0; i < rcategs; i++)  {
    for (j = 0; j < categs; j++) 
      free(tbl[i][j]);
    free(tbl[i]);
  }
  free(tbl);
}


static void inittable()
{
  /* Define a lookup table. Precompute values and print them out in tables */
  long i, j;
  double sumrates;
  
  tbl = (valrec ***) Malloc( rcategs * sizeof(valrec **));
  for (i = 0; i < rcategs; i++) {
    tbl[i] = (valrec **) Malloc( categs*sizeof(valrec *));
    for (j = 0; j < categs; j++)
      tbl[i][j] = (valrec *) Malloc( sizeof(valrec));
  }

  for (i = 0; i < rcategs; i++) {
    for (j = 0; j < categs; j++) {
      tbl[i][j]->rat = rrate[i]*rate[j];
      tbl[i][j]->ratxi = tbl[i][j]->rat * xi;
      tbl[i][j]->ratxv = tbl[i][j]->rat * xv;

      /* Allocate assuming bifurcations, will be changed later if
         neccesarry (i.e. there's a user tree) */
      tbl[i][j]->ww   = (double *) Malloc( 2 * sizeof (double));
      tbl[i][j]->zz   = (double *) Malloc( 2 * sizeof (double));
      tbl[i][j]->wwzz = (double *) Malloc( 2 * sizeof (double));
      tbl[i][j]->vvzz = (double *) Malloc( 2 * sizeof (double));
    }
  }
  sumrates = 0.0;
  for (i = 0; i < endsite; i++) {
    for (j = 0; j < rcategs; j++)
      sumrates += aliasweight[i] * probcat[j] 
        * tbl[j][category[alias[i] - 1] - 1]->rat;
  }
  sumrates /= (double)sites;
  for (i = 0; i < rcategs; i++)
    for (j = 0; j < categs; j++) {
      tbl[i][j]->rat /= sumrates;
      tbl[i][j]->ratxi /= sumrates;
      tbl[i][j]->ratxv /= sumrates;
    }

  if(jumb > 1)
    return;

  if (gama || invar) {
    fprintf(outfile, "\nDiscrete approximation to gamma distributed rates\n");
    fprintf(outfile,
    " Coefficient of variation of rates = %f  (alpha = %f)\n", cv, alpha);
  }
  if (rcategs > 1) {
    fprintf(outfile, "\nState in HMM    Rate of change    Probability\n\n");
    for (i = 0; i < rcategs; i++)
      if (probcat[i] < 0.0001)
        fprintf(outfile, "%9ld%16.3f%20.6f\n", i+1, rrate[i], probcat[i]);
      else if (probcat[i] < 0.001)
          fprintf(outfile, "%9ld%16.3f%19.5f\n", i+1, rrate[i], probcat[i]);
        else if (probcat[i] < 0.01)
            fprintf(outfile, "%9ld%16.3f%18.4f\n", i+1, rrate[i], probcat[i]);
          else
            fprintf(outfile, "%9ld%16.3f%17.3f\n", i+1, rrate[i], probcat[i]);
    putc('\n', outfile);
    if (auto_) {
      fprintf(outfile,
     "Expected length of a patch of sites having the same rate = %8.3f\n",
             1/lambda);
      putc('\n', outfile);
    }
  }
  if (categs > 1) {
    fprintf(outfile, "\nSite category   Rate of change\n\n");
    for (i = 0; i < categs; i++)
      fprintf(outfile, "%9ld%16.3f\n", i+1, rate[i]);
    fprintf(outfile, "\n\n");
  }
}  /* inittable */


static void alloc_nvd(long num_sibs, nuview_data *local_nvd)
{
  /* Allocate blocks of memory appropriate for the number of siblings
     a given node has */
  local_nvd->yy     = (double *) Malloc( num_sibs * sizeof (double));
  local_nvd->wwzz   = (double *) Malloc( num_sibs * sizeof (double));
  local_nvd->vvzz   = (double *) Malloc( num_sibs * sizeof (double));
  local_nvd->vzsumr = (double *) Malloc( num_sibs * sizeof (double));
  local_nvd->vzsumy = (double *) Malloc( num_sibs * sizeof (double));
  local_nvd->sum    = (double *) Malloc( num_sibs * sizeof (double));
  local_nvd->sumr   = (double *) Malloc( num_sibs * sizeof (double));
  local_nvd->sumy   = (double *) Malloc( num_sibs * sizeof (double));
  local_nvd->xx     = (sitelike *) Malloc( num_sibs * sizeof (sitelike));
}  /* alloc_nvd */


static void free_nvd(nuview_data *local_nvd)
{
  /* The natural complement to the alloc version */
  free (local_nvd->yy);
  free (local_nvd->wwzz);
  free (local_nvd->vvzz);
  free (local_nvd->vzsumr);
  free (local_nvd->vzsumy);
  free (local_nvd->sum);
  free (local_nvd->sumr);
  free (local_nvd->sumy);
  free (local_nvd->xx);
}  /* free_nvd */


static node *invalid_descendant_view(node *p)
{  
  /* Useful as an assertion - Traverses all descendants of p looking for
   * uninitialized views.  Returns NULL if none exist, otherwise a pointer to
   * the first one found. */

  node *q = NULL;
  node *s = NULL;

  if ( p == NULL || p->tip )
    return NULL;

  for ( q = p->next; q != p; q = q->next ) {
    s = invalid_descendant_view(q->back);
    if ( s != NULL )
      return s;
  }

  return s;
} /* invalid_descendant_view */


static boolean nuview(node *p)
{
  /* Recursively update summary data for subtree rooted at p. Returns true if
   * view has changed. */

  long i, j, k, l, num_sibs = 0, sib_index;
  nuview_data *local_nvd;
  node *q;
  node *sib_ptr, *sib_back_ptr;
  sitelike p_xx;
  double lw;
  double correction;
  double maxx;

  assert(p != NULL);
  if ( p == NULL ) return false;
  if ( p->tip ) return false; /* Tips do not need to be initialized */

  for (q = p->next; q != p; q = q->next) {
    num_sibs++;
    if ( q->back != NULL && !q->tip) {
      if ( nuview(q->back) )
        p->initialized = false;
    }
  }

  if ( p->initialized )
    return false;

  /* At this point, all views downstream should be initialized.
   * If not, we have a problem. */
  assert( invalid_descendant_view(p) == NULL );

  /* Allocate the structure and blocks therein for variables used in
     this function */
  local_nvd = (nuview_data *) Malloc( sizeof (nuview_data));
  alloc_nvd (num_sibs, local_nvd);

  /* Loop 1: makes assignments to tbl based on some combination of
     what's already in tbl and the children's value of v */
  sib_ptr = p;
  for (sib_index=0; sib_index < num_sibs; sib_index++) {
    sib_ptr      = sib_ptr->next;
    sib_back_ptr = sib_ptr->back;
    
    if (sib_back_ptr != NULL)
      lw = -fabs(p->tyme - sib_back_ptr->tyme);
    else
      lw = 0.0;

    for (i = 0; i < rcategs; i++)
      for (j = 0; j < categs; j++) {
        tbl[i][j]->ww[sib_index]   = exp(tbl[i][j]->ratxi * lw);
        tbl[i][j]->zz[sib_index]   = exp(tbl[i][j]->ratxv * lw);
        tbl[i][j]->wwzz[sib_index] = tbl[i][j]->ww[sib_index] * tbl[i][j]->zz[sib_index];
        tbl[i][j]->vvzz[sib_index] = (1.0 - tbl[i][j]->ww[sib_index]) *
                                                                tbl[i][j]->zz[sib_index];
      }
  }

  /* Loop 2: */
  for (i = 0; i < endsite; i++) {
    correction = 0; 
    maxx = 0;
    k = category[alias[i]-1] - 1;
    for (j = 0; j < rcategs; j++) {

      /* Loop 2.1 */
      sib_ptr = p;
      for (sib_index=0; sib_index < num_sibs; sib_index++) {
        sib_ptr         = sib_ptr->next;
        sib_back_ptr    = sib_ptr->back;
        

        local_nvd->wwzz[sib_index] = tbl[j][k]->wwzz[sib_index];
        local_nvd->vvzz[sib_index] = tbl[j][k]->vvzz[sib_index];
        local_nvd->yy[sib_index]   = 1.0 - tbl[j][k]->zz[sib_index];
        if (sib_back_ptr != NULL) {
          memcpy(local_nvd->xx[sib_index],
               sib_back_ptr->x[i][j],
               sizeof(sitelike));
          if ( j == 0)
            correction += sib_back_ptr->underflows[i];
        }
        else {
          local_nvd->xx[sib_index][0] = 1.0;
          local_nvd->xx[sib_index][(long)C - (long)A] = 1.0;
          local_nvd->xx[sib_index][(long)G - (long)A] = 1.0;
          local_nvd->xx[sib_index][(long)T - (long)A] = 1.0;
        }
      }

      /* Loop 2.2 */
      for (sib_index=0; sib_index < num_sibs; sib_index++) {
        local_nvd->sum[sib_index] =
          local_nvd->yy[sib_index] *
          (freqa * local_nvd->xx[sib_index][(long)A] +
           freqc * local_nvd->xx[sib_index][(long)C] +
           freqg * local_nvd->xx[sib_index][(long)G] +
           freqt * local_nvd->xx[sib_index][(long)T]);
        local_nvd->sumr[sib_index] =
          freqar * local_nvd->xx[sib_index][(long)A] +
          freqgr * local_nvd->xx[sib_index][(long)G];
        local_nvd->sumy[sib_index] =
          freqcy * local_nvd->xx[sib_index][(long)C] +
          freqty * local_nvd->xx[sib_index][(long)T];
        local_nvd->vzsumr[sib_index] =
          local_nvd->vvzz[sib_index] * local_nvd->sumr[sib_index];
        local_nvd->vzsumy[sib_index] =
          local_nvd->vvzz[sib_index] * local_nvd->sumy[sib_index];
      }

      /* Initialize to one, multiply incremental values for every
         sibling a node has */
      p_xx[(long)A] = 1 ;
      p_xx[(long)C] = 1 ; 
      p_xx[(long)G] = 1 ;
      p_xx[(long)T] = 1 ;

      for (sib_index=0; sib_index < num_sibs; sib_index++) {
        p_xx[(long)A] *=
          local_nvd->sum[sib_index] +
          local_nvd->wwzz[sib_index] *
          local_nvd->xx[sib_index][(long)A] +
          local_nvd->vzsumr[sib_index];
        p_xx[(long)C] *=
          local_nvd->sum[sib_index] +
          local_nvd->wwzz[sib_index] *
          local_nvd->xx[sib_index][(long)C] +
          local_nvd->vzsumy[sib_index];
        p_xx[(long)G] *=
          local_nvd->sum[sib_index] +
          local_nvd->wwzz[sib_index] *
          local_nvd->xx[sib_index][(long)G] +
          local_nvd->vzsumr[sib_index];
        p_xx[(long)T] *=
          local_nvd->sum[sib_index] +
          local_nvd->wwzz[sib_index] *
          local_nvd->xx[sib_index][(long)T] +
          local_nvd->vzsumy[sib_index];
      }

      for ( l = 0 ; l < ((long)T - (long)A + 1 ) ; l++ ) {
        if (  p_xx[l] > maxx )
          maxx = p_xx[l];
      }

      /* And the final point of this whole function: */
      memcpy(p->x[i][j], p_xx, sizeof(sitelike));
    }
    p->underflows[i] = 0;
    if ( maxx < MIN_DOUBLE)
      fix_x(p, i, maxx,rcategs);
    p->underflows[i] += correction;
  }

  free_nvd (local_nvd);
  free (local_nvd);

  p->initialized = true;
  return true;
}  /* nuview */


static double dnamlk_evaluate(node *p)
{ /* Evaluate and return the log likelihood of the current tree
   * as seen from the branch from p to p->back. If p is the root node,
   * the first child branch is used instead. Views are updated as needed. */
  contribarr tterm;
  static contribarr like, nulike, clai;
  double sum, sum2, sumc=0, y, lz, y1, z1zz, z1yy, 
                prod12, prod1, prod2, prod3, sumterm, lterm;
  long i, j, k, lai;
  node *q, *r;
  double *x1, *x2;              /* pointers to sitelike elements in node->x */
  sum = 0.0;

  assert( all_tymes_valid(curtree.root, 0.0, false) );

  /* Root node has no branch, so use branch to first child */
  if (p == curtree.root)
    p = p->next;

  r = p;
  q = p->back;
  nuview(r);
  nuview(q);
  y = fabs(r->tyme - q->tyme);
  lz = -y;

  for (i = 0; i < rcategs; i++)
    for (j = 0; j < categs; j++) {
    tbl[i][j]->orig_zz = exp(tbl[i][j]->ratxi * lz);
    tbl[i][j]->z1 = exp(tbl[i][j]->ratxv * lz);
    tbl[i][j]->z1zz = tbl[i][j]->z1 * tbl[i][j]->orig_zz;
    tbl[i][j]->z1yy = tbl[i][j]->z1 - tbl[i][j]->z1zz;
  }
  for (i = 0; i < endsite; i++) {
    k = category[alias[i]-1] - 1;
    for (j = 0; j < rcategs; j++) {
      if (y > 0.0) {
        y1 = 1.0 - tbl[j][k]->z1;
        z1zz = tbl[j][k]->z1zz;
        z1yy = tbl[j][k]->z1yy;
      } else {
        y1 = 0.0;
        z1zz = 1.0;
        z1yy = 0.0;
      }
      x1 = r->x[i][j];
      prod1 = freqa * x1[0] + freqc * x1[(long)C - (long)A] +
             freqg * x1[(long)G - (long)A] + freqt * x1[(long)T - (long)A];
      x2 = q->x[i][j];
      prod2 = freqa * x2[0] + freqc * x2[(long)C - (long)A] +
             freqg * x2[(long)G - (long)A] + freqt * x2[(long)T - (long)A];
      prod3 = (x1[0] * freqa + x1[(long)G - (long)A] * freqg) *
              (x2[0] * freqar + x2[(long)G - (long)A] * freqgr) +
          (x1[(long)C - (long)A] * freqc + x1[(long)T - (long)A] * freqt) *
         (x2[(long)C - (long)A] * freqcy + x2[(long)T - (long)A] * freqty);
      prod12 = freqa * x1[0] * x2[0] +
               freqc * x1[(long)C - (long)A] * x2[(long)C - (long)A] +
               freqg * x1[(long)G - (long)A] * x2[(long)G - (long)A] +
               freqt * x1[(long)T - (long)A] * x2[(long)T - (long)A];
      tterm[j] = z1zz * prod12 + z1yy * prod3 + y1 * prod1 * prod2;
    }
    sumterm = 0.0;
    for (j = 0; j < rcategs; j++)
      sumterm += probcat[j] * tterm[j];
    lterm = log(sumterm) + p->underflows[i] + q->underflows[i];
    for (j = 0; j < rcategs; j++)
      clai[j] = tterm[j] / sumterm;
    memcpy(contribution[i], clai, sizeof(contribarr));
    if (!auto_ && usertree && (which <= shimotrees))
      l0gf[which - 1][i] = lterm;
    sum += aliasweight[i] * lterm;
  }
  if (auto_) {
    for (j = 0; j < rcategs; j++)
      like[j] = 1.0;
    for (i = 0; i < sites; i++) {
      sumc = 0.0;
      for (k = 0; k < rcategs; k++)
        sumc += probcat[k] * like[k];
      sumc *= lambda;
      if ((ally[i] > 0) && (location[ally[i]-1] > 0)) {
        lai = location[ally[i] - 1];
        memcpy(clai, contribution[lai - 1], sizeof(contribarr));
        for (j = 0; j < rcategs; j++)
          nulike[j] = ((1.0 - lambda) * like[j] + sumc) * clai[j];
      } else {
        for (j = 0; j < rcategs; j++)
          nulike[j] = ((1.0 - lambda) * like[j] + sumc);
      }
      memcpy(like, nulike, sizeof(contribarr));
    }
    sum2 = 0.0;
    for (i = 0; i < rcategs; i++)
      sum2 += probcat[i] * like[i];
    sum += log(sum2);
  }
  curtree.likelihood = sum;
  if (auto_ || !usertree)
    return sum;
  if(which <= shimotrees)
    l0gl[which - 1] = sum;
  if (which == 1) {
    maxwhich = 1;
    maxlogl = sum;
    return sum;
  }
  if (sum > maxlogl) {
    maxwhich = which;
    maxlogl = sum;
  }
  return sum;
}  /* dnamlk_evaluate */


static boolean update(node *p)
{
  /* Conditionally optimize tyme at a node. Return true if successful. */

  if (p == NULL)
    return false;

  if ( (!usertree) || (usertree && !lngths) )
    return makenewv(p);
  return false;
}  /* update */


static boolean smooth(node *p)
{
  node *q = NULL;
  boolean success;

  if (p == NULL)
    return false;
  if (p->tip)
    return false;

  /* optimize tyme here */
  success = update(p);

  if (smoothit || polishing) {
    for (q = p->next; q != p; q = q->next) {
      /* smooth subtrees */
      success = smooth(q->back) || success;
      /* optimize tyme again after each subtree */
      success = update(p) || success;  
    }
  }

  return success;
}  /* smooth */


static void restoradd(node *below, node *newtip, node *newfork, double prevtyme)
{
  /* restore "new" tip and fork to place "below".  restore tymes */
  /* assumes bifurcation */

  hookup(newfork, below->back);
  hookup(newfork->next, below);
  hookup(newtip, newfork->next->next);
  curtree.nodep[newfork->index-1] = newfork;
  setnodetymes(newfork,prevtyme);
} /* restoradd */


static void dnamlk_add(node *below, node *newtip, node *newfork)
{
  /* inserts the nodes newfork and its descendant, newtip, into the tree. */
  long i;
  node *p;
  node *above;
  double newtyme;
  boolean success;

  assert( all_tymes_valid(curtree.root, 0.98*MIN_BRANCH_LENGTH, false) );
  assert( floating_fork(newfork) );
  assert( newtip->back == NULL );

  /* Get parent nodelets */
  below = pnode(&curtree, below);
  newfork = pnode(&curtree, newfork);
  newtip = pnode(&curtree, newtip);

  /* Join above node to newfork */
  if (below->back == NULL)
    newfork->back = NULL;
  else {
    above = below->back;
    /* unhookup(below, above); */
    hookup(newfork, above);
  }

  /* Join below to newfork->next->next */
  hookup(below, newfork->next->next);
  /* Join newtip to newfork->next */
  hookup(newfork->next, newtip);

  /* Move root if inserting there */
  if (curtree.root == below)
    curtree.root = newfork;

  /* p = child with lowest tyme */
  p = newtip->tyme < below->tyme ? newtip : below;

  /* If not at root, set newfork tyme to average below/above */
  if (newfork->back != NULL) {
    if (p->tyme > newfork->back->tyme)
      newtyme = (p->tyme + newfork->back->tyme) / 2.0;
    else
      newtyme = p->tyme - INSERT_MIN_TYME;

    if (p->tyme - newtyme < MIN_BRANCH_LENGTH)
      newtyme = p->tyme - MIN_BRANCH_LENGTH;
    setnodetymes(newfork, newtyme);
    /* Now move from newfork to root, setting parent tymes older than children 
     * by at least MIN_BRANCH_LENGTH */
    p = newfork;
    while (p != curtree.root) {
      if (p->back->tyme > p->tyme - MIN_BRANCH_LENGTH)
        setnodetymes(p->back, p->tyme - MIN_BRANCH_LENGTH);
      else
        break;
      /* get parent node */
      p = pnode(&curtree, p->back);
    }
  } else { /* root == newfork */
    /* make root 2x older */
    setnodetymes(newfork, p->tyme - 2*INSERT_MIN_TYME);
  }

  assert( all_tymes_valid(curtree.root, 0.98*MIN_BRANCH_LENGTH, false) );
  
  /* Adjust branch lengths throughout */
  for ( i = 1; i < smoothings; i++ ) {
    success = smooth(newfork);
    success = smooth(newfork->back) || success;
    if ( !success ) break;
  }
}  /* dnamlk_add */


static void dnamlk_re_move(node **item, node **fork, boolean tempadd)
{
  /* removes nodes *item and its parent (returned in *fork), from the tree.
    the new descendant of fork's ancestor is made to be fork's descendant other
    than item. Item must point to node*, but *fork is not read */
  node *p, *q;
  long i;
  boolean success;

  if ((*item)->back == NULL) {
    *fork = NULL;
    return;
  }
  *item = curtree.nodep[(*item)->index-1];
  *fork = curtree.nodep[(*item)->back->index - 1];
  if (curtree.root == *fork) {
    if (*item == (*fork)->next->back)
      curtree.root = (*fork)->next->next->back;
    else
      curtree.root = (*fork)->next->back;
  }
  p = (*item)->back->next->back;
  q = (*item)->back->next->next->back;
  if (p != NULL)
    p->back = q;
  if (q != NULL)
    q->back = p;
  (*fork)->back = NULL;
  p = (*fork)->next;
  while (p != *fork) {
    p->back = NULL;
    p = p->next;
  }
  (*item)->back = NULL;
  inittrav(p);
  inittrav(q);
  if (tempadd)
    return;

  for ( i = 1; i <= smoothings; i++ ) {
    success = smooth(q);
    if ( smoothit )
      success = smooth(q->back) || success;
    if ( !success ) break;
  }
}  /* dnamlk_re_move */


static void tryadd(node *p, node **item, node **nufork)
{  /* temporarily adds one fork and one tip to the tree.
    if the location where they are added yields greater
    likelihood than other locations tested up to that
    time, then keeps that location as there */

  if ( !global2 ) 
      save_tymes(&curtree,tymes);
  dnamlk_add(p, *item, *nufork);
  like = dnamlk_evaluate(p);
  if (lastsp) {
      if (like >= bestree.likelihood || bestree.likelihood == UNDEFINED)  {
        copy_(&curtree, &bestree, nonodes, rcategs);
        if ( global2 ) 
          /* To be restored in maketree() */
            save_tymes(&curtree,tymes);
      }
  }
  if (like > bestyet || bestyet == UNDEFINED) {
    bestyet = like;
    there = p;
  }
  dnamlk_re_move(item, nufork, true);
  if ( !global2 ) {
      restore_tymes(&curtree,tymes);
  }
}  /* tryadd */


static void addpreorder(node *p, node *item_, node *nufork_, boolean contin,
                        boolean continagain)
{
  /* Traverse tree, adding item at different locations until we
   * find a better likelihood. Afterwards, global 'there' will be
   * set to the best add location, or will be left alone if no
   * better could be found. */

  node *item, *nufork;

  item = item_;
  nufork = nufork_;
  if (p == NULL)
    return;
  tryadd(p, &item, &nufork);
  contin = continagain;
  if ((!p->tip) && contin) {
    /* assumes bifurcation (OK) */
    addpreorder(p->next->back, item, nufork, contin, continagain);
    addpreorder(p->next->next->back, item, nufork, contin, continagain);
  }
}  /* addpreorder */


static boolean tryrearr(node *p)
{
  /* evaluates one rearrangement of the tree.
    if the new tree has greater likelihood than the old
    keeps the new tree and returns true.
    otherwise, restores the old tree and returns false. */

  node *forknode;       /* parent fork of p */
  node *frombelow;      /* other child of forknode */
  node *whereto;        /* parent fork of forknode */

  double oldlike;       /* likelihood before rearrangement */
  double prevtyme;      /* forknode->tyme before rearrange */
  double like_delta;    /* improvement in likelihood */
  boolean wasonleft;    /* true if p first child of forknode */

  if (p == curtree.root)
    return false;
  /* forknode = parent fork of p */
  forknode = curtree.nodep[p->back->index - 1];
  if (forknode == curtree.root)
    return false;
  oldlike = bestyet;
  prevtyme = forknode->tyme;
  /* assumes bifurcation (OK) */
  /* frombelow = other child of forknode (not p) */
  if (forknode->next->back == p) {
    frombelow = forknode->next->next->back;
    wasonleft = true;
  }
  else {
    frombelow = forknode->next->back;
    wasonleft = false;
  }
  whereto = curtree.nodep[forknode->back->index - 1];

  /* remove forknode and p */
  dnamlk_re_move(&p, &forknode, true);

  /* add p and forknode as parent of whereto */
  dnamlk_add(whereto, p, forknode);

  like = dnamlk_evaluate(p);
  like_delta = like - oldlike;
  if ( like_delta < LIKE_EPSILON && oldlike != UNDEFINED) {
    dnamlk_re_move(&p, &forknode, true);
    restoradd(frombelow, p, forknode, prevtyme);
    if (wasonleft && (forknode->next->next->back == p)) {
       hookup (forknode->next->back, forknode->next->next);
       hookup (forknode->next, p);
    }
    curtree.likelihood = oldlike;

    /* assumes bifurcation (OK) */
    inittrav(forknode);
    inittrav(forknode->next);
    inittrav(forknode->next->next);

    return false;
  } else {
    bestyet = like;
  }
  return true;
}  /* tryrearr */


static boolean repreorder(node *p)
{
  /* traverses a binary tree, calling function tryrearr
    at a node before calling tryrearr at its descendants.
    Returns true the first time rearrangement increases the
    tree's likelihood. */

  if (p == NULL)
    return false;
  if ( !tryrearr(p) )
    return false;
  if (p->tip)
    return true;
  /* assumes bifurcation */
  if ( !repreorder(p->next->back) )
    return false;
  if ( !repreorder(p->next->next->back) )
    return false;
  
  return true;
}  /* repreorder */


static void rearrange(node **r)
{
  /* traverses the tree (preorder), finding any local
    rearrangement which increases the likelihood.
    if traversal succeeds in increasing the tree's
    likelihood, function rearrange runs traversal again */
  while ( repreorder(*r) )
    /* continue */;

}  /* rearrange */


static void initdnamlnode(node **p, node **grbg, node *q, long len, long nodei,
                     long *ntips, long *parens, initops whichinit,
                     pointarray treenode, pointarray nodep, Char *str, Char *ch,
                     FILE *intree)
{
  /* Initializes each node as it is read from user tree by treeread().
   * whichinit specifies which type of initialization is to be done.
   */
  boolean minusread;
  double valyew, divisor;

  switch (whichinit) {
  case bottom:
    gnu(grbg, p);
    (*p)->index = nodei;
    (*p)->tip = false;
    malloc_pheno((*p), endsite, rcategs);
    nodep[(*p)->index - 1] = (*p);
    break;
  case nonbottom:
    gnu(grbg, p);
    malloc_pheno(*p, endsite, rcategs);
    (*p)->index = nodei;
    break;
  case tip:
    match_names_to_data (str, nodep, p, spp);
    break;
  case iter:
    (*p)->initialized = false;
    /* Initial branch lengths start at 0.0. tymetrav() enforces
     * MIN_BRANCH_LENGTH */
    (*p)->v = 0.0;
    (*p)->iter = true;
    if ((*p)->back != NULL)
      (*p)->back->iter = true;
    break;
  case length:
    processlength(&valyew, &divisor, ch, &minusread, intree, parens);
    (*p)->v = valyew / divisor / fracchange;
    (*p)->iter = false;
    if ((*p)->back != NULL) {
      (*p)->back->v = (*p)->v;
      (*p)->back->iter = false;
    }
    break;
  case hslength:
    break;
  case hsnolength:
    if (usertree && lengthsopt && lngths) {
      printf("Warning: one or more lengths not defined in user tree number %ld.\n", which);
      printf("DNAMLK will attempt to optimize all branch lengths.\n\n");
      lngths = false;
    }
    break;
  case treewt:
    break;
  case unittrwt:
    break;
  }
} /* initdnamlnode */


static double tymetrav(node *p)
{
  /* Recursively convert branch lengths to node tymes. Returns the maximum
   * branch length p's parent can have, which is p->tyme - max(p->v,
   * MIN_BRANCH_LENGTH) */

  node *q;
  double xmax;
  double x;

  xmax = 0.0;
  if (!p->tip) {
    for (q = p->next; q != p; q = q->next) {
      x = tymetrav(q->back);
      if (xmax > x)
        xmax = x;
    }
  } else {
    x = 0.0;
  }

  setnodetymes(p,xmax);

  if (p->v < MIN_BRANCH_LENGTH)
    return xmax - MIN_BRANCH_LENGTH;
  else 
    return xmax - p->v;
}  /* tymetrav */


static void reconstr(node *p, long n)
{
  /* reconstruct and print out base at site n+1 at node p */
  long i, j, k, m, first, second, num_sibs;
  double f, sum, xx[4];
  node *q;

  if ((ally[n] == 0) || (location[ally[n]-1] == 0))
    putc('.', outfile);
  else {
    j = location[ally[n]-1] - 1;
    for (i = 0; i < 4; i++) {
      f = p->x[j][mx-1][i];
      num_sibs = count_sibs(p);
      q = p;
      for (k = 0; k < num_sibs; k++) {
        q = q->next;
        f *= q->x[j][mx-1][i];
      }
      f = sqrt(f);
      xx[i] = f;
    }
    xx[0] *= freqa;
    xx[1] *= freqc;
    xx[2] *= freqg;
    xx[3] *= freqt;
    sum = xx[0]+xx[1]+xx[2]+xx[3];
    for (i = 0; i < 4; i++)
      xx[i] /= sum;
    first = 0;
    for (i = 1; i < 4; i++)
      if (xx [i] > xx[first])
        first = i;
    if (first == 0)
      second = 1;
    else
      second = 0;
    for (i = 0; i < 4; i++)
      if ((i != first) && (xx[i] > xx[second]))
        second = i;
    m = 1 << first;
    if (xx[first] < 0.4999995)
      m = m + (1 << second);
    if (xx[first] > 0.95)
      putc(toupper(basechar[m - 1]), outfile);
    else
      putc(basechar[m - 1], outfile);
    if (rctgry && rcategs > 1)
      mx = mp[n][mx - 1];    
    else
      mx = 1;
  }
} /* reconstr */


static void rectrav(node *p, long m, long n)
{
  /* print out segment of reconstructed sequence for one branch */
  long num_sibs, i;
  node *sib_ptr;

  putc(' ', outfile);
  if (p->tip) {
    for (i = 0; i < nmlngth; i++)
      putc(nayme[p->index-1][i], outfile);
  } else
    fprintf(outfile, "%4ld      ", p->index - spp);
  fprintf(outfile, "  ");
  mx = mx0;
  for (i = m; i <= n; i++) {
    if ((i % 10 == 0) && (i != m))
      putc(' ', outfile);
    if (p->tip)
      putc(y[p->index-1][i], outfile);
    else
      reconstr(p, i);
  }
  putc('\n', outfile);
  if (!p->tip) {
    num_sibs = count_sibs(p);
    sib_ptr = p;
    for (i = 0; i < num_sibs; i++) {
      sib_ptr = sib_ptr->next;
      rectrav(sib_ptr->back, m, n);
    }
  }
  mx1 = mx;
}  /* rectrav */


static void summarize(FILE *fp)
{
  long i, j, mm;
  double mode, sum;
  double like[maxcategs], nulike[maxcategs];
  double **marginal;

  mp = (long **)Malloc(sites * sizeof(long *));
  for (i = 0; i <= sites-1; ++i)
    mp[i] = (long *)Malloc(sizeof(long)*rcategs);
  fprintf(fp, "\nLn Likelihood = %11.5f\n\n", curtree.likelihood);
  fprintf(fp, " Ancestor      Node      Node Height     Length\n");
  fprintf(fp, " --------      ----      ---- ------     ------\n");
  mlk_describe(fp, &curtree, fracchange);
  putc('\n', fp);
  if (rctgry && rcategs > 1) {
    for (i = 0; i < rcategs; i++)
      like[i] = 1.0;
    for (i = sites - 1; i >= 0; i--) {
      sum = 0.0;
      for (j = 0; j < rcategs; j++) {
        nulike[j] = (lambda1 + lambda * probcat[j]) * like[j];
        mp[i][j] = j + 1;
        for (k = 1; k <= rcategs; k++) {
          if (k != j + 1) {
            if (lambda * probcat[k - 1] * like[k - 1] > nulike[j]) {
              nulike[j] = lambda * probcat[k - 1] * like[k - 1];
              mp[i][j] = k;
            }
          }
        }
        if ((ally[i] > 0) && (location[ally[i]-1] > 0))
          nulike[j] *= contribution[location[ally[i] - 1] - 1][j];
        sum += nulike[j];
      }
      for (j = 0; j < rcategs; j++)
        nulike[j] /= sum;
      memcpy(like, nulike, rcategs * sizeof(double));
    }
    mode = 0.0;
    mx = 1;
    for (i = 1; i <= rcategs; i++) {
      if (probcat[i - 1] * like[i - 1] > mode) {
        mx = i;
        mode = probcat[i - 1] * like[i - 1];
      }
    }
    mx0 = mx;
    fprintf(fp,
 "Combination of categories that contributes the most to the likelihood:\n\n");
    for (i = 1; i <= nmlngth + 3; i++)
      putc(' ', fp);
    for (i = 1; i <= sites; i++) {
      fprintf(fp, "%ld", mx);
      if (i % 10 == 0)
        putc(' ', fp);
      if (i % 60 == 0 && i != sites) {
        putc('\n', fp);
        for (j = 1; j <= nmlngth + 3; j++)
          putc(' ', fp);
      }
      mx = mp[i - 1][mx - 1];
    }
    fprintf(fp, "\n\n");
    marginal = (double **) Malloc( sites*sizeof(double *));
    for (i = 0; i < sites; i++)
      marginal[i] = (double *) Malloc( rcategs*sizeof(double));
    for (i = 0; i < rcategs; i++)
      like[i] = 1.0;
    for (i = sites - 1; i >= 0; i--) {
      sum = 0.0;
      for (j = 0; j < rcategs; j++) {
        nulike[j] = (lambda1 + lambda * probcat[j]) * like[j];
        for (k = 1; k <= rcategs; k++) {
          if (k != j + 1)
              nulike[j] += lambda * probcat[k - 1] * like[k - 1];
        }
        if ((ally[i] > 0) && (location[ally[i]-1] > 0))
          nulike[j] *= contribution[location[ally[i] - 1] - 1][j];
        sum += nulike[j];
      }
      for (j = 0; j < rcategs; j++) {
        nulike[j] /= sum;
        marginal[i][j] = nulike[j];
      }
      memcpy(like, nulike, rcategs * sizeof(double));
    }
    for (i = 0; i < rcategs; i++)
      like[i] = 1.0;
    for (i = 0; i < sites; i++) {
      sum = 0.0;
      for (j = 0; j < rcategs; j++) {
        nulike[j] = (lambda1 + lambda * probcat[j]) * like[j];
        for (k = 1; k <= rcategs; k++) {
          if (k != j + 1)
              nulike[j] += lambda * probcat[k - 1] * like[k - 1];
        }
        marginal[i][j] *= like[j] * probcat[j];
        sum += nulike[j];
      }
      for (j = 0; j < rcategs; j++)
        nulike[j] /= sum;
      memcpy(like, nulike, rcategs * sizeof(double));
      sum = 0.0;
      for (j = 0; j < rcategs; j++)
        sum += marginal[i][j];
      for (j = 0; j < rcategs; j++)
        marginal[i][j] /= sum;
    }
    fprintf( fp, "Most probable category at each site if > 0.95 probability "
        "(\".\" otherwise)\n\n" );
    for (i = 1; i <= nmlngth + 3; i++)
      putc(' ', fp);
    for (i = 0; i < sites; i++) {
      mm = 0;
      sum = 0.0;
      for (j = 0; j < rcategs; j++)
        if (marginal[i][j] > sum) {
          sum = marginal[i][j];
          mm = j;
        }
        if (sum >= 0.95)
        fprintf(fp, "%ld", mm+1);
      else
        putc('.', fp);
      if ((i+1) % 60 == 0) {
        if (i != 0) {
          putc('\n', fp);
          for (j = 1; j <= nmlngth + 3; j++)
            putc(' ', fp);
        }
      }
      else if ((i+1) % 10 == 0)
        putc(' ', fp);
    }
    putc('\n', fp);
    for (i = 0; i < sites; i++)
      free(marginal[i]);
    free(marginal);
  }
  putc('\n', fp);
  putc('\n', fp);
  if (hypstate) {
    fprintf(fp, "Probable sequences at interior nodes:\n\n");
    fprintf(fp, "  node       ");
    for (i = 0; (i < 13) && (i < ((sites + (sites-1)/10 - 39) / 2)); i++)
      putc(' ', fp);
    fprintf(fp, "Reconstructed sequence (caps if > 0.95)\n\n");
    if (!rctgry || (rcategs == 1))
      mx0 = 1;
    for (i = 0; i < sites; i += 60) {
      k = i + 59;
      if (k >= sites)
        k = sites - 1;
      rectrav(curtree.root, i, k);
      putc('\n', fp);
      mx0 = mx1;
    }
  }
  for (i = 0; i < sites; ++i)
    free(mp[i]);
  free(mp);
}  /* summarize */


static void dnamlk_treeout(node *p)
{
  /* write out file with representation of final tree */
  node *sib_ptr;
  long i, n, w, num_sibs;
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
    sib_ptr = p;
    num_sibs = count_sibs(p);
    putc('(', outtree);
    col++;

    for (i=0; i < (num_sibs - 1); i++) {
      sib_ptr = sib_ptr->next;
      dnamlk_treeout(sib_ptr->back);
      putc(',', outtree);
      col++;
      if (col > 55) {
        putc('\n', outtree);
        col = 0;
      }
    }
    sib_ptr = sib_ptr->next;
    dnamlk_treeout(sib_ptr->back);
    putc(')', outtree);
    col++;
  }
  if (p == curtree.root) {
    fprintf(outtree, ";\n");
    return;
  }
  x = fracchange * (p->tyme - curtree.nodep[p->back->index - 1]->tyme);
  if (x > 0.0)
    w = (long)(0.4342944822 * log(x));
  else if (x == 0.0)
    w = 0;
  else
    w = (long)(0.4342944822 * log(-x)) + 1;
  if (w < 0)
    w = 0;
  fprintf(outtree, ":%*.5f", (int)(w + 7), x);
  col += w + 8;
}  /* dnamlk_treeout */


static void init_tymes(node *p, double minlength)
{
  /* Set all node tymes closest to the tips but with no branches shorter than
   * minlength */

  long i, num_sibs;
  node *sib_ptr, *sib_back_ptr;

  /* traverse to set up times in subtrees */
  if (p->tip)
    return;

  sib_ptr = p;
  num_sibs = count_sibs(p);
  for (i=0 ; i < num_sibs; i++) {
    sib_ptr      = sib_ptr->next;
    sib_back_ptr = sib_ptr->back;
    init_tymes(sib_back_ptr, minlength);
  }

  /* set time at this node */
  setnodetymes(p, min_child_tyme(p) - minlength);

}  /* init_tymes */


static void treevaluate()
{
  /* evaluate likelihood of tree, after iterating branch lengths */
  long i;

  if ( !usertree || (usertree && !lngths) ) {
    polishing = true;
    smoothit = true;
    for (i = 0; i < smoothings; ) {
      if ( !smooth(curtree.root) )
        i++;
    }
  }
  dnamlk_evaluate(curtree.root);
}  /* treevaluate */


static void maketree()
{
  /* constructs a binary tree from the pointers in curtree.nodep,
     adds each node at location which yields highest likelihood
     then rearranges the tree for greatest likelihood */

  long i, j, numtrees;
  node *item, *nufork, *dummy, *q, *root=NULL;
  boolean succeded, dummy_haslengths, dummy_first, goteof;
  long max_nonodes;        /* Maximum number of nodes required to
                            * express all species in a bifurcating tree
                            * */
  long nextnode;
  pointarray dummy_treenode=NULL;
  double oldbest;
  node *tmp;

  inittable();  

  if (!usertree) {
    for (i = 1; i <= spp; i++)
      enterorder[i - 1] = i;
    if (jumble)
      randumize(seed, enterorder);
    curtree.root = curtree.nodep[spp];
    curtree.root->back = NULL;
    for (i = 0; i < spp; i++)
       curtree.nodep[i]->back = NULL;
    for (i = spp; i < nonodes; i++) {
       q = curtree.nodep[i];
       q->back = NULL;
       while ((q = q->next) != curtree.nodep[i])
         q->back = NULL;
    }
    polishing = false;
    dnamlk_add(curtree.nodep[enterorder[0] - 1], curtree.nodep[enterorder[1] - 1],
        curtree.nodep[spp]);
    if (progress) {
      printf("\nAdding species:\n");
      writename(0, 2, enterorder);
#ifdef WIN32
      phyFillScreenColor();
#endif
    }
    lastsp = false;
    smoothit = false;
    for (i = 3; i <= spp; i++) {
      bestree.likelihood = UNDEFINED;
      bestyet = UNDEFINED;
      there = curtree.root;
      item = curtree.nodep[enterorder[i - 1] - 1];
      nufork = curtree.nodep[spp + i - 2];
      lastsp = (i == spp);
      addpreorder(curtree.root, item, nufork, true, true);
      dnamlk_add(there, item, nufork);
      like = dnamlk_evaluate(curtree.root);
      copy_(&curtree, &bestree, nonodes, rcategs);
      rearrange(&curtree.root);
      if (curtree.likelihood > bestree.likelihood) {
        copy_(&curtree, &bestree, nonodes, rcategs);
      }
      if (progress) {
        writename(i - 1, 1, enterorder);
#ifdef WIN32
        phyFillScreenColor();
#endif
      }
      if (lastsp && global) {
        /* perform global rearrangements */
        if (progress) {
          printf("Doing global rearrangements\n");
          printf("  !");
          for (j = 1; j <= nonodes; j++)
            if ( j % (( nonodes / 72 ) + 1 ) == 0 )
              putchar('-');
          printf("!\n");
        }
        global2 = false;
        do {
          succeded = false;
          if (progress)
            printf("   ");
          /* FIXME: tymes gets clobbered by tryadd() */
          /* save_tymes(&curtree, tymes); */
          for (j = 0; j < nonodes; j++) {
            oldbest = bestree.likelihood;
            bestyet = UNDEFINED;
            item = curtree.nodep[j];
            if (item != curtree.root) {
              nufork = pnode(&curtree, item->back); /* parent fork */
              
              if (nufork != curtree.root) {  
                tmp = nufork->next->back;
                if (tmp == item) 
                    tmp = nufork->next->next->back; 
                    /* can't figure out why we never get here */
              }
              else {
                  if (nufork->next->back != item)
                      tmp  = nufork->next->back;
                  else tmp = nufork->next->next->back;
              } /* if we add item at tmp we have done nothing */
              assert( all_tymes_valid(curtree.root, 0.98*MIN_BRANCH_LENGTH, false) );
              dnamlk_re_move(&item, &nufork, false);
              /* there = curtree.root; */
              there = tmp;
              addpreorder(curtree.root, item, nufork, true, true);
              if ( tmp != there && bestree.likelihood > oldbest)
                succeded = true;
              dnamlk_add(there, item, nufork);
              if (global2)
                restore_tymes(&curtree,tymes);
            }
            if (progress) {
              if ( j % (( nonodes / 72 ) + 1 ) == 0 )
                putchar('.');
              fflush(stdout);
            }
          } 
          if (progress)
            putchar('\n');
        } while ( succeded );
      }
    }
    if (njumble > 1 && lastsp) {
      for (i = 0; i < spp; i++ )
        dnamlk_re_move(&curtree.nodep[i], &dummy, false);
      if (jumb == 1 || bestree2.likelihood < bestree.likelihood)
        copy_(&bestree, &bestree2, nonodes, rcategs);
    }
    if (jumb == njumble) {
      if (njumble > 1)
        copy_(&bestree2, &curtree, nonodes, rcategs);
      else copy_(&bestree, &curtree, nonodes, rcategs);
      fprintf(outfile, "\n\n");
      treevaluate();
      curtree.likelihood = dnamlk_evaluate(curtree.root);
      if (treeprint)
        mlk_printree(outfile, &curtree);
      summarize(outfile);
      if (trout) {
        col = 0;
        dnamlk_treeout(curtree.root);
      }
    } 
  } else { /* if ( usertree ) */
    /* Open in binary: ftell() is broken for UNIX line-endings under WIN32 */
    openfile(&intree, INTREE, "input tree file", "rb", progname, intreename);
    inittable_for_usertree (intree);
    numtrees = countsemic(&intree);
    if(numtrees > MAXSHIMOTREES)
      shimotrees = MAXSHIMOTREES;
    else
      shimotrees = numtrees;
    if (numtrees > 2)
      initseed(&inseed, &inseed0, seed);
    l0gl = (double *)Malloc(shimotrees * sizeof(double));
    l0gf = (double **)Malloc(shimotrees * sizeof(double *));
    for (i=0; i < shimotrees; ++i)
      l0gf[i] = (double *)Malloc(endsite * sizeof(double));
    if (treeprint) {
      fprintf(outfile, "User-defined tree");
      if (numtrees > 1)
        putc('s', outfile);
      fprintf(outfile, ":\n\n");
    }
    fprintf(outfile, "\n\n");
    which = 1;
    max_nonodes = nonodes;
    while (which <= numtrees) {
      
      /* These initializations required each time through the loop
         since multiple trees require re-initialization */
      dummy_haslengths = true;
      nextnode         = 0;
      dummy_first      = true;
      goteof           = false;
      lngths           = lengthsopt;
      nonodes          = max_nonodes;

      treeread(intree, &root, dummy_treenode, &goteof, &dummy_first,
               curtree.nodep, &nextnode, &dummy_haslengths, &grbg, 
               initdnamlnode, false, nonodes);

      if (goteof && (which <= numtrees)) {
        /* if we hit the end of the file prematurely */
        printf ("\n");
        printf ("ERROR: trees missing at end of file.\n");
        printf ("\tExpected number of trees:\t\t%ld\n", numtrees);
        printf ("\tNumber of trees actually in file:\t%ld.\n\n", which - 1);
        exxit(-1);
      }

      nonodes = nextnode;
      root = curtree.nodep[root->index - 1];
      curtree.root = root;

      if (lngths)
        tymetrav(curtree.root);
      else
        init_tymes(curtree.root, initialv);

      treevaluate();
      if (treeprint)
        mlk_printree(outfile, &curtree);
      summarize(outfile);
      if (trout) {
        col = 0;
        dnamlk_treeout(curtree.root);
      }
      if(which < numtrees){
        freex_notip(nonodes, curtree.nodep);
        gdispose(curtree.root, &grbg, curtree.nodep);
      }
      which++;
    }

    FClose(intree);
    if (!auto_ && numtrees > 1 && weightsum > 1 )
      standev2(numtrees, maxwhich, 0, endsite, maxlogl, l0gl, l0gf,
               aliasweight, seed);
  }
  
  if (jumb == njumble) {
    if (progress) {
      printf("\nOutput written to file \"%s\"\n\n", outfilename);
      if (trout)
        printf("Tree also written onto file \"%s\"\n\n", outtreename);
    }
    free(contribution);
    freex(nonodes, curtree.nodep);
    if (!usertree) {
      freex(nonodes, bestree.nodep);
      if (njumble > 1)
        freex(nonodes, bestree2.nodep);
    }
  }
  free(root);
  freetable();
} /* maketree */

/*?? Dnaml has a clean-up function for freeing memory, closing files, etc.
     Put one here too? */

int main(int argc, Char *argv[])
{  /* DNA Maximum Likelihood with molecular clock */
  /* Initialize mlclock.c */
  mlclock_init(&curtree, &dnamlk_evaluate);

#ifdef MAC
  argc = 1;                /* macsetup("Dnamlk", "Dnamlk");        */
  argv[0] = "Dnamlk";
#endif
  init(argc,argv);
  progname = argv[0];
  openfile(&infile, INFILE, "input file", "r", argv[0], infilename);
  openfile(&outfile, OUTFILE, "output file", "w", argv[0], outfilename);

  ibmpc = IBMCRT;
  ansi = ANSICRT;
  datasets = 1;
  mulsets = false;
  firstset = true;
  
  doinit();

  ttratio0 = ttratio;

  /* Open output tree, categories, and weights files if needed */
  if ( trout )
    openfile(&outtree, OUTTREE, "output tree file", "w", argv[0], outtreename);
  if ( ctgry )
    openfile(&catfile, CATFILE, "categories file", "r", argv[0], catfilename);
  if ( weights || justwts )
   openfile(&weightfile, WEIGHTFILE, "weights file", "r", argv[0], weightfilename);

  /* Data set loop */
  for ( ith = 1; ith <= datasets; ith++ ) {
    ttratio = ttratio0;
    if ( datasets > 1 ) {
      fprintf(outfile, "Data set # %ld:\n\n", ith);
      if ( progress )
        printf("\nData set # %ld:\n", ith);
    }
    getinput();
    if ( ith == 1 )
      firstset = false;
    
    /* Jumble loop */
    if (usertree)
      maketree();
    else
      for (jumb = 1; jumb <= njumble; jumb++)
        maketree();
  }

  /* Close files */
  FClose(infile);  
  FClose(outfile);
  if ( trout )
    FClose(outtree);
  if ( ctgry )
    FClose(catfile);
  if ( weights || justwts )
    FClose(weightfile);

#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(outtreename);
#endif
  printf("Done.\n\n");
#ifdef WIN32
  phyRestoreConsoleAttributes();
#endif
  return 0;
}  /* DNA Maximum Likelihood with molecular clock */
