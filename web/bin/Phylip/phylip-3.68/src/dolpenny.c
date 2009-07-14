#include "phylip.h"
#include "disc.h"
#include "dollo.h"

/* version 3.6. (c) Copyright 1993-2004 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define maxtrees        1000  /* maximum number of trees to be printed out   */
#define often           100   /* how often to notify how many trees examined */
#define many            1000  /* how many multiples of howoften before stop  */

typedef long *treenumbers;
typedef double *valptr;
typedef long *placeptr;

#ifndef OLDC
/* function prototypes */
void   getoptions(void);
void   allocrest(void);
void   doinit(void);
void   inputoptions(void);
void   doinput(void);
void   preorder(node *);
void   evaluate(node *);
void   addtraverse(node *, node *, node *, placeptr, valptr, long *);
void   addit(long);
void   describe(void);
void   maketree(void);
/* function prototypes */
#endif

Char infilename[FNMLNGTH], outfilename[FNMLNGTH], outtreename[FNMLNGTH],  weightfilename[FNMLNGTH], ancfilename[FNMLNGTH];
node *root;
long howmanny, howoften, col, msets, ith;
boolean weights, thresh, ancvar, questions, dollo, simple, trout,
               progress, treeprint, stepbox, ancseq,
               mulsets, firstset, justwts;
boolean *ancone, *anczero, *ancone0, *anczero0;
pointptr treenode;   /* pointers to all nodes in tree */
double fracdone, fracinc;
double threshold;
double *threshwt;
boolean *added;
Char *guess;
steptr numsteps, numsone, numszero;
gbit *garbage;
long **bestorders, **bestrees;

/* Local variables for maketree, propagated globally for C version: */
long examined, mults;
boolean firsttime, done;
double like, bestyet;
treenumbers current, order;
long fullset;
bitptr zeroanc, oneanc;
bitptr stps;


void getoptions()
{
  /* interactively set options */
  long loopcount, loopcount2;
  Char ch, ch2;
  boolean done;

  fprintf(outfile, "\nPenny algorithm for Dollo or polymorphism");
  fprintf(outfile, " parsimony, version %s\n",VERSION);
  fprintf(outfile, " branch-and-bound to find all");
  fprintf(outfile, " most parsimonious trees\n\n");
  howoften = often;
  howmanny = many;
  simple = true;
  thresh = false;
  threshold = spp;
  trout = true;
  weights = false;
  justwts = false;
  ancvar = false;
  dollo = true;
  printdata = false;
  progress = true;
  treeprint = true;
  stepbox = false;
  ancseq = false;
  loopcount = 0;
  do {
    cleerhome();
    printf("\nPenny algorithm for Dollo or polymorphism parsimony,");
    printf(" version %s\n",VERSION);
    printf(" branch-and-bound to find all most parsimonious trees\n\n");
    printf("Settings for this run:\n");
    printf("  P                     Parsimony method?  %s\n",
           (dollo ? "Dollo" : "Polymorphism"));
    printf("  H        How many groups of %4ld trees:%6ld\n",howoften,howmanny);
    printf("  F        How often to report, in trees:%5ld\n",howoften);
    printf("  S           Branch and bound is simple?  %s\n",
           (simple ? "Yes" : "No.  Reconsiders order of species"));
    printf("  T              Use Threshold parsimony?");
    if (thresh)
      printf("  Yes, count steps up to%4.1f per char.\n", threshold);
    else
      printf("  No, use ordinary parsimony\n");
    printf("  A                 Use ancestral states?  %s\n",
           (ancvar ? "Yes" : "No"));
    printf("  W                       Sites weighted?  %s\n",
           (weights ? "Yes" : "No"));
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
        printf("  Yes, %2ld %s\n", msets,
               (justwts ? "sets of weights" : "data sets"));
    else
      printf("  No\n");
    printf("  0   Terminal type (IBM PC, ANSI, none)?  %s\n",
           (ibmpc ? "IBM PC" : ansi  ? "ANSI" : "(none)"));
    printf("  1    Print out the data at start of run  %s\n",
           (printdata ? "Yes" : "No"));
    printf("  2  Print indications of progress of run  %s\n",
           (progress ? "Yes" : "no"));
    printf("  3                        Print out tree  %s\n",
           (treeprint ? "Yes": "No"));
    printf("  4     Print out steps in each character  %s\n",
           (stepbox ? "Yes" : "No"));
    printf("  5     Print states at all nodes of tree  %s\n",
           (ancseq ? "Yes" : "No"));
    printf("  6       Write out trees onto tree file?  %s\n",
           (trout ? "Yes" : "No"));
    if(weights && justwts){
        printf(
         "WARNING:  W option and Multiple Weights options are both on.  ");
        printf(
         "The W menu option is unnecessary and has no additional effect. \n");
    }
    printf("\nAre these settings correct? (type Y or the letter for one to change)\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    fflush(stdout);
    scanf("%c%*[^\n]", &ch);
    getchar();
    uppercase(&ch);
    done = (ch == 'Y');
    if (!done) {
      if (strchr("WHMSTAPF1234560",ch) != NULL){
        switch (ch) {

        case 'W':
        weights = !weights;
        break;

        case 'H':
          inithowmany(&howmanny, howoften);
          break;

        case 'F':
          inithowoften(&howoften);
          break;

        case 'A':
          ancvar = !ancvar;
          break;

        case 'P':
          dollo = !dollo;
          break;

        case 'S':
          simple = !simple;
          break;

        case 'T':
          thresh = !thresh;
          if (thresh)
            initthreshold(&threshold);
          break;

        case 'M':
          mulsets = !mulsets;
          if (mulsets){
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
            justweights(&msets);
          else
            initdatasets(&msets);
          }
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
          stepbox = !stepbox;
          break;

        case '5':
          ancseq = !ancseq;
          break;

        case '6':
          trout = !trout;
          break;
        }
      } else
        printf("Not a possible option!\n");
    }
    countup(&loopcount, 100);
  } while (!done);
}  /* getoptions */


void allocrest()
{
  long i;

  extras = (long *)Malloc(chars*sizeof(long));
  weight = (long *)Malloc(chars*sizeof(long));
  threshwt = (double *)Malloc(chars*sizeof(double));
  guess = (Char *)Malloc(chars*sizeof(Char));
  numsteps = (long *)Malloc(chars*sizeof(long));
  numszero = (long *)Malloc(chars*sizeof(long));
  numsone = (long *)Malloc(chars*sizeof(long));
  bestorders = (long **)Malloc(maxtrees*sizeof(long *));
  bestrees = (long **)Malloc(maxtrees*sizeof(long *));
  for (i = 1; i <= maxtrees; i++) {
    bestorders[i - 1] = (long *)Malloc(spp*sizeof(long));
    bestrees[i - 1] = (long *)Malloc(spp*sizeof(long));
  }
  current = (treenumbers)Malloc(spp*sizeof(long));
  order = (treenumbers)Malloc(spp*sizeof(long));
  nayme = (naym *)Malloc(spp*sizeof(naym));
  added = (boolean *)Malloc(nonodes*sizeof(boolean));
  ancone = (boolean *)Malloc(chars*sizeof(boolean));
  anczero = (boolean *)Malloc(chars*sizeof(boolean));
  ancone0 = (boolean *)Malloc(chars*sizeof(boolean));
  anczero0 = (boolean *)Malloc(chars*sizeof(boolean));
  zeroanc = (bitptr)Malloc(words*sizeof(long));
  oneanc = (bitptr)Malloc(words*sizeof(long));
}  /* allocrest */


void doinit()
{
  /* initializes variables */

  inputnumbers(&spp, &chars, &nonodes, 1);
  words = chars / bits + 1;
  getoptions();
  if (printdata)
    fprintf(outfile, "%2ld species, %3ld characters\n", spp, chars);
  alloctree(&treenode);
  setuptree(treenode);
  allocrest();
}  /* doinit */


void inputoptions()
{
  /* input the information on the options */
  long i;

  if(justwts){
      if(firstset){
          scan_eoln(infile);
          if (ancvar) {
              inputancestors(anczero0, ancone0);
          }
      }
      for (i = 0; i < (chars); i++)
          weight[i] = 1;
      inputweights(chars, weight, &weights);
  }
  else {
      if (!firstset)
          samenumsp(&chars, ith);
      scan_eoln(infile);
      for (i = 0; i < (chars); i++)
          weight[i] = 1;
      if (ancvar)
          inputancestors(anczero0, ancone0);
      if (weights)
          inputweights(chars, weight, &weights);
  }
  for (i = 0; i < (chars); i++) {
    if (!ancvar) {
      anczero[i] = true;
      ancone[i] = false;
    } else {
      anczero[i] = anczero0[i];
      ancone[i] = ancone0[i];
    }
  }
  questions = false;
  if (!thresh)
    threshold = spp;
  for (i = 0; i < (chars); i++) {
    questions = (questions || (ancone[i] && anczero[i]));
    threshwt[i] = threshold * weight[i];
  }
}  /* inputoptions */


void doinput()
{
  /* reads the input data */
  inputoptions();
  if(!justwts || firstset)
      inputdata(treenode, dollo, printdata, outfile);
}  /* doinput */


void preorder(node *p)
{
  /* go back up tree setting up and counting interior node
     states */
  long i;

  if (!p->tip) {
    correct(p, fullset, dollo, zeroanc, treenode);
    preorder(p->next->back);
    preorder(p->next->next->back);
  }
  if (p->back == NULL)
    return;
  if (dollo) {
    for (i = 0; i < (words); i++)
      stps[i] = (treenode[p->back->index - 1]->stateone[i] & p->statezero[i] &
                 zeroanc[i]) |
                (treenode[p->back->index - 1]->statezero[i] & p->stateone[i] &
                 fullset & (~zeroanc[i]));
  } else {
    for (i = 0; i < (words); i++)
      stps[i] = treenode[p->back->index - 1]->stateone[i] &
                treenode[p->back->index - 1]->statezero[i] & p->stateone[i] &
                p->statezero[i];
  }
  count(stps, zeroanc, numszero, numsone);
}  /* preorder */


void evaluate(node *r)
{
  /* Determines the number of losses or polymorphisms needed
     for a tree. This is the minimum number needed to evolve
     chars on this tree */
  long i, stepnum, smaller;
  double sum;

  sum = 0.0;
  for (i = 0; i < (chars); i++) {
    numszero[i] = 0;
    numsone[i] = 0;
  }
  for (i = 0; i < (words); i++)
    zeroanc[i] = fullset;
  postorder(r);
  preorder(r);
  for (i = 0; i < (words); i++)
    zeroanc[i] = 0;
  postorder(r);
  preorder(r);
  for (i = 0; i < (chars); i++) {
    smaller = spp * weight[i];
    numsteps[i] = smaller;
    if (anczero[i]) {
      numsteps[i] = numszero[i];
      smaller = numszero[i];
    }
    if (ancone[i] && numsone[i] < smaller)
      numsteps[i] = numsone[i];
    stepnum = numsteps[i] + extras[i];
    if (stepnum <= threshwt[i])
      sum += stepnum;
    else
      sum += threshwt[i];
    guess[i] = '?';
    if (!ancone[i] || (anczero[i] && numszero[i] < numsone[i]))
      guess[i] = '0';
    else if (!anczero[i] || (ancone[i] && numsone[i] < numszero[i]))
      guess[i] = '1';
  }
  if (examined == 0 && mults == 0)
    bestyet = -1.0;
  like = sum;
}  /* evaluate */


void addtraverse(node *a, node *b, node *c, placeptr place,
                valptr valyew, long *n)
{
  /* traverse all places to add b */
  if (done)
    return;
  add(a, b, c, &root, treenode);
  (*n)++;
  evaluate(root);
  examined++;
  if (examined == howoften) {
    examined = 0;
    mults++;
    if (mults == howmanny)
      done = true;
    if (progress) {
      printf("%6ld",mults);
      if (bestyet >= 0)
        printf("%18.5f", bestyet);
      else
        printf("         -        ");
      printf("%17ld%20.2f\n", nextree - 1, fracdone * 100);
#ifdef WIN32
      phyFillScreenColor();
#endif
    }
  }
  valyew[*n - 1] = like;
  place[*n - 1] = a->index;
  re_move(&b, &c, &root, treenode);
  if (!a->tip) {
    addtraverse(a->next->back, b, c, place,valyew,n);
    addtraverse(a->next->next->back, b, c, place,valyew,n);
  }
}  /* addtraverse */


void addit(long m)
{
  /* adds the species one by one, recursively */
  long n;
  valptr valyew;
  placeptr place;
  long i, j, n1, besttoadd = 0;
  valptr bestval;
  placeptr bestplace;
  double oldfrac, oldfdone, sum, bestsum;

  valyew = (valptr)Malloc(nonodes*sizeof(double));
  bestval = (valptr)Malloc(nonodes*sizeof(double));
  place = (placeptr)Malloc(nonodes*sizeof(long));
  bestplace = (placeptr)Malloc(nonodes*sizeof(long));
  if (simple && !firsttime) {
    n = 0;
    added[order[m - 1] - 1] = true;
    addtraverse(root, treenode[order[m - 1] - 1], treenode[spp + m - 2],
                place, valyew, &n);
    besttoadd = order[m - 1];
    memcpy(bestplace, place, nonodes*sizeof(long));
    memcpy(bestval, valyew, nonodes*sizeof(double));
  } else {
    bestsum = -1.0;
    for (i = 1; i <= (spp); i++) {
      if (!added[i - 1]) {
        n = 0;
        added[i - 1] = true;
        addtraverse(root, treenode[i - 1], treenode[spp + m - 2],
                    place, valyew, &n);
        added[i - 1] = false;
        sum = 0.0;
        for (j = 0; j < n; j++)
          sum += valyew[j];
        if (sum > bestsum) {
          bestsum = sum;
          besttoadd = i;
          memcpy(bestplace, place, nonodes*sizeof(long));
          memcpy(bestval, valyew, nonodes*sizeof(double));
        }
      }
    }
  }
  order[m - 1] = besttoadd;
  memcpy(place, bestplace, nonodes*sizeof(long));
  memcpy(valyew, bestval, nonodes*sizeof(double));
  shellsort(valyew, place, n);
  oldfrac = fracinc;
  oldfdone = fracdone;
  n1 = 0;
  for (i = 0; i < (n); i++) {
    if (valyew[i] <= bestyet || bestyet < 0.0)
      n1++;
  }
  if (n1 > 0)
    fracinc /= n1;
  for (i = 0; i < n; i++) {
    if (valyew[i] <=bestyet ||bestyet < 0.0) {
      current[m - 1] = place[i];
      add(treenode[place[i] - 1], treenode[besttoadd - 1],
          treenode[spp + m - 2], &root, treenode);
      added[besttoadd - 1] = true;
      if (m < spp)
        addit(m + 1);
      else {
        if (valyew[i] < bestyet || bestyet < 0.0) {
          nextree = 1;
          bestyet = valyew[i];
        }
        if (nextree <= maxtrees) {
          memcpy(bestorders[nextree - 1], order,
                 spp*sizeof(long));
          memcpy(bestrees[nextree - 1], current,
                 spp*sizeof(long));
        }
        nextree++;
        firsttime = false;
      }
      re_move(&treenode[besttoadd - 1], &treenode[spp + m - 2], &root, treenode);
      added[besttoadd - 1] = false;
    }
    fracdone += fracinc;
  }
  fracinc = oldfrac;
  fracdone = oldfdone;
  free(valyew);
  free(bestval);
  free(place);
  free(bestplace);
}  /* addit */


void describe()
{
  /* prints ancestors, steps and table of numbers of steps in
     each character */
  if (stepbox) {
    putc('\n', outfile);
    writesteps(weights, dollo, numsteps);
  }
  if (questions)
    guesstates(guess);
  if (ancseq) {
    hypstates(fullset, dollo, guess, treenode, root, garbage, zeroanc, oneanc);
    putc('\n', outfile);
  }
  putc('\n', outfile);
  if (trout) {
    col = 0;
    treeout(root, nextree, &col, root);
  }
}  /* describe */


void maketree()
{
  /* tree construction recursively by branch and bound */
  long i, j, k;
  node *dummy;

  fullset = (1L << (bits + 1)) - (1L << 1);
  if (progress) {
    printf("\nHow many\n");
    printf("trees looked                                       Approximate\n");
    printf("at so far      Length of        How many           percentage\n");
    printf("(multiples     shortest tree    trees this long    searched\n");
    printf("of %4ld):      found so far     found so far       so far\n",
           howoften);
    printf("----------     ------------     ------------       ------------\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
  }
  done = false;
  mults = 0;
  examined = 0;
  nextree = 1;
  root = treenode[0];
  firsttime = true;
  for (i = 0; i < (spp); i++)
    added[i] = false;
  added[0] = true;
  order[0] = 1;
  k = 2;
  fracdone = 0.0;
  fracinc = 1.0;
  bestyet = -1.0;
  stps = (bitptr)Malloc(words*sizeof(long));
  addit(k);
  if (done) {
    if (progress) {
      printf("Search broken off!  Not guaranteed to\n");
      printf(" have found the most parsimonious trees.\n");
    }
    if (treeprint) {
      fprintf(outfile, "Search broken off!  Not guaranteed to\n");
      fprintf(outfile, " have found the most parsimonious\n");
      fprintf(outfile, " trees, but here is what we found:\n");
    }
  }
  if (treeprint) {
    fprintf(outfile, "\nrequires a total of %18.3f\n\n", bestyet);
    if (nextree == 2)
      fprintf(outfile, "One most parsimonious tree found:\n");
    else
      fprintf(outfile, "%5ld trees in all found\n", nextree - 1);
  }
  if (nextree > maxtrees + 1) {
    if (treeprint)
      fprintf(outfile, "here are the first%4ld of them\n", (long)maxtrees);
    nextree = maxtrees + 1;
  }
  if (treeprint)
    putc('\n', outfile);
  for (i = 0; i < (spp); i++)
    added[i] = true;
  for (i = 0; i <= (nextree - 2); i++) {
    for (j = k; j <= (spp); j++)
      add(treenode[bestrees[i][j - 1] - 1],
          treenode[bestorders[i][j - 1] - 1], treenode[spp + j - 2],
          &root, treenode);
    evaluate(root);
    printree(1.0, treeprint, root);
    describe();
    for (j = k - 1; j < (spp); j++)
      re_move(&treenode[bestorders[i][j] - 1], &dummy, &root, treenode);
  }
  if (progress) {
    printf("\nOutput written to file \"%s\"\n\n", outfilename);
    if (trout)
      printf("Trees also written onto file \"%s\"\n\n", outtreename);
  }
  free(stps);
  if (ancseq)
    freegarbage(&garbage);
}  /* maketree */


int main(int argc, Char *argv[])
{ /* branch-and-bound method for Dollo, polymorphism parsimony */
  /* Reads in the number of species, number of characters,
     options and data.  Then finds all most parsimonious trees */
#ifdef MAC
  argc = 1;                /* macsetup("Dolpenny","");                */
  argv[0] = "Dolpenny";
#endif
  init(argc, argv);
  openfile(&infile,INFILE,"input file", "r",argv[0],infilename);
  openfile(&outfile,OUTFILE,"output file", "w",argv[0],outfilename);
  ibmpc = IBMCRT;
  ansi = ANSICRT;
  garbage = NULL;
  mulsets = false;
  msets = 1;
  firstset = true;
  bits = 8*sizeof(long) - 1;
  doinit();
  if (weights || justwts)
    openfile(&weightfile,WEIGHTFILE,"weights file","r",argv[0],weightfilename);
  if (trout)
    openfile(&outtree,OUTTREE,"output tree file", "w",argv[0],outtreename);
  if(ancvar)
      openfile(&ancfile,ANCFILE,"ancestors file", "r",argv[0],ancfilename);
      fprintf(outfile,"%s parsimony method\n\n",dollo ? "Dollo" : "Polymorphism");  
  for (ith = 1; ith <= msets; ith++) {
    doinput();
    if (msets > 1 && !justwts) {
      fprintf(outfile, "Data set # %ld:\n\n",ith);
      if (progress)
        printf("\nData set # %ld:\n",ith);
    } 
    if (justwts){
        fprintf(outfile, "Weights set # %ld:\n\n", ith);
        if (progress)
            printf("\nWeights set # %ld:\n\n", ith);
    }
    if (printdata){
        if (weights || justwts)
            printweights(outfile, 0, chars, weight, "Characters");
        if (ancvar)
            printancestors(outfile, anczero, ancone);
    }
    if (ith == 1)
      firstset = false;
    maketree();
  }

  FClose(infile);
  FClose(outfile);
  FClose(outtree);
#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(outtreename);
#endif
#ifdef WIN32
  phyRestoreConsoleAttributes();
#endif
  return 0;
}  /* branch-and-bound method for Dollo, polymorphism parsimony */
