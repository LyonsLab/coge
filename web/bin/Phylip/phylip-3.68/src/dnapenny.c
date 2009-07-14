#include "phylip.h"
#include "seq.h"

/* version 3.6. (c) Copyright 1993-2004 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define maxtrees        1000  /* maximum number of trees to be printed out   */
#define often           100   /* how often to notify how many trees examined */
#define many            1000  /* how many multiples of howoften before stop  */

typedef node **pointptr;
typedef long *treenumbers;
typedef double *valptr;
typedef long *placeptr;

#ifndef OLDC
/* function prototypes */
void   getoptions(void);
void   allocrest(void);
void   doinit(void);
void   makeweights(void);
void   doinput(void);
void   supplement(node *);
void   evaluate(node *);
void   addtraverse(node *, node *, node *, long *, long *, valptr, placeptr);
void   addit(long );
void   dnapenny_reroot(node *);
void   describe(void);
void   maketree(void);
void reallocchars(void);
/* function prototypes */
#endif


Char infilename[FNMLNGTH],outfilename[FNMLNGTH],outtreename[FNMLNGTH], weightfilename[FNMLNGTH];
node *root, *p;
long *zeros=NULL;
long chars, howmanny, howoften, col, msets, ith;
boolean weights, thresh, simple, trout, progress, stepbox, ancseq,
               mulsets, firstset, justwts;
double threshold;
steptr oldweight;
pointptr treenode;   /* pointers to all nodes in tree */
double fracdone, fracinc;
boolean *added;
gbases *garbage;
node **grbg;
Char basechar[]="ACMGRSVTWYHKDBNO???????????????";

/* Variables for maketree, propagated globally for C version: */
long examined, mults;
boolean firsttime, done, recompute;
double like, bestyet;
treenumbers *bestorders, *bestrees;
treenumbers current, order;
long *threshwt;
baseptr nothing;
node *temp, *temp1;
long suppno[] =
  { 0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,4};

long suppset[] =          /* this was previously a function. */
{                                  /* in C, it doesn't need to be.    */
   1 << ((long)A),
   1 << ((long)C),
   1 << ((long)G),
   1 << ((long)T),
   1 << ((long)O),
   (1 << ((long)A)) | (1 << ((long)C)),
   (1 << ((long)A)) | (1 << ((long)G)),
   (1 << ((long)A)) | (1 << ((long)T)),
   (1 << ((long)A)) | (1 << ((long)O)),
   (1 << ((long)C)) | (1 << ((long)G)),
   (1 << ((long)C)) | (1 << ((long)T)),
   (1 << ((long)C)) | (1 << ((long)O)),
   (1 << ((long)G)) | (1 << ((long)T)),
   (1 << ((long)G)) | (1 << ((long)O)),
   (1 << ((long)T)) | (1 << ((long)O)),
   (1 << ((long)A)) | (1 << ((long)C)) | (1 << ((long)G)),
   (1 << ((long)A)) | (1 << ((long)C)) | (1 << ((long)T)),
   (1 << ((long)A)) | (1 << ((long)C)) | (1 << ((long)O)),
   (1 << ((long)A)) | (1 << ((long)G)) | (1 << ((long)T)),
   (1 << ((long)A)) | (1 << ((long)G)) | (1 << ((long)O)),
   (1 << ((long)A)) | (1 << ((long)T)) | (1 << ((long)O)),
   (1 << ((long)C)) | (1 << ((long)G)) | (1 << ((long)T)),
   (1 << ((long)C)) | (1 << ((long)G)) | (1 << ((long)O)),
   (1 << ((long)C)) | (1 << ((long)T)) | (1 << ((long)O)),
   (1 << ((long)G)) | (1 << ((long)T)) | (1 << ((long)O)),
   (1 << ((long)A))|(1 << ((long)C))|(1 << ((long)G))|(1 << ((long)T)),
   (1 << ((long)A))|(1 << ((long)C))|(1 << ((long)G))|(1 << ((long)O)),
   (1 << ((long)A))|(1 << ((long)C))|(1 << ((long)T))|(1 << ((long)O)),
   (1 << ((long)A))|(1 << ((long)G))|(1 << ((long)T))|(1 << ((long)O)),
   (1 << ((long)C))|(1 << ((long)G))|(1 << ((long)T)) | (1 << ((long)O)),
   (1 << ((long)A))|(1 << ((long)C))|(1 << ((long)G)) | (1 << ((long)T)) | (1 << ((long)O))};


void getoptions()
{
  /* interactively set options */
  long loopcount, loopcount2;
  Char ch, ch2;

  fprintf(outfile, "\nPenny algorithm for DNA, version %s\n",VERSION);
  fprintf(outfile, " branch-and-bound to find all");
  fprintf(outfile, " most parsimonious trees\n\n");
  howoften = often;
  howmanny = many;
  outgrno = 1;
  outgropt = false;
  simple = true;
  thresh = false;
  threshold = spp;
  trout = true;
  weights = false;
  justwts = false;
  printdata = false;
  dotdiff = true;
  progress = true;
  treeprint = true;
  stepbox = false;
  ancseq = false;
  interleaved = true;
  loopcount = 0;
  for (;;) {
    cleerhome();
    printf("\nPenny algorithm for DNA, version %s\n",VERSION);
    printf(" branch-and-bound to find all most parsimonious trees\n\n");
    printf("Settings for this run:\n");
    printf("  H        How many groups of %4ld trees:%6ld\n", howoften, 
        howmanny);
    printf("  F        How often to report, in trees:  %4ld\n", howoften);
    printf("  S           Branch and bound is simple?  %s\n",
           (simple ?  "Yes" : "No. reconsiders order of species"));
    printf("  O                        Outgroup root?  %s%3ld\n",
           (outgropt ? "Yes, at sequence number" :
                       "No, use as outgroup species"),outgrno);
    printf("  T              Use Threshold parsimony?");
    if (thresh)
      printf("  Yes, count steps up to%4.1f per site\n", threshold);
    else
      printf("  No, use ordinary parsimony\n");
    printf("  W                       Sites weighted?  %s\n",
           (weights ? "Yes" : "No"));
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld %s\n", msets,
              (justwts ? "sets of weights" : "data sets"));
    else
      printf("  No\n");
    printf("  I          Input sequences interleaved?  %s\n",
           (interleaved ? "Yes" : "No, sequential"));
    printf("  0   Terminal type (IBM PC, ANSI, none)?  %s\n",
           (ibmpc ? "IBM PC" : ansi ? "ANSI" : "(none)"));
    printf("  1    Print out the data at start of run  %s\n",
           (printdata ? "Yes" : "No"));
    printf("  2  Print indications of progress of run  %s\n",
           (progress ? "Yes" : "No"));
    printf("  3                        Print out tree  %s\n",
           (treeprint ? "Yes" : "No"));
    printf("  4          Print out steps in each site  %s\n",
           (stepbox ? "Yes" : "No" ));
    printf("  5  Print sequences at all nodes of tree  %s\n",
           (ancseq ? "Yes" : "No"));
    printf("  6       Write out trees onto tree file?  %s\n",
           (trout ? "Yes" : "No"));
    if(weights && justwts){
      printf("WARNING:  W option and Multiple Weights options are both on.  ");
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
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
    if (ch == 'Y')
      break;
    uppercase(&ch);
    if ((strchr("WHMSOFTI1234560",ch)) != NULL){
      switch (ch) {
        
      case 'H':
        inithowmany(&howmanny, howoften);
        break;

      case 'W':
        weights = !weights;
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
       
      case 'F':
        inithowoften(&howoften);
        break;
        
      case 'S':
        simple = !simple;
        break;

      case 'O':
        outgropt = !outgropt;
        if (outgropt)
          initoutgroup(&outgrno, spp);
        else
           outgrno = 1;
        break;
        
      case 'T':
        thresh = !thresh;
        if (thresh)
          initthreshold(&threshold);
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
    countup(&loopcount, 100);
  }
}  /* getoptions */

void allocrest()
{
  long i;

  y = (Char **)Malloc(spp*sizeof(Char *));
  for (i = 0; i < spp; i++)
    y[i] = (Char *)Malloc(chars*sizeof(Char));
  weight = (long *)Malloc(chars*sizeof(long));
  oldweight = (long *)Malloc(chars*sizeof(long));
  alias = (steptr)Malloc(chars*sizeof(long));
  ally = (steptr)Malloc(chars*sizeof(long));
  location = (steptr)Malloc(chars*sizeof(long));
  nayme = (naym *)Malloc(spp*sizeof(naym));
  bestorders = (treenumbers *)Malloc(maxtrees*sizeof(treenumbers));
  for (i = 1; i <= maxtrees; i++)
    bestorders[i - 1] = (treenumbers)Malloc(spp*sizeof(long));
  bestrees = (treenumbers *)Malloc(maxtrees*sizeof(treenumbers));
  for (i = 1; i <= maxtrees; i++)
    bestrees[i - 1] = (treenumbers)Malloc(spp*sizeof(long));
  current = (treenumbers)Malloc(spp*sizeof(long));
  order = (treenumbers)Malloc(spp*sizeof(long));
  added = (boolean *)Malloc(nonodes*sizeof(boolean));
}  /* allocrest */

void reallocchars(void) 
{/* The amount of chars can change between runs 
    this function reallocates all the variables 
    whose size depends on the amount of chars */
  long i;

  for (i = 0; i < spp; i++) {
    free(y[i]);
    y[i] = (Char *)Malloc(chars*sizeof(Char));
  }
  
  free(weight);
  free(oldweight);
  free(alias);
  free(ally);
  free(location);
  
  weight = (long *)Malloc(chars*sizeof(long));
  oldweight = (long *)Malloc(chars*sizeof(long));
  alias = (steptr)Malloc(chars*sizeof(long));
  ally = (steptr)Malloc(chars*sizeof(long));
  location = (steptr)Malloc(chars*sizeof(long));
} /* reallocchars */

void doinit()
{
  /* initializes variables */
  inputnumbers(&spp, &chars, &nonodes, 1);
  getoptions();
  if (printdata)
    fprintf(outfile, "%2ld species, %3ld  sites\n", spp, chars);
  alloctree(&treenode, nonodes, false);
  allocrest();
}  /* doinit */

void makeweights()
{
  /* make up weights vector to avoid duplicate computations */
  long i;

  for (i = 1; i <= chars; i++) {
    alias[i - 1] = i;
    oldweight[i - 1] = weight[i - 1];
    ally[i - 1] = i;
  }
  sitesort(chars, weight);
  sitecombine(chars);
  sitescrunch(chars);
  endsite = 0;
  for (i = 1; i <= chars; i++) {
    if (ally[i - 1] == i)
      endsite++;
  }
  for (i = 1; i <= endsite; i++)
    location[alias[i - 1] - 1] = i;
  if (!thresh)
    threshold = spp;
  threshwt = (long *)Malloc(endsite*sizeof(long));
  for (i = 0; i < endsite; i++) {
    weight[i] *= 10;
    threshwt[i] = (long)(threshold * weight[i] + 0.5);
  }
  if ( zeros  != NULL )
    free(zeros);
  zeros = (long *)Malloc(endsite*sizeof(long)); /*in makeweights()*/
  for (i = 0; i < endsite; i++)
    zeros[i] = 0;
}  /* makeweights */


void doinput()
{
  /* reads the input data */
  long i;

  if (justwts) {
    if (firstset)
      inputdata(chars);
    for (i = 0; i < chars; i++)
      weight[i] = 1;
    inputweights(chars, weight, &weights);
    if (justwts) {
      fprintf(outfile, "\n\nWeights set # %ld:\n\n", ith);
      if (progress)
        printf("\nWeights set # %ld:\n\n", ith);
    }
    if (printdata)
      printweights(outfile, 0, chars, weight, "Sites");
  } else {
    if (!firstset){
      samenumsp(&chars, ith);
      reallocchars();
    }
    inputdata(chars);
    for (i = 0; i < chars; i++)
      weight[i] = 1;
    if (weights) {
      inputweights(chars, weight, &weights);
      if (printdata)
        printweights(outfile, 0, chars, weight, "Sites");
    }
  }

  makeweights();
  makevalues(treenode, zeros, false);
  alloctemp(&temp, zeros, endsite);
  alloctemp(&temp1, zeros, endsite);
}  /* doinput */


void supplement(node *r)
{
  /* determine minimum number of steps more which will
     be added when rest of species are put in tree */
  long i, j, k, sum, sumall=0, sumadded=0;
  boolean doneadded, allhave, addedhave, has;
  long supps;

  for (i = 0; i < endsite; i++) {
    sum = 3;
    j = 1;
    doneadded = false;
    do {
      allhave = true;
      addedhave = true;
      supps = suppset[j-1];
      for (k = 0; k < spp; k++) {
        has = ((treenode[k]->base[i] & supps) != 0);
        if (added[k] && !doneadded)
          addedhave = (addedhave && has);
        allhave = (allhave && has);
      }
      if (allhave)
        sumall = suppno[j - 1];
      if (addedhave)
        sumadded = suppno[j - 1];
      doneadded = (doneadded || addedhave);
      j++;
    } while (!(j > 31 || (allhave && doneadded)));
    if (addedhave && allhave)
      sum = sumall - sumadded;
    r->numsteps[i] += sum * weight[i];
  }
}  /* supplement */


void evaluate(node *r)
{
  /* determines the number of steps needed for a tree. this is
     the minimum number of steps needed to evolve sequences on
     this tree */
  long i, steps;
  double sum;

  sum = 0.0;
  supplement(r);
  for (i = 0; i < endsite; i++) {
    steps = r->numsteps[i];
    if ((long)steps <= threshwt[i])
      sum += steps;
    else
      sum += threshwt[i];
  }
  if (examined == 0 && mults == 0)
    bestyet = -1.0;
  like = sum;
}  /* evaluate */


void addtraverse(node *p, node *item, node *fork, long *m, long *n,
                        valptr valyew, placeptr place)
{
  /* traverse all places to add item */
  if (done)
    return;
  if (*m <= 2 || (p != root && p != root->next->back)) {
    if (p == root)
      fillin(temp, item, p);
    else {
      fillin(temp1, item, p);
      fillin(temp, temp1, p->back);
    }
    (*n)++;
    evaluate(temp);
    examined++;
    if (examined == howoften) {
      examined = 0;
      mults++;
      if (mults == howmanny)
        done = true;
      if (progress) {
        printf("%7ld", mults);
        if (bestyet >= 0)
          printf("%16.1f", bestyet / 10.0);
        else
          printf("            -   ");
        printf("%17ld%20.2f\n", nextree - 1, fracdone * 100);
#ifdef WIN32
        phyFillScreenColor();
#endif
      }
    }
    valyew[(*n) - 1] = like;
    place[(*n) - 1] = p->index;
  }
  if (!p->tip) {
    addtraverse(p->next->back, item, fork, m,n,valyew,place);
    addtraverse(p->next->next->back, item, fork,m,n,valyew,place);
  }
}  /* addtraverse */


void addit(long m)
{
  /* adds the species one by one, recursively */
  long n;
  valptr valyew;
  placeptr place;
  long i, j, n1, besttoadd=0;
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
    addtraverse(root, treenode[order[m - 1] - 1],
                treenode[spp + m - 2], &m,&n,valyew,place);
    besttoadd = order[m - 1];
    memcpy(bestplace, place, nonodes*sizeof(long));
    memcpy(bestval, valyew, nonodes*sizeof(double));
  } else {
    bestsum = -1.0;
    for (i = 1; i <= spp; i++) {
      if (!added[i - 1]) {
        n = 0;
        added[i - 1] = true;
        addtraverse(root, treenode[i - 1], treenode[spp + m - 2],
                    &m,&n,valyew,place);
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
  for (i = 0; i < n; i++) {
    if (valyew[i] <= bestyet || bestyet < 0.0)
      n1++;
  }
  if (n1 > 0)
    fracinc /= n1;
  for (i = 0; i < n; i++) {
    if (valyew[i] <= bestyet || bestyet < 0.0) {
      current[m - 1] = place[i];
      recompute = (m < spp);
      add(treenode[place[i] - 1], treenode[besttoadd - 1],
          treenode[spp + m - 2], &root, recompute, treenode, grbg, zeros);
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
      recompute = (m < spp);
      re_move(treenode[besttoadd - 1], &treenode[spp + m - 2], &root,
                recompute, treenode, grbg, zeros);
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


void dnapenny_reroot(node *outgroup)
{
  /* reorients tree, putting outgroup in desired position. */
  node *p, *q, *newbottom, *oldbottom;

  if (outgroup->back->index == root->index)
    return;
  newbottom = outgroup->back;
  p = treenode[newbottom->index - 1]->back;
  while (p->index != root->index) {
    oldbottom = treenode[p->index - 1];
    treenode[p->index - 1] = p;
    p = oldbottom->back;
  }
  p = root->next;
  q = root->next->next;
  p->back->back = q->back;
  q->back->back = p->back;
  p->back = outgroup;
  q->back = outgroup->back;
  outgroup->back->back = root->next->next;
  outgroup->back = root->next;
  treenode[newbottom->index - 1] = newbottom;
}  /* dnapenny_reroot */


void describe()
{
  /* prints ancestors, steps and table of numbers of steps in
     each site */

  if (stepbox)
    writesteps(chars, weights, oldweight, root);
  if (ancseq) {
    hypstates(chars, root, treenode, &garbage, basechar);
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

  if (progress) {
    printf("\nHow many\n");
    printf("trees looked                                       Approximate\n");
    printf("at so far      Length of        How many           percentage\n");
    printf("(multiples     shortest tree    trees this short   searched\n");
    printf("of %4ld):      found so far     found so far       so far\n",
           howoften);
    printf("----------     ------------     ------------       ------------\n");
  }
#ifdef WIN32
  phyFillScreenColor();
#endif
  done = false;
  mults = 0;
  examined = 0;
  nextree = 1;
  root = treenode[0];
  firsttime = true;
  for (i = 0; i < spp; i++)
    added[i] = false;
  added[0] = true;
  order[0] = 1;
  k = 2;
  fracdone = 0.0;
  fracinc = 1.0;
  bestyet = -1.0;
  recompute = true;
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
    fprintf(outfile, "\nrequires a total of %18.3f\n\n", bestyet / 10.0);
    if (nextree == 2)
      fprintf(outfile, "One most parsimonious tree found:\n");
    else
      fprintf(outfile, "%6ld trees in all found\n", nextree - 1);
  }
  if (nextree > maxtrees + 1) {
    if (treeprint)
      fprintf(outfile, "here are the first%4ld of them\n", (long)maxtrees);
    nextree = maxtrees + 1;
  }
  if (treeprint)
    putc('\n', outfile);
  for (i = 0; i < spp; i++)
    added[i] = true;
  for (i = 0; i <= nextree - 2; i++) {
    root = treenode[0];
    for (j = k; j <= spp; j++)
      add(treenode[bestrees[i][j - 1] - 1],
          treenode[bestorders[i][j - 1] - 1], treenode[spp + j - 2],
          &root, recompute, treenode, grbg, zeros);
    dnapenny_reroot(treenode[outgrno - 1]);
    postorder(root);
    evaluate(root);
    printree(root, 1.0);
    describe();
    for (j = k - 1; j < spp; j++)
      re_move(treenode[bestorders[i][j] - 1], &dummy, &root,
                recompute, treenode, grbg, zeros);
  }
  if (progress) {
    printf("\nOutput written to file \"%s\"\n\n", outfilename);
    if (trout)
      printf("Trees also written onto file \"%s\"\n\n", outtreename);
  }
  freetemp(&temp);
  freetemp(&temp1);
  if (ancseq)
    freegarbage(&garbage);
}  /* maketree */


int main(int argc, Char *argv[])
{  /* Penny's branch-and-bound method for DNA sequences */
#ifdef MAC
   argc = 1;                /* macsetup("Dnapenny","");        */
   argv[0] = "Dnapenny";
#endif
  init(argc, argv);

  /* Reads in the number of species, number of characters,
     options and data.  Then finds all most parsimonious trees */
  openfile(&infile,INFILE,"input file", "r",argv[0],infilename);
  openfile(&outfile,OUTFILE,"output file", "w",argv[0],outfilename);

  ibmpc = IBMCRT;
  ansi = ANSICRT;
  mulsets = false;
  garbage = NULL;
  msets = 1;
  firstset = true;
  doinit();
  if (weights || justwts)
    openfile(&weightfile,WEIGHTFILE,"weights file","r",argv[0],weightfilename);
  if (trout)
    openfile(&outtree,OUTTREE,"output tree file", "w",argv[0],outtreename);
  for (ith = 1; ith <= msets; ith++) {
    doinput();
    if (ith == 1)
      firstset = false;
    if (msets > 1 && !justwts) {
      fprintf(outfile, "\nData set # %ld:\n",ith);
      if (progress)
        printf("\nData set # %ld:\n",ith);
    }
    maketree();
    free(threshwt);
    freenodes(nonodes,treenode);
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
}  /* Penny's branch-and-bound method for DNA sequences */
