
#include "phylip.h"
#include "seq.h"

/* version 3.6. (c) Copyright 1993-2004 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define maxtrees        100   /* maximum number of tied trees stored     */

typedef boolean *boolptr;

#ifndef OLDC
/* function prototypes */
void   getoptions(void);
void   allocrest(void);
void   doinit(void);
void   initdnacompnode(node **, node **, node *, long, long, long *,
               long *, initops, pointarray, pointarray, Char *, Char *, FILE *);
void   makeweights(void);
void   doinput(void);
void   mincomp(long );
void   evaluate(node *);
void   localsavetree(void);

void   tryadd(node *, node *, node *);
void   addpreorder(node *, node *, node *);
void   tryrearr(node *, boolean *);
void   repreorder(node *, boolean *);
void   rearrange(node **);
void   describe(void);
void   initboolnames(node *, boolean *);
void   maketree(void);
void   freerest(void);
void   standev3(long, long, long, double, double *, long **, longer);
void   reallocchars(void);
/* function prototypes */
#endif


Char infilename[FNMLNGTH], outfilename[FNMLNGTH], intreename[FNMLNGTH], outtreename[FNMLNGTH], weightfilename[FNMLNGTH];
node *root, *p;
long chars, col, ith, njumble, jumb, msets;
long inseed, inseed0;
boolean jumble, usertree, trout, weights,
               progress, stepbox, ancseq, firstset, mulsets, justwts;
steptr oldweight, necsteps;
pointarray treenode;   /* pointers to all nodes in tree */
long *enterorder;
Char basechar[32]="ACMGRSVTWYHKDBNO???????????????";
bestelm *bestrees;
boolean dummy;
longer seed;
gbases *garbage;
Char ch;
Char progname[20];
long *zeros;

/* Local variables for maketree, propogated globally for C version: */
  long maxwhich;
  double like, maxsteps, bestyet, bestlike, bstlike2;
  boolean lastrearr, recompute;
  double nsteps[maxuser];
  long **fsteps;
  node *there;
  long *place;
  boolptr in_tree;
  baseptr nothing;
  node *temp, *temp1;
  node *grbg;

void getoptions()
{
  /* interactively set options */
  long loopcount, loopcount2;
  Char ch, ch2;

  fprintf(outfile, "\nDNA compatibility algorithm, version %s\n\n",VERSION);
  putchar('\n');
  jumble = false;
  njumble = 1;
  outgrno = 1;
  outgropt = false;
  trout = true;
  usertree = false;
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
    printf("\nDNA compatibility algorithm, version %s\n\n",VERSION);
    printf("Settings for this run:\n");
    printf("  U                 Search for best tree?  %s\n",
           (usertree ? "No, use user trees in input file" : "Yes"));
    if (!usertree) {
      printf("  J   Randomize input order of sequences?");
      if (jumble) {
        printf(
         "  Yes (seed =%8ld,%3ld times)\n", inseed0, njumble);
      }
      else
        printf("  No. Use input order\n");
    }
    printf("  O                        Outgroup root?");
    if (outgropt)
      printf("  Yes, at sequence number%3ld\n", outgrno);
    else
      printf("  No, use as outgroup species%3ld\n", outgrno);
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
           ibmpc ? "IBM PC" : ansi  ? "ANSI"   : "(none)");
    printf("  1    Print out the data at start of run  %s\n",
           (printdata ? "Yes" : "No"));
    printf("  2  Print indications of progress of run  %s\n",
           (progress ? "Yes" : "No"));
    printf("  3                        Print out tree  %s\n",
           (treeprint ? "Yes" : "No"));
    printf("  4  Print steps & compatibility at sites  %s\n",
           (stepbox ? "Yes" : "No"));
    printf("  5  Print sequences at all nodes of tree  %s\n",
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
    if (ch == 'Y')
      break;
    if (((!usertree) && (strchr("WJOTUMI1234560", ch) != NULL))
        || (usertree && ((strchr("WOTUMI1234560", ch) != NULL)))){
      switch (ch) {

      case 'J':
        jumble = !jumble;
        if (jumble)
          initjumble(&inseed, &inseed0, seed, &njumble);
        else njumble = 1;
        break;

      case 'O':
        outgropt = !outgropt;
        if (outgropt)
          initoutgroup(&outgrno, spp);
        break;

      case 'U':
        usertree = !usertree;
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
          if (!jumble) {
            jumble = true;
            initjumble(&inseed, &inseed0, seed, &njumble);
          }
        }
        break;

      case 'I':
        interleaved = !interleaved;
        break;

      case 'W':
        weights = !weights;
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
  free(enterorder);
  free(necsteps);
  free(alias);
  free(ally);
  free(location);
  free(in_tree);

  weight = (steptr)Malloc(chars*sizeof(long));
  oldweight = (steptr)Malloc(chars*sizeof(long));
  enterorder = (long *)Malloc(spp*sizeof(long));
  necsteps = (steptr)Malloc(chars*sizeof(long));
  alias = (steptr)Malloc(chars*sizeof(long));
  ally = (steptr)Malloc(chars*sizeof(long));
  location = (steptr)Malloc(chars*sizeof(long));
  in_tree = (boolptr)Malloc(chars*sizeof(boolean));
}


void allocrest()
{
  long i;

  y = (Char **)Malloc(spp*sizeof(Char *));
  for (i = 0; i < spp; i++)
    y[i] = (Char *)Malloc(chars*sizeof(Char));
  bestrees = (bestelm *) Malloc(maxtrees*sizeof(bestelm));
  for (i = 1; i <= maxtrees; i++)
    bestrees[i - 1].btree = (long *)Malloc(nonodes*sizeof(long));
  nayme = (naym *)Malloc(spp*sizeof(naym));
  weight = (steptr)Malloc(chars*sizeof(long));
  oldweight = (steptr)Malloc(chars*sizeof(long));
  enterorder = (long *)Malloc(spp*sizeof(long));
  necsteps = (steptr)Malloc(chars*sizeof(long));
  alias = (steptr)Malloc(chars*sizeof(long));
  ally = (steptr)Malloc(chars*sizeof(long));
  location = (steptr)Malloc(chars*sizeof(long));
  place = (long *)Malloc((2*spp-1)*sizeof(long));
  in_tree = (boolptr)Malloc(spp*sizeof(boolean));
}  /* allocrest */


void doinit()
{
  /* initializes variables */

  inputnumbers(&spp, &chars, &nonodes, 1);
  getoptions();
  if (printdata)
    fprintf(outfile, "%2ld species, %3ld  sites\n", spp, chars);
  alloctree(&treenode, nonodes, usertree);
  allocrest();
}  /* doinit */


void initdnacompnode(node **p, node **grbg, node *q, long len, long nodei,
                     long *ntips, long *parens, initops whichinit,
                     pointarray treenode, pointarray nodep, Char *str, Char *ch,
                     FILE *intree)
{
  /* initializes a node */
  boolean minusread;
  double valyew, divisor;


  switch (whichinit) {
  case bottom:
    gnutreenode(grbg, p, nodei, endsite, zeros);
    treenode[nodei - 1] = *p;
    break;
  case nonbottom:
    gnutreenode(grbg, p, nodei, endsite, zeros);
    break;
  case tip:
    match_names_to_data (str, treenode, p, spp);
    break;
  case length:
    processlength(&valyew, &divisor, ch, &minusread, intree, parens);
    /* process and discard lengths */
  default:
    break;
  }
} /* initdnacompnode */


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
  zeros = (long *)Malloc(endsite*sizeof(long));
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
  makevalues(treenode, zeros, usertree);
  allocnode(&temp, zeros, endsite);
  allocnode(&temp1, zeros, endsite);
}  /* doinput */


void mincomp(long n)
{
  /* computes for each site the minimum number of steps
    necessary to accomodate those species already
    in the analysis, adding in species n */
  long i, j, k, l, m;
  bases b;
  long s;
  boolean allowable, deleted;

  in_tree[n - 1] = true;
  for (i = 0; i < endsite; i++)
    necsteps[i] = 3;
  for (m = 0; m <= 31; m++) {
    s = 0;
    l = -1;
    k = m;
    for (b = A; (long)b <= (long)O; b = (bases)((long)b + 1)) {
      if ((k & 1) == 1) {
        s |= 1L << ((long)b);
        l++;
      }
      k /= 2;
    }
    for (j = 0; j < endsite; j++) {
      allowable = true;
      i = 1;
      while (allowable && i <= spp) {
        if (in_tree[i - 1] && treenode[i - 1]->base[j] != 0) {
          if ((treenode[i - 1]->base[j] & s) == 0)
            allowable = false;
        }
        i++;
      }
      if (allowable) {
        if (l < necsteps[j])
          necsteps[j] = l;
      }
    }
  }
  for (j = 0; j < endsite; j++) {
    deleted = false;
    for (i = 0; i < spp; i++) {
      if (in_tree[i] && treenode[i]->base[j] == 0)
        deleted = true;
    }
    if (deleted)
      necsteps[j]++;
  }
  for (i = 0; i < endsite; i++)
    necsteps[i] *= weight[i];
}  /* mincomp */


void evaluate(node *r)
{
  /* determines the number of steps needed for a tree. this is
    the minimum number of steps needed to evolve sequences on
    this tree */
  long i, term;
  double sum;

   sum = 0.0;
   for (i = 0; i < endsite; i++) {
    if (r->numsteps[i] == necsteps[i])
      term = weight[i];
    else
      term = 0;
    sum += term;
    if (usertree && which <= maxuser)
      fsteps[which - 1][i] = term;
  }
  if (usertree && which <= maxuser) {
    nsteps[which - 1] = sum;
    if (which == 1) {
      maxwhich = 1;
      maxsteps = sum;
    } else if (sum > maxsteps) {
      maxwhich = which;
      maxsteps = sum;
    }
  }
  like = sum;
}  /* evaluate */


void localsavetree()
{
  /* record in place where each species has to be
    added to reconstruct this tree */
  long i, j;
  node *p;
  boolean done;

  reroot(treenode[outgrno - 1], root);
  savetraverse(root);
  for (i = 0; i < nonodes; i++)
    place[i] = 0;
  place[root->index - 1] = 1;
  for (i = 1; i <= spp; i++) {
    p = treenode[i - 1];
    while (place[p->index - 1] == 0) {
      place[p->index - 1] = i;
      while (!p->bottom)
        p = p->next;
      p = p->back;
    }
    if (i > 1) {
      place[i - 1] = place[p->index - 1];
      j = place[p->index - 1];
      done = false;
      while (!done) {
        place[p->index - 1] = spp + i - 1;
        while (!p->bottom)
          p = p->next;
        p = p->back;
        done = (p == NULL);
        if (!done)
          done = (place[p->index - 1] != j);
      }
    }
  }
}  /* localsavetree */


void tryadd(node *p, node *item, node *nufork)
{
  /* temporarily adds one fork and one tip to the tree.
    if the location where they are added yields greater
    "likelihood" than other locations tested up to that
    time, then keeps that location as there */
  long pos;
  boolean found;
  node *rute, *q;

  if (p == root)
    fillin(temp, item, p);
  else {
    fillin(temp1, item, p);
    fillin(temp, temp1, p->back);
  }
  evaluate(temp);
  if (lastrearr) {
    if (like < bestlike) {
      if (item == nufork->next->next->back) {
        q = nufork->next;
        nufork->next = nufork->next->next;
        nufork->next->next = q;
        q->next = nufork;
      }
    } else if (like >= bstlike2) {
      recompute = false;
      add(p, item, nufork, &root, recompute, treenode, &grbg, zeros);
      rute = root->next->back;
      localsavetree();
      reroot(rute, root);
      if (like > bstlike2) {
        bestlike = bstlike2 = like;
        pos = 1;
        nextree = 1;
        addtree(pos, &nextree, dummy, place, bestrees);
      } else {
        pos = 0;
        findtree(&found, &pos, nextree, place, bestrees);
        if (!found) {
          if (nextree <= maxtrees)
            addtree(pos, &nextree, dummy, place, bestrees);
        }
      }
      re_move(item, &nufork, &root, recompute, treenode, &grbg, zeros);
      recompute = true;
    }
  }
  if (like > bestyet) {
    bestyet = like;
    there = p;
  }
}  /* tryadd */


void addpreorder(node *p, node *item, node *nufork)
{
  /* traverses a binary tree, calling PROCEDURE tryadd
    at a node before calling tryadd at its descendants */

  if (p == NULL)
    return;
  tryadd(p, item, nufork);
  if (!p->tip) {
    addpreorder(p->next->back, item, nufork);
    addpreorder(p->next->next->back, item, nufork);
  }
}  /* addpreorder */


void tryrearr(node *p, boolean *success)
{
  /* evaluates one rearrangement of the tree.
    if the new tree has greater "likelihood" than the old
    one sets success := TRUE and keeps the new tree.
    otherwise, restores the old tree */
  node *frombelow, *whereto, *forknode, *q;
  double oldlike;

  if (p->back == NULL)
    return;
  forknode = treenode[p->back->index - 1];
  if (forknode->back == NULL)
    return;
  oldlike = bestyet;
  if (p->back->next->next == forknode)
    frombelow = forknode->next->next->back;
  else
    frombelow = forknode->next->back;
  whereto = treenode[forknode->back->index - 1];
  if (whereto->next->back == forknode)
    q = whereto->next->next->back;
  else
    q = whereto->next->back;
  fillin(temp1, frombelow, q);
  fillin(temp, temp1, p);
  fillin(temp1, temp, whereto->back);
  evaluate(temp1);
  if (like <= oldlike + LIKE_EPSILON) {
    if (p != forknode->next->next->back)
      return;
    q = forknode->next;
    forknode->next = forknode->next->next;
    forknode->next->next = q;
    q->next = forknode;
    return;
  }
  recompute = false;
  re_move(p, &forknode, &root, recompute, treenode, &grbg, zeros);
  fillin(whereto, whereto->next->back, whereto->next->next->back);
  recompute = true;
  add(whereto, p, forknode, &root, recompute, treenode, &grbg, zeros);
  *success = true;
  bestyet = like;
}  /* tryrearr */


void repreorder(node *p, boolean *success)
{
  /* traverses a binary tree, calling PROCEDURE tryrearr
    at a node before calling tryrearr at its descendants */
  if (p == NULL)
    return;
  tryrearr(p,success);
  if (!p->tip) {
    repreorder(p->next->back,success);
    repreorder(p->next->next->back,success);
  }
}  /* repreorder */


void rearrange(node **r)
{
  /* traverses the tree (preorder), finding any local
    rearrangement which decreases the number of steps.
    if traversal succeeds in increasing the tree's
    "likelihood", PROCEDURE rearrange runs traversal again */
  boolean success=true;

  while (success) {
    success = false;
    repreorder(*r,&success);
  }
}  /* rearrange */


void describe()
{
  /* prints ancestors, steps and table of numbers of steps in
    each site and table of compatibilities */
  long i, j, k;

  if (treeprint) {
    fprintf(outfile, "\ntotal number of compatible sites is ");
    fprintf(outfile, "%10.1f\n", like);
  }
  if (stepbox) {
    writesteps(chars, weights, oldweight, root);
    fprintf(outfile,
            "\n compatibility (Y or N) of each site with this tree:\n\n");
    fprintf(outfile, "      ");
    for (i = 0; i <= 9; i++)
      fprintf(outfile, "%ld", i);
    fprintf(outfile, "\n     *----------\n");
    for (i = 0; i <= (chars / 10); i++) {
      putc(' ', outfile);
      fprintf(outfile, "%3ld !", i * 10);
      for (j = 0; j <= 9; j++) {
        k = i * 10 + j;
        if (k > 0 && k <= chars) {
          if (root->numsteps[location[ally[k - 1] - 1] - 1] ==
              necsteps[location[ally[k - 1] - 1] - 1]) {
            if (oldweight[k - 1] > 0)
              putc('Y', outfile);
            else
              putc('y', outfile);
          } else {
            if (oldweight[k - 1] > 0)
              putc('N', outfile);
            else
              putc('n', outfile);
          }
        } else
          putc(' ', outfile);
      }
      putc('\n', outfile);
    }
  }
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


void initboolnames(node *p, boolean *names)
{
  /* sets BOOLEANs that indicate tips */
  node *q;

  if (p->tip) {
    names[p->index - 1] = true;
    return;
  }
  q = p->next;
  while (q != p) {
    initboolnames(q->back, names);
    q = q->next;
  }
}  /* initboolnames */


void standev3(long chars, long numtrees, long maxwhich, double maxsteps,
              double *nsteps, long **fsteps, longer seed)
{  /* do paired sites test (KHT or SH) on user-defined trees */
  long i, j, k;
  double wt, sumw, sum, sum2, sd;
  double temp;
  double **covar, *P, *f, *r;

#define SAMPLES 1000
  if (numtrees == 2) {
    fprintf(outfile, "Kishino-Hasegawa-Templeton test\n\n");
    fprintf(outfile, "Tree    Compatible  Difference  Its S.D.");
    fprintf(outfile, "   Significantly worse?\n\n");
    which = 1;
    while (which <= numtrees) {
      fprintf(outfile, "%3ld  %11.1f", which, nsteps[which - 1]);
      if (maxwhich == which)
        fprintf(outfile, "  <------ best\n");
      else {
        sumw = 0.0;
        sum = 0.0;
        sum2 = 0.0;
        for (i = 0; i < chars; i++) {
          if (weight[i] > 0) {
            wt = weight[i];
            sumw += wt;
            temp = (fsteps[maxwhich - 1][i] - fsteps[which - 1][i]);
            sum += temp;
            sum2 += temp * temp / wt;
          }
        }
        sd = sqrt(sumw / (sumw - 1.0) * (sum2 - sum * sum /sumw));
        fprintf(outfile, " %10.1f %11.4f",
                (maxsteps-nsteps[which - 1]), sd);
        if (sum > 1.95996 * sd)
          fprintf(outfile, "           Yes\n");
        else
          fprintf(outfile, "           No\n");
      }
      which++;
    }
    fprintf(outfile, "\n\n");
  } else {           /* Shimodaira-Hasegawa test using normal approximation */
    if(numtrees > MAXSHIMOTREES){
      fprintf(outfile, "Shimodaira-Hasegawa test on first %d of %ld trees\n\n"
              , MAXSHIMOTREES, numtrees);
      numtrees = MAXSHIMOTREES;
    } else {
      fprintf(outfile, "Shimodaira-Hasegawa test\n\n");
    }
    covar = (double **)Malloc(numtrees*sizeof(double *));  
    sumw = 0.0;
    for (i = 0; i < chars; i++)
      sumw += weight[i];
    for (i = 0; i < numtrees; i++)
      covar[i] = (double *)Malloc(numtrees*sizeof(double));  
    for (i = 0; i < numtrees; i++) {        /* compute covariances of trees */
      sum = nsteps[i]/sumw;
      for (j = 0; j <=i; j++) {
        sum2 = nsteps[j]/sumw;
        temp = 0.0;
        for (k = 0; k < chars; k++) {
          if (weight[k] > 0) {
            wt = weight[k];
            temp = temp + (fsteps[i][k]- wt * sum)
                          * (fsteps[j][k]- wt * sum2) / wt;
          }
        }
        covar[i][j] = temp;
        if (i != j)
          covar[j][i] = temp;
      }
    }
    for (i = 0; i < numtrees; i++) { /* in-place Cholesky decomposition
                                        of trees x trees covariance matrix */
      sum = 0.0;
      for (j = 0; j <= i-1; j++)
        sum = sum + covar[i][j] * covar[i][j];
      if (sqrt(covar[i][i] <= sum))
        temp = 0.0;
      else
        temp = sqrt(covar[i][i] - sum);
      covar[i][i] = temp;
      for (j = i+1; j < numtrees; j++) {
        sum = 0.0;
        for (k = 0; k < i; k++)
          sum = sum + covar[i][k] * covar[j][k];
        if (fabs(temp) < 1.0E-12)
          covar[j][i] = 0.0;
        else
          covar[j][i] = (covar[j][i] - sum)/temp;
      }
    }
    f = (double *)Malloc(numtrees*sizeof(double)); /* resamples sums */
    P = (double *)Malloc(numtrees*sizeof(double)); /* vector of P's of trees */
    r = (double *)Malloc(numtrees*sizeof(double)); /* store normal variates */
    for (i = 0; i < numtrees; i++)
      P[i] = 0.0;
    sum2 = nsteps[0];             /* sum2 will be largest # of compat. sites */
    for (i = 1; i < numtrees; i++)
      if (sum2 < nsteps[i])
        sum2 = nsteps[i];
    for (i = 1; i <= SAMPLES; i++) {           /* loop over resampled trees */
      for (j = 0; j < numtrees; j++)          /* draw Normal variates */
        r[j] = normrand(seed);
      for (j = 0; j < numtrees; j++) {        /* compute vectors */
        sum = 0.0;
        for (k = 0; k <= j; k++)
          sum += covar[j][k]*r[k];
        f[j] = sum;
      }
      sum = f[1];
      for (j = 1; j < numtrees; j++)          /* get max of vector */
        if (f[j] > sum)
          sum = f[j];
      for (j = 0; j < numtrees; j++)          /* accumulate P's */
        if (sum2-nsteps[j] < sum-f[j])
          P[j] += 1.0/SAMPLES;
    }
    fprintf(outfile, "Tree   Compatible  Difference   P value");
    fprintf(outfile, "   Significantly worse?\n\n");
    for (i = 0; i < numtrees; i++) {
      fprintf(outfile, "%3ld  %10.1f", i+1, nsteps[i]);
      if ((maxwhich-1) == i)
        fprintf(outfile, "  <------ best\n");
      else {
        fprintf(outfile, " %10.1f  %10.3f", sum2-nsteps[i], P[i]);
        if (P[i] < 0.05)
          fprintf(outfile, "           Yes\n");
        else
          fprintf(outfile, "           No\n");
      }
    }
  fprintf(outfile, "\n");
  free(P);             /* free the variables we Malloc'ed */
  free(f);
  free(r);
  for (i = 0; i < numtrees; i++)
    free(covar[i]);
  free(covar);
  }
}  /* standev */


void maketree()
{
  /* constructs a binary tree from the pointers in treenode.
    adds each node at location which yields highest "likelihood"
    then rearranges the tree for greatest "likelihood" */
  long i, j, numtrees, nextnode;
  boolean firsttree, goteof, haslengths;
  double gotlike;
  node *item, *nufork, *dummy;
  pointarray nodep;
  boolean *names;

  if (!usertree) {
    recompute = true;
    for (i = 0; i < spp; i++)
      in_tree[i] = false;
    for (i = 1; i <= spp; i++)
      enterorder[i - 1] = i;
    if (jumble)
      randumize(seed, enterorder);
    root = treenode[enterorder[0] - 1];
    add(treenode[enterorder[0] - 1], treenode[enterorder[1] - 1],
        treenode[spp], &root, recompute, treenode, &grbg, zeros);
    if (progress) {
      printf("Adding species:\n");
      writename(0, 2, enterorder);
#ifdef WIN32
      phyFillScreenColor();
#endif
    }
    in_tree[0] = true;
    in_tree[1] = true;
    lastrearr = false;
    for (i = 3; i <= spp; i++) {
      mincomp(i);
      bestyet = -350.0 * spp * chars;
      item = treenode[enterorder[i - 1] - 1];
      nufork = treenode[spp + i - 2];
      there = root;
      addpreorder(root, item, nufork);
      add(there, item, nufork, &root, recompute, treenode, &grbg, zeros);
      like = bestyet;
      rearrange(&root);
      if (progress) {
        writename(i - 1, 1, enterorder);
#ifdef WIN32
        phyFillScreenColor();
#endif
      }
      lastrearr = (i == spp);
      if (lastrearr) {
        if (progress) {
          printf("\nDoing global rearrangements\n");
          printf("  !");
          for (j = 1; j <= nonodes; j++)
            if ( j % (( nonodes / 72 ) + 1 ) == 0 )
              putchar('-');
          printf("!\n");
#ifdef WIN32
          phyFillScreenColor();
#endif
        }
        bestlike = bestyet;
        if (jumb == 1) {
          bstlike2 = bestlike;
          nextree = 1;
        }
        do {
          if (progress)
            printf("   ");
          gotlike = bestlike;
          for (j = 0; j < nonodes; j++) {
            bestyet = -10.0 * spp * chars;
            item = treenode[j];
            there = root;
            if (item != root) {
              re_move(item, &nufork, &root, recompute, treenode, &grbg, zeros);
              there = root;
              addpreorder(root, item, nufork);
              add(there, item, nufork, &root, recompute, treenode, &grbg, zeros);
            }
            if (progress) {
              if ( j % (( nonodes / 72 ) + 1 ) == 0 )
                putchar('.');
              fflush(stdout);
            }
          }
          if (progress)
            putchar('\n');
        } while (bestlike > gotlike);
      }
    }
    if (progress)
      putchar('\n');
    for (i = spp - 1; i >= 1; i--)
      re_move(treenode[i], &dummy, &root, recompute, treenode, &grbg, zeros);
    if (jumb == njumble) {
      if (treeprint) {
        putc('\n', outfile);
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
      recompute = false;
      for (i = 0; i <= (nextree - 2); i++) {
        root = treenode[0];
        add(treenode[0], treenode[1], treenode[spp], &root, recompute,
              treenode, &grbg, zeros);
        for (j = 3; j <= spp; j++)
          add(treenode[bestrees[i].btree[j - 1] - 1], treenode[j - 1],
            treenode[spp + j - 2], &root, recompute, treenode, &grbg, zeros);
        reroot(treenode[outgrno - 1], root);
        postorder(root);
        evaluate(root);
        printree(root, 1.0);
        describe();
        for (j = 1; j < spp; j++)
          re_move(treenode[j], &dummy, &root, recompute, treenode, &grbg, zeros);
      }
    }
  } else {
    /* Open in binary: ftell() is broken for UNIX line-endings under WIN32 */
    openfile(&intree, INTREE, "input tree file", "rb", progname, intreename);
    numtrees = countsemic(&intree);
    if (numtrees > 2)
      initseed(&inseed, &inseed0, seed);
    if (treeprint) {
      fprintf(outfile, "User-defined tree");
      if (numtrees > 1)
        putc('s', outfile);
      fprintf(outfile, ":\n");
    }
    fsteps = (long **)Malloc(maxuser*sizeof(long *));
    for (j = 1; j <= maxuser; j++)
      fsteps[j - 1] = (long *)Malloc(endsite*sizeof(long));
    names = (boolean *)Malloc(spp*sizeof(boolean));
    nodep = NULL;
    maxsteps = 0.0;
    which = 1;
    while (which <= numtrees) {
      firsttree = true;
      nextnode = 0;
      haslengths = true;
      treeread(intree, &root, treenode, &goteof, &firsttree, nodep, &nextnode,
          &haslengths, &grbg, initdnacompnode,false,nonodes);
      for (j = 0; j < spp; j++)
        names[j] = false;
      initboolnames(root, names);
      for (j = 0; j < spp; j++)
        in_tree[j] = names[j];
      j = 1;
      while (!in_tree[j - 1])
        j++;
      mincomp(j);
      if (outgropt)
        reroot(treenode[outgrno - 1], root);
      postorder(root);
      evaluate(root);
      printree(root, 1.0);
      describe();
      which++;
    }
    FClose(intree);
    putc('\n', outfile);
    if (numtrees > 1 && chars > 1 ) {
      standev3(chars, numtrees, maxwhich, maxsteps, nsteps, fsteps, seed);
    }
    for (j = 1; j <= maxuser; j++)
      free(fsteps[j - 1]);
    free(fsteps);
    free(names);
  }
  if (jumb == njumble) {
    if (progress) {
      printf("\nOutput written to file \"%s\"\n\n", outfilename);
      if (trout)
        printf("Trees also written onto file \"%s\"\n\n", outtreename);
    }
  }
}  /* maketree */


void freerest()
{
  if (!usertree) {
    freenode(&temp);
    freenode(&temp1);
  }
  freegrbg(&grbg);
  if (ancseq)
    freegarbage(&garbage);
  free(zeros);
  freenodes(nonodes, treenode);
}  /*  freerest */


int main(int argc, Char *argv[])
{  /* DNA compatibility by uphill search */
  /* reads in spp, chars, and the data. Then calls maketree to
    construct the tree */
#ifdef MAC
  argc = 1;                /* macsetup("Dnacomp","");        */
  argv[0]="Dnacomp";
#endif
  init(argc, argv);
  openfile(&infile,INFILE,"input file", "r",argv[0],infilename);
  openfile(&outfile,OUTFILE,"output file", "w",argv[0],outfilename);
  mulsets = false;
  garbage = NULL;
  grbg = NULL;
  ibmpc = IBMCRT;
  ansi = ANSICRT;
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
      fprintf(outfile, "Data set # %ld:\n\n", ith);
      if (progress)
        printf("Data set # %ld:\n\n", ith);
    }
    for (jumb = 1; jumb <= njumble; jumb++)
      maketree();
    freerest();
  }
  freetree(nonodes, treenode);
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
  exxit(0);
  return 0;
}  /* DNA compatibility by uphill search */


