/* version 3.6. (c) Copyright 1993-2004 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#include <float.h>

#include "phylip.h"
#include "dist.h"

#define epsilonk         0.000001   /* a very small but not too small number */

#ifndef OLDC
/* function prototypes */
void   getoptions(void);
void   doinit(void);
void   inputoptions(void);
void   getinput(void);
void   input_data(void);
void   add(node *, node *, node *);
void   re_move(node **, node **);
void   scrunchtraverse(node *, node **, double *);
void   combine(node *, node *);
void   scrunch(node *);

void   secondtraverse(node *, node *, node *, node *, long, long,
                long , double *);
void   firstraverse(node *, node *, double *);
void   sumtraverse(node *, double *);
void   evaluate(node *);
void   tryadd(node *, node **, node **);
void   addpreorder(node *, node *, node *);
void   tryrearr(node *, node **, boolean *);
void   repreorder(node *, node **, boolean *);
void   rearrange(node **);

void   dtraverse(node *);
void   describe(void);
void   copynode(node *, node *);
void   copy_(tree *, tree *);
void   maketree(void);
/* function prototypes */
#endif



Char infilename[FNMLNGTH], outfilename[FNMLNGTH], intreename[FNMLNGTH], outtreename[FNMLNGTH];
long nonodes, numtrees, col, datasets, ith, njumble, jumb;
/*   numtrees is used by usertree option part of maketree */
long inseed;
tree curtree, bestree;   /* pointers to all nodes in tree */
boolean minev, jumble, usertree, lower, upper, negallowed, replicates, trout,
        printdata, progress, treeprint, mulsets, firstset;
longer seed;
double power;
long *enterorder;
/* Local variables for maketree, propagated globally for C version: */
  long examined;
  double like, bestyet;
  node *there;
  boolean *names;
  Char ch;
  char *progname;
double trweight; /* to make treeread happy */
boolean goteof, haslengths, lengths;  /* ditto ... */


void getoptions()
{
  /* interactively set options */
  long inseed0, loopcount;
  Char ch;

  minev = false;
  jumble = false;
  njumble = 1;
  lower = false;
  negallowed = false;
  power = 2.0;
  replicates = false;
  upper = false;
  usertree = false;
  trout = true;
  printdata = false;
  progress = true;
  treeprint = true;
  loopcount = 0;
  for(;;) {
    cleerhome();
    printf("\nFitch-Margoliash method ");
    printf("with contemporary tips, version %s\n\n",VERSION);
    printf("Settings for this run:\n");
    printf("  D      Method (F-M, Minimum Evolution)?  %s\n",
           (minev ? "Minimum Evolution" : "Fitch-Margoliash"));
    printf("  U                 Search for best tree?  %s\n",
           usertree ? "No, use user trees in input file" : "Yes");
    printf("  P                                Power?%9.5f\n",power);
    printf("  -      Negative branch lengths allowed?  %s\n",
           (negallowed ? "Yes" : "No"));
    printf("  L         Lower-triangular data matrix?  %s\n",
           (lower ? "Yes" : "No"));
    printf("  R         Upper-triangular data matrix?  %s\n",
           (upper ? "Yes" : "No"));
    printf("  S                        Subreplicates?  %s\n",
           (replicates ? "Yes" : "No"));
    if (!usertree) {
      printf("  J     Randomize input order of species?");
      if (jumble)
            printf("  Yes (seed =%8ld,%3ld times)\n", inseed0, njumble);
      else
        printf("  No. Use input order\n");
    }
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld sets\n", datasets);
    else
      printf("  No\n");
    printf("  0   Terminal type (IBM PC, ANSI, none)?  %s\n",
           (ibmpc ? "IBM PC" : ansi  ? "ANSI" : "(none)"));
    printf("  1    Print out the data at start of run  %s\n",
           (printdata ? "Yes" : "No"));
    printf("  2  Print indications of progress of run  %s\n",
           (progress ? "Yes" : "No"));
    printf("  3                        Print out tree  %s\n",
    (treeprint ? "Yes" : "No"));
    printf("  4       Write out trees onto tree file?  %s\n",
          (trout ? "Yes" : "No"));
    printf("\n  Y to accept these or type the letter for one to change\n");
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
    if (((!usertree) && (strchr("DJUP-LRSM12340", ch) != NULL))
        || (usertree && ((strchr("DUP-LRSM12340", ch) != NULL)))){
      switch (ch) {

      case 'D':
        minev = !minev;
        if (!negallowed)
          negallowed = true;
        break;

      case '-':
        negallowed = !negallowed;
        break;

      case 'J':
        jumble = !jumble;
        if (jumble)
          initjumble(&inseed, &inseed0, seed, &njumble);
        else njumble = 1;
        break;

      case 'L':
        lower = !lower;
        break;

      case 'P':
        initpower(&power);
        break;

      case 'R':
        upper = !upper;
        break;

      case 'S':
        replicates = !replicates;
        break;

      case 'U':
        usertree = !usertree;
        break;

      case 'M':
        mulsets = !mulsets;
        if (mulsets)
          initdatasets(&datasets);
        jumble = true;
        if (jumble)
          initseed(&inseed, &inseed0, seed);
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
  if (upper && lower) {
    printf("ERROR: Data matrix cannot be both uppeR and Lower triangular\n");
    exxit(-1);
  }
}  /* getoptions */


void doinit()
{
  /* initializes variables */

  inputnumbers2(&spp, &nonodes, 1);
  getoptions();
  alloctree(&curtree.nodep, nonodes);
  allocd(nonodes, curtree.nodep);
  allocw(nonodes, curtree.nodep);
  if (!usertree && njumble > 1) {
    alloctree(&bestree.nodep, nonodes);
    allocd(nonodes, bestree.nodep);
    allocw(nonodes, bestree.nodep);
  }
  nayme = (naym *)Malloc(spp*sizeof(naym));
  enterorder = (long *)Malloc(spp*sizeof(long));
}  /* doinit */


void inputoptions()
{
  /* print options information */
  if (!firstset)
    samenumsp2(ith);
  fprintf(outfile, "\nFitch-Margoliash method ");
  fprintf(outfile, "with contemporary tips, version %s\n\n",VERSION);
  if (minev)
    fprintf(outfile, "Minimum evolution method option\n\n");
  fprintf(outfile, "                  __ __             2\n");
  fprintf(outfile, "                  \\  \\   (Obs - Exp)\n");
  fprintf(outfile, "Sum of squares =  /_ /_  ------------\n");
  fprintf(outfile, "                               ");
  if (power == (long)power)
    fprintf(outfile, "%2ld\n", (long)power);
  else
    fprintf(outfile, "%4.1f\n", power);
  fprintf(outfile, "                   i  j      Obs\n\n");
  fprintf(outfile, "negative branch lengths");
  if (!negallowed)
    fprintf(outfile, " not");
  fprintf(outfile, " allowed\n\n");
}  /* inputoptions */


void getinput()
{
  /* reads the input data */
  inputoptions();
}  /* getinput */


void input_data()
{
  /* read in distance matrix */
  long i, j, k, columns, n;
  boolean skipit, skipother;
  double x;
  columns = replicates ? 4 : 6;
  if (printdata) {
    fprintf(outfile, "\nName                       Distances");
    if (replicates)
      fprintf(outfile, " (replicates)");
    fprintf(outfile, "\n----                       ---------");
    if (replicates)
      fprintf(outfile, "-------------");
    fprintf(outfile, "\n\n");
  }
  setuptree(&curtree, nonodes);
  if (!usertree && njumble > 1)
    setuptree(&bestree, nonodes);
  for (i = 0; i < (spp); i++) {
    curtree.nodep[i]->d[i] = 0.0;
    curtree.nodep[i]->w[i] = 0.0;
    curtree.nodep[i]->weight = 0.0;
    scan_eoln(infile);
    initname(i);
    for (j = 1; j <= (spp); j++) {
      skipit = ((lower && j >= i + 1) || (upper && j <= i + 1));
      skipother = ((lower && i + 1 >= j) || (upper && i + 1 <= j));
      if (!skipit) {
        if (eoln(infile))
          scan_eoln(infile);
        fscanf(infile, "%lf", &x);
        curtree.nodep[i]->d[j - 1] = x;
        if (replicates) {
          if (eoln(infile)) 
            scan_eoln(infile);
          fscanf(infile, "%ld", &n);
        } else
          n = 1;
        if (n > 0 && x < 0) {
          printf("NEGATIVE DISTANCE BETWEEN SPECIES%5ld AND %5ld\n",
                 i + 1, j);
          exxit(-1);
        }
        curtree.nodep[i]->w[j - 1] = n;
        if (skipother) {
          curtree.nodep[j - 1]->d[i] = curtree.nodep[i]->d[j - 1];
          curtree.nodep[j - 1]->w[i] = curtree.nodep[i]->w[j - 1];
        }
        if ((i == j) && (fabs(curtree.nodep[i-1]->d[j-1]) > 0.000000001)) {
       printf("\nERROR: diagonal element of row %ld of distance matrix ", i+2);
          printf("is not zero.\n");
          printf("       Is it a distance matrix?\n\n");
          exxit(-1);        
        }
        if ((j < i) && (fabs(curtree.nodep[i]->d[j-1]-curtree.nodep[j-1]->d[i])
             > 0.000000001)) {
          printf("ERROR: distance matrix is not symmetric:\n");
          printf("       (%ld,%ld) element and (%ld,%ld) element are unequal.\n",
            i+1, j+1, j+1, i+1);
          printf("       They are %10.6f and %10.6f, respectively.\n",
                  curtree.nodep[i]->d[j-1], curtree.nodep[j]->d[i-1]);
          printf("       Is it a distance matrix?\n\n");
          exxit(-1);
        }
      }
    }
  }
  scan_eoln(infile);
  if (printdata) {
    for (i = 0; i < (spp); i++) {
      for (j = 0; j < nmlngth; j++)
        putc(nayme[i][j], outfile);
      putc(' ', outfile);
      for (j = 1; j <= (spp); j++) {
        fprintf(outfile, "%10.5f", curtree.nodep[i]->d[j - 1]);
        if (replicates)
          fprintf(outfile, " (%3ld)", (long)curtree.nodep[i]->w[j - 1]);
        if (j % columns == 0 && j < spp) {
          putc('\n', outfile);
          for (k = 1; k <= nmlngth + 1; k++)
            putc(' ', outfile);
        }
      }
      putc('\n', outfile);
    }
    putc('\n', outfile);
  }
  for (i = 0; i < (spp); i++) {
    for (j = 0; j < (spp); j++) {
      if (i + 1 != j + 1) {
        if (curtree.nodep[i]->d[j] < epsilonk)
          curtree.nodep[i]->d[j] = epsilonk;
        curtree.nodep[i]->w[j] /= exp(power * log(curtree.nodep[i]->d[j]));
      }
    }
  }
}  /* inputdata */


void add(node *below, node *newtip, node *newfork)
{
  /* inserts the nodes newfork and its left descendant, newtip,
     to the tree.  below becomes newfork's right descendant */
  if (below != curtree.nodep[below->index - 1])
    below = curtree.nodep[below->index - 1];
  if (below->back != NULL)
    below->back->back = newfork;
  newfork->back = below->back;
  below->back = newfork->next->next;
  newfork->next->next->back = below;
  newfork->next->back = newtip;
  newtip->back = newfork->next;
  if (curtree.root == below)
    curtree.root = newfork;
  curtree.root->back = NULL;
}  /* add */


void re_move(node **item, node **fork)
{
  /* removes nodes item and its ancestor, fork, from the tree.
     the new descendant of fork's ancestor is made to be
     fork's second descendant (other than item).  Also
     returns pointers to the deleted nodes, item and fork */
  node *p, *q;

  if ((*item)->back == NULL) {
    *fork = NULL;
    return;
  }
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
}  /* remove */


void scrunchtraverse(node *u, node **closest, double *tmax)
{
  /* traverse to find closest node to the current one */
  if (!u->sametime) {
    if (u->t > *tmax) {
      *closest = u;
      *tmax = u->t;
    }
    return;
  }
  u->t = curtree.nodep[u->back->index - 1]->t;
  if (!u->tip) {
    scrunchtraverse(u->next->back, closest,tmax);
    scrunchtraverse(u->next->next->back, closest,tmax);
  }
}  /* scrunchtraverse */


void combine(node *a, node *b)
{
  /* put node b into the set having the same time as a */
  if (a->weight + b->weight <= 0.0)
    a->t = 0.0;
  else
    a->t = (a->t * a->weight + b->t * b->weight) / (a->weight + b->weight);
  a->weight += b->weight;
  b->sametime = true;
}  /* combine */


void scrunch(node *s)
{
  /* see if nodes can be combined to prevent negative lengths */
  double tmax;
  node *closest;
  boolean found;

  closest = NULL;
  tmax = -1.0;
  do {
    if (!s->tip) {
      scrunchtraverse(s->next->back, &closest,&tmax);
      scrunchtraverse(s->next->next->back, &closest,&tmax);
    }
    found = (tmax > s->t);
    if (found)
      combine(s, closest);
    tmax = -1.0;
  } while (found);
}  /* scrunch */


void secondtraverse(node *a, node *q, node *u, node *v, long i, long j,
                        long k, double *sum)
{
  /* recalculate distances, add to sum */
  long l;
  double wil, wjl, wkl, wli, wlj, wlk, TEMP;

  if (!(a->processed || a->tip)) {
    secondtraverse(a->next->back, q,u,v,i,j,k,sum);
    secondtraverse(a->next->next->back, q,u,v,i,j,k,sum);
    return;
  }
  if (!(a != q && a->processed))
    return;
  l = a->index;
  wil = u->w[l - 1];
  wjl = v->w[l - 1];
  wkl = wil + wjl;
  wli = a->w[i - 1];
  wlj = a->w[j - 1];
  wlk = wli + wlj;
  q->w[l - 1] = wkl;
  a->w[k - 1] = wlk;
  if (wkl <= 0.0)
    q->d[l - 1] = 0.0;
  else
    q->d[l - 1] = (wil * u->d[l - 1] + wjl * v->d[l - 1]) / wkl;
  if (wlk <= 0.0)
    a->d[k - 1] = 0.0;
  else
    a->d[k - 1] = (wli * a->d[i - 1] + wlj * a->d[j - 1]) / wlk;
  if (minev)
    return;
  if (wkl > 0.0) {
    TEMP = u->d[l - 1] - v->d[l - 1];
    (*sum) += wil * wjl / wkl * (TEMP * TEMP);
  }
  if (wlk > 0.0) {
    TEMP = a->d[i - 1] - a->d[j - 1];
    (*sum) += wli * wlj / wlk * (TEMP * TEMP);
  }
}  /* secondtraverse */


void firstraverse(node *q_, node *r, double *sum)
{  /* firsttraverse                              */
   /* go through tree calculating branch lengths */
  node *q;
  long i, j, k;
  node *u, *v;

  q = q_;
  if (q == NULL)
    return;
  q->sametime = false;
  if (!q->tip) {
    firstraverse(q->next->back, r,sum);
    firstraverse(q->next->next->back, r,sum);
  }
  q->processed = true;
  if (q->tip)
    return;
  u = q->next->back;
  v = q->next->next->back;
  i = u->index;
  j = v->index;
  k = q->index;
  if (u->w[j - 1] + v->w[i - 1] <= 0.0)
    q->t = 0.0;
  else
    q->t = (u->w[j - 1] * u->d[j - 1] +  v->w[i - 1] * v->d[i - 1]) /
             (2.0 * (u->w[j - 1] + v->w[i - 1]));
  q->weight = u->weight + v->weight + u->w[j - 1] + v->w[i - 1];
  if (!negallowed)
    scrunch(q);
  u->v = q->t - u->t;
  v->v = q->t - v->t;
  u->back->v = u->v;
  v->back->v = v->v;
  secondtraverse(r,q,u,v,i,j,k,sum);
}  /* firstraverse */


void sumtraverse(node *q, double *sum)
{
  /* traverse to finish computation of sum of squares */
  long i, j;
  node *u, *v;
  double TEMP, TEMP1;

  if (minev && (q != curtree.root))
    *sum += q->v;
  if (q->tip)
    return;
  sumtraverse(q->next->back, sum);
  sumtraverse(q->next->next->back, sum);
  if (!minev) {
    u = q->next->back;
    v = q->next->next->back;
    i = u->index;
    j = v->index;
    TEMP = u->d[j - 1] - 2.0 * q->t;
    TEMP1 = v->d[i - 1] - 2.0 * q->t;
    (*sum) += u->w[j - 1] * (TEMP * TEMP) + v->w[i - 1] * (TEMP1 * TEMP1);
  }
}  /* sumtraverse */


void evaluate(node *r)
{
  /* fill in times and evaluate sum of squares for tree */
  double sum;
  long i;
  sum = 0.0;
  for (i = 0; i < (nonodes); i++)
    curtree.nodep[i]->processed = curtree.nodep[i]->tip;
  firstraverse(r, r,&sum);
  sumtraverse(r, &sum);
  examined++;
  if (replicates && (lower || upper))
    sum /= 2;
  like = -sum;
}  /* evaluate */


void tryadd(node *p, node **item, node **nufork)
{
  /* temporarily adds one fork and one tip to the tree.
     if the location where they are added yields greater
     "likelihood" than other locations tested up to that
     time, then keeps that location as there */
  add(p, *item, *nufork);
  evaluate(curtree.root);
  if (like > bestyet) {
    bestyet = like;
    there = p;
  }
  re_move(item, nufork);
}  /* tryadd */


void addpreorder(node *p, node *item, node *nufork)
{
  /* traverses a binary tree, calling PROCEDURE tryadd
     at a node before calling tryadd at its descendants */
  if (p == NULL)
    return;
  tryadd(p, &item,&nufork);
  if (!p->tip) {
    addpreorder(p->next->back, item, nufork);
    addpreorder(p->next->next->back, item, nufork);
  }
}  /* addpreorder */


void tryrearr(node *p, node **r, boolean *success)
{
  /* evaluates one rearrangement of the tree.
     if the new tree has greater "likelihood" than the old
     one sets success := TRUE and keeps the new tree.
     otherwise, restores the old tree */
  node *frombelow, *whereto, *forknode;
  double oldlike;

  if (p->back == NULL)
    return;
  forknode = curtree.nodep[p->back->index - 1];
  if (forknode->back == NULL)
    return;
  oldlike = like;
  if (p->back->next->next == forknode)
    frombelow = forknode->next->next->back;
  else
    frombelow = forknode->next->back;
  whereto = forknode->back;
  re_move(&p, &forknode);
  add(whereto, p, forknode);
  if ((*r)->back != NULL)
    *r = curtree.nodep[(*r)->back->index - 1];
  evaluate(*r);
  if (like - oldlike > LIKE_EPSILON) {
    bestyet = like;
    *success = true;
    return;
  }
  re_move(&p, &forknode);
  add(frombelow, p, forknode);
  if ((*r)->back != NULL)
    *r = curtree.nodep[(*r)->back->index - 1];
  like = oldlike;
}  /* tryrearr */


void repreorder(node *p, node **r, boolean *success)
{
  /* traverses a binary tree, calling PROCEDURE tryrearr
     at a node before calling tryrearr at its descendants */
  if (p == NULL)
    return;
  tryrearr(p,r,success);
  if (!p->tip) {
    repreorder(p->next->back,r,success);
    repreorder(p->next->next->back,r,success);
  }
}  /* repreorder */


void rearrange(node **r_)
{
  /* traverses the tree (preorder), finding any local
     rearrangement which decreases the number of steps.
     if traversal succeeds in increasing the tree's
     "likelihood", PROCEDURE rearrange runs traversal again */
  node **r;
  boolean success;

  r = r_;
  success = true;
  while (success) {
    success = false;
    repreorder(*r,r,&success);
  }
}  /* rearrange */


void dtraverse(node *q)
{
  /* print table of lengths etc. */
  long i;

  if (!q->tip)
    dtraverse(q->next->back);
  if (q->back != NULL) {
    fprintf(outfile, "%4ld   ", q->back->index - spp);
    if (q->index <= spp) {
      for (i = 0; i < nmlngth; i++)
        putc(nayme[q->index - 1][i], outfile);
    } else
      fprintf(outfile, "%4ld      ", q->index - spp);
    fprintf(outfile, "%13.5f", curtree.nodep[q->back->index - 1]->t - q->t);
    q->v = curtree.nodep[q->back->index - 1]->t - q->t;
    q->back->v = q->v;
    fprintf(outfile, "%16.5f\n", curtree.root->t - q->t);
  }
  if (!q->tip)
    dtraverse(q->next->next->back);
}  /* dtraverse */


void describe()
{
  /* prints table of lengths, times, sum of squares, etc. */
  long i, j;
  double totalnum;
  double TEMP;

  if (!minev)
    fprintf(outfile, "\nSum of squares = %10.3f\n\n", -like);
  else
    fprintf(outfile, "Sum of branch lengths = %10.3f\n\n", -like);
  if ((fabs(power - 2) < 0.01) && !minev) {
    totalnum = 0.0;
    for (i = 0; i < (spp); i++) {
      for (j = 0; j < (spp); j++) {
        if (i + 1 != j + 1 && curtree.nodep[i]->d[j] > 0.0) {
          TEMP = curtree.nodep[i]->d[j];
          totalnum += curtree.nodep[i]->w[j] * (TEMP * TEMP);
        }
      }
    }
    totalnum -= 2;
    if (replicates && (lower || upper))
      totalnum /= 2;
    fprintf(outfile, "Average percent standard deviation =");
    fprintf(outfile, "%10.5f\n\n", 100 * sqrt(-(like / totalnum)));
  }
  fprintf(outfile, "From     To            Length          Height\n");
  fprintf(outfile, "----     --            ------          ------\n\n");
  dtraverse(curtree.root);
  putc('\n', outfile);
  if (trout) {
    col = 0;
    treeoutr(curtree.root,&col,&curtree);
  }
}  /* describe */


void copynode(node *c, node *d)
{
  /* make a copy of a node */

  memcpy(d->d, c->d, nonodes*sizeof(double));
  memcpy(d->w, c->w, nonodes*sizeof(double));
  d->t = c->t;
  d->sametime = c->sametime;
  d->weight = c->weight;
  d->processed = c->processed;
  d->xcoord = c->xcoord;
  d->ycoord = c->ycoord;
  d->ymin = c->ymin;
  d->ymax = c->ymax;
}  /* copynode */


void copy_(tree *a, tree *b)
{
  /* make a copy of a tree */
  long i, j=0;
  node *p, *q;

  for (i = 0; i < spp; i++) {
    copynode(a->nodep[i], b->nodep[i]);
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
  for (i = spp; i < nonodes; i++) {
    p = a->nodep[i];
    q = b->nodep[i];
    for (j = 1; j <= 3; j++) {
      copynode(p, q);
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
  b->root = a->root;
}  /* copy */


void maketree()
{
  /* constructs a binary tree from the pointers in curtree.nodep.
     adds each node at location which yields highest "likelihood"
     then rearranges the tree for greatest "likelihood" */
  long i, j, which;
  double bestlike, bstlike2=0, gotlike;
  boolean lastrearr;
  node *item, *nufork;

  if (!usertree) {
    if (jumb == 1) {
      input_data();
      examined = 0;
    }
    for (i = 1; i <= (spp); i++)
      enterorder[i - 1] = i;
    if (jumble)
      randumize(seed, enterorder);
    curtree.root = curtree.nodep[enterorder[0] - 1];
    add(curtree.nodep[enterorder[0] - 1], curtree.nodep[enterorder[1] - 1],
        curtree.nodep[spp]);
    if (progress) {
      printf("Adding species:\n");
      writename(0, 2, enterorder);
#ifdef WIN32
      phyFillScreenColor();
#endif
    }
    for (i = 3; i <= (spp); i++) {
      bestyet = -DBL_MAX;
      item = curtree.nodep[enterorder[i - 1] - 1];
      nufork = curtree.nodep[spp + i - 2];
      addpreorder(curtree.root, item, nufork);
      add(there, item, nufork);
      like = bestyet;
      rearrange(&curtree.root);
      evaluate(curtree.root);
      examined--;
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
          for (j = 1; j <= (nonodes); j++)
            if ( j % (( nonodes / 72 ) + 1 ) == 0 )
              putchar('-');
          printf("!\n");
#ifdef WIN32
          phyFillScreenColor();
#endif
        }
        bestlike = bestyet;
        do {
          gotlike = bestlike;
          if (progress)
            printf("   ");
          for (j = 0; j < (nonodes); j++) {
            there = curtree.root;
            bestyet = -DBL_MAX;
            item = curtree.nodep[j];
            if (item != curtree.root) {
              re_move(&item, &nufork);
              there = curtree.root;
              addpreorder(curtree.root, item, nufork);
              add(there, item, nufork);
            }
            if (progress) {
              if ( j % (( nonodes / 72 ) + 1 ) == 0 )
                putchar('.');
              fflush(stdout);
            }
          }
          if (progress) {
            putchar('\n');
#ifdef WIN32
            phyFillScreenColor();
#endif
          }
        } while (bestlike > gotlike);
        if (njumble > 1) {
          if (jumb == 1 || (jumb > 1 && bestlike > bstlike2)) {
            copy_(&curtree, &bestree);
            bstlike2 = bestlike;
          }
        }
      }
    }
    if (njumble == jumb) {
      if (njumble > 1)
        copy_(&bestree, &curtree);
      evaluate(curtree.root);
      printree(curtree.root, treeprint, false, true);
      describe();
    }
  } else {
    input_data();
    /* Open in binary: ftell() is broken for UNIX line-endings under WIN32 */
    openfile(&intree,INTREE,"input tree file","rb",progname,intreename);
    numtrees = countsemic(&intree);
    if (treeprint)
      fprintf(outfile, "\n\nUser-defined trees:\n\n");
    names = (boolean *)Malloc(spp*sizeof(boolean));
    which = 1;
    while (which <= numtrees ) {
      treeread2 (intree, &curtree.root, curtree.nodep, lengths, &trweight, 
                      &goteof, &haslengths, &spp,false,nonodes);
      if (curtree.root->back) {
        printf("Error:  Kitsch cannot read unrooted user trees\n");
        exxit(-1);
      }
      evaluate(curtree.root);
      printree(curtree.root, treeprint, false, true);
      describe();
      which++;
    }
    FClose(intree);
    free(names);
  }
  if (jumb == njumble && progress) {
    printf("\nOutput written to file \"%s\"\n\n", outfilename);
    if (trout)
      printf("Tree also written onto file \"%s\"\n", outtreename);
    putchar('\n');
  }
}  /* maketree */


int main(int argc, Char *argv[])
{  /* Fitch-Margoliash criterion with contemporary tips */
#ifdef MAC
  argc = 1;                /* macsetup("Kitsch","");        */
  argv[0] = "Kitsch";
#endif
  init(argc,argv);
  /* reads in spp, options, and the data, then calls maketree to
     construct the tree */
  progname = argv[0];
  openfile(&infile,INFILE,"input file","r",argv[0],infilename);
  openfile(&outfile,OUTFILE,"output file","w",argv[0],outfilename);

  ibmpc = IBMCRT;
  ansi = ANSICRT;
  mulsets = false;
  firstset = true;
  datasets = 1;
  doinit();
  openfile(&outtree,OUTTREE,"output tree file","w",argv[0],outtreename);
  for (ith = 1; ith <= datasets; ith++) {
    if (datasets > 1) {
      fprintf(outfile, "\nData set # %ld:\n",ith);
      if (progress)
        printf("\nData set # %ld:\n",ith);
    }
    getinput();
    for (jumb = 1; jumb <= njumble; jumb++)
      maketree();
    firstset = false;
    if (eoln(infile) && (ith < datasets))
      scan_eoln(infile);
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
  printf("Done.\n\n");
  return 0;
}  /* Fitch-Margoliash criterion with contemporary tips */
