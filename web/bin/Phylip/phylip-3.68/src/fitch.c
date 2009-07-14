
#include "phylip.h"
#include "dist.h"
#include "float.h"

/* version 3.6. (c) Copyright 1993-2004 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define zsmoothings     10    /* number of zero-branch correction iterations */
#define epsilonf        0.000001   /* a very small but not too small number  */
#define delta           0.0001      /* a not quite so small number */
#define MAXNUMTREES   100000000 /* a number bigger than conceivable numtrees */


#ifndef OLDC
/* function prototypes */
void   getoptions(void);
void   allocrest(void);
void   doinit(void);
void   inputoptions(void);
void   fitch_getinput(void);
void   secondtraverse(node *, double , long *, double *);
void   firsttraverse(node *, long *, double *);
double evaluate(tree *);
void   nudists(node *, node *);
void   makedists(node *);

void   makebigv(node *);
void   correctv(node *);
void   alter(node *, node *);
void   nuview(node *);
void   update(node *);
void   smooth(node *);
void   filltraverse(node *, node *, boolean);
void   fillin(node *, node *, boolean);
void   insert_(node *, node *, boolean);
void   copynode(node *, node *);

void   copy_(tree *, tree *);
void   setuptipf(long, tree *);
void   buildnewtip(long , tree *, long);
void   buildsimpletree(tree *, long);
void   addtraverse(node *, node *, boolean, long *, boolean *);
void   re_move(node **, node **);
void   rearrange(node *, long *, long *, boolean *);
void   describe(node *);
void   summarize(long);
void   nodeinit(node *);
void   initrav(node *);
void   treevaluate(void);
void   maketree(void);
void   globrearrange(long* numtrees,boolean* succeeded);
/* function prototypes */
#endif



Char infilename[FNMLNGTH], outfilename[FNMLNGTH], intreename[FNMLNGTH], outtreename[FNMLNGTH];
long nonodes2, outgrno, nums, col, datasets, ith, njumble, jumb=0;
long inseed;
vector *x;
intvector *reps;
boolean minev, global, jumble, lengths, usertree, lower, upper, negallowed,
        outgropt, replicates, trout, printdata, progress, treeprint,
        mulsets, firstset;
double power;
double trweight; /* to make treeread happy */
boolean goteof, haslengths;  /* ditto ... */
boolean first; /* ditto ... */
node *addwhere;

longer seed;
long *enterorder;
tree curtree, priortree, bestree, bestree2;
Char ch;
char *progname;



void getoptions()
{
  /* interactively set options */
  long inseed0=0, loopcount;
  Char ch;
  boolean done=false;

  putchar('\n');
  minev = false;
  global = false;
  jumble = false;
  njumble = 1;
  lengths = false;
  lower = false;
  negallowed = false;
  outgrno = 1;
  outgropt = false;
  power = 2.0;
  replicates = false;
  trout = true;
  upper = false;
  usertree = false;
  printdata = false;
  progress = true;
  treeprint = true;
  loopcount = 0;
  do {
    cleerhome();
    printf("\nFitch-Margoliash method version %s\n\n",VERSION);
    printf("Settings for this run:\n");
    printf("  D      Method (F-M, Minimum Evolution)?  %s\n",
             (minev ? "Minimum Evolution" : "Fitch-Margoliash"));
    printf("  U                 Search for best tree?  %s\n",
           (usertree ? "No, use user trees in input file" : "Yes"));
    if (usertree) {
      printf("  N          Use lengths from user trees?  %s\n",
             (lengths ? "Yes" : "No"));
    }
    printf("  P                                Power?%9.5f\n",power);
    printf("  -      Negative branch lengths allowed?  %s\n",
           negallowed ? "Yes" : "No");
    printf("  O                        Outgroup root?");
    if (outgropt)
      printf("  Yes, at species number%3ld\n", outgrno);
    else
      printf("  No, use as outgroup species%3ld\n", outgrno);
    printf("  L         Lower-triangular data matrix?");
    if (lower)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  R         Upper-triangular data matrix?");
    if (upper)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  S                        Subreplicates?");
    if (replicates)
      printf("  Yes\n");
    else
      printf("  No\n");
    if (!usertree) {
      printf("  G                Global rearrangements?");
      if (global)
        printf("  Yes\n");
      else
        printf("  No\n");
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
    printf(
   "\n  Y to accept these or type the letter for one to change\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    fflush(stdout);
    scanf("%c%*[^\n]", &ch);
    getchar();
    uppercase(&ch);
    done = (ch == 'Y');
    if (!done) {
      if (((!usertree) && (strchr("DJOUNPG-LRSM01234", ch) != NULL))
          || (usertree && ((strchr("DOUNPG-LRSM01234", ch) != NULL)))){
        switch (ch) {

        case 'D':
          minev = !minev;
          if (minev && (!negallowed))
            negallowed = true;
          break;

        case '-':
          negallowed = !negallowed;
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
          lower = !lower;
          break;

         case 'N':
           lengths = !lengths;
           break;

        case 'O':
          outgropt = !outgropt;
          if (outgropt)
            initoutgroup(&outgrno, spp);
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
    }
    countup(&loopcount, 100);
  } while (!done);
  if (lower && upper) {
    printf("ERROR: Data matrix cannot be both uppeR and Lower triangular\n");
    exxit(-1);
  }
}  /* getoptions */


void allocrest()
{
  long i;

  x = (vector *)Malloc(spp*sizeof(vector));
  reps = (intvector *)Malloc(spp*sizeof(intvector));
  for (i=0;i<spp;++i){
    x[i]=(vector)Malloc(nonodes2 * sizeof(double));
    reps[i]=(intvector)Malloc(spp * sizeof(long));
  }
  nayme = (naym *)Malloc(spp*sizeof(naym));
  enterorder = (long *)Malloc(spp*sizeof(long));
}


void doinit()
{
  /* initializes variables */

  inputnumbers2(&spp, &nonodes2, 1);
  getoptions();
  if ( !usertree )
    nonodes2--;
  alloctree(&curtree.nodep, nonodes2);
  allocd(nonodes2, curtree.nodep);
  allocw(nonodes2, curtree.nodep);
  if (!usertree) {
    alloctree(&bestree.nodep, nonodes2);
    allocd(nonodes2, bestree.nodep);
    allocw(nonodes2, bestree.nodep);
    alloctree(&priortree.nodep, nonodes2);
    allocd(nonodes2, priortree.nodep);
    allocw(nonodes2, priortree.nodep);
    if (njumble > 1) {
      alloctree(&bestree2.nodep, nonodes2);
      allocd(nonodes2, bestree2.nodep);
      allocw(nonodes2, bestree2.nodep);
    }
  }
  allocrest();
}  /* doinit */


void inputoptions()
{
  /* print options information */
  if (!firstset)
    samenumsp2(ith);
  fprintf(outfile, "\nFitch-Margoliash method version %s\n\n",VERSION);
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
  fprintf(outfile, "Negative branch lengths ");
  if (!negallowed)
    fprintf(outfile, "not ");
  fprintf(outfile, "allowed\n\n");
  if (global)
    fprintf(outfile, "global optimization\n\n");
}  /* inputoptions */


void fitch_getinput()
{
  /* reads the input data */
  inputoptions();
}  /* fitch_getinput */


void secondtraverse(node *q, double y, long *nx, double *sum)
{
  /* from each of those places go back to all others */
   /* nx comes from firsttraverse */
   /* sum comes from evaluate via firsttraverse */
  double z=0.0, TEMP=0.0;

  z = y + q->v;
  if (q->tip) {
    TEMP = q->d[(*nx) - 1] - z;
    *sum += q->w[(*nx) - 1] * (TEMP * TEMP);
  } else {
    secondtraverse(q->next->back, z, nx, sum);
    secondtraverse(q->next->next->back, z, nx,sum);
  }
}  /* secondtraverse */


void firsttraverse(node *p, long *nx, double *sum)
{
  /* go through tree calculating branch lengths */
  if (minev && (p != curtree.start))
    *sum += p->v;
  if (p->tip) {
    if (!minev) {
      *nx = p->index;
      secondtraverse(p->back, 0.0, nx, sum);
      }
  } else {
    firsttraverse(p->next->back, nx,sum);
    firsttraverse(p->next->next->back, nx,sum);
  }
}  /* firsttraverse */


double evaluate(tree *t)
{
  double sum=0.0;
  long nx=0;
  /* evaluate likelihood of a tree */
  firsttraverse(t->start->back ,&nx, &sum);
  firsttraverse(t->start, &nx, &sum);
  if ((!minev) && replicates && (lower || upper))
    sum /= 2;
  t->likelihood = -sum;
  return (-sum);
}  /* evaluate */


void nudists(node *x, node *y)
{
  /* compute distance between an interior node and tips */
  long nq=0, nr=0, nx=0, ny=0;
  double dil=0, djl=0, wil=0, wjl=0, vi=0, vj=0;
  node *qprime, *rprime;

  qprime = x->next;
  rprime = qprime->next->back;
  qprime = qprime->back;
  ny = y->index;
  dil = qprime->d[ny - 1];
  djl = rprime->d[ny - 1];
  wil = qprime->w[ny - 1];
  wjl = rprime->w[ny - 1];
  vi = qprime->v;
  vj = rprime->v;
  x->w[ny - 1] = wil + wjl;
  if (wil + wjl <= 0.0)
    x->d[ny - 1] = 0.0;
  else
    x->d[ny - 1] = ((dil - vi) * wil + (djl - vj) * wjl) / (wil + wjl);
  nx = x->index;
  nq = qprime->index;
  nr = rprime->index;
  dil = y->d[nq - 1];
  djl = y->d[nr - 1];
  wil = y->w[nq - 1];
  wjl = y->w[nr - 1];
  y->w[nx - 1] = wil + wjl;
  if (wil + wjl <= 0.0)
    y->d[nx - 1] = 0.0;
  else
    y->d[nx - 1] = ((dil - vi) * wil + (djl - vj) * wjl) / (wil + wjl);
}  /* nudists */


void makedists(node *p)
{
  /* compute distances among three neighbors of a node */
  long i=0, nr=0, ns=0;
  node *q, *r, *s;

  r = p->back;
  nr = r->index;
  for (i = 1; i <= 3; i++) {
    q = p->next;
    s = q->back;
    ns = s->index;
    if (s->w[nr - 1] + r->w[ns - 1] <= 0.0)
      p->dist = 0.0;
    else
      p->dist = (s->w[nr - 1] * s->d[nr - 1] + r->w[ns - 1] * r->d[ns - 1]) /
                (s->w[nr - 1] + r->w[ns - 1]);
    p = q;
    r = s;
    nr = ns;
  }
}  /* makedists */


void makebigv(node *p)
{
  /* make new branch length */
  long i=0;
  node *temp, *q, *r;

  q = p->next;
  r = q->next;
  for (i = 1; i <= 3; i++) {
    if (p->iter) {
      p->v = (p->dist + r->dist - q->dist) / 2.0;
      p->back->v = p->v;
    }
    temp = p;
    p = q;
    q = r;
    r = temp;
  }
}  /* makebigv */


void correctv(node *p)
{
  /* iterate branch lengths if some are to be zero */
  node *q, *r, *temp;
  long i=0, j=0, n=0, nq=0, nr=0, ntemp=0;
  double wq=0.0, wr=0.0;

  q = p->next;
  r = q->next;
  n = p->back->index;
  nq = q->back->index;
  nr = r->back->index;
  for (i = 1; i <= zsmoothings; i++) {
    for (j = 1; j <= 3; j++) {
      if (p->iter) {
        wr = r->back->w[n - 1] + p->back->w[nr - 1];
        wq = q->back->w[n - 1] + p->back->w[nq - 1];
        if (wr + wq <= 0.0 && !negallowed)
          p->v = 0.0;
        else
          p->v = ((p->dist - q->v) * wq + (r->dist - r->v) * wr) / (wr + wq);
        if (p->v < 0 && !negallowed)
          p->v = 0.0;
        p->back->v = p->v;
      }
      temp = p;
      p = q;
      q = r;
      r = temp;
      ntemp = n;
      n = nq;
      nq = nr;
      nr = ntemp;
    }
  }
}  /* correctv */


void alter(node *x, node *y)
{
  /* traverse updating these views */
  nudists(x, y);
  if (!y->tip) {
    alter(x, y->next->back);
    alter(x, y->next->next->back);
  }
}  /* alter */


void nuview(node *p)
{
  /* renew information about subtrees */
  long i=0;
  node *q, *r, *pprime, *temp;

  q = p->next;
  r = q->next;
  for (i = 1; i <= 3; i++) {
    temp = p;
    pprime = p->back;
    alter(p, pprime);
    p = q;
    q = r;
    r = temp;
  }
}  /* nuview */


void update(node *p)
{
  /* update branch lengths around a node */

  if (p->tip)
    return;
  makedists(p);
  if (p->iter || p->next->iter || p->next->next->iter) {
    makebigv(p);
    correctv(p);
  }
  nuview(p);
}  /* update */


void smooth(node *p)
{
  /* go through tree getting new branch lengths and views */
  if (p->tip)
    return;
  update(p);
  smooth(p->next->back);
  smooth(p->next->next->back);
}  /* smooth */


void filltraverse(node *pb, node *qb, boolean contin)
{
  if (qb->tip)
    return;
  if (contin) {
    filltraverse(pb, qb->next->back,contin);
    filltraverse(pb, qb->next->next->back,contin);
    nudists(qb, pb);
    return;
  }
  if (!qb->next->back->tip)
    nudists(qb->next->back, pb);
  if (!qb->next->next->back->tip)
    nudists(qb->next->next->back, pb);
}  /* filltraverse */


void fillin(node *pa, node *qa, boolean contin)
{
  if (!pa->tip) {
    fillin(pa->next->back, qa, contin);
    fillin(pa->next->next->back, qa, contin);
  }
  filltraverse(pa, qa, contin);
}  /* fillin */


void insert_(node *p, node *q, boolean contin_)
{
  /* put p and q together and iterate info. on resulting tree */
  double x=0.0, oldlike;
  hookup(p->next->next, q->back);
  hookup(p->next, q);
  x = q->v / 2.0;
  p->v = 0.0;
  p->back->v = 0.0;
  p->next->v = x;
  p->next->back->v = x;
  p->next->next->back->v = x;
  p->next->next->v = x;
  fillin(p->back, p, contin_);
  evaluate(&curtree);
  do {
    oldlike = curtree.likelihood;
    smooth(p);
    smooth(p->back);
    evaluate(&curtree);
  } while (fabs(curtree.likelihood - oldlike) > delta);
}  /* insert_ */


void copynode(node *c, node *d)
{
  /* make a copy of a node */

  memcpy(d->d, c->d, nonodes2*sizeof(double));
  memcpy(d->w, c->w, nonodes2*sizeof(double));
  d->v = c->v;
  d->iter = c->iter;
  d->dist = c->dist;
  d->xcoord = c->xcoord;
  d->ycoord = c->ycoord;
  d->ymin = c->ymin;
  d->ymax = c->ymax;
}  /* copynode */


void copy_(tree *a, tree *b)
{
  /* make copy of a tree a to tree b */
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
  for (i = spp; i < nonodes2; i++) {
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
  b->likelihood = a->likelihood;
  b->start = a->start;
}  /* copy_ */


void setuptipf(long m, tree *t)
{
  /* initialize branch lengths and views in a tip */
  long i=0;
  intvector n=(long *)Malloc(spp * sizeof(long)); 
  node *WITH;

  WITH = t->nodep[m - 1];
  memcpy(WITH->d, x[m - 1], (nonodes2 * sizeof(double)));
  memcpy(n, reps[m - 1], (spp * sizeof(long)));
  for (i = 0; i < spp; i++) {
    if (i + 1 != m && n[i] > 0) {
      if (WITH->d[i] < epsilonf)
        WITH->d[i] = epsilonf;
      WITH->w[i] = n[i] / exp(power * log(WITH->d[i]));
    } else {
      WITH->w[i] = 0.0;
      WITH->d[i] = 0.0;
    }
  }
  for (i = spp; i < nonodes2; i++) {
    WITH->w[i] = 1.0;
    WITH->d[i] = 0.0;
  }
  WITH->index = m;
  if (WITH->iter) WITH->v = 0.0;
  free(n);
}  /* setuptipf */


void buildnewtip(long m, tree *t, long nextsp)
{
  /* initialize and hook up a new tip */
  node *p;
  setuptipf(m, t);
  p = t->nodep[nextsp + spp - 3];
  hookup(t->nodep[m - 1], p);
}  /* buildnewtip */


void buildsimpletree(tree *t, long nextsp)
{
  /* make and initialize a three-species tree */
  curtree.start=curtree.nodep[enterorder[0] - 1]; 
  setuptipf(enterorder[0], t);
  setuptipf(enterorder[1], t);
  hookup(t->nodep[enterorder[0] - 1], t->nodep[enterorder[1] - 1]);
  buildnewtip(enterorder[2], t, nextsp);
  insert_(t->nodep[enterorder[2] - 1]->back, t->nodep[enterorder[0] - 1],
          false);
}  /* buildsimpletree */


void addtraverse(node *p, node *q, boolean contin, long *numtrees,
                        boolean *succeeded)
{
 /* traverse through a tree, finding best place to add p */
  insert_(p, q, true);
  (*numtrees)++;
  if (evaluate(&curtree) > (bestree.likelihood + 
                            epsilonf * fabs(bestree.likelihood))){
    copy_(&curtree, &bestree);
    addwhere = q;
    (*succeeded)=true;
  }
  copy_(&priortree, &curtree);
  if (!q->tip && contin) {
    addtraverse(p, q->next->back, contin,numtrees,succeeded);
    addtraverse(p, q->next->next->back, contin,numtrees,succeeded);
  }
}  /* addtraverse */


void re_move(node **p, node **q)
{
  /* re_move p and record in q where it was */
  *q = (*p)->next->back;
  hookup(*q, (*p)->next->next->back);
  (*p)->next->back = NULL;
  (*p)->next->next->back = NULL;
  update(*q);
  update((*q)->back);
}  /* re_move */


void globrearrange(long* numtrees,boolean* succeeded) 
{
  /* does global rearrangements */
  tree globtree;
  tree oldtree;
  int i,j,k,num_sibs,num_sibs2;
  node *where,*sib_ptr,*sib_ptr2;
  double oldbestyet = curtree.likelihood;
  int success = false;
 
  alloctree(&globtree.nodep,nonodes2);
  alloctree(&oldtree.nodep,nonodes2);
  setuptree(&globtree,nonodes2);
  setuptree(&oldtree,nonodes2);
  allocd(nonodes2, globtree.nodep);
  allocd(nonodes2, oldtree.nodep);
  allocw(nonodes2, globtree.nodep);
  allocw(nonodes2, oldtree.nodep);
  copy_(&curtree,&globtree);
  copy_(&curtree,&oldtree);
  for ( i = spp ; i < nonodes2 ; i++ ) {
    num_sibs = count_sibs(curtree.nodep[i]);
    sib_ptr  = curtree.nodep[i];
    if ( (i - spp) % (( nonodes2 / 72 ) + 1 ) == 0 )
      putchar('.');
    fflush(stdout);
    for ( j = 0 ; j <= num_sibs ; j++ ) {
      re_move(&sib_ptr,&where);
      copy_(&curtree,&priortree);
      
      if (where->tip) {
        copy_(&oldtree,&curtree);
        copy_(&oldtree,&bestree);
        sib_ptr=sib_ptr->next;
        continue;
      }
      else num_sibs2 = count_sibs(where);
      sib_ptr2 = where;
      for ( k = 0 ; k < num_sibs2 ; k++ ) {
        addwhere = NULL;
        addtraverse(sib_ptr,sib_ptr2->back,true,numtrees,succeeded);
        if ( addwhere && where != addwhere && where->back != addwhere
              && bestree.likelihood > globtree.likelihood) {
            copy_(&bestree,&globtree);
            success = true;
        }
        sib_ptr2 = sib_ptr2->next;
      } 
      copy_(&oldtree,&curtree);
      copy_(&oldtree,&bestree);
      sib_ptr = sib_ptr->next;
    }
  }
  copy_(&globtree,&curtree);
  copy_(&globtree,&bestree);
  if (success && globtree.likelihood > oldbestyet)  {
    *succeeded = true;
  }
  else  {
    *succeeded = false;
  }
  freed(nonodes2, globtree.nodep);
  freed(nonodes2, oldtree.nodep);
  freew(nonodes2, globtree.nodep);
  freew(nonodes2, oldtree.nodep);
  freetree(&globtree.nodep,nonodes2);
  freetree(&oldtree.nodep,nonodes2);
}
 

void rearrange(node *p, long *numtrees, long *nextsp, boolean *succeeded)
{
  node *q, *r;
  if (!p->tip && !p->back->tip) {
    r = p->next->next;
    re_move(&r, &q);
    copy_(&curtree, &priortree);
    addtraverse(r, q->next->back, false, numtrees,succeeded);
    addtraverse(r, q->next->next->back, false, numtrees,succeeded);
    copy_(&bestree, &curtree);
    if (global && ((*nextsp) == spp)) {
      putchar('.');
      fflush(stdout);
    }
  }
  if (!p->tip) {
    rearrange(p->next->back, numtrees,nextsp,succeeded);
    rearrange(p->next->next->back, numtrees,nextsp,succeeded);
  }
}  /* rearrange */


void describe(node *p)
{
  /* print out information for one branch */
  long i=0;
  node *q;

  q = p->back;
  fprintf(outfile, "%4ld          ", q->index - spp);
  if (p->tip) {
    for (i = 0; i < nmlngth; i++)
      putc(nayme[p->index - 1][i], outfile);
  } else
    fprintf(outfile, "%4ld      ", p->index - spp);
  fprintf(outfile, "%15.5f\n", q->v);
  if (!p->tip) {
    describe(p->next->back);
    describe(p->next->next->back);
  }
}  /* describe */


void summarize(long numtrees)
{
  /* print out branch lengths etc. */
  long i, j, totalnum;

  fprintf(outfile, "\nremember:");
  if (outgropt)
    fprintf(outfile, " (although rooted by outgroup)");
  fprintf(outfile, " this is an unrooted tree!\n\n");
  if (!minev)
    fprintf(outfile, "Sum of squares = %11.5f\n\n", -curtree.likelihood);
  else
    fprintf(outfile, "Sum of branch lengths = %11.5f\n\n", -curtree.likelihood);
  if ((power == 2.0) && !minev) {
    totalnum = 0;
    for (i = 1; i <= nums; i++) {
      for (j = 1; j <= nums; j++) {
        if (i != j)
          totalnum += reps[i - 1][j - 1];
      }
    }
    fprintf(outfile, "Average percent standard deviation = ");
    fprintf(outfile, "%11.5f\n\n",
            100 * sqrt(-curtree.likelihood / (totalnum - 2)));
  }
  fprintf(outfile, "Between        And            Length\n");
  fprintf(outfile, "-------        ---            ------\n");
  describe(curtree.start->next->back);
  describe(curtree.start->next->next->back);
  describe(curtree.start->back);
  fprintf(outfile, "\n\n");
  if (trout) {
    col = 0;
    treeout(curtree.start, &col, 0.43429445222, true,
              curtree.start);
  }
}  /* summarize */


void nodeinit(node *p)
{
  /* initialize a node */
  long i, j;

  for (i = 1; i <= 3; i++) {
    for (j = 0; j < nonodes2; j++) {
      p->w[j] = 1.0;
      p->d[j] = 0.0;
    }
    p = p->next;
  }
  if ((!lengths) || p->iter)
    p->v = 1.0;
  if ((!lengths) || p->back->iter)
    p->back->v = 1.0;
}  /* nodeinit */


void initrav(node *p)
{
  /* traverse to initialize */
  if (p->tip)
    return;
  nodeinit(p);
  initrav(p->next->back);
  initrav(p->next->next->back);
}  /* initrav */

void treevaluate()
{
  /* evaluate user-defined tree, iterating branch lengths */
  long i;
  double oldlike;

  for (i = 1; i <= spp; i++)
    setuptipf(i, &curtree);
  unroot(&curtree,nonodes2);

  initrav(curtree.start);
  if (curtree.start->back != NULL) {
    initrav(curtree.start->back);
    evaluate(&curtree);
    do {
      oldlike = curtree.likelihood;
      smooth(curtree.start);
      evaluate(&curtree);
    } while (fabs(curtree.likelihood - oldlike) > delta);
  }
  evaluate(&curtree);
}  /* treevaluate */


void maketree()
{
  /* contruct the tree */
  long nextsp,numtrees;
  boolean succeeded=false;
  long i, j, which;

  if (usertree) {
    inputdata(replicates, printdata, lower, upper, x, reps);
    setuptree(&curtree, nonodes2);
    for (which = 1; which <= spp; which++)
      setuptipf(which, &curtree);
    if (eoln(infile))
      scan_eoln(infile);
    /* Open in binary: ftell() is broken for UNIX line-endings under WIN32 */
    openfile(&intree,INTREE,"input tree file","rb",progname,intreename);
    numtrees = countsemic(&intree);
    if (numtrees > MAXNUMTREES) {
      printf("\nERROR: number of input trees is read incorrectly from %s\n",
        intreename);
      exxit(-1);
    }
    if (treeprint) {
      fprintf(outfile, "User-defined tree");
      if (numtrees > 1)
        putc('s', outfile);
      fprintf(outfile, ":\n\n");
    }
    first = true;
    which = 1;
    while (which <= numtrees) {
      treeread2 (intree, &curtree.start, curtree.nodep,
        lengths, &trweight, &goteof, &haslengths, &spp,false,nonodes2);
      nums = spp;
      curtree.start = curtree.nodep[outgrno - 1]->back;
      treevaluate();
      printree(curtree.start, treeprint, false, false);
      summarize(numtrees);
      clear_connections(&curtree,nonodes2);
      which++;
    }
    FClose(intree);
  } else {
    if (jumb == 1) {
      inputdata(replicates, printdata, lower, upper, x, reps);
      setuptree(&curtree, nonodes2);
      setuptree(&priortree, nonodes2);
      setuptree(&bestree, nonodes2);
      if (njumble > 1) setuptree(&bestree2, nonodes2);
    }
    for (i = 1; i <= spp; i++)
      enterorder[i - 1] = i;
    if (jumble)
      randumize(seed, enterorder);
    nextsp = 3;
    buildsimpletree(&curtree, nextsp);
    curtree.start = curtree.nodep[enterorder[0] - 1]->back;
    if (jumb == 1) numtrees = 1;
    nextsp = 4;
    if (progress) {
      printf("Adding species:\n");
      writename(0, 3, enterorder);
#ifdef WIN32
      phyFillScreenColor();
#endif
    }
    while (nextsp <= spp) {
      nums = nextsp;
      buildnewtip(enterorder[nextsp - 1], &curtree, nextsp);
      copy_(&curtree, &priortree);
      bestree.likelihood = -DBL_MAX;
      curtree.start = curtree.nodep[enterorder[0] - 1]->back;
      addtraverse(curtree.nodep[enterorder[nextsp - 1] - 1]->back,
                  curtree.start, true, &numtrees,&succeeded);
      copy_(&bestree, &curtree);
      if (progress) {
        writename(nextsp  - 1, 1, enterorder);
#ifdef WIN32
        phyFillScreenColor();
#endif
      }
      if (global && nextsp == spp) {
        if (progress) {
          printf("Doing global rearrangements\n");
          printf("  !");
          for (j = spp; j < nonodes2; j++)
            if ( (j - spp) % (( nonodes2 / 72 ) + 1 ) == 0 )
              putchar('-');
          printf("!\n");
          printf("   ");
        }
      }
      succeeded = true;
      while (succeeded) {
        succeeded = false;
        curtree.start = curtree.nodep[enterorder[0] - 1]->back;
        if (nextsp == spp  && global)
          globrearrange (&numtrees,&succeeded);
        else{
          rearrange(curtree.start,&numtrees,&nextsp,&succeeded);
        }
        if (global && ((nextsp) == spp) && progress)
          printf("\n   ");
      }
      if (global && nextsp == spp) {
        putc('\n', outfile);
        if (progress)
          putchar('\n');
      }
      if (njumble > 1) {
        if (jumb == 1 && nextsp == spp)
          copy_(&bestree, &bestree2);
        else if (nextsp == spp) {
          if (bestree2.likelihood < bestree.likelihood)
            copy_(&bestree, &bestree2);
        }
      }
      if (nextsp == spp && jumb == njumble) {
        if (njumble > 1) copy_(&bestree2, &curtree);
        curtree.start = curtree.nodep[outgrno - 1]->back;
        printree(curtree.start, treeprint, true, false);
        summarize(numtrees);
      }
      nextsp++;
    }
  }
  if (jumb == njumble && progress) {
    printf("\nOutput written to file \"%s\"\n\n", outfilename);
    if (trout) {
      printf("Tree also written onto file \"%s\"\n", outtreename);
      putchar('\n');
    }
  }
}  /* maketree */


int main(int argc, Char *argv[])
{
  int i;
#ifdef MAC
  argc = 1;                /* macsetup("Fitch","");        */
  argv[0]="Fitch";
#endif
  init(argc,argv);
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
  for (i=0;i<spp;++i){
    enterorder[i]=0;}
  for (ith = 1; ith <= datasets; ith++) {
    if (datasets > 1) {
      fprintf(outfile, "Data set # %ld:\n\n",ith);
      if (progress)
        printf("\nData set # %ld:\n\n",ith);
    }
    fitch_getinput();
    for (jumb = 1; jumb <= njumble; jumb++)
        maketree();
    firstset = false;
    if (eoln(infile) && (ith < datasets))
      scan_eoln(infile);
  }
  if (trout)
    FClose(outtree);
  FClose(outfile);
  FClose(infile);
#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(outtreename);
#endif
  printf("Done.\n\n");
#ifdef WIN32
  phyRestoreConsoleAttributes();
#endif
  return 0;
}
