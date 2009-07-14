/* version 3.6. (c) Copyright 1993-2004 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#include <float.h>

#include "phylip.h"
#include "cont.h"


#define epsilon1        0.000001   /* small number */
#define epsilon2        0.02   /* not such a small number */

#define smoothings      4   /* number of passes through smoothing algorithm */

#define over            60


#ifndef OLDC
/* function prototypes */
void   getoptions(void);
void   allocrest(void);
void   doinit(void);
void   getalleles(void);
void   inputdata(void);
void   transformgfs(void);
void   getinput(void);
void   sumlikely(node *, node *, double *);
double evaluate(tree *);
double distance(node *, node *);
void   makedists(node *);

void   makebigv(node *, boolean *);
void   correctv(node *);
void   littlev(node *);
void   nuview(node *);
void   update(node *);
void   smooth(node *);
void   insert_(node *, node *);
void   copynode(node *, node *);
void   copy_(tree *, tree *);
void   inittip(long, tree *);

void   buildnewtip(long, tree *, long);
void   buildsimpletree(tree *);
void   addtraverse(node *, node *, boolean);
void   re_move(node **, node **);
void   rearrange(node *);
void   coordinates(node *, double, long *, double *);
void   drawline(long, double);
void   printree(void);
void   treeout(node *);
void   describe(node *, double, double);

void   summarize(void);
void   nodeinit(node *);
void   initrav(node *);
void   treevaluate(void);
void   maketree(void);
void   globrearrange(void);
/* function prototypes */
#endif


Char infilename[FNMLNGTH], outfilename[FNMLNGTH],
      intreename[FNMLNGTH], outtreename[FNMLNGTH];
long nonodes2, loci, totalleles, df, outgrno, col,
      datasets, ith, njumble, jumb=0;
long inseed, inseed0;
long *alleles, *locus, *weight;
phenotype3 *x;
boolean all, contchars, global, jumble, lengths, outgropt, trout,
               usertree, printdata, progress, treeprint, mulsets, firstset;
longer seed;
long *enterorder;
tree curtree, priortree, bestree, bestree2;
long nextsp, numtrees, which, maxwhich, shimotrees;
 /* From maketree, propagated to global */
boolean succeeded;
double maxlogl;
double l0gl[MAXSHIMOTREES];
double *pbar, *sqrtp, *l0gf[MAXSHIMOTREES];
Char ch;
char *progname;
double trweight;   /* added to make treeread happy */
boolean goteof;
boolean haslengths;   /* end of ones added to make treeread happy */
node *addwhere;


void getoptions()
{ /* interactively set options */
  long inseed0, loopcount;
  Char ch;
  boolean done;

  fprintf(outfile, "\nContinuous character Maximum Likelihood");
  fprintf(outfile, " method version %s\n\n",VERSION);
  putchar('\n');
  global = false;
  jumble = false;
  njumble = 1;
  lengths = false;
  outgrno = 1;
  outgropt = false;
  all = false;
  contchars = false;
  trout = true;
  usertree = false;
  printdata = false;
  progress = true;
  treeprint = true;
  loopcount = 0;
  do {
    cleerhome();
    printf("\nContinuous character Maximum Likelihood");
    printf(" method version %s\n\n",VERSION);
    printf("Settings for this run:\n");
    printf("  U                       Search for best tree?  %s\n",
           (usertree ? "No, use user trees in input" : "Yes"));
    if (usertree) {
      printf("  L                Use lengths from user trees?%s\n",
             (lengths ? "  Yes" : "  No"));
    }
    printf("  C  Gene frequencies or continuous characters?  %s\n",
           (contchars ? "Continuous characters" : "Gene frequencies"));
    if (!contchars)
      printf("  A   Input file has all alleles at each locus?  %s\n",
             (all ? "Yes" : "No, one allele missing at each"));
    printf("  O                              Outgroup root?  %s %ld\n",
           (outgropt ? "Yes, at species number" :
                       "No, use as outgroup species"),outgrno);
    if (!usertree) {
      printf("  G                      Global rearrangements?  %s\n",
             (global ? "Yes" : "No"));
      printf("  J           Randomize input order of species?");
      if (jumble)
        printf("  Yes (seed=%8ld,%3ld times)\n", inseed0, njumble);
      else
        printf("  No. Use input order\n");
    }
    printf("  M                 Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld sets\n", datasets);
    else
      printf("  No\n");
    printf("  0         Terminal type (IBM PC, ANSI, none)?  %s\n",
           ibmpc ? "IBM PC" : ansi  ? "ANSI" : "(none)");
    printf("  1          Print out the data at start of run  %s\n",
           (printdata ? "Yes" : "No"));
    printf("  2        Print indications of progress of run  %s\n",
           (progress ? "Yes" : "No"));
    printf("  3                              Print out tree  %s\n",
           (treeprint ? "Yes" : "No"));
    printf("  4             Write out trees onto tree file?  %s\n",
           (trout ? "Yes" : "No"));
    printf("\n  Y to accept these or type the letter for one to change\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    fflush(stdout);
    scanf("%c%*[^\n]", &ch);
    getchar();
    uppercase(&ch);
    done = (ch == 'Y');
    if (!done) {
      if (((!usertree) && (strchr("JOUGACM12340", ch) != NULL))
          || (usertree && ((strchr("LOUACM12340", ch) != NULL)))){
        switch (ch) {

        case 'A':
          if (!contchars)
            all = !all;
          break;

        case 'C':
          contchars = !contchars;
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
           lengths = !lengths;
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
          if (mulsets)
            initdatasets(&datasets);
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
}  /* getoptions */


void allocrest()
{  /* allocate arrays for number of alleles, the data coordinates, names etc */

  alleles = (long *)Malloc(loci*sizeof(long));
  if (contchars)
    locus = (long *)Malloc(loci*sizeof(long));
  x = (phenotype3 *)Malloc(spp*sizeof(phenotype3));
  nayme = (naym *)Malloc(spp*sizeof(naym));
  enterorder = (long *)Malloc(spp*sizeof(long));
}  /* allocrest */


void doinit()
{ /* initializes variables */

  inputnumbers(&spp, &loci, &nonodes2, 1);
  getoptions();
  if(!usertree)
    nonodes2--;
  if (printdata)
    fprintf(outfile, "\n%4ld Populations, %4ld Loci\n", spp, loci);
  alloctree(&curtree.nodep, nonodes2);
  if (!usertree) {
    alloctree(&bestree.nodep, nonodes2);
    alloctree(&priortree.nodep, nonodes2);
    if (njumble > 1) {
      alloctree(&bestree2.nodep, nonodes2);
    }
  }
  allocrest();
}  /* doinit */


void getalleles()
{ /* set up number of alleles at loci */
  long i, j, m;

  if (!firstset)
    samenumsp(&loci, ith);
  if (contchars ) {
    totalleles = loci;
    for (i = 1; i <= loci; i++) {
      locus[i - 1] = i;
      alleles[i - 1] = 1;
    }
    df = loci;
  } else {
    totalleles = 0;
    scan_eoln(infile);
    if (printdata) {
      fprintf(outfile, "\nNumbers of alleles at the loci:\n");
      fprintf(outfile, "------- -- ------- -- --- -----\n\n");
    }
    for (i = 1; i <= loci; i++) {
      if (eoln(infile)) 
        scan_eoln(infile);
      if (fscanf(infile, "%ld", &alleles[i - 1]) != 1) {
        printf("ERROR: Unable to read number of alleles at locus %ld\n", i);
        exxit(-1);
      }
      if (alleles[i - 1] <= 0) {
    printf("ERROR: Bad number of alleles: %ld at locus %ld\n", alleles[i-1], i);
        exxit(-1);
      }
      totalleles += alleles[i - 1];
      if (printdata)
        fprintf(outfile, "%4ld", alleles[i - 1]);
    }
    locus = (long *)Malloc(totalleles*sizeof(long));
    m = 0;
    for (i = 1; i <= loci; i++) {
      for (j = 0; j < alleles[i - 1]; j++)
        locus[m+j] = i;
      m += alleles[i - 1];
    }
    df = totalleles - loci;
  }
  allocview(&curtree, nonodes2, totalleles);
  if (!usertree) {
    allocview(&bestree, nonodes2, totalleles);
    allocview(&priortree, nonodes2, totalleles);
    if (njumble > 1)
      allocview(&bestree2, nonodes2, totalleles);
  }
  for (i = 0; i < spp; i++)
    x[i] = (phenotype3)Malloc(totalleles*sizeof(double));
  pbar = (double *)Malloc(totalleles*sizeof(double));
  if (usertree)
    for (i = 0; i < MAXSHIMOTREES; i++)
      l0gf[i] = (double *)Malloc(totalleles*sizeof(double));
  if (printdata)
    putc('\n', outfile);
}  /* getalleles */


void inputdata()
{ /* read species data */
  long i, j, k, l, m, m0, n, p;
  double sum;

  if (printdata) {
    fprintf(outfile, "\nName");
    if (contchars)
      fprintf(outfile, "                       Phenotypes\n");
    else
      fprintf(outfile, "                 Gene Frequencies\n");
    fprintf(outfile, "----");
    if (contchars)
      fprintf(outfile, "                       ----------\n");
    else
      fprintf(outfile, "                 ---- -----------\n");
    putc('\n', outfile);
    if (!contchars) {
      for (j = 1; j <= nmlngth - 8; j++)
        putc(' ', outfile);
      fprintf(outfile, "locus:");
      p = 1;
      for (j = 1; j <= loci; j++) {
        if (all)
          n = alleles[j - 1];
        else
          n = alleles[j - 1] - 1;
        for (k = 1; k <= n; k++) {
          fprintf(outfile, "%10ld", j);
          if (p % 6 == 0 && (all || p < df)) {
            putc('\n', outfile);
            for (l = 1; l <= nmlngth - 2; l++)
              putc(' ', outfile);
          }
          p++;
        }
      }
      fprintf(outfile, "\n\n");
    }
  }
  for (i = 0; i < spp; i++) {
    scan_eoln(infile);
    initname(i);
    if (printdata)
      for (j = 0; j < nmlngth; j++)
        putc(nayme[i][j], outfile);
    m = 1;
    p = 1;
    for (j = 1; j <= loci; j++) {
      m0 = m;
      sum = 0.0;
      if (contchars)
        n = 1;
      else if (all)
        n = alleles[j - 1];
      else
        n = alleles[j - 1] - 1;
      for (k = 1; k <= n; k++) {
        if (eoln(infile)) 
          scan_eoln(infile);
        if (fscanf(infile, "%lf", &x[i][m - 1]) != 1) {
          printf("ERROR: unable to read allele frequency"
              "for species %ld, locus %ld\n", i+1, j);
          exxit(-1);
        }
        sum += x[i][m - 1];
        if (!contchars && x[i][m - 1] < 0.0) {
          printf("\n\nERROR: locus %ld in species %ld: an allele", j, i+1);
          printf(" frequency is negative\n");
          exxit(-1);
        }
        if (printdata) {
          fprintf(outfile, "%10.5f", x[i][m - 1]);
          if (p % 6 == 0 && (all || p < df)) {
            putc('\n', outfile);
            for (l = 1; l <= nmlngth; l++)
              putc(' ', outfile);
          }
        }
        p++;
        m++;
      }
      if (all && !contchars) {
        if (fabs(sum - 1.0) > epsilon2) {
          printf(
      "\n\nERROR: Locus %ld in species %ld: frequencies do not add up to 1\n",
                  j, i + 1);
          printf("\nFrequencies are:\n");
          for (l = m0; l <= m-3; l++)
            printf("%f+", x[i][l]);
          printf("%f = %f\n\n", x[i][m-2], sum);
          exxit(-1);
        } else {
          for (l = 0; l <= m-2; l++)
            x[i][l] /= sum;
        }
      }
      if (!all && !contchars) {
        x[i][m-1] = 1.0 - sum;
        if (x[i][m-1] < 0.0) {
          if (x[i][m-1] > -epsilon2) {
            for (l = 0; l <= m-2; l++)
              x[i][l] /= sum;
            x[i][m-1] = 0.0;
          } else {
            printf("\n\nERROR: Locus %ld in species %ld: ", j, i + 1);
            printf("frequencies add up to more than 1\n");
            printf("\nFrequencies are:\n");
            for (l = m0-1; l <= m-3; l++)
              printf("%f+", x[i][l]);
            printf("%f = %f\n\n", x[i][m-2], sum);
            exxit(-1);
          }
        }
        m++;
      }
    }
    if (printdata)
      putc('\n', outfile);
  }
  scan_eoln(infile);
  if (printdata)
    putc('\n', outfile);
}  /* inputdata */


void transformgfs()
{ /* do stereographic projection transformation on gene frequencies to
    get variables that come closer to independent Brownian motions */
  long i, j, k, l, m, n, maxalleles;
  double f, sum;
  double *sumprod, *sqrtp, *pbar;
  phenotype3 *c;

  sumprod = (double *)Malloc(loci*sizeof(double));
  sqrtp = (double *)Malloc(totalleles*sizeof(double));
  pbar = (double *)Malloc(totalleles*sizeof(double));
  for (i = 0; i < totalleles; i++) {  /* get mean gene frequencies */
    pbar[i] = 0.0;
    for (j = 0; j < spp; j++)
      pbar[i] += x[j][i];
    pbar[i] /= spp;
    if (pbar[i] == 0.0)
      sqrtp[i] = 0.0;
    else
      sqrtp[i] = sqrt(pbar[i]);
  }
  for (i = 0; i < spp; i++) {
    for (j = 0; j < loci; j++)    /* for each locus, sum of root(p*x) */
      sumprod[j] = 0.0;
    for (j = 0; j < totalleles; j++)
      if ((pbar[j]*x[i][j]) >= 0.0)
        sumprod[locus[j]-1] += sqrtp[j]*sqrt(x[i][j]);
    for (j = 0; j < totalleles; j++) { /* the projection to tangent plane */
      f = (1.0 + sumprod[locus[j]-1])/2.0;
      if (x[i][j] == 0.0)
        x[i][j] = (2.0/f - 1.0)*sqrtp[j];
      else
        x[i][j] = (1.0/f)*sqrt(x[i][j]) + (1.0/f - 1.0)*sqrtp[j];
    }
  }
  maxalleles = 0;
  for (i = 0; i < loci; i++)
    if (alleles[i] > maxalleles)
      maxalleles = alleles[i];
  c = (phenotype3 *)Malloc(maxalleles*sizeof(phenotype3));
  for (i = 0; i < maxalleles; i++)  /* enough room for any locus's contrasts */
    c[i] = (double *)Malloc(maxalleles*sizeof(double));
  m = 0;        
  for (j = 0; j < loci; j++) {            /* do this for each locus */
    for (k = 0; k < alleles[j]-1; k++) {  /* one fewer than # of alleles */
      c[k][0] = 1.0;
      for (l = 0; l < k; l++) {          /* for contrasts 1 to k make it ... */
        sum = 0.0;
        for (n = 0; n <= l; n++)
          sum += c[k][n]*c[l][n];
        if (fabs(c[l][l+1]) > 0.000000001) /* ... orthogonal to those ones */
          c[k][l+1] = -sum / c[l][l+1];    /* set coeff to make orthogonal */
        else
          c[k][l+1] = 1.0;
      }
      sum = 0.0;
      for (l = 0; l <= k; l++) /* make it orthogonal to vector of sqrtp's */
        sum += c[k][l]*sqrtp[m+l];
      if (sqrtp[m+k+1] > 0.0000000001)
        c[k][k+1] = - sum / sqrtp[m+k+1];  /* ... setting last coeff */
      else {
        for (l = 0; l <= k; l++)
          c[k][l] = 0.0;
        c[k][k+1] = 1.0;
      }
      sum = 0.0;
      for (l = 0; l <= k+1; l++)
        sum += c[k][l]*c[k][l];
      sum = sqrt(sum);
      for (l = 0; l <= k+1; l++)
        if (sum > 0.0000000001)
          c[k][l] /= sum;
    }
    for (i = 0; i < spp; i++) {     /* the orthonormal axes in the plane */
      for (l = 0; l < alleles[j]-1; l++) { /* compute the l-th one */
        c[maxalleles-1][l] = 0.0;  /* temporarily store it ... */
        for (n = 0; n <= l+1; n++) 
          c[maxalleles-1][l] += c[l][n]*x[i][m+n];
      }
      for (l = 0; l < alleles[j]-1; l++)
        x[i][m+l] = c[maxalleles-1][l];   /* replace the gene freqs by it */
    }
    m += alleles[j];
  }
  for (i = 0; i < maxalleles; i++)
    free(c[i]);
  free(c);
  free(sumprod);
  free(sqrtp);
  free(pbar);
} /* transformgfs */


void getinput()
{ /* reads the input data */
  getalleles();
  inputdata();
  if (!contchars) {
    transformgfs();
  }
}  /* getinput */


void sumlikely(node *p, node *q, double *sum)
{ /* sum contribution to likelihood over forks in tree */
  long i, j, m;
  double term, sumsq, vee;
  double temp;

  if (!p->tip)
    sumlikely(p->next->back, p->next->next->back, sum);
  if (!q->tip)
    sumlikely(q->next->back, q->next->next->back, sum);
  if (p->back == q)
    vee = p->v;
  else
    vee = p->v + q->v;
  vee += p->deltav + q->deltav;
  if (vee <= 1.0e-10) {
    printf("ERROR: check for two identical species ");
    printf("and eliminate one from the data\n");
    exxit(-1);
  }
  sumsq = 0.0;
  if (usertree && which <= MAXSHIMOTREES) {
    for (i = 0; i < loci; i++)
      l0gf[which - 1][i] += (1 - alleles[i]) * log(vee) / 2.0;
  }
  if (contchars) {
    m = 0;
    for (i = 0; i < loci; i++) {
      temp = p->view[i] - q->view[i];
      term = temp * temp;
      if (usertree && which <= MAXSHIMOTREES)
        l0gf[which - 1][i] -= term / (2.0 * vee);
      sumsq += term;
    }
  } else {
    m = 0;
    for (i = 0; i < loci; i++) {
      for (j = 1; j < alleles[i]; j++) {
        temp = p->view[m+j-1] - q->view[m+j-1];
        term = temp * temp;
        if (usertree && which <= MAXSHIMOTREES)
          l0gf[which - 1][i] -= term / (2.0 * vee);
        sumsq += term;
      }
      m += alleles[i];
    }
  }
  (*sum) += df * log(vee) / -2.0 - sumsq / (2.0 * vee);
}  /* sumlikely */


double evaluate(tree *t)
{ /* evaluate likelihood of a tree */
  long i;
  double sum;

  sum = 0.0;
  if (usertree && which <= MAXSHIMOTREES) {
    for (i = 0; i < loci; i++)
      l0gf[which - 1][i] = 0.0;
  }
  sumlikely(t->start->back, t->start, &sum);
  if (usertree && which <= MAXSHIMOTREES) {
    l0gl[which - 1] = sum;
    if (which == 1) {
      maxwhich = 1;
      maxlogl = sum;
    } else if (sum > maxlogl) {
      maxwhich = which;
      maxlogl = sum;
    }
  }
  t->likelihood = sum;
  return sum;
}  /* evaluate */


double distance(node *p, node *q)
{ /* distance between two nodes */
  long i, j, m;
  double sum, temp;

  sum = 0.0;
  if (!contchars) {
    m = 0;
    for (i = 0; i < loci; i++) {
      for (j = 0; j < alleles[i]-1; j++) {
        temp = p->view[m+j] - q->view[m+j];
        sum += temp * temp;
      }
      m += alleles[i];
    }
  }
  else {
    for (i = 0; i < totalleles; i++) {
      temp = p->view[i] - q->view[i];
      sum += temp * temp;
    }
  }
  return sum;
}  /* distance */


void makedists(node *p)
{ /* compute distances among three neighbors of a node */
  long i;
  node *q;

  for (i = 1; i <= 3; i++) {
    q = p->next;
    p->dist = distance(p->back, q->back);
    p = q;
  }
}  /* makedists */


void makebigv(node *p, boolean *negatives)
{ /* make new branch length */
  long i;
  node *temp, *q, *r;

  q = p->next;
  r = q->next;
  *negatives = false;
  for (i = 1; i <= 3; i++) {
    p->bigv = p->v + p->back->deltav;
    if (p->iter) {
      p->bigv = (p->dist + r->dist - q->dist) / (df * 2);
      p->back->bigv = p->bigv;
      if (p->bigv < p->back->deltav)
        *negatives = true;
    }
    temp = p;
    p = q;
    q = r;
    r = temp;
  }
}  /* makebigv */


void correctv(node *p)
{ /* iterate branch lengths if some are to be zero */
  node *q, *r, *temp;
  long i, j;
  double f1, f2, vtot;

  q = p->next;
  r = q->next;
  for (i = 1; i <= smoothings; i++) {
    for (j = 1; j <= 3; j++) {
      vtot = q->bigv + r->bigv;
      if (vtot > 0.0)
        f1 = q->bigv / vtot;
      else
        f1 = 0.5;
      f2 = 1.0 - f1;
      p->bigv = (f1 * r->dist + f2 * p->dist - f1 * f2 * q->dist) / df;
      p->bigv -= vtot * f1 * f2;
      if (p->bigv < p->back->deltav)
        p->bigv = p->back->deltav;
      p->back->bigv = p->bigv;
      temp = p;
      p = q;
      q = r;
      r = temp;
    }
  }
}  /* correctv */


void littlev(node *p)
{ /* remove part of it that belongs to other barnches */
  long i;

  for (i = 1; i <= 3; i++) {
    if (p->iter)
      p->v = p->bigv - p->back->deltav;
    if (p->back->iter)
      p->back->v = p->v;
    p = p->next;
  }
}  /* littlev */


void nuview(node *p)
{ /* renew information about subtrees */
  long i, j, k, m;
  node *q, *r, *a, *b, *temp;
  double v1, v2, vtot, f1, f2;

  q = p->next;
  r = q->next;
  for (i = 1; i <= 3; i++) {
    a = q->back;
    b = r->back;
    v1 = q->bigv;
    v2 = r->bigv;
    vtot = v1 + v2;
    if (vtot > 0.0)
      f1 = v2 / vtot;
    else
      f1 = 0.5;
    f2 = 1.0 - f1;
    m = 0;
    for (j = 0; j < loci; j++) {
      for (k = 1; k <= alleles[j]; k++)
        p->view[m+k-1] = f1 * a->view[m+k-1] + f2 * b->view[m+k-1];
      m += alleles[j];
    }
    p->deltav = v1 * f1;
    temp = p;
    p = q;
    q = r;
    r = temp;
  }
}  /* nuview */


void update(node *p)
{ /* update branch lengths around a node */
  boolean negatives;

  if (p->tip)
    return;
  makedists(p);
  makebigv(p,&negatives);
  if (negatives)
    correctv(p);
  littlev(p);
  nuview(p);
}  /* update */


void smooth(node *p)
{ /* go through tree getting new branch lengths and views */
  if (p->tip)
    return;
  update(p);
  smooth(p->next->back);
  smooth(p->next->next->back);
}  /* smooth */


void insert_(node *p, node *q)
{ /* put p and q together and iterate info. on resulting tree */
  long i;

  hookup(p->next->next, q->back);
  hookup(p->next, q);
  for (i = 1; i <= smoothings; i++) {
    smooth(p);
    smooth(p->back);
  }
}  /* insert_ */


void copynode(node *c, node *d)
{ /* make a copy of a node */
  memcpy(d->view, c->view, totalleles*sizeof(double));
  d->v = c->v;
  d->iter = c->iter;
  d->deltav = c->deltav;
  d->bigv = c->bigv;
  d->dist = c->dist;
  d->xcoord = c->xcoord;
  d->ycoord = c->ycoord;
  d->ymin = c->ymin;
  d->ymax = c->ymax;
}  /* copynode */


void copy_(tree *a, tree *b)
{ /* make a copy of tree a to tree b */
  long i, j;
  node *p, *q;

  for (i = 0; i < spp; i++) {
    copynode(a->nodep[i], b->nodep[i]);
    if (a->nodep[i]->back) {
      if (a->nodep[i]->back == a->nodep[a->nodep[i]->back->index - 1])
        b->nodep[i]->back = b->nodep[a->nodep[i]->back->index - 1];
      else if (a->nodep[i]->back == a->nodep[a->nodep[i]->back->index - 1]->next)
        b->nodep[i]->back = b->nodep[a->nodep[i]->back->index - 1]->next;
      else
        b->nodep[i]->back = b->nodep[a->nodep[i]->back->index - 1]->next->next;
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


void inittip(long m, tree *t)
{ /* initialize branch lengths and views in a tip */
  node *tmp;

  tmp = t->nodep[m - 1];
  memcpy(tmp->view, x[m - 1], totalleles*sizeof(double));
  tmp->deltav = 0.0;
  tmp->v = 0.0;
}  /* inittip */


void buildnewtip(long m, tree *t, long nextsp)
{ /* initialize and hook up a new tip */
  node *p;

  inittip(m, t);
  p = t->nodep[nextsp + spp - 3];
  hookup(t->nodep[m - 1], p);
}  /* buildnewtip */


void buildsimpletree(tree *t)
{ /* make and initialize a three-species tree */
  inittip(enterorder[0], t);
  inittip(enterorder[1], t);
  hookup(t->nodep[enterorder[0] - 1], t->nodep[enterorder[1] - 1]);
  buildnewtip(enterorder[2], t, nextsp);
  insert_(t->nodep[enterorder[2] - 1]->back, t->nodep[enterorder[0] - 1]);
}  /* buildsimpletree */


void addtraverse(node *p, node *q, boolean contin)
{ /* traverse through a tree, finding best place to add p */
  insert_(p, q);
  numtrees++;
  if (evaluate(&curtree) > bestree.likelihood) {
    copy_(&curtree, &bestree);
    addwhere = q;
  }
  copy_(&priortree, &curtree);
  if (!q->tip && contin) {
    addtraverse(p, q->next->back, contin);
    addtraverse(p, q->next->next->back, contin);
  }
}  /* addtraverse */


void re_move(node **p, node **q)
{ /* remove p and record in q where it was */
  *q = (*p)->next->back;
  hookup(*q, (*p)->next->next->back);
  (*p)->next->back = NULL;
  (*p)->next->next->back = NULL;
  update(*q);
  update((*q)->back);
}  /* re_move */


void globrearrange() 
{ /* does global rearrangements */
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
  allocview(&oldtree, nonodes2, totalleles);
  allocview(&globtree, nonodes2, totalleles);
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
        sib_ptr = sib_ptr->next;
        continue;
      }
      else num_sibs2 = count_sibs(where);
      sib_ptr2 = where;
      for ( k = 0 ; k < num_sibs2 ; k++ ) {
        addwhere = NULL;
        addtraverse(sib_ptr,sib_ptr2->back,true);
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
    succeeded = true;
  }
  else  {
    succeeded = false;
  }
  freeview(&oldtree, nonodes2);
  freeview(&globtree, nonodes2);
  freetree(&globtree.nodep,nonodes2);
  freetree(&oldtree.nodep,nonodes2);
}

void rearrange(node *p)
{ /* rearranges the tree locally */
  node *q, *r;

  if (!p->tip && !p->back->tip) {
    r = p->next->next;
    re_move(&r, &q );
    copy_(&curtree, &priortree);
    addtraverse(r, q->next->back, false);
    addtraverse(r, q->next->next->back, false);
    copy_(&bestree, &curtree);
  }
  if (!p->tip) {
    rearrange(p->next->back);
    rearrange(p->next->next->back);
  }
}  /* rearrange */


void coordinates(node *p, double lengthsum, long *tipy, double *tipmax)
{ /* establishes coordinates of nodes */
  node *q, *first, *last;

  if (p->tip) {
    p->xcoord = lengthsum;
    p->ycoord = *tipy;
    p->ymin = *tipy;
    p->ymax = *tipy;
    (*tipy) += down;
    if (lengthsum > (*tipmax))
      (*tipmax) = lengthsum;
    return;
  }
  q = p->next;
  do {
    coordinates(q->back, lengthsum + q->v, tipy,tipmax);
    q = q->next;
  } while ((p == curtree.start || p != q) &&
           (p != curtree.start || p->next != q));
  first = p->next->back;
  q = p;
  while (q->next != p)
    q = q->next;
  last = q->back;
  p->xcoord = lengthsum;
  if (p == curtree.start)
    p->ycoord = p->next->next->back->ycoord;
  else
    p->ycoord = (first->ycoord + last->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* coordinates */


void drawline(long i, double scale)
{ /* draws one row of the tree diagram by moving up tree */
  node *p, *q;
  long n, j;
  boolean extra;
  node *r, *first = NULL, *last = NULL;
  boolean done;

  p = curtree.start;
  q = curtree.start;
  extra = false;
  if (i == (long)p->ycoord && p == curtree.start) {
    if (p->index - spp >= 10)
      fprintf(outfile, " %2ld", p->index - spp);
    else
      fprintf(outfile, "  %ld", p->index - spp);
    extra = true;
  } else
    fprintf(outfile, "  ");
  do {
    if (!p->tip) {
      r = p->next;
      done = false;
      do {
        if (i >= (long)r->back->ymin && i <= (long)r->back->ymax) {
          q = r->back;
          done = true;
        }
        r = r->next;
      } while (!(done || (p != curtree.start && r == p) ||
                 (p == curtree.start && r == p->next)));
      first = p->next->back;
      r = p;
      while (r->next != p)
        r = r->next;
      last = r->back;
      if (p == curtree.start)
        last = p->back;
    }
    done = (p->tip || p == q);
    n = (long)(scale * (q->xcoord - p->xcoord) + 0.5);
    if (n < 3 && !q->tip)
      n = 3;
    if (extra) {
      n--;
      extra = false;
    }
    if ((long)q->ycoord == i && !done) {
      if ((long)p->ycoord != (long)q->ycoord)
        putc('+', outfile);
      else
        putc('-', outfile);
      if (!q->tip) {
        for (j = 1; j <= n - 2; j++)
          putc('-', outfile);
        if (q->index - spp >= 10)
          fprintf(outfile, "%2ld", q->index - spp);
        else
          fprintf(outfile, "-%ld", q->index - spp);
        extra = true;
      } else {
        for (j = 1; j < n; j++)
          putc('-', outfile);
      }
    } else if (!p->tip) {
      if ((long)last->ycoord > i && (long)first->ycoord < i
           && i != (long)p->ycoord) {
        putc('!', outfile);
        for (j = 1; j < n; j++)
          putc(' ', outfile);
      } else {
        for (j = 1; j <= n; j++)
          putc(' ', outfile);
      }
    } else {
      for (j = 1; j <= n; j++)
        putc(' ', outfile);
    }
    if (q != p)
      p = q;
  } while (!done);
  if ((long)p->ycoord == i && p->tip) {
    for (j = 0; j < nmlngth; j++)
      putc(nayme[p->index - 1][j], outfile);
  }
  putc('\n', outfile);
}  /* drawline */


void printree()
{ /* prints out diagram of the tree */
  long i;
  long tipy;
  double tipmax,scale;

  if (!treeprint)
    return;
  putc('\n', outfile);
  tipy = 1;
  tipmax = 0.0;
  coordinates(curtree.start, 0.0, &tipy,&tipmax);
  scale = over / (tipmax + 0.0001);
  for (i = 1; i <= (tipy - down); i++)
    drawline(i,scale);
  putc('\n', outfile);
}  /* printree */


void treeout(node *p)
{ /* write out file with representation of final tree */
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
    treeout(p->next->back);
    putc(',', outtree);
    col++;
    if (col > 55) {
      putc('\n', outtree);
      col = 0;
    }
    treeout(p->next->next->back);
    if (p == curtree.start) {
      putc(',', outtree);
      col++;
      if (col > 45) {
        putc('\n', outtree);
        col = 0;
      }
      treeout(p->back);
    }
    putc(')', outtree);
    col++;
  }
  x = p->v;
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
    fprintf(outtree, ":%*.8f", (int)w + 7, x);
    col += w + 8;
  }
}  /* treeout */


void describe(node *p, double chilow, double chihigh)
{ /* print out information for one branch */
  long i;
  node *q;
  double bigv, delta;

  q = p->back;
  fprintf(outfile, "%3ld       ", q->index - spp);
  if (p->tip) {
    for (i = 0; i < nmlngth; i++)
      putc(nayme[p->index - 1][i], outfile);
  } else
    fprintf(outfile, "%4ld      ", p->index - spp);
  fprintf(outfile, "%15.8f", q->v);
  delta = p->deltav + p->back->deltav;
  bigv = p->v + delta;
  if (p->iter)
     fprintf(outfile, "   (%12.8f,%12.8f)",
             chilow * bigv - delta, chihigh * bigv - delta);
  fprintf(outfile, "\n");
  if (!p->tip) {
    describe(p->next->back, chilow,chihigh);
    describe(p->next->next->back, chilow,chihigh);
  }
}  /* describe */


void summarize(void)
{ /* print out branch lengths etc. */
  double chilow,chihigh;

  fprintf(outfile, "\nremember: ");
  if (outgropt)
    fprintf(outfile, "(although rooted by outgroup) ");
  fprintf(outfile, "this is an unrooted tree!\n\n");
  fprintf(outfile, "Ln Likelihood = %11.5f\n", curtree.likelihood);
  if (df == 1) {
    chilow = 0.000982;
    chihigh = 5.02389;
  } else if (df == 2) {
    chilow = 0.05064;
    chihigh = 7.3777;
  } else {
    chilow = 1.0 - 2.0 / (df * 9);
    chihigh = chilow;
    chilow -= 1.95996 * sqrt(2.0 / (df * 9));
    chihigh += 1.95996 * sqrt(2.0 / (df * 9));
    chilow *= chilow * chilow;
    chihigh *= chihigh * chihigh;
  }
  fprintf(outfile, "\nBetween     And             Length");
  fprintf(outfile, "      Approx. Confidence Limits\n");
  fprintf(outfile, "-------     ---             ------");
  fprintf(outfile, "      ------- ---------- ------\n");
  describe(curtree.start->next->back, chilow,chihigh);
  describe(curtree.start->next->next->back, chilow,chihigh);
  describe(curtree.start->back, chilow, chihigh);
  fprintf(outfile, "\n\n");
  if (trout) {
    col = 0;
    treeout(curtree.start);
  }
}  /* summarize */


void nodeinit(node *p)
{ /* initialize a node */
  node *q, *r;
  long i, j, m;

  if (p->tip)
    return;
  q = p->next->back;
  r = p->next->next->back;
  nodeinit(q);
  nodeinit(r);
  m = 0;
  for (i = 0; i < loci; i++) {
    for (j = 1; j < alleles[i]; j++)
      p->view[m+j-1] = 0.5 * q->view[m+j-1] + 0.5 * r->view[m+j-1];
    m += alleles[i];
  }
  if ((!lengths) || p->iter)
    p->v = 0.1;
  if ((!lengths) || p->back->iter)
    p->back->v = 0.1;
}  /* nodeinit */


void initrav(node *p)
{ /* traverse to initialize */
  node* q; 

  if (p->tip)
    nodeinit(p->back);
  else {
    q = p->next;
    while ( q != p)  {
      initrav(q->back);
      q = q->next;
    }
  }
}  /* initrav */


void treevaluate()
{ /* evaluate user-defined tree, iterating branch lengths */
  long i;
  
  unroot(&curtree,nonodes2);
  initrav(curtree.start);
  initrav(curtree.start->back);
  for (i = 1; i <= smoothings * 4; i++)
    smooth(curtree.start);
  evaluate(&curtree);
}  /* treevaluate */


void maketree()
{ /* construct the tree */
  long i;

  if (usertree) {
    /* Open in binary: ftell() is broken for UNIX line-endings under WIN32 */
    openfile(&intree,INTREE,"input tree file", "rb",progname,intreename);
    numtrees = countsemic(&intree);
    if(numtrees > MAXSHIMOTREES)
      shimotrees = MAXSHIMOTREES;
    else
      shimotrees = numtrees;
    if (numtrees > 2)
      initseed(&inseed, &inseed0, seed);
    if (treeprint) {
      fprintf(outfile, "User-defined tree");
      if (numtrees > 1)
        putc('s', outfile);
      putc('\n', outfile);
    }
    setuptree(&curtree, nonodes2);
    for (which = 1; which <= spp; which++)
      inittip(which, &curtree);
    which = 1;
    while (which <= numtrees) {
      for (i = 0 ; i < nonodes2 ; i++) {
        if ( i > spp) {
          /* must do this since not all nodes may be used if an 
             unrooted tree is read in after a rooted one */
          curtree.nodep[i]->back = NULL;
          curtree.nodep[i]->next->back = NULL;
          curtree.nodep[i]->next->next->back = NULL;
        }
        else curtree.nodep[i]->back = NULL;
      }
      treeread2 (intree, &curtree.start, curtree.nodep,
        lengths, &trweight, &goteof, &haslengths, &spp,false,nonodes2);
      curtree.start = curtree.nodep[outgrno - 1]->back;
      treevaluate();
      printree();
      summarize();
      which++;
    }
    FClose(intree);
    if (numtrees > 1 && loci > 1 ) {
      weight = (long *)Malloc(loci*sizeof(long));
      for (i = 0; i < loci; i++)
        weight[i] = 1;
      standev2(numtrees, maxwhich, 0, loci-1, maxlogl, l0gl, l0gf, seed);
      free(weight);
      fprintf(outfile, "\n\n");
    }
  } else { /* if ( !usertree ) */
    if (jumb == 1) {
      setuptree(&curtree, nonodes2);
      setuptree(&priortree, nonodes2);
      setuptree(&bestree, nonodes2);
      if (njumble > 1) 
        setuptree(&bestree2, nonodes2);
    }
    for (i = 1; i <= spp; i++)
      enterorder[i - 1] = i;
    if (jumble)
      randumize(seed, enterorder);
    nextsp = 3;
    buildsimpletree(&curtree);
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
      buildnewtip(enterorder[nextsp - 1], &curtree, nextsp);
      copy_(&curtree, &priortree);
      bestree.likelihood = -DBL_MAX;
      addtraverse(curtree.nodep[enterorder[nextsp - 1] - 1]->back,
                  curtree.start, true );
      copy_(&bestree, &curtree);
      if (progress) {
        writename(nextsp - 1, 1, enterorder);
#ifdef WIN32
        phyFillScreenColor();
#endif
      }
      if (global && nextsp == spp) {
        if (progress) {
          printf("\nDoing global rearrangements\n");
          printf("  !");
          for (i = 1; i <= spp - 2; i++)
            if ( (i - spp) % (( nonodes2 / 72 ) + 1 ) == 0 )
              putchar('-');
          printf("!\n");
          printf("   ");
        }
      }
      succeeded = true;
      while (succeeded) {
        succeeded = false;
        if ( global && nextsp == spp )
          globrearrange();
        else 
          rearrange(curtree.start);
        if (global && nextsp == spp)
          putc('\n', outfile);
      }
      if (global && nextsp == spp && progress)
        putchar('\n');
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
        printree();
        summarize();
      }
      nextsp++;
    }
  }
  if ( jumb < njumble)
    return;
  if (progress) {
    printf("\n\nOutput written to file \"%s\"\n\n", outfilename);
    if (trout)
      printf("Tree also written onto file \"%s\"\n\n", outtreename);
  }
  freeview(&curtree, nonodes2);
  if (!usertree) {
    freeview(&bestree, nonodes2);
    freeview(&priortree, nonodes2);
  }
  for (i = 0; i < spp; i++)
    free(x[i]);
  if (!contchars) {
    free(locus);
    free(pbar);
  }
}  /* maketree */


int main(int argc, Char *argv[])
{  /* main program */
  long i;

#ifdef MAC
  argc = 1;                /* macsetup("Contml","");                */
  argv[0] = "Contml";
#endif
  init(argc, argv);
  progname = argv[0];
  openfile(&infile,INFILE,"input file", "r",argv[0],infilename);
  openfile(&outfile,OUTFILE,"output file", "w",argv[0],outfilename);
  ibmpc = IBMCRT;
  ansi = ANSICRT;
  mulsets = false;
  firstset = true;
  datasets = 1;
  doinit();
  if (trout)
    openfile(&outtree,OUTTREE,"output tree file", "w",argv[0],outtreename);
  for (ith = 1; ith <= datasets; ith++) {
    getinput();
    if (ith == 1)
      firstset = false;
    if (datasets > 1) {
      fprintf(outfile, "Data set # %ld:\n\n", ith);
      if (progress)
        printf("\nData set # %ld:\n", ith);
    }
    for (jumb = 1; jumb <= njumble; jumb++)
      maketree();
    if (usertree)
      for (i = 0; i < MAXSHIMOTREES; i++)
        free(l0gf[i]);
  }
  FClose(outfile);
  FClose(outtree);
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

