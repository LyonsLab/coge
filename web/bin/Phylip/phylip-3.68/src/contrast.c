
/* version 3.6. (c) Copyright 1993-2004 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#include "phylip.h"
#include "cont.h"



#ifndef OLDC
/* function prototypes */
void   getoptions(void);
void   getdata(void);
void   allocrest(void);
void   doinit(void);
void   contwithin(void);
void   contbetween(node *, node *);
void   nuview(node *);
void   makecontrasts(node *);
void   writecontrasts(void);

void   regressions(void);
double logdet(double **);
void   invert(double **);
void   initcovars(boolean);
double normdiff(boolean);
void   matcopy(double **, double **);
void   newcovars(boolean);
void   printcovariances(boolean);
void   emiterate(boolean);
void   initcontrastnode(node **, node **, node *, long, long, long *,
             long *, initops, pointarray, pointarray, Char *, Char *, FILE *);

void   maketree(void);
/* function prototypes */
#endif


Char infilename[FNMLNGTH], outfilename[FNMLNGTH], intreename[FNMLNGTH];
long nonodes, chars, numtrees;
long *sample, contnum;
phenotype3 **x, **cntrast, *ssqcont;
double **vara, **vare, **oldvara, **oldvare,
       **Bax, **Bex, **temp1, **temp2, **temp3;
double logL, logLvara, logLnovara;
boolean nophylo, printdata, progress, reg, mulsets,
  varywithin, writecont, bifurcating;
Char ch;
long contno;
node *grbg;

/* Local variables for maketree, propagated globally for c version: */
tree curtree;

/* Variables declared just to make treeread happy */
boolean haslengths, goteof, first;
double trweight;


void getoptions()
{
  /* interactively set options */
  long loopcount, loopcount2;
  Char ch;
  boolean done, done1;

  mulsets = false;
  nophylo = false;
  printdata = false;
  progress = true;
  varywithin = false;
  writecont = false;
  loopcount = 0;
  do {
    cleerhome();
    printf("\nContinuous character comparative analysis, version %s\n\n",
             VERSION);
    printf("Settings for this run:\n");
    printf("  W        Within-population variation in data?");
    if (varywithin)
      printf("  Yes, multiple individuals\n");
    else {
      printf("  No, species values are means\n");
      printf("  R     Print out correlations and regressions?  %s\n",
             (reg ? "Yes" : "No"));
    }
    if (varywithin) {
      printf("  A      LRT test of no phylogenetic component?");
      if (nophylo)
        printf("  Yes, with and without VarA\n");
      else
        printf("  No, just assume it is there\n");
    }
    if (!varywithin)
      printf("  C                        Print out contrasts?  %s\n",
               (writecont? "Yes" : "No"));
    printf("  M                     Analyze multiple trees?");
    if (mulsets)
      printf("  Yes, %2ld trees\n", numtrees);
    else
      printf("  No\n");
    printf("  0         Terminal type (IBM PC, ANSI, none)?  %s\n",
           ibmpc ? "IBM PC"  :
           ansi  ? "ANSI"    : "(none)");
    printf("  1          Print out the data at start of run  %s\n",
           (printdata ? "Yes" : "No"));
    printf("  2        Print indications of progress of run  %s\n",
           (progress ? "Yes" : "No"));
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
    done = (ch == 'Y');
    if (!done) {
      if (strchr("RAMWC120", ch) != NULL) {
        switch (ch) {

        case 'R':
          reg = !reg;
          break;

        case 'A':
          nophylo = !nophylo;
          break;

        case 'M':
          mulsets = !mulsets;
          if (mulsets) {
            loopcount2 = 0;
            do {
              printf("How many trees?\n");
#ifdef WIN32
              phyFillScreenColor();
#endif
              fflush(stdout);
              scanf("%ld%*[^\n]", &numtrees);
              getchar();
              done1 = (numtrees >= 1);
              if (!done1)
                printf("BAD TREES NUMBER:  it must be greater than 1\n");
              countup(&loopcount2, 10);
            } while (done1 != true);
          }
          break;

        case 'C':
          writecont = !writecont;
          break;

        case 'W':
          varywithin = !varywithin;
          if (!nophylo)
            nophylo = true;
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
        }
      } else
        printf("Not a possible option!\n");
    }
    countup(&loopcount, 100);
  } while (!done);
}  /* getoptions */


void getdata()
{
  /* read species data */
  long i, j, k, l;

  if (printdata) {
    fprintf(outfile,
       "\nContinuous character contrasts analysis, version %s\n\n",VERSION);
    fprintf(outfile, "%4ld Populations, %4ld Characters\n\n", spp, chars);
    fprintf(outfile, "Name");
    fprintf(outfile, "                       Phenotypes\n");
    fprintf(outfile, "----");
    fprintf(outfile, "                       ----------\n\n");
  }
  x = (phenotype3 **)Malloc((long)spp*sizeof(phenotype3 *));
  cntrast = (phenotype3 **)Malloc((long)spp*sizeof(phenotype3 *));
  ssqcont = (phenotype3 *)Malloc((long)spp*sizeof(phenotype3 *));
  contnum = spp-1;
  for (i = 0; i < spp; i++) {
    scan_eoln(infile);
    initname(i);
    if (varywithin) {
      if (fscanf(infile, "%ld", &sample[i]) == 1) {
        contnum += sample[i]-1;
        scan_eoln(infile);
      }
      else {
        printf("Error reading number of individuals for species %ld\n", i+1);
        exxit(-1);
      }
    }
    else sample[i] = 1;
    if (printdata)
      for(j = 0; j < nmlngth; j++)
        putc(nayme[i][j], outfile);
    x[i] = (phenotype3 *)Malloc((long)sample[i]*sizeof(phenotype3));
    cntrast[i] = (phenotype3 *)Malloc((long)(sample[i]*sizeof(phenotype3)));
    ssqcont[i] = (double *)Malloc((long)(sample[i]*sizeof(double)));
    for (k = 0; k <= sample[i]-1; k++) {
      x[i][k] = (phenotype3)Malloc((long)chars*sizeof(double));
      cntrast[i][k] = (phenotype3)Malloc((long)chars*sizeof(double));
      for (j = 1; j <= chars; j++) {
        if (eoln(infile)) 
          scan_eoln(infile);
        if (fscanf(infile, "%lf", &x[i][k][j - 1]) != 1) {
          printf("Error in input file at species %ld\n", i+1);
          exxit(-1);
        }
        if (printdata) {
          fprintf(outfile, " %9.5f", x[i][k][j - 1]);
          if (j % 6 == 0) {
            putc('\n', outfile);
            for (l = 1; l <= nmlngth; l++)
              putc(' ', outfile);
          }
        }
      }
    }
    if (printdata)
      putc('\n', outfile);
  }
  scan_eoln(infile);
  if (printdata)
    putc('\n', outfile);
}  /* getdata */


void allocrest()
{
  long i;

  /* otherwise if individual variation, these are allocated in getdata  */
  sample = (long *)Malloc((long)spp*sizeof(long));
  nayme = (naym *)Malloc((long)spp*sizeof(naym));
  vara = (double **)Malloc((long)chars*sizeof(double *));
  oldvara = (double **)Malloc((long)chars*sizeof(double *));
  vare = (double **)Malloc((long)chars*sizeof(double *));
  oldvare = (double **)Malloc((long)chars*sizeof(double *));
  Bax = (double **)Malloc((long)chars*sizeof(double *));
  Bex = (double **)Malloc((long)chars*sizeof(double *));
  temp1 = (double **)Malloc((long)chars*sizeof(double *));
  temp2 = (double **)Malloc((long)chars*sizeof(double *));
  temp3 = (double **)Malloc((long)chars*sizeof(double *));
  for (i = 0; i < chars; i++) {
     vara[i] = (double *)Malloc((long)chars*sizeof(double));
     oldvara[i] = (double *)Malloc((long)chars*sizeof(double));
     vare[i] = (double *)Malloc((long)chars*sizeof(double));
     oldvare[i] = (double *)Malloc((long)chars*sizeof(double));
     Bax[i] = (double *)Malloc((long)chars*sizeof(double));
     Bex[i] = (double *)Malloc((long)chars*sizeof(double));
     temp1[i] = (double *)Malloc((long)chars*sizeof(double));
     temp2[i] = (double *)Malloc((long)chars*sizeof(double));
     temp3[i] = (double *)Malloc((long)chars*sizeof(double));
  }
}  /* allocrest */


void doinit()
{
  /* initializes variables */

  inputnumbers(&spp, &chars, &nonodes, 1);
  getoptions();
  allocrest();
}  /* doinit */


void contwithin()
{
  /* compute the within-species contrasts, if any */
  long i, j, k;
  double *sumphen;

  sumphen = (double *)Malloc((long)chars*sizeof(double));
  for (i = 0; i <= spp-1 ; i++) {
    for (j = 0; j < chars; j++)
      sumphen[j] = 0.0;
    for (k = 0; k <= (sample[i]-1); k++) {
      for (j = 0; j < chars; j++) {
        if (k > 0)
          cntrast[i][k][j]
              = (sumphen[j] - k*x[i][k][j])/sqrt((double)(k*(k+1)));
        sumphen[j] += x[i][k][j];
        if (k == (sample[i]-1))
          curtree.nodep[i]->view[j] = sumphen[j]/sample[i];
        x[i][0][j] = sumphen[j]/sample[i];
      }
      if (k == 0)
        curtree.nodep[i]->ssq = 1.0/sample[i]; /* sum of squares for sp. i */
      else
        ssqcont[i][k] = 1.0;   /* if a within contrast */
    }
  }
  contno = 1;
}  /* contwithin */


void contbetween(node *p, node *q)
{
  /* compute one contrast */
  long i;
  double v1, v2;

  if ((p->v < 0.0) || (q->v < 0.0)) {
    printf("\nERROR: input tree has a negative branch length,");
    printf(" which is not allowed.\n\n");
    exxit(-1);
  }
  for (i = 0; i < chars; i++)
    cntrast[contno - 1][0][i] = (p->view[i] - q->view[i])/sqrt(p->ssq+q->ssq);
  v1 = q->v + q->deltav;
  if (p->back != q)
    v2 = p->v + p->deltav;
  else v2 = p->deltav;
  ssqcont[contno - 1][0] = (v1 + v2)/(p->ssq + q->ssq);
     /* this is really the variance of the contrast */
  contno++;
}  /* contbetween */


void nuview(node *p)
{
  /* renew information about subtrees */
  long j;
  node *q, *r;
  double v1, v2, vtot, f1, f2;

  q = p->next->back;
  r = p->next->next->back;
  v1 = q->v + q->deltav;
  v2 = r->v + r->deltav;
  vtot = v1 + v2;
  if (vtot > 0.0)
    f1 = v2 / vtot;
  else
    f1 = 0.5;
  f2 = 1.0 - f1;
  for (j = 0; j < chars; j++)
    p->view[j] = f1 * q->view[j] + f2 * r->view[j];
  p->deltav = v1 * f1;
  p->ssq = f1*f1*q->ssq + f2*f2*r->ssq;
}  /* nuview */


void makecontrasts(node *p)
{
  /* compute the contrasts, recursively */
  if (p->tip)
    return;
  makecontrasts(p->next->back);
  makecontrasts(p->next->next->back);
  nuview(p);
  contbetween(p->next->back, p->next->next->back);
}  /* makecontrasts */


void writecontrasts()
{
  /* write out the contrasts */
  long i, j;

  if (printdata || reg) {
    fprintf(outfile, "\nContrasts (columns are different characters)\n");
    fprintf(outfile, "--------- -------- --- --------- -----------\n\n");
  }
  for (i = 0; i <= contno - 2; i++) {
    for (j = 0; j < chars; j++)
      fprintf(outfile, " %9.5f", cntrast[i][0][j]/sqrt(ssqcont[i][0]));
    putc('\n', outfile);
  }
}  /* writecontrasts */


void regressions()
{
  /* compute regressions and correlations among contrasts */
  long i, j, k;
  double **sumprod;

  sumprod = (double **)Malloc((long)chars*sizeof(double *));
  for (i = 0; i < chars; i++) {
    sumprod[i] = (double *)Malloc((long)chars*sizeof(double));
    for (j = 0; j < chars; j++)
      sumprod[i][j] = 0.0;
  }
    for (i = 0; i <= contno - 2; i++) {
    for (j = 0; j < chars; j++) {
      for (k = 0; k < chars; k++)
        sumprod[j][k] += cntrast[i][0][j] * cntrast[i][0][k] / ssqcont[i][0];
    }
  }
  fprintf(outfile, "\nCovariance matrix\n");
  fprintf(outfile, "---------- ------\n\n");
  for (i = 0; i < chars; i++) {
    for (j = 0; j < chars; j++)
      sumprod[i][j] /= contno - 1;
  }
  for (i = 0; i < chars; i++) {
    for (j = 0; j < chars; j++)
      fprintf(outfile, " %9.4f", sumprod[i][j]);
    putc('\n', outfile);
  }
  fprintf(outfile, "\nRegressions (columns on rows)\n");
  fprintf(outfile, "----------- -------- -- -----\n\n");
  for (i = 0; i < chars; i++) {
    for (j = 0; j < chars; j++)
      fprintf(outfile, " %9.4f", sumprod[i][j] / sumprod[i][i]);
    putc('\n', outfile);
  }
  fprintf(outfile, "\nCorrelations\n");
  fprintf(outfile, "------------\n\n");
  for (i = 0; i < chars; i++) {
    for (j = 0; j < chars; j++)
      fprintf(outfile, " %9.4f",
              sumprod[i][j] / sqrt(sumprod[i][i] * sumprod[j][j]));
    putc('\n', outfile);
  }
  for (i = 0; i < chars; i++)
    free(sumprod[i]);
  free(sumprod);
}  /* regressions */


double logdet(double **a)
{
  /* Gauss-Jordan log determinant calculation.
     in place, overwriting previous contents of a.  On exit,
      matrix a contains the inverse. Works only for positive definite A */
  long i, j, k;
  double temp, sum;

  sum = 0.0;
  for (i = 0; i < chars; i++) {
    if (fabs(a[i][i]) < 1.0E-37) {
       printf("ERROR: tried to invert singular matrix.\n");
       exxit(-1);
    }
    sum += log(a[i][i]);
    temp = 1.0 / a[i][i];
    a[i][i] = 1.0;
    for (j = 0; j < chars; j++)
      a[i][j] *= temp;
    for (j = 0; j < chars; j++) {
      if (j != i) {
        temp = a[j][i];
        a[j][i] = 0.0;
        for (k = 0; k < chars; k++)
          a[j][k] -= temp * a[i][k];
      }
    }
  }
  return(sum);
}  /* logdet */


void invert(double **a)
{
  /* Gauss-Jordan reduction -- invert chars x chars matrix a
     in place, overwriting previous contents of a.  On exit,
      matrix a contains the inverse.*/
  long i, j, k;
  double temp;

  for (i = 0; i < chars; i++) {
    if (fabs(a[i][i]) < 1.0E-37) {
       printf("ERROR: tried to invert singular matrix.\n");
       exxit(-1);
    }
    temp = 1.0 / a[i][i];
    a[i][i] = 1.0;
    for (j = 0; j < chars; j++)
      a[i][j] *= temp;
    for (j = 0; j < chars; j++) {
      if (j != i) {
        temp = a[j][i];
        a[j][i] = 0.0;
        for (k = 0; k < chars; k++)
          a[j][k] -= temp * a[i][k];
      }
    }
  }
}  /*invert*/


void initcovars(boolean novara)
{
/* Initialize covariance estimates */
   long i, j, k, l, contswithin;

   /* zero the matrices */
   for (i = 0; i < chars; i++)
     for (j = 0; j < chars; j++) {
       vara[i][j] = 0.0;
       vare[i][j] = 0.0;
   }
   /* estimate VE from within contrasts -- unbiasedly */
   contswithin = 0;
   for (i = 0; i < spp; i++) {
     for (j = 1; j < sample[i]; j++) {
       contswithin++;
       for (k = 0; k < chars; k++)
         for (l = 0; l < chars; l++)
           vare[k][l] += cntrast[i][j][k]*cntrast[i][j][l];
     }
   }
   /* estimate VA from between contrasts -- biasedly: does not take out VE */
   if (!novara) {   /* leave VarA = 0 if no A component assumed present */
     for (i = 0; i < spp-1; i++) {
       for (j = 0; j < chars; j++)
         for (k = 0; k < chars; k++)
           if (ssqcont[i][0] <= 0.0)
             vara[j][k] += cntrast[i][0][j]*cntrast[i][0][k];
           else
             vara[j][k] += cntrast[i][0][j]*cntrast[i][0][k]
                           / ((long)(spp-1)*ssqcont[i][0]);
     }
   }
   for (k = 0; k < chars; k++)
     for (l = 0; l < chars; l++)
       if (contswithin > 0)
         vare[k][l] /= contswithin;
       else {
         if (!novara) {
           vara[k][l] = 0.5 * vara[k][l];
           vare[k][l] = vara[k][l];
         }
       }
}  /* initcovars */


double normdiff(boolean novara)
{
/* Get relative norm of difference between old, new covariances */
   double s;
   long i, j;

   s = 0.0;
   for (i = 0; i < chars; i++)
     for (j = 0; j < chars; j++) {
        if (!novara) {
          if (fabs(oldvara[i][j]) <= 0.00000001)
            s += vara[i][j];
          else
            s += fabs(vara[i][j]/oldvara[i][j]-1.0);
        }
        if (fabs(oldvare[i][j]) <= 0.00000001)
          s += vare[i][j];
        else
          s += fabs(vare[i][j]/oldvare[i][j]-1.0);
     }
   return s/((double)(chars*chars));
}  /* normdiff */


void matcopy(double **a, double **b)
{
/* Copy matrices chars x chars: a to b */
   long i;
   
   for (i = 0; i < chars; i++) {
      memcpy(b[i], a[i], chars*sizeof(double));
   }
}  /* matcopy */


void newcovars(boolean novara)
{
/* one EM update of covariances, compute old likelihood too */
  long i, j, k, l, m;
  double sum, sum2, sum3, sqssq;

  if (!novara)
    matcopy(vara, oldvara); 
  matcopy(vare, oldvare); 
  sum2 = 0.0;         /* log likelihood of old parameters accumulates here */
  for (i = 0; i < chars; i++)                    /* zero out vara and vare */
    for (j = 0; j < chars; j++) {
      if (!novara)
        vara[i][j] = 0.0;
      vare[i][j] = 0.0;
    }
  for (i = 0; i < spp-1; i++) {            /* accumulate over contrasts ... */
    if (i <= spp-2) {       /* E(aa'|x) and E(ee'|x) for "between" contrasts */
      sqssq = sqrt(ssqcont[i][0]);
      for (k = 0; k < chars; k++)       /* compute (dA+E) for this contrast */
        for (l = 0; l < chars; l++)
          if (!novara)
            temp1[k][l] = ssqcont[i][0] * oldvara[k][l] + oldvare[k][l]; 
          else
            temp1[k][l] = oldvare[k][l]; 
      matcopy(temp1, temp2);
      invert(temp2);                               /* compute (dA+E)^(-1)  */
      /* sum of - x (dA+E)^(-1) x'/2 for old A, E */
      for (k = 0; k < chars; k++)
        for (l = 0; l < chars; l++)
          sum2 -= cntrast[i][0][k]*temp2[k][l]*cntrast[i][0][l]/2.0;
      matcopy(temp1, temp3);
      sum2 -= 0.5 * logdet(temp3);             /* log determinant term too */
      if (!novara) {
        for (k = 0; k < chars; k++)
          for (l = 0; l < chars; l++) {
            sum = 0.0;
            for (j = 0; j < chars; j++)
              sum += temp2[k][j] * sqssq * oldvara[j][l];
            Bax[k][l] = sum;           /*  Bax = (dA+E)^(-1) * sqrt(d) * A */
          }
      }
      for (k = 0; k < chars; k++)
        for (l = 0; l < chars; l++) {
          sum = 0.0;
          for (j = 0; j < chars; j++)
            sum += temp2[k][j] * oldvare[j][l];
          Bex[k][l] = sum;                       /*  Bex = (dA+E)^(-1) * E */
        }
      if (!novara) {
        for (k = 0; k < chars; k++)
          for (l = 0; l < chars; l++) {
             sum = 0.0;
             for (m = 0; m < chars; m++)
               sum += Bax[m][k] * (cntrast[i][0][m]*cntrast[i][0][l]
                                     -temp1[m][l]);
             temp2[k][l] = sum;                  /*  Bax'*(xx'-(dA+E)) ... */
          }
        for (k = 0; k < chars; k++)
          for (l = 0; l < chars; l++) {
             sum = 0.0;
             for (m = 0; m < chars; m++)
               sum += temp2[k][m] * Bax[m][l];
             vara[k][l] += sum;                          /*   ... * Bax  */
          }
      }
      for (k = 0; k < chars; k++)
        for (l = 0; l < chars; l++) {
           sum = 0.0;
           for (m = 0; m < chars; m++)
             sum += Bex[m][k] * (cntrast[i][0][m]*cntrast[i][0][l]
                                   -temp1[m][l]);
           temp2[k][l] = sum;                   /*  Bex'*(xx'-(dA+E)) ... */
        }
      for (k = 0; k < chars; k++)
        for (l = 0; l < chars; l++) {
           sum = 0.0;
           for (m = 0; m < chars; m++)
             sum += temp2[k][m] * Bex[m][l];
           vare[k][l] += sum;                            /*   ... * Bex  */
        }
    }
  }
  matcopy(oldvare, temp2);
  invert(temp2);                          /* get E^(-1) */
  matcopy(oldvare, temp3);
  sum3 = 0.5 * logdet(temp3);             /* get 1/2 log det(E) */
  for (i = 0; i < spp; i++) {
    if (sample[i] > 1) {
      for (j = 1; j < sample[i]; j++) {   /* E(aa'|x) (invisibly) and
                                             E(ee'|x) for within contrasts */
        for (k = 0; k < chars; k++)
          for (l = 0; l < chars; l++) {
            vare[k][l] += cntrast[i][j][k] * cntrast[i][j][l] - oldvare[k][l];
            sum2 -= cntrast[i][j][k] * temp2[k][l] * cntrast[i][j][l] / 2.0;
                                  /* accumulate - x*E^(-1)*x'/2 for old E */
          }
        sum2 -= sum3;                          /* log determinant term too */
      }
    }
  }
  for (i = 0; i < chars; i++)        /* complete EM by dividing by denom ... */
    for (j = 0; j < chars; j++) {    /* ... and adding old VA, VE */
      if (!novara) {
        vara[i][j] /= (double)contnum;
        vara[i][j] += oldvara[i][j];
      }
      vare[i][j] /= (double)contnum;
      vare[i][j] += oldvare[i][j];
    }
  logL = sum2;                       /* log likelihood for old values */
}  /* newcovars */


void printcovariances(boolean novara)
{
/* print out ML covariances and regressions in the error-covariance case */
   long i, j;

   fprintf(outfile, "\n\n");
   if (novara)
     fprintf(outfile, "Estimates when VarA is not in the model\n\n");
   else
     fprintf(outfile, "Estimates when VarA is in the model\n\n");
   if (!novara) {
     fprintf(outfile, "Estimate of VarA\n");
     fprintf(outfile, "-------- -- ----\n");
     fprintf(outfile, "\n");
     for (i = 0; i < chars; i++) {
       for (j = 0; j < chars; j++)
         fprintf(outfile, " %12.6f ", vara[i][j]);
       fprintf(outfile, "\n");
     }
     fprintf(outfile, "\n");
   }
   fprintf(outfile, "Estimate of VarE\n");
   fprintf(outfile, "-------- -- ----\n");
   fprintf(outfile, "\n");
   for (i = 0; i < chars; i++) {
     for (j = 0; j < chars; j++)
       fprintf(outfile, " %12.6f ", vare[i][j]);
     fprintf(outfile, "\n");
   }
   fprintf(outfile, "\n");
   if (!novara) {
     fprintf(outfile, "VarA Regressions (columns on rows)\n");
     fprintf(outfile, "---- ----------- -------- -- -----\n\n");
     for (i = 0; i < chars; i++) {
       for (j = 0; j < chars; j++)
         fprintf(outfile, " %9.4f", vara[i][j] / vara[i][i]);
       putc('\n', outfile);
     }
     fprintf(outfile, "\n");
     fprintf(outfile, "VarA Correlations\n");
     fprintf(outfile, "---- ------------\n\n");
     for (i = 0; i < chars; i++) {
       for (j = 0; j < chars; j++)
         fprintf(outfile, " %9.4f",
                 vara[i][j] / sqrt(vara[i][i] * vara[j][j]));
       putc('\n', outfile);
     }
     fprintf(outfile, "\n");
   }
   fprintf(outfile, "VarE Regressions (columns on rows)\n");
   fprintf(outfile, "---- ----------- -------- -- -----\n\n");
   for (i = 0; i < chars; i++) {
     for (j = 0; j < chars; j++)
       fprintf(outfile, " %9.4f", vare[i][j] / vare[i][i]);
     putc('\n', outfile);
   }
   fprintf(outfile, "\n");
   fprintf(outfile, "\nVarE Correlations\n");
   fprintf(outfile, "---- ------------\n\n");
   for (i = 0; i < chars; i++) {
     for (j = 0; j < chars; j++)
       fprintf(outfile, " %9.4f",
               vare[i][j] / sqrt(vare[i][i] * vare[j][j]));
     putc('\n', outfile);
   }
   fprintf(outfile, "\n\n");
} /* printcovariances */


void emiterate(boolean novara)
{
/* EM iteration of error and phylogenetic covariances */
   /* How to handle missing values? */
   long its;
   double relnorm;

   initcovars(novara);
   its = 1;
   do {
     newcovars(novara);
     relnorm = normdiff(novara);
     if (its % 100 == 0)
      printf("Iteration no. %ld:  ln L = %f, Norm = %f\n", its, logL, relnorm);
     its++;
   } while ((relnorm > 0.00001) && (its < 10000));
   if (its == 10000) {
     printf("\nWARNING: Iterations did not converge.");
     printf("  Results may be unreliable.\n");
   }
} /* emiterate */


void initcontrastnode(node **p, node **grbg, node *q, long len,
                        long nodei, long *ntips, long *parens, initops whichinit,
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
    nodep[(*p)->index - 1] = (*p);
    (*p)->view = (phenotype3)Malloc((long)chars*sizeof(double));
    break;
  case nonbottom:
    gnu(grbg, p);
    (*p)->index = nodei;
    (*p)->view = (phenotype3)Malloc((long)chars*sizeof(double));
    break;
  case tip:
    match_names_to_data (str, nodep, p, spp);
    (*p)->view = (phenotype3)Malloc((long)chars*sizeof(double));
    (*p)->deltav = 0.0;
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
  default:        /* cases of hslength,iter,hsnolength,treewt,unittrwt*/
    break;        /* not handled                                            */
  }
} /* initcontrastnode */


void maketree()
{
  /* set up the tree and use it */
  long which, nextnode;
  node *q, *r;

  alloctree(&curtree.nodep, nonodes);
  setuptree(&curtree, nonodes);
  which = 1;
  while (which <= numtrees) {
    if ((printdata || reg) && numtrees > 1) {
      fprintf(outfile, "Tree number%4ld\n", which);
      fprintf(outfile, "==== ====== ====\n\n");
    }
    nextnode = 0;
    treeread (intree, &curtree.start, curtree.nodep, &goteof, &first, 
                    curtree.nodep, &nextnode, &haslengths, &grbg, 
                    initcontrastnode,false,nonodes);
    q = curtree.start;
    r = curtree.start;
    while (!(q->next == curtree.start))
      q = q->next;
    q->next = curtree.start->next;
    curtree.start = q;
    chuck(&grbg, r);
    curtree.nodep[spp] = q;
    bifurcating = (curtree.start->next->next == curtree.start);
    contwithin();
    makecontrasts(curtree.start);
    if (!bifurcating) {
      makecontrasts(curtree.start->back);
      contbetween(curtree.start, curtree.start->back);
    }
    if (!varywithin) {
      if (writecont)
        writecontrasts();
      if (reg)
        regressions();
      putc('\n', outfile);
      }
    else {
      emiterate(false);
      printcovariances(false);
      if (nophylo) {
        logLvara = logL;
        emiterate(nophylo);
        printcovariances(nophylo);
        logLnovara = logL;
        fprintf(outfile, "\n\n\n    Likelihood Ratio Test");
        fprintf(outfile,  " of no VarA component\n");
        fprintf(outfile, "    ---------- ----- ----");
        fprintf(outfile,  " -- -- ---- ---------\n\n");
        fprintf(outfile, "      Log likelihood with varA     = %13.5f,",
                          logLvara);
        fprintf(outfile, "  %ld parameters\n\n", chars*(chars+1));
        fprintf(outfile, "      Log likelihood without varA  = %13.5f,",
                          logLnovara);
        fprintf(outfile, "  %ld parameters\n\n", chars*(chars+1)/2);
        fprintf(outfile, "                     difference    = %13.5f\n\n",
                          logLvara-logLnovara);
        fprintf(outfile, "      Chi-square value = %13.5f, ",
                  2.0*(logLvara-logLnovara));
        fprintf(outfile, "  %ld  degrees of freedom\n\n", chars*(chars+1)/2);
      }
    }
    which++;
  }
  if (progress)
    printf("\nOutput written to file \"%s\"\n\n", outfilename);
}  /* maketree */


int main(int argc, Char *argv[])
{  /* main program */
#ifdef MAC
  argc = 1;                /* macsetup("Contrast","Contrast");                */
  argv[0] = "Contrast";
#endif
  init(argc, argv);
  openfile(&infile,INFILE,"input data","r",argv[0],infilename);
  /* Open in binary: ftell() is broken for UNIX line-endings under WIN32 */
  openfile(&intree,INTREE,"input tree", "rb",argv[0],intreename);
  openfile(&outfile,OUTFILE,"output", "w",argv[0],outfilename);
  ibmpc = IBMCRT;
  ansi = ANSICRT;
  reg = true;
  numtrees = 1;
  doinit();
  getdata();
  maketree();
  FClose(infile);
  FClose(outfile);
  FClose(intree);
  printf("Done.\n\n");
#ifdef WIN32
  phyRestoreConsoleAttributes();
#endif
  return 0; 
}
