
#include "phylip.h"
#include "seq.h"

/* version 3.6. (c) Copyright 1994-2004 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define initialv        0.1     /* starting value of branch length          */
#define iterationsr     20      /* how many Newton iterations per distance  */


#ifndef OLDC
/* function prototypes */
void restdist_inputnumbers(void);
void getoptions(void);
void allocrest(void);
void doinit(void);
void inputoptions(void);
void restdist_inputdata(void);
void restdist_sitesort(void);
void restdist_sitecombine(void);
void makeweights(void);
void makev(long, long, double *);
void makedists(void);
void writedists(void);
void getinput(void);
void reallocsites(void);
/* function prototypes */
#endif


Char infilename[FNMLNGTH], outfilename[FNMLNGTH]; 
long sites, weightsum, datasets, ith;
boolean  restsites, neili, gama, weights, lower,
           progress, mulsets, firstset;
double ttratio, fracchange, cvi, sitelength, xi, xv;
double **d;
steptr aliasweight;
char *progname;
Char ch;

void restdist_inputnumbers()
{
  /* read and print out numbers of species and sites */
  fscanf(infile, "%ld%ld", &spp, &sites);
}  /* restdist_inputnumbers */


void getoptions()
{
  /* interactively set options */
  long loopcount, loopcount2;
  Char ch;

  putchar('\n');
  sitelength = 6.0;
  neili = false;
  gama = false;
  cvi = 0.0;
  weights = false;
  lower = false;
  printdata = false;
  progress = true;
  restsites = true;
  interleaved = true;
  ttratio = 2.0;
  loopcount = 0;
  for (;;) {
    cleerhome();
    printf("\nRestriction site or fragment distances, ");
    printf("version %s\n\n",VERSION);
    printf("Settings for this run:\n");
    printf("  R           Restriction sites or fragments?  %s\n",
           (restsites ? "Sites" : "Fragments"));
    printf("  N        Original or modified Nei/Li model?  %s\n",
           (neili ? "Original" : "Modified"));
    if (!neili) {
      printf("  G  Gamma distribution of rates among sites?");
      if (!gama)
        printf("  No\n");
      else
        printf("  Yes\n");
      printf("  T            Transition/transversion ratio?  %f\n", ttratio);
    }
    printf("  S                              Site length? %4.1f\n",sitelength);
    printf("  L                  Form of distance matrix?  %s\n",
           (lower ? "Lower-triangular" : "Square"));
    printf("  M               Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld sets\n", datasets);
    else
      printf("  No\n");
    printf("  I              Input sequences interleaved?  %s\n",
           (interleaved ? "Yes" : "No, sequential"));
    printf("  0       Terminal type (IBM PC, ANSI, none)?  %s\n",
           ibmpc ? "IBM PC" : ansi  ? "ANSI" : "(none)");
    printf("  1       Print out the data at start of run?  %s\n",
           (printdata ? "Yes" : "No"));
    printf("  2     Print indications of progress of run?  %s\n",
           (progress ? "Yes" : "No"));
    printf("\n  Y to accept these or type the letter for one to change\n");
    fflush(stdout);
    scanf("%c%*[^\n]", &ch);
    getchar();
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
    if (ch == 'Y')
      break;
    if (strchr("RDNGTSLMI012",ch) != NULL){
      switch (ch) {
        
      case 'R':
        restsites = !restsites;
        break;
        
      case 'G':
        if (!neili) 
          gama = !gama;
        break;

      case 'N':
        neili = !neili;
        break;

      case 'T':
        if (!neili)
          initratio(&ttratio);
        break;

      case 'S':
        loopcount2 = 0;
        do {
          printf("New Sitelength?\n");
          fflush(stdout);
          scanf("%lf%*[^\n]", &sitelength);
          getchar();
          if (sitelength < 1.0)
            printf("BAD RESTRICTION SITE LENGTH: %f\n", sitelength); 
          countup(&loopcount2, 10);
        } while (sitelength < 1.0);
        break;
        
      case 'L':
        lower = !lower;
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
        
      }
    } else
      printf("Not a possible option!\n");
    countup(&loopcount, 100);
  }
  if (gama) {
    loopcount = 0;
    do {
      printf(
"\nCoefficient of variation of substitution rate among sites (must be positive)?\n");
      fflush(stdout);
      scanf("%lf%*[^\n]", &cvi);
      getchar();
      countup(&loopcount, 100);
    } while (cvi <= 0.0);
    cvi = 1.0 / (cvi * cvi);
    printf("\n");
  }
  xi = (ttratio - 0.5)/(ttratio + 0.5);
  xv = 1.0 - xi;
  fracchange = xi*0.5 + xv*0.75;
}  /* getoptions */


void reallocsites() 
{
  long i;

  for (i = 0; i < spp; i++){
    free(y[i]);
    y[i] = (Char *)Malloc(sites*sizeof(Char));
  }
  
  free(weight);
  free(alias);
  free(aliasweight);

  weight = (steptr)Malloc((sites+1)*sizeof(long));
  alias = (steptr)Malloc((sites+1)*sizeof(long));
  aliasweight = (steptr)Malloc((sites+1)*sizeof(long));
  makeweights();
}

void allocrest()
{
  long i;

  y = (Char **)Malloc(spp*sizeof(Char *));
  for (i = 0; i < spp; i++)
    y[i] = (Char *)Malloc(sites*sizeof(Char));
  nayme = (naym *)Malloc(spp*sizeof(naym));
  weight = (steptr)Malloc((sites+1)*sizeof(long));
  alias = (steptr)Malloc((sites+1)*sizeof(long));
  aliasweight = (steptr)Malloc((sites+1)*sizeof(long));
  d = (double **)Malloc(spp*sizeof(double *));
  for (i = 0; i < spp; i++)
    d[i] = (double*)Malloc(spp*sizeof(double));
}  /* allocrest */


void doinit()
{
  /* initializes variables */
  restdist_inputnumbers();
  getoptions();
  if (printdata)
    fprintf(outfile, "\n %4ld Species, %4ld Sites\n", spp, sites);
  allocrest();
}  /* doinit */


void inputoptions()
{
  /* read the options information */
  Char ch;
  long i, extranum, cursp, curst;

  if (!firstset) {
    if (eoln(infile)) 
      scan_eoln(infile);
    fscanf(infile, "%ld%ld", &cursp, &curst);
    if (cursp != spp) {
      printf("\nERROR: INCONSISTENT NUMBER OF SPECIES IN DATA SET %4ld\n",
             ith);
      exxit(-1);
    }
    sites = curst;
    reallocsites();
  }
  for (i = 1; i <= sites; i++)
    weight[i] = 1;
  weightsum = sites;
  extranum = 0;
  fscanf(infile, "%*[ 0-9]");
  readoptions(&extranum, "W");
  for (i = 1; i <= extranum; i++) {
      matchoptions(&ch, "W");
      inputweights2(1, sites+1, &weightsum, weight, &weights, "RESTDIST");
  }
}  /* inputoptions */


void restdist_inputdata()
{
  /* read the species and sites data */
  long i, j, k, l, sitesread = 0, sitesnew = 0;
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
            printf(" ERROR -- Bad symbol %c",ch);
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
}  /* restdist_inputdata */


void restdist_sitesort()
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
}  /* restdist_sitesort */


void restdist_sitecombine()
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
}  /* restdist_sitecombine */


void makeweights()
{
  /* make up weights vector to avoid duplicate computations */
  long i;

  for (i = 1; i <= sites; i++) {
    alias[i] = i;
    aliasweight[i] = weight[i];
  }
  restdist_sitesort();
  restdist_sitecombine();
  sitescrunch2(sites + 1, 2, 3, aliasweight);
  for (i = 1; i <= sites; i++) {
    weight[i] = aliasweight[i];
    if (weight[i] > 0)
      endsite = i;
  }
  weight[0] = 1;
}  /* makeweights */


void makev(long m, long n, double *v)
{
  /* compute one distance */
  long i, ii, it, numerator, denominator;
  double f, g=0, h, p1, p2, p3, q1, pp, tt, delta, vv;

  numerator = 0;
  denominator = 0;
  for (i = 0; i < endsite; i++) {
    ii = alias[i + 1];
    if ((y[m-1][ii-1] == '+') || 
        (y[n-1][ii-1] == '+')) {
      denominator += weight[i + 1];
      if ((y[m-1][ii-1] == '+') && (y[n-1][ii-1] == '+')) {
        numerator += weight[i + 1];
      }
    }
  }
  f = 2*numerator/(double)(denominator+numerator);
  if (restsites) {
    if (exp(-sitelength*1.38629436) > f) {
      printf("\nERROR: Infinite distance between ");
      printf(" species %3ld and %3ld\n", m, n);
      exxit(-1);
    }
  }
  if (!restsites) {
    if (!neili) {
      f = (sqrt(f*(f+8.0))-f)/2.0;
    }
    else {
      g = initialv;
      delta = g;
      it = 0;
      while (fabs(delta) > 0.00002 && it < iterationsr) {
        it++;
        h = g;
        g = exp(0.25*log(f * (3-2*g)));
        delta = g - h;
      }
    }
  }
  if ((!restsites) && neili)
    vv = - (2.0/sitelength) * log(g);
  else {
    if (neili && restsites) {
      pp = exp(log(f)/(2*sitelength));
      vv = -(3.0/2.0)*log((4.0/3.0)*pp - (1.0/3.0));      
    } else {
      pp = exp(log(f)/sitelength); 
      delta = initialv;
      tt = delta;
      it = 0;
      while (fabs(delta) > 0.000001 && it < iterationsr) {
        it++;
        if (gama) {
          p1 = exp(-cvi * log(1 + tt / cvi));
          p2 = exp(-cvi * log(1 + xv * tt / cvi))
            - exp(-cvi * log(1 + tt / cvi));
          p3 = 1.0 - exp(-cvi * log(1 + xv * tt / cvi));
        } else {
          p1 = exp(-tt);
          p2 = exp(-xv * tt) - exp(-tt);
          p3 = 1.0 - exp(-xv * tt);
        }
        q1 = p1 + p2 / 2.0 + p3 / 4.0;
        g = q1 - pp;
        if (g < 0.0)
          delta = fabs(delta) / -2.0;
        else
          delta = fabs(delta);
        tt += delta;
      }
      vv = fracchange * tt;
    }
  }
  *v = fabs(vv);
}  /* makev */


void makedists()
{
  /* compute distance matrix */
  long i, j;
  double v;

  if (progress)
    printf("Distances calculated for species\n");
  for (i = 0; i < spp; i++)
    d[i][i] = 0.0;
  for (i = 1; i < spp; i++) {
    if (progress) {
      printf("    ");
      for (j = 0; j < nmlngth; j++)
        putchar(nayme[i - 1][j]);
      printf("   ");
    }
    for (j = i + 1; j <= spp; j++) {
      makev(i, j, &v);
      d[i - 1][j - 1] = v;
      d[j - 1][i - 1] = v;
      if (progress)
        putchar('.');
    }
    if (progress)
      putchar('\n');
  }
  if (progress) {
    printf("    ");
    for (j = 0; j < nmlngth; j++)
      putchar(nayme[spp - 1][j]);
    putchar('\n');
  }
}  /* makedists */


void writedists()
{
  /* write out distances */
  long i, j, k;

  if (!printdata)
    fprintf(outfile, "%5ld\n", spp);
  for (i = 0; i < spp; i++) {
    for (j = 0; j < nmlngth; j++)
      putc(nayme[i][j], outfile);
    if (lower)
      k = i;
    else
      k = spp;
    for (j = 1; j <= k; j++) {
      if (d[i][j-1] < 100.0)
        fprintf(outfile, "%10.6f", d[i][j-1]);
      else if (d[i][j-1] < 1000.0)
        fprintf(outfile, " %10.6f", d[i][j-1]);
        else
          fprintf(outfile, " %11.6f", d[i][j-1]);
      if ((j + 1) % 7 == 0 && j < k)
        putc('\n', outfile);
    }
    putc('\n', outfile);
  }
  if (progress)
    printf("\nDistances written to file \"%s\"\n\n", outfilename);
}  /* writedists */


void getinput()
{
  /* reads the input data */
  inputoptions();
  restdist_inputdata();
  makeweights();
}  /* getinput */


int main(int argc, Char *argv[])
{  /* distances from restriction sites or fragments */
#ifdef MAC
  argc = 1;                /* macsetup("Restdist",""); */
  argv[0] = "Restdist";
#endif
  init(argc,argv);
  progname = argv[0];
  openfile(&infile,INFILE,"input data file","r",argv[0],infilename);
  openfile(&outfile,OUTFILE,"output file","w",argv[0],outfilename);
  ibmpc = IBMCRT;
  ansi = ANSICRT;
  mulsets = false;
  datasets = 1;
  firstset = true;
  doinit();
  for (ith = 1; ith <= datasets; ith++) {
    getinput();
    if (ith == 1)
      firstset = false;
    if (datasets > 1 && progress)
      printf("\nData set # %ld:\n\n",ith);
    makedists();
    writedists();
  }
  FClose(infile);
  FClose(outfile);
#ifdef MAC
  fixmacfile(outfilename);
#endif

  printf("Done.\n\n");
#ifdef WIN32
  phyRestoreConsoleAttributes();
#endif
  return 0;
}  /* distances from restriction sites or fragments */
