#include "phylip.h"
#include "disc.h"

/* version 3.6. (c) Copyright 1993-2004 by the University of Washington.
   Written by Joseph Felsenstein, Jerry Shurman, Hisashi Horino,
   Akiko Fuseki, Sean Lamont, and Andrew Keeffe.  Permission is granted
   to copy and use this program provided no fee is charged for it and
   provided that this copyright notice is not removed. */

#define FormWide        80   /* width of outfile page */

typedef boolean *aPtr;
typedef long *SpPtr, *ChPtr;

typedef struct vecrec {
  aPtr vec;
  struct vecrec *next;
} vecrec;

typedef vecrec **aDataPtr;
typedef vecrec **Matrix;

#ifndef OLDC
/* function prototypes */
void clique_gnu(vecrec **);
void clique_chuck(vecrec *);
void nunode(node **);
void getoptions(void);
void clique_setuptree(void);
void allocrest(void);
void doinit(void);
void clique_inputancestors(void);
void clique_printancestors(void);
void clique_inputfactors(void);

void inputoptions(void);
void clique_inputdata(void);
boolean Compatible(long, long);
void SetUp(vecrec **);
void Intersect(boolean *, boolean *, boolean *);
long CountStates(boolean *);
void Gen1(long , long, boolean *, boolean *, boolean *);
boolean Ingroupstate(long );
void makeset(void);
void Init(long *, long *, long *, aPtr);

void ChSort(long *, long *, long);
void PrintClique(boolean *);
void bigsubset(long *, long);
void recontraverse(node **, long *, long, long);
void reconstruct(long, long);
void reroot(node *);
void clique_coordinates(node *, long *, long);
void clique_drawline(long);
void clique_printree(void);
void DoAll(boolean *, boolean *, boolean *, long);

void Gen2(long, long, boolean *, boolean *, boolean *);
void GetMaxCliques(vecrec **);
void reallocchars(void);
/* function prototypes */
#endif



Char infilename[FNMLNGTH], outfilename[FNMLNGTH], outtreename[FNMLNGTH], ancfilename[FNMLNGTH], factfilename[FNMLNGTH], weightfilename[FNMLNGTH];
long ActualChars, Cliqmin, outgrno,
       col, ith, msets, setsz;
boolean ancvar, Clmin, Factors, outgropt, trout, weights, noroot, justwts,
               printcomp, progress, treeprint, mulsets, firstset;
long nodes;

aPtr ancone;
Char *Factor;
long *ActChar, *oldweight;
aDataPtr Data;
Matrix Comp;            /* the character compatibility matrix      */
node *root;
long **grouping;
pointptr treenode;   /* pointers to all nodes in tree              */
vecrec *garbage;

/* these variables are to DoAll in the pascal Version. */
 aPtr aChars;
 boolean *Rarer;
 long n, MaxChars;
 SpPtr SpOrder;
 ChPtr ChOrder;

/* variables for GetMaxCliques: */
 vecrec **Comp2;
 long tcount;
 aPtr Temp, Processed, Rarer2;


void clique_gnu(vecrec **p)
{  /* this and the following are do-it-yourself garbage collectors.
     Make a new node or pull one off the garbage list */

 if (garbage != NULL) {
    *p = garbage;
    garbage = garbage->next;
  } else {
    *p = (vecrec *)Malloc((long)sizeof(vecrec));
    (*p)->vec = (aPtr)Malloc((long)chars*sizeof(boolean));
  }
  (*p)->next = NULL;
}  /* clique_gnu */


void clique_chuck(vecrec *p)
{  /* collect garbage on p -- put it on front of garbage list */

 p->next = garbage;
  garbage = p;
}  /* clique_chuck */


void nunode(node **p)
{  /* replacement for NEW */
  *p = (node *)Malloc((long)sizeof(node));
  (*p)->next = NULL;
  (*p)->tip = false;
}  /* nunode */


void getoptions(void)
{
  /* interactively set options */
  long loopcount, loopcount2;
  Char ch, ch2;
  boolean done;

  fprintf(outfile, "\nLargest clique program, version %s\n\n",VERSION);
  putchar('\n');
  ancvar = false;
  Clmin = false;
  Factors = false;
  outgrno = 1;
  outgropt = false;
  trout = true;
  weights = false;
  justwts = false;
  printdata = false;
  printcomp = false;
  progress = true;
  treeprint = true;
  loopcount = 0;
  do {
    cleerhome();
    printf("\nLargest clique program, version %s\n\n",VERSION);
    printf("Settings for this run:\n");
    printf("  A   Use ancestral states in input file?  %s\n",
           (ancvar ? "Yes" : "No"));
    printf("  F              Use factors information?  %s\n",
           Factors ? "Yes" : "No");
    printf("  W                       Sites weighted?  %s\n",
           (weights ? "Yes" : "No"));
    printf("  C          Specify minimum clique size?");
    if (Clmin)
      printf("  Yes, at size%3ld\n", Cliqmin);
    else
      printf("  No\n");
    printf("  O                        Outgroup root?  %s%3ld\n",
           (outgropt ? "Yes, at species number" :
                       "No, use as outgroup species"),outgrno);
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld %s\n", msets,
             (justwts ? "sets of weights" : "data sets"));
    else
      printf("  No\n");
    printf("  0   Terminal type (IBM PC, ANSI, none)?  %s\n",
           ibmpc ? "IBM PC" : ansi  ? "ANSI"   : "(none)");
    printf("  1    Print out the data at start of run  %s\n",
           (printdata ? "Yes" : "No"));
    printf("  2  Print indications of progress of run  %s\n",
           (progress ? "Yes" : "No"));
    printf("  3        Print out compatibility matrix  %s\n",
           (printcomp ? "Yes" : "No"));
    printf("  4                        Print out tree  %s\n",
           (treeprint ? "Yes" : "No"));
    printf("  5       Write out trees onto tree file?  %s\n",
           (trout ? "Yes" : "No"));
    if(weights && justwts){
        printf(
         "WARNING:  W option and Multiple Weights options are both on.  ");
        printf(
         "The W menu option is unnecessary and has no additional effect. \n");
    }
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
      if (strchr("OACMFW012345",ch) != NULL){
        switch (ch) {

        case 'A':
          ancvar = !ancvar;
          break;

        case 'F':
          Factors = !Factors;
          break;

        case 'W':
          weights = !weights;
          break;

        case 'C':
          Clmin = !Clmin;
          if (Clmin) {
            loopcount2 = 0;
            do {
              printf("Minimum clique size:\n");
#ifdef WIN32
              phyFillScreenColor();
#endif
              fflush(stdout);
              scanf("%ld%*[^\n]", &Cliqmin);
              getchar();
              countup(&loopcount2, 10);
            } while (Cliqmin < 0);
          }
          break;

        case 'O':
          outgropt = !outgropt;
          if (outgropt)
            initoutgroup(&outgrno, spp);
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
          printcomp = !printcomp;
          break;

        case '4':
          treeprint = !treeprint;
          break;

        case '5':
          trout = !trout;
          break;
        }
      } else
        printf("Not a possible option!\n");
      countup(&loopcount, 100);
    }
  } while (!done);
}  /* getoptions */


void clique_setuptree(void)
{
  /* initialization of tree pointers, variables */
  long i;

  treenode = (pointptr)Malloc((long)spp*sizeof(node *));
  for (i = 0; i < spp; i++) {
    treenode[i] = (node *)Malloc((long)sizeof(node));
    treenode[i]->next = NULL;
    treenode[i]->back = NULL;
    treenode[i]->index = i + 1;
    treenode[i]->tip = false;
  }
}  /* clique_setuptree */


void reallocchars()
{
  long i;
  Comp = (Matrix)Malloc((long)chars*sizeof(vecrec *));
  for (i = 0; i < (chars); i++)
    clique_gnu(&Comp[i]);
  ancone = (aPtr)Malloc((long)chars*sizeof(boolean));
  Factor = (Char *)Malloc((long)chars*sizeof(Char));
  ActChar = (long *)Malloc((long)chars*sizeof(long));
  oldweight = (long *)Malloc((long)chars*sizeof(long));
  weight = (long *)Malloc((long)chars*sizeof(long));
  ActualChars = chars;
  for (i = 1; i <= (chars); i++)
    ActChar[i - 1] = i;
}

void allocrest(void)
{
  long i;

  Data = (aDataPtr)Malloc((long)spp*sizeof(vecrec *));
  for (i = 0; i < (spp); i++)
    clique_gnu(&Data[i]);
  Comp = (Matrix)Malloc((long)chars*sizeof(vecrec *));
  for (i = 0; i < (chars); i++)
    clique_gnu(&Comp[i]);
  setsz = (long)ceil(((double)spp+1.0)/(double)SETBITS);
  ancone = (aPtr)Malloc((long)chars*sizeof(boolean));
  Factor = (Char *)Malloc((long)chars*sizeof(Char));
  ActChar = (long *)Malloc((long)chars*sizeof(long));
  oldweight = (long *)Malloc((long)chars*sizeof(long));
  weight = (long *)Malloc((long)chars*sizeof(long));
  nayme = (naym *)Malloc((long)spp*sizeof(naym));
}  /* allocrest */


void doinit(void)
{
  /* initializes variables */

  inputnumbersold(&spp, &chars, &nonodes, 1);
  getoptions();
  if (printdata)
    fprintf(outfile, "%2ld species, %3ld  characters\n", spp, chars);
  clique_setuptree();
  allocrest();
}  /* doinit */


void clique_inputancestors(void)
{
  /* reads the ancestral states for each character */
  long i;
  Char ch;

  for (i = 0; i < (chars); i++) {
    do {
      if (eoln(ancfile)) 
        scan_eoln(ancfile);
      ch = gettc(ancfile);
      if (ch == '\n')
        ch = ' ';
    } while (ch == ' ');
    switch (ch) {
    
    case '1':
      ancone[i] = true;
      break;

    case '0':
      ancone[i] = false;
      break;
    
    default:
      printf("BAD ANCESTOR STATE: %c AT CHARACTER %4ld\n", ch, i + 1);
      exxit(-1);
    }
  }
  scan_eoln(ancfile);
}  /* clique_inputancestors */


void clique_printancestors(void)
{
  /* print out list of ancestral states */
  long i;

  fprintf(outfile, "Ancestral states:\n");
  for (i = 1; i <= nmlngth + 2; i++)
    putc(' ', outfile);
  for (i = 1; i <= (chars); i++) {
    newline(outfile, i, 55, (long)nmlngth + 1);
    if (ancone[i - 1])
      putc('1', outfile);
    else
      putc('0', outfile);
    if (i % 5 == 0)
      putc(' ', outfile);
  }
  fprintf(outfile, "\n\n");
}  /* clique_printancestors */


void clique_inputfactors(void)
{
  /* reads the factor symbols */
  long i;

  ActualChars = 1;
  for (i = 1; i <= (chars); i++) {
    if (eoln(factfile)) 
      scan_eoln(factfile);
    Factor[i - 1] = gettc(factfile);
    if (i > 1) {
      if (Factor[i - 1] != Factor[i - 2])
        ActualChars++;
    }
    ActChar[i - 1] = ActualChars;
  }
  scan_eoln(factfile);
}  /* clique_inputfactors */


void inputoptions(void)
{
  /* reads the species names and character data */
  long i;
  if(justwts){
    if (!firstset)
      samenumsp(&chars, ith);
    if(firstset){
      ActualChars = chars;
      for (i = 1; i <= (chars); i++)
        ActChar[i - 1] = i;
      scan_eoln(infile);
    } else reallocchars();

    for (i = 0; i < (chars); i++)
      oldweight[i] = 1;
    inputweights(chars, oldweight, &weights);
    if(firstset && ancvar)
      clique_inputancestors();
    if(firstset && Factors)
      clique_inputfactors();
    if (printdata)
      printweights(outfile, 0, ActualChars, oldweight, "Characters");
    if (Factors)
      printfactors(outfile, chars, Factor, "");
    if (firstset && ancvar && printdata)
      clique_printancestors();
    noroot = !(outgropt || ancvar);
  } else {
    ActualChars = chars;
    for (i = 1; i <= (chars); i++)
      ActChar[i - 1] = i;
    for (i = 0; i < (chars); i++)
      oldweight[i] = 1;
    scan_eoln(infile);
    if(weights)
      inputweights(chars, oldweight, &weights);
    if(ancvar)
      clique_inputancestors();
    if(Factors)
      clique_inputfactors();
    if (weights && printdata)
      printweights(outfile, 0, ActualChars, oldweight, "Characters");
    if (Factors)
      printfactors(outfile, chars, Factor, "");
    if (ancvar && printdata)
      clique_printancestors();
    noroot = !(outgropt || ancvar);
  }
} /* inputoptions */


void clique_inputdata(void)
{
  long i, j;
  Char ch;

  j = chars / 2 + (chars / 5 - 1) / 2 - 5;
  if (j < 0)
    j = 0;
  if (j > 27)
    j = 27;
  if (printdata) {
    fprintf(outfile, "Species  ");
    for (i = 1; i <= j; i++)
      putc(' ', outfile);
    fprintf(outfile, "Character states\n");
    fprintf(outfile, "-------  ");
    for (i = 1; i <= j; i++)
      putc(' ', outfile);
    fprintf(outfile, "--------- ------\n\n");
  }
  for (i = 0; i < (spp); i++) {
    initname(i);
    if (printdata)
      for (j = 0; j < nmlngth; j++)
        putc(nayme[i][j], outfile);
    if (printdata)
      fprintf(outfile, "  ");
    for (j = 1; j <= (chars); j++) {
      do {
        if (eoln(infile)) 
          scan_eoln(infile);
        ch = gettc(infile);
      } while (ch == ' ' || ch == '\t');
      if (printdata) {
        putc(ch, outfile);
        newline(outfile, j, 55, (long)nmlngth + 1);
        if (j % 5 == 0)
          putc(' ', outfile);
      }
      if (ch != '0' && ch != '1') {
        printf("\n\nERROR: Bad character state: %c (not 0 or 1)", ch);
        printf(" at character %ld of species %ld\n\n", j, i + 1);
        exxit(-1);
      }
      Data[i]->vec[j - 1] = (ch == '1');
    }
    scan_eoln(infile);
    if (printdata)
      putc('\n', outfile);
  }
  putc('\n', outfile);
  for (i = 0; i < (chars); i++) {
    if (i + 1 == 1 || !Factors)
      weight[i] = oldweight[i];
    else if (Factor[i] != Factor[i - 1])
      weight[ActChar[i] - 1] = oldweight[i];
  }
}  /* clique_inputdata */


boolean Compatible(long ch1, long ch2)
{
  /* TRUE if two characters ch1 < ch2 are compatible */
  long i, j, k;
  boolean Compt, Done1, Done2;
  boolean Info[4];

  Compt = true;
  j = 1;
  while (ch1 > ActChar[j - 1])
    j++;
  Done1 = (ch1 != ActChar[j - 1]);
  while (!Done1) {
    k = j;
    while (ch2 > ActChar[k - 1])
      k++;
    Done2 = (ch2 != ActChar[k - 1]);
    while (!Done2) {
      for (i = 0; i <= 3; i++)
        Info[i] = false;
      if (ancvar) {
        if (ancone[j - 1] && ancone[k - 1])
          Info[0] = true;
        else if (ancone[j - 1] && !ancone[k - 1])
          Info[1] = true;
        else if (!ancone[j - 1] && ancone[k - 1])
          Info[2] = true;
        else
          Info[3] = true;
      }
      for (i = 0; i < (spp); i++) {
        if (Data[i]->vec[j - 1] && Data[i]->vec[k - 1])
          Info[0] = true;
        else if (Data[i]->vec[j - 1] && !Data[i]->vec[k - 1])
          Info[1] = true;
        else if (!Data[i]->vec[j - 1] && Data[i]->vec[k - 1])
          Info[2] = true;
        else
          Info[3] = true;
      }
      Compt = (Compt && !(Info[0] && Info[1] && Info[2] && Info[3]));
      k++;
      Done2 = (k > chars);
      if (!Done2)
        Done2 = (ch2 != ActChar[k - 1]);
    }
    j++;
    Done1 = (j > chars);
    if (!Done1)
      Done1 = (ch1 != ActChar[j - 1]);
  }
  return Compt;
}  /* Compatible */


void SetUp(vecrec **Comp)
{
  /* sets up the compatibility matrix */
  long i, j;

  if (printcomp) {
    if (Factors)
      fprintf(outfile, "      (For original multistate characters)\n");
    fprintf(outfile, "Character Compatibility Matrix (1 if compatible)\n");
    fprintf(outfile, "--------- ------------- ------ -- -- -----------\n\n");
  }
  for (i = 0; i < (ActualChars); i++) {
    if (printcomp) {
      for (j = 1; j <= ((48 - ActualChars) / 2); j++)
        putc(' ', outfile);
      for (j = 1; j < i + 1; j++) {
        if (Comp[i]->vec[j - 1])
          putc('1', outfile);
        else
          putc('.', outfile);
        newline(outfile, j, 70, (long)nmlngth + 1);
      }
    }
    Comp[i]->vec[i] = true;
    if (printcomp)
      putc('1', outfile);
    for (j = i + 1; j < (ActualChars); j++) {
      Comp[i]->vec[j] = Compatible(i + 1, j + 1);
      if (printcomp) {
        if (Comp[i]->vec[j])
          putc('1', outfile);
        else
          putc('.', outfile);
      }
      Comp[j]->vec[i] = Comp[i]->vec[j];
    }
    if (printcomp)
      putc('\n', outfile);
  }
  putc('\n', outfile);
}  /* SetUp */


void Intersect(boolean *V1, boolean *V2, boolean *V3)
{
  /* takes the logical intersection V1 AND V2 */
  long i;

  for (i = 0; i < (ActualChars); i++)
    V3[i] = (V1[i] && V2[i]);
}  /* Intersect */


long CountStates(boolean *V)
{
  /* counts the 1's in V */
  long i, TempCount;

  TempCount = 0;
  for (i = 0; i < (ActualChars); i++) {
    if (V[i])
      TempCount += weight[i];
  }
  return TempCount;
}  /* CountStates */


void Gen1(long i, long CurSize, boolean *aChars, boolean *Candidates,
                        boolean *Excluded)
{
  /* finds largest size cliques and prints them out */
  long CurSize2, j, k, Actual, Possible;
  boolean Futile;
  vecrec *Chars2, *Cands2, *Excl2, *Cprime, *Exprime;

  clique_gnu(&Chars2);
  clique_gnu(&Cands2);
  clique_gnu(&Excl2);
  clique_gnu(&Cprime);
  clique_gnu(&Exprime);
  CurSize2 = CurSize;
  memcpy(Chars2->vec, aChars, chars*sizeof(boolean));
  memcpy(Cands2->vec, Candidates, chars*sizeof(boolean));
  memcpy(Excl2->vec, Excluded, chars*sizeof(boolean));
  j = i;
  while (j <= ActualChars) {
    if (Cands2->vec[j - 1]) {
      Chars2->vec[j - 1] = true;
      Cands2->vec[j - 1] = false;
      CurSize2 += weight[j - 1];
      Possible = CountStates(Cands2->vec);
      Intersect(Cands2->vec, Comp2[j - 1]->vec, Cprime->vec);
      Actual = CountStates(Cprime->vec);
      Intersect(Excl2->vec, Comp2[j - 1]->vec, Exprime->vec);
      Futile = false;
      for (k = 0; k <= j - 2; k++) {
        if (Exprime->vec[k] && !Futile) {
          Intersect(Cprime->vec, Comp2[k]->vec, Temp);
          Futile = (CountStates(Temp) == Actual);
        }
      }
      if (CurSize2 + Actual >= Cliqmin && !Futile) {
        if (Actual > 0)
          Gen1(j + 1,CurSize2,Chars2->vec,Cprime->vec,Exprime->vec);
        else if (CurSize2 > Cliqmin) {
          Cliqmin = CurSize2;
          if (tcount >= 0)
            tcount = 1;
        } else if (CurSize2 == Cliqmin)
          tcount++;
      }
      if (Possible > Actual) {
        Chars2->vec[j - 1] = false;
        Excl2->vec[j - 1] = true;
        CurSize2 -= weight[j - 1];
      } else
        j = ActualChars;
    }
    j++;
  }
  clique_chuck(Chars2);
  clique_chuck(Cands2);
  clique_chuck(Excl2);
  clique_chuck(Cprime);
  clique_chuck(Exprime);
}  /* Gen1 */


boolean Ingroupstate(long i)
{
  /* the ingroup state for the i-th character */
  boolean outstate;

  if (noroot) {
    outstate = Data[0]->vec[i - 1];
    return (!outstate);
  }
  if (ancvar)
    outstate = ancone[i - 1];
  else
    outstate = Data[outgrno - 1]->vec[i - 1];
  return (!outstate);
}  /* Ingroupstate */


void makeset(void)
{
  /* make up set of species for given set of characters */
  long i, j, k, m;
  boolean instate;
  long *st;

  st = (long *)Malloc(setsz*sizeof(long));
  n = 0;
  for (i = 0; i < (MaxChars); i++) {
    for (j = 0; j < setsz; j++)
      st[j] = 0;
    instate = Ingroupstate(ChOrder[i]);
    for (j = 0; j < (spp); j++) {
      if (Data[SpOrder[j] - 1]->vec[ChOrder[i] - 1] == instate) {
        m = (long)(SpOrder[j]/SETBITS);
        st[m] = ((long)st[m]) | (1L << (SpOrder[j] % SETBITS));
      }
    }
    memcpy(grouping[++n - 1], st, setsz*sizeof(long));
  }
  for (i = 0; i < (spp); i++) {
    k = (long)(SpOrder[i]/SETBITS);
    grouping[++n - 1][k] = 1L << (SpOrder[i] % SETBITS);
  }
  free(st);
}  /* makeset */


void Init(long *ChOrder, long *Count, long *MaxChars, aPtr aChars)
{
  /* initialize vectors and character count */
  long i, j, temp;
  boolean instate;

  *MaxChars = 0;
  for (i = 1; i <= (chars); i++) {
    if (aChars[ActChar[i - 1] - 1]) {
      (*MaxChars)++;
      ChOrder[*MaxChars - 1] = i;
      instate = Ingroupstate(i);
      temp = 0;
      for (j = 0; j < (spp); j++) {
        if (Data[j]->vec[i - 1] == instate)
          temp++;
      }
      Count[i - 1] = temp;
    }
  }
}  /*Init */


void ChSort(long *ChOrder, long *Count, long MaxChars)
{
  /* sorts the characters by number of ingroup states */
  long j, temp;
  boolean ordered;

  ordered = false;
  while (!ordered) {
    ordered = true;
    for (j = 1; j < MaxChars; j++) {
      if (Count[ChOrder[j - 1] - 1] < Count[ChOrder[j] - 1]) {
        ordered = false;
        temp = ChOrder[j - 1];
        ChOrder[j - 1] = ChOrder[j];
        ChOrder[j] = temp;
      }
    }
  }
}  /* ChSort */


void PrintClique(boolean *aChars)
{
  /* prints the characters in a clique */
  long i, j;
  fprintf(outfile, "\n\n");
  if (Factors) {
    fprintf(outfile, "Actual Characters: (");
    j = 0;
    for (i = 1; i <= (ActualChars); i++) {
      if (aChars[i - 1]) {
        fprintf(outfile, "%3ld", i);
        j++;
        newline(outfile, j, (long)((FormWide - 22) / 3),
                                (long)nmlngth + 1);
      }
    }
    fprintf(outfile, ")\n");
  }
  if (Factors)
    fprintf(outfile, "Binary ");
  fprintf(outfile, "Characters: (");
  j = 0;
  for (i = 1; i <= (chars); i++) {
    if (aChars[ActChar[i - 1] - 1]) {
      fprintf(outfile, "%3ld", i);
      j++;
      if (Factors)
        newline(outfile, j, (long)((FormWide - 22) / 3), 
                                (long)nmlngth + 1);
      else
        newline(outfile, j, (long)((FormWide - 15) / 3), 
                                (long)nmlngth + 1);
    }
  }
  fprintf(outfile, ")\n\n");
}  /* PrintClique */


void bigsubset(long *st, long n)
{
  /* find a maximal subset of st among the groupings */
  long i, j;
  long *su;
  boolean max, same;

  su = (long *)Malloc(setsz*sizeof(long));
  for (i = 0; i < setsz; i++)
    su[i] = 0;
  for (i = 0; i < n; i++) {
    max = true;
    for (j = 0; j < setsz; j++)
      if ((grouping[i][j] & ~st[j]) != 0)
       max = false;
    if (max) {
      same = true;
      for (j = 0; j < setsz; j++)
        if (grouping[i][j] != st[j])
          same = false;
      if (!same) {
        for (j = 0; j < setsz; j++)
          if ((su[j] & ~grouping[i][j]) != 0)
            max = false;
        if (max) {
          same = true;
          for (j = 0; j < setsz; j++)
            if (grouping[i][j] != su[j])
              same = false;
          if (!same)
            memcpy(su, grouping[i], setsz*sizeof(long));
        }
      }
    }
  }
  memcpy(st, su, setsz*sizeof(long));
  free(su);
}  /* bigsubset */


void recontraverse(node **p, long *st, long n, long MaxChars)
{
  /* traverse to reconstruct the tree from the characters */
  long i, j, k, maxpos;
  long *tempset, *st2;
  boolean found, zero, zero2, same;
  node *q;

  j = k = 0;
  for (i = 1; i <= (spp); i++) {
    if (((1L << (i % SETBITS)) & st[(long)(i / SETBITS)]) != 0) {
      k++;
      j = i;
    }
  }
  if (k == 1) {
    *p = treenode[j - 1];
    (*p)->tip = true;
    (*p)->index = j;
    return;
  }
  nunode(p);
  (*p)->index = 0;
  tempset = (long*)Malloc(setsz*sizeof(long));
  memcpy(tempset, st, setsz*sizeof(long));
  q = *p;
  zero = true;
  for (i = 0; i < setsz; i++)
    if (tempset[i] != 0)
      zero = false;
  if (!zero)
    bigsubset(tempset, n);
  zero = true;
  zero2 = true;
  for (i = 0; i < setsz; i++)
    if (st[i] != 0)
      zero = false;
  if (!zero) {
    for (i = 0; i < setsz; i++)
      if (tempset[i] != 0)
        zero2 = false;
  }
  st2 = (long *)Malloc(setsz*sizeof(long));
  memcpy(st2, st, setsz*sizeof(long));
  while (!zero2) {
    nunode(&q->next);
    q = q->next;
    recontraverse(&q->back, tempset, n,MaxChars);
    i = 1;
    maxpos = 0;
    while (i <= MaxChars) {
      same = true;
      for (j = 0; j < setsz; j++)
        if (grouping[i - 1][j] != tempset[j])
          same = false;
      if (same)
        maxpos = i;
      i++;
    }
    q->back->maxpos = maxpos;
    q->back->back = q;
    for (j = 0; j < setsz; j++)
      st2[j] &= ~tempset[j];
    memcpy(tempset, st2, setsz*sizeof(long));
    found = false;
    i = 1;
    while (!found && i <= n) {
      same = true;
      for (j = 0; j < setsz; j++)
        if (grouping[i - 1][j] != tempset[j])
          same = false;
      if (same)
        found = true;
      else
        i++;
    }
    zero = true;
    for (j = 0; j < setsz; j++)
      if (tempset[j] != 0)
        zero = false;
    if (!zero && !found)
      bigsubset(tempset, n);
    zero = true;
    zero2 = true;
    for (j = 0; j < setsz; j++)
      if (st2[j] != 0)
        zero = false;
    if (!zero)
      for (j = 0; j < setsz; j++)
        if (tempset[j] != 0)
          zero2 = false;
}
  q->next = *p;
  free(tempset);
  free(st2);
}  /* recontraverse */


void reconstruct(long n, long MaxChars)
{  /* reconstruct tree from the subsets */
  long i;
  long *s;
  s = (long *)Malloc(setsz*sizeof(long));
  for (i = 0; i < setsz; i++) {
    if (i+1 == setsz) {
      s[i] = 1L << ((spp % SETBITS) + 1);
      if (setsz > 1)
        s[i] -= 1;
      else s[i] -= 1L << 1;
    }
    else if (i == 0) {
      if (setsz > 1)
        s[i] = ~0L - 1;
    }
    else {
      if (setsz > 2)
        s[i] = ~0L;
    }
  }
  recontraverse(&root,s,n,MaxChars);
  free(s);
}  /* reconstruct */


void reroot(node *outgroup)
{
  /* reorients tree, putting outgroup in desired position. */
  long i;
  boolean nroot;
  node *p, *q;

  nroot = false;
  p = root->next;
  while (p != root) {
    if (outgroup->back == p) {
      nroot = true;
      p = root;
    } else
      p = p->next;
  }
  if (nroot)
    return;
  p = root;
  i = 0;
  while (p->next != root) {
    p = p->next;
    i++;
  }
  if (i == 2) {
    root->next->back->back = p->back;
    p->back->back = root->next->back;
    q = root->next;
  } else {
    p->next = root->next;
    nunode(&root->next);
    q = root->next;
    nunode(&q->next);
    p = q->next;
    p->next = root;
    q->tip = false;
    p->tip = false;
  }
  q->back = outgroup;
  p->back = outgroup->back;
  outgroup->back->back = p;
  outgroup->back = q;
}  /* reroot */


void clique_coordinates(node *p, long *tipy, long MaxChars)
{
  /* establishes coordinates of nodes */
  node *q, *first, *last;
  long maxx;

  if (p->tip) {
    p->xcoord = 0;
    p->ycoord = *tipy;
    p->ymin = *tipy;
    p->ymax = *tipy;
    (*tipy) += down;
    return;
  }
  q = p->next;
  maxx = 0;
  while (q != p) {
    clique_coordinates(q->back, tipy, MaxChars);
    if (!q->back->tip) {
      if (q->back->xcoord > maxx)
        maxx = q->back->xcoord;
    }
    q = q->next;
  }
  first = p->next->back;
  q = p;
  while (q->next != p)
    q = q->next;
  last = q->back;
  p->xcoord = (MaxChars - p->maxpos) * 3 - 2;
  if (p == root)
    p->xcoord += 2;
  p->ycoord = (first->ycoord + last->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* clique_coordinates */


void clique_drawline(long i)
{
  /* draws one row of the tree diagram by moving up tree */
  node *p, *q;
  long n, m, j, k, l, sumlocpos, size, locpos, branchpos;
  long *poslist;
  boolean extra, done, plus, found, same;
  node *r, *first = NULL, *last = NULL;

  poslist = (long *)Malloc((long)(spp + MaxChars)*sizeof(long));
  branchpos = 0;
  p = root;
  q = root;
  fprintf(outfile, "  ");
  extra = false;
  plus = false;
  do {
    if (!p->tip) {
      found = false;
      r = p->next;
      while (r != p && !found) {
        if (i >= r->back->ymin && i <= r->back->ymax) {
          q = r->back;
          found = true;
        } else
          r = r->next;
      }
      first = p->next->back;
      r = p;
      while (r->next != p)
        r = r->next;
      last = r->back;
    }
    done = (p->tip || p == q);
    n = p->xcoord - q->xcoord;
    m = n;
    if (extra) {
      n--;
      extra = false;
    }
    if ((long)q->ycoord == i && !done) {
      if (!q->tip) {
        putc('+', outfile);
        plus = true;
        j = 1;
        for (k = 1; k <= (q->maxpos); k++) {
          same = true;
          for (l = 0; l < setsz; l++)
            if (grouping[k - 1][l] != grouping[q->maxpos - 1][l])
              same = false;
          if (same) {
            poslist[j - 1] = k;
            j++;
          }
        }
        size = j - 1;
        if (size == 0) {
          for (k = 1; k < n; k++)
            putc('-', outfile);
          sumlocpos = n;
        } else {
          sumlocpos = 0;
          j = 1;
          while (j <= size) {
            locpos = poslist[j - 1] * 3;
            if (j != 1)
              locpos -= poslist[j - 2] * 3;
            else
              locpos -= branchpos;
            for (k = 1; k < locpos; k++)
              putc('-', outfile);
            if (Rarer[ChOrder[poslist[j - 1] - 1] - 1])
              putc('1', outfile);
            else
              putc('0', outfile);
            sumlocpos += locpos;
            j++;
          }
          for (j = sumlocpos + 1; j < n; j++)
            putc('-', outfile);
          putc('+', outfile);
          if (m > 0)
            branchpos += m;
          extra = true;
        }
      } else {
        if (!plus) {
          putc('+', outfile);
          plus = false;
        } else
          n++;
        j = 1;
        for (k = 1; k <= (q->maxpos); k++) {
          same = true;
          for (l = 0; l < setsz; l++)
             if (grouping[k - 1][l] != grouping[q->maxpos - 1][l])
              same = false;
          if (same) {
            poslist[j - 1] = k;
            j++;
          }
        }
        size = j - 1;
        if (size == 0) {
          for (k = 1; k <= n; k++)
            putc('-', outfile);
          sumlocpos = n;
        } else {
          sumlocpos = 0;
          j = 1;
          while (j <= size) {
            locpos = poslist[j - 1] * 3;
            if (j != 1)
              locpos -= poslist[j - 2] * 3;
            else
              locpos -= branchpos;
            for (k = 1; k < locpos; k++)
              putc('-', outfile);
            if (Rarer[ChOrder[poslist[j - 1] - 1] - 1])
              putc('1', outfile);
            else
              putc('0', outfile);
            sumlocpos += locpos;
            j++;
          }
          for (j = sumlocpos + 1; j <= n; j++)
            putc('-', outfile);
          if (m > 0)
            branchpos += m;
        }
        putc('-', outfile);
      }
    } else if (!p->tip && (long)last->ycoord > i && (long)first->ycoord < i &&
               (i != (long)p->ycoord || p == root)) {
      putc('!', outfile);
      for (j = 1; j < n; j++)
        putc(' ', outfile);
      plus = false;
      if (m > 0)
        branchpos += m;
    } else {
      for (j = 1; j <= n; j++)
        putc(' ', outfile);
      plus = false;
      if (m > 0)
        branchpos += m;
    }
    if (q != p)
      p = q;
  } while (!done);
  if (p->ycoord == i && p->tip) {
    for (j = 0; j < nmlngth; j++)
      putc(nayme[p->index - 1][j], outfile);
}
  putc('\n', outfile);
  free(poslist);
}  /* clique_drawline */


void clique_printree(void)
{
  /* prints out diagram of the tree */
  long tipy, i;

  if (!treeprint)
    return;
  tipy = 1;
  clique_coordinates(root, &tipy, MaxChars);
  fprintf(outfile, "\n  Tree and");
  if (Factors)
    fprintf(outfile, " binary");
  fprintf(outfile, " characters:\n\n");
  fprintf(outfile, "   ");
  for (i = 0; i < (MaxChars); i++)
    fprintf(outfile, "%3ld", ChOrder[i]);
  fprintf(outfile, "\n   ");
  for (i = 0; i < (MaxChars); i++) {
    if (Rarer[ChOrder[i] - 1])
      fprintf(outfile, "%3c", '1');
    else
      fprintf(outfile, "%3c", '0');
  }
  fprintf(outfile, "\n\n");
  for (i = 1; i <= (tipy - down); i++)
    clique_drawline(i);
  fprintf(outfile, "\nremember: this is an unrooted tree!\n\n");
}  /* clique_printree */


void DoAll(boolean *Chars_,boolean *Processed,boolean *Rarer_,long tcount)
{
  /* print out a clique and its tree */
  long i, j;
  ChPtr Count;

  aChars = (aPtr)Malloc((long)chars*sizeof(boolean));
  SpOrder = (SpPtr)Malloc((long)spp*sizeof(long));
  ChOrder = (ChPtr)Malloc((long)chars*sizeof(long));
  Count = (ChPtr)Malloc((long)chars*sizeof(long));
  memcpy(aChars, Chars_, chars*sizeof(boolean));
  Rarer = Rarer_;
  Init(ChOrder, Count, &MaxChars, aChars);
  ChSort(ChOrder, Count, MaxChars);
  for (i = 1; i <= (spp); i++)
    SpOrder[i - 1] = i;
  for (i = 1; i <= (chars); i++) {
    if (aChars[ActChar[i - 1] - 1]) {
      if (!Processed[ActChar[i - 1] - 1]) {
        Rarer[i - 1] = Ingroupstate(i);
        Processed[ActChar[i - 1] - 1] = true;
      }
    }
  }
  PrintClique(aChars);
  grouping = (long **)Malloc((long)(spp + MaxChars)*sizeof(long *));
  for (i = 0; i < spp + MaxChars; i++) {
    grouping[i] = (long *)Malloc(setsz*sizeof(long));
    for (j = 0; j < setsz; j++)
      grouping[i][j] = 0;
  }
  makeset();
  clique_setuptree();
  reconstruct(n,MaxChars);
  if (noroot)
    reroot(treenode[outgrno - 1]);
  clique_printree();
  if (trout) {
    col = 0;
    treeout(root, tcount+1, &col, root);
  }
  free(SpOrder);
  free(ChOrder);
  free(Count);
  for (i = 0; i < spp + MaxChars; i++)
    free(grouping[i]);
  free(grouping);
}  /* DoAll */


void Gen2(long i, long CurSize, boolean *aChars, boolean *Candidates,
                        boolean *Excluded)
{
  /* finds largest size cliques and prints them out */
  long CurSize2, j, k, Actual, Possible;
  boolean Futile;
  vecrec *Chars2, *Cands2, *Excl2, *Cprime, *Exprime;

  clique_gnu(&Chars2);
  clique_gnu(&Cands2);
  clique_gnu(&Excl2);
  clique_gnu(&Cprime);
  clique_gnu(&Exprime);
  CurSize2 = CurSize;
  memcpy(Chars2->vec, aChars, chars*sizeof(boolean));
  memcpy(Cands2->vec, Candidates, chars*sizeof(boolean));
  memcpy(Excl2->vec, Excluded, chars*sizeof(boolean));
  j = i;
  while (j <= ActualChars) {
    if (Cands2->vec[j - 1]) {
      Chars2->vec[j - 1] = true;
      Cands2->vec[j - 1] = false;
      CurSize2 += weight[j - 1];
      Possible = CountStates(Cands2->vec);
      Intersect(Cands2->vec, Comp2[j - 1]->vec, Cprime->vec);
      Actual = CountStates(Cprime->vec);
      Intersect(Excl2->vec, Comp2[j - 1]->vec, Exprime->vec);
      Futile = false;
      for (k = 0; k <= j - 2; k++) {
        if (Exprime->vec[k] && !Futile) {
          Intersect(Cprime->vec, Comp2[k]->vec, Temp);
          Futile = (CountStates(Temp) == Actual);
        }
      }
      if (CurSize2 + Actual >= Cliqmin && !Futile) {
        if (Actual > 0)
          Gen2(j + 1,CurSize2,Chars2->vec,Cprime->vec,Exprime->vec);
        else
          DoAll(Chars2->vec,Processed,Rarer2,tcount);
      }
      if (Possible > Actual) {
        Chars2->vec[j - 1] = false;
        Excl2->vec[j - 1] = true;
        CurSize2 -= weight[j - 1];
      } else
        j = ActualChars;
    }
    j++;
  }
  clique_chuck(Chars2);
  clique_chuck(Cands2);
  clique_chuck(Excl2);
  clique_chuck(Cprime);
  clique_chuck(Exprime);
}  /* Gen2 */


void GetMaxCliques(vecrec **Comp_)
{
  /* recursively generates the largest cliques */
  long i;
  aPtr aChars, Candidates, Excluded;

  Temp = (aPtr)Malloc((long)chars*sizeof(boolean));
  Processed = (aPtr)Malloc((long)chars*sizeof(boolean));
  Rarer2 = (aPtr)Malloc((long)chars*sizeof(boolean));
  aChars = (aPtr)Malloc((long)chars*sizeof(boolean));
  Candidates = (aPtr)Malloc((long)chars*sizeof(boolean));
  Excluded = (aPtr)Malloc((long)chars*sizeof(boolean));
  Comp2 = Comp_;
  putc('\n', outfile);
  if (Clmin) {
    fprintf(outfile, "Cliques with at least%3ld characters\n", Cliqmin);
    fprintf(outfile, "------- ---- -- ----- -- ----------\n");
  } else {
    Cliqmin = 0;
    fprintf(outfile, "Largest Cliques\n");
    fprintf(outfile, "------- -------\n");
    for (i = 0; i < (ActualChars); i++) {
      aChars[i] = false;
      Excluded[i] = false;
      Candidates[i] = true;
    }
    tcount = 0;
    Gen1(1, 0, aChars, Candidates, Excluded);
  }
  for (i = 0; i < (ActualChars); i++) {
    aChars[i] = false;
    Candidates[i] = true;
    Processed[i] = false;
    Excluded[i] = false;
  }
  Gen2(1, 0, aChars, Candidates, Excluded);
  putc('\n', outfile);
  free(Temp);
  free(Processed);
  free(Rarer2);
  free(aChars);
  free(Candidates);
  free(Excluded);
}  /* GetMaxCliques */


int main(int argc, Char *argv[])
{  /* Main Program */
#ifdef MAC
   argc = 1;                /* macsetup("Clique","Clique");                */
   argv[0] = "Clique";
#endif
  init(argc, argv);
  openfile(&infile,INFILE,"input file", "r",argv[0],infilename);
  openfile(&outfile,OUTFILE,"output file", "w",argv[0],outfilename);
  openfile(&outtree,OUTTREE,"output tree file", "w",argv[0],outtreename);
  ibmpc = IBMCRT;
  ansi = ANSICRT;
  mulsets = false;
  firstset = true;
  msets = 1;
  doinit();
  if(ancvar)
    openfile(&ancfile,ANCFILE,"ancestors file", "r",argv[0],ancfilename);
  if(Factors)
      openfile(&factfile,FACTFILE,"factors file", "r",argv[0],factfilename);
  if (weights || justwts)
    openfile(&weightfile,WEIGHTFILE,"weights file","r",argv[0],weightfilename);
  for (ith = 1; ith <= (msets); ith++) {
    inputoptions();
    if(!justwts || firstset)
      clique_inputdata();
    firstset = false;
    SetUp(Comp);
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
    GetMaxCliques(Comp);
    if (progress) {
      printf("\nOutput written to file \"%s\"\n",outfilename);
      if (trout)
        printf("\nTree");
        if (tcount > 1)
          printf("s");
        printf(" written on file \"%s\"\n\n", outtreename);
    }
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
}

