
#include "phylip.h"
#include "seq.h"

/* version 3.6. (c) Copyright 1993-2004 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define maxtrees        100   /* maximum number of tied trees stored */

typedef enum {
  universal, ciliate, mito, vertmito, flymito, yeastmito
} codetype;

/* nodes will form a binary tree */

typedef struct gseq {
  seqptr seq;
  struct gseq *next;
} gseq;

#ifndef OLDC
/* function prototypes */
void   protgnu(gseq **);
void   protchuck(gseq *);
void   code(void);
void   setup(void);
void   getoptions(void);
void   protalloctree(void);
void   allocrest(void);
void   doinit(void);
void   protinputdata(void);

void   protmakevalues(void);
void   doinput(void);
void   protfillin(node *, node *, node *);
void   protpreorder(node *);
void   protadd(node *, node *, node *);
void   protre_move(node **, node **);
void   evaluate(node *);
void   protpostorder(node *);
void   protreroot(node *);
void   protsavetraverse(node *, long *, boolean *);

void   protsavetree(long *, boolean *);
void   tryadd(node *, node **, node **);
void   addpreorder(node *, node *, node *);
void   tryrearr(node *, boolean *);
void   repreorder(node *, boolean *);
void   rearrange(node **);
void   protgetch(Char *);
void   protaddelement(node **, long *, long *, boolean *);
void   prottreeread(void);
void   protancestset(long *, long *, long *, long *, long *);
                
void   prothyprint(long , long , boolean *, node *, boolean *, boolean *);
void   prothyptrav(node *, sitearray *, long, long, long *, boolean *,
                sitearray);
void   prothypstates(long *);
void   describe(void);
void   maketree(void);
void   reallocnode(node* p);
void   reallocchars(void);
/* function prototypes */
#endif



Char infilename[FNMLNGTH], outfilename[FNMLNGTH], intreename[FNMLNGTH], outtreename[FNMLNGTH], weightfilename[FNMLNGTH];
node *root;
long chars, col, msets, ith, njumble, jumb;
/*   chars = number of sites in actual sequences */
long inseed, inseed0;
boolean jumble, usertree, weights, thresh, trout, progress, stepbox,
    justwts, ancseq, mulsets, firstset;
codetype whichcode;
long fullset, fulldel;
pointarray treenode;   /* pointers to all nodes in tree */
double threshold;
steptr threshwt;
longer seed;
long *enterorder;
sitearray translate[(long)quest - (long)ala + 1];
aas trans[4][4][4];
long **fsteps;
bestelm *bestrees;
boolean dummy;
gseq *garbage;
node *temp, *temp1;
Char ch;
aas tmpa;
char *progname;

/* Local variables for maketree, propagated globally for c version: */
long minwhich;
double like, bestyet, bestlike, minsteps, bstlike2;
boolean lastrearr, recompute;
node *there;
double nsteps[maxuser];
long *place;
boolean *names;


void protgnu(gseq **p)
{
  /* this and the following are do-it-yourself garbage collectors.
     Make a new node or pull one off the garbage list */
  if (garbage != NULL) {
    *p = garbage;
    free((*p)->seq);
    (*p)->seq = (seqptr)Malloc(chars*sizeof(sitearray));
    garbage = garbage->next;
  } else {
    *p = (gseq *)Malloc(sizeof(gseq));
    (*p)->seq = (seqptr)Malloc(chars*sizeof(sitearray));
  }
  (*p)->next = NULL;
}  /* protgnu */


void protchuck(gseq *p)
{
  /* collect garbage on p -- put it on front of garbage list */
  p->next = garbage;
  garbage = p;
}  /* protchuck */


void code()
{
  /* make up table of the code 1 = u, 2 = c, 3 = a, 4 = g */
  trans[0][0][0] = phe;
  trans[0][0][1] = phe;
  trans[0][0][2] = leu;
  trans[0][0][3] = leu;
  trans[0][1][0] = ser1;
  trans[0][1][1] = ser1;
  trans[0][1][2] = ser1;
  trans[0][1][3] = ser1;
  trans[0][2][0] = tyr;
  trans[0][2][1] = tyr;
  trans[0][2][2] = stop;
  trans[0][2][3] = stop;
  trans[0][3][0] = cys;
  trans[0][3][1] = cys;
  trans[0][3][2] = stop;
  trans[0][3][3] = trp;
  trans[1][0][0] = leu;
  trans[1][0][1] = leu;
  trans[1][0][2] = leu;
  trans[1][0][3] = leu;
  trans[1][1][0] = pro;
  trans[1][1][1] = pro;
  trans[1][1][2] = pro;
  trans[1][1][3] = pro;
  trans[1][2][0] = his;
  trans[1][2][1] = his;
  trans[1][2][2] = gln;
  trans[1][2][3] = gln;
  trans[1][3][0] = arg;
  trans[1][3][1] = arg;
  trans[1][3][2] = arg;
  trans[1][3][3] = arg;
  trans[2][0][0] = ileu;
  trans[2][0][1] = ileu;
  trans[2][0][2] = ileu;
  trans[2][0][3] = met;
  trans[2][1][0] = thr;
  trans[2][1][1] = thr;
  trans[2][1][2] = thr;
  trans[2][1][3] = thr;
  trans[2][2][0] = asn;
  trans[2][2][1] = asn;
  trans[2][2][2] = lys;
  trans[2][2][3] = lys;
  trans[2][3][0] = ser2;
  trans[2][3][1] = ser2;
  trans[2][3][2] = arg;
  trans[2][3][3] = arg;
  trans[3][0][0] = val;
  trans[3][0][1] = val;
  trans[3][0][2] = val;
  trans[3][0][3] = val;
  trans[3][1][0] = ala;
  trans[3][1][1] = ala;
  trans[3][1][2] = ala;
  trans[3][1][3] = ala;
  trans[3][2][0] = asp;
  trans[3][2][1] = asp;
  trans[3][2][2] = glu;
  trans[3][2][3] = glu;
  trans[3][3][0] = gly;
  trans[3][3][1] = gly;
  trans[3][3][2] = gly;
  trans[3][3][3] = gly;
  if (whichcode == mito)
    trans[0][3][2] = trp;
  if (whichcode == vertmito) {
    trans[0][3][2] = trp;
    trans[2][3][2] = stop;
    trans[2][3][3] = stop;
    trans[2][0][2] = met;
  }
  if (whichcode == flymito) {
    trans[0][3][2] = trp;
    trans[2][0][2] = met;
    trans[2][3][2] = ser2;
  }
  if (whichcode == yeastmito) {
    trans[0][3][2] = trp;
    trans[1][0][2] = thr;
    trans[2][0][2] = met;
  }
} /* code */


void setup()
{
  /* set up set table to get aasets from aas */
  aas a, b;
  long i, j, k, l, s;

  for (a = ala; (long)a <= (long)stop; a = (aas)((long)a + 1)) {
    translate[(long)a - (long)ala][0] = 1L << ((long)a);
    translate[(long)a - (long)ala][1] = 1L << ((long)a);
  }
  for (i = 0; i <= 3; i++) {
    for (j = 0; j <= 3; j++) {
      for (k = 0; k <= 3; k++) {
        for (l = 0; l <= 3; l++) {
          translate[(long)trans[i][j][k]][1] |= (1L << (long)trans[l][j][k]);
          translate[(long)trans[i][j][k]][1] |= (1L << (long)trans[i][l][k]);
          translate[(long)trans[i][j][k]][1] |= (1L << (long)trans[i][j][l]);
        }
      }
    }
  }
  translate[(long)del - (long)ala][1] = 1L << ((long)del);
  fulldel = (1L << ((long)stop + 1)) - (1L << ((long)ala));
  fullset = fulldel & (~(1L << ((long)del)));
  translate[(long)asx - (long)ala][0]
    = (1L << ((long)asn)) | (1L << ((long)asp));
  translate[(long)glx - (long)ala][0]
    = (1L << ((long)gln)) | (1L << ((long)glu));
  translate[(long)ser - (long)ala][0]
    = (1L << ((long)ser1)) | (1L << ((long)ser2));
  translate[(long)unk - (long)ala][0] = fullset;
  translate[(long)quest - (long)ala][0] = fulldel;
  translate[(long)asx - (long)ala][1] = translate[(long)asn - (long)ala][1]
                                       | translate[(long)asp - (long)ala][1];
  translate[(long)glx - (long)ala][1] = translate[(long)gln - (long)ala][1]
                                       | translate[(long)glu - (long)ala][1];
  translate[(long)ser - (long)ala][1] = translate[(long)ser1 - (long)ala][1]
                                       | translate[(long)ser2 - (long)ala][1];
  translate[(long)unk - (long)ala][1] = fullset;
  translate[(long)quest - (long)ala][1] = fulldel;
  for (a = ala; (long)a <= (long)quest; a = (aas)((long)a + 1)) {
    s = 0;
    for (b = ala; (long)b <= (long)stop; b = (aas)((long)b + 1)) {
      if (((1L << ((long)b)) & translate[(long)a - (long)ala][1]) != 0)
        s |= translate[(long)b - (long)ala][1];
    }
    translate[(long)a - (long)ala][2] = s;
  }
}  /* setup */


void getoptions()
{
  /* interactively set options */
  long loopcount, loopcount2;
  Char ch, ch2;

  fprintf(outfile, "\nProtein parsimony algorithm, version %s\n\n",VERSION);
  putchar('\n');
  jumble = false;
  njumble = 1;
  outgrno = 1;
  outgropt = false;
  thresh = false;
  trout = true;
  usertree = false;
  weights = false;
  whichcode = universal;
  printdata = false;
  progress = true;
  treeprint = true;
  stepbox = false;
  ancseq = false;
  dotdiff = true;
  interleaved = true;
  loopcount = 0;
  for (;;) {
    cleerhome();
    printf("\nProtein parsimony algorithm, version %s\n\n",VERSION);
    printf("Setting for this run:\n");
    printf("  U                 Search for best tree?  %s\n",
           (usertree ? "No, use user trees in input file" : "Yes"));
    if (!usertree) {
      printf("  J   Randomize input order of sequences?");
      if (jumble)
        printf("  Yes (seed =%8ld,%3ld times)\n", inseed0, njumble);
      else
        printf("  No. Use input order\n");
    }
    printf("  O                        Outgroup root?");
    if (outgropt)
      printf("  Yes, at sequence number%3ld\n", outgrno);
    else
      printf("  No, use as outgroup species%3ld\n", outgrno);
    printf("  T              Use Threshold parsimony?");
    if (thresh)
      printf("  Yes, count steps up to%4.1f per site\n", threshold);
    else
      printf("  No, use ordinary parsimony\n");
    printf("  C               Use which genetic code?  %s\n",
      (whichcode == universal) ? "Universal"                  :
      (whichcode == ciliate)   ? "Ciliate"                    :
      (whichcode == mito)      ? "Universal mitochondrial"    :
      (whichcode == vertmito)  ? "Vertebrate mitochondrial"   :
      (whichcode == flymito)   ? "Fly mitochondrial"          :
      (whichcode == yeastmito) ? "Yeast mitochondrial"        : "");
    printf("  W                       Sites weighted?  %s\n",
           (weights ? "Yes" : "No"));
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
        printf("  Yes, %2ld %s\n", msets,
               (justwts ? "sets of weights" : "data sets"));
    else
      printf("  No\n");
    printf("  I          Input sequences interleaved?  %s\n",
           (interleaved ? "Yes" : "No"));
    printf("  0   Terminal type (IBM PC, ANSI, none)?  %s\n",
           (ibmpc ? "IBM PC" : ansi ? "ANSI" : "(none)"));
    printf("  1    Print out the data at start of run  %s\n",
           (printdata ? "Yes" : "No"));
    printf("  2  Print indications of progress of run  %s\n",
           (progress ? "Yes" : "No"));
    printf("  3                        Print out tree  %s\n",
           (treeprint ? "Yes" : "No"));
    printf("  4          Print out steps in each site  %s\n",
           (stepbox ? "Yes" : "No"));
    printf("  5  Print sequences at all nodes of tree  %s\n",
           (ancseq ? "Yes" : "No"));
    if (ancseq || printdata)
      printf("  .  Use dot-differencing to display them  %s\n",
              dotdiff ? "Yes" : "No");
    printf("  6       Write out trees onto tree file?  %s\n",
           (trout ? "Yes" : "No"));
    if(weights && justwts){
        printf(
         "WARNING:  W option and Multiple Weights options are both on.  ");
        printf(
         "The W menu option is unnecessary and has no additional effect. \n");
    }
    printf(
   "\nAre these settings correct? (type Y or the letter for one to change)\n");
    fflush(stdout);
    scanf("%c%*[^\n]", &ch);
    getchar();
    uppercase(&ch);
    if (ch == 'Y')
      break;
    if (((!usertree) && (strchr("WCJOTUMI12345.60", ch) != NULL))
        || (usertree && ((strchr("WCOTUMI12345.60", ch) != NULL)))){
      switch (ch) {
        
      case 'J':
        jumble = !jumble;
        if (jumble)
          initjumble(&inseed, &inseed0, seed, &njumble);
        else njumble = 1;
        break;
      
      case 'W':
        weights = !weights;
        break;
  
      case 'O':
        outgropt = !outgropt;
        if (outgropt)
          initoutgroup(&outgrno, spp);
        else outgrno = 1;
        break;
        
      case 'T':
        thresh = !thresh;
        if (thresh)
          initthreshold(&threshold);
        break;

      case 'C':
        printf("\nWhich genetic code?\n");
        printf(" type         for\n\n");
        printf("   U           Universal\n");
        printf("   M           Mitochondrial\n");
        printf("   V           Vertebrate mitochondrial\n");
        printf("   F           Fly mitochondrial\n");
        printf("   Y           Yeast mitochondrial\n\n");
        loopcount2 = 0;
        do {
          printf("type U, M, V, F, or Y\n");
          fflush(stdout);
          scanf("%c%*[^\n]", &ch);
          getchar();
          if (ch == '\n')
            ch = ' ';
          uppercase(&ch);
          countup(&loopcount2, 10);
        } while (ch != 'U' && ch != 'M' && ch != 'V'
                  && ch != 'F' && ch != 'Y');
        switch (ch) {

        case 'U':
          whichcode = universal;
          break;

        case 'M':
          whichcode = mito;
          break;

        case 'V':
          whichcode = vertmito;
          break;

        case 'F':
          whichcode = flymito;
          break;

        case 'Y':
          whichcode = yeastmito;
          break;
        }
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
        
      case 'U':
        usertree = !usertree;
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

      case '.':
        dotdiff = !dotdiff;
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


void protalloctree()
{ /* allocate treenode dynamically */
  long i, j;
  node *p, *q;

  treenode = (pointarray)Malloc(nonodes*sizeof(node *));
  for (i = 0; i < (spp); i++) {
    treenode[i] = (node *)Malloc(sizeof(node));
    treenode[i]->numsteps = (steptr)Malloc(chars*sizeof(long));
    treenode[i]->siteset = (seqptr)Malloc(chars*sizeof(sitearray));
    treenode[i]->seq = (aas *)Malloc(chars*sizeof(aas));
  }
  for (i = spp; i < (nonodes); i++) {
    q = NULL;
    for (j = 1; j <= 3; j++) {
      p = (node *)Malloc(sizeof(node));
      p->numsteps = (steptr)Malloc(chars*sizeof(long));
      p->siteset = (seqptr)Malloc(chars*sizeof(sitearray));
      p->seq = (aas *)Malloc(chars*sizeof(aas));
      p->next = q;
      q = p;
    }
    p->next->next->next = p;
    treenode[i] = p;
  }
}  /* protalloctree */


void reallocnode(node* p) 
{
  free(p->numsteps);
  free(p->siteset);
  free(p->seq);
  p->numsteps = (steptr)Malloc(chars*sizeof(long));
  p->siteset = (seqptr)Malloc(chars*sizeof(sitearray));
  p->seq = (aas *)Malloc(chars*sizeof(aas));
}


void reallocchars(void) 
{ /* reallocates variables that are dependand on the number of chars
   * do we need to reallocate the garbage list too? */
  long i;
  node *p;

  if (usertree)
    for (i = 0; i < maxuser; i++) {
      free(fsteps[i]);
      fsteps[i] = (long *)Malloc(chars*sizeof(long));
    }
  
  for (i = 0; i < nonodes; i++) {
    reallocnode(treenode[i]);  
    if (i >= spp) {
      p=treenode[i]->next;
      while (p != treenode[i])  {
        reallocnode(p); 
        p = p->next;
      }
    }
  }

  free(weight);
  free(threshwt);
  free(temp->numsteps);
  free(temp->siteset);
  free(temp->seq);
  free(temp1->numsteps);
  free(temp1->siteset); 
  free(temp1->seq);

  weight = (steptr)Malloc(chars*sizeof(long));
  threshwt = (steptr)Malloc(chars*sizeof(long));
  temp->numsteps = (steptr)Malloc(chars*sizeof(long));
  temp->siteset = (seqptr)Malloc(chars*sizeof(sitearray));
  temp->seq = (aas *)Malloc(chars*sizeof(aas));
  temp1->numsteps = (steptr)Malloc(chars*sizeof(long));
  temp1->siteset = (seqptr)Malloc(chars*sizeof(sitearray));
  temp1->seq = (aas *)Malloc(chars*sizeof(aas));
}


void allocrest()
{ /* allocate remaining global arrays and variables dynamically */
  long i;

  if (usertree) {
    fsteps = (long **)Malloc(maxuser*sizeof(long *));
    for (i = 0; i < maxuser; i++)
      fsteps[i] = (long *)Malloc(chars*sizeof(long));
  }
  bestrees = (bestelm *)Malloc(maxtrees*sizeof(bestelm));
  for (i = 1; i <= maxtrees; i++)
    bestrees[i - 1].btree = (long *)Malloc(spp*sizeof(long));
  nayme = (naym *)Malloc(spp*sizeof(naym));
  enterorder = (long *)Malloc(spp*sizeof(long));
  place = (long *)Malloc(nonodes*sizeof(long));
  weight = (steptr)Malloc(chars*sizeof(long));
  threshwt = (steptr)Malloc(chars*sizeof(long));
  temp = (node *)Malloc(sizeof(node));
  temp->numsteps = (steptr)Malloc(chars*sizeof(long));
  temp->siteset = (seqptr)Malloc(chars*sizeof(sitearray));
  temp->seq = (aas *)Malloc(chars*sizeof(aas));
  temp1 = (node *)Malloc(sizeof(node));
  temp1->numsteps = (steptr)Malloc(chars*sizeof(long));
  temp1->siteset = (seqptr)Malloc(chars*sizeof(sitearray));
  temp1->seq = (aas *)Malloc(chars*sizeof(aas));
}  /* allocrest */


void doinit()
{
  /* initializes variables */

  inputnumbers(&spp, &chars, &nonodes, 1);
  getoptions();
  if (printdata)
    fprintf(outfile, "%2ld species, %3ld  sites\n\n", spp, chars);
  protalloctree();
  allocrest();
}  /* doinit*/


void protinputdata()
{
  /* input the names and sequences for each species */
  long i, j, k, l, aasread, aasnew = 0;
  Char charstate;
  boolean allread, done;
  aas aa;   /* temporary amino acid for input */

  if (printdata)
    headings(chars, "Sequences", "---------");
  aasread = 0;
  allread = false;
  while (!(allread)) {
    /* eat white space -- if the separator line has spaces on it*/
    do {
      charstate = gettc(infile);
    } while (charstate == ' ' || charstate == '\t');
    ungetc(charstate, infile);
    if (eoln(infile)) {
      scan_eoln(infile);
    }
    i = 1;
    while (i <= spp) {
      if ((interleaved && aasread == 0) || !interleaved)
        initname(i - 1);
      j = interleaved ? aasread : 0;
      done = false;
      while (!done && !eoff(infile)) {
        if (interleaved)
          done = true;
        while (j < chars && !(eoln(infile) || eoff(infile))) {
          charstate = gettc(infile);
          if (charstate == '\n' || charstate == '\t')
            charstate = ' ';
          if (charstate == ' ' || (charstate >= '0' && charstate <= '9'))
            continue;
          uppercase(&charstate);
          if ((!isalpha(charstate) && charstate != '?' &&
               charstate != '-' && charstate != '*') || charstate == 'J' ||
              charstate == 'O' || charstate == 'U') {
            printf("WARNING -- BAD AMINO ACID:%c",charstate);
            printf(" AT POSITION%5ld OF SPECIES %3ld\n",j,i);
            exxit(-1);
          }
          j++;
          aa = (charstate == 'A') ?  ala :
               (charstate == 'B') ?  asx :
               (charstate == 'C') ?  cys :
               (charstate == 'D') ?  asp :
               (charstate == 'E') ?  glu :
               (charstate == 'F') ?  phe :
               (charstate == 'G') ?  gly : aa;
          aa = (charstate == 'H') ?  his :
               (charstate == 'I') ? ileu :
               (charstate == 'K') ?  lys :
               (charstate == 'L') ?  leu :
               (charstate == 'M') ?  met :
               (charstate == 'N') ?  asn :
               (charstate == 'P') ?  pro :
               (charstate == 'Q') ?  gln :
               (charstate == 'R') ?  arg : aa;
          aa = (charstate == 'S') ?  ser :
               (charstate == 'T') ?  thr :
               (charstate == 'V') ?  val :
               (charstate == 'W') ?  trp :
               (charstate == 'X') ?  unk :
               (charstate == 'Y') ?  tyr :
               (charstate == 'Z') ?  glx :
               (charstate == '*') ? stop :
               (charstate == '?') ? quest:
               (charstate == '-') ? del  :  aa;

          treenode[i - 1]->seq[j - 1] = aa;
          memcpy(treenode[i - 1]->siteset[j - 1],
                 translate[(long)aa - (long)ala], sizeof(sitearray));
        }
        if (interleaved)
          continue;
        if (j < chars) 
          scan_eoln(infile);
        else if (j == chars)
          done = true;
      }
      if (interleaved && i == 1)
        aasnew = j;
      scan_eoln(infile);
      if ((interleaved && j != aasnew) || ((!interleaved) && j != chars)){
        printf("ERROR: SEQUENCES OUT OF ALIGNMENT\n");
        exxit(-1);}
      i++;
    }
    if (interleaved) {
      aasread = aasnew;
      allread = (aasread == chars);
    } else
      allread = (i > spp);
  }
  if (printdata) {
    for (i = 1; i <= ((chars - 1) / 60 + 1); i++) {
      for (j = 1; j <= (spp); j++) {
        for (k = 0; k < nmlngth; k++)
          putc(nayme[j - 1][k], outfile);
        fprintf(outfile, "   ");
        l = i * 60;
        if (l > chars)
          l = chars;
        for (k = (i - 1) * 60 + 1; k <= l; k++) {
          if (j > 1 && treenode[j - 1]->seq[k - 1] == treenode[0]->seq[k - 1])
            charstate = '.';
          else {
            tmpa = treenode[j-1]->seq[k-1];
              charstate =  (tmpa == ala) ? 'A' :
                           (tmpa == asx) ? 'B' :
                           (tmpa == cys) ? 'C' :
                           (tmpa == asp) ? 'D' :
                           (tmpa == glu) ? 'E' :
                           (tmpa == phe) ? 'F' :
                           (tmpa == gly) ? 'G' :
                           (tmpa == his) ? 'H' :
                           (tmpa ==ileu) ? 'I' :
                           (tmpa == lys) ? 'K' :
                           (tmpa == leu) ? 'L' : charstate;
              charstate =  (tmpa == met) ? 'M' :
                           (tmpa == asn) ? 'N' :
                           (tmpa == pro) ? 'P' :
                           (tmpa == gln) ? 'Q' :
                           (tmpa == arg) ? 'R' :
                           (tmpa == ser) ? 'S' :
                           (tmpa ==ser1) ? 'S' :
                           (tmpa ==ser2) ? 'S' : charstate;
              charstate =  (tmpa == thr) ? 'T' :
                           (tmpa == val) ? 'V' :
                           (tmpa == trp) ? 'W' :
                           (tmpa == unk) ? 'X' :
                           (tmpa == tyr) ? 'Y' :
                           (tmpa == glx) ? 'Z' :
                           (tmpa == del) ? '-' :
                           (tmpa ==stop) ? '*' :
                           (tmpa==quest) ? '?' : charstate;
        }
          putc(charstate, outfile);
          if (k % 10 == 0 && k % 60 != 0)
            putc(' ', outfile);
        }
        putc('\n', outfile);
      }
      putc('\n', outfile);
    }
    putc('\n', outfile);
  }
  putc('\n', outfile);
}  /* protinputdata */


void protmakevalues()
{
  /* set up fractional likelihoods at tips */
  long i, j;
  node *p;

  for (i = 1; i <= nonodes; i++) {
    treenode[i - 1]->back = NULL;
    treenode[i - 1]->tip = (i <= spp);
    treenode[i - 1]->index = i;
    for (j = 0; j < (chars); j++)
      treenode[i - 1]->numsteps[j] = 0;
    if (i > spp) {
      p = treenode[i - 1]->next;
      while (p != treenode[i - 1]) {
        p->back = NULL;
        p->tip = false;
        p->index = i;
        for (j = 0; j < (chars); j++)
          p->numsteps[j] = 0;
        p = p->next;
      }
    }
  }
}  /* protmakevalues */


void doinput()
{
  /* reads the input data */
  long i;

  if (justwts) {
    if (firstset)
      protinputdata();
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
    for (i = 0; i < chars; i++)
      weight[i] = 1;
    if (weights) {
      inputweights(chars, weight, &weights);
    }
    if (weights)
      printweights(outfile, 0, chars, weight, "Sites");
    protinputdata();
  }
  if(!thresh)
    threshold = spp * 3.0;
  for(i = 0 ; i < (chars) ; i++){
    weight[i]*=10;
    threshwt[i] = (long)(threshold * weight[i] + 0.5);      
  }

  protmakevalues();
}  /* doinput */


void protfillin(node *p, node *left, node *rt)
{
  /* sets up for each node in the tree the aa set for site m
     at that point and counts the changes.  The program
     spends much of its time in this function */
  boolean counted, done;
  aas aa;
  long s = 0;
  sitearray ls, rs, qs;
  long i, j, m, n;

  for (m = 0; m < chars; m++) {
    if (left != NULL)
      memcpy(ls, left->siteset[m], sizeof(sitearray));
    if (rt != NULL)
      memcpy(rs, rt->siteset[m], sizeof(sitearray));
    if (left == NULL) {
      n = rt->numsteps[m];
      memcpy(qs, rs, sizeof(sitearray));
    }
    else if (rt == NULL) {
      n = left->numsteps[m];
      memcpy(qs, ls, sizeof(sitearray));
    }
    else {
      n = left->numsteps[m] + rt->numsteps[m];
      if ((ls[0] == rs[0]) && (ls[1] == rs[1]) && (ls[2] == rs[2])) {
        qs[0] = ls[0];
        qs[1] = ls[1];
        qs[2] = ls[2];
      }
      else {
        counted = false;
        for (i = 0; (!counted) && (i <= 3); i++) {
          switch (i) {
 
            case 0:
              s = ls[0] & rs[0];
              break;

            case 1:
              s = (ls[0] & rs[1]) | (ls[1] & rs[0]);
              break;

            case 2:
              s = (ls[0] & rs[2]) | (ls[1] & rs[1]) | (ls[2] & rs[0]);
              break;

            case 3:
              s = ls[0] | (ls[1] & rs[2]) | (ls[2] & rs[1]) | rs[0];
              break;

          }
          if (s != 0) {
            qs[0] = s;
            counted = true;
          } else
              n += weight[m];
        }
        switch (i) {
          case 1:
            qs[1] = qs[0] | (ls[0] & rs[1]) | (ls[1] & rs[0]);
            qs[2] = qs[1] | (ls[0] & rs[2]) | (ls[1] & rs[1]) | (ls[2] & rs[0]);
            break;
          case 2:
            qs[1] = qs[0] | (ls[0] & rs[2]) | (ls[1] & rs[1]) | (ls[2] & rs[0]);
            qs[2] = qs[1] | ls[0] | (ls[1] & rs[2]) | (ls[2] & rs[1]) | rs[0];
            break;
          case 3:
            qs[1] = qs[0] | ls[0] | (ls[1] & rs[2]) | (ls[2] & rs[1]) | rs[0];
            qs[2] = qs[1] | ls[1] | (ls[2] & rs[2]) | rs[1];
            break;
          case 4:
            qs[1] = qs[0] | ls[1] | (ls[2] & rs[2]) | rs[1];
            qs[2] = qs[1] | ls[2] | rs[2];
            break;
        }
        for (aa = ala; (long)aa <= (long)stop; aa = (aas)((long)aa + 1)) {
          done = false;
          for (i = 0; (!done) && (i <= 1); i++) {
            if (((1L << ((long)aa)) & qs[i]) != 0) {
              for (j = i+1; j <= 2; j++)
                qs[j] |= translate[(long)aa - (long)ala][j-i];
              done = true;
            }
          }
        }
      }
    }
    p->numsteps[m] = n;
    memcpy(p->siteset[m], qs, sizeof(sitearray));
  }
}  /* protfillin */


void protpreorder(node *p)
{
  /* recompute number of steps in preorder taking both ancestoral and
     descendent steps into account */
  if (p != NULL && !p->tip) {
    protfillin (p->next, p->next->next->back, p->back);
    protfillin (p->next->next, p->back, p->next->back);
    protpreorder (p->next->back);
    protpreorder (p->next->next->back);
  }
} /* protpreorder */


void protadd(node *below, node *newtip, node *newfork)
{
  /* inserts the nodes newfork and its left descendant, newtip,
     to the tree.  below becomes newfork's right descendant */

  if (below != treenode[below->index - 1])
    below = treenode[below->index - 1];
  if (below->back != NULL)
    below->back->back = newfork;
  newfork->back = below->back;
  below->back = newfork->next->next;
  newfork->next->next->back = below;
  newfork->next->back = newtip;
  newtip->back = newfork->next;
  if (root == below)
    root = newfork;
  root->back = NULL;

  if (recompute) {
    protfillin (newfork, newfork->next->back, newfork->next->next->back);
    protpreorder(newfork);
    if (newfork != root)
      protpreorder(newfork->back);
  }
}  /* protadd */


void protre_move(node **item, node **fork)
{
  /* removes nodes item and its ancestor, fork, from the tree.
     the new descendant of fork's ancestor is made to be
     fork's second descendant (other than item).  Also
     returns pointers to the deleted nodes, item and fork */
  node *p, *q, *other;

  if ((*item)->back == NULL) {
    *fork = NULL;
    return;
  }
  *fork = treenode[(*item)->back->index - 1];
  if ((*item) == (*fork)->next->back)
    other = (*fork)->next->next->back;
  else other = (*fork)->next->back;
  if (root == *fork)
    root = other;
  p = (*item)->back->next->back;
  q = (*item)->back->next->next->back;
  if (p != NULL) p->back = q;
  if (q != NULL) q->back = p;
  (*fork)->back = NULL;
  p = (*fork)->next;
  do {
    p->back = NULL;
    p = p->next;
  } while (p != (*fork));
  (*item)->back = NULL;
  if (recompute) {
    protpreorder(other);
    if (other != root) protpreorder(other->back);
  }
}  /* protre_move */


void evaluate(node *r)
{
  /* determines the number of steps needed for a tree. this is the
     minimum number of steps needed to evolve sequences on this tree */
  long i, steps, term;
  double sum;

  sum = 0.0;
  for (i = 0; i < (chars); i++) {
    steps = r->numsteps[i];
    if (steps <= threshwt[i])
      term = steps;
    else
      term = threshwt[i];
    sum += term;
    if (usertree && which <= maxuser)
      fsteps[which - 1][i] = term;
  }
  if (usertree && which <= maxuser) {
    nsteps[which - 1] = sum;
    if (which == 1) {
      minwhich = 1;
      minsteps = sum;
    } else if (sum < minsteps) {
      minwhich = which;
      minsteps = sum;
    }
  }
  like = -sum;
}  /* evaluate */


void protpostorder(node *p)
{
  /* traverses a binary tree, calling PROCEDURE fillin at a
     node's descendants before calling fillin at the node */
  if (p->tip)
    return;
  protpostorder(p->next->back);
  protpostorder(p->next->next->back);
  protfillin(p, p->next->back, p->next->next->back);
}  /* protpostorder */


void protreroot(node *outgroup)
{
  /* reorients tree, putting outgroup in desired position. */
  node *p, *q;

  if (outgroup->back->index == root->index)
    return;
  p = root->next;
  q = root->next->next;
  p->back->back = q->back;
  q->back->back = p->back;
  p->back = outgroup;
  q->back = outgroup->back;
  outgroup->back->back = q;
  outgroup->back = p;
}  /* protreroot */


void protsavetraverse(node *p, long *pos, boolean *found)
{
  /* sets BOOLEANs that indicate which way is down */
  p->bottom = true;
  if (p->tip)
    return;
  p->next->bottom = false;
  protsavetraverse(p->next->back, pos,found);
  p->next->next->bottom = false;
  protsavetraverse(p->next->next->back, pos,found);
}  /* protsavetraverse */


void protsavetree(long *pos, boolean *found)
{
  /* record in place where each species has to be
     added to reconstruct this tree */
  long i, j;
  node *p;
  boolean done;

  protreroot(treenode[outgrno - 1]);
  protsavetraverse(root, pos,found);
  for (i = 0; i < (nonodes); i++)
    place[i] = 0;
  place[root->index - 1] = 1;
  for (i = 1; i <= (spp); i++) {
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
}  /* protsavetree */


void tryadd(node *p, node **item, node **nufork)
{
  /* temporarily adds one fork and one tip to the tree.
     if the location where they are added yields greater
     "likelihood" than other locations tested up to that
     time, then keeps that location as there */
  long pos;
  boolean found;
  node *rute, *q;

  if (p == root)
    protfillin(temp, *item, p);
  else {
    protfillin(temp1, *item, p);
    protfillin(temp, temp1, p->back);
  }
  evaluate(temp);
  if (lastrearr) {
    if (like < bestlike) {
      if ((*item) == (*nufork)->next->next->back) {
        q = (*nufork)->next;
        (*nufork)->next = (*nufork)->next->next;
        (*nufork)->next->next = q;
        q->next = (*nufork);
      }
    }
    else if (like >= bstlike2) {
      recompute = false;
      protadd(p, (*item), (*nufork));
      rute = root->next->back;
      protsavetree(&pos,&found);
      protreroot(rute);
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
      protre_move (item, nufork);
      recompute = true;
    }
  }
  if (like >= bestyet) {
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
  tryadd(p, &item,&nufork);
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
  protfillin(temp1, frombelow, q);
  protfillin(temp, temp1, p);
  protfillin(temp1, temp, whereto->back);
  evaluate(temp1);
  if (like - oldlike < LIKE_EPSILON) {
    if (p == forknode->next->next->back) {
      q = forknode->next;
      forknode->next = forknode->next->next;
      forknode->next->next = q;
      q->next = forknode;
    }
  }
  else {
    recompute = false;
    protre_move(&p, &forknode);
    protfillin(whereto, whereto->next->back, whereto->next->next->back);
    recompute = true;
    protadd(whereto, p, forknode);
    *success = true;
    bestyet = like;
  }
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
  boolean success = true;
  while (success) {
    success = false;
    repreorder(*r, &success);
  }
}  /* rearrange */


void protgetch(Char *c)
{
  /* get next nonblank character */
  do {
    if (eoln(intree))
      scan_eoln(intree);
    *c = gettc(intree);
    if (*c == '\n' || *c == '\t')
      *c = ' ';
  } while (!(*c != ' ' || eoff(intree)));
}  /* protgetch */


void protaddelement(node **p,long *nextnode,long *lparens,boolean *names)
{
  /* recursive procedure adds nodes to user-defined tree */
  node *q;
  long i, n;
  boolean found;
  Char str[nmlngth];

  protgetch(&ch);
  
  if (ch == '(' ) {
    if ((*lparens) >= spp - 1) {
      printf("\nERROR IN USER TREE: TOO MANY LEFT PARENTHESES\n");
      exxit(-1);
    }
    (*nextnode)++;
    (*lparens)++;
    q = treenode[(*nextnode) - 1];
    protaddelement(&q->next->back, nextnode,lparens,names);
    q->next->back->back = q->next;
    findch(',', &ch, which);
    protaddelement(&q->next->next->back, nextnode,lparens,names);
    q->next->next->back->back = q->next->next;
    findch(')', &ch, which);
    *p = q;
    return;
  }
  for (i = 0; i < nmlngth; i++)
    str[i] = ' ';
  n = 1;
  do {
    if (ch == '_')
      ch = ' ';
    str[n - 1] = ch;
    if (eoln(intree))
      scan_eoln(intree);
    ch = gettc(intree);
    n++;
  } while (ch != ',' && ch != ')' && ch != ':' && n <= nmlngth);
  n = 1;
  do {
    found = true;
    for (i = 0; i < nmlngth; i++)
      found = (found && ((str[i] == nayme[n - 1][i]) ||
                         ((nayme[n - 1][i] == '_') && (str[i] == ' '))));
    if (found) {
      if (names[n - 1] == false) {
        *p = treenode[n - 1];
        names[n - 1] = true;
      } else {
        printf("\nERROR IN USER TREE: DUPLICATE NAME FOUND -- ");
        for (i = 0; i < nmlngth; i++)
          putchar(nayme[n - 1][i]);
        putchar('\n');
        exxit(-1);
      }
    } else
      n++;
  } while (!(n > spp || found));
  if (n <= spp)
    return;
  printf("CANNOT FIND SPECIES: ");
  for (i = 0; i < nmlngth; i++)
    putchar(str[i]);
  putchar('\n');
}  /* protaddelement */


void prottreeread()
{
  /* read in user-defined tree and set it up */
  long nextnode, lparens, i;

  root = treenode[spp];
  nextnode = spp;
  root->back = NULL;
  names = (boolean *)Malloc(spp*sizeof(boolean));
  for (i = 0; i < (spp); i++)
    names[i] = false;
  lparens = 0;
  protaddelement(&root, &nextnode,&lparens,names);
  if (ch == '[') {
    do
      ch = gettc(intree);
    while (ch != ']');
    ch = gettc(intree);
  }
  findch(';', &ch, which);
  scan_eoln(intree);
  free(names);
}  /* prottreeread */


void protancestset(long *a, long *b, long *c, long *d, long *k)
{
  /* sets up the aa set array. */
  aas aa;
  long s, sa, sb;
  long i, j, m, n;
  boolean counted;

  counted = false;
  *k = 0;
  for (i = 0; i <= 5; i++) {
    if (*k < 3) {
      s = 0;
      if (i > 3)
        n = i - 3;
      else
        n = 0;
      for (j = n; j <= (i - n); j++) {
        if (j < 3)
          sa = a[j];
        else
          sa = fullset;
        for (m = n; m <= (i - j - n); m++) {
          if (m < 3)
            sb = sa & b[m];
          else
            sb = sa;
          if (i - j - m < 3)
            sb &= c[i - j - m];
          s |= sb;
        }
      }
      if (counted || s != 0) {
        d[*k] = s;
        (*k)++;
        counted = true;
      }
    }
  }
  for (i = 0; i <= 1; i++) {
    for (aa = ala; (long)aa <= (long)stop; aa = (aas)((long)aa + 1)) {
      if (((1L << ((long)aa)) & d[i]) != 0) {
        for (j = i + 1; j <= 2; j++)
          d[j] |= translate[(long)aa - (long)ala][j - i];
      }
    }
  }
}  /* protancestset */


void prothyprint(long b1, long b2, boolean *bottom, node *r,
                        boolean *nonzero, boolean *maybe)
{
  /* print out states in sites b1 through b2 at node */
  long i;
  boolean dot;
  Char ch = 0;
  aas aa;

  if (*bottom) {
    if (!outgropt)
      fprintf(outfile, "      ");
    else
      fprintf(outfile, "root  ");
  } else
    fprintf(outfile, "%3ld   ", r->back->index - spp);
  if (r->tip) {
    for (i = 0; i < nmlngth; i++)
      putc(nayme[r->index - 1][i], outfile);
  } else
    fprintf(outfile, "%4ld      ", r->index - spp);
  if (*bottom)
    fprintf(outfile, "          ");
  else if (*nonzero)
    fprintf(outfile, "   yes    ");
  else if (*maybe)
    fprintf(outfile, "  maybe   ");
  else
    fprintf(outfile, "   no     ");
  for (i = b1 - 1; i < b2; i++) {
    aa = r->seq[i];
    switch (aa) {

    case ala:
      ch = 'A';
      break;

    case asx:
      ch = 'B';
      break;

    case cys:
      ch = 'C';
      break;

    case asp:
      ch = 'D';
      break;

    case glu:
      ch = 'E';
      break;

    case phe:
      ch = 'F';
      break;

    case gly:
      ch = 'G';
      break;

    case his:
      ch = 'H';
      break;

    case ileu:
      ch = 'I';
      break;

    case lys:
      ch = 'K';
      break;

    case leu:
      ch = 'L';
      break;

    case met:
      ch = 'M';
      break;

    case asn:
      ch = 'N';
      break;

    case pro:
      ch = 'P';
      break;

    case gln:
      ch = 'Q';
      break;

    case arg:
      ch = 'R';
      break;

    case ser:
      ch = 'S';
      break;

    case ser1:
      ch = 'S';
      break;

    case ser2:
      ch = 'S';
      break;

    case thr:
      ch = 'T';
      break;

    case trp:
      ch = 'W';
      break;

    case tyr:
      ch = 'Y';
      break;

    case val:
      ch = 'V';
      break;

    case glx:
      ch = 'Z';
      break;

    case del:
      ch = '-';
      break;

    case stop:
      ch = '*';
      break;

    case unk:
      ch = 'X';
      break;

    case quest:
      ch = '?';
      break;
    }
    if (!(*bottom) && dotdiff)
      dot = (r->siteset[i] [0] == treenode[r->back->index - 1]->siteset[i][0]
             || ((r->siteset[i][0] &
                  (~((1L << ((long)ser1)) | (1L << ((long)ser2)) |
                                         (1L << ((long)ser))))) == 0 &&
                (treenode[r->back->index - 1]->siteset[i] [0] &
                  (~((1L << ((long)ser1)) | (1L << ((long)ser2)) |
                                         (1L << ((long)ser))))) == 0));
    else
      dot = false;
    if (dot)
      putc('.', outfile);
    else
      putc(ch, outfile);
    if ((i + 1) % 10 == 0)
      putc(' ', outfile);
  }
  putc('\n', outfile);
}  /* prothyprint */


void prothyptrav(node *r, sitearray *hypset, long b1, long b2, long *k,
                        boolean *bottom, sitearray nothing)
{
  boolean maybe, nonzero;
  long i;
  aas aa;
  long anc = 0, hset;
  gseq *ancset, *temparray;

  protgnu(&ancset);
  protgnu(&temparray);
  maybe = false;
  nonzero = false;
  for (i = b1 - 1; i < b2; i++) {
    if (!r->tip) {
      protancestset(hypset[i], r->next->back->siteset[i],
                r->next->next->back->siteset[i], temparray->seq[i], k);
      memcpy(r->siteset[i], temparray->seq[i], sizeof(sitearray));
    }
    if (!(*bottom))
      anc = treenode[r->back->index - 1]->siteset[i][0];
    if (!r->tip) {
      hset = r->siteset[i][0];
      r->seq[i] = quest;
      for (aa = ala; (long)aa <= (long)stop; aa = (aas)((long)aa + 1)) {
        if (hset == 1L << ((long)aa))
          r->seq[i] = aa;
      }
      if (hset == ((1L << ((long)asn)) | (1L << ((long)asp))))
        r->seq[i] = asx;
      if (hset == ((1L << ((long)gln)) | (1L << ((long)gly))))
        r->seq[i] = glx;
      if (hset == ((1L << ((long)ser1)) | (1L << ((long)ser2))))
        r->seq[i] = ser;
      if (hset == fullset)
        r->seq[i] = unk;
    }
    nonzero = (nonzero || (r->siteset[i][0] & anc) == 0);
    maybe = (maybe || r->siteset[i][0] != anc);
  }
  prothyprint(b1, b2,bottom,r,&nonzero,&maybe);
  *bottom = false;
  if (!r->tip) {
    memcpy(temparray->seq, r->next->back->siteset, chars*sizeof(sitearray));
    for (i = b1 - 1; i < b2; i++)
      protancestset(hypset[i], r->next->next->back->siteset[i], nothing,
                ancset->seq[i], k);
    prothyptrav(r->next->back, ancset->seq, b1, b2,k,bottom,nothing );
    for (i = b1 - 1; i < b2; i++)
      protancestset(hypset[i], temparray->seq[i], nothing, ancset->seq[i],k);
    prothyptrav(r->next->next->back, ancset->seq, b1, b2, k,bottom,nothing);
  }
  protchuck(temparray);
  protchuck(ancset);
}  /* prothyptrav */


void prothypstates(long *k)
{
  /* fill in and describe states at interior nodes */
  boolean bottom;
  sitearray nothing;
  long i, n;
  seqptr hypset;

  fprintf(outfile, "\nFrom    To     Any Steps?    State at upper node\n");
  fprintf(outfile, "                             ");
  fprintf(outfile, "( . means same as in the node below it on tree)\n\n");
  memcpy(nothing, translate[(long)quest - (long)ala], sizeof(sitearray));
  hypset = (seqptr)Malloc(chars*sizeof(sitearray));
  for (i = 0; i < (chars); i++)
    memcpy(hypset[i], nothing, sizeof(sitearray));
  bottom = true;
  for (i = 1; i <= ((chars - 1) / 40 + 1); i++) {
    putc('\n', outfile);
    n = i * 40;
    if (n > chars)
      n = chars;
    bottom = true;
    prothyptrav(root, hypset, i * 40 - 39, n, k,&bottom,nothing);
  }
  free(hypset);
}  /* prothypstates */


void describe()
{
  /* prints ancestors, steps and table of numbers of steps in
     each site */
  long i,j,k;

  if (treeprint)
    fprintf(outfile, "\nrequires a total of %10.3f\n", like / -10);
  if (stepbox) {
    putc('\n', outfile);
    if (weights)
      fprintf(outfile, "weighted ");
    fprintf(outfile, "steps in each position:\n");
    fprintf(outfile, "      ");
    for (i = 0; i <= 9; i++)
      fprintf(outfile, "%4ld", i);
    fprintf(outfile, "\n     *-----------------------------------------\n");
    for (i = 0; i <= (chars / 10); i++) {
      fprintf(outfile, "%5ld", i * 10);
      putc('!', outfile);
      for (j = 0; j <= 9; j++) {
        k = i * 10 + j;
        if (k == 0 || k > chars)
          fprintf(outfile, "    ");
        else
          fprintf(outfile, "%4ld", root->numsteps[k - 1] / 10);
      }
      putc('\n', outfile);
    }
  }
  if (ancseq) {
    prothypstates(&k);
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
  /* constructs a binary tree from the pointers in treenode.
     adds each node at location which yields highest "likelihood"
     then rearranges the tree for greatest "likelihood" */
  long i, j, numtrees;
  double gotlike;
  node *item, *nufork, *dummy;

  if (!usertree) {
    for (i = 1; i <= (spp); i++)
      enterorder[i - 1] = i;
    if (jumble)
      randumize(seed, enterorder);
    root = treenode[enterorder[0] - 1];
    recompute = true;
    protadd(treenode[enterorder[0] - 1], treenode[enterorder[1] - 1],
        treenode[spp]);
    if (progress) {
      printf("\nAdding species:\n");
      writename(0, 2, enterorder);
    }
    lastrearr = false;
    for (i = 3; i <= (spp); i++) {
      bestyet = -30.0*spp*chars;
      there = root;
      item = treenode[enterorder[i - 1] - 1];
      nufork = treenode[spp + i - 2];
      addpreorder(root, item, nufork);
      protadd(there, item, nufork);
      like = bestyet;
      rearrange(&root);
      if (progress)
        writename(i - 1, 1, enterorder);
      lastrearr = (i == spp);
      if (lastrearr) {
        if (progress) {
          printf("\nDoing global rearrangements\n");
          printf("  !");
          for (j = 1; j <= nonodes; j++)
            if ( j % (( nonodes / 72 ) + 1 ) == 0 )
              putchar('-');
          printf("!\n");
        }
        bestlike = bestyet;
        if (jumb == 1) {
          bstlike2 = bestlike = -30.0*spp*chars;
          nextree = 1;
        }
        do {
          if (progress)
            printf("   ");
          gotlike = bestlike;
          for (j = 0; j < (nonodes); j++) {
            bestyet = -30.0*spp*chars;
            item = treenode[j];
            if (item != root) {
              nufork = treenode[treenode[j]->back->index - 1];
              protre_move(&item, &nufork);
              there = root;
              addpreorder(root, item, nufork);
              protadd(there, item, nufork);
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
      protre_move(&treenode[i], &dummy);
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
        protadd(treenode[0], treenode[1], treenode[spp]);
        for (j = 3; j <= (spp); j++)
          protadd(treenode[bestrees[i].btree[j - 1] - 1], treenode[j - 1],
              treenode[spp + j - 2]);
        protreroot(treenode[outgrno - 1]);
        protpostorder(root);
        evaluate(root);
        printree(root, 1.0);
        describe();
        for (j = 1; j < (spp); j++)
          protre_move(&treenode[j], &dummy);
      }
    }
  } else {
    /* Open in binary: ftell() is broken for UNIX line-endings under WIN32 */
    openfile(&intree,INTREE,"input tree file", "rb",progname,intreename);
    numtrees = countsemic(&intree);
    if (treeprint) {
      fprintf(outfile, "User-defined tree");
      if (numtrees > 1)
        putc('s', outfile);
      fprintf(outfile, ":\n\n\n\n");
    }
    which = 1;
    while (which <= numtrees) {
      prottreeread();
      if (outgropt)
        protreroot(treenode[outgrno - 1]);
      protpostorder(root);
      evaluate(root);
      printree(root, 1.0);
      describe();
      which++;
    }
    printf("\n");
    FClose(intree);
    putc('\n', outfile);
    if (numtrees > 1 && chars > 1 )
      standev(chars, numtrees, minwhich, minsteps, nsteps, fsteps, seed);
  }
  if (jumb == njumble && progress) {
    printf("Output written to file \"%s\"\n\n", outfilename);
    if (trout)
      printf("Trees also written onto file \"%s\"\n\n", outtreename);
  }
}  /* maketree */


int main(int argc, Char *argv[])
{  /* Protein parsimony by uphill search */
#ifdef MAC
  argc = 1;         /* macsetup("Protpars","");                */
  argv[0] = "Protpars";
#endif
  init(argc,argv);
  progname = argv[0];
  openfile(&infile,INFILE,"input file", "r",argv[0],infilename);
  openfile(&outfile,OUTFILE,"output file", "w",argv[0],outfilename);

  ibmpc = IBMCRT;
  ansi = ANSICRT;
  garbage = NULL;
  mulsets = false;
  msets = 1;
  firstset = true;
  code();
  setup();
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
      fprintf(outfile, "Data set # %ld:\n\n",ith);
      if (progress)
        printf("Data set # %ld:\n\n",ith);
    }
    for (jumb = 1; jumb <= njumble; jumb++)
      maketree();
  }
  FClose(infile);
  FClose(outfile);
  FClose(outtree);
#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(outtreename);
#endif
  printf("Done.\n\n");

#ifdef WIN32
  phyRestoreConsoleAttributes();
#endif
  return 0;
}  /* Protein parsimony by uphill search */
