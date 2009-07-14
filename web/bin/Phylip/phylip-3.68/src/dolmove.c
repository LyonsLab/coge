#include "phylip.h"
#include "moves.h"
#include "disc.h"
#include "dollo.h"

/* version 3.6. (c) Copyright 1993-2004 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define overr           4
#define which           1

typedef enum {
  horiz, vert, up, overt, upcorner, downcorner, onne, zerro, question, polym
} chartype;

typedef enum {  rearr, flipp, reroott, none } rearrtype;

typedef enum {
  arb, use, spec
} howtree;

#ifndef OLDC
/* function prototypes */
void   getoptions(void);
void   inputoptions(void);
void   allocrest(void);
void   doinput(void);
void   configure(void);
void   prefix(chartype);
void   postfix(chartype);
void   makechar(chartype);
void   dolmove_correct(node *);
void   dolmove_count(node *);

void   preorder(node *);
void   evaluate(node *);
void   reroot(node *);
void   dolmove_hyptrav(node *);
void   dolmove_hypstates(void);
void   grwrite(chartype, long, long *);
void   dolmove_drawline(long);
void   dolmove_printree(void);
void   arbitree(void);
void   yourtree(void);
                
void   initdolmovenode(node **, node **, node *, long, long, long *,
        long *, initops, pointarray, pointarray, Char *, Char *, FILE *);
void   buildtree(void);
void   rearrange(void);
void   tryadd(node *, node **, node **, double *);
void   addpreorder(node *, node *, node *, double *);
void   try(void);
void   undo(void);
void   treewrite(boolean);
void   clade(void);
void   flip(void);

void   changeoutgroup(void);
void   redisplay(void);
void   treeconstruct(void);
/* function prototypes */
#endif

Char infilename[FNMLNGTH],intreename[FNMLNGTH],outtreename[FNMLNGTH], ancfilename[FNMLNGTH], factfilename[FNMLNGTH], weightfilename[FNMLNGTH];
node *root;
long outgrno, col, screenlines, screenwidth, scrollinc,treelines,
        leftedge,topedge,vmargin,hscroll,vscroll,farthest;
/* outgrno indicates outgroup */
boolean weights, thresh, ancvar, questions, dollo, factors,
               waswritten;
boolean *ancone, *anczero, *ancone0, *anczero0;
Char *factor;
pointptr treenode;   /* pointers to all nodes in tree */
double threshold;
double *threshwt;
unsigned char cha[10];
boolean reversed[10];
boolean graphic[10];
howtree how;
char *progname;
char ch;

/* Variables for treeread */
boolean usertree, goteof, firsttree, haslengths;
pointarray nodep;
node *grbg;
long *zeros;

/* Local variables for treeconstruct, propagated globally for c version: */
long dispchar, dispword, dispbit, atwhat, what, fromwhere, towhere,
  oldoutgrno, compatible;
double like, bestyet, gotlike;
Char *guess;
boolean display, newtree, changed, subtree, written, oldwritten, restoring,
  wasleft, oldleft, earlytree;
boolean *in_tree;
steptr numsteps;
long fullset;
bitptr zeroanc, oneanc;
node *nuroot;
rearrtype lastop;
steptr numsone, numszero;
boolean *names;


void getoptions()
{
  /* interactively set options */
  long loopcount;
  Char ch;
  boolean done, gotopt;
  char input[100];

  how = arb;
  usertree = false;
  goteof = false;
  thresh = false;
  threshold = spp;
  weights = false;
  ancvar = false;
  factors = false;
  dollo = true;
  loopcount = 0;
  do {
    cleerhome();
    printf("\nInteractive Dollo or polymorphism parsimony,");
    printf(" version %s\n\n",VERSION);
    printf("Settings for this run:\n");
    printf("  P                        Parsimony method?");
    printf("  %s\n",(dollo ? "Dollo" : "Polymorphism"));
    printf("  A                    Use ancestral states?  %s\n",
           ancvar ? "Yes" : "No");
    printf("  F                 Use factors information?  %s\n",
           factors ? "Yes" : "No");
    printf("  W                          Sites weighted?  %s\n",
           (weights ? "Yes" : "No"));
    printf("  T                 Use Threshold parsimony?");
    if (thresh)
      printf("  Yes, count steps up to%4.1f\n", threshold);
    else
      printf("  No, use ordinary parsimony\n");
    printf("  A      Use ancestral states in input file?");
    printf("  %s\n",(ancvar ? "Yes" : "No"));
    printf("  U Initial tree (arbitrary, user, specify)?");
    printf("  %s\n",(how == arb) ? "Arbitrary" :
                    (how == use) ? "User tree from tree file" :
                                   "Tree you specify");
    printf("  0      Graphics type (IBM PC, ANSI, none)?  %s\n",
           ibmpc ? "IBM PC" : ansi  ? "ANSI"   : "(none)");
    printf("  L               Number of lines on screen?%4ld\n",screenlines);
    printf("  S                Width of terminal screen?%4ld\n",screenwidth);
    printf(
      "\n\nAre these settings correct? (type Y or the letter for one to change)\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    getstryng(input);
    ch = input[0];
    uppercase(&ch);
    done = (ch == 'Y');
    gotopt = (strchr("SFTPULA0W",ch) != NULL) ? true : false;
    if (gotopt) {
      switch (ch) {

      case 'A':
        ancvar = !ancvar;
        break;

      case 'F':
        factors = !factors;
        break;

      case 'W':
        weights = !weights;
        break;
        
      case 'P':
        dollo = !dollo;
        break;

      case 'T':
        thresh = !thresh;
        if (thresh)
          initthreshold(&threshold);
        break;

      case 'U':
        if (how == arb)
          how = use;
        else if (how == use)
            how = spec;
          else
            how = arb;
        break;

      case '0':
        initterminal(&ibmpc, &ansi);
        break;

      case 'L':
        initnumlines(&screenlines);
        break;

      case 'S':
        screenwidth = readlong("Width of terminal screen (in characters)?\n");
      }
    }
    else  
      printf("Not a possible option!\n");
    countup(&loopcount, 100);
  } while (!done);
  if (scrollinc < screenwidth / 2.0)
    hscroll = scrollinc;
  else
    hscroll = screenwidth / 2;
  if (scrollinc < screenlines / 2.0)
    vscroll = scrollinc;
  else
    vscroll = screenlines / 2;
}  /* getoptions */


void inputoptions()
{
  /* input the information on the options */
  long i;

  scan_eoln(infile);
  for (i = 0; i < (chars); i++)
      weight[i] = 1;
  if (ancvar)
      inputancestors(anczero0, ancone0);
  if (factors)
      inputfactors(chars, factor, &factors);
  if (weights)
      inputweights(chars, weight, &weights);
  putchar('\n');
  if (weights)
    printweights(stdout, 0, chars, weight, "Characters");
  if (factors)
    printfactors(stdout, chars, factor, "");
  for (i = 0; i < (chars); i++) {
    if (!ancvar) {
      anczero[i] = true;
      ancone[i] = false;
    } else {
      anczero[i] = anczero0[i];
      ancone[i] = ancone0[i];
    }
  }
  if (ancvar)
    printancestors(stdout, anczero, ancone);
  if (!thresh)
    threshold = spp;
  questions = false;
  for (i = 0; i < (chars); i++) {
    questions = (questions || (ancone[i] && anczero[i]));
    threshwt[i] = threshold * weight[i];
  }
}  /* inputoptions */


void allocrest()
{
  nayme = (naym *)Malloc(spp*sizeof(naym));
  in_tree = (boolean *)Malloc(nonodes*sizeof(boolean));
  extras = (steptr)Malloc(chars*sizeof(long));
  weight = (steptr)Malloc(chars*sizeof(long));
  numszero = (steptr)Malloc(chars*sizeof(long));
  numsone = (steptr)Malloc(chars*sizeof(long));
  threshwt = (double *)Malloc(chars*sizeof(double));
  factor = (Char *)Malloc(chars*sizeof(Char));
  ancone = (boolean *)Malloc(chars*sizeof(boolean));
  anczero = (boolean *)Malloc(chars*sizeof(boolean));
  ancone0 = (boolean *)Malloc(chars*sizeof(boolean));
  anczero0 = (boolean *)Malloc(chars*sizeof(boolean));
  zeroanc = (bitptr)Malloc(words*sizeof(long));
  oneanc = (bitptr)Malloc(words*sizeof(long));
}  /* allocrest */


void doinput()
{
  /* reads the input data */

  inputnumbers(&spp, &chars, &nonodes, 1);
  words = chars / bits + 1;
  printf("%2ld species, %3ld characters\n", spp, chars);
  printf("\nReading input file ...\n\n");
  getoptions();
  if (weights)
      openfile(&weightfile,WEIGHTFILE,"weights file","r",progname,weightfilename);
  if(ancvar)
      openfile(&ancfile,ANCFILE,"ancestors file", "r",progname,ancfilename);
  if(factors)
      openfile(&factfile,FACTFILE,"factors file", "r",progname,factfilename);

  alloctree(&treenode);
  setuptree(treenode);
  allocrest();
  inputoptions();
  inputdata(treenode, dollo, false, stdout);
}  /* doinput */


void configure()
{
  /* configure to machine -- set up special characters */
  chartype a;

  for (a = horiz; (long)a <= (long)polym; a = (chartype)((long)a + 1))
    reversed[(long)a] = false;
  for (a = horiz; (long)a <= (long)polym; a = (chartype)((long)a + 1))
    graphic[(long)a] = false;
  if (ibmpc) {
    cha[(long)horiz] = 205;
    graphic[(long)horiz] = true;
    cha[(long)vert] = 186;
    graphic[(long)vert] = true;
    cha[(long)up] = 186;
    graphic[(long)up] = true;
    cha[(long)overt] = 205;
    graphic[(long)overt] = true;
    cha[(long)onne] = 219;
    reversed[(long)onne] = true;
    cha[(long)zerro] = 176;
    graphic[(long)zerro] = true;
    cha[(long)question] = 178;   /* or try CHR(177) */
    cha[(long)polym] = '\001';
    reversed[(long)polym] = true;
    cha[(long)upcorner] = 200;
    graphic[(long)upcorner] = true;
    cha[(long)downcorner] = 201;
    graphic[(long)downcorner] = true;
    graphic[(long)question] = true;
    return;
  }
  if (ansi) {
    cha[(long)onne] = ' ';
    reversed[(long)onne] = true;
    cha[(long)horiz] = cha[(long)onne];
    reversed[(long)horiz] = true;
    cha[(long)vert] = cha[(long)onne];
    reversed[(long)vert] = true;
    cha[(long)up] = 'x';
    graphic[(long)up] = true;
    cha[(long)overt] = 'q';
    graphic[(long)overt] = true;
    cha[(long)zerro] = 'a';
    graphic[(long)zerro] = true;
    reversed[(long)zerro] = true;
    cha[(long)question] = '?';
    cha[(long)question] = '?';
    reversed[(long)question] = true;
    cha[(long)polym] = '%';
    reversed[(long)polym] = true;
    cha[(long)upcorner] = 'm';
    graphic[(long)upcorner] = true;
    cha[(long)downcorner] = 'l';
    graphic[(long)downcorner] = true;
    return;
  }
  cha[(long)horiz] = '=';
  cha[(long)vert] = ' ';
  cha[(long)up] = '!';
  cha[(long)overt] = '-';
  cha[(long)onne] = '*';
  cha[(long)zerro] = '=';
  cha[(long)question] = '.';
  cha[(long)polym] = '%';
  cha[(long)upcorner] = '`';
  cha[(long)downcorner] = ',';
}  /* configure */


void prefix(chartype a)
{
  /* give prefix appropriate for this character */
  if (reversed[(long)a])
    prereverse(ansi);
  if (graphic[(long)a])
    pregraph(ansi);
}  /* prefix */


void postfix(chartype a)
{
  /* give postfix appropriate for this character */
  if (reversed[(long)a])
    postreverse(ansi);
  if (graphic[(long)a])
    postgraph(ansi);
}  /* postfix */


void makechar(chartype a)
{
  /* print out a character with appropriate prefix or postfix */
  prefix(a);
  putchar(cha[(long)a]);
  postfix(a);
}  /* makechar */


void dolmove_correct(node *p)
{  /* get final states for intermediate nodes */
  long i;
  long z0, z1, s0, s1, temp;

  if (p->tip)
    return;
  for (i = 0; i < (words); i++) {
    if (p->back == NULL) {
      s0 = zeroanc[i];
      s1 = oneanc[i];
    } else {
      s0 = treenode[p->back->index - 1]->statezero[i];
      s1 = treenode[p->back->index - 1]->stateone[i];
    }
    z0 = (s0 & p->statezero[i]) |
         (p->next->back->statezero[i] & p->next->next->back->statezero[i]);
    z1 = (s1 & p->stateone[i]) |
         (p->next->back->stateone[i] & p->next->next->back->stateone[i]);
    if (dollo) {
      temp = z0 & (~(zeroanc[i] & z1));
      z1 &= ~(oneanc[i] & z0);
      z0 = temp;
    }
    temp = fullset & (~z0) & (~z1);
    p->statezero[i] = z0 | (temp & s0 & (~s1));
    p->stateone[i] = z1 | (temp & s1 & (~s0));
  }
}  /* dolmove_correct */


void dolmove_count(node *p)
{
  /* counts the number of steps in a fork of the tree.
    The program spends much of its time in this PROCEDURE */
  long i, j, l;
  bitptr steps;

  steps = (bitptr)Malloc(words*sizeof(long));
  if (dollo) {
    for (i = 0; i < (words); i++)
      steps[i] = (treenode[p->back->index - 1]->stateone[i] &
                  p->statezero[i] & zeroanc[i]) |
                 (treenode[p->back->index - 1]->statezero[i] &
                  p->stateone[i] & oneanc[i]);
  } else {
    for (i = 0; i < (words); i++)
      steps[i] = treenode[p->back->index - 1]->stateone[i] &
                 treenode[p->back->index - 1]->statezero[i] & p->stateone[i] &
                 p->statezero[i];
  }
  j = 1;
  l = 0;
  for (i = 0; i < (chars); i++) {
    l++;
    if (l > bits) {
      l = 1;
      j++;
    }
    if (((1L << l) & steps[j - 1]) != 0) {
      if (((1L << l) & zeroanc[j - 1]) != 0)
        numszero[i] += weight[i];
      else
        numsone[i] += weight[i];
    }
  }
  free(steps);
}  /* dolmove_count */


void preorder(node *p)
{
  /* go back up tree setting up and counting interior node
    states */

  if (!p->tip) {
    dolmove_correct(p);
    preorder(p->next->back);
    preorder(p->next->next->back);
  }
  if (p->back != NULL)
    dolmove_count(p);
}  /* preorder */


void evaluate(node *r)
{
  /* Determines the number of losses or polymorphisms needed for a tree.
     This is the minimum number needed to evolve chars on this tree */
  long i, stepnum, smaller;
  double sum;
  boolean nextcompat, thiscompat, done;

  sum = 0.0;
  for (i = 0; i < (chars); i++) {
    numszero[i] = 0;
    numsone[i] = 0;
  }
  for (i = 0; i < (words); i++) {
    zeroanc[i] = fullset;
    oneanc[i] = 0;
  }
  compatible = 0;
  nextcompat = true;
  postorder(r);
  preorder(r);
  for (i = 0; i < (words); i++) {
    zeroanc[i] = 0;
    oneanc[i] = fullset;
  }
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
    thiscompat = (stepnum <= weight[i]);
    if (factors) {
      done = (i + 1 == chars);
      if (!done)
        done = (factor[i + 1] != factor[i]);
      nextcompat = (nextcompat && thiscompat);
      if (done) {
        if (nextcompat)
          compatible += weight[i];
        nextcompat = true;
      }
    } else if (thiscompat)
      compatible += weight[i];
    guess[i] = '?';
    if (!ancone[i] ||
        (anczero[i] && numszero[i] < numsone[i]))
      guess[i] = '0';
    else if (!anczero[i] ||
             (ancone[i] && numsone[i] < numszero[i]))
      guess[i] = '1';
  }
  like = -sum;
}  /* evaluate */


void reroot(node *outgroup)
{
  /* reorients tree, putting outgroup in desired position. */
  node *p, *q, *newbottom, *oldbottom;
  boolean onleft;

  if (outgroup->back->index == root->index)
    return;
  newbottom = outgroup->back;
  p = treenode[newbottom->index - 1]->back;
  while (p->index != root->index) {
    oldbottom = treenode[p->index - 1];
    treenode[p->index - 1] = p;
    p = oldbottom->back;
  }
  onleft = (p == root->next);
  if (restoring)
    if (!onleft && wasleft){
      p = root->next->next;
      q = root->next;
    } else {
      p = root->next;
      q = root->next->next;
    }
  else {
    if (onleft)
      oldoutgrno = root->next->next->back->index;
    else
      oldoutgrno = root->next->back->index;
    wasleft = onleft;
    p = root->next;
    q = root->next->next;
  }
  p->back->back = q->back;
  q->back->back = p->back;
  p->back = outgroup;
  q->back = outgroup->back;
  if (restoring) {
    if (!onleft && wasleft) {
      outgroup->back->back = root->next;
      outgroup->back = root->next->next;
    } else {
      outgroup->back->back = root->next->next;
      outgroup->back = root->next;
    }
  } else {
    outgroup->back->back = root->next->next;
    outgroup->back = root->next;
  }
  treenode[newbottom->index - 1] = newbottom;
}  /* reroot */


void dolmove_hyptrav(node *r)
{
  /* compute states at interior nodes for one character */
  if (!r->tip)
    dolmove_correct(r);
  if (((1L << dispbit) & r->stateone[dispword - 1]) != 0) {
    if (((1L << dispbit) & r->statezero[dispword - 1]) != 0) {
      if (dollo)
        r->state = '?';
      else
        r->state = 'P';
    } else
      r->state = '1';
  } else {
    if (((1L << dispbit) & r->statezero[dispword - 1]) != 0)
      r->state = '0';
    else
      r->state = '?';
  }
  if (!r->tip) {
    dolmove_hyptrav(r->next->back);
    dolmove_hyptrav(r->next->next->back);
  }
}  /* dolmove_hyptrav */


void dolmove_hypstates()
{
  /* fill in and describe states at interior nodes */
  long i, j, k;

  for (i = 0; i < (words); i++) {
    zeroanc[i] = 0;
    oneanc[i] = 0;
  }
  for (i = 0; i < (chars); i++) {
    j = i / bits + 1;
    k = i % bits + 1;
    if (guess[i] == '0')
      zeroanc[j - 1] = ((long)zeroanc[j - 1]) | (1L << k);
    if (guess[i] == '1')
      oneanc[j - 1] = ((long)oneanc[j - 1]) | (1L << k);
  }
  filltrav(root);
  dolmove_hyptrav(root);
}  /* dolmove_hypstates */


void grwrite(chartype c, long num, long *pos)
{
  int i;

  prefix(c);
  for (i = 1; i <= num; i++) {
    if ((*pos) >= leftedge && (*pos) - leftedge + 1 < screenwidth)
      putchar(cha[(long)c]);
    (*pos)++;
  }
  postfix(c);
}  /* grwrite */

 
void dolmove_drawline(long i)
{
  /* draws one row of the tree diagram by moving up tree */
  node *p, *q, *r, *first =NULL, *last =NULL;
  long n, j, pos;
  boolean extra, done;
  Char s, cc;
  chartype c, d;

  pos = 1;
  p = nuroot;
  q = nuroot;
  extra = false;
  if (i == (long)p->ycoord && (p == root || subtree)) {
    c = overt;
    if (display) {
      if (p == root)
        cc = guess[dispchar - 1];
      else
        cc = p->state;
      switch (cc) {

      case '1':
        c = onne;
        break;

      case '0':
        c = zerro;
        break;

      case '?':
        c = question;
        break;

      case 'P':
        c = polym;
        break;
      }
    }
    
    if ((subtree))
      stwrite("Subtree:", 8, &pos, leftedge, screenwidth);
    if (p->index >= 100)
      nnwrite(p->index, 3, &pos, leftedge, screenwidth);
    else if (p->index >= 10) {
      grwrite(c, 1, &pos);
      nnwrite(p->index, 2, &pos, leftedge, screenwidth);
    } else {
      grwrite(c, 2, &pos);
      nnwrite(p->index, 1, &pos, leftedge, screenwidth);
    }
    extra = true;
  } else {
    if (subtree)
      stwrite("          ", 10, &pos, leftedge, screenwidth);
    else
      stwrite("  ", 2, &pos, leftedge, screenwidth);
  }
  do {
    if (!p->tip) {
      r = p->next;
      done = false;
      do {
        if (i >= r->back->ymin && i <= r->back->ymax) {
          q = r->back;
          done = true;
        }
        r = r->next;
      } while (!(done || r == p));
      first = p->next->back;
      r = p->next;
      while (r->next != p)
        r = r->next;
      last = r->back;
    }
    done = (p == q);
    n = (long)p->xcoord - (long)q->xcoord;
    if (n < 3 && !q->tip)
      n = 3;
    if (extra) {
      n--;
      extra = false;
    }
    if ((long)q->ycoord == i && !done) {
      if ((long)q->ycoord > (long)p->ycoord)
        d = upcorner;
      else
        d = downcorner;
      c = overt;
      s = q->state;
      if (s == 'P' && p->state != 'P')
        s = p->state;
      if (display) {
        switch (s) {

        case '1':
          c = onne;
          break;

        case '0':
          c = zerro;
          break;

        case '?':
          c = question;
          break;

        case 'P':
          c = polym;
          break;
        }
        d = c;
      }
      if (n > 1) {
        grwrite(d, 1, &pos);
        grwrite(c, n - 3, &pos);
      }
      if (q->index >= 100)
        nnwrite(q->index, 3, &pos, leftedge, screenwidth);
      else if (q->index >= 10) {
        grwrite(c, 1, &pos);
        nnwrite(q->index, 2, &pos, leftedge, screenwidth);
      } else {
        grwrite(c, 2, &pos);
        nnwrite(q->index, 1, &pos, leftedge, screenwidth);
      }
      extra = true;
    } else if (!q->tip) {
      if ((long)last->ycoord > i && (long)first->ycoord < i
        && i != (long)p->ycoord) {
        c = up;
        if (i < (long)p->ycoord)
          s = p->next->back->state;
        else
          s = p->next->next->back->state;
        if (s == 'P' && p->state != 'P')
          s = p->state;
        if (display) {
          switch (s) {

          case '1':
            c = onne;
            break;

          case '0':
            c = zerro;
            break;

          case '?':
            c = question;
            break;

          case 'P':
            c = polym;
            break;
          }
        }
        grwrite(c, 1, &pos);
        chwrite(' ', n - 1, &pos, leftedge, screenwidth);
      } else
        chwrite(' ', n, &pos, leftedge, screenwidth);
    } else
      chwrite(' ', n, &pos, leftedge, screenwidth);
    if (p != q)
      p = q;
  } while (!done);
  if ((long)p->ycoord == i && p->tip) {
    n = 0;
    for (j = 1; j <= nmlngth; j++) {
      if (nayme[p->index - 1][j - 1] != '\0')
        n = j;
    }
    chwrite(':', 1, &pos, leftedge, screenwidth);
    for (j = 0; j < n; j++)
      chwrite(nayme[p->index - 1][j], 1, &pos, leftedge, screenwidth);
  }
  putchar('\n');
}  /* dolmove_drawline */


void dolmove_printree()
{
  /* prints out diagram of the tree */
  long tipy;
  long i, dow;

  if (!subtree)
    nuroot = root;
  if (changed || newtree)
    evaluate(root);
  if (display)
    dolmove_hypstates();
#ifdef WIN32
  if(ibmpc || ansi){
    phyClearScreen();
  } else {
    printf("\n");
  }
#else
  if (ansi || ibmpc)
    printf("\033[2J\033[H");
  else
    putchar('\n');
#endif
  tipy = 1;
  dow = down;
  if (spp * dow > screenlines && !subtree)
    dow--;
  printf("(unrooted)");
  if (display) {
    printf(" ");
    makechar(onne);
    printf(":1 ");
    makechar(question);
    printf(":? ");
    makechar(zerro);
    printf(":0 ");
    makechar(polym);
    printf(":0/1");
  } else
    printf("                    ");
  if (!earlytree) {
    printf("%10.1f Steps", -like);
  }
  if (display)
    printf("  SITE%4ld", dispchar);
  else
    printf("         ");
  if (!earlytree) {
    printf("  %3ld chars compatible\n", compatible);
  }

  printf("%-20s",dollo ? "Dollo" : "Polymorphism");
  
  if (changed && !earlytree) {
    if (-like < bestyet) {
      printf("     BEST YET!");
      bestyet = -like;
    } else if (fabs(-like - bestyet) < 0.000001)
      printf("     (as good as best)");
    else {
      if (-like < gotlike)
        printf("     better");
      else if (-like > gotlike)
        printf("     worse!");
    }
  }
  printf("\n");
  farthest = 0;
  coordinates(nuroot, &tipy, 1.5, &farthest);
  vmargin = 5;
  treelines = tipy - dow;
  if (topedge != 1){
        printf("** %ld lines above screen **\n", topedge - 1);
    vmargin++;}
  if ((treelines - topedge + 1) > (screenlines - vmargin))
    vmargin++;
  for (i = 1; i <= treelines; i++) {
    if (i >= topedge && i < topedge + screenlines - vmargin)
      dolmove_drawline(i);
  }
  if ((treelines - topedge + 1) > (screenlines - vmargin)) 
    printf("** %ld lines below screen **\n",
           treelines - (topedge - 1 + screenlines - vmargin));
  if (treelines - topedge + vmargin + 1 < screenlines)
    putchar('\n');
  gotlike = -like;
  changed = false;
}  /* dolmove_printree */


void arbitree()
{
  long i;

  root = treenode[0];
  add2(treenode[0], treenode[1], treenode[spp], &root, restoring, wasleft,
         treenode);
  for (i = 3; i <= (spp); i++)
    add2(treenode[spp + i - 3], treenode[i - 1], treenode[spp + i - 2], &root,
      restoring, wasleft, treenode);
  for (i = 0; i < (nonodes); i++)
    in_tree[i] = true;
}  /* arbitree */


void yourtree()
{
  long i, j;
  boolean ok;

  root = treenode[0];
  add2(treenode[0], treenode[1], treenode[spp], &root, restoring, wasleft,
         treenode);
  i = 2;
  do {
    i++;
    dolmove_printree();
    printf("Add species%3ld: ", i);
    for (j = 0; j < nmlngth; j++)
      putchar(nayme[i - 1][j]);
    do {
      printf("\nbefore node (type number): ");
      inpnum(&j, &ok);
      ok = (ok && ((j >= 1 && j < i) || (j > spp && j < spp + i - 1)));
      if (!ok)
        printf("Impossible number. Please try again:\n");
    } while (!ok);
    add2(treenode[j - 1], treenode[i - 1], treenode[spp + i - 2], &root,
      restoring, wasleft, treenode);
  } while (i != spp);
  for (i = 0; i < (nonodes); i++)
    in_tree[i] = true;
}  /* yourtree */


void initdolmovenode(node **p, node **grbg, node *q, long len, long nodei,
                        long *ntips, long *parens, initops whichinit,
                        pointarray treenode, pointarray nodep, Char *str, Char *ch,
                        FILE *intree)
{
  /* initializes a node */
  /* LM 7/27  I added this function and the commented lines around */
  /* treeread() to get the program running, but all 4 move programs*/
  /* are improperly integrated into the v4.0 support files.  As is */
  /* this is a patchwork function                                   */
  boolean minusread;
  double valyew, divisor;

  switch (whichinit) {
  case bottom:
    gnutreenode(grbg, p, nodei, chars, zeros);
    treenode[nodei - 1] = *p;
    break;
  case nonbottom:
    gnutreenode(grbg, p, nodei, chars, zeros);
    break;
  case tip:
    match_names_to_data (str, treenode, p, spp);
    break;
  case length:
    processlength(&valyew, &divisor, ch, &minusread, intree, parens);
    /* process lengths and discard */
  default:      /*cases hslength,hsnolength,treewt,unittrwt,iter,*/
    break;      /*length should never occur                      */
  }
} /* initdolmovenode */


void buildtree()
{
  long i, j, nextnode;
  node *p;

  changed = false;
  newtree = false;
  switch (how) {

  case arb:
    arbitree();
    break;

  case use:
    /* Open in binary: ftell() is broken for UNIX line-endings under WIN32 */
    openfile(&intree,INTREE,"input tree file", "rb",progname,intreename);
    names = (boolean *)Malloc(spp*sizeof(boolean));
    firsttree = true;                       /**/
    nodep = NULL;                           /**/
    nextnode = 0;                           /**/
    haslengths = 0;                         /**/
    zeros = (long *)Malloc(chars*sizeof(long));         /**/
    for (i = 0; i < chars; i++)             /**/
      zeros[i] = 0;                         /**/
    treeread(intree, &root, treenode, &goteof, &firsttree, nodep, &nextnode,
        &haslengths, &grbg, initdolmovenode,false,nonodes);
    for (i = spp; i < (nonodes); i++) {
      p = treenode[i];
      for (j = 1; j <= 3; j++) {
        p->stateone = (bitptr)Malloc(words*sizeof(long));
        p->statezero = (bitptr)Malloc(words*sizeof(long));
        p = p->next;
      }
    } /* debug: see comment at initdolmovenode() */

    /*treeread(which, ch, &root, treenode, names);*/
    for (i = 0; i < (spp); i++)
      in_tree[i] = names[i];
    free(names);
    FClose(intree);
    break;

  case spec:
    yourtree();
    break;
  }
  outgrno = root->next->back->index;
  if (in_tree[outgrno - 1])
    reroot(treenode[outgrno - 1]);
}  /* buildtree */


void rearrange()
{
  long i, j;
  boolean ok1, ok2;
  node *p, *q;

  printf("Remove everything to the right of which node? ");
  inpnum(&i, &ok1);
  ok1 = (ok1 && i >= 1 && i < spp * 2 && i != root->index);
  if (ok1) {
    printf("Add before which node? ");
    inpnum(&j, &ok2);
    ok2 = (ok2 && j >= 1 && j < spp * 2);
    if (ok2) {
      ok2 = (treenode[j - 1] != treenode[treenode[i - 1]->back->index - 1]);
      p = treenode[j - 1];
      while (p != root) {
        ok2 = (ok2 && p != treenode[i - 1]);
        p = treenode[p->back->index - 1];
      }
      if (ok1 && ok2) {
        what = i;
        q = treenode[treenode[i - 1]->back->index - 1];
        if (q->next->back->index == i)
          fromwhere = q->next->next->back->index;
        else
          fromwhere = q->next->back->index;
        towhere = j;
        re_move2(&treenode[i - 1], &q, &root, &wasleft, treenode);
        add2(treenode[j - 1], treenode[i - 1], q, &root, restoring, wasleft, treenode);
      }
      lastop = rearr;
    }
  }
  changed = (ok1 && ok2);
  dolmove_printree();
  if (!(ok1 && ok2))
    printf("Not a possible rearrangement.  Try again: \n");
  else {
    oldwritten = written;
    written = false;
  }
}  /* rearrange */


void tryadd(node *p, node **item, node **nufork, double *place)
{
  /* temporarily adds one fork and one tip to the tree.
    Records scores in ARRAY place */
  add2(p, *item, *nufork, &root, restoring, wasleft, treenode);
  evaluate(root);
  place[p->index - 1] = -like;
  re_move2(item, nufork, &root, &wasleft, treenode);
}  /* tryadd */


void addpreorder(node *p, node *item_, node *nufork_, double *place)
{
  /* traverses a binary tree, calling PROCEDURE tryadd
    at a node before calling tryadd at its descendants */
  node *item, *nufork;


  item = item_;
  nufork = nufork_;
  if (p == NULL)
    return;
  tryadd(p, &item,&nufork,place);
  if (!p->tip) {
    addpreorder(p->next->back, item,nufork,place);
    addpreorder(p->next->next->back,item,nufork,place);
  }
}  /* addpreorder */


void try()
{
  /* Remove node, try it in all possible places */
  double *place;
  long i, j, oldcompat;
  double current;
  node *q, *dummy, *rute;
  boolean tied, better, ok;

  printf("Try other positions for which node? ");
  inpnum(&i, &ok);
  if (!(ok && i >= 1 && i <= nonodes && i != root->index)) {
    printf("Not a possible choice! ");
    return;
  }
  printf("WAIT ...\n");
  place = (double *)Malloc(nonodes*sizeof(double));
  for (j = 0; j < (nonodes); j++)
    place[j] = -1.0;
  evaluate(root);
  current = -like;
  oldcompat = compatible;
  what = i;
  q = treenode[treenode[i - 1]->back->index - 1];
  if (q->next->back->index == i)
    fromwhere = q->next->next->back->index;
  else
    fromwhere = q->next->back->index;
  rute = root;
  if (root->index == treenode[i - 1]->back->index) {
    if (treenode[treenode[i - 1]->back->index - 1]->next->back == treenode[i - 1])
      rute = treenode[treenode[i - 1]->back->index - 1]->next->next->back;
    else
      rute = treenode[treenode[i - 1]->back->index - 1]->next->back;
  }
  re_move2(&treenode[i - 1], &dummy, &root, &wasleft, treenode);
  oldleft = wasleft;
  root = rute;
  addpreorder(root, treenode[i - 1], dummy, place);
  wasleft = oldleft;
  restoring = true;
  add2(treenode[fromwhere - 1], treenode[what - 1], dummy, &root,
    restoring, wasleft, treenode);
  like = -current;
  compatible = oldcompat;
  restoring = false;
  better = false;
  printf("       BETTER: ");
  for (j = 1; j <= (nonodes); j++) {
    if (place[j - 1] < current && place[j - 1] >= 0.0) {
      printf("%3ld:%6.2f", j, place[j - 1]);
      better = true;
    }
  }
  if (!better)
    printf(" NONE");
  printf("\n       TIED:    ");
  tied = false;
  for (j = 1; j <= (nonodes); j++) {
    if (fabs(place[j - 1] - current) < 1.0e-6 && j != fromwhere) {
      if (j < 10)
        printf("%2ld", j);
      else
        printf("%3ld", j);
      tied = true;
    }
  }
  if (tied)
    printf(":%6.2f\n", current);
  else
    printf("NONE\n");
  changed = true;
  free(place);
}  /* try */

void undo()
{
  /* restore to tree before last rearrangement */
  long temp;
  boolean btemp;
  node *q;

  switch (lastop) {

  case rearr:
    restoring = true;
    oldleft = wasleft;
    re_move2(&treenode[what - 1], &q, &root, &wasleft, treenode);
    btemp = wasleft;
    wasleft = oldleft;
    add2(treenode[fromwhere - 1], treenode[what - 1], q, &root,
      restoring, wasleft, treenode);
    wasleft = btemp;
    restoring = false;
    temp = fromwhere;
    fromwhere = towhere;
    towhere = temp;
    changed = true;
    break;

  case flipp:
    q = treenode[atwhat - 1]->next->back;
    treenode[atwhat - 1]->next->back =
      treenode[atwhat - 1]->next->next->back;
    treenode[atwhat - 1]->next->next->back = q;
    treenode[atwhat - 1]->next->back->back = treenode[atwhat - 1]->next;
    treenode[atwhat - 1]->next->next->back->back =
      treenode[atwhat - 1]->next->next;
    break;

  case reroott:
    restoring = true;
    temp = oldoutgrno;
    oldoutgrno = outgrno;
    outgrno = temp;
    reroot(treenode[outgrno - 1]);
    restoring = false;
    break;

  case none:
    /* blank case */
    break;
  }
  dolmove_printree();
  if (lastop == none) {
    printf("No operation to undo! ");
    return;
  }
  btemp = oldwritten;
  oldwritten = written;
  written = btemp;
}  /* undo */


void treewrite(boolean done)
{
  /* write out tree to a file */
  Char ch;

  treeoptions(waswritten, &ch, &outtree, outtreename, progname);
  if (!done)
    dolmove_printree();
  if (waswritten && ch == 'N')
    return;
  col = 0;
  treeout(root, 1, &col, root);
  printf("\nTree written to file \"%s\"\n\n", outtreename);
  waswritten = true;
  written = true;
  FClose(outtree);
#ifdef MAC
  fixmacfile(outtreename);
#endif
}  /* treewrite */


void clade()
{
  /* pick a subtree and show only that on screen */
  long i;
  boolean ok;

  printf("Select subtree rooted at which node (0 for whole tree)? ");
  inpnum(&i, &ok);
  ok = (ok && ((unsigned)i) <= ((unsigned)nonodes));
  if (ok) {
    subtree = (i > 0);
    if (subtree)
      nuroot = treenode[i - 1];
    else
      nuroot = root;
  }
  dolmove_printree();
  if (!ok)
    printf("Not possible to use this node. ");
}  /* clade */


void flip()
{
  /* flip at a node left-right */
  long i;
  boolean ok;
  node *p;

  printf("Flip branches at which node? ");
  inpnum(&i, &ok);
  ok = (ok && i > spp && i <= nonodes);
  if (ok) {
    p = treenode[i - 1]->next->back;
    treenode[i - 1]->next->back = treenode[i - 1]->next->next->back;
    treenode[i - 1]->next->next->back = p;
    treenode[i - 1]->next->back->back = treenode[i - 1]->next;
    treenode[i - 1]->next->next->back->back = treenode[i - 1]->next->next;
    atwhat = i;
    lastop = flipp;
  }
  dolmove_printree();
  if (ok) {
    oldwritten = written;
    written = false;
    return;
  }
  if (i >= 1 && i <= spp)
    printf("Can't flip there. ");
  else
    printf("No such node. ");
}  /* flip */


void changeoutgroup()
{
  long i;
  boolean ok;

  oldoutgrno = outgrno;
  do {
    printf("Which node should be the new outgroup? ");
    inpnum(&i, &ok);
    ok = (ok && in_tree[i - 1] && i >= 1 && i <= nonodes &&
          i != root->index);
    if (ok)
      outgrno = i;
  } while (!ok);
  if (in_tree[outgrno - 1])
    reroot(treenode[outgrno - 1]);
  changed = true;
  lastop = reroott;
  dolmove_printree();
  oldwritten = written;
  written = false;
}  /* changeoutgroup */


void redisplay()
{
  boolean done;
  char input[100];

  done = false;
  waswritten = false;
  do {
    printf("NEXT? (Options: R # + - S . T U W O F H J K L C ? X Q) ");
    printf("(? for Help) ");
#ifdef WIN32
    phyFillScreenColor();
#endif
    getstryng(input);
    ch = input[0];
    uppercase(&ch);
    if (strchr("RSWH#.O?+TFX-UCQHJKL",ch) != NULL){
      switch (ch) {

      case 'R':
        rearrange();
        break;

      case '#':
        nextinc(&dispchar, &dispword, &dispbit, chars, bits, &display,
                  numsteps, weight);
        dolmove_printree();
        break;

      case '+':
        nextchar(&dispchar, &dispword, &dispbit, chars, bits, &display);
        dolmove_printree();
        break;

      case '-':
        prevchar(&dispchar, &dispword, &dispbit, chars, bits, &display);
        dolmove_printree();
        break;

      case 'S':
        show(&dispchar, &dispword, &dispbit, chars, bits, &display);
        dolmove_printree();
        break;

      case '.':
        dolmove_printree();
        break;

      case 'T':
        try();
        break;

      case 'U':
        undo();
        break;

      case 'W':
        treewrite(done);
        break;

      case 'O':
        changeoutgroup();
        break;

      case 'F':
        flip();
        break;

      case 'H':
        window(left, &leftedge, &topedge, hscroll, vscroll, treelines,
                 screenlines, screenwidth, farthest, subtree);
        dolmove_printree();
        break;

      case 'J':
        window(downn, &leftedge, &topedge, hscroll, vscroll, treelines,
                 screenlines, screenwidth, farthest, subtree);
        dolmove_printree();
        break;

      case 'K':
        window(upp, &leftedge, &topedge, hscroll, vscroll, treelines,
                 screenlines, screenwidth, farthest, subtree);
        dolmove_printree();
        break;

      case 'L':
        window(right, &leftedge, &topedge, hscroll, vscroll, treelines,
                 screenlines, screenwidth, farthest, subtree);
        dolmove_printree();
        break;

      case 'C':
        clade();
        break;

      case '?':
        help("character");
        dolmove_printree();
        break;

      case 'X':
        done = true;
        break;

      case 'Q':
        done = true;
        break;
      }
    }
  } while (!done);
  if (!written) {
    do {
      printf("Do you want to write out the tree to a file? (Y or N) ");
#ifdef WIN32
      phyFillScreenColor();
#endif
      getstryng(input);
      ch = input[0];
    } while (ch != 'Y' && ch != 'y' && ch != 'N' && ch != 'n');
  }
  if (ch == 'Y' || ch == 'y')
    treewrite(done);
}  /* redisplay */


void treeconstruct()
{
  /* constructs a binary tree from the pointers in treenode. */

  restoring = false;
  subtree = false;
  display = false;
  dispchar = 0;
  fullset = (1L << (bits + 1)) - (1L << 1);
  guess = (Char *)Malloc(chars*sizeof(Char));
  numsteps = (steptr)Malloc(chars*sizeof(long));
  earlytree = true;
  buildtree();
  waswritten = false;
  printf("\nComputing steps needed for compatibility in sites ...\n\n");
  newtree = true;
  earlytree = false;
  dolmove_printree();
  bestyet = -like;
  gotlike = -like;
  lastop = none;
  newtree = false;
  written = false;
  lastop = none;
  redisplay();
}  /* treeconstruct */


int main(int argc, Char *argv[])
{ /* Interactive Dollo/polymorphism parsimony */
  /* reads in spp, chars, and the data. Then calls treeconstruct to
    construct the tree and query the user */
#ifdef MAC
  argc = 1;                /* macsetup("Dolmove","");                */
  argv[0] = "Dolmove";
#endif
  init(argc, argv);
  progname = argv[0];
  strcpy(infilename,INFILE);
  strcpy(outtreename,OUTTREE);
  strcpy(intreename,INTREE);

  openfile(&infile,infilename,"input file", "r",argv[0],infilename);

  screenlines = 24;
  scrollinc   = 20;
  screenwidth = 80;
  topedge     = 1;
  leftedge    = 1;
  ibmpc = IBMCRT;
  ansi = ANSICRT;
  root = NULL;
  bits = 8*sizeof(long) - 1;
  doinput();
  configure();
  treeconstruct();
  if (waswritten) {
    FClose(outtree);
#ifdef MAC
    fixmacfile(outtreename);
#endif
  }
  FClose(infile);
#ifdef WIN32
  phyRestoreConsoleAttributes();
#endif
  return 0;
}  /* Interactive Dollo/polymorphism parsimony */

