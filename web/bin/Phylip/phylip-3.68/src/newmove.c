
#include "phylip.h"
#include "disc.h"
#include "moves.h"
#include "wagner.h"

/* version 3.6. (c) Copyright 1993-2002 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define overr           4
#define which           1
#define maxsz           999   /* size of pointer array for the undo trees */
                              /* this can be large without eating memory */
typedef struct treeset_t 
{
  
  node *root;
  pointarray nodep;
  pointarray treenode;
  long nonodes;
  boolean waswritten, hasmult, haslengths, nolengths, initialized;
}
treeset_t;

treeset_t treesets[2];

node **treeone, **treetwo;


typedef enum {
  horiz, vert, up, overt, upcorner, midcorner, downcorner, 
  onne, zerro, question
} chartype;

typedef enum {
  arb, use, spec
} howtree;

typedef enum {
  rearr, flipp, reroott, none
} rearrtype;

typedef enum {beforenode, atnode} movet;

movet fromtype;

#ifndef OLDC
/*function prototypes */
void   getoptions(void);
void   inputoptions(void);
void   allocrest(void);
void   doinput(void);
void   configure(void);
void   prefix(chartype);
void   postfix(chartype);
void   makechar(chartype);
void   move_fillin(node *);
void   move_postorder(node *);

void   evaluate(node *);
void   reroot(node *);
void   move_filltrav(node *);
void   move_hyptrav(node *);
void   move_hypstates(void);
void   grwrite(chartype, long, long *);
void   move_drawline(long, long);
void   move_printree(void);
void   arbitree(void);
void   yourtree(void);

void   initmovenode(node **, node **, node *, long, long, long *, long *,
        initops, pointarray, pointarray, Char *, Char *, FILE *);
void   buildtree(void);
void   rearrange(void);
void   tryadd(node *, node *, node *, double *);
void   addpreorder(node *, node *, node *, double *);
void   try(void);
void   undo(void);
void   treewrite(boolean);
void   clade(void);
void   flip(long);

void   changeoutgroup(void);
void   redisplay(void);
void   treeconstruct(void);
void   maketriad(node **, long);
/*function prototypes */
#endif
 
char infilename[FNMLNGTH],intreename[FNMLNGTH],outtreename[FNMLNGTH], weightfilename[FNMLNGTH], ancfilename[FNMLNGTH], mixfilename[FNMLNGTH], factfilename[FNMLNGTH];
node *root;
long outgrno, screenlines, col, treelines, leftedge, topedge,
  vmargin, hscroll, vscroll, scrollinc, screenwidth, farthest, whichtree, othertree;
/* outgrno indicates outgroup */
boolean weights, thresh, outgropt, ancvar, questions, allsokal,
               allwagner, mixture, factors, noroot,  waswritten;
boolean *ancone, *anczero, *ancone0, *anczero0;
Char *factor;
pointptr treenode;   /* pointers to all nodes in tree */
double threshold;
double *threshwt;
bitptr wagner, wagner0;
unsigned char che[9];
boolean reversed[9];
boolean graphic[9];
howtree how;
gbit *garbage;
char* progname;
Char ch;

/* Variables for treeread */
boolean usertree, goteof, firsttree, haslengths;
pointarray nodep;
node *grbg;
long *zeros;

/* Local variables for treeconstruct, propagated globally for C vesion: */
long dispchar, dispword, dispbit, atwhat, what, fromwhere, towhere,
     oldoutgrno, compatible;
double like, bestyet, gotlike;
Char *guess;
boolean display, newtree, changed, subtree, written, oldwritten, restoring,
  wasleft, oldleft, earlytree;
boolean *in_tree;
steptr numsteps, numsone, numszero;
long fullset;
bitptr steps, zeroanc, oneanc;
node *nuroot;
rearrtype lastop;
boolean *names;
steptr steps_a;


void maketriad(node **p, long index)
{
  /* Initiate an internal node with stubs for two children */
  long i, j;
  node *q;
  q = NULL;
  for (i = 1; i <= 3; i++) {
    gnu(&grbg, p);
    (*p)->index = index;
    (*p)->hasname = false;
    (*p)->haslength = false;
    (*p)->deleted=false;
    (*p)->deadend=false;
    (*p)->onebranch=false;
    (*p)->onebranchhaslength=false;
    if(!(*p)->statezero)
      (*p)->statezero = (bitptr)Malloc(words*sizeof(long));
    if(!(*p)->stateone)
      (*p)->stateone = (bitptr)Malloc(words*sizeof(long));
    if(!(*p)->zerosteps)
      (*p)->zerosteps = (bitptr)Malloc(words*sizeof(long));
    if(!(*p)->onesteps)
      (*p)->onesteps = (bitptr)Malloc(words*sizeof(long));
    for (j=0;j<MAXNCH;j++)
      (*p)->nayme[j] = '\0';
    (*p)->next = q;
    q = *p;
  }
  (*p)->next->next->next = *p;
  q = (*p)->next;
  while (*p != q) {
    (*p)->back = NULL;
    (*p)->tip = false;
    *p = (*p)->next;
  }
  treenode[index - 1] = *p;
}  /* maketriad */


void move_copynode(node *fromnode, node *tonode)
{
  
  /* Copy the contents of a node from fromnode to tonode. */
  int i = 0;
  
  /*
  printf("copynode: fromnode->numnuc = %ld, tonode->numnuc = %ld\n",
    fromnode->numnuc,tonode->numnuc);
  */
  if (fromnode->statezero != NULL)
    memcpy(tonode->statezero, fromnode->statezero, words*sizeof(long));

  if (fromnode->stateone != NULL)
    memcpy(tonode->stateone, fromnode->stateone, words*sizeof(long));

  if (fromnode->zerosteps != NULL)
    memcpy(tonode->zerosteps, fromnode->zerosteps, words*sizeof(long));

  if (fromnode->onesteps != NULL)
    memcpy(tonode->onesteps, fromnode->onesteps, words*sizeof(long));


  tonode->state   = fromnode->state;
  tonode->index   = fromnode->index;
  tonode->tip     = fromnode->tip;
  for (i=0;i<MAXNCH;i++)
    tonode->nayme[i] = fromnode->nayme[i];
}
/* move_copynode */


node *copytrav(node *p)
{
  /* Traverse the tree from p on down, copying nodes to the other tree */
  node *q, *newnode, *newnextnode, *temp;
  gnu(&grbg, &newnode);
  if(!newnode->statezero)
    newnode->statezero = (bitptr)Malloc(words*sizeof(long));
  if(!newnode->stateone)
    newnode->stateone = (bitptr)Malloc(words*sizeof(long));
  if(!newnode->zerosteps)
    newnode->zerosteps = (bitptr)Malloc(words*sizeof(long));
  if(!newnode->onesteps)
    newnode->onesteps = (bitptr)Malloc(words*sizeof(long));
  move_copynode(p,newnode);

  if (treenode[p->index-1] == p)
    treesets[othertree].treenode[p->index-1] = newnode;

  /* if this is a tip, return now */
  if (p->tip)
    return newnode;

  /* go around the ring, copying as we go */

  q = p->next;
  gnu(&grbg, &newnextnode);

  if(!newnextnode->statezero)
    newnextnode->statezero = (bitptr)Malloc(words*sizeof(long));
  if(!newnextnode->stateone)
    newnextnode->stateone = (bitptr)Malloc(words*sizeof(long));
  if(!newnextnode->zerosteps)
    newnextnode->zerosteps = (bitptr)Malloc(words*sizeof(long));
  if(!newnextnode->onesteps)
    newnextnode->onesteps = (bitptr)Malloc(words*sizeof(long));
  move_copynode(q, newnextnode);
  newnode->next = newnextnode;
  do {
    newnextnode->back = copytrav(q->back);
    newnextnode->back->back = newnextnode;
    q = q->next;
    if (q == p)
      newnextnode->next = newnode;
    else {
      temp = newnextnode;
      gnu(&grbg, &newnextnode);
      if(!newnextnode->statezero)
        newnextnode->statezero = (bitptr)Malloc(words*sizeof(long));
      if(!newnextnode->stateone)
        newnextnode->stateone = (bitptr)Malloc(words*sizeof(long));
      if(!newnextnode->zerosteps)
        newnextnode->zerosteps = (bitptr)Malloc(words*sizeof(long));
      if(!newnextnode->onesteps)
        newnextnode->onesteps = (bitptr)Malloc(words*sizeof(long));
      move_copynode(q, newnextnode);
      temp->next = newnextnode;
    }
  } while (q != p);
  return newnode;
}
/* copytrav */

void chucktree(node *p)
{
  /* recursively run through a tree and chuck all of its nodes,
     putting them on the garbage list */
  int i, numNodes = 1;
  node *q, *r;

  /* base case -- tip */
  if(p->tip){
    chuck(&grbg, p);
    return;
  }

  /* recursively callchuck tree on all decendants */
  q = p->next;
  while(q != p){
    chucktree(q->back);
    numNodes++;
    q = q->next;
  }

  /* now chuck all sub-nodes in the node ring */
  for(i=0 ; i < numNodes ; i++){
    r = q->next;
    chuck(&grbg, q);
    q = r;
  }
}
/* chucktree */

void copytree()
{
  /* Make a complete copy of the current tree for undo purposes */
  if (whichtree == 1)
    othertree = 0;
  else
    othertree = 1;
  
  if(treesets[othertree].root){
    chucktree(treesets[othertree].root);
  }
  
  treesets[othertree].root = copytrav(root);
  treesets[othertree].nonodes = nonodes;
  treesets[othertree].waswritten = waswritten;
  treesets[othertree].initialized = true;
}
/* copytree */



void getoptions()
{
  /* interactively set options */
  long loopcount;
  Char ch;
  boolean done, gotopt;

  how = arb;
  usertree = false;
  goteof = false;
  outgrno = 1;
  outgropt = false;
  thresh = false;
  threshold = spp;
  weights = false;
  ancvar = false;
  allsokal = false;
  allwagner = true;
  mixture = false;
  factors = false;
  loopcount = 0;
  do {
    printf((ansi || ibmpc) ? "\033[2J\033[H" : "\n");
    printf("\n\nInteractive mixed parsimony algorithm, version %s\n\n",
           VERSION);
    printf("Settings for this run:\n");
    printf("  X                         Use Mixed method?  %s\n",
           mixture ? "Yes" : "No");
    printf("  P                         Parsimony method?  %s\n",
           (allwagner && !mixture)   ? "Wagner"       :
           (!(allwagner || mixture)) ? "Camin-Sokal"  : "(methods in mixture)");
        
    printf("  A                     Use ancestral states?  %s\n",
           ancvar ? "Yes" : "No");
    printf("  F                  Use factors information?  %s\n",
           factors ? "Yes" : "No");
    printf("  O                            Outgroup root?  %s %3ld\n",
           outgropt ? "Yes, at species number" : "No, use as outgroup species",
           outgrno);
    printf("  W                           Sites weighted?  %s\n",
           (weights ? "Yes" : "No"));
    printf("  T                  Use Threshold parsimony?");
    if (thresh)
      printf("  Yes, count steps up to%4.1f\n", threshold);
    else
      printf("  No, use ordinary parsimony\n");
    printf("  U  Initial tree (arbitrary, user, specify)?  %s\n",
           (how == arb) ? "Arbitrary"                   :
           (how == use) ? "User tree from tree file"    : "Tree you specify");
    printf("  0       Graphics type (IBM PC, ANSI, none)?  %s\n",
           ibmpc ? "IBM PC" : ansi  ? "ANSI" : "(none)");
    printf("  S                 Width of terminal screen?");
    printf("%4ld\n", screenwidth);
    printf("  L                Number of lines on screen?%4ld",screenlines);
    printf("\n\nAre these settings correct?");
    printf(" (type Y or the letter for one to change)\n");
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
    gotopt = (strchr("SFOTXPAU0WL",ch) != NULL) ? true : false;
    if (gotopt) {
      switch (ch) {

      case 'F':
        factors = !factors;
        break;

      case 'X':
        mixture = !mixture;
        break;

      case 'W':
          weights = !weights;
          break;

      case 'P':
        allwagner = !allwagner;
        break;

      case 'A':
        ancvar = !ancvar;
        break;

      case 'O':
        outgropt = !outgropt;
        if (outgropt)
          initoutgroup(&outgrno, spp);
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

      case 'S':
        screenwidth= readlong("Width of terminal screen (in characters)?\n");
        break;

      case 'L':
        initnumlines(&screenlines);
        break;
      }
    }
    if (!(gotopt || done))
      printf("Not a possible option!\n");
    countup(&loopcount, 100);
  } while (!done);
  allsokal = (!allwagner && !mixture);
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
    if (factors) {
      factor = (Char *)Malloc(chars*sizeof(Char));
      inputfactors(chars, factor, &factors);
    }
    if (mixture)
        inputmixture(wagner0);
    if (weights)
      inputweights(chars, weight, &weights);
  putchar('\n');
  if (weights)
    printweights(stdout, 0, chars, weight, "Characters");
  for (i = 0; i < (words); i++) {
    if (mixture)
      wagner[i] = wagner0[i];
    else if (allsokal)
      wagner[i] = 0;
    else
      wagner[i] = (1L << (bits + 1)) - (1L << 1);
  }
  if (allsokal && !mixture)
    printf("Camin-Sokal parsimony method\n\n");
  if (allwagner && !mixture)
    printf("Wagner parsimony method\n\n");
  if (mixture)
    printmixture(stdout, wagner);
  for (i = 0; i < (chars); i++) {
    if (!ancvar) {
      anczero[i] = true;
      ancone[i] = (((1L << (i % bits + 1)) & wagner[i / bits]) != 0);
    } else {
      anczero[i] = anczero0[i];
      ancone[i] = ancone0[i];
    }
  }
  if (factors)
    printfactors(stdout, chars, factor, "");
  if (ancvar)
    printancestors(stdout, anczero, ancone);
  noroot = true;
  questions = false;
  for (i = 0; i < (chars); i++) {
    if (weight[i] > 0) {
      noroot = (noroot && ancone[i] && anczero[i] &&
          ((((1L << (i % bits + 1)) & wagner[i / bits]) != 0)
            || threshold <= 2.0));
    }
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
  numsteps = (steptr)Malloc(chars*sizeof(long));
  numsone = (steptr)Malloc(chars*sizeof(long));
  numszero = (steptr)Malloc(chars*sizeof(long));
  threshwt = (double *)Malloc(chars*sizeof(double));
  guess = (Char *)Malloc(chars*sizeof(Char));
  ancone = (boolean *)Malloc(chars*sizeof(boolean));
  anczero = (boolean *)Malloc(chars*sizeof(boolean));
  ancone0 = (boolean *)Malloc(chars*sizeof(boolean));
  anczero0 = (boolean *)Malloc(chars*sizeof(boolean));
  wagner = (bitptr)Malloc(words*sizeof(long));
  wagner0 = (bitptr)Malloc(words*sizeof(long));
  steps = (bitptr)Malloc(words*sizeof(long));
  zeroanc = (bitptr)Malloc(words*sizeof(long));
  oneanc = (bitptr)Malloc(words*sizeof(long));
  steps_a = (steptr)Malloc(chars*sizeof(long));
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
  if(mixture)
      openfile(&mixfile,MIXFILE,"mixture file", "r",progname,mixfilename);
  if(factors)
      openfile(&factfile,FACTFILE,"factors file", "r",progname,factfilename);

  alloctree(&treenode);
  setuptree(treenode);
  allocrest();
  inputoptions();
  inputdata(treenode, true, false, stdout);
}  /* doinput */


void configure()
{
  /* configure to machine -- set up special characters */
  chartype a;

  for (a = horiz; (long)a <= (long)question; a = (chartype)((long)a + 1))
    reversed[(long)a] = false;
  for (a = horiz; (long)a <= (long)question; a = (chartype)((long)a + 1))
    graphic[(long)a] = false;
  if (ibmpc) {
    che[(long)horiz] = 205;
    graphic[(long)horiz] = true;
    che[(long)vert] = 186;
    graphic[(long)vert] = true;
    che[(long)up] = 186;
    graphic[(long)up] = true;
    che[(long)overt] = 205;
    graphic[(long)overt] = true;
    che[(long)onne] = 219;
    reversed[(long)onne] = true;
    che[(long)zerro] = 176;
    graphic[(long)zerro] = true;
    che[(long)question] = 178;   /* or try CHR(177) */
    graphic[(long)question] = true;
    che[(long)upcorner] = 200;
    graphic[(long)upcorner] = true;
    che[(long)midcorner] = 204;
    graphic[(long)midcorner] = true;
    che[(long)downcorner] = 201;
    graphic[(long)downcorner] = true;
    return;
  }
  if (ansi) {
    che[(long)onne] = ' ';
    reversed[(long)onne] = true;
    che[(long)horiz] = che[(long)onne];
    reversed[(long)horiz] = true;
    che[(long)vert] = che[(long)onne];
    reversed[(long)vert] = true;
    che[(long)up] = 'x';
    graphic[(long)up] = true;
    che[(long)overt] = 'q';
    graphic[(long)overt] = true;
    che[(long)zerro] = 'a';
    graphic[(long)zerro] = true;
    reversed[(long)zerro] = true;
    che[(long)question] = '?';
    reversed[(long)question] = true;
    che[(long)upcorner] = 'm';
    graphic[(long)upcorner] = true;
    che[(long)midcorner] = 't';
    graphic[(long)midcorner] = true;
    che[(long)downcorner] = 'l';
    graphic[(long)downcorner] = true;
    return;
  }
  che[(long)horiz] = '=';
  che[(long)vert] = ' ';
  che[(long)up] = '!';
  che[(long)overt] = '-';
  che[(long)onne] = '*';
  che[(long)zerro] = '=';
  che[(long)question] = '.';
  che[(long)upcorner] = '`';
  che[(long)midcorner] = '+';
  che[(long)downcorner] = ',';
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
  putchar(che[(long)a]);
  postfix(a);
}  /* makechar */


void add_at(node *below, node *newtip, node *newfork)
{
  /* inserts the nodes newfork and its left descendant, newtip,
    to the tree.  below becomes newfork's right descendant */
  node *leftdesc, *rtdesc;

  if (below != treenode[below->index - 1])
    below = treenode[below->index - 1];

  if (newfork == NULL) {
    nonodes++;
    maketriad (&newfork, nonodes);
  }
  if (below->back != NULL) {
    below->back->back = newfork;
  }
  newfork->back = below->back;
  leftdesc = newtip;
  rtdesc = below;
  rtdesc->back = newfork->next->next;
  newfork->next->next->back = rtdesc;
  newfork->next->back = leftdesc;
  leftdesc->back = newfork->next;
  if (root == below)
    root = newfork;
  root->back = NULL;
}  /* add_at */


void add_before(node *atnode, node *newtip)
{
  /* inserts the node newtip together with its ancestral fork
     into the tree next to the node atnode. */
  node *q;

  if (atnode != treenode[atnode->index - 1])
    atnode = treenode[atnode->index - 1];
  q = treenode[newtip->index-1]->back;
  if (q != NULL) {
    q = treenode[q->index-1];
    if (newtip == q->next->next->back) {
      q->next->back = newtip;
      newtip->back = q->next;
      q->next->next->back = NULL;
    }
  }
  if (newtip->back != NULL) {
    add_at(atnode, newtip, treenode[newtip->back->index-1]);
  } else {
    add_at(atnode, newtip, NULL);
  }
}  /* add_before */


void add_child(node *parent, node *newchild)
{
  /* adds the node newchild into the tree as the last child of parent */

  int i;
  node *newnode, *q;

  if (parent != treenode[parent->index - 1])
    parent = treenode[parent->index - 1];
  gnu(&grbg, &newnode);
  newnode->tip = false;
  newnode->deleted=false;
  newnode->deadend=false;
  newnode->onebranch=false;
  newnode->onebranchhaslength=false;
  for (i=0;i<MAXNCH;i++)
    newnode->nayme[i] = '\0';
  newnode->index = parent->index;
  if(!newnode->statezero)
    newnode->statezero = (bitptr)Malloc(words*sizeof(long));
  if(!newnode->stateone)
    newnode->stateone = (bitptr)Malloc(words*sizeof(long));
  if(!newnode->zerosteps)
    newnode->zerosteps = (bitptr)Malloc(words*sizeof(long));
  if(!newnode->onesteps)
    newnode->onesteps = (bitptr)Malloc(words*sizeof(long));
  q = parent;
  do {
    q = q->next;
  } while (q->next != parent);
  newnode->next = parent;
  q->next = newnode;
  newnode->back = newchild;
  newchild->back = newnode;
  if (newchild->haslength) {
    newnode->length = newchild->length;
    newnode->haslength = true;
  } else
    newnode->haslength = false;
}  /* add_child */


void dnamove_add(node *below, node *newtip, node *newfork)
{
  /* inserts the nodes newfork and its left descendant, newtip,
     to the tree.  below becomes newfork's right descendant */
  boolean putleft;
  node *leftdesc, *rtdesc;

  if (below != treenode[below->index - 1])
    below = treenode[below->index - 1];
  if (below->back != NULL)
    below->back->back = newfork;
  newfork->back = below->back;
  putleft = true;
  if (restoring)
    putleft = wasleft;
  if (putleft) {
    leftdesc = newtip;
    rtdesc = below;
  } else {
    leftdesc = below;
    rtdesc = newtip;
  }
  rtdesc->back = newfork->next->next;
  newfork->next->next->back = rtdesc;
  newfork->next->back = leftdesc;
  leftdesc->back = newfork->next;
  if (root == below)
    root = newfork;
  root->back = NULL;
  newfork->numdesc = 2;
}  /* dnamove_add */


void move_re_move(node **item, node **fork)
{
/* Removes node item from the tree.  If item has one sibling,
   removes its ancestor, fork, from the tree as well and attach
   item's sib to fork's ancestor.  In this case, it returns a pointer
   to the removed fork node which is still attached to item.
*/
  node *p=NULL, *q;
  int nodecount;

  if ((*item)->back == NULL) {
    *fork = NULL;
    return;
  }
  *fork = treenode[(*item)->back->index - 1];
  nodecount = 0;
  if ((*fork)->next->back == *item)
    p = *fork;
  q = (*fork)->next;
  do {
    nodecount++;
    if (q->next->back == *item)
      p = q;
    q = q->next;
  } while (*fork != q);

  if (nodecount > 2)
  {
    fromtype = atnode;
    p->next = (*item)->back->next;
    chuck(&grbg, (*item)->back);
    (*item)->back = NULL;
    *fork = NULL;
  } else {
    /* traditional (binary tree) remove code */
    if (*item == (*fork)->next->back) {
      if (root == *fork)
     root = (*fork)->next->next->back;
    } else {
      if (root == *fork)
     root = (*fork)->next->back;
    }
    fromtype = beforenode;
    /* stitch nodes together, leaving out item */
    p = (*item)->back->next->back;
    q = (*item)->back->next->next->back;
    if (p != NULL)
      p->back = q;
    if (q != NULL)
      q->back = p;
    if (haslengths) {
      if (p != NULL && q != NULL) {
     p->length += q->length;
     q->length = p->length;
      } else
     (*item)->length = (*fork)->next->length + (*fork)->next->next->length;
    }
    (*fork)->back = NULL;
    p = (*fork)->next;
    while (p != *fork) {
      p->back = NULL;
      p = p->next;
    }
    (*item)->back = NULL;
  } /* endif nodecount > 2 else */
}  /* dnamove_re_move */


void new_move_fillin3(node *p)
{
  node *q;
  long n = 0;
  long n0, n1;
  long i;
  
  /* Establish number of children of this node */
  q = p->next;
  do {
    q = q->next;
    n++;
  } while (q != p);

  for (i = 0; i < (chars); i++) {
    n0 = 0;
    n1 = 0;
    q = p->next;
    do {
      if (q->next->back->statezero[i]) {
        n0++;
      }
      if (q->next->back->stateone[i]) {
        n1++;
      }
      q = q->next;
    } while (q != p);
    if (n1 > n0) {
      p->stateone[i] = 1;
    } else {
      if (n0 > n1) {
        p->statezero[i] = 1;
      } else { /* n0, n1 equal */
        p->statezero[i] = 1;
        p->stateone[i] = 1;
      }
    }

    if ((n - n1) <= (n - n0)) {
      steps[i] = (n - n1);
    } else {
      steps[i] = (n - n0);
    }
  }
  
} /* new_move_fillin */

void new_move_fillin(node *p)
{
  node *q;
  long n = 0;
  long n0, n1;
  long nzerosteps, nonesteps;
  long i, j, k;
  long wa, za, oa;
  boolean wa_2[chars];
  boolean za_2[chars];
  boolean oa_2[chars];
  
  for (i=0;i<chars;i++) {
    j = i / bits + 1;
    k = i % bits + 1;
    printf("i = %d, j = %d, k = %d\n",i,j,k);
    
    wa_2[i] = ((long)wagner[j - 1]) & (1L << k);
    za_2[i] = ((long)zeroanc[j - 1]) & (1L << k);
    oa_2[i] = ((long)oneanc[j - 1]) & (1L << k);
  }
  /**/

  /* Establish number of children of this node */
  q = p->next;
  do {
    q = q->next;
    n++;
  } while (q != p);

  for (i = 0; i < (chars); i++) {
    j = i / bits + 1;
    k = i % bits + 1;
    wa = wa_2[i];
    za = za_2[i];
    oa = oa_2[i];

    if (wa) {
      /*xx*/  printf("wagner true, i = %d\n",i);
            
      //------------- wagner[i] if
      /*p->statezero[i] = 0;*/
      p->statezero[j] &= ~(1L << k); /*xx*/
      
      /*p->stateone[i] = 0;*/
      p->stateone[j] &= ~(1L << k);
      
      n0 = 0;
      n1 = 0;
      q = p->next;
        
      do {
#if 0
        if (q->back->statezero[i]) {
          n0++;
        }
        if (q->back->stateone[i]) {
          n1++;
        }
#endif
        if ((q->back->statezero[j]) & (1L << k)) {
          n0++;
        }
        if ((q->back->stateone[j]) & (1L << k)) {
          n1++;
        }
        q = q->next;
      } while (q != p);
      if (n1 > n0) {
        /*p->stateone[i] = 1;*/
        p->stateone[j] |= (1L << k);
      } else {
        if (n0 > n1) {
          /*p->statezero[i] = 1;*/
          p->statezero[j] |= (1L << k);
        } else { /* n0, n1 equal */
          /*p->statezero[i] = 1;
            p->stateone[i] = 1;*/
          p->statezero[j] |= (1L << k);
          p->stateone[j] |= (1L << k);
        }
      }

      /*xx work here */      
      if ((n - n1) <= (n - n0)) {
        steps[i] = (n - n1);
      } else {
        steps[i] = (n - n0);
      }
      /*xx*/printf("new_move_fillin: steps[%d] = %d\n",i,steps[i]);
    } else {
      
      /*xx*/  printf("wagner false, i = %d\n",i);
    
    //---------------- wagner[i] else

      nzerosteps = ((n1 > 0) ? n1 : 0);
      nonesteps = ((n0 > 0) ? n0 : 0);
      if (za)
        p->zerosteps[i] = nzerosteps;

      if (oa)
        p->onesteps[i] = nonesteps;

      if (za && oa) {
        p->zerosteps[i] = ((nzerosteps < nonesteps) ? nzerosteps : nonesteps);
        p->onesteps[i] = p->zerosteps[i];
      }
    }

    //-----------------
  }
  puts("new_move_fillin ending ----");
  
} /* new_move_fillin */


void move_fillin(node *p)
{
  /* Sets up for each node in the tree two statesets.
     stateone and statezero are the sets of character
     states that must be 1 or must be 0, respectively,
     in a most parsimonious reconstruction, based on the
     information at or above this node.   Note that this
     state assignment may change based on information further
     down the tree.   If a character is in both sets it is in
     state "P".   If in neither, it is "?". */
  long i;
  long l0, l1, r0, r1, st, wa, za, oa;

  for (i = 0; i < (words); i++) {
    l0 = p->next->back->statezero[i];
    l1 = p->next->back->stateone[i];
    r0 = p->next->next->back->statezero[i];
    r1 = p->next->next->back->stateone[i];
    wa = wagner[i];
    za = zeroanc[i];
    oa = oneanc[i];
    st = (l1 & r0) | (l0 & r1);
    steps[i] = st;
    printf("move_fillin: steps[%d] = %d\n",i,steps[i]);/*xx*/
    p->stateone[i] = (l1 | r1) & (~(st & (wa | za)));
    p->statezero[i] = (l0 | r0) & (~(st & (wa | oa)));
  }
}  /* move_fillin */



void move_postorder(node *p)
{
  /* traverses a binary tree, calling function fillin at a
     node's descendants before calling fillin at the node */
  node *q;
  
  if (p->tip)
    return;

  q = p->next;
  do {
    move_postorder(q->back);
    q = q->next;
  } while (p != q);

  /*  move_fillin(p);*/
      new_move_fillin(p);
  
  count(steps, zeroanc, numszero, numsone);
}  /* move_postorder */


void evaluate(node *r)
{
  /* Determines the number of steps needed for a tree.
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
  move_postorder(r);
  count(r->stateone, zeroanc, numszero, numsone);
  for (i = 0; i < (words); i++) {
    zeroanc[i] = 0;
    oneanc[i] = fullset;
  }
  move_postorder(r);
  count(r->statezero, zeroanc, numszero, numsone);
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


void move_reroot(node *outgroup)
{
  /* Reorient tree so that outgroup is by itself on the left of the root */ 
  node *p, *q, *r;
  long nodecount = 0;
  double templen;

  if(outgroup->back->index == root->index)
    return;

  q = root->next;
  do {                    /* when this loop exits, p points to the internal */
    p = q;                /* node to the right of root */
    nodecount++;
    q = p->next;
  } while (q != root);
  r = p;

  /* reorient nodep array

     The nodep array must point to the ring member of each ring
     that is closest to the root.  The while loop changes the ring member
     pointed to by treenode[] for those nodes that will have their
     orientation changed by the reroot operation.
  */
  p = outgroup->back;
  while (p->index != root->index) {
    q = treenode[p->index - 1]->back;
    treenode[p->index - 1] = p;
    p = q;
  }
  if (nodecount > 2)
    treenode[p->index - 1] = p;

  /* If nodecount > 2, the current node ring to which root is pointing
     will remain in place and root will point somewhere else. */
  /* detach root from old location */
  if (nodecount > 2) {
    r->next = root->next;
    root->next = NULL;
    nonodes++;
    maketriad(&root, nonodes);

    if (haslengths) {
      /* root->haslength remains false, or else treeout() will generate
         a bogus extra length */
      root->next->haslength = true;
      root->next->next->haslength = true;
    }
  } else { /* if (nodecount > 2) else */
    q = root->next;
    q->back->back = r->back;
    r->back->back = q->back;

    if (haslengths) {
      r->back->length = r->back->length + q->back->length;
      q->back->length = r->back->length;
    }
  } /* if (nodecount > 2) endif */

  /* tie root into new location */
  root->next->back = outgroup;
  root->next->next->back = outgroup->back;
  outgroup->back->back = root->next->next;
  outgroup->back = root->next;

  /* place root equidistant between left child (outgroup) and
     right child by deviding outgroup's length */
  if (haslengths) {
    templen = outgroup->length / 2.0;
    outgroup->length = templen;
    outgroup->back->length = templen;
    root->next->next->length = templen;
    root->next->next->back->length = templen;
  }
} /* move_reroot */


void old_reroot(node *outgroup)
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
}  /* old_reroot */


void move_filltrav(node *r)
{
  /* traverse to fill in interior node states */
  node *q;
  
  if (r->tip)
    return;
  q = r->next;
  do {
    move_filltrav(q->back);
    q = q->next;
  } while ( q != r);
  
  move_fillin(r);
}  /* move_filltrav */


void move_hyptrav(node *r)
{
  /* compute states at one interior node */
  long i;
  boolean bottom;
  long l0, l1, r0, r1, s0, s1, a0, a1, temp, wa;
  gbit *zerobelow = NULL, *onebelow = NULL;

  disc_gnu(&zerobelow, &garbage);
  disc_gnu(&onebelow, &garbage);
  bottom = (r->back == NULL);
  if (bottom) {
    memcpy(zerobelow->bits_, zeroanc, words*sizeof(long));
    memcpy(onebelow->bits_, oneanc, words*sizeof(long));
  } else {
    memcpy(zerobelow->bits_, treenode[r->back->index - 1]->statezero,
           words*sizeof(long));
    memcpy(onebelow->bits_, treenode[r->back->index - 1]->stateone,
           words*sizeof(long));
  }
  for (i = 0; i < (words); i++) {
    s0 = r->statezero[i];
    s1 = r->stateone[i];
    a0 = zerobelow->bits_[i];
    a1 = onebelow->bits_[i];
    if (!r->tip) {
      wa = wagner[i];
      l0 = r->next->back->statezero[i];
      l1 = r->next->back->stateone[i];
      r0 = r->next->next->back->statezero[i];
      r1 = r->next->next->back->stateone[i];
      s0 = (wa & ((a0 & l0) | (a0 & r0) | (l0 & r0))) |
           (fullset & (~wa) & s0);
      s1 = (wa & ((a1 & l1) | (a1 & r1) | (l1 & r1))) |
           (fullset & (~wa) & s1);
      temp = fullset & (~(s0 | s1 | l1 | l0 | r1 | r0));
      s0 |= temp & a0;
      s1 |= temp & a1;
      r->statezero[i] = s0;
      r->stateone[i] = s1;
    }
  }
  if (((1L << dispbit) & r->stateone[dispword - 1]) != 0) {
    if (((1L << dispbit) & r->statezero[dispword - 1]) != 0)
      r->state = '?';
    else
      r->state = '1';
  } else {
    if (((1L << dispbit) & r->statezero[dispword - 1]) != 0)
      r->state = '0';
    else
      r->state = '?';
  }
  if (!r->tip) {
    move_hyptrav(r->next->back);
    move_hyptrav(r->next->next->back);
  }
  disc_chuck(zerobelow, &garbage);
  disc_chuck(onebelow, &garbage);
}  /* move_hyptrav */


void move_hypstates()
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
  move_filltrav(root);
  move_hyptrav(root);
}  /* move_hypstates */


void grwrite(chartype c, long num, long *pos)
{
  long i;

  prefix(c);
  for (i = 1; i <= num; i++) {
    if ((*pos) >= leftedge && (*pos) - leftedge + 1 < screenwidth)
      putchar(che[(long)c]);
    (*pos)++;
  }
  postfix(c);
}  /* grwrite */


void move_drawline(long i, long lastline)
{
  /* draws one row of the tree diagram by moving up tree */
  node *p, *q, *r, *first =NULL, *last =NULL;
  long n, j, pos;
  boolean extra, done;
  Char st;
  chartype c, d;

  pos = 1;
  p = nuroot;
  q = nuroot;
  extra = false;
  if (i == (long)p->ycoord && (p == root || subtree)) {
    extra = true;
    c = overt;
    if (display) {
      switch (p->state) {

      case '1':
        c = onne;
        break;

      case '0':
        c = zerro;
        break;

      case '?':
        c = question;
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
      c = overt;
      if (q == first)
        d = downcorner;
      
      else if (q == last)
        d = upcorner;
      
      else if ((long)q->ycoord == (long)p->ycoord)
        d = c;
      
      else
        d = midcorner;
      


      if (display) {
        switch (q->state) {

        case '1':
          c = onne;
          break;

        case '0':
          c = zerro;
          break;

        case '?':
          c = question;
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
          st = p->next->back->state;
        else
          st = p->next->next->back->state;
        if (display) {
          switch (st) {

          case '1':
            c = onne;
            break;

          case '0':
            c = zerro;
            break;

          case '?':
            c = question;
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
}  /* move_drawline */


void move_printree()
{
  /* prints out diagram of the tree */
  long tipy, i, dow;

  if (!subtree)
    nuroot = root;
  if (changed || newtree)
    evaluate(root);
  if (display)
    move_hypstates();
  if (ansi || ibmpc)
    printf("\033[2J\033[H");
  else
    putchar('\n');
  tipy = 1;
  dow = down;
  if (spp * dow > screenlines && !subtree) {
    dow--;
  }
  if (noroot)
    printf("(unrooted)");
  if (display) {
    printf(" ");
    makechar(onne);
    printf(":1 ");
    makechar(question);
    printf(":? ");
    makechar(zerro);
    printf(":0 ");
  } else
    printf("             ");
  if (!earlytree) {
    printf("%10.1f Steps", -like);
  }
  if (display)
    printf("  CHAR%3ld", dispchar);
  else
    printf("         ");
  if (!earlytree) {
    printf("  %3ld chars compatible\n", compatible);
  }
  printf("                            ");
  if (changed && !earlytree) {
    if (-like < bestyet) {
      printf("     BEST YET!");
      bestyet = -like;
    } else if (fabs(-like - bestyet) < 0.000001)
      printf("  (as good as best)");
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
  vmargin = 4;
  treelines = tipy - dow;
  if (topedge != 1) {
    printf("** %ld lines above screen **\n", topedge - 1);
    vmargin++;
  }

  if ((treelines - topedge + 1) > (screenlines - vmargin))
    vmargin++;
  for (i = 1; i <= treelines; i++) {
    if (i >= topedge && i < topedge + screenlines - vmargin)
      move_drawline(i, treelines);
  }

  if ((treelines - topedge + 1) > (screenlines - vmargin)) {
    printf("** %ld", treelines - (topedge - 1 + screenlines - vmargin));
    printf(" lines below screen **\n");
  }
  if (treelines - topedge + vmargin + 1 < screenlines)
    putchar('\n');
  gotlike = -like;
  changed = false;
}  /* move_printree */


void arbitree()
{
  long i;

  root = treenode[0];
  add2(treenode[0], treenode[1], treenode[spp], &root,
    restoring, wasleft, treenode);
  for (i = 3; i <= (spp); i++)
    add2(treenode[spp+ i - 3], treenode[i - 1], treenode[spp + i - 2],
      &root, restoring, wasleft, treenode);
  for (i = 0; i < (nonodes); i++)
    in_tree[i] = true;
}  /* arbitree */


void yourtree()
{
  long i, j;
  boolean ok;

  root = treenode[0];
  dnamove_add(treenode[0], treenode[1], treenode[spp]);
  i = 2;
  do {
    i++;
    move_printree();
    printf("Add species%3ld: ", i);
    for (j = 0; j < nmlngth; j++)
      putchar(nayme[i - 1][j]);
    do {
      printf("\n at or before which node (type number): ");
      inpnum(&j, &ok);
      ok = (ok && ((j >= 1 && j < i) || (j > spp && j < spp + i - 1)));
      if (!ok)
        printf("Impossible number. Please try again:\n");
    } while (!ok);

    if (j >= i) {   /* has user chosen a non-tip? if so, offer choice */
      do {
        printf(" Insert at node (A) or before node (B)? ");
#ifdef WIN32
        phyFillScreenColor();
#endif
        fflush(stdout);
        scanf("%c%*[^\n]", &ch);
        getchar();
        if (ch == '\n')
          ch = ' ';
        ch = isupper(ch) ? ch : toupper(ch);
      } while (ch != 'A' && ch != 'B');
    }
    else ch = 'B';   /* if user has chosen a tip, set Before */

    if (j != 0) {
      if (ch == 'A') {
        if (!treenode[j - 1]->tip) {
          add_child(treenode[j - 1], treenode[i - 1]);
        }
      } else {
        printf("dnamove_add(below %ld, newtip %ld, newfork %ld)\n",j-1,i-1,spp+i-2);
        dnamove_add(treenode[j - 1], treenode[i - 1], treenode[spp + i - 2]);
      } /* endif (before or at node) */
    }


  } while (i != spp);
}  /* yourtree */


void initmovenode(node **p, node **grbg, node *q, long len, long nodei,
                        long *ntips, long *parens, initops whichinit,
                        pointarray treenode, pointarray nodep, Char *str, Char *ch,
                        FILE *intree)
{
  /* initializes a node */
  /* LM 7/27  I added this function and the commented lines around */
  /* treeread() to get the program running, but all 4 move programs*/
  /* are improperly integrated into the v4.0 support files.  As is */
  /* this is a patchwork function              */
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
    break;
  }
} /* initmovenode */


void buildtree()
{
  long i, j, nextnode;
  node *p;

  treeone = (node **)Malloc(maxsz*sizeof(node *));
  treetwo = (node **)Malloc(maxsz*sizeof(node *));
  treesets[othertree].treenode = treetwo;
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
    firsttree = true;                                            /**/
    nodep = NULL;                                       /**/
    nextnode = 0;                                           /**/
    haslengths = 0;                                         /**/
    zeros = (long *)Malloc(chars*sizeof(long));         /**/
    for (i = 0; i < chars; i++)                           /**/
      zeros[i] = 0;                                         /**/
    treesets[whichtree].nodep = nodep;
    treeread(intree, &root, treenode, &goteof, &firsttree,
                nodep, &nextnode, &haslengths,
                &grbg, initmovenode); /*debug*/
    for (i = spp; i < (nonodes); i++) {
      p = treenode[i];
      for (j = 1; j <= 3; j++) {
        p->stateone = (bitptr)Malloc(words*sizeof(long));
        p->statezero = (bitptr)Malloc(words*sizeof(long));
        p->onesteps = (bitptr)Malloc(words*sizeof(long));
        p->zerosteps = (bitptr)Malloc(words*sizeof(long));
        p = p->next;
      }
    } /* debug: see comment at initmovenode() */

   /*  treeread(which, ch, &root, treenode, names);*/
    for (i = 0; i < (spp); i++)
      in_tree[i] = names[i];
    free(names);
    FClose(intree);
    break;

  case spec:
    yourtree();
    break;
  }
  if (!outgropt)
    outgrno = root->next->back->index;
  if (outgropt && in_tree[outgrno - 1])
    move_reroot(treenode[outgrno - 1]);
}  /* buildtree */


void consolidatetree(long index)
{
  node *start, *r, *q;
  int i;

  i = 0;

  start = treenode[index - 1];
  q = start->next;

  while (q != start) {
    r = q;
    q = q->next;
    chuck(&grbg, r);
  } 
  chuck(&grbg, q);
  
  i = index;

  while (i <= nonodes) {
 
    r = treenode[i - 1];
    if (!(r->tip))
       r->index--;
    if (!(r->tip)) {
      q = r->next;
      do {
        q->index--;
        q = q->next;
      } while (r != q && q != NULL);
    }
    treenode[i - 1] = treenode[i];
    i++;
  }

  nonodes--;
} /* consolidatetree */


void rearrange()
{
  long i, j, maxinput;
  boolean ok1, ok2;
  node *p, *q;
  char ch;

  printf("Remove everything to the right of which node? ");
  inpnum(&i, &ok1);
  ok1 = (ok1 && i >= 1 && i <= (spp * 2 - 1) && i != root->index);
  if (ok1)
    ok1 = !treenode[i - 1]->deleted;
  if (ok1) {
    printf("Add at or before which node? ");
    inpnum(&j, &ok2);
    ok2 = (ok2 && j >= 1 && j <= (spp * 2 - 1));
    if (ok2) {
      if (j != root->index)
        ok2 = !treenode[treenode[j - 1]->back->index - 1]->deleted;
    }
    if (ok2) {
/*xx This edit says "j must not be i's parent."
     Is this necessary anymore?   */
    /*  ok2 = (nodep[j - 1] != nodep[nodep[i - 1]->back->index - 1]);*/
      p = treenode[j - 1];
      /* make sure that j is not a descendent of i */
      while (p != root) {
        ok2 = (ok2 && p != treenode[i - 1]);
        p = treenode[p->back->index - 1];
      }
      if (ok1 && ok2) {
        maxinput = 1;
        do {
          printf("Insert at node (A) or before node (B)? ");
#ifdef WIN32
          phyFillScreenColor();
#endif
          fflush(stdout);
          scanf("%c%*[^\n]", &ch);
          getchar();
          if (ch == '\n')
            ch = ' ';
          ch = isupper(ch) ? ch : toupper(ch);
          maxinput++;
          if (maxinput == 100) {
            printf("ERROR: too many tries at choosing option\n");
            exxit(-1);
          }
        } while (ch != 'A' && ch != 'B');
        if (ch == 'A') {
          if (!(treenode[j - 1]->deleted) && !treenode[j - 1]->tip) {
            changed = true;
            copytree();
            move_re_move(&treenode[i - 1], &q);
            add_child(treenode[j - 1], treenode[i - 1]);
            if (fromtype == beforenode)
              consolidatetree(q->index);
          } else
            ok2 = false;
        } else {
          if (j != root->index) { /* can't insert at root */
            changed = true;
            copytree();
            move_re_move(&treenode[i - 1], &q);
            if (q != NULL) {
              treenode[q->index-1]->next->back = treenode[i-1];
              treenode[i-1]->back = treenode[q->index-1]->next;
            }
            add_before(treenode[j - 1], treenode[i - 1]);
          } else
            ok2 = false;
        } /* endif (before or at node) */
      } /* endif (ok to do move) */
    } /* endif (destination node ok) */
  } /* endif (from node ok) */
  move_printree();
  if (!(ok1 && ok2))
    printf("Not a possible rearrangement.  Try again: \n");
  else {
    written = false;
  }
}  /* rearrange */


void old_rearrange()
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
        changed = true;
        copytree();
        what = i;
        q = treenode[treenode[i - 1]->back->index - 1];
        if (q->next->back->index == i)
          fromwhere = q->next->next->back->index;
        else
          fromwhere = q->next->back->index;
        towhere = j;
        re_move2(&treenode[i - 1], &q, &root, &wasleft, treenode);
        add2(treenode[j - 1], treenode[i - 1], q, &root,
          restoring, wasleft, treenode);
      }
      lastop = rearr;
    }
  }
  move_printree();
  if (!(ok1 && ok2))
    printf("Not a possible rearrangement.   Try again: ");
  else {
    oldwritten =written;
    written = false;
  }
}  /* old_rearrange */


void tryadd(node *p, node *item, node *nufork, double *place)
{
  /* temporarily adds one fork and one tip to the tree.
     Records scores in array place */
  add2(p, item, nufork, &root, restoring, wasleft, treenode);
  evaluate(root);
  place[p->index - 1] = -like;
  re_move2(&item, &nufork, &root, &wasleft, treenode);
}  /* tryadd */


void addpreorder(node *p, node *item, node *nufork, double *place)
{
  /* traverses a binary tree, calling function tryadd
     at a node before calling tryadd at its descendants */
  if (p == NULL)
    return;
  tryadd(p,item,nufork,place);
  if (!p->tip) {
    addpreorder(p->next->back, item, nufork, place);
    addpreorder(p->next->next->back, item, nufork, place);
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
  copytree();
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
  wasleft =oldleft;
  restoring = true;
  add2(treenode[fromwhere - 1], treenode[what - 1],dummy, &root,
    restoring, wasleft, treenode);
  like = -current;
  compatible = oldcompat;
  restoring = false;
  better = false;
  printf("          BETTER: ");
  for (j = 1; j <= (nonodes); j++) {
    if (place[j - 1] < current && place[j - 1] >= 0.0) {
      printf("%3ld:%6.2f", j, place[j - 1]);
      better = true;
    }
  }
  if (!better)
    printf(" NONE");
  printf("\n          TIED:    ");
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
  boolean btemp;

  /* don't undo to an uninitialized tree */
  if (!treesets[othertree].initialized) {
    move_printree();
    printf("Nothing to undo.\n");
    return;
  }

  treesets[whichtree].root = root;
  treesets[whichtree].treenode = treenode;
  treesets[whichtree].nonodes = nonodes;
  treesets[whichtree].waswritten = waswritten;
  treesets[whichtree].initialized = true;

  whichtree = othertree;

  root = treesets[whichtree].root;
  treenode = treesets[whichtree].treenode;
  nonodes = treesets[whichtree].nonodes;
  waswritten = treesets[whichtree].waswritten;

  if (othertree == 0)
    othertree = 1;
  else
    othertree = 0;

  changed = true;
  move_printree();
  btemp = oldwritten;
  oldwritten = written;
  written = btemp;
}  /* undo */

#if 0
void old_undo()
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
    add2(treenode[fromwhere - 1], treenode[what - 1],q, &root,
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
  move_printree();
  if (lastop == none) {
    printf("No operation to undo! \n");
    return;
  }
  btemp = oldwritten;
  oldwritten = written;
  written = btemp;
}  /* old_undo */
#endif

void treewrite(boolean done)
{
  /* write out tree to a file */
  Char ch;

  treeoptions(waswritten, &ch, &outtree, outtreename, progname);
  if (!done)
    move_printree();
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
}
  /* treewrite */


void clade()
{
  /* pick a subtree and show only that on screen */
  long i;
  boolean ok;

  printf("Select subtree rooted at which node (0 for whole tree)? ");
  inpnum(&i, &ok);
  ok = (ok && (unsigned)i <= nonodes);
  if (ok) {
    subtree = (i > 0);
    if (subtree)
      nuroot = treenode[i - 1];
    else
      nuroot = root;
  }
  move_printree();
  if (!ok)
    printf("Not possible to use this node. ");
}  /* clade */


void oldflip()
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
  move_printree();
  if (ok) {
    oldwritten = written;
    written = false;
    return;
  }
  if (i >= 1 && i <= spp)
    printf("Can't flip there. ");
  else
    printf("No such node. ");
}  /* oldflip */


void fliptrav(node *p, boolean recurse)
{
  node *q, *temp, *r =NULL, *rprev =NULL, *l, *lprev;
  boolean lprevflag;
  int nodecount, loopcount, i;

  if (p->tip)
    return;

  q = p->next;
  l = q;
  lprev = p;
  nodecount = 0;

  do {
    nodecount++;
    if (q->next->next == p) {
      rprev = q;
      r = q->next;
    }
    q = q->next;
  } while (p != q);

  if (nodecount == 1)
    return;
  loopcount = nodecount / 2;

  for (i=0; i<loopcount; i++) {
    lprev->next = r;
    rprev->next = l;
    temp = r->next;
    r->next = l->next;
    l->next = temp;
    if (i < (loopcount - 1)) {
      lprevflag = false;
      q = p->next;
      do {
        if (q == lprev->next && !lprevflag) {
          lprev = q;
          l = q->next;
          lprevflag = true;
        }
        if (q->next == rprev) {
          rprev = q;
          r = q->next;
        }
        q = q->next;
      } while (p != q);
    }
  }
  if (recurse) {
    q = p->next;
    do {
      fliptrav(q->back, true);
      q = q->next;
    } while (p != q);
  }
}  /* fliptrav */


void flip(long atnode)
{
  /* flip at a node left-right */
  long i;
  boolean ok;

  if (atnode == 0) {
    printf("Flip branches at which node? ");
    inpnum(&i, &ok);
    ok = (ok && i > spp && i <= nonodes);
  } else {
    i = atnode;
    ok = true;
  }
  if (ok) {
    copytree();
    fliptrav(treenode[i - 1], true);
  }
  if (atnode == 0)
    move_printree();
  if (ok) {
    written = false;
    return;
  }
  if ((i >= 1 && i <= spp) ||
      (i > spp && i <= nonodes))
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

/*xx may have to remove intree[] references */
  if (in_tree[outgrno - 1])
    move_reroot(treenode[outgrno - 1]);
  changed = true;
  lastop = reroott;
  move_printree();
  oldwritten = written;
  written = false;
}  /* changeoutgroup */


void redisplay()
{
  boolean done=false;
  waswritten = false;
  do {
    printf("\nNEXT? (Options: R # + - S . T U W O F H J K L C ? X Q) ");
    printf("(? for Help) ");
#ifdef WIN32
    phyFillScreenColor();
#endif
    fflush(stdout);
    scanf("%c%*[^\n]", &ch);
    getchar();
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
    if (strchr("R#+-S.TUWOFHJKLC?XQ",ch) != NULL){
      switch (ch) {

      case 'R':
        rearrange();
        break;

      case '#':
        nextinc(&dispchar, &dispword, &dispbit, chars, bits, &display,
                  numsteps, weight);
        move_printree();
        break;

      case '+':
        nextchar(&dispchar, &dispword, &dispbit, chars, bits, &display);
        move_printree();
        break;

      case '-':
        prevchar(&dispchar, &dispword, &dispbit, chars, bits, &display);
        move_printree();
        break;

      case 'S':
        show(&dispchar, &dispword, &dispbit, chars, bits, &display);
        move_printree();
        break;

      case '.':
        move_printree();
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
        flip(0);
        break;

      case 'H':
        window(left, &leftedge, &topedge, hscroll, vscroll, treelines,
                 screenlines, screenwidth, farthest, subtree);
        move_printree();
        break;

      case 'J':
        window(downn, &leftedge, &topedge, hscroll, vscroll, treelines,
                 screenlines, screenwidth, farthest, subtree);
        move_printree();
        break;

      case 'K':
        window(upp, &leftedge, &topedge, hscroll, vscroll, treelines,
                 screenlines, screenwidth, farthest, subtree);
        move_printree();
        break;

      case 'L':
        window(right, &leftedge, &topedge, hscroll, vscroll, treelines,
                 screenlines, screenwidth, farthest, subtree);
        move_printree();
        break;

      case 'C':
        clade();
        break;

      case '?':
        help("character");
        move_printree();
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
      fflush(stdout);
      scanf("%c%*[^\n]", &ch);
      getchar();
      if (ch == '\n')
        ch = ' ';
    } while (ch != 'Y' && ch != 'y' && ch != 'N' && ch != 'n');
  }
  if (ch == 'Y' || ch == 'y')
    treewrite(done);
}  /* redisplay */


void treeconstruct()
{
  /* constructs a binary tree from the pointers in treenode. */

  int i;
  restoring = false;
  subtree = false;
  display = false;
  dispchar = 0;
  fullset = (1L << (bits + 1)) - (1L << 1);
  earlytree = true;
  buildtree();

  /* get an accurate value for nonodes by finding out where the nodes
     really stop
  */
  for (i=0;i<nonodes;i++) {
    
    if (treenode[i]==NULL)
      break;
    
  }
  nonodes = i;

  waswritten = false;
  printf("\nComputing steps needed for compatibility in characters...\n\n");
  newtree = true;
  earlytree = false;
  move_printree();
  bestyet = -like;
  gotlike = -like;
  lastop = none;
  newtree = false;
  written = false;
  redisplay();
}  /* treeconstruct */


int main(int argc, Char *argv[])
{  /* Interactive mixed parsimony                                      */
   /* reads in spp, chars, and the data. Then calls treeconstruct to   */
   /*  construct the tree and query the user                           */
#ifdef MAC
  argc = 1;                /* macsetup("Move","");                */
  argv[0] = "Move";
#endif
  init(argc, argv);
  progname = argv[0];
  strcpy(infilename,INFILE);
  strcpy(intreename,INTREE);
  strcpy(outtreename,OUTTREE);
  openfile(&infile,infilename,"input file", "r",argv[0],infilename);

  whichtree = 0;
  othertree = 1;
  treesets[whichtree].initialized = false;
  treesets[othertree].initialized = false;
  screenlines = 24;
  scrollinc = 20;
  screenwidth = 80;
  topedge = 1;
  leftedge = 1;
  ibmpc = IBMCRT;
  ansi = ANSICRT;
  root = NULL;
  bits = 8*sizeof(long) - 1;
  doinput();
  configure();
  treeconstruct();
  FClose(outtree);
#ifdef MAC
  fixmacfile(outtreename);
#endif
#ifdef WIN32
  phyRestoreConsoleAttributes();
#endif
  return 0;
}  /* Interactive mixed parsimony */

