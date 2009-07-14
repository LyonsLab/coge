
#include "phylip.h"
#include "moves.h"
#include "seq.h"

/* version 3.6. (c) Copyright 1993-2004 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#define maxsz           999   /* size of pointer array for the undo trees */
                              /* this can be large without eating memory */

typedef struct treeset_t {
  node *root;
  pointarray nodep;
  pointarray treenode;
  long nonodes;
  boolean waswritten, hasmult, haslengths, nolengths, initialized;
} treeset_t;

treeset_t treesets[2];

node **treeone, **treetwo;

typedef enum {
  horiz, vert, up, overt, upcorner, midcorner, downcorner, aa, cc, gg, tt, question
  } chartype;

typedef enum {
  rearr, flipp, reroott, none
  } rearrtype;

typedef struct gbase2 {
  baseptr2 base2;
  struct gbase2 *next;
} gbase2;

typedef enum {
  arb, use, spec
  } howtree;

typedef enum {beforenode, atnode} movet;

movet fromtype;

typedef node **pointptr;

#ifndef OLDC
/* function prototypes */
void dnamove_gnu(gbases **);
void dnamove_chuck(gbases *);
void getoptions(void);
void inputoptions(void);
void allocrest(void);
void doinput(void);
void configure(void);
void prefix(chartype);
void postfix(chartype);
       
void makechar(chartype);
void dnamove_add(node *, node *, node *);
void dnamove_re_move(node **, node **);
void evaluate(node *);
void dnamove_reroot(node *);
void firstrav(node *, long);
void dnamove_hyptrav(node *, long *, long, boolean *);

void grwrite(chartype, long, long *);
void dnamove_drawline(long);
void dnamove_printree(void);
void arbitree(void);
void yourtree(void);
void initdnamovenode(node **, node **, node *, long, long, long *,
        long *, initops, pointarray, pointarray, Char *, Char *,
        FILE *);
void buildtree(void);
void setorder(void);
void mincomp(void);
       
void rearrange(void);
void dnamove_nextinc(void);
void dnamove_nextchar(void);
void dnamove_prevchar(void);
void dnamove_show(void);
void tryadd(node *, node **, node **, double *);
void addpreorder(node *, node *, node *, double *);
void try(void);
void undo(void);
void treewrite(boolean);
       
void clade(void);
void flip(long);
void changeoutgroup(void);
void redisplay(void);
void treeconstruct(void);
void maketriad(node **, long);
void newdnamove_hyptrav(node *, long *, long, long, boolean,
                        pointarray);
void prepare_node(node *p);
void dnamove_copynode(node *fromnode, node *tonode);
node *copytrav(node *p);
void numdesctrav(node *p);
void chucktree(node *p);
void copytree(void);
void makeweights(void);
void add_at(node *below, node *newtip, node *newfork);
void add_before(node *atnode, node *newtip);
void add_child(node *parent, node *newchild);
void newdnamove_hypstates(long chars, node *root, pointarray treenode);
void consolidatetree(long index);
void fliptrav(node *p, boolean recurse);
/* function prototypes */
#endif



char infilename[FNMLNGTH],intreename[FNMLNGTH],outtreename[FNMLNGTH], weightfilename[FNMLNGTH];
node *root;
long chars, screenlines, col, treelines, leftedge, topedge, vmargin,
   hscroll, vscroll, scrollinc, screenwidth, farthest, whichtree, othertree;
boolean weights, thresh, waswritten;
boolean usertree, goteof, firsttree, haslengths;  /*treeread variables*/
pointarray nodep;                                  /*treeread variables*/
node *grbg = NULL;                                  /*treeread variables*/
long *zeros;                                          /*treeread variables*/
pointptr treenode;   /* pointers to all nodes in tree */
double threshold;
double *threshwt;
boolean reversed[(long)question - (long)horiz + 1];
boolean graphic[(long)question - (long)horiz + 1];
unsigned char chh[(long)question - (long)horiz + 1];
howtree how;
gbases *garbage;
char *progname;

/* Local variables for treeconstruct, propogated global for C version: */

long dispchar, atwhat, what, fromwhere, towhere, oldoutgrno, compatible;
double like, bestyet, gotlike;
boolean display, newtree, changed, subtree, written, oldwritten, restoring,
  wasleft, oldleft, earlytree;
steptr necsteps;
boolean *in_tree;
long sett[31];
steptr numsteps;
node *nuroot;
rearrtype lastop;
Char  ch;
boolean *names;


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
    if(!(*p)->base)
      (*p)->base = (baseptr)Malloc(chars*sizeof(long));
    if(!(*p)->numnuc)
      (*p)->numnuc = (nucarray *)Malloc(endsite*sizeof(nucarray));
    if(!(*p)->numsteps)
      (*p)->numsteps = (steptr)Malloc(endsite*sizeof(long));
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


void prepare_node(node *p) {
/* This function allocates the base, numnuc and numsteps arrays for
   a node.  Because a node can change roles between tip, internal and
   ring member, all nodes need to have these in case they are used.
*/   
  p->base = (baseptr)Malloc(chars*sizeof(long));
  p->numnuc = (nucarray *)Malloc(endsite*sizeof(nucarray));
  p->numsteps = (steptr)Malloc(endsite*sizeof(long));
} /* prepare_tip */


void dnamove_gnu(gbases **p)
{
  /* this and the following are do-it-yourself garbage collectors.
     Make a new node or pull one off the garbage list */
  if (garbage != NULL) {
    *p = garbage;
    garbage = garbage->next;
  } else {
    *p = (gbases *)Malloc(sizeof(gbases));
    (*p)->base = (baseptr2)Malloc(chars*sizeof(long));
  }
  (*p)->next = NULL;
}  /* dnamove_gnu */


void dnamove_chuck(gbases *p)
{
  /* collect garbage on p -- put it on front of garbage list */
  p->next = garbage;
  garbage = p;
}  /* dnamove_chuck */


void dnamove_copynode(node *fromnode, node *tonode)
{
  /* Copy the contents of a node from fromnode to tonode. */
  int i = 0;
/*
  printf("copynode: fromnode = %d, tonode = %d\n",
    fromnode->index,tonode->index);
  printf("copynode: fromnode->base = %ld, tonode->base = %ld\n",
    fromnode->base,tonode->base);
*/
  memcpy(tonode->base, fromnode->base, chars*sizeof(long));
/*
  printf("copynode: fromnode->numnuc = %ld, tonode->numnuc = %ld\n",
    fromnode->numnuc,tonode->numnuc);
*/
  if (fromnode->numnuc != NULL)
    memcpy(tonode->numnuc, fromnode->numnuc, endsite*sizeof(nucarray));

  if (fromnode->numsteps != NULL)
    memcpy(tonode->numsteps, fromnode->numsteps, endsite*sizeof(long));

  tonode->numdesc = fromnode->numdesc;
  tonode->state   = fromnode->state;
  tonode->index   = fromnode->index;
  tonode->tip     = fromnode->tip;
  for (i=0;i<MAXNCH;i++)
    tonode->nayme[i] = fromnode->nayme[i]; 

} /* dnamove_copynode */


node *copytrav(node *p)
{
  /* Traverse the tree from p on down, copying nodes to the other tree */
  node *q, *newnode, *newnextnode, *temp;
  gnu(&grbg, &newnode);
  if(!newnode->base)
    newnode->base = (baseptr)Malloc(chars*sizeof(long));
  if(!newnode->numnuc)
    newnode->numnuc = (nucarray *)Malloc(endsite*sizeof(nucarray));
  if(!newnode->numsteps)
    newnode->numsteps = (steptr)Malloc(endsite*sizeof(long));
  dnamove_copynode(p,newnode);

  if (treenode[p->index-1] == p)
    treesets[othertree].treenode[p->index-1] = newnode;

  /* if this is a tip, return now */
  if (p->tip)
    return newnode;

  /* go around the ring, copying as we go */

  q = p->next;
  gnu(&grbg, &newnextnode);
  if(!newnextnode->base)
    newnextnode->base = (baseptr)Malloc(chars*sizeof(long));
  if(!newnextnode->numnuc)
    newnextnode->numnuc = (nucarray *)Malloc(endsite*sizeof(nucarray));
  if(!newnextnode->numsteps)
    newnextnode->numsteps = (steptr)Malloc(endsite*sizeof(long));
  dnamove_copynode(q, newnextnode);
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
      if(!newnextnode->base)
        newnextnode->base = (baseptr)Malloc(chars*sizeof(long));
      if(!newnextnode->numnuc)
        newnextnode->numnuc = (nucarray *)Malloc(endsite*sizeof(nucarray));
      if(!newnextnode->numsteps)
        newnextnode->numsteps = (steptr)Malloc(endsite*sizeof(long));
      dnamove_copynode(q, newnextnode);
      temp->next = newnextnode;
    }
  } while (q != p);
  return newnode;
} /* copytrav */


void numdesctrav(node *p)
{
  node *q;
  long childcount = 0;

  if (p->tip) {
    p->numdesc = 0;
    return;
  }

  q = p->next;

  do {
    numdesctrav(q->back);
    childcount++;
    q = q->next;
  } while (q != p);

  p->numdesc = childcount;
} /* numdesctrav */


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
} /* chucktree */


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

} /* copytree */


void getoptions()
{
  /* interactively set options */
  Char ch;
  boolean done, gotopt;
  long loopcount;

  how = arb;
  usertree = false;
  goteof = false;
  outgrno = 1;
  outgropt = false;
  thresh = false;
  weights = false;
  interleaved = true;
  loopcount = 0;
  do {
    cleerhome();
    printf("\nInteractive DNA parsimony, version %s\n\n",VERSION);
    printf("Settings for this run:\n");
    printf("  O                             Outgroup root?");
    if (outgropt)
      printf("  Yes, at sequence number%3ld\n", outgrno);
    else
      printf("  No, use as outgroup species%3ld\n", outgrno);
    printf("  W                            Sites weighted?  %s\n",
           (weights ? "Yes" : "No"));
    printf("  T                   Use Threshold parsimony?");
    if (thresh)
      printf("  Yes, count up to%4.1f per site\n", threshold);
    else
      printf("  No, use ordinary parsimony\n");
    printf("  I               Input sequences interleaved?  %s\n",
           (interleaved ? "Yes" : "No, sequential"));

    printf("  U   Initial tree (arbitrary, user, specify)?  %s\n",
           (how == arb) ? "Arbitrary"                :
           (how == use) ? "User tree from tree file" : "Tree you specify");
    printf("  0        Graphics type (IBM PC, ANSI, none)?  %s\n",
           ibmpc ? "IBM PC" : ansi  ? "ANSI"     : "(none)");
    printf("  S                  Width of terminal screen?");
    printf("%4ld\n", screenwidth);
    printf("  L                 Number of lines on screen?%4ld\n",screenlines);
    printf("\nAre these settings correct? ");
    printf("(type Y or the letter for one to change)\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    fflush(stdout);
    scanf("%c%*[^\n]", &ch);
    getchar();
    uppercase(&ch);
    done = (ch == 'Y');
    gotopt = (strchr("SOTIU0WL",ch) != NULL) ? true : false;
    if (gotopt) {
      switch (ch) {
        
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
        
        case 'I':
          interleaved = !interleaved;
          break;

        case 'W':
          weights = !weights;
          break;

        case 'U':
          if (how == arb){
            how = use;
            usertree = 1;}
          else if (how == use){
            how = spec;
            usertree = 0;}
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

  for (i = 0; i < (chars); i++)
    weight[i] = 1;
  if (weights){
      inputweights(chars, weight, &weights);
      printweights(stdout, 0, chars, weight, "Sites");
  }
  if (!thresh)
    threshold = spp;
  for (i = 0; i < (chars); i++)
    threshwt[i] = threshold * weight[i];
}  /* inputoptions */


void allocrest()
{
  long i;

  nayme = (naym *)Malloc(spp*sizeof(naym));
  in_tree = (boolean *)Malloc(nonodes*sizeof(boolean));
  weight = (steptr)Malloc(chars*sizeof(long));
  numsteps = (steptr)Malloc(chars*sizeof(long));
  necsteps = (steptr)Malloc(chars*sizeof(long));
  threshwt = (double *)Malloc(chars*sizeof(double));
  alias = (long *)Malloc(chars*sizeof(long));     /* from dnapars */
  ally = (long *)Malloc(chars*sizeof(long));      /* from dnapars */
  y = (Char **)Malloc(spp*sizeof(Char *));        /* from dnapars */
  for (i = 0; i < spp; i++)                       /* from dnapars */
    y[i] = (Char *)Malloc(chars*sizeof(Char));    /* from dnapars */
  location = (long *)Malloc(chars*sizeof(long));  /* from dnapars */
}  /* allocrest */

void makeweights()
{
  /* make up weights vector to avoid duplicate computations */
  long i;

  for (i = 1; i <= chars; i++) {
    alias[i - 1] = i;
    ally[i - 1] = i;
  }
  endsite = 0;
  for (i = 1; i <= chars; i++) { 
    if (ally[i - 1] == i) 
      endsite++; 
  } 

  for (i = 1; i <= endsite; i++)
    location[alias[i - 1] - 1] = i;
  if (!thresh)
    threshold = spp;
  zeros = (long *)Malloc(endsite*sizeof(long));
  for (i = 0; i < endsite; i++)
    zeros[i] = 0;
}  /* makeweights */


void doinput()
{
  /* reads the input data */
  inputnumbers(&spp, &chars, &nonodes, 1);
  printf("%2ld species, %3ld  sites\n", spp, chars);
  getoptions();
  printf("\nReading input file ...\n\n");
  if (weights)
    openfile(&weightfile,WEIGHTFILE,"weights file","r",progname,weightfilename);
  allocrest();
  inputoptions();
  alloctree(&treenode, nonodes, usertree);
  setuptree(treenode, nonodes, usertree);
  inputdata(chars);
  makeweights();
  makevalues(treenode, zeros, usertree);
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
    chh[(long)horiz] = 205;
    graphic[(long)horiz] = true;
    chh[(long)vert] = 186;
    graphic[(long)vert] = true;
    chh[(long)up] = 186;
    graphic[(long)up] = true;
    chh[(long)overt] = 205;
    graphic[(long)overt] = true;
    chh[(long)upcorner] = 200;
    graphic[(long)upcorner] = true;
    chh[(long)midcorner] = 204;
    graphic[(long)midcorner] = true;
    chh[(long)downcorner] = 201;
    graphic[(long)downcorner] = true;
    chh[(long)aa] = 176;
    chh[(long)cc] = 178;
    chh[(long)gg] = 177;
    chh[(long)tt] = 219;
    chh[(long)question] = '\001';
    return;
  }
  if (ansi) {
    chh[(long)horiz] = ' ';
    reversed[(long)horiz] = true;
    chh[(long)vert] = chh[(long)horiz];
    reversed[(long)vert] = true;
    chh[(long)up] = 'x';
    graphic[(long)up] = true;
    chh[(long)overt] = 'q';
    graphic[(long)overt] = true;
    chh[(long)upcorner] = 'm';
    graphic[(long)upcorner] = true;
    chh[(long)midcorner] = 't';
    graphic[(long)midcorner] = true;
    chh[(long)downcorner] = 'l';
    graphic[(long)downcorner] = true;
    chh[(long)aa] = 'a';
    reversed[(long)aa] = true;
    chh[(long)cc] = 'c';
    reversed[(long)cc] = true;
    chh[(long)gg] = 'g';
    reversed[(long)gg] = true;
    chh[(long)tt] = 't';
    reversed[(long)tt] = true;
    chh[(long)question] = '?';
    reversed[(long)question] = true;
    return;
  }
  chh[(long)horiz] = '=';
  chh[(long)vert] = ' ';
  chh[(long)up] = '!';
  chh[(long)upcorner] = '`';
  chh[(long)midcorner] = '+';
  chh[(long)downcorner] = ',';
  chh[(long)overt] = '-';
  chh[(long)aa] = 'a';
  chh[(long)cc] = 'c';
  chh[(long)gg] = 'g';
  chh[(long)tt] = 't';
  chh[(long)question] = '.';
}  /* configure */


void prefix(chartype a)
{
  /* give prefix appropriate for this character */
  if (reversed[(long)a])
    prereverse(ansi);
  if (graphic[(long)a])
    pregraph2(ansi);
}  /* prefix */


void postfix(chartype a)
{
  /* give postfix appropriate for this character */
  if (reversed[(long)a])
    postreverse(ansi);
  if (graphic[(long)a])
    postgraph2(ansi);
}  /* postfix */


void makechar(chartype a)
{
  /* print out a character with appropriate prefix or postfix */
  prefix(a);
  putchar(chh[(long)a]);
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
  if(!newnode->base)
    newnode->base = (baseptr)Malloc(chars*sizeof(long));
  if(!newnode->numnuc)
    newnode->numnuc = (nucarray *)Malloc(endsite*sizeof(nucarray));
  if(!newnode->numsteps)
    newnode->numsteps = (steptr)Malloc(endsite*sizeof(long));
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


void dnamove_re_move(node **item, node **fork)
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


void evaluate(node *r)
{
  /* determines the number of steps needed for a tree. this is
     the minimum number of steps needed to evolve sequences on
     this tree */
  long i, steps;
  double sum;

  compatible = 0;
  sum = 0.0;
  for (i = 0; i < (chars); i++)
    numsteps[i] = 0;
  /* set numdesc at each node to reflect current number of descendants */
  numdesctrav(root);
  postorder(r);
  for (i = 0; i < endsite; i++) {
    steps = r->numsteps[i];
    if (steps <= threshwt[i]) {
      sum += steps;
    } else {
      sum += threshwt[i];
    }
    if (steps <= necsteps[i] && !earlytree)
      compatible += weight[i];
  }
  like = -sum;
}  /* evaluate */


void dnamove_reroot(node *outgroup)
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
} /* dnamove_reroot */


void newdnamove_hyptrav(node *r_, long *hypset_, long b1, long b2, boolean bottom_,
                        pointarray treenode)
{
  /*  compute, print out states at one interior node */
  struct LOC_hyptrav Vars;
  long i, j, k;
  long largest;
  gbases *ancset;
  nucarray *tempnuc;
  node *p, *q;

  Vars.bottom = bottom_;
  Vars.r = r_;
  Vars.hypset = hypset_;
  dnamove_gnu(&ancset);
  tempnuc = (nucarray *)Malloc(endsite*sizeof(nucarray));
  Vars.maybe = false;
  Vars.nonzero = false;
  if (!Vars.r->tip)
    zeronumnuc(Vars.r, endsite);

  for (i = b1 - 1; i < b2; i++) {
    j = location[ally[i] - 1];
    Vars.anc = Vars.hypset[j - 1];
    if (!Vars.r->tip) {
      p = Vars.r->next;
      for (k = (long)A; k <= (long)O; k++)
        if (Vars.anc & (1 << k))
          Vars.r->numnuc[j - 1][k]++;
      do {
        for (k = (long)A; k <= (long)O; k++)
          if (p->back->base[j - 1] & (1 << k))
            Vars.r->numnuc[j - 1][k]++;
        p = p->next;
      } while (p != Vars.r);
      largest = getlargest(Vars.r->numnuc[j - 1]);
      Vars.tempset = 0;
      for (k = (long)A; k <= (long)O; k++) {
        if (Vars.r->numnuc[j - 1][k] == largest)
          Vars.tempset |= (1 << k);
      }
      Vars.r->base[j - 1] = Vars.tempset;
    }
    if (!Vars.bottom)
      Vars.anc = treenode[Vars.r->back->index - 1]->base[j - 1];
    Vars.nonzero = (Vars.nonzero || (Vars.r->base[j - 1] & Vars.anc) == 0);
    Vars.maybe = (Vars.maybe || Vars.r->base[j - 1] != Vars.anc);
  }

  j = location[ally[dispchar - 1] - 1];
  Vars.tempset = Vars.r->base[j - 1];
  Vars.anc = Vars.hypset[j - 1];
  if (!Vars.bottom)
     Vars.anc = treenode[Vars.r->back->index - 1]->base[j - 1];

  r_->state = '?';
  if (Vars.tempset == (1 << A))
    r_->state = 'A';
  if (Vars.tempset == (1 << C))
    r_->state = 'C';
  if (Vars.tempset == (1 << G))
    r_->state = 'G';
  if (Vars.tempset == (1 << T))
    r_->state = 'T';

  Vars.bottom = false;
  if (!Vars.r->tip) {
    memcpy(tempnuc, Vars.r->numnuc, endsite*sizeof(nucarray));
    q = Vars.r->next;
    do {
      memcpy(Vars.r->numnuc, tempnuc, endsite*sizeof(nucarray));
      for (i = b1 - 1; i < b2; i++) {
        j = location[ally[i] - 1];
        for (k = (long)A; k <= (long)O; k++)
          if (q->back->base[j - 1] & (1 << k))
            Vars.r->numnuc[j - 1][k]--;
        largest = getlargest(Vars.r->numnuc[j - 1]);
        ancset->base[j - 1] = 0;
        for (k = (long)A; k <= (long)O; k++)
          if (Vars.r->numnuc[j - 1][k] == largest)
            ancset->base[j - 1] |= (1 << k);
        if (!Vars.bottom)
          Vars.anc = ancset->base[j - 1];
      }
      newdnamove_hyptrav(q->back, ancset->base, b1, b2, Vars.bottom,
                treenode);
      q = q->next;
    } while (q != Vars.r);
  }
  dnamove_chuck(ancset);
}  /* newdnamove_hyptrav */


void newdnamove_hypstates(long chars, node *root, pointarray treenode)
{
  /* fill in and describe states at interior nodes */
  /* used in dnacomp, dnapars, & dnapenny */
  long i, n;
  baseptr nothing;

  /* garbage is passed along without usage to newdnamove_hyptrav,
     which also does not use it. */
  nothing = (baseptr)Malloc(endsite*sizeof(long));
  for (i = 0; i < endsite; i++)
    nothing[i] = 0;
  for (i = 1; i <= ((chars - 1) / 40 + 1); i++) {
    putc('\n', outfile);
    n = i * 40;
    if (n > chars)
      n = chars;
    newdnamove_hyptrav(root, nothing, i * 40 - 39, n, true, treenode);
  }
  free(nothing);
}  /* newdnamove_hypstates */


void grwrite(chartype c, long num, long *pos)
{
  long i;

  prefix(c);
  for (i = 1; i <= num; i++) {
    if ((*pos) >= leftedge && (*pos) - leftedge + 1 < screenwidth)
      putchar(chh[(long)c]);
    (*pos)++;
  }
  postfix(c);
}  /* grwrite */


void dnamove_drawline(long i)
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
  if (i == p->ycoord && (p == root || subtree)) {
    extra = true;
    c = overt;
    if (display) {
      switch (p->state) {
        
      case 'A':
        c = aa;
        break;
        
      case 'C':
        c = cc;
        break;
        
      case 'G':
        c = gg;
        break;
        
      case 'T':
        c = tt;
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
    if ((subtree))
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
    n = p->xcoord - q->xcoord;
    if (n < 3 && !q->tip)
      n = 3;
    if (extra) {
      n--;
      extra = false;
    }
    if (q->ycoord == i && !done) {
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
        
        case 'A':
          c = aa;
          break;
        
        case 'C':
          c = cc;
          break;
        
        case 'G':
          c = gg;
          break;
        
        case 'T':
          c = tt;
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
      if (last->ycoord > i && first->ycoord < i && i != p->ycoord) {
        c = up;
        if (i < p->ycoord)
          st = p->next->back->state;
        else
          st = p->next->next->back->state;
        if (display) {
          switch (st) {
        
          case 'A':
            c = aa;
            break;
        
          case 'C':
            c = cc;
            break;
        
          case 'G':
            c = gg;
            break;
        
          case 'T':
            c = tt;
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
  if (p->ycoord == i && p->tip) {
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
}  /* dnamove_drawline */

void dnamove_printree()
{
  /* prints out diagram of the tree */
  long tipy;
  long i, dow;

  if (!subtree)
    nuroot = root;
  if (changed || newtree)
    evaluate(root);
  if (display) {
    outfile = stdout;
    newdnamove_hypstates(chars, root, treenode);
  }
#ifdef WIN32
  if(ibmpc || ansi)
    phyClearScreen();
  else
    printf("\n");
#else
  printf((ansi || ibmpc) ? "\033[2J\033[H" : "\n");
#endif
  tipy = 1;
  dow = down;
  if (spp * dow > screenlines && !subtree)
    dow--;

  printf("  (unrooted)");
  if (display) {
    printf(" ");
    makechar(aa);
    printf(":A ");
    makechar(cc);
    printf(":C ");
    makechar(gg);
    printf(":G ");
    makechar(tt);
    printf(":T ");
    makechar(question);
    printf(":?");
  } else
    printf("                    ");
  if (!earlytree) {
    printf("%10.1f Steps", -like);
  }
  if (display)
    printf(" SITE%4ld", dispchar);
  else
    printf("         ");
  if (!earlytree) {
    printf("  %3ld sites compatible\n", compatible);
  }

  printf("                            ");
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
      dnamove_drawline(i);
  }
  if ((treelines - topedge + 1) > (screenlines - vmargin)) {
    printf("** %ld", treelines - (topedge - 1 + screenlines - vmargin));
    printf(" lines below screen **\n");
  }
  if (treelines - topedge + vmargin + 1 < screenlines)
    putchar('\n');
  gotlike = -like;
  changed = false;
}  /* dnamove_printree */


void arbitree()
{
  long i;
  root = treenode[0];
  dnamove_add(treenode[0], treenode[1], treenode[spp]);
  for (i = 3; i <= (spp); i++) {
    dnamove_add(treenode[spp + i - 3], treenode[i - 1], treenode[spp + i - 2]);
  }
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
    dnamove_printree();
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


void initdnamovenode(node **p, node **grbg, node *q, long len, long nodei,
                        long *ntips, long *parens, initops whichinit,
                        pointarray treenode, pointarray nodep, Char *str, Char *ch,
                        FILE *intree)
{
  /* initializes a node */
  /* LM 7/27  I added this function and the commented lines around */
  /* treeread() to get the program running, but all 4 move programs*/
  /* are improperly integrated into the v4.0 support files.  As is */
  /* endsite = chars and this is a patchwork function                   */
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
    /* process lengths and discard */
  default:      /*cases hslength,hsnolength,treewt,unittrwt,iter,*/
    break;
  }
} /* initdnamovenode */


void buildtree()
{
  long i, nextnode;
  node *p;
  long j;

  treeone = (node **)Malloc(maxsz*sizeof(node *));
  treetwo = (node **)Malloc(maxsz*sizeof(node *));
  treesets[othertree].treenode = treetwo;
  changed = false;
  newtree = false;
  switch (how) {

  case arb:
    treesets[othertree].treenode = treetwo;
    arbitree();
    break;

  case use:
    /* Open in binary: ftell() is broken for UNIX line-endings under WIN32 */
    openfile(&intree,intreename,"input tree file", "rb",progname,intreename);
    names = (boolean *)Malloc(spp*sizeof(boolean));
    firsttree = true; 
    nodep = NULL;
    nextnode = 0;
    haslengths = 0;
    for (i = 0; i < endsite; i++)
      zeros[i] = 0;
    treesets[whichtree].nodep = nodep;
    treeread(intree, &root, treenode, &goteof, &firsttree, nodep, &nextnode, 
                    &haslengths, &grbg, initdnamovenode,true,nonodes); 
    for (i = spp; i < (nextnode); i++) {
      p = treenode[i];
      for (j = 1; j <= 3; j++) {
        p->base = (baseptr2)Malloc(chars*sizeof(long));
        p = p->next;
      } 
    } /* debug: see comment at initdnamovenode() */

    free(names);
    FClose(intree);
    break;

  case spec:
    treesets[othertree].treenode = treetwo;
    yourtree();
    break;
  }
  if (!outgropt)
    outgrno = root->next->back->index;
  if (outgropt)
    dnamove_reroot(treenode[outgrno - 1]);
}  /* buildtree */


void setorder()
{
  /* sets in order of number of members */
  sett[0] = 1L << ((long)A);
  sett[1] = 1L << ((long)C);
  sett[2] = 1L << ((long)G);
  sett[3] = 1L << ((long)T);
  sett[4] = 1L << ((long)O);
  sett[5] = (1L << ((long)A)) | (1L << ((long)C));
  sett[6] = (1L << ((long)A)) | (1L << ((long)G));
  sett[7] = (1L << ((long)A)) | (1L << ((long)T));
  sett[8] = (1L << ((long)A)) | (1L << ((long)O));
  sett[9] = (1L << ((long)C)) | (1L << ((long)G));
  sett[10] = (1L << ((long)C)) | (1L << ((long)T));
  sett[11] = (1L << ((long)C)) | (1L << ((long)O));
  sett[12] = (1L << ((long)G)) | (1L << ((long)T));
  sett[13] = (1L << ((long)G)) | (1L << ((long)O));
  sett[14] = (1L << ((long)T)) | (1L << ((long)O));
  sett[15] = (1L << ((long)A)) | (1L << ((long)C)) | (1L << ((long)G));
  sett[16] = (1L << ((long)A)) | (1L << ((long)C)) | (1L << ((long)T));
  sett[17] = (1L << ((long)A)) | (1L << ((long)C)) | (1L << ((long)O));
  sett[18] = (1L << ((long)A)) | (1L << ((long)G)) | (1L << ((long)T));
  sett[19] = (1L << ((long)A)) | (1L << ((long)G)) | (1L << ((long)O));
  sett[20] = (1L << ((long)A)) | (1L << ((long)T)) | (1L << ((long)O));
  sett[21] = (1L << ((long)C)) | (1L << ((long)G)) | (1L << ((long)T));
  sett[22] = (1L << ((long)C)) | (1L << ((long)G)) | (1L << ((long)O));
  sett[23] = (1L << ((long)C)) | (1L << ((long)T)) | (1L << ((long)O));
  sett[24] = (1L << ((long)G)) | (1L << ((long)T)) | (1L << ((long)O));
  sett[25] = (1L << ((long)A)) | (1L << ((long)C)) | (1L << ((long)G)) |
    (1L << ((long)T));
  sett[26] = (1L << ((long)A)) | (1L << ((long)C)) | (1L << ((long)G)) |
    (1L << ((long)O));
  sett[27] = (1L << ((long)A)) | (1L << ((long)C)) | (1L << ((long)T)) |
    (1L << ((long)O));
  sett[28] = (1L << ((long)A)) | (1L << ((long)G)) | (1L << ((long)T)) |
    (1L << ((long)O));
  sett[29] = (1L << ((long)C)) | (1L << ((long)G)) | (1L << ((long)T)) |
    (1L << ((long)O));
  sett[30] = (1L << ((long)A)) | (1L << ((long)C)) | (1L << ((long)G)) |
    (1L << ((long)T)) | (1L << ((long)O));
}  /* setorder */


void mincomp()
{
  /* computes for each site the minimum number of steps
     necessary to accomodate those species already
     in the analysis */
  long i, j, k;
  boolean done;

  for (i = 0; i < (chars); i++) {
    done = false;
    j = 0;
    while (!done) {
      j++;
      done = true;
      k = 1;
      do {
        if (k < nonodes)
          done = (done && (treenode[k - 1]->base[i] & sett[j - 1]) != 0);
        k++;
      } while (k <= spp && done);
    }
    if (j == 31)
      necsteps[i] = 4;
    if (j <= 30)
      necsteps[i] = 3;
    if (j <= 25)
      necsteps[i] = 2;
    if (j <= 15)
      necsteps[i] = 1;
    if (j <= 5)
      necsteps[i] = 0;
    necsteps[i] *= weight[i];
  }
}  /* mincomp */


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
            dnamove_re_move(&treenode[i - 1], &q);
            add_child(treenode[j - 1], treenode[i - 1]);
            if (fromtype == beforenode)
              consolidatetree(q->index);
          } else
            ok2 = false;
        } else {
          if (j != root->index) { /* can't insert at root */
            changed = true;
            copytree();
            dnamove_re_move(&treenode[i - 1], &q);
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
  dnamove_printree();
  if (!(ok1 && ok2))
    printf("Not a possible rearrangement.  Try again: \n");
  else {
    written = false;
  }
}  /* rearrange */


void dnamove_nextinc()
{
  /* show next incompatible site */
  long disp0;
  boolean done;

  display = true;
  disp0 = dispchar;
  done = false;
  do {
    dispchar++;
    if (dispchar > chars) {
      dispchar = 1;
      done = (disp0 == 0);
    }
  } while (!(necsteps[dispchar - 1] != numsteps[dispchar - 1] ||
             dispchar == disp0 || done));
  dnamove_printree();
}  /* dnamove_nextinc */

void dnamove_nextchar()
{
  /* show next site */
  display = true;
  dispchar++;
  if (dispchar > chars)
    dispchar = 1;
  dnamove_printree();
}  /* dnamove_nextchar */

void dnamove_prevchar()
{
  /* show previous site */
  display = true;
  dispchar--;
  if (dispchar < 1)
    dispchar = chars;
  dnamove_printree();
}  /* dnamove_prevchar */

void dnamove_show()
{
  long i;
  boolean ok;

  do {
    printf("SHOW: (Character number or 0 to see none)? ");
    inpnum(&i, &ok);
    ok = (ok && (i == 0 || (i >= 1 && i <= chars)));
    if (ok && i != 0) {
      display = true;
      dispchar = i;
    }
    if (ok && i == 0)
      display = false;
  } while (!ok);
  dnamove_printree();
}  /* dnamove_show */


void tryadd(node *p, node **item, node **nufork, double *place)
{
  /* temporarily adds one fork and one tip to the tree.
     Records scores in ARRAY place */
  dnamove_add(p, *item, *nufork);
  evaluate(root);
  place[p->index - 1] = -like;
  dnamove_re_move(item, nufork);
}  /* tryadd */


void addpreorder(node *p, node *item_, node *nufork_, double *place)
{
  /* traverses a binary tree, calling PROCEDURE tryadd
     at a node before calling tryadd at its descendants */
  node *item, *nufork, *q;

  item = item_;
  nufork = nufork_;
  if (p == NULL)
    return;
  tryadd(p,&item,&nufork,place);
  if (!p->tip) {
    q = p->next;
    do {
      addpreorder(q->back, item,nufork,place);
      q = q->next;
    } while (q != p);
  }
}  /* addpreorder */


void try()
{
  /* Remove node, try it in all possible places */
  double *place;
  long i, j, oldcompat, saveparent;
  double current;
  node *q, *dummy, *rute;
  boolean tied, better, ok, madenode;

  madenode = false;
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
/* q = ring base of i's parent */
  q = treenode[treenode[i - 1]->back->index - 1];
  saveparent = q->index;
/* if i is a left child, fromwhere = index of right sibling (binary) */
/* if i is a right child, fromwhere = index of left sibling (binary) */
  if (q->next->back->index == i)
    fromwhere = q->next->next->back->index;
  else
    fromwhere = q->next->back->index;
  rute = root;
 
/* if root is i's parent ... */
  if (q->next->next->next == q) {
    if (root == treenode[treenode[i - 1]->back->index - 1]) {
      /* if i is left child then rute becomes right child,
         and vice-versa */
      if (treenode[treenode[i - 1]->back->index - 1]->next->back == treenode[i - 1])
        rute = treenode[treenode[i - 1]->back->index - 1]->next->next->back;
      else
        rute = treenode[treenode[i - 1]->back->index - 1]->next->back;
    }
  }

/* Remove i and perhaps its parent node from the tree.  If i is part of a
   multifurcation, *dummy will come back null.  If so, make a new internal
   node to be i's parent as it is inserted in various places around the
   tree.
*/
  dnamove_re_move(&treenode[i - 1], &dummy);
  if (dummy == NULL) {
    madenode = true;
    nonodes++;
    maketriad(&dummy, nonodes);
  }
  oldleft = wasleft;                                                
  root = rute;
  addpreorder(root, treenode[i - 1], dummy, place);
  wasleft = oldleft;
  restoring = true;
  if (madenode) {
    add_child(treenode[saveparent - 1], treenode[i - 1]);
    nonodes--;
  } else
    dnamove_add(treenode[fromwhere - 1], treenode[what - 1], q);
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
  boolean btemp;

  /* don't undo to an uninitialized tree */
  if (!treesets[othertree].initialized) {
    dnamove_printree();
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
  dnamove_printree();
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
    dnamove_printree();
  if (waswritten && ch != 'A' && ch != 'R')
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
  ok = (ok && ((unsigned)i) <= ((unsigned)nonodes));
  if (ok) {
    subtree = (i > 0);
    if (subtree)
      nuroot = treenode[i - 1];
    else
      nuroot = root;
  }
  dnamove_printree();
  if (!ok)
    printf("Not possible to use this node. ");
}  /* clade */


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
    dnamove_printree();
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
    ok = (ok && i >= 1 && i <= nonodes &&
          i != root->index);
    if (ok)
      outgrno = i;
  } while (!ok);

  copytree();
  dnamove_reroot(treenode[outgrno - 1]);
  changed = true;
  lastop = reroott;
  dnamove_printree();
  oldwritten = written;
  written = false;
}  /* changeoutgroup */


void redisplay()
{
  boolean done = false;
  waswritten = false;
  do {
    printf("NEXT? (Options: R # + - S . T U W O F H J K L C ? X Q) ");
    printf("(? for Help) ");
#ifdef WIN32
    phyFillScreenColor();
#endif
    fflush(stdout);
    scanf("%c%*[^\n]", &ch);
    getchar();
    uppercase(&ch); 
    if (strchr("HJKLCFORSTUXQ+#-.W?",ch) != NULL){
      switch (ch) {
        
      case 'R':
        rearrange();
        break;
        
      case '#':
        dnamove_nextinc();
        break;
        
      case '+':
        dnamove_nextchar();
        break;
        
      case '-':
        dnamove_prevchar();
        break;
        
      case 'S':
        dnamove_show();
        break;
        
      case '.':
        dnamove_printree();
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
        dnamove_printree();
        break;

      case 'J':
        window(downn, &leftedge, &topedge, hscroll, vscroll, treelines,
                 screenlines, screenwidth, farthest, subtree);
        dnamove_printree();
        break;

      case 'K':
        window(upp, &leftedge, &topedge, hscroll, vscroll, treelines,
                 screenlines, screenwidth, farthest, subtree);
        dnamove_printree();
        break;

      case 'L':
        window(right, &leftedge, &topedge, hscroll, vscroll, treelines,
                 screenlines, screenwidth, farthest, subtree);
        dnamove_printree();
        break;

      case 'C':
        clade();
        break;
        
      case '?':
        help("site");
        dnamove_printree();
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
  if (written)
    return;
  do {
    printf("Do you want to write out the tree to a file? (Y or N) ");
#ifdef WIN32
    phyFillScreenColor();
#endif
    fflush(stdout);
    scanf("%c%*[^\n]", &ch);
    getchar();
    if (ch == 'Y' || ch == 'y')
      treewrite(done);
  } while (ch != 'Y' && ch != 'y' && ch != 'N' && ch != 'n');
}  /* redisplay */


void treeconstruct()
{
  /* constructs a binary tree from the pointers in treenode. */
  int i;
  restoring = false;
  subtree = false;
  display = false;
  dispchar = 0;
  earlytree = true;
  waswritten = false;
  buildtree();

  /* get an accurate value for nonodes by finding out where the nodes
     really stop
   */
  for (i=0;i<nonodes;i++) {
    if (treenode[i]==NULL)
      break;
  }
  nonodes = i;

  printf("\nComputing steps needed for compatibility in sites ...\n\n");
  setorder();
  mincomp();
  newtree = true;
  earlytree = false;
  dnamove_printree();
  bestyet = -like;
  gotlike = -like;
  lastop = none;
  newtree = false;
  written = false;
  redisplay();
}  /* treeconstruct */


int main(int argc, Char *argv[])
{ /* Interactive DNA parsimony */
  /* reads in spp, chars, and the data. Then calls treeconstruct to
     construct the tree and query the user */
#ifdef MAC
  argc = 1;                /* macsetup("Dnamove","");        */
  argv[0] = "Dnamove";
#endif
  init(argc, argv);
  progname = argv[0];
  strcpy(infilename,INFILE);
  strcpy(intreename,INTREE);
  strcpy(outtreename,OUTTREE);

  openfile(&infile,infilename,"input file", "r",argv[0],infilename);
  openfile(&outtree,outtreename,"output file", "w",argv[0],outtreename);

  whichtree = 0;
  othertree = 1;
  treesets[whichtree].initialized = false;
  treesets[othertree].initialized = false;
  garbage = NULL;
  screenlines = 24;
  scrollinc = 20;
  screenwidth = 80;
  printdata = 0;
  topedge = 1;
  leftedge = 1;
  ibmpc = IBMCRT;
  ansi = ANSICRT;
  doinput();
  configure();
  treeconstruct();
  FClose(infile);
  FClose(outtree);
#ifdef MAC
  fixmacfile(outtreename);
#endif
#ifdef WIN32
  phyRestoreConsoleAttributes();
#endif
  return 0;
}  /* Interactive DNA parsimony */
