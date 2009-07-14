
#include "phylip.h"
#include "moves.h"

/* version 3.6. (c) Copyright 1993-2008 by the University of Washington.
   Written by Joseph Felsenstein and Andrew Keeffe.  Permission is granted to
   copy and use this program provided no fee is charged for it and provided
   that this copyright notice is not removed. */

/* maximum number of species               */
#define maxsp           5000

/* size of pointer array.  >= 2*maxsp - 1  */
/* (this can be large without eating memory */
#define maxsz           9999

#define overr           4
#define which           1


typedef enum {valid, remoov, quit} reslttype;

typedef enum {
  horiz, vert, up, updel, ch_over, upcorner, midcorner, downcorner, aa, cc, 
  gg, tt, deleted
} chartype;

typedef struct treeset_t {
  node *root;
  pointarray nodep;
  long nonodes;
  boolean waswritten, hasmult, haslengths, nolengths, initialized;
} treeset_t;

treeset_t treesets[2];
treeset_t simplifiedtree;

typedef enum { arb, use, spec } howtree;

typedef enum {beforenode, atnode} movet;

movet fromtype;


#ifndef OLDC
/* function prototypes */
void   initretreenode(node **, node **, node *, long, long, long *,
    long *, initops, pointarray, pointarray, Char *, Char *, FILE *);
void   gdispose(node *);
void   maketriad(node **, long);
void   maketip(node **, long);
void   copynode(node *, node *);
node  *copytrav(node *);
void   copytree(void);
void   getoptions(void);
void   configure(void);
void   prefix(chartype);

void   postfix(chartype);
void   ltrav(node *, boolean *);
boolean ifhaslengths(void);
void   add_at(node *, node *, node *);
void   add_before(node *, node *);
void   add_child(node *, node *);
void   re_move(node **, node **);
void   reroot(node *);
void   ltrav_(node *, double, double, double *, long *, long *);

void   precoord(node *, boolean *, double *, long *);
void   coordinates(node *, double, long *, long *, double *);
void   flatcoordinates(node *, long *);
void   grwrite(chartype, long, long *);
void   drawline(long, node *, boolean *);
void   printree(void);
void   togglelengths(void);
void   arbitree(void);
void   yourtree(void);
void   buildtree(void);
void   unbuildtree(void);
void   retree_help(void);
void   consolidatetree(long);
void   rearrange(void);
boolean any_deleted(node *);
void   fliptrav(node *, boolean);
void   flip(long);
void   transpose(long);
void   ifdeltrav(node *, boolean *);
double oltrav(node *);
void   outlength(void);
void   midpoint(void);
void   deltrav(node *, boolean );
void   reg_del(node *, boolean);
boolean isdeleted(long);
void   deletebranch(void);
void   restorebranch(void);
void   del_or_restore(void);
void   undo(void);
void   treetrav(node *);
void   simcopynode(node *, node *);
node  *simcopytrav(node *);
void   simcopytree(void);
void   writebranchlength(double);
void   treeout(node *, boolean, double, long);
void   maketemptriad(node **, long);
void   roottreeout(boolean *);
void   notrootedtorooted(void);
void   rootedtonotrooted(void);
void   treewrite(boolean *);
void   retree_window(adjwindow);
void   getlength(double *, reslttype *, boolean *);
void   changelength(void);
void   changename(void);
void   clade(void);
void   changeoutgroup(void);
void   redisplay(void);
void   treeconstruct(void);
void   fill_del(node*p);
/* function prototypes */
#endif


node *root, *garbage;

long nonodes, outgrno, screenwidth, vscreenwidth,
     screenlines, col, treenumber, leftedge, topedge, treelines,
     hscroll, vscroll, scrollinc, whichtree, othertree,
     numtrees, treesread;

double     trweight;
boolean    waswritten, onfirsttree, hasmult, haslengths,
           nolengths, nexus, xmltree;

node **treeone, **treetwo;
pointarray nodep;           /* pointers to all nodes in current tree */

node *grbg;

boolean    reversed[14];
boolean    graphic[14];
unsigned char       cch[14];
howtree    how;

char       intreename[FNMLNGTH], outtreename[FNMLNGTH];

boolean    subtree, written, readnext;
node      *nuroot;
Char      ch;

boolean delarray[maxsz];


void initretreenode(node **p, node **grbg, node *q, long len,
    long nodei, long *ntips, long *parens, initops whichinit,
    pointarray treenode, pointarray nodep, Char *str,
    Char *ch, FILE *intree)
{
  /* initializes a node */
  long i;
  boolean minusread;
  double valyew, divisor;

  switch (whichinit) {

    case bottom:
      gnu(grbg, p);
      (*p)->index = nodei;
      (*p)->tip = false;
      (*p)->deleted=false;
      (*p)->deadend=false;
      (*p)->onebranch=false;
      (*p)->onebranchhaslength=false;
      for (i=0;i<MAXNCH;i++)
        (*p)->nayme[i] = '\0';
      nodep[(*p)->index - 1] = (*p);
      break;

    case nonbottom:
      gnu(grbg, p);
      (*p)->index = nodei;
      break;

    case hslength:
      if ((*p)->back) {
        (*p)->back->back = *p;
        (*p)->haslength = (*p)->back->haslength;
        if ((*p)->haslength)
          (*p)->length = (*p)->back->length;
      }
      break;

    case tip:
      (*ntips)++;
      gnu(grbg, p);
      nodep[(*ntips) - 1] = *p;
      (*p)->index = *ntips;
      (*p)->tip = true;
      (*p)->hasname = true;
      strncpy ((*p)->nayme, str, MAXNCH);
      break;

    case length:
      (*p)->haslength = true;
      if ((*p)->back != NULL)
        (*p)->back->haslength = (*p)->haslength;
      processlength(&valyew, &divisor, ch,
          &minusread, intree, parens);
      if (!minusread)
        (*p)->length = valyew / divisor;
      else
        (*p)->length = 0.0;
      (*p)->back  = q;
      if (haslengths && q != NULL) {
        (*p)->back->haslength = (*p)->haslength;
        (*p)->back->length = (*p)->length;
      }
      break;

    case hsnolength:
      haslengths = (haslengths && q == NULL);
      (*p)->haslength = false;
      (*p)->back  = q;
      break;

    default:        /*cases iter, treewt, unttrwt         */
      break;        /*should not occur                */

  }
} /* initretreenode */


void gdispose(node *p)
{
  /* go through tree throwing away nodes */
  node *q, *r;

  if (p->tip)
    return;
  q = p->next;
  while (q != p) {
    gdispose(q->back);
    q->tip = false;
    q->hasname = false;
    q->haslength = false;
    r = q;
    q = q->next;
    chuck(&grbg, r);
  }
  q->tip = false;
  q->hasname = false;
  q->haslength = false;
  chuck(&grbg, q);
}  /*  gdispose */


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
  nodep[index - 1] = *p;
}  /* maketriad */


void maketip(node **p, long index)
{
  /*  Initiate a tip node */
  gnu(&grbg, p);
  (*p)->index = index;
  (*p)->tip = true;
  (*p)->hasname = false;
  (*p)->haslength = false;
  nodep[index - 1] = *p;
}  /* maketip */


void copynode(node *fromnode, node *tonode)
{
  /* Copy the contents of a node from fromnode to tonode. */
  int i;

  tonode->index   = fromnode->index;
  tonode->deleted = fromnode->deleted;
  tonode->tip     = fromnode->tip;
  tonode->hasname = fromnode->hasname;
  if (fromnode->hasname)
    for (i=0;i<MAXNCH;i++)
      tonode->nayme[i] = fromnode->nayme[i]; 
  tonode->haslength = fromnode->haslength;
  if (fromnode->haslength)
    tonode->length = fromnode->length;

} /* copynode */


node *copytrav(node *p)
{
  /* Traverse the tree from p on down, copying nodes to the other tree */
  node *q, *newnode, *newnextnode, *temp;

  gnu(&grbg, &newnode);
  copynode(p,newnode);

  if (nodep[p->index-1] == p)
    treesets[othertree].nodep[p->index-1] = newnode;

  /* if this is a tip, return now */
  if (p->tip)
    return newnode;

  /* go around the ring, copying as we go */

  q = p->next;
  gnu(&grbg, &newnextnode);
  copynode(q, newnextnode);
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
      copynode(q, newnextnode);
      temp->next = newnextnode;
    }
  } while (q != p);
  return newnode;
} /* copytrav */


void copytree()
{
  /* Make a complete copy of the current tree for undo purposes */
  if (whichtree == 1)
    othertree = 0;
  else
    othertree = 1;

  treesets[othertree].root = copytrav(root);
  treesets[othertree].nonodes = nonodes;
  treesets[othertree].waswritten = waswritten;
  treesets[othertree].hasmult = hasmult;
  treesets[othertree].haslengths = haslengths;
  treesets[othertree].nolengths = nolengths;
  treesets[othertree].initialized = true;

} /* copytree */


void getoptions()
{
  /* interactively set options */
  long loopcount;
  Char ch;
  boolean done, gotopt;

  how = use;
  outgrno = 1;
  loopcount = 0;
  onfirsttree = true;
  do {
    cleerhome();
    printf("\nTree Rearrangement, version %s\n\n",VERSION);
    printf("Settings for this run:\n");
    printf("  U          Initial tree (arbitrary, user, specify)?");
    if (how == arb)
      printf("  Arbitrary\n");
    else if (how == use)
      printf("  User tree from tree file\n");
    else
      printf("  Tree you specify\n");
    printf("  N   Format to write out trees (PHYLIP, Nexus, XML)?");
    if (nexus)
      printf("  Nexus\n");
    else {
      if (xmltree)
        printf("  XML\n");
      else
        printf("  PHYLIP\n");
    }
    printf("  0               Graphics type (IBM PC, ANSI, none)?");
    if (ibmpc)
      printf("  IBM PC\n");
    if (ansi )
      printf("  ANSI\n");
    if (!(ibmpc || ansi))
      printf("  (none)\n");
    printf("  W       Width of terminal screen, of plotting area?");
    printf("%4ld, %2ld\n", screenwidth, vscreenwidth);
    printf("  L                        Number of lines on screen?");
    printf("%4ld\n", screenlines);
    printf("\nAre these settings correct?");
    printf(" (type Y or the letter for one to change)\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    fflush(stdout);
    scanf("%c%*[^\n]", &ch);
    getchar();
    if (ch == '\n')
      ch = ' ';
    ch = (isupper(ch)) ? ch : toupper(ch);
    done = (ch == 'Y');
    gotopt = (ch == 'U' || ch == 'N' || ch == '0' || ch == 'L' || ch == 'W');
    if (gotopt) {
      switch (ch) {

        case 'U':
          if (how == arb)
            how = use;
          else if (how == use)
            how = spec;
          else
            how = arb;
          break;

        case 'N':
          if (nexus) {
            nexus = false;
            xmltree = true;
          }
          else if (xmltree)
            xmltree = false;
          else
            nexus = true;
          break;

        case '0':
          initterminal(&ibmpc, &ansi);
          break;

        case 'L':
          initnumlines(&screenlines);
          break;

        case 'W':
          screenwidth= readlong("Width of terminal screen (in characters)?\n");
          vscreenwidth=readlong("Width of plotting area (in characters)?\n");
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


void configure()
{
  /* configure to machine -- set up special characters */
  chartype a;

  for (a = horiz; (long)a <= (long)deleted; a = (chartype)((long)a + 1))
    reversed[(long)a] = false;
  for (a = horiz; (long)a <= (long)deleted; a = (chartype)((long)a + 1))
    graphic[(long)a] = false;
  cch[(long)deleted] = '.';
  cch[(long)updel] = ':';
  if (ibmpc) {
    cch[(long)horiz] = '>';
    cch[(long)vert] = 186;
    graphic[(long)vert] = true;
    cch[(long)up] = 186;
    graphic[(long)up] = true;
    cch[(long)ch_over] = 205;
    graphic[(long)ch_over] = true;
    cch[(long)upcorner] = 200;
    graphic[(long)upcorner] = true;
    cch[(long)midcorner] = 204;
    graphic[(long)midcorner] = true;
    cch[(long)downcorner] = 201;
    graphic[(long)downcorner] = true;
    return;
  }
  if (ansi) {
    cch[(long)horiz] = '>';
    cch[(long)vert] = cch[(long)horiz];
    reversed[(long)vert] = true;
    cch[(long)up] = 'x';
    graphic[(long)up] = true;
    cch[(long)ch_over] = 'q';
    graphic[(long)ch_over] = true;
    cch[(long)upcorner] = 'm';
    graphic[(long)upcorner] = true;
    cch[(long)midcorner] = 't';
    graphic[(long)midcorner] = true;
    cch[(long)downcorner] = 'l';
    graphic[(long)downcorner] = true;
    return;
  }
  cch[(long)horiz] = '>';
  cch[(long)vert] = ' ';
  cch[(long)up] = '!';
  cch[(long)upcorner] = '`';
  cch[(long)midcorner] = '+';
  cch[(long)downcorner] = ',';
  cch[(long)ch_over] = '-';
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


void ltrav(node *p, boolean *localhl)
{
  /* Traversal function for ifhaslengths() */
  node *q;

  if (p->tip) {
    (*localhl) = ((*localhl) && p->haslength);
    return;
  }
  q = p->next;
  do {
    (*localhl) = ((*localhl) && q->haslength);
    if ((*localhl))
      ltrav(q->back, localhl);
    q = q->next;
  } while (p != q);
}  /* ltrav */


boolean ifhaslengths()
{
  /* return true if every branch in tree has a length */
  boolean localhl;
  localhl = true;
  ltrav(root, &localhl);
  return localhl;
}  /* ifhaslengths */


void add_at(node *below, node *newtip, node *newfork)
{
  /* inserts the nodes newfork and its left descendant, newtip,
     to the tree.  below becomes newfork's right descendant */
  node *leftdesc, *rtdesc;
  double length;

  if (below != nodep[below->index - 1])
    below = nodep[below->index - 1];

  if (newfork == NULL) {
    nonodes++;
    maketriad (&newfork, nonodes);
    if (haslengths) {
      newfork->haslength = true;
      newfork->next->haslength = true;
      newfork->next->next->haslength = true;
    }
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
  if (!haslengths)
    return;
  if (newfork->back != NULL) {
    length = newfork->back->length / 2.0;
    newfork->length = length;
    newfork->back->length = length;
    below->length = length;
    below->back->length = length;
  } else {
    length = newtip->length / 2.0;
    newtip->length = length;
    newtip->back->length = length;
    below->length = length;
    below->back->length = length;
    below->haslength = true;
  }
  newtip->back->length = newtip->length;
}  /* add_at */


void add_before(node *atnode, node *newtip)
{
  /* inserts the node newtip together with its ancestral fork
     into the tree next to the node atnode. */
  /*xx ?? debug what to do if no ancestral node -- have to create one */
  /*xx this case is handled by add_at.  However, add_at does not account for
    when there is more than one sibling for the relocated newtip */
  node *q;

  if (atnode != nodep[atnode->index - 1])
    atnode = nodep[atnode->index - 1];
  q = nodep[newtip->index-1]->back;
  if (q != NULL) {
    q = nodep[q->index-1];
    if (newtip == q->next->next->back) {
      q->next->back = newtip;
      newtip->back = q->next;
      q->next->next->back = NULL;
    }
  }
  if (newtip->back != NULL) {
    add_at(atnode, newtip, nodep[newtip->back->index-1]);
  } else {
    add_at(atnode, newtip, NULL);
  }
}  /* add_before */


void add_child(node *parent, node *newchild)
{
  /* adds the node newchild into the tree as the last child of parent */

  int i;
  node *newnode, *q;

  if (parent != nodep[parent->index - 1])
    parent = nodep[parent->index - 1];
  gnu(&grbg, &newnode);
  newnode->tip = false;
  newnode->deleted=false;
  newnode->deadend=false;
  newnode->onebranch=false;
  newnode->onebranchhaslength=false;
  for (i=0;i<MAXNCH;i++)
    newnode->nayme[i] = '\0';
  newnode->index = parent->index;
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


void re_move(node **item, node **fork)
{
  /* Removes node item from the tree.  If item has one sibling,
     removes its ancestor, fork, from the tree as well and attach
     item's sib to fork's ancestor.  In this case, it returns a pointer
     to the removed fork node which is still attached to item.
     */
  node *p =NULL, *q;
  int nodecount;

  if ((*item)->back == NULL) {
    *fork = NULL;
    return;
  }
  *fork = nodep[(*item)->back->index - 1];
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
    /*xx*/ *fork = NULL;
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
}  /* re_move */


void reroot(node *outgroup)
{
  /* Reorient tree so that outgroup is by itself on the left of the root */ 
  node *p, *q, *r;
  long nodecount = 0;
  double templen;

  q = root->next;
  do {                    /* when this loop exits, p points to the internal */
    p = q;                /* node to the right of root */
    nodecount++;
    q = p->next;
  } while (q != root);
  r = p;

  /* There is no point in proceeding if 
     1. outgroup is a child of root, and
     2. the tree bifurcates at the root.
     */
  if((outgroup->back->index == root->index) && !(nodecount > 2))
    return;

  /* reorient nodep array

     The nodep array must point to the ring member of each ring
     that is closest to the root.  The while loop changes the ring member
     pointed to by nodep[] for those nodes that will have their
     orientation changed by the reroot operation.
     */
  p = outgroup->back;
  while (p->index != root->index) {
    q = nodep[p->index - 1]->back;
    nodep[p->index - 1] = p;
    p = q;
  }
  if (nodecount > 2)
    nodep[p->index - 1] = p;

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
     right child by dividing outgroup's length */
  if (haslengths) {
    templen = outgroup->length / 2.0;
    outgroup->length = templen;
    outgroup->back->length = templen;
    root->next->next->length = templen;
    root->next->next->back->length = templen;
  }
} /* reroot */


void ltrav_(node *p, double lengthsum, double lmin, double *tipmax,
    long *across, long *maxchar)
{
  node *q;
  long rchar, nl;
  double sublength;

  if (p->tip) {
    if (lengthsum > (*tipmax))
      (*tipmax) = lengthsum;
    if (lmin == 0.0)
      return;
    rchar = (long)(lengthsum / (*tipmax) * (*across) + 0.5);

    nl = strlen(nodep[p->index - 1]->nayme);
    if (rchar + nl > (*maxchar))
      (*across) = (*maxchar) - (long)(nl * (*tipmax) / lengthsum + 0.5);
    return;
  }
  q = p->next;
  do {
    if (q->length >= lmin)
      sublength = q->length;
    else
      sublength = lmin;
    ltrav_(q->back, lengthsum + sublength, lmin, tipmax, across, maxchar);
    q = q->next;
  } while (p != q);
}  /* ltrav */


void precoord(node *nuroot,boolean *subtree,double *tipmax,long *across)
{
  /* set tipmax and across so that tree is scaled to screenwidth */

  double oldtipmax, minimum;
  long i, maxchar;

  (*tipmax) = 0.0;
  if ((*subtree))
    maxchar = vscreenwidth - 13;
  else
    maxchar = vscreenwidth - 5;
  (*across) = maxchar;
  ltrav_(nuroot, 0.0, 0.0, tipmax, across, &maxchar);
  i = 0;
  do {
    oldtipmax = (*tipmax);
    minimum = 3.0 / (*across) * (*tipmax);
    ltrav_(nuroot, 0.0, minimum, tipmax, across, &maxchar);
    i++;
  } while (fabs((*tipmax) - oldtipmax) > 0.01 * oldtipmax && i <= 40);
}  /* precoord */


void coordinates(node *p, double lengthsum, long *across, long *tipy,
    double *tipmax)
{
  /* establishes coordinates of nodes for display with lengths */
  node *q, *first, *last;

  if (p->tip) {
    p->xcoord = (long)((*across) * lengthsum / (*tipmax) + 0.5);
    p->ycoord = (*tipy);
    p->ymin   = (*tipy);
    p->ymax   = (*tipy);
    (*tipy)  += down;
    return;
  }
  q = p->next;
  do {
    coordinates(q->back, lengthsum + q->length, across, tipy, tipmax);
    q = q->next;
  } while (p != q);
  first = p->next->back;
  q = p;
  while (q->next != p)
    q = q->next;
  last = q->back;
  p->xcoord = (long)((*across) * lengthsum / (*tipmax) + 0.5);
  if (p == root) {
    if (root->next->next->next == root)
      p->ycoord = (first->ycoord + last->ycoord) / 2;
    else
      p->ycoord = p->next->next->back->ycoord;
  }
  else
    p->ycoord = (first->ycoord + last->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* coordinates */


void flatcoordinates(node *p, long *tipy)
{
  /* establishes coordinates of nodes for display without lengths */
  node *q, *first, *last;

  if (p->tip) {
    p->xcoord = 0;
    p->ycoord = (*tipy);
    p->ymin   = (*tipy);
    p->ymax   = (*tipy);
    (*tipy) += down;
    return;
  }
  q = p->next;
  do {
    flatcoordinates(q->back, tipy);
    q = q->next;
  } while (p != q);
  first = p->next->back;
  q = p->next;
  while (q->next != p)
    q = q->next;
  last = q->back;
  p->xcoord = (last->ymax - first->ymin) * 3 / 2;
  if (p == root) {
    if (root->next->next->next == root)
      p->ycoord = (first->ycoord + last->ycoord) / 2;
    else
      p->ycoord = p->next->next->back->ycoord;
  }
  else
    p->ycoord = (first->ycoord + last->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* flatcoordinates */


void grwrite(chartype c, long num, long *pos)
{
  long i;

  prefix(c);
  for (i = 1; i <= num; i++) {
    if ((*pos) >= leftedge && (*pos) - leftedge + 1 < screenwidth)
      putchar(cch[(long)c]);
    (*pos)++;
  }
  postfix(c);
}  /* grwrite */


void drawline(long i, node *nuroot, boolean *subtree)
{
  /* draws one row of the tree diagram by moving up tree */
  long pos;
  node *p, *q, *r, *s, *first =NULL, *last =NULL;
  long n, j;
  long up_nondel, down_nondel;
  boolean extra, done;
  chartype c, d;
  pos = 1;
  p = nuroot;
  q = nuroot;

  extra = false;
  if (i == (long)p->ycoord && (p == root || (*subtree))) {
    c = ch_over;
    if ((*subtree))
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
    if ((*subtree))
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
    if (haslengths && !nolengths)
      n = (long)(q->xcoord - p->xcoord);
    else
      n = (long)(p->xcoord - q->xcoord);
    if (n < 3 && !q->tip)
      n = 3;
    if (extra) {
      n--;
      extra = false;
    }
    if ((long)q->ycoord == i && !done) {
      c = ch_over;
      if (!haslengths && !q->haslength)
        c = horiz;
      if (q->deleted)
        c = deleted;
      if (q == first)
        d = downcorner;
      else if (q == last)
        d = upcorner;
      else if ((long)q->ycoord == (long)p->ycoord)
        d = c;
      else
        d = midcorner;
      if (n > 1 || q->tip) {
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
      if ((long)last->ycoord > i && (long)first->ycoord < i &&
          i != (long)p->ycoord) {
        c = up;
        if(p->deleted)
          c = updel;
        if (!p->tip) {
          up_nondel = 0;
          down_nondel = 0;
          r = p->next;
          do {
            s = r->back;
            if ((long)s->ycoord < (long)p->ycoord && !s->deleted)
              up_nondel = (long)s->ycoord;
            if (s->ycoord > p->ycoord && !s->deleted && 
                (down_nondel == 0))
              down_nondel = (long)s->ycoord;
            if (i < (long)p->ycoord && s->deleted && i > (long)s->ycoord)
              c = updel;
            if (i > (long)p->ycoord && s->deleted && i < (long)s->ycoord)
              c = updel;
            r = r->next;
          } while (r != p);

          if ((up_nondel != 0) && i < (long)p->ycoord && i > up_nondel)
            c = up;
          if ((down_nondel != 0) && i > (long)p->ycoord && i < down_nondel)
            c = up;
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
    if (p->hasname) {
      n = 0;
      for (j = 1; j <= MAXNCH; j++) {
        if (nodep[p->index - 1]->nayme[j - 1] != '\0')
          n = j;
      }
      chwrite(':', 1, &pos, leftedge, screenwidth);
      for (j = 0; j < n; j++)
        chwrite(nodep[p->index - 1]->nayme[j],
            1, &pos, leftedge, screenwidth);
    }
  }
  putchar('\n');
}  /* drawline */


void printree()
{
  /* prints out diagram of the tree */
  long across;
  long tipy;
  double tipmax;
  long i, dow, vmargin;

  haslengths = ifhaslengths();
  if (!subtree)
    nuroot = root;
  cleerhome();
  tipy = 1;
  dow = down;
  if (spp * dow > screenlines && !subtree) {
    dow--;
  }
  if (haslengths && !nolengths) {
    precoord(nuroot, &subtree, &tipmax, &across);
    /* protect coordinates() from div/0 errors if user decides to
       examine a tip as a subtree */
    if (tipmax == 0)
      tipmax = 0.01;
    coordinates(nuroot, 0.0, &across, &tipy, &tipmax);
  } else
    flatcoordinates(nuroot, &tipy);
  vmargin = 2;
  treelines = tipy - dow;
  if (topedge != 1) {
    printf("** %ld lines above screen **\n", topedge - 1);
    vmargin++;
  }
  if ((treelines - topedge + 1) > (screenlines - vmargin))
    vmargin++;
  for (i = 1; i <= treelines; i++) {
    if (i >= topedge && i < topedge + screenlines - vmargin)
      drawline(i, nuroot,&subtree);
  }
  if (leftedge > 1)
    printf("** %ld characters to left of screen ", leftedge);
  if ((treelines - topedge + 1) > (screenlines - vmargin)) {
    printf("** %ld", treelines - (topedge - 1 + screenlines - vmargin));
    printf(" lines below screen **\n");
  }
  if (treelines - topedge + vmargin + 1 < screenlines)
    putchar('\n');
}  /* printree */


void togglelengths()
{
  nolengths = !nolengths;
  printree();
}  /* togglengths */


void arbitree()
{
  long i, maxinput;
  node *newtip, *newfork;

  maxinput = 1;
  do {
    spp = readlong("How many species?\n");
    maxinput++;
    if (maxinput == 100) {
      printf("ERROR: too many tries at choosing species\n");
      exxit(-1);
    }
  } while (spp <= 0);
  nonodes = spp * 2 - 1;
  maketip(&root, 1);
  maketip(&newtip, 2);
  maketriad(&newfork, spp + 1);
  add_at(root, newtip, newfork);
  for (i = 3; i <= spp; i++) {
    maketip(&newtip, i);
    maketriad(&newfork, spp + i - 1);
    add_at(nodep[spp + i - 3], newtip, newfork);
  }
}  /* arbitree */


void yourtree()
{
  long uniquearray[maxsz];
  long uniqueindex = 0;
  long i, j, k, k_max, maxinput;
  boolean ok, done;
  node *newtip, *newfork;
  Char ch;

  uniquearray[0] = 0;
  spp = 2;
  nonodes = spp * 2 - 1;
  maketip(&root, 1);
  maketip(&newtip, 2);
  maketriad(&newfork, spp + 3);
  add_at(root, newtip, newfork);
  i = 2;
  maxinput = 1;
  k_max = 5;
  do {
    i++;
    printree();
    printf("Enter 0 to stop building tree.\n");
    printf("Add species%3ld", i);
    do {
      printf("\n at or before which node (type number): ");
      inpnum(&j, &ok);
      ok = (ok && ((unsigned long)j < i || (j > spp + 2 && j < spp + i + 1)));
      if (!ok)
        printf("Impossible number. Please try again:\n");
      maxinput++;
      if (maxinput == 100) {
        printf("ERROR: too many tries at choosing number\n");
        exxit(-1);
      }
    } while (!ok);
    maxinput = 1;
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
        maxinput++;
        if (maxinput == 100) {
          printf("ERROR: too many tries at choosing option\n");
          exxit(-1);
        }
      } while (ch != 'A' && ch != 'B');
    }
    else ch = 'B';   /* if user has chosen a tip, set Before */

    if (j != 0) {
      if (ch == 'A') {
        if (!nodep[j - 1]->tip) {
          maketip(&newtip, i);
          add_child(nodep[j - 1], nodep[i - 1]);
        }
      } else {
        maketip(&newtip, i);
        maketriad(&newfork, spp + i + 1);
        nodep[i-1]->back = newfork;
        newfork->back = nodep[i-1];
        add_before(nodep[j - 1], nodep[i - 1]);
      } /* endif (before or at node) */
    }

    done = (j == 0);
    if (!done) {
      if (ch == 'B')
        k = spp * 2 + 3;
      else
        k = spp * 2 + 2;

      k_max = k;
      do {
        if (nodep[k - 2] != NULL) {
          nodep[k - 1] = nodep[k - 2];
          nodep[k - 1]->index = k;
          nodep[k - 1]->next->index = k;
          nodep[k - 1]->next->next->index = k;
        }
        k--;
      } while (k != spp + 3); 
      if (j > spp + 1)
        j++;
      spp++;
      nonodes = spp * 2 - 1;
    }
  } while (!done);

  for (i = spp + 1; i <= k_max; i++) {
    if ((nodep[i - 1] != nodep[i]) && (nodep[i - 1] != NULL)) {
      uniquearray[uniqueindex++] = i;
      uniquearray[uniqueindex] = 0;
    }
  }

  for ( i = 0; uniquearray[i] != 0; i++) {
    nodep[spp + i] = nodep[uniquearray[i] - 1];
    nodep[spp + i]->index = spp + i + 1;
    nodep[spp + i]->next->index = spp + i + 1;
    nodep[spp + i]->next->next->index = spp + i + 1;
  }
  for (i = spp + uniqueindex; i <= k_max; i++)
    nodep[i] = NULL;

  nonodes = spp * 2 - 1;
}  /* yourtree */


void buildtree()
{
  /* variables needed to be passed to treeread() */
  long    nextnode   = 0;
  pointarray dummy_treenode=NULL;  /* Ignore what happens to this */
  boolean goteof     = false;
  boolean haslengths = false;
  boolean firsttree;
  node *p, *q;
  long nodecount = 0;


  /* These assignments moved from treeconstruct -- they seem to happen
     only here. */
  /*xx treeone & treetwo assignments should probably happen in 
    treeconstruct.  Memory leak if user reads multiple trees. */
  treeone = (node **)Malloc(maxsz*sizeof(node *));
  treetwo = (node **)Malloc(maxsz*sizeof(node *));
  simplifiedtree.nodep = (node **)Malloc(maxsz*sizeof(node *));
  subtree     = false;
  topedge     = 1;
  leftedge    = 1;
  switch (how) {

    case arb:
      nodep = treeone;
      treesets[othertree].nodep = treetwo;
      arbitree();
      break;

    case use:
      printf("\nReading tree file ...\n\n");

      if (!readnext) {
        /* This is the first time through here, act accordingly */
        firsttree = true;
        /* Open in binary: ftell() is broken for UNIX line-endings under WIN32 */
        openfile(&intree,INTREE,"input tree file", "rb","retree",intreename);
        numtrees = countsemic(&intree);
        treesread = 0;
      } else {
        /* This isn't the first time through here ... */
        firsttree = false;
      }
      allocate_nodep(&nodep, &intree, &spp);
      treesets[whichtree].nodep = nodep;

      if (firsttree)
        nayme = (naym *)Malloc(spp*sizeof(naym));
      treeread(intree, &root, dummy_treenode, &goteof, &firsttree, nodep, 
          &nextnode, &haslengths, &grbg, initretreenode,true,-1);
      nonodes = nextnode;
      treesread++;
      treesets[othertree].nodep = treetwo;
      break;

    case spec:
      nodep = treeone;
      treesets[othertree].nodep = treetwo;
      yourtree();
      break;
  }

  q = root->next;
  do {
    p = q;
    nodecount++;
    q = p->next;
  } while (q != root);

  outgrno = root->next->back->index;
  if(!(nodecount > 2)) {
    reroot(nodep[outgrno - 1]);
  }
}  /* buildtree */


void unbuildtree()
{
  /* throw all nodes of the tree onto the garbage heap */
  long i;

  gdispose(root);
  for (i = 0; i < nonodes; i++)
    nodep[i] = NULL;
}  /* unbuildtree */


void retree_help()
{
  /* display help information */
  char tmp[100];
  printf("\n\n . Redisplay the same tree again\n");
  if (haslengths) {
    printf(" = Redisplay the same tree with");
    if (!nolengths)
      printf("out/with");
    else
      printf("/without");
    printf(" lengths\n");
  }
  printf(" U Undo the most recent change in the tree\n");
  printf(" W Write tree to a file\n");
  printf(" + Read next tree from file (may blow up if none is there)\n");
  printf("\n");
  printf(" R Rearrange a tree by moving a node or group\n");
  printf(" O select an Outgroup for the tree\n");
  if (haslengths)
    printf(" M Midpoint root the tree\n");
  printf(" T Transpose immediate branches at a node\n");
  printf(" F Flip (rotate) subtree at a node\n");
  printf(" D Delete or restore nodes\n");
  printf(" B Change or specify the length of a branch\n");
  printf(" N Change or specify the name(s) of tip(s)\n");
  printf("\n");
  printf(" H Move viewing window to the left\n");
  printf(" J Move viewing window downward\n");
  printf(" K Move viewing window upward\n");
  printf(" L Move viewing window to the right\n");
  printf(" C show only one Clade (subtree) (might be useful if tree is ");
  printf("too big)\n");
  printf(" ? Help (this screen)\n");
  printf(" Q (Quit) Exit from program\n");
  printf(" X Exit from program\n\n");
  printf(" TO CONTINUE, PRESS ON THE Return OR Enter KEY");
  getstryng(tmp);
  printree();
}  /* retree_help */


void consolidatetree(long index)
{
  node *start, *r, *q;
  int i;

  start = nodep[index - 1];
  q = start->next;

  while (q != start) {
    r = q;
    q = q->next;
    chuck(&grbg, r);
  } 
  chuck(&grbg, q);

  i = index;
  while (nodep[i-1] != NULL) {
    r = nodep[i - 1];
    if (!(r->tip))
      r->index--;
    if (!(r->tip)) {
      q = r->next;
      do {
        q->index--;
        q = q->next;
      } while (r != q && q != NULL);
    }
    nodep[i - 1] = nodep[i];
    i++;
  }

  nonodes--;
} /* consolidatetree */


void rearrange()
{
  long i, j, maxinput;
  boolean ok;
  node *p, *q;
  char ch;

  printf("Remove everything to the right of which node? ");
  inpnum(&i, &ok);
  if ( ok == false ) {
    /* fall through */
  } else if ( i < 1 || i > spp*2 - 1 ) {
    /* i is not in range */
    ok = false;
  } else if (i == root->index ) {
    /* i is root */
    ok = false;
  } else if ( nodep[i-1]->deleted ) {
    /* i has been deleted */
    ok = false;
  } else {
    printf("Add at or before which node? ");
    inpnum(&j, &ok);
    if ( ok == false ) {
      /* fall through */
    } else if ( j < 1 || j > spp*2 - 1 ) {
      /* j is not in range */
      ok = false;
    } else if ( nodep[j-1]->deleted ) {
      /* j has been deleted */
      ok = false;
    } else if (j != root->index && nodep[nodep[j-1]->back->index - 1]->deleted ) {
      /* parent of j has been deleted */
      ok = false;
    } else if ( nodep[j-1] == nodep[nodep[i-1]->back->index -1] ) {
      /* i is j's parent */
      ok = false;
    } else {
      /* make sure that j is not a descendant of i */
      for ( p = nodep[j-1]; p != root; p = nodep[p->back->index - 1] ) {
        if ( p == nodep[i-1] ) {
          ok = false;
          break;
        }
      }
      if ( ok ) {
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
          ch = toupper(ch);
          maxinput++;
          if (maxinput == 100) {
            printf("ERROR: Input failed too many times.\n");
            exxit(-1);
          }
        } while (ch != 'A' && ch != 'B');

        if (ch == 'A') {
          if ( nodep[j - 1]->deleted || nodep[j - 1]->tip ) {
            /* If j is a tip or has been deleted */
            ok = false;
          } else if ( nodep[j-1] == nodep[nodep[i-1]->back->index -1] ) {
            /* If j is i's parent */
            ok = false;
          } else {
            copytree();
            re_move(&nodep[i - 1], &q);
            add_child(nodep[j - 1], nodep[i - 1]);
            if (fromtype == beforenode)
              consolidatetree(q->index);
          }
        } else { /* ch == 'B' */
          if (j == root->index) { /* can't insert at root */
            ok = false;
          } else {
            copytree();
            printf("Insert before node %ld\n",j);
            re_move(&nodep[i - 1], &q);
            if (q != NULL) {
              nodep[q->index-1]->next->back = nodep[i-1];
              nodep[i-1]->back = nodep[q->index-1]->next;
            }
            add_before(nodep[j - 1], nodep[i - 1]);
          }
        } /* endif (before or at node) */
      } /* endif (ok to do move) */
    } /* endif (destination node ok) */
  } /* endif (from node ok) */

  printree();

  if ( !ok )
    printf("Not a possible rearrangement.  Try again: \n");
  else {
    written = false;
  }
}  /* rearrange */


boolean any_deleted(node *p)
{
  /* return true if there are any deleted branches from branch on down */
  boolean localdl;
  localdl = false;
  ifdeltrav(p, &localdl);
  return localdl;
}  /* any_deleted */


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
    if (ok)
      ok = !any_deleted(nodep[i - 1]);
  } else {
    i = atnode;
    ok = true;
  }
  if (ok) {
    copytree();
    fliptrav(nodep[i - 1], true);
  }
  if (atnode == 0)
    printree();
  if (ok) {
    written = false;
    return;
  }
  if ((i >= 1 && i <= spp) ||
      (i > spp && i <= nonodes && any_deleted(nodep[i - 1])))
    printf("Can't flip there. ");
  else
    printf("No such node. ");
}  /* flip */


void transpose(long atnode)
{
  /* flip at a node left-right */
  long i;
  boolean ok;

  if (atnode == 0) {
    printf("Transpose branches at which node? ");
    inpnum(&i, &ok);
    ok = (ok && i > spp && i <= nonodes);
    if (ok)
      ok = !nodep[i - 1]->deleted;
  } else {
    i = atnode;
    ok = true;
  }
  if (ok) {
    copytree();
    fliptrav(nodep[i - 1], false);
  }
  if (atnode == 0)
    printree();
  if (ok) {
    written = false;
    return;
  }
  if ((i >= 1 && i <= spp) ||
      (i > spp && i <= nonodes && nodep[i - 1]->deleted))
    printf("Can't transpose there. ");
  else
    printf("No such node. ");
}  /* transpose */


void ifdeltrav(node *p, boolean *localdl)
{
  node *q;

  if (*localdl)
    return;

  if (p->tip) {
    (*localdl) = ((*localdl) || p->deleted);
    return;
  }
  q = p->next;
  do {
    (*localdl) = ((*localdl) || q->deleted);
    ifdeltrav(q->back, localdl);
    q = q->next;
  } while (p != q);
}  /* ifdeltrav */


double oltrav(node *p)
{
  node *q;
  double maxlen, templen;
  if (p->deleted)
    return 0.0;
  if (p->tip) {
    p->beyond = 0.0;
    return 0.0;
  } else {
    q = p->next;
    maxlen = 0;
    do {
      templen = q->back->deleted ? 0.0 : q->length + oltrav(q->back);
      maxlen = (maxlen > templen) ? maxlen : templen;
      q->beyond = templen;
      q = q->next;
    } while (p != q);
    p->beyond = maxlen;
    return (maxlen);
  }
}  /* oltrav */


void outlength()
{
  /* compute the farthest combined length out from each node */
  oltrav(root);
}  /* outlength */


void midpoint()
{
  /* midpoint root the tree */
  double balance, greatlen, lesslen, grlen, maxlen;
  node *maxnode, *grnode, *lsnode =NULL;
  boolean ok = true;
  boolean changed = false;
  node *p, *q;
  long nodecount = 0;
  boolean multi = false;

  copytree();
  p = root;
  outlength();
  q = p->next;
  greatlen = 0;
  grnode = q->back;
  lesslen = 0;

  q = root->next;
  do {
    p = q;
    nodecount++;
    q = p->next;
  } while (q != root);
  if (nodecount > 2)
    multi = true;

  /* Find the two greatest lengths reaching from root to tips.
     Also find the lengths and node pointers of the first nodes in the 
     direction of those two greatest lengths.  */
  p = root;
  q = root->next;
  do {
    if (greatlen <= q->beyond) {
      lesslen = greatlen;
      lsnode = grnode;
      greatlen = q->beyond;
      grnode = q->back;
    }
    if ((greatlen > q->beyond) && (q->beyond >= lesslen)) {
      lesslen = q->beyond;
      lsnode = q->back;
    }
    q = q->next;
  } while (p != q);

  /* If we don't have two non-deleted nodes to balance between then
     we can't midpoint root the tree */
  if (grnode->deleted || lsnode->deleted || grnode == lsnode)
    ok = false;
  balance = greatlen - (greatlen + lesslen) / 2.0;
  grlen = grnode->length;

  while ((balance - grlen > 1e-10) && ok) {
    /* First, find the most distant immediate child of grnode
       and reroot to it. */
    p = grnode;
    q = p->next;
    maxlen = 0;
    maxnode = q->back;
    do {
      if (maxlen <= q->beyond) {
        maxlen = q->beyond;
        maxnode = q->back;
      }
      q = q->next;
    } while (p != q);
    reroot(maxnode);
    changed = true;

    /* Reassess the situation, using the same "find the two greatest
       lengths" code as occurs before the while loop.  If another reroot
       is necessary, this while loop will repeat. */
    p = root;
    outlength();
    q = p->next;
    greatlen = 0;
    grnode = q->back;
    lesslen = 0;
    do {
      if (greatlen <= q->beyond) {
        lesslen = greatlen;
        lsnode = grnode;
        greatlen = q->beyond;
        grnode = q->back;
      }
      if ((greatlen > q->beyond) && (q->beyond >= lesslen)) {
        lesslen = q->beyond;
        lsnode = q->back;
      }
      q = q->next;
    } while (p != q);
    if (grnode->deleted || lsnode->deleted || grnode == lsnode)
      ok = false;
    balance = greatlen - (greatlen + lesslen) / 2.0;
    grlen = grnode->length;
  }; /* end of while ((balance > grlen) && ok) */

  if (ok) {
    /*xx the following ignores deleted nodes */
    /*   this may be ok because deleted nodes are omitted from length calculations */
    if (multi) {
      reroot(grnode); /*xx need length corrections */

      p = root;
      outlength();
      q = p->next;
      greatlen = 0;
      grnode = q->back;
      lesslen = 0;

      do {
        if (greatlen <= q->beyond) {
          lesslen = greatlen;
          lsnode = grnode;
          greatlen = q->beyond;
          grnode = q->back;
        }
        if ((greatlen > q->beyond) && (q->beyond >= lesslen)) {
          lesslen = q->beyond;
          lsnode = q->back;
        }
        q = q->next;
      } while (p != q);
      balance = greatlen - (greatlen + lesslen) / 2.0;
    }
    grnode->length -= balance;
    if (((grnode->length) < 0.0) && (grnode->length > -1.0e-10))
      grnode->length = 0.0;
    grnode->back->length = grnode->length;
    lsnode->length += balance;
    if (((lsnode->length) < 0.0) && (lsnode->length > -1.0e-10))
      lsnode->length = 0.0;
    lsnode->back->length = lsnode->length;
  }
  printree();
  if (ok) {
    if (any_deleted(root))
      printf("Deleted nodes were not used in midpoint calculations.\n");
  }
  else {
    printf("Can't perform midpoint because of deleted branches.\n");
    if (changed) {
      undo();
      printf("Tree restored to original state.  Undo information lost.\n");
    }
  }
} /* midpoint */


void deltrav(node *p, boolean value)
{
  /* register p and p's children as deleted or extant, depending on value */
  node *q;

  p->deleted = value;

  if (p->tip)
    return;

  q = p->next;
  do {
    deltrav(q->back, value);
    q = q->next;
  } while (p != q);
}  /* deltrav */


void fill_del(node*p)
{
  int alldell;
  node *q = p;

  if ( p->next == NULL) return;

  q=p->next;
  while ( q != p) {
    fill_del(q->back); 
    q=q->next;
  }

  alldell = 1;

  q=p->next;
  while ( q != p) {
    if ( !q->back->deleted ) {
      alldell = 0;
    }
    q=q->next;
  }

  p->deleted = alldell;

}



void reg_del(node *delp, boolean value)
{
  /* register delp and all of delp's children as deleted */
  deltrav(delp, value);
  fill_del(root);
}  /* reg_del */


boolean isdeleted(long nodenum)
{
  /* true if nodenum is a node number in a deleted branch */
  return(nodep[nodenum - 1]->deleted);
} /* isdeleted */


void deletebranch()
{
  /* delete a node */
  long i;
  boolean ok1;

  printf("Delete everything to the right of which node? ");
  inpnum(&i, &ok1);
  ok1 = (ok1 && i >= 1 && i <= nonodes && i != root->index && !isdeleted(i));
  if (ok1) {
    copytree();
    reg_del(nodep[i - 1],true);
  }
  printree();
  if (!ok1)
    printf("Not a possible deletion.  Try again.\n");
  else {
    written = false;
  }
}  /* deletebranch */


void restorebranch()
{
  /* restore deleted branches */
  long i;
  boolean ok1;

  printf("Restore everything to the right of which node? ");
  inpnum(&i, &ok1);
  ok1 = (ok1 && i >= 1 && i < spp * 2 && isdeleted(i) && 
      ( i == root->index || !nodep[nodep[i - 1]->back->index - 1]->deleted));

  if (ok1) {
    reg_del(nodep[i - 1],false);
  }
  printree();
  if (!ok1)
    printf("Not a possible restoration.  Try again: \n");
  else {
    written = false;
  }
} /* restorebranch */


void del_or_restore()
{
  /* delete or restore a branch */
  long maxinput;
  Char ch;

  if (any_deleted(root)) {
    maxinput = 1;
    do {
      printf("Enter D to delete a branch\n");
      printf("OR enter R to restore a branch: ");
#ifdef WIN32
      phyFillScreenColor();
#endif
      fflush(stdout);
      scanf("%c%*[^\n]", &ch);
      getchar();
      if (ch == '\n')
        ch = ' ';
      ch = (isupper(ch)) ? ch : toupper(ch);
      maxinput++;
      if (maxinput == 100) {
        printf("ERROR: too many tries at choosing option\n");
        exxit(-1);
      }
    } while (ch != 'D' && ch != 'R');
    if (ch == 'R')
      restorebranch();
    else
      deletebranch();
  }
  else
    deletebranch();
} /* del_or_restore */


void undo()
{
  /* don't undo to an uninitialized tree */
  if (!treesets[othertree].initialized) {
    printree();
    printf("Nothing to undo.\n");
    return;
  }

  treesets[whichtree].root = root;
  treesets[whichtree].nodep = nodep;
  treesets[whichtree].nonodes = nonodes;
  treesets[whichtree].waswritten = waswritten;
  treesets[whichtree].hasmult = hasmult;
  treesets[whichtree].haslengths = haslengths;
  treesets[whichtree].nolengths = nolengths;
  treesets[whichtree].initialized = true;

  whichtree = othertree;

  root = treesets[whichtree].root;
  nodep = treesets[whichtree].nodep;
  nonodes = treesets[whichtree].nonodes;
  waswritten = treesets[whichtree].waswritten;
  hasmult = treesets[whichtree].hasmult;
  haslengths = treesets[whichtree].haslengths;
  nolengths = treesets[whichtree].nolengths;

  if (othertree == 0)
    othertree = 1;
  else
    othertree = 0;

  printree();
}  /* undo */

/* 
   These attributes of nodes in the tree are modified by treetrav()
   in preparation for writing a tree to disk.

   boolean deadend      This node is not deleted but all of its children
   are, so this node will be treated as such when
   the tree is written or displayed.

   boolean onebranch    This node has only one valid child, so that this
   node will not be written and its child will be 
   written as a child of its grandparent with the
   appropriate summing of lengths.

   nodep *onebranchnode
   Used if onebranch is true.  Onebranchnode points
   to the one valid child.  This child may be one or
   more generations down from the current node.

   double onebranchlength
   Used if onebranch is true.  Onebranchlength is
   the length from the current node to the valid 
   child.
   */

void treetrav(node *p)
{
  long branchcount = 0;
  node *q, *onebranchp =NULL;

  /* Count the non-deleted branches hanging off of this node into branchcount.
     If there is only one such branch, onebranchp points to that branch. */

  if (p->tip)
    return;

  q = p->next;
  do {
    if (!q->back->deleted) {
      if (!q->back->tip)
        treetrav(q->back);

      if (!q->back->deadend && !q->back->deleted) {
        branchcount++;
        onebranchp = q->back;
      }
    }
    q = q->next;
  } while (p != q);

  if (branchcount == 0)
    p->deadend = true;
  else
    p->deadend = false;
  p->onebranch = false;
  if (branchcount == 1 && onebranchp->tip) {
    p->onebranch = true;
    p->onebranchnode = onebranchp;
    p->onebranchhaslength = (p->haslength || (p == root))
      && onebranchp->haslength;
    if (p->onebranchhaslength)
      p->onebranchlength = onebranchp->length + p->length;
  }
  if (branchcount == 1 && !onebranchp->tip) {
    p->onebranch = true;
    if (onebranchp->onebranch) {
      p->onebranchnode = onebranchp->onebranchnode;
      p->onebranchhaslength = (p->haslength || (p == root))
        && onebranchp->onebranchhaslength;
      if (p->onebranchhaslength)
        p->onebranchlength = onebranchp->onebranchlength + p->length;
    } else {
      p->onebranchnode = onebranchp;
      p->onebranchhaslength = p->haslength && onebranchp->haslength;
      if (p->onebranchhaslength)
        p->onebranchlength = onebranchp->length + p->length;
    }
  }
} /* treetrav */


void simcopynode(node *fromnode, node *tonode)
{
  /* Copy the contents of a node from fromnode to tonode. */
  int i;

  tonode->index   = fromnode->index;
  tonode->deleted = fromnode->deleted;
  tonode->tip     = fromnode->tip;
  tonode->hasname = fromnode->hasname;
  if (fromnode->hasname)
    for (i=0;i<MAXNCH;i++)
      tonode->nayme[i] = fromnode->nayme[i]; 
  tonode->haslength = fromnode->haslength;
  if (fromnode->haslength)
    tonode->length = fromnode->length;

} /* simcopynode */


node *simcopytrav(node *p)
{
  /* Traverse the tree from p on down, copying nodes to the other tree */
  node *q, *newnode, *newnextnode, *temp;
  long lastnodeidx = 0;

  gnu(&grbg, &newnode);
  simcopynode(p, newnode);

  if (nodep[p->index - 1] == p)
    simplifiedtree.nodep[p->index - 1] = newnode;

  /* if this is a tip, return now */
  if (p->tip)
    return newnode;
  if (p->onebranch && p->onebranchnode->tip) {
    simcopynode(p->onebranchnode, newnode);
    if (p->onebranchhaslength)
      newnode->length = p->onebranchlength;
    return newnode;
  } else if (p->onebranch && !p->onebranchnode->tip) {
    /* recurse down p->onebranchnode */
    p->onebranchnode->length = p->onebranchlength;
    p->onebranchnode->haslength = p->onebranchnode->haslength;
    return simcopytrav(p->onebranchnode);
  } else {
    /* Multiple non-deleted branch case:  go round the node recursing
       down the branches. Don't go down deleted branches or dead ends. */
    q = p->next;
    while (q != p) {
      if (!q->back->deleted && !q->back->deadend)
        lastnodeidx = q->back->index;
      q = q->next;
    }

    q = p->next;
    gnu(&grbg, &newnextnode);
    simcopynode(q, newnextnode);
    newnode->next = newnextnode;
    do {
      /* If branch is deleted or is a dead end, do not recurse
         down the branch. */
      if (!q->back->deleted && !q->back->deadend) {
        newnextnode->back = simcopytrav(q->back);
        newnextnode->back->back = newnextnode;
        q = q->next;
        if (newnextnode->back->index == lastnodeidx) {
          newnextnode->next = newnode;
          break;
        }
        if (q == p) {
          newnextnode->next = newnode;
        } else {
          temp = newnextnode;
          gnu(&grbg, &newnextnode);
          simcopynode(q, newnextnode);
          temp->next = newnextnode;
        }
      } else { /*xx this else and q=q->next are experimental 
                 (seems to be working) */
        q = q->next;
      }

    } while (q != p);
  }
  return newnode;
}  /* simcopytrav */


void simcopytree()
{
  /* Make a simplified copy of the current tree for rooting/unrooting
     on output.  Deleted notes are removed and lengths are consolidated. */

  simplifiedtree.root = simcopytrav(root);
  /*xx If there are deleted nodes, nonodes will be different.
    However, nonodes is not used in the simplified tree. */
  simplifiedtree.nonodes = nonodes;
  simplifiedtree.waswritten = waswritten;
  simplifiedtree.hasmult = hasmult;
  simplifiedtree.haslengths = haslengths;
  simplifiedtree.nolengths = nolengths;
  simplifiedtree.initialized = true;
} /* simcopytree */


void writebranchlength(double x)
{
  long w;

  /* write branch length onto output file, keeping track of what
     column of line you are in, and writing to correct precision */

  if (x > 0.0)
    w = (long)(0.43429448222 * log(x));
  else if (x == 0.0)
    w = 0;
  else
    w = (long)(0.43429448222 * log(-x)) + 1;
  if (w < 0)
    w = 0;
  if ((long)(100000*x) == 100000*(long)x) {
    if (!xmltree)
      putc(':', outtree);
    fprintf(outtree, "%*.1f", (int)(w + 2), x);
    col += w + 3;
  }
  else {
    if ((long)(100000*x) == 10000*(long)(10*x)) {
      if (!xmltree)
        putc(':', outtree);
      fprintf(outtree, "%*.1f", (int)(w + 3), x);
      col += w + 4;
    }
    else {
      if ((long)(100000*x) == 1000*(long)(100*x)) {
        if (!xmltree)
          putc(':', outtree);
        fprintf(outtree, "%*.2f", (int)(w + 4), x);
        col += w + 5;
      }
      else {
        if ((long)(100000*x) == 100*(long)(1000*x)) {
          if (!xmltree)
            putc(':', outtree);
          fprintf(outtree, "%*.3f", (int)(w + 5), x);
          col += w + 6;
        }
        else {
          if ((long)(100000*x) == 10*(long)(10000*x)) {
            if (!xmltree)
              putc(':', outtree);
            fprintf(outtree, "%*.4f", (int)(w + 6), x);
            col += w + 7;
          }
          else {
            if (!xmltree)
              putc(':', outtree);
            fprintf(outtree, "%*.5f", (int)(w + 7), x);
            col += w + 8;
          }
        }
      }
    }
  }
} /* writebranchlength */


void treeout(node *p, boolean writeparens, double addlength, long indent)
{
  /* write out file with representation of final tree */
  long i, n, lastnodeidx = 0;
  Char c;
  double x;
  boolean comma;
  node *q;

  /* If this is a tip or there are no non-deleted branches from this node,
     render this node as a tip (write its name).  */

  if (p == root) {
    if (xmltree)
      indent = 0;
    else
      indent = 0;
    if (xmltree) {
      fprintf(outtree, "<phylogeny>");  /* assumes no length at root! */
    } else putc('(', outtree);
  }
  if (p->tip) {
    if (p->hasname) {
      n = 0;
      for (i = 1; i <= MAXNCH; i++) {
        if ((nodep[p->index - 1]->nayme[i - 1] != '\0')
            && (nodep[p->index - 1]->nayme[i - 1] != ' '))
          n = i;
      }
      indent += 2;
      if (xmltree) {
        putc('\n', outtree);
        for (i = 1; i <= indent; i++)
          putc(' ', outtree);
        fprintf(outtree, "<clade");
        if (p->haslength) {
          fprintf(outtree, " length=\"");
          x = p->length;
          writebranchlength(x);
          fprintf(outtree,"\"");
        }
        putc('>', outtree);
        fprintf(outtree, "<name>");
      }
      for (i = 0; i < n; i++) {
        c = nodep[p->index - 1]->nayme[i];
        if (c == ' ')
          c = '_';
        putc(c, outtree);
      }
      col += n;
      if (xmltree)
        fprintf(outtree, "</name></clade>");
    }
  } else if (p->onebranch && p->onebranchnode->tip) { 
    if (p->onebranchnode->hasname) {
      n = 0;
      for (i = 1; i <= MAXNCH; i++) {
        if ((nodep[p->index - 1]->nayme[i - 1] != '\0')
            && (nodep[p->index - 1]->nayme[i - 1] != ' '))
          n = i;
        indent += 2;
        if (xmltree) {
          putc('\n', outtree);
          for (i = 1; i <= indent; i++)
            putc(' ', outtree);
          fprintf(outtree, "<clade");
          if ((p->haslength && writeparens) || p->onebranch) {
            if (!(p->onebranch && !p->onebranchhaslength)) {
              fprintf(outtree, " length=");
              if (p->onebranch)
                x = p->onebranchlength;
              else
                x = p->length;
              x += addlength;
              writebranchlength(x);
            }
            fprintf(outtree, "<name>");
          }
        }
        for (i = 0; i < n; i++) {
          c = p->onebranchnode->nayme[i];
          if (c == '_')
            c = ' ';
          putc(c, outtree);
        }
        col += n;
        if (xmltree)
          fprintf(outtree, "</name></clade>");
      }
    }
  } else if (p->onebranch && !p->onebranchnode->tip) {
    treeout(p->onebranchnode, true, 0.0, indent);
  } else {
    /* Multiple non-deleted branch case:  go round the node
       recursing down the branches.  */
    if (xmltree) {
      putc('\n', outtree);
      indent += 2;
      for (i = 1; i <= indent; i++)
        putc(' ', outtree);
      if (p == root)
        fprintf(outtree, "<clade>");
    }
    if (p != root) {
      if (xmltree) {
        fprintf(outtree, "<clade");
        if ((p->haslength && writeparens) || p->onebranch) {
          if (!(p->onebranch && !p->onebranchhaslength)) {
            fprintf(outtree, " length=\"");
            if (p->onebranch)
              x = p->onebranchlength;
            else
              x = p->length;
            x += addlength;
            writebranchlength(x);
          }
          fprintf(outtree, "\">");
        }
        else fprintf(outtree, ">");
      } else putc('(', outtree); 
    }
    (col)++;
    q = p->next; 
    while (q != p) {
      if (!q->back->deleted && !q->back->deadend)
        lastnodeidx = q->back->index;
      q = q->next;
    }
    q = p->next; 
    while (q != p) {
      comma = true;
      /* If branch is deleted or is a dead end, do not recurse
         down the branch and do not write a comma afterwards.
         */
      if (!q->back->deleted && !q->back->deadend)
        treeout(q->back, true, 0.0, indent);
      else
        comma = false;
      if (q->back->index == lastnodeidx)
        comma = false;
      q = q->next;
      if (q == p)
        break;
      if ((q->next == p) && (q->back->deleted || q->back->deadend))
        break;
      if (comma && !xmltree)
        putc(',', outtree);
      (col)++;
      if ((!xmltree) && col > 65) {
        putc('\n', outtree);
        col = 0;
      }
    }
    /* The right paren ')' closes off this level of recursion. */
    if (p != root) {
      if (xmltree) {
        fprintf(outtree, "\n");
        for (i = 1; i <= indent; i++)
          putc(' ', outtree);
      }
      if (xmltree) { 
        fprintf(outtree, "</clade>");
      } else putc(')', outtree);
    }
    (col)++;
  }
  if (!xmltree)
    if ((p->haslength && writeparens) || p->onebranch) {
      if (!(p->onebranch && !p->onebranchhaslength)) {
        if (p->onebranch)
          x = p->onebranchlength;
        else
          x = p->length;
        x += addlength;
        writebranchlength(x);
      }
    }
  if (p == root) {
    if (xmltree) {
      fprintf(outtree, "\n  </clade>\n</phylogeny>\n");
    } else putc(')', outtree); 
  }
}  /* treeout */


void maketemptriad(node **p, long index)
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
}  /* maketemptriad */


void roottreeout(boolean *userwantsrooted)
{
  /* write out file with representation of final tree */
  long trnum, trnumwide;
  boolean treeisrooted = false;

  treetrav(root);
  simcopytree(); /* Prepare a copy of the going tree without deleted branches */
  treesets[whichtree].root = root;        /* Store the current root */

  if (nexus) {
    trnum = treenumber;
    trnumwide = 1;
    while (trnum >= 10) {
      trnum /= 10;
      trnumwide++;
    }
    fprintf(outtree, "TREE PHYLIP_%*ld = ", (int)trnumwide, treenumber);
    if (!(*userwantsrooted))
      fprintf(outtree, "[&U] ");
    else
      fprintf(outtree, "[&R] ");
    col += 15;
  }
  root = simplifiedtree.root;                /* Point root at simplified tree */
  root->haslength = false;              /* Root should not have a length */
  if (root->tip)
    treeisrooted = true;
  else {
    if (root->next->next->next == root)
      treeisrooted = true;
    else 
      treeisrooted = false;
  }
  if (*userwantsrooted && !treeisrooted)
    notrootedtorooted();
  if (!(*userwantsrooted) && treeisrooted)
    rootedtonotrooted();
  if ((*userwantsrooted && treeisrooted) ||
      (!(*userwantsrooted) && !treeisrooted)) {
    treeout(root,true,0.0, 0);
  }
  root = treesets[whichtree].root;     /* Point root at original (real) tree */
  if (!xmltree) {
    if (hasmult)
      fprintf(outtree, "[%6.4f];\n", trweight);
    else
      fprintf(outtree, ";\n");
  }
}  /* roottreeout */


void notrootedtorooted()
{
  node *newbase, *temp;

  /* root halfway along leftmost branch of unrooted tree */
  /* create a new triad for the new base */
  maketemptriad(&newbase,nonodes+1);

  /* Take left branch and make it the left branch of newbase */
  newbase->next->back = root->next->back;
  newbase->next->next->back = root;
  /* If needed, divide length between left and right branches */
  if (newbase->next->back->haslength) {
    newbase->next->back->length /= 2.0;
    newbase->next->next->back->length =
      newbase->next->back->length;
    newbase->next->next->back->haslength = true;
  }
  /* remove leftmost ring node from old base ring */
  temp = root->next->next;
  chuck(&grbg, root->next);
  root->next = temp;
  /* point root at new base and write the tree */
  root = newbase;
  treeout(root,true,0.0, 0);
  /* (since tree mods are to simplified tree and will not be used
     for general purpose tree editing, much initialization can be
     skipped.) */
} /* notrootedtorooted */


void rootedtonotrooted()
{
  node *q, *r, *temp, *newbase;
  boolean sumhaslength = false;
  double sumlength = 0;

  /* Use the leftmost non-tip immediate descendant of the root,
     root at that, write a multifurcation with that as the base.
     If both descendants are tips, write tree as is. */
  root = simplifiedtree.root;
  /* first, search for leftmost non-tip immediate descendent of root */
  q = root->next->back;
  r = root->next->next->back;
  if (q->tip && r->tip) {
    treeout(root,true,0.0, 0);
  } else if (!(q->tip)) {
    /* allocate new base pointer */
    gnu(&grbg,&newbase);
    newbase->next = q->next;
    q->next = newbase;
    q->back = r;
    r->back = q;
    if (q->haslength && r->haslength) {
      sumlength = q->length + r->length;
      sumhaslength = true;
    }
    if (sumhaslength) {
      q->length = sumlength;
      q->back->length = sumlength;
    } else {
      q->haslength = false;
      r->haslength = false;
    }
    chuck(&grbg, root->next->next);
    chuck(&grbg, root->next);
    chuck(&grbg, root);
    root = newbase;
    treeout(root, true, 0.0, 0);
  } else if (q-tip && !(r->tip)) {
    temp = r;
    do {
      temp = temp->next;
    } while (temp->next != r);
    gnu(&grbg,&newbase);
    newbase->next = temp->next;
    temp->next = newbase;      
    q->back = r;
    r->back = q;
    if (q->haslength && r->haslength) {
      sumlength = q->length + r->length;
      sumhaslength = true;
    }
    if (sumhaslength) {
      q->length = sumlength;
      q->back->length = sumlength;
    } else {
      q->haslength = false;
      r->haslength = false;
    }
    chuck(&grbg, root->next->next);
    chuck(&grbg, root->next);
    chuck(&grbg, root);
    root = newbase;
    treeout(root, true, 0.0, 0);
  }
} /* rootedtonotrooted */


void treewrite(boolean *done)
{
  /* write out tree to a file */
  long maxinput;
  boolean rooted;

  if ( root->deleted ) {
    printf("Cannot write tree because every branch in the tree is deleted\n");
    return;
  }
  openfile(&outtree,OUTTREE,"output tree file","w","retree",outtreename);
  if (nexus && onfirsttree) {
    fprintf(outtree, "#NEXUS\n");
    fprintf(outtree, "BEGIN TREES\n");
    fprintf(outtree, "TRANSLATE;\n"); /* MacClade needs this */
  }
  if (xmltree && onfirsttree) {
    fprintf(outtree, "<phylogenies>\n");
  }
  onfirsttree = false;
  maxinput = 1;
  do {
    printf("Enter R if the tree is to be rooted\n");
    printf("OR enter U if the tree is to be unrooted: ");
#ifdef WIN32
    phyFillScreenColor();
#endif
    fflush(stdout);
    scanf("%c%*[^\n]", &ch);
    getchar();
    if (ch == '\n')
      ch = ' ';
    ch = (isupper(ch)) ? ch : toupper(ch);
    maxinput++;
    if (maxinput == 100) {
      printf("ERROR: too many tries at choosing option\n");
      exxit(-1);
    }
  } while (ch != 'R' && ch != 'U');
  col = 0;
  rooted = (ch == 'R');
  roottreeout(&rooted);
  treenumber++;
  printf("\nTree written to file \"%s\"\n\n", outtreename);
  waswritten = true;
  written = true;
  if (!(*done))
    printree();
  FClose(outtree);
}  /* treewrite */


void retree_window(adjwindow action)
{
  /* move viewing window of tree */
  switch (action) {

    case left:
      if (leftedge != 1)
        leftedge -= hscroll;
      break;

    case downn:
      /* The 'topedge + 3' is needed to allow downward scrolling
         when part of the tree is above the screen and only 1 or 2 lines
         are below it. */
      if (treelines - topedge + 3 >= screenlines)
        topedge += vscroll;
      break;

    case upp:
      if (topedge != 1)
        topedge -= vscroll;
      break;

    case right:
      if (leftedge < vscreenwidth+2)
        if (hscroll > leftedge - vscreenwidth + 1)
          leftedge = vscreenwidth;
        else
          leftedge += hscroll;
      break;
  }
  printree();
}  /* retree_window */


void getlength(double *length, reslttype *reslt, boolean *hslngth)
{
  long maxinput;
  double valyew;
  char tmp[100];

  valyew = 0.0;
  maxinput = 1;
  do {
    printf("\nEnter the new branch length\n");
    printf("OR enter U to leave the length unchanged\n");
    if (*hslngth)
      printf("OR enter R to remove the length from this branch: \n");
    getstryng(tmp);

    if (tmp[0] == 'u' || tmp[0] == 'U'){
      *reslt = quit;
      break;
    }
    else if (tmp[0] == 'r' || tmp[0] == 'R') {
      (*reslt) = remoov;
      break;}
    else if (sscanf(tmp,"%lf",&valyew) == 1){
      (*reslt) = valid;
      break;}
      maxinput++;
      if (maxinput == 100) {
        printf("ERROR: too many tries at choosing option\n");
        exxit(-1);
      }
  } while (1);
  (*length) = valyew;
}  /* getlength */


void changelength()
{
  /* change or specify the length of a tip */
  boolean hslngth;
  boolean ok;
  long i, w, maxinput;
  double length, x;
  Char ch;
  reslttype reslt;
  node *p;

  maxinput = 1;
  do {
    printf("Specify length of which branch (0 = all branches)? ");
    inpnum(&i, &ok);
    ok = (ok && (unsigned long)i <= nonodes);
    if (ok && (i != 0))
      ok = (ok && !nodep[i - 1]->deleted);
    if (i == 0)
      ok = (nodep[i - 1] != root);
    maxinput++;
    if (maxinput == 100) {
      printf("ERROR: too many tries at choosing option\n");
      exxit(-1);
    }
  } while (!ok);
  if (i != 0) {
    p = nodep[i - 1];
    putchar('\n');
    if (p->haslength) {
      x = p->length;
      if (x > 0.0)
        w = (long)(0.43429448222 * log(x));
      else if (x == 0.0)
        w = 0;
      else
        w = (long)(0.43429448222 * log(-x)) + 1;
      if (w < 0)
        w = 0;
      printf("The current length of this branch is %*.5f\n", (int)(w + 7), x);
    } else
      printf("This branch does not have a length\n");
    hslngth = p->haslength;
    getlength(&length, &reslt, &hslngth);
    switch (reslt) {

      case valid:
        copytree();

        p->length = length;
        p->haslength = true;
        if (p->back != NULL) {
          p->back->length = length;
          p->back->haslength = true;
        }
        break;

      case remoov:
        copytree();

        p->haslength = false;
        if (p->back != NULL)
          p->back->haslength = false;
        break;

      case quit:
        /* blank case */
        break;
    }
  } else {
    printf("\n (this operation cannot be undone)\n");
    maxinput = 1;
    do {
      printf("\n   enter U to leave the lengths unchanged\n");
      printf("OR enter R to remove the lengths from all branches: \n");
#ifdef WIN32
      phyFillScreenColor();
#endif
      fflush(stdout);
      scanf("%c%*[^\n]", &ch);
      getchar();
      if (ch == '\n')
        ch = ' ';
      maxinput++;
      if (maxinput == 100) {
        printf("ERROR: too many tries at choosing option\n");
        exxit(-1);
      }
    } while (ch != 'U' && ch != 'u' && ch != 'R' && ch != 'r');
    if (ch == 'R' || ch == 'r') {
      copytree();
      for (i = 0; i < spp; i++)
        nodep[i]->haslength = false;
      for (i = spp; i < nonodes; i++) {
        if (nodep[i] != NULL) {
          nodep[i]->haslength = false;
          nodep[i]->next->haslength = false;
          nodep[i]->next->next->haslength = false;
        }
      }
    }
  }
  printree();
}  /* changelength */


void changename()
{
  /* change or specify the name of a tip */
  boolean ok;
  long i, n, tipno;
  char tipname[100];

  for(;;) {
    for(;;) {
      printf("Specify name of which tip? (enter its number or 0 to quit): ");
      inpnum(&i, &ok);
      if (i > 0 && ((unsigned long)i <= spp) && ok)
        if (!nodep[i - 1]->deleted) {
          tipno = i;
          break;
        }
      if (i == 0) {
        tipno = 0;
        break;
      }
    }
    if (tipno == 0)
      break;
    if (nodep[tipno - 1]->hasname) {
      n = 0;

      /* this is valid because names are padded out to MAXNCH with nulls */
      for (i = 1; i <= MAXNCH; i++) {
        if (nodep[tipno - 1]->nayme[i - 1] != '\0')
          n = i;
      }
      printf("The current name of tip %ld is \"", tipno);
      for (i = 0; i < n; i++)
        putchar(nodep[tipno - 1]->nayme[i]);
      printf("\"\n");
    }
    copytree();
    for (i = 0; i < MAXNCH; i++)
      nodep[tipno - 1]->nayme[i] = ' ';
    printf("Enter new tip name: ");
    i = 1;
    getstryng(tipname);
    strncpy(nodep[tipno-1]->nayme,tipname,MAXNCH);
    nodep[tipno - 1]->hasname = true;
    printree();
  }
  printree();
}  /* changename */


void clade()
{
  /* pick a subtree and show only that on screen */
  long i;
  boolean ok;

  printf("Select subtree rooted at which node (0 for whole tree)? ");
  inpnum(&i, &ok);
  ok = (ok && (unsigned long)i <= nonodes);
  if (ok) {
    subtree = (i > 0);
    if (subtree)
      nuroot = nodep[i - 1];
    else
      nuroot = root;
  }
  printree();
  if (!ok)
    printf("Not possible to use this node. ");
}  /* clade */


void changeoutgroup()
{
  long i, maxinput;
  boolean ok;

  maxinput = 1;
  do {
    printf("Which node should be the new outgroup? ");
    inpnum(&i, &ok);
    ok = (ok && i >= 1 && i <= nonodes && i != root->index);
    if (ok)
      ok = (ok && !nodep[i - 1]->deleted);
    if (ok)
      ok = !nodep[nodep[i - 1]->back->index - 1]->deleted;
    if (ok)
      outgrno = i;
    maxinput++;
    if (maxinput == 100) {
      printf("ERROR: too many tries at choosing option\n");
      exxit(-1);
    }
  } while (!ok);
  copytree();
  reroot(nodep[outgrno - 1]);
  printree();
  written = false;
}  /* changeoutgroup */


void redisplay()
{
  long maxinput;
  boolean done;
  char    ch;

  done = false;
  maxinput = 1;
  do {
    printf("\nNEXT? (Options: R . ");
    if (haslengths)
      printf("= ");
    printf("U W O ");
    if (haslengths)
      printf("M ");
    printf("T F D B N H J K L C + ? X Q) (? for Help) ");
#ifdef WIN32
    phyFillScreenColor();
#endif
    fflush(stdout);
    scanf("%c%*[^\n]", &ch);
    getchar();
    if (ch == '\n')
      ch = ' ';
    ch = isupper(ch) ? ch : toupper(ch);
    if (ch == 'C' || ch == 'F' || ch == 'O' || ch == 'R' ||
        ch == 'U' || ch == 'X' || ch == 'Q' || ch == '.' ||
        ch == 'W' || ch == 'B' || ch == 'N' || ch == '?' ||
        ch == 'H' || ch == 'J' || ch == 'K' || ch == 'L' ||
        ch == '+' || ch == 'T' || ch == 'D' ||
        (haslengths && ch == 'M') ||
        (haslengths && ch == '=')) {

      switch (ch) {

        case 'R':
          rearrange();
          break;

        case '.':
          printree();
          break;

        case '=':
          togglelengths();
          break;

        case 'U':
          undo();
          break;

        case 'W':
          treewrite(&done);
          break;

        case 'O':
          changeoutgroup();
          break;

        case 'M':
          midpoint();
          break;

        case 'T':
          transpose(0);
          break;

        case 'F':
          flip(0);
          break;

        case 'C':
          clade();
          break;

        case 'D':
          del_or_restore();
          break;

        case 'B':
          changelength();
          break;

        case 'N':
          changename();
          break;

        case 'H':
          retree_window(left);
          break;

        case 'J':
          retree_window(downn);
          break;

        case 'K':
          retree_window(upp);
          break;

        case 'L':
          retree_window(right);
          break;

        case '?':
          retree_help();
          break;

        case '+':
          if (treesread <numtrees) {
            readnext = true;
            done = true;
          } else {
            printf("No more trees to read in input file.\n");
          }
          break;

        case 'X':
          done = true;
          break;

        case 'Q':
          done = true;
          break;
      }
    } else {
      printf("Not a possible option!\n");
      maxinput++;
      if (maxinput == 100) {
        printf("ERROR: too many tries at choosing option\n");
        exxit(-1);
      }
    }
  } while (!done);
  if (!written) {
    maxinput = 1;
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
      if (ch == 'Y' || ch == 'y')
        treewrite(&done);
      maxinput++;
      if (maxinput == 100) {
        printf("ERROR: too many tries at choosing option\n");
        exxit(-1);
      }
    } while (ch != 'Y' && ch != 'y' && ch != 'N' && ch != 'n');
  }
}  /* redisplay */


void treeconstruct()
{
  /* constructs a tree from the pointers in nodep. */

  waswritten = false;
  readnext = false;
  do {
    whichtree = 0;
    othertree = 1;
    treesets[whichtree].initialized = false;
    treesets[othertree].initialized = false;
    /*xx    nodep = treesets[whichtree].nodep;   debug */
    root = treesets[whichtree].root;
    subtree = false;
    topedge = 1;
    leftedge = 1;
    buildtree();
    readnext = false;
    written = false;
    waswritten = false;
    printree();
    redisplay();
    if (readnext)
      unbuildtree();
  } while (readnext);
}  /* treeconstruct */


int main(int argc, Char *argv[])
{
  /* reads in spp. Then calls treeconstruct to                              *
     construct the tree and query the user                                  */
  int i;

#ifdef MAC
  argc = 1;                /* macsetup("Retree","");        */
  argv[0] = "Retree";
#endif
  init(argc,argv);
  /*  treesets[0].nodep = treeone;
      treesets[1].nodep = treetwo;*/
  nexus     = false;
  nolengths = false;
  grbg = NULL;
  scrollinc = 20;
  screenlines = 24;
  screenwidth = 80;
  vscreenwidth = 80;
  ibmpc = IBMCRT;
  ansi  = ANSICRT;
  for(i = 0; i < maxsz; i++)
    delarray[i] = false;
  treenumber = 1;
  getoptions();
  configure();
  treeconstruct();
  if (waswritten) {
    outtree = fopen(outtreename, "a");
    if (xmltree)
      fprintf(outtree, "</phylogenies>\n");
    else if (nexus)
      fprintf(outtree, "END;\n");
  }
  FClose(intree);
  FClose(outtree);
#ifdef MAC
  fixmacfile(outtreename);
#endif
#ifdef WIN32
  phyRestoreConsoleAttributes();
#endif
  return 0;
}  /* Retree */

