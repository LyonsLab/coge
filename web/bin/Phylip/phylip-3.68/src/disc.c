#include "phylip.h"
#include "disc.h"

/* version 3.6. (c) Copyright 1993-2004 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

long chars, nonodes, nextree, which;
/*  nonodes = number of nodes in tree                                        *
 *  chars = number of binary characters                                      *
 *  words = number of words needed to represent characters of one organism   */
steptr weight, extras;
boolean printdata;


void inputdata(pointptr treenode,boolean dollo,boolean printdata,FILE *outfile)
{
  /* input the names and character state data for species */
  /* used in Dollop, Dolpenny, Dolmove, & Move */
  long i, j, l;
  char k;
  Char charstate;
  /* possible states are '0', '1', 'P', 'B', and '?' */

  if (printdata)
    headings(chars, "Characters", "----------");
  for (i = 0; i < (chars); i++)
    extras[i] = 0;
  for (i = 1; i <= spp; i++) {
    initname(i-1);
    if (printdata) {
      for (j = 0; j < nmlngth; j++)
        putc(nayme[i - 1][j], outfile);
      fprintf(outfile, "   ");
    }
    for (j = 0; j < (words); j++) {
      treenode[i - 1]->stateone[j] = 0;
      treenode[i - 1]->statezero[j] = 0;
    }
    for (j = 1; j <= (chars); j++) {
      k = (j - 1) % bits + 1;
      l = (j - 1) / bits + 1;
      do {
        if (eoln(infile)) 
          scan_eoln(infile);
        charstate = gettc(infile);
      } while (charstate == ' ' || charstate == '\t');
      if (charstate == 'b')
        charstate = 'B';
      if (charstate == 'p')
        charstate = 'P';
      if (charstate != '0' && charstate != '1' && charstate != '?' &&
          charstate != 'P' && charstate != 'B') {
        printf("\n\nERROR: Bad character state: %c ",charstate);
        printf("at character %ld of species %ld\n\n", j, i);
        exxit(-1);
      }
      if (printdata) {
        newline(outfile, j, 55, nmlngth + 3);
        putc(charstate, outfile);
        if (j % 5 == 0)
          putc(' ', outfile);
      }
      if (charstate == '1')
        treenode[i - 1]->stateone[l - 1] =
          ((long)treenode[i - 1]->stateone[l - 1]) | (1L << k);
      if (charstate == '0')
        treenode[i - 1]->statezero[l - 1] =
          ((long)treenode[i - 1]->statezero[l - 1]) | (1L << k);
      if (charstate == 'P' || charstate == 'B') {
        if (dollo)
          extras[j - 1] += weight[j - 1];
        else {
          treenode[i - 1]->stateone[l - 1] =
            ((long)treenode[i - 1]->stateone[l - 1]) | (1L << k);
          treenode[i - 1]->statezero[l - 1] =
            ((long)treenode[i - 1]->statezero[l - 1]) | (1L << k);
        }
      }
    }
    scan_eoln(infile);
    if (printdata)
      putc('\n', outfile);
  }
  if (printdata)
    fprintf(outfile, "\n\n");
}  /* inputdata */


void inputdata2(pointptr2 treenode)
{
  /* input the names and character state data for species */
  /* used in Mix & Penny */
  long i, j, l;
  char k;
  Char charstate;
  /* possible states are '0', '1', 'P', 'B', and '?' */

  if (printdata)
    headings(chars, "Characters", "----------");
  for (i = 0; i < (chars); i++)
    extras[i] = 0;
  for (i = 1; i <= spp; i++) {
    initname(i-1);
    if (printdata) {
      for (j = 0; j < nmlngth; j++)
        putc(nayme[i - 1][j], outfile);
    }
    fprintf(outfile, "   ");
    for (j = 0; j < (words); j++) {
      treenode[i - 1]->fulstte1[j] = 0;
      treenode[i - 1]->fulstte0[j] = 0;
      treenode[i - 1]->empstte1[j] = 0;
      treenode[i - 1]->empstte0[j] = 0;
    }
    for (j = 1; j <= (chars); j++) {
      k = (j - 1) % bits + 1;
      l = (j - 1) / bits + 1;
      do {
        if (eoln(infile)) 
          scan_eoln(infile);
        charstate = gettc(infile);
        if (charstate == '\n' || charstate == '\t')
          charstate = ' ';
      } while (charstate == ' ');
      if (charstate == 'b')          charstate = 'B';
      if (charstate == 'p')          charstate = 'P';
      if (charstate != '0' && charstate != '1' && charstate != '?' &&
          charstate != 'P' && charstate != 'B') {
        printf("\n\nERROR: Bad character state: %c ",charstate);
        printf("at character %ld of species %ld\n\n", j, i);
        exxit(-1);
      }
      if (printdata) {
        newline(outfile, j, 55, nmlngth + 3);
        putc(charstate, outfile);
        if (j % 5 == 0)
          putc(' ', outfile);
      }
      if (charstate == '1') {
        treenode[i-1]->fulstte1[l-1] =
          ((long)treenode[i-1]->fulstte1[l-1]) | (1L << k);
        treenode[i-1]->empstte1[l-1] =
          treenode[i-1]->fulstte1[l-1];
      }
      if (charstate == '0') {
        treenode[i-1]->fulstte0[l-1] =
          ((long)treenode[i-1]->fulstte0[l-1]) | (1L << k);
        treenode[i-1]->empstte0[l-1] =
          treenode[i-1]->fulstte0[l-1];
      }
      if (charstate == 'P' || charstate == 'B')
        extras[j-1] += weight[j-1];
    }
    scan_eoln(infile);
    if (printdata)
      putc('\n', outfile);
  }
  fprintf(outfile, "\n\n");
}  /* inputdata2 */


void alloctree(pointptr *treenode)
{
  /* allocate tree nodes dynamically */
  /* used in dollop, dolmove, dolpenny, & move */
  long i, j;
  node *p, *q;

  (*treenode) = (pointptr)Malloc(nonodes*sizeof(node *));
  for (i = 0; i < (spp); i++) {
    (*treenode)[i] = (node *)Malloc(sizeof(node));
    (*treenode)[i]->stateone = (bitptr)Malloc(words*sizeof(long));
    (*treenode)[i]->statezero = (bitptr)Malloc(words*sizeof(long));
  }
  for (i = spp; i < (nonodes); i++) {
    q = NULL;
    for (j = 1; j <= 3; j++) {
      p = (node *)Malloc(sizeof(node));
      p->stateone = (bitptr)Malloc(words*sizeof(long));
      p->statezero = (bitptr)Malloc(words*sizeof(long));
      p->next = q;
      q = p;
    }
    p->next->next->next = p;
    (*treenode)[i] = p;
  }
}  /* alloctree */


void alloctree2(pointptr2 *treenode)
{
  /* allocate tree nodes dynamically */
  /* used in mix & penny */
  long i, j;
  node2 *p, *q;

  (*treenode) = (pointptr2)Malloc(nonodes*sizeof(node2 *));
  for (i = 0; i < (spp); i++) {
    (*treenode)[i] = (node2 *)Malloc(sizeof(node2));
    (*treenode)[i]->fulstte1 = (bitptr)Malloc(words*sizeof(long));
    (*treenode)[i]->fulstte0 = (bitptr)Malloc(words*sizeof(long));
    (*treenode)[i]->empstte1 = (bitptr)Malloc(words*sizeof(long));
    (*treenode)[i]->empstte0 = (bitptr)Malloc(words*sizeof(long));
    (*treenode)[i]->fulsteps = (bitptr)Malloc(words*sizeof(long));
    (*treenode)[i]->empsteps = (bitptr)Malloc(words*sizeof(long));
  }
  for (i = spp; i < (nonodes); i++) {
    q = NULL;
    for (j = 1; j <= 3; j++) {
      p = (node2 *)Malloc(sizeof(node2));
      p->fulstte1 = (bitptr)Malloc(words*sizeof(long));
      p->fulstte0 = (bitptr)Malloc(words*sizeof(long));
      p->empstte1 = (bitptr)Malloc(words*sizeof(long));
      p->empstte0 = (bitptr)Malloc(words*sizeof(long));
      p->fulsteps = (bitptr)Malloc(words*sizeof(long));
      p->empsteps = (bitptr)Malloc(words*sizeof(long));
      p->next = q;
      q = p;
    }
    p->next->next->next = p;
    (*treenode)[i] = p;
  }
} /* alloctree2 */


void setuptree(pointptr treenode)
{
  /* initialize tree nodes */
  /* used in dollop, dolmove, dolpenny, & move */
  long i;
  node *p;

  for (i = 1; i <= (nonodes); i++) {
    treenode[i-1]->back = NULL;
    treenode[i-1]->tip = (i <= spp);
    treenode[i-1]->index = i;
    if (i > spp) {
      p = treenode[i-1]->next;
      while (p != treenode[i-1]) {
        p->back = NULL;
        p->tip = false;
        p->index = i;
        p = p->next;
      }
    }
  }
} /* setuptree */


void setuptree2(pointptr2 treenode)
{
  /* initialize tree nodes */
  /* used in mix & penny */
  long i;
  node2 *p;

  for (i = 1; i <= (nonodes); i++) {
    treenode[i-1]->back = NULL;
    treenode[i-1]->tip = (i <= spp);
    treenode[i-1]->index = i;
    if (i > spp) {
      p = treenode[i-1]->next;
      while (p != treenode[i-1]) {
        p->back = NULL;
        p->tip = false;
        p->index = i;
        p = p->next;
      }
    }
  }
} /* setuptree2 */


void inputancestors(boolean *anczero0, boolean *ancone0)
{
  /* reads the ancestral states for each character */
  /* used in dollop, dolmove, dolpenny, mix, move, & penny */
  long i;
  Char ch;

  for (i = 0; i < (chars); i++) {
    anczero0[i] = true;
    ancone0[i] = true;
    do {
      if (eoln(ancfile))
        scan_eoln(ancfile);
      ch = gettc(ancfile);
      if (ch == '\n')
        ch = ' ';
    } while (ch == ' ');
    if (ch == 'p')
      ch = 'P';
    if (ch == 'b')
      ch = 'B';
    if (strchr("10PB?",ch) != NULL){
      anczero0[i] = (ch == '1') ? false : anczero0[i];
      ancone0[i] = (ch == '0') ? false : ancone0[i];
    } else {
      printf("BAD ANCESTOR STATE: %cAT CHARACTER %4ld\n", ch, i + 1);
      exxit(-1);
    }
  }
  scan_eoln(ancfile);
}  /* inputancestorsnew */


void printancestors(FILE *filename, boolean *anczero, boolean *ancone)
{
  /* print out list of ancestral states */
  /* used in dollop, dolmove, dolpenny, mix, move, & penny */
  long i;

  fprintf(filename, "    Ancestral states:\n");
  for (i = 1; i <= nmlngth + 3; i++)
    putc(' ', filename);
  for (i = 1; i <= (chars); i++) {
    newline(filename, i, 55, nmlngth + 3);
    if (ancone[i-1] && anczero[i-1])
      putc('?', filename);
    else if (ancone[i-1])
      putc('1', filename);
    else
      putc('0', filename);
    if (i % 5 == 0)
      putc(' ', filename);
  }
  fprintf(filename, "\n\n");
}  /* printancestor */


void add(node *below, node *newtip, node *newfork, node **root,
                        pointptr treenode)
{
  /* inserts the nodes newfork and its left descendant, newtip,
     to the tree.  below becomes newfork's right descendant.
     The global variable root is also updated */
  /* used in dollop & dolpenny */

  if (below != treenode[below->index - 1])
    below = treenode[below->index - 1];
  if (below->back != NULL)
    below->back->back = newfork;
  newfork->back = below->back;
  below->back = newfork->next->next;
  newfork->next->next->back = below;
  newfork->next->back = newtip;
  newtip->back = newfork->next;
  if (*root == below)
    *root = newfork;
}  /* add */


void add2(node *below, node *newtip, node *newfork, node **root,
                        boolean restoring, boolean wasleft, pointptr treenode)
{
  /* inserts the nodes newfork and its left descendant, newtip,
     to the tree.   below becomes newfork's right descendant */
  /* used in move & dolmove */
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
  if (*root == below)
    *root = newfork;
  (*root)->back = NULL;
}  /* add2 */


void add3(node2 *below, node2 *newtip, node2 *newfork, node2 **root,
                        pointptr2 treenode)
{
  /* inserts the nodes newfork and its left descendant, newtip,
     to the tree.  below becomes newfork's right descendant.
     The global variable root is also updated */
  /* used in mix & penny */
  node2 *p;

  if (below != treenode[below->index - 1])
    below = treenode[below->index - 1];
  if (below->back != NULL)
    below->back->back = newfork;
  newfork->back = below->back;
  below->back = newfork->next->next;
  newfork->next->next->back = below;
  newfork->next->back = newtip;
  newtip->back = newfork->next;
  if (*root == below)
    *root = newfork;
  (*root)->back = NULL;
  p = newfork;
  do {
    p->visited = false;
    p = p->back;
    if (p != NULL) p = treenode[p->index - 1];
  } while (p != NULL);
}  /* add3 */


void re_move(node **item, node **fork, node **root, pointptr treenode)
{
  /* removes nodes item and its ancestor, fork, from the tree.
     the new descendant of fork's ancestor is made to be
     fork's second descendant (other than item).  Also
     returns pointers to the deleted nodes, item and fork.
     The global variable root is also updated */
  /* used in dollop & dolpenny */
  node *p, *q;

  if ((*item)->back == NULL) {
    *fork = NULL;
    return;
  }
  *fork = treenode[(*item)->back->index - 1];
  if (*root == *fork) {
    if (*item == (*fork)->next->back)
      *root = (*fork)->next->next->back;
    else
      *root = (*fork)->next->back;
  }
  p = (*item)->back->next->back;
  q = (*item)->back->next->next->back;
  if (p != NULL)
    p->back = q;
  if (q != NULL)
    q->back = p;
  (*fork)->back = NULL;
  p = (*fork)->next;
  while (p != *fork) {
    p->back = NULL;
    p = p->next;
  }
  (*item)->back = NULL;
}  /* re_move */


void re_move2(node **item, node **fork, node **root, boolean *wasleft,
                        pointptr treenode)
{
  /* removes nodes item and its ancestor, fork, from the tree.
     the new descendant of fork's ancestor is made to be
     fork's second descendant (other than item).   Also
     returns pointers to the deleted nodes, item and fork */
  /* used in move & dolmove */
  node *p, *q;

  if ((*item)->back == NULL) {
    *fork = NULL;
    return;
  }
  *fork = treenode[(*item)->back->index - 1];
  if (*item == (*fork)->next->back) {
    if (*root == *fork)
      *root = (*fork)->next->next->back;
    (*wasleft) = true;
  } else {
    if (*root == *fork)
      *root = (*fork)->next->back;
    (*wasleft) = false;
  }
  p = (*item)->back->next->back;
  q = (*item)->back->next->next->back;
  if (p != NULL)
    p->back = q;
  if (q != NULL)
    q->back = p;
  (*fork)->back = NULL;
  p = (*fork)->next;
  while (p != *fork) {
    p->back = NULL;
    p = p->next;
  }
  (*item)->back = NULL;
}  /* re_move2 */


void re_move3(node2 **item, node2 **fork, node2 **root, pointptr2 treenode)
{
  /* removes nodes item and its ancestor, fork, from the tree.
     the new descendant of fork's ancestor is made to be
     fork's second descendant (other than item).  Also
     returns pointers to the deleted nodes, item and fork.
     The global variable *root is also updated */
  /* used in mix & penny */
  node2 *p, *q;

  if ((*item)->back == NULL) {
    *fork = NULL;
    return;
  }
  *fork = treenode[(*item)->back->index - 1];
  if (*root == *fork) {
    if (*item == (*fork)->next->back)
      *root = (*fork)->next->next->back;
    else
      *root = (*fork)->next->back;
  }
  p = (*item)->back->next->back;
  q = (*item)->back->next->next->back;
  if (p != NULL)
    p->back = q;
  if (q != NULL)
    q->back = p;
  q = (*fork)->back;
  (*fork)->back = NULL;
  p = (*fork)->next;
  while (p != *fork) {
    p->back = NULL;
    p = p->next;
  }
  (*item)->back = NULL;
  if (q != NULL) q = treenode[q->index - 1];
  while (q != NULL) {
    q-> visited = false;
    q = q->back;
    if (q != NULL) q = treenode[q->index - 1];
  }
}  /* re_move3 */


void coordinates(node *p, long *tipy, double f, long *fartemp)
{
  /* establishes coordinates of nodes */
  /* used in dollop, dolpenny, dolmove, & move */
  node *q, *first, *last;

  if (p->tip) {
    p->xcoord = 0;
    p->ycoord = *tipy;
    p->ymin = *tipy;
    p->ymax = *tipy;
    *tipy += down;
    return;
  }
  q = p->next;
  do {
    coordinates(q->back, tipy, f, fartemp);
    q = q->next;
  } while (p != q);
  first = p->next->back;
  q = p->next;
  while (q->next != p)
    q = q->next;
  last = q->back;
  p->xcoord = (last->ymax - first->ymin) * f;
  p->ycoord = (first->ycoord + last->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
  if (p->xcoord > *fartemp)
    *fartemp = p->xcoord;
}  /* coordinates */

void coordinates2(node2 *p, long *tipy)
{
  /* establishes coordinates2 of nodes */
  node2 *q, *first, *last;

  if (p->tip) {
    p->xcoord = 0;
    p->ycoord = *tipy;
    p->ymin = *tipy;
    p->ymax = *tipy;
    (*tipy) += down;
    return;
  }
  q = p->next;
  do {
    coordinates2(q->back, tipy);
    q = q->next;
  } while (p != q);
  first = p->next->back;
  q = p->next;
  while (q->next != p)
    q = q->next;
  last = q->back;
  p->xcoord = last->ymax - first->ymin;
  p->ycoord = (first->ycoord + last->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* coordinates2 */


void treeout(node *p, long nextree, long *col, node *root)
{
  /* write out file with representation of final tree */
  /* used in dollop, dolmove, dolpenny, & move */
  long i, n;
  Char c;
  node *q;

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
    *col += n;
  } else {
    q = p->next;
    putc('(', outtree);
    (*col)++;
    while (q != p) {
      treeout(q->back, nextree, col, root);
      q = q->next;
      if (q == p)
        break;
      putc(',', outtree);
      (*col)++;
      if (*col > 65) {
        putc('\n', outtree);
        *col = 0;
      }
    }
    putc(')', outtree);
    (*col)++;
  }
  if (p != root)
    return;
  if (nextree > 2)
    fprintf(outtree, "[%6.4f];\n", 1.0 / (nextree - 1));
  else
    fprintf(outtree, ";\n");
}  /* treeout */

void treeout2(node2 *p, long *col, node2 *root)
{
  /* write out file with representation of final tree */
  /* used in mix & penny */
  long i, n;
  Char c;

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
    *col += n;
  } else {
    putc('(', outtree);
    (*col)++;
    treeout2(p->next->back, col, root);
    putc(',', outtree);
    (*col)++;
    if (*col > 65) {
      putc('\n', outtree);
      *col = 0;
    }
    treeout2(p->next->next->back, col, root);
    putc(')', outtree);
    (*col)++;
  }
  if (p != root)
    return;
  if (nextree > 2)
    fprintf(outtree, "[%6.4f];\n", 1.0 / (nextree - 1));
  else
    fprintf(outtree, ";\n");
}  /* treeout2 */


void standev(long numtrees, long minwhich, double minsteps,
                        double *nsteps, double **fsteps, longer seed)
{  /* paired sites tests (KHT or SH) on user-defined trees */
   /* used in pars */
  long i, j, k;
  double wt, sumw, sum, sum2, sd;
  double temp;
  double **covar, *P, *f, *r;

#define SAMPLES 1000
  if (numtrees > maxuser) {
    printf("TOO MANY USER-DEFINED TREES");
    printf("  test only performed in the first %ld of them\n", (long)maxuser);
  } else
  if (numtrees == 2) {
    fprintf(outfile, "Kishino-Hasegawa-Templeton test\n\n");
    fprintf(outfile, "Tree    Steps   Diff Steps   Its S.D.");
    fprintf(outfile, "   Significantly worse?\n\n");
    which = 1;
    while (which <= numtrees) {
      fprintf(outfile, "%3ld%10.1f", which, nsteps[which - 1]);
      if (minwhich == which)
        fprintf(outfile, "  <------ best\n");
      else {
        sumw = 0.0;
        sum = 0.0;
        sum2 = 0.0;
        for (i = 0; i < chars; i++) {
          if (weight[i] > 0) {
            wt = weight[i];
            sumw += wt;
            temp = (fsteps[which - 1][i] - fsteps[minwhich - 1][i]) / 10.0;
            sum += wt*temp;
            sum2 += wt * temp * temp;
          }
        }
        temp = sum / sumw;
        sd = sqrt(sumw / (sumw - 1.0) * (sum2 - temp * temp));
        fprintf(outfile, "%10.1f%12.4f",
                (nsteps[which - 1] - minsteps) / 10, sd);
        if (sum > 1.95996 * sd)
          fprintf(outfile, "           Yes\n");
        else
          fprintf(outfile, "           No\n");
      }
      which++;
    }
    fprintf(outfile, "\n\n");
  } else {           /* Shimodaira-Hasegawa test using normal approximation */
    if(numtrees > MAXSHIMOTREES){
      fprintf(outfile, "Shimodaira-Hasegawa test on first %d of %ld trees\n\n"
              , MAXSHIMOTREES, numtrees);
      numtrees = MAXSHIMOTREES;
    } else {
      fprintf(outfile, "Shimodaira-Hasegawa test\n\n");
    }
    covar = (double **)Malloc(numtrees*sizeof(double *));  
    for (i = 0; i < numtrees; i++)
      covar[i] = (double *)Malloc(numtrees*sizeof(double));  
    sumw = 0.0;
    for (i = 0; i < chars; i++)
      sumw += weight[i];
    for (i = 0; i < numtrees; i++) {        /* compute covariances of trees */
      sum = nsteps[i]/(10.0*sumw);
      for (j = 0; j <=i; j++) {
        sum2 = nsteps[j]/(10.0*sumw);
        temp = 0.0;
        for (k = 0; k < chars; k++) {
          if (weight[k] > 0)
            temp = temp + weight[k]*(fsteps[i][k]/10.0-sum)
                                   *(fsteps[j][k]/10.0-sum2);
        }
        covar[i][j] = temp;
        if (i != j)
          covar[j][i] = temp;
      }
    }
    for (i = 0; i < numtrees; i++) { /* in-place Cholesky decomposition
                                        of trees x trees covariance matrix */
      sum = 0.0;
      for (j = 0; j <= i-1; j++)
        sum = sum + covar[i][j] * covar[i][j];
      if (covar[i][i]-sum <= 0.0)
        temp = 0.0;
      else 
        temp = sqrt(covar[i][i] - sum);
      covar[i][i] = temp;
      for (j = i+1; j < numtrees; j++) {
        sum = 0.0;
        for (k = 0; k < i; k++)
          sum = sum + covar[i][k] * covar[j][k];
        if (fabs(temp) < 1.0E-12)
          covar[j][i] = 0.0;
        else 
          covar[j][i] = (covar[j][i] - sum)/temp;
      }
    }
    f = (double *)Malloc(numtrees*sizeof(double)); /* resampled sums */
    P = (double *)Malloc(numtrees*sizeof(double)); /* vector of P's of trees */
    r = (double *)Malloc(numtrees*sizeof(double)); /* store Normal variates */
    for (i = 0; i < numtrees; i++)
      P[i] = 0.0;
    sum2 = nsteps[0];               /* sum2 will be smallest # of steps */
    for (i = 1; i < numtrees; i++)
      if (sum2 > nsteps[i])
        sum2 = nsteps[i];
    for (i = 1; i <= SAMPLES; i++) {          /* loop over resampled trees */
      for (j = 0; j < numtrees; j++)          /* draw Normal variates */
        r[j] = normrand(seed);
      for (j = 0; j < numtrees; j++) {        /* compute vectors */
        sum = 0.0;
        for (k = 0; k <= j; k++)
          sum += covar[j][k]*r[k];
        f[j] = sum;
      }
      sum = f[1];
      for (j = 1; j < numtrees; j++)          /* get min of vector */
        if (f[j] < sum)
          sum = f[j];
      for (j = 0; j < numtrees; j++)          /* accumulate P's */
        if (nsteps[j]-sum2 < f[j] - sum)
          P[j] += 1.0/SAMPLES;
    }
    fprintf(outfile, "Tree    Steps   Diff Steps   P value");
    fprintf(outfile, "   Significantly worse?\n\n");
    for (i = 0; i < numtrees; i++) {
      fprintf(outfile, "%3ld%10.1f", i+1, nsteps[i]);
      if ((minwhich-1) == i)
        fprintf(outfile, "  <------ best\n");
      else {
        fprintf(outfile, "  %9.1f %10.3f", nsteps[i]-sum2, P[i]);
        if (P[i] < 0.05)
          fprintf(outfile, "           Yes\n");
        else
          fprintf(outfile, "           No\n");
      }
    }
  fprintf(outfile, "\n");
  free(P);             /* free the variables we Malloc'ed */
  free(f);
  free(r);
  for (i = 0; i < numtrees; i++)
    free(covar[i]);
  free(covar);
  }
}  /* standev */


void guesstates(Char *guess)
{ /* write best guesses of ancestral states */
  /* used in dollop, dolpenny, mix, & penny */
  long i, j;

  fprintf(outfile, "best guesses of ancestral states:\n");
  fprintf(outfile, "      ");
  for (i = 0; i <= 9; i++)
    fprintf(outfile, "%2ld", i);
  fprintf(outfile, "\n     *--------------------\n");
  for (i = 0; i <= (chars / 10); i++) {
    fprintf(outfile, "%5ld!", i * 10);
    for (j = 0; j <= 9; j++) {
      if (i * 10 + j == 0 || i * 10 + j > chars)
        fprintf(outfile, "  ");
      else
        fprintf(outfile, " %c", guess[i * 10 + j - 1]);
    }
    putc('\n', outfile);
  }
  putc('\n', outfile);
}  /* guesstates */


void freegarbage(gbit **garbage)
{
  /* used in dollop, dolpenny, mix, & penny */
  gbit *p;

  while (*garbage) {
    p = *garbage;
    *garbage = (*garbage)->next;
    free(p->bits_);
    free(p);
  }
}  /* freegarbage */


void disc_gnu(gbit **p, gbit **grbg)
{
  /* this is a do-it-yourself garbage collectors for move
     Make a new node or pull one off the garbage list */

  if (*grbg != NULL) {
    *p = *grbg;
    *grbg = (*grbg)->next;
  } else {
    *p = (gbit *)Malloc(sizeof(gbit));
    (*p)->bits_ = (bitptr)Malloc(words*sizeof(long));
  }
  (*p)->next       = NULL;
}  /* disc_gnu */


void disc_chuck(gbit *p, gbit **grbg)
{
  /* collect garbage on p -- put it on front of garbage list */
  p->next = *grbg;
  *grbg = p;
}  /* disc_chuck */

