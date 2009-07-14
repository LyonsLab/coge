
#include "phylip.h"
#include "disc.h"
#include "wagner.h"

/* version 3.6. (c) Copyright 1993-2004 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


void inputmixture(bitptr wagner0)
{
  /* input mixture of methods */
  /* used in mix, move, & penny */
  long i, j, k;
  Char ch;
  boolean wag;

  for (i = 0; i < (words); i++)
    wagner0[i] = 0;
  j = 0;
  k = 1;
  for (i = 1; i <= (chars); i++) {
    do {
      if (eoln(mixfile))
        scan_eoln(mixfile);
      ch = gettc(mixfile);
      if (ch == '\n')
        ch = ' ';
    } while (ch == ' ');
    uppercase(&ch);
    wag = false;
    if (ch == 'W' || ch == '?')
      wag = true;
    else if (ch == 'S' || ch == 'C')
      wag = false;
    else {
      printf("BAD METHOD: %c\n", ch);
      exxit(-1);
    }
    if (wag)
      wagner0[k - 1] = (long)wagner0[k - 1] | (1L << j);
    j++;
    if (j > bits) {
      j = 1;
      k++;
    }
  }
  scan_eoln(mixfile);
}  /* inputmixture */


void printmixture(FILE *filename, bitptr wagner)
{
  /* print out list of parsimony methods */
  /* used in mix, move, & penny */
  long i, k, l;

  fprintf(filename, "Parsimony methods:\n");
  l = 0;
  k = 1;
  for (i = 1; i <= nmlngth + 3; i++)
    putc(' ', filename);
  for (i = 1; i <= (chars); i++) {
    newline(filename, i, 55, nmlngth + 3);
    l++;
    if (l > bits) {
      l = 1;
      k++;
    }
    if (((1L << l) & wagner[k - 1]) != 0)
      putc('W', filename);
    else
      putc('S', filename);
    if (i % 5 == 0)
      putc(' ', filename);
  }
  fprintf(filename, "\n\n");
}  /* printmixture */


void fillin(node2 *p,long fullset,boolean full,bitptr wagner,bitptr zeroanc)
{
  /* Sets up for each node in the tree two statesets.
     stateone and statezero are the sets of character
     states that must be 1 or must be 0, respectively,
     in a most parsimonious reconstruction, based on the
     information at or above this node.  Note that this
     state assignment may change based on information further
     down the tree.  If a character is in both sets it is in
     state "P".  If in neither, it is "?". */
  long i;
  long l0, l1, r0, r1, st, wa, za;

  for (i = 0; i < (words); i++) {
    if (full) {
      l0 = p->next->back->fulstte0[i];
      l1 = p->next->back->fulstte1[i];
      r0 = p->next->next->back->fulstte0[i];
      r1 = p->next->next->back->fulstte1[i];
    }
    else {
      l0 = p->next->back->empstte0[i];
      l1 = p->next->back->empstte1[i];
      r0 = p->next->next->back->empstte0[i];
      r1 = p->next->next->back->empstte1[i];
    }
    st = (l1 & r0) | (l0 & r1);
    wa = wagner[i];
    za = zeroanc[i];
    if (full) {
      p->fulstte1[i] = (l1 | r1) & (~(st & (wa | za)));
      p->fulstte0[i] = (l0 | r0) & (~(st & (wa | (fullset & (~za)))));
      p->fulsteps[i] = st;
    }
    else {
      p->empstte1[i] = (l1 | r1) & (~(st & (wa | za)));
      p->empstte0[i] = (l0 | r0) & (~(st & (wa | (fullset & (~za)))));
      p->empsteps[i] = st;
    }
  }
}  /* fillin */


void count(long *stps, bitptr zeroanc, steptr numszero, steptr numsone)
{
  /* counts the number of steps in a fork of the tree.
     The program spends much of its time in this PROCEDURE */
  /* used in mix & penny */
  long i, j, l;

  j = 1;
  l = 0;
  for (i = 0; i < (chars); i++) {
    l++;
    if (l > bits) {
      l = 1;
      j++;
    }
    if (((1L << l) & stps[j - 1]) != 0) {
      if (((1L << l) & zeroanc[j - 1]) != 0)
        numszero[i] += weight[i];
      else
        numsone[i] += weight[i];
    }
  }
}  /* count */


void postorder(node2 *p, long fullset, boolean full, bitptr wagner,
                        bitptr zeroanc)
{
  /* traverses a binary tree, calling PROCEDURE fillin at a
     node's descendants before calling fillin at the node2 */
  /* used in mix & penny */
  if (p->tip)
    return;
  postorder(p->next->back, fullset, full, wagner, zeroanc);
  postorder(p->next->next->back, fullset, full, wagner, zeroanc);
  if (!p->visited) {
    fillin(p, fullset, full, wagner, zeroanc);
    if (!full) p->visited = true;
  }
}  /* postorder */

void cpostorder(node2 *p, boolean full, bitptr zeroanc, steptr numszero,
                        steptr numsone)
{
  /* traverses a binary tree, calling PROCEDURE count at a
     node's descendants before calling count at the node2 */
  /* used in mix & penny */
  if (p->tip)
    return;
  cpostorder(p->next->back, full, zeroanc, numszero, numsone);
  cpostorder(p->next->next->back, full, zeroanc, numszero, numsone);
  if (full)
    count(p->fulsteps, zeroanc, numszero, numsone);
  else
    count(p->empsteps, zeroanc, numszero, numsone);
}  /* cpostorder */


void filltrav(node2 *r, long fullset, boolean full, bitptr wagner,
                        bitptr zeroanc)
{
  /* traverse to fill in interior node states */
  if (r->tip)
    return;
  filltrav(r->next->back, fullset, full, wagner, zeroanc);
  filltrav(r->next->next->back, fullset, full, wagner, zeroanc);
  fillin(r, fullset, full, wagner, zeroanc);
}  /* filltrav */


void hyprint(struct htrav_vars2 *htrav, boolean unknown, boolean noroot,
                        boolean didreroot, bitptr wagner, Char *guess)
{
  /* print out states at node2 */
  long i, j, k;
  char l;
  boolean dot, a0, a1, s0, s1;

  if (htrav->bottom) {
    if (noroot && !didreroot)
      fprintf(outfile, "       ");
    else
      fprintf(outfile, "root   ");
  } else
    fprintf(outfile, "%3ld    ", htrav->r->back->index - spp);
  if (htrav->r->tip) {
    for (i = 0; i < nmlngth; i++)
      putc(nayme[htrav->r->index - 1][i], outfile);
  } else
    fprintf(outfile, "%4ld      ", htrav->r->index - spp);
  if (htrav->bottom && noroot && !didreroot)
    fprintf(outfile, "          ");
  else if (htrav->nonzero)
    fprintf(outfile, "   yes    ");
  else if (unknown)
    fprintf(outfile, "    ?     ");
  else if (htrav->maybe)
    fprintf(outfile, "  maybe   ");
  else
    fprintf(outfile, "   no     ");
  for (j = 1; j <= (chars); j++) {
    newline(outfile, j, 40, nmlngth + 17);
    k = (j - 1) / bits + 1;
    l = (j - 1) % bits + 1;
    dot = (((1L << l) & wagner[k - 1]) == 0 && guess[j - 1] == '?');
    s0 = (((1L << l) & htrav->r->empstte0[k - 1]) != 0);
    s1 = (((1L << l) & htrav->r->empstte1[k - 1]) != 0);
    a0 = (((1L << l) & htrav->zerobelow->bits_[k - 1]) != 0);
    a1 = (((1L << l) & htrav->onebelow->bits_[k - 1]) != 0);
    dot = (dot || ((!htrav->bottom || !noroot || didreroot) && a1 == s1 &&
                   a0 == s0));
    if (dot)
      putc('.', outfile);
    else {
      if (s0)
        putc('0', outfile);
      else if (s1)
        putc('1', outfile);
      else
        putc('?', outfile);
    }
    if (j % 5 == 0)
      putc(' ', outfile);
  }
  putc('\n', outfile);
}  /* hyprint */


void hyptrav(node2 *r_, boolean unknown, bitptr dohyp, long fullset,
        boolean noroot, boolean didreroot, bitptr wagner,
        bitptr zeroanc, bitptr oneanc, pointptr2 treenode, Char *guess,
        gbit *garbage)
{
  /* compute, print out states at one interior node2 */
  /* used in mix & penny */
  struct htrav_vars2 vars;
  long i;
  long l0, l1, r0, r1, s0, s1, a0, a1, temp, dh, wa;

  vars.r = r_;
  disc_gnu(&vars.zerobelow, &garbage);
  disc_gnu(&vars.onebelow, &garbage);
  vars.bottom = (vars.r->back == NULL);
  vars.maybe = false;
  vars.nonzero = false;
  if (vars.bottom) {
    memcpy(vars.zerobelow->bits_, zeroanc, words*sizeof(long));
    memcpy(vars.onebelow->bits_, oneanc, words*sizeof(long));
  } else {
    memcpy(vars.zerobelow->bits_, treenode[vars.r->back->index - 1]->empstte0,
           words*sizeof(long));
    memcpy(vars.onebelow->bits_, treenode[vars.r->back->index - 1]->empstte1,
           words*sizeof(long));
  }
  for (i = 0; i < (words); i++) {
    dh = dohyp[i];
    s0 = vars.r->empstte0[i];
    s1 = vars.r->empstte1[i];
    a0 = vars.zerobelow->bits_[i];
    a1 = vars.onebelow->bits_[i];
    if (!vars.r->tip) {
      wa = wagner[i];
      l0 = vars.r->next->back->empstte0[i];
      l1 = vars.r->next->back->empstte1[i];
      r0 = vars.r->next->next->back->empstte0[i];
      r1 = vars.r->next->next->back->empstte1[i];
      s0 = (wa & ((a0 & l0) | (a0 & r0) | (l0 & r0))) |
           (dh & fullset & (~wa) & s0);
      s1 = (wa & ((a1 & l1) | (a1 & r1) | (l1 & r1))) |
           (dh & fullset & (~wa) & s1);
      temp = fullset & (~(s0 | s1 | l1 | l0 | r1 | r0));
      s0 |= temp & a0;
      s1 |= temp & a1;
      vars.r->empstte0[i] = s0;
      vars.r->empstte1[i] = s1;
    }
    vars.maybe = (vars.maybe || (dh & (s0 | s1)) != (a0 | a1));
    vars.nonzero = (vars.nonzero || ((s1 & a0) | (s0 & a1)) != 0);
  }
  hyprint(&vars,unknown, noroot, didreroot, wagner, guess);
  if (!vars.r->tip) {                
    hyptrav(vars.r->next->back,unknown,dohyp, fullset, noroot,didreroot,
               wagner, zeroanc, oneanc, treenode, guess, garbage);
    hyptrav(vars.r->next->next->back, unknown,dohyp, fullset, noroot,
               didreroot, wagner, zeroanc, oneanc, treenode, guess, garbage);
  }
  disc_chuck(vars.zerobelow, &garbage);
  disc_chuck(vars.onebelow, &garbage);
}  /* hyptrav */


void hypstates(long fullset, boolean full, boolean noroot, boolean didreroot,
                        node2 *root, bitptr wagner, bitptr zeroanc, bitptr oneanc,
                        pointptr2 treenode, Char *guess, gbit *garbage)
{
  /* fill in and describe states at interior nodes */
  /* used in mix & penny */
  boolean unknown;
  bitptr dohyp;
  long i, j, k;

  for (i = 0; i < (words); i++) {
    zeroanc[i] = 0;
    oneanc[i] = 0;
  }
  unknown = false;
  for (i = 0; i < (chars); i++) {
    j = i / bits + 1;
    k = i % bits + 1;
    if (guess[i] == '0')
      zeroanc[j - 1] = ((long)zeroanc[j - 1]) | (1L << k);
    if (guess[i] == '1')
      oneanc[j - 1] = ((long)oneanc[j - 1]) | (1L << k);
    unknown = (unknown || ((((1L << k) & wagner[j - 1]) == 0) &&
                  guess[i] == '?'));
  }
  dohyp = (bitptr)Malloc(words*sizeof(long));
  for (i = 0; i < (words); i++)
    dohyp[i] = wagner[i] | zeroanc[i] | oneanc[i];
  filltrav(root, fullset, full, wagner, zeroanc);
  fprintf(outfile, "From    To     Any Steps?    ");
  fprintf(outfile, "State at upper node\n");
  fprintf(outfile, "                             ");
  fprintf(outfile, "( . means same as in the node below it on tree)\n\n");
  hyptrav(root,unknown,dohyp, fullset, noroot, didreroot, wagner,
                zeroanc, oneanc, treenode, guess, garbage);
  free(dohyp);
}  /* hypstates */


void drawline(long i, double scale, node2 *root)
{
  /* draws one row of the tree diagram by moving up tree */
  node2 *p, *q, *r, *first =NULL, *last =NULL;
  long n, j;
  boolean extra, done;

  p = root;
  q = root;
  extra = false;
  if (i == p->ycoord && p == root) {
    if (p->index - spp >= 10)
      fprintf(outfile, "-%2ld", p->index - spp);
    else
      fprintf(outfile, "--%ld", p->index - spp);
    extra = true;
  } else
    fprintf(outfile, "  ");
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
    n = (long)(scale * (p->xcoord - q->xcoord) + 0.5);
    if (n < 3 && !q->tip)
      n = 3;
    if (extra) {
      n--;
      extra = false;
    }
    if (q->ycoord == i && !done) {
      putc('+', outfile);
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
      if (last->ycoord > i && first->ycoord < i && i != p->ycoord) {
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
    if (p != q)
      p = q;
  } while (!done);
  if (p->ycoord == i && p->tip) {
    for (j = 0; j < nmlngth; j++)
      putc(nayme[p->index - 1][j], outfile);
  }
  putc('\n', outfile);
}  /* drawline */


void printree(boolean treeprint,boolean noroot,boolean didreroot,node2 *root)
{
  /* prints out diagram of the tree */
  /* used in mix & penny */
  long tipy, i;
  double scale;

  putc('\n', outfile);
  if (!treeprint)
    return;
  putc('\n', outfile);
  tipy = 1;
  coordinates2(root, &tipy);
  scale = 1.5;
  putc('\n', outfile);
  for (i = 1; i <= (tipy - down); i++)
    drawline(i, scale, root);
  if (noroot) {
    fprintf(outfile, "\n  remember:");
    if (didreroot)
      fprintf(outfile, " (although rooted by outgroup)");
    fprintf(outfile, " this is an unrooted tree!\n");
  }
  putc('\n', outfile);
}  /* printree */


void writesteps(boolean weights, steptr numsteps)
{ /* write number of steps */
  /* used in mix & penny */
  long i, j, k;

  if (weights)
    fprintf(outfile, "weighted ");
  fprintf(outfile, "steps in each character:\n");
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
      fprintf(outfile, "%4ld", numsteps[k - 1] + extras[k - 1]);
    }
    putc('\n', outfile);
  }
  putc('\n', outfile);
}  /* writesteps */

