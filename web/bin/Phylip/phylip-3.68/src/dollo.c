#include "phylip.h"
#include "disc.h"
#include "dollo.h"

/* version 3.6. (c) Copyright 1993-2004 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

void correct(node *p, long fullset, boolean dollo,
                        bitptr zeroanc, pointptr treenode)
{  /* get final states for intermediate nodes */
  long i;
  long z0, z1, s0, s1, temp;

  if (p->tip)
    return;
  for (i = 0; i < (words); i++) {
    if (p->back == NULL) {
      s0 = zeroanc[i];
      s1 = fullset & (~zeroanc[i]);
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
      z1 &= ~(fullset & (~zeroanc[i]) & z0);
      z0 = temp;
    }
    temp = fullset & (~z0) & (~z1);
    p->statezero[i] = z0 | (temp & s0 & (~s1));
    p->stateone[i] = z1 | (temp & s1 & (~s0));
  }
}  /* correct */

void fillin(node *p)
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

  for (i = 0; i < words; i++) {
    p->stateone[i] = p->next->back->stateone[i] |
                     p->next->next->back->stateone[i];
    p->statezero[i] = p->next->back->statezero[i] |
                      p->next->next->back->statezero[i];
  }
}  /* fillin */

void postorder(node *p)
{
  /* traverses a binary tree, calling PROCEDURE fillin at a
     node's descendants before calling fillin at the node */
  /* used in dollop, dolmove, & move */
  if (p->tip)
    return;
  postorder(p->next->back);
  postorder(p->next->next->back);
  fillin(p);
}  /* postorder */

void count(long *stps, bitptr zeroanc, steptr numszero, steptr numsone)
{
  /* counts the number of steps in a branch of the tree.
     The program spends much of its time in this PROCEDURE */
  /* used in dolpenny & move */
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

void filltrav(node *r)
{
  /* traverse to fill in interior node states */
  if (r->tip)
    return;
  filltrav(r->next->back);
  filltrav(r->next->next->back);
  fillin(r);
}  /* filltrav */


void hyprint(struct htrav_vars *Hyptrav, boolean *unknown,
                        bitptr dohyp, Char *guess)
{
  /* print out states at node */
  long i, j, k;
  char l;
  boolean dot, a0, a1, s0, s1;

  if (Hyptrav->bottom)
    fprintf(outfile, "root   ");
  else
    fprintf(outfile, "%3ld    ", Hyptrav->r->back->index - spp);
  if (Hyptrav->r->tip) {
    for (i = 0; i < nmlngth; i++)
      putc(nayme[Hyptrav->r->index - 1][i], outfile);
  } else
    fprintf(outfile, "%4ld      ", Hyptrav->r->index - spp);
  if (Hyptrav->nonzero)
    fprintf(outfile, "   yes    ");
  else if (*unknown)
    fprintf(outfile, "    ?     ");
  else
    fprintf(outfile, "   no     ");
  for (j = 1; j <= (chars); j++) {
    newline(outfile, j, 40, nmlngth + 17);
    k = (j - 1) / bits + 1;
    l = (j - 1) % bits + 1;
    dot = (((1L << l) & dohyp[k - 1]) == 0 && guess[j - 1] == '?');
    s0 = (((1L << l) & Hyptrav->r->statezero[k - 1]) != 0);
    s1 = (((1L << l) & Hyptrav->r->stateone[k - 1]) != 0);
    a0 = (((1L << l) & Hyptrav->zerobelow->bits_[k - 1]) != 0);
    a1 = (((1L << l) & Hyptrav->onebelow->bits_[k - 1]) != 0);
    dot = (dot || (a1 == s1 && a0 == s0));
    if (dot)
      putc('.', outfile);
    else {
      if (s0) {
        if (s1)
          putc('P', outfile);
        else
          putc('0', outfile);
      } else if (s1)
        putc('1', outfile);
      else
        putc('?', outfile);
    }
    if (j % 5 == 0)
      putc(' ', outfile);
  }
  putc('\n', outfile);
}  /* hyprint */

void hyptrav(node *r_, boolean *unknown, bitptr dohyp, long fullset,
        boolean dollo, Char *guess, pointptr treenode, gbit *garbage,
        bitptr zeroanc, bitptr oneanc)
{
  /*  compute, print out states at one interior node */
  struct htrav_vars HypVars;
  long i;

  HypVars.r = r_;
  disc_gnu(&HypVars.zerobelow, &garbage);
  disc_gnu(&HypVars.onebelow, &garbage);
  if (!HypVars.r->tip)
    correct(HypVars.r, fullset, dollo, zeroanc, treenode);
  HypVars.bottom = (HypVars.r->back == NULL);
  HypVars.nonzero = false;
  if (HypVars.bottom) {
    memcpy(HypVars.zerobelow->bits_, zeroanc, words*sizeof(long));
    memcpy(HypVars.onebelow->bits_, oneanc, words*sizeof(long));
  } else {
    memcpy(HypVars.zerobelow->bits_,
          treenode[HypVars.r->back->index - 1]->statezero, words*sizeof(long));
    memcpy(HypVars.onebelow->bits_,
          treenode[HypVars.r->back->index - 1]->stateone, words*sizeof(long));
  }
  for (i = 0; i < (words); i++)
    HypVars.nonzero = (HypVars.nonzero ||
                      ((HypVars.r->stateone[i] & HypVars.zerobelow->bits_[i])
                      | (HypVars.r->statezero[i]
                        & HypVars.onebelow->bits_[i])) != 0);
  hyprint(&HypVars,unknown,dohyp, guess);
  if (!HypVars.r->tip) {
    hyptrav(HypVars.r->next->back, unknown,dohyp, fullset, dollo, guess,
              treenode, garbage, zeroanc, oneanc);
    hyptrav(HypVars.r->next->next->back, unknown,dohyp, fullset, dollo,
              guess, treenode, garbage, zeroanc, oneanc);
  }
  disc_chuck(HypVars.zerobelow, &garbage);
  disc_chuck(HypVars.onebelow, &garbage);
}  /* hyptrav */

void hypstates(long fullset, boolean dollo, Char *guess, pointptr treenode,
                        node *root, gbit *garbage, bitptr zeroanc, bitptr oneanc)
{
  /* fill in and describe states at interior nodes */
  /* used in dollop & dolpenny */
  boolean unknown  = false;
  bitptr dohyp;
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
    unknown = (unknown || guess[i] == '?');
  }
  dohyp = (bitptr)Malloc(words*sizeof(long));
  for (i = 0; i < words; i++)
    dohyp[i] = zeroanc[i] | oneanc[i];
  filltrav(root);
  fprintf(outfile, "From    To     Any Steps?");
  fprintf(outfile, "    State at upper node\n");
  fprintf(outfile, "                            ");
  fprintf(outfile, " ( . means same as in the node below it on tree)\n\n");
  hyptrav(root, &unknown,dohyp, fullset, dollo, guess, treenode, garbage,
            zeroanc, oneanc);
  free(dohyp);
}  /* hypstates */


void drawline(long i, double scale, node *root)
{
  /* draws one row of the tree diagram by moving up tree */
  node *p, *q, *r, *first =NULL, *last =NULL;
  long n, j;
  boolean extra, done;

  p = root;
  q = root;
  extra = false;
  if (i == (long)p->ycoord && p == root) {
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
    if ((long)q->ycoord == i && !done) {
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
      if ((long)last->ycoord > i && (long)first->ycoord < i
         && i != (long)p->ycoord) {
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
  if ((long)p->ycoord == i && p->tip) {
    for (j = 0; j < nmlngth; j++)
      putc(nayme[p->index - 1][j], outfile);
  }
  putc('\n', outfile);
}  /* drawline */

void printree(double f, boolean treeprint, node *root)
{
  /* prints out diagram of the tree */
  /* used in dollop & dolpenny */
  long i, tipy, dummy;
  double scale;

  putc('\n', outfile);
  if (!treeprint)
    return;
  putc('\n', outfile);
  tipy = 1;
  dummy = 0;
  coordinates(root, &tipy, f, &dummy);
  scale = 1.5;
  putc('\n', outfile);
  for (i = 1; i <= (tipy - down); i++)
    drawline(i, scale, root);
  putc('\n', outfile);
}  /* printree */


void writesteps(boolean weights, boolean dollo, steptr numsteps)
{ /* write number of steps */
  /* used in dollop & dolpenny */
  long i, j, k;

  if (weights)
    fprintf(outfile, "weighted");
  if (dollo)
    fprintf(outfile, " reversions ");
  else
    fprintf(outfile, " polymorphisms ");
  fprintf(outfile, "in each character:\n");
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

