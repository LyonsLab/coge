#include "printree.h"

static void mlk_drawline(FILE *fp, tree *t, long i, double scale)
{
  /* draws one row of the tree diagram by moving up tree */
  node *p, *q, *r, *first =NULL, *last =NULL;
  long n, j;
  boolean extra, done;

  p = t->root;
  q = t->root;
  extra = false;
  if ((long)(p->ycoord) == i) {
    if (p->index - spp >= 10)
      fprintf(fp, "-%2ld", p->index - spp);
    else
      fprintf(fp, "--%ld", p->index - spp);
    extra = true;
  } else
    fprintf(fp, "  ");
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
    n = (long)(scale * ((long)(p->xcoord) - (long)(q->xcoord)) + 0.5);
    if (n < 3 && !q->tip)
      n = 3;
    if (extra) {
      n--;
      extra = false;
    }
    if ((long)(q->ycoord) == i && !done) {
      if (p->ycoord != q->ycoord)
        putc('+', fp);
      else
        putc('-', fp);
      if (!q->tip) {
        for (j = 1; j <= n - 2; j++)
          putc('-', fp);
        if (q->index - spp >= 10)
          fprintf(fp, "%2ld", q->index - spp);
        else
          fprintf(fp, "-%ld", q->index - spp);
        extra = true;
      } else {
        for (j = 1; j < n; j++)
          putc('-', fp);
      }
    } else if (!p->tip) {
      if ((long)(last->ycoord) > i && (long)(first->ycoord) < i && 
           i != (long)(p->ycoord)) {
        putc('!', fp);
        for (j = 1; j < n; j++)
          putc(' ', fp);
      } else {
        for (j = 1; j <= n; j++)
          putc(' ', fp);
      }
    } else {
      for (j = 1; j <= n; j++)
        putc(' ', fp);
    }
    if (p != q)
      p = q;
  } while (!done);
  if ((long)(p->ycoord) == i && p->tip) {
    for (j = 0; j < nmlngth; j++)
      putc(nayme[p->index - 1][j], fp);
  }
  putc('\n', fp);
}  /* mlk_drawline */


static void mlk_coordinates(node *p, long *tipy)
{
  /* establishes coordinates of nodes */
  node *q, *first, *last, *pp1 =NULL, *pp2 =NULL;
  long num_sibs, p1, p2, i;

  if (p->tip) {
    p->xcoord = 0;
    p->ycoord = (*tipy);
    p->ymin   = (*tipy);
    p->ymax   = (*tipy);
    (*tipy)  += down;
    return;
  }
  q = p->next;
  do {
    mlk_coordinates(q->back, tipy);
    q = q->next;
  } while (p != q);
  num_sibs = count_sibs(p);
  p1 = (long)((num_sibs+1)/2.0);
  p2 = (long)((num_sibs+2)/2.0);
  i = 1;
  q = p->next;
  first  = q->back;
  do {
    if (i == p1) pp1 = q->back;
    if (i == p2) pp2 = q->back;
    last = q->back;
    q = q->next;
    i++;
  } while (q != p);
  p->xcoord = (long)(0.5 - over * p->tyme);
  p->ycoord = (pp1->ycoord + pp2->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* mlk_coordinates */

void mlk_printree(FILE *fp, tree *t)
{
 /* prints out diagram of the tree */
  long tipy;
  double scale;
  long i;
  node *p;

  assert(fp != NULL);

  putc('\n', fp);
  tipy = 1;
  mlk_coordinates(t->root, &tipy);
  p = t->root;
  while (!p->tip)
    p = p->next->back;
  scale = 1.0 / (long)(p->tyme - t->root->tyme + 1.000);
  putc('\n', fp);
  for (i = 1; i <= tipy - down; i++)
    mlk_drawline(fp, t, i, scale);
  putc('\n', fp);
}  /* dnamlk_printree */


static void describe_r(FILE *fp, tree *t, node *p, double fracchange)
{
  long i, num_sibs;
  node *sib_ptr, *sib_back_ptr;
  double v;

  if (p == t->root)
    fprintf(fp, " root         ");
  else
    fprintf(fp, "%4ld          ", p->back->index - spp);
  if (p->tip) {
    for (i = 0; i < nmlngth; i++)
      putc(nayme[p->index - 1][i], fp);
  } else
    fprintf(fp, "%4ld      ", p->index - spp);
  if (p != t->root) {
    fprintf(fp, "%11.5f", fracchange * (p->tyme - t->root->tyme));
    v = fracchange * (p->tyme - t->nodep[p->back->index - 1]->tyme);
    fprintf(fp, "%13.5f", v);
  }
  putc('\n', fp);
  if (!p->tip) {

    sib_ptr = p;
    num_sibs = count_sibs(p);
    for (i=0 ; i < num_sibs; i++) {
      sib_ptr      = sib_ptr->next;
      sib_back_ptr = sib_ptr->back;
      describe_r(fp, t, sib_back_ptr, fracchange);
    }
  }
}  /* describe */


void mlk_describe(FILE *fp, tree *t, double fracchange)
{
  describe_r(fp, t, t->root, fracchange);
}
  
