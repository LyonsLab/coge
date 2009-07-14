
/* version 3.6. (c) Copyright 1993-2000 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

/*
    disc.h: included in mix, move, penny, dollop, dolmove, dolpenny,
            & clique
*/


/* node and pointptr used in Dollop, Dolmove, Dolpenny, Move, & Clique */

typedef node **pointptr;

/* node and pointptr used in Mix & Penny */

typedef struct node2 {         /* describes a tip species or an ancestor   */
  struct node2 *next, *back;
  long index;
  boolean tip, bottom,visited;/* present species are tips of tree         */
  bitptr fulstte1, fulstte0;  /* see in PROCEDURE fillin                  */
  bitptr empstte1, empstte0;  /* see in PROCEDURE fillin                  */
  bitptr fulsteps,empsteps;
  long xcoord, ycoord, ymin;  /* used by printree                        */
  long ymax;
} node2;

typedef node2 **pointptr2;

typedef struct gbit {
  bitptr bits_;
  struct gbit *next;
} gbit;

typedef struct htrav_vars {
  node *r;
  boolean bottom, nonzero;
  gbit *zerobelow, *onebelow;
} htrav_vars;

typedef struct htrav_vars2 {
  node2 *r;
  boolean bottom, maybe, nonzero;
  gbit *zerobelow, *onebelow;
} htrav_vars2;


extern long chars, nonodes,  nextree, which;
/*  nonodes = number of nodes in tree                                        *
 *  chars = number of binary characters                                      */
extern steptr weight, extras;
extern boolean printdata;

#ifndef OLDC
/*function prototypes*/
void inputdata(pointptr, boolean, boolean, FILE *);
void inputdata2(pointptr2);
void alloctree(pointptr *);
void alloctree2(pointptr2 *);
void setuptree(pointptr);
void setuptree2(pointptr2);
void inputancestors(boolean *, boolean *);
void inputancestorsnew(boolean *, boolean *);
void printancestors(FILE *, boolean *, boolean *);
void add(node *, node *, node *, node **, pointptr);
void add2(node *, node *, node *, node **, boolean, boolean, pointptr);

void add3(node2 *, node2 *, node2 *, node2 **, pointptr2);
void re_move(node **, node **, node **, pointptr);
void re_move2(node **, node **, node **, boolean *, pointptr);
void re_move3(node2 **, node2 **, node2 **, pointptr2);
void coordinates(node *, long *, double , long *);
void coordinates2(node2 *, long *);
void treeout(node *, long, long *, node *);
void treeout2(node2 *, long *, node2 *);
void standev(long, long, double, double *, double **, longer);
void guesstates(Char *);

void freegarbage(gbit **);
void disc_gnu(gbit **, gbit **);
void disc_chuck(gbit *, gbit **);
/*function prototypes*/
#endif
