
/* version 3.6. (c) Copyright 1993-2000 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

/*
    dist.h: included in fitch, kitsch, & neighbor
*/

#define over            60


typedef long *intvector;

typedef node **pointptr;

#ifndef OLDC
/*function prototypes*/
void alloctree(pointptr *, long);
void freetree(pointptr *, long);
void allocd(long, pointptr);
void freed(long, pointptr);
void allocw(long, pointptr);
void freew(long, pointptr);
void setuptree(tree *, long);
void inputdata(boolean, boolean, boolean, boolean, vector *, intvector *);
void coordinates(node *, double, long *, double *, node *, boolean);
void drawline(long, double, node *, boolean);
void printree(node *, boolean, boolean, boolean);
void treeoutr(node *, long *, tree *);
void treeout(node *, long *, double, boolean, node *);
/*function prototypes*/
#endif

