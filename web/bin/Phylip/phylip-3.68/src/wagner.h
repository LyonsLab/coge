
/* version 3.6. (c) Copyright 1993-2000 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

/*
  wagner.h: included in move, mix & penny 
*/

#ifndef OLDC
/* function prototypes */
void   inputmixture(bitptr);
void   inputmixturenew(bitptr);
void   printmixture(FILE *, bitptr);
void   fillin(node2 *,long, boolean, bitptr, bitptr);
void   count(long *, bitptr, steptr, steptr);
void   postorder(node2 *, long, boolean, bitptr, bitptr);
void   cpostorder(node2 *, boolean, bitptr, steptr, steptr);
void   filltrav(node2 *, long, boolean, bitptr, bitptr);
void   hyprint(struct htrav_vars2 *,boolean,boolean,boolean,bitptr,Char *);
void   hyptrav(node2 *, boolean, bitptr, long, boolean, boolean, bitptr,
                bitptr, bitptr, pointptr2, Char *, gbit *);
void   hypstates(long, boolean, boolean, boolean, node2 *, bitptr, bitptr,
                bitptr, pointptr2, Char *, gbit *);

void   drawline(long, double, node2 *);
void   printree(boolean, boolean, boolean, node2 *);
void   writesteps(boolean, steptr);
/* function prototypes */
#endif

