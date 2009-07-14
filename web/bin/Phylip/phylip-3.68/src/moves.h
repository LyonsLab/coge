
/*
  moves.h: included in dnamove, move, dolmove, & retree
*/

typedef enum {
  left, downn, upp, right
} adjwindow;

#ifndef OLDC
/* function prototypes */
void   inpnum(long *, boolean *);
void   prereverse(boolean);
void   postreverse(boolean);
void   chwrite(Char, long, long *, long, long);
void   nnwrite(long, long, long *, long, long);
void   stwrite(const char *,long,long *,long,long);
void   help(const char *);
void   treeoptions(boolean, Char *, FILE **, Char *, Char *);
void   window(adjwindow, long *, long *, long, long, long, long, long,
                long, boolean);
void   pregraph(boolean);

void   pregraph2(boolean);
void   postgraph(boolean);
void   postgraph2(boolean);
void   nextinc(long *, long *, long *, long, long, boolean *, steptr,
                steptr);
void   nextchar(long *, long *, long *, long, long, boolean *);
void   prevchar(long *, long *, long *, long, long, boolean *);
void   show(long *, long *, long *, long, long, boolean *);
/* function prototypes */
#endif

