
#include "phylip.h"

#ifdef X
#include <X11/Xlib.h>
#include <X11/Intrinsic.h>
#include <X11/StringDefs.h>
#include <X11/Xaw/Label.h>
#include <X11/Xatom.h>
#include <X11/Shell.h>
#include <X11/Xaw/SimpleMenu.h>
#include <X11/Xaw/Paned.h>
#include <X11/Xaw/Box.h>
#include <X11/Xaw/MenuButton.h>
#include <X11/Xaw/SmeBSB.h>
#include <X11/Xaw/Form.h>
#include <X11/Xaw/Dialog.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#endif

#ifdef MAC
#include "interface.h"
#endif

#include "math.h"  
#define  DEFAULT_STRIPE_HEIGHT 20

#define minus           '-'
#define stripewidth     3000L
#define maxstripedepth  3500
#define fontsize        3800
#define pi              3.1415926535897932384626433
#define ebcdic          EBCDIC
#define segments        40
#define xstart          10
#define ystart          35
#define LF              10
#define CR              13
#define escape  (ebcdic ?  '\'' :  '\033')
#define null  '\000'
#define AFMDIR "/usr/lib/transcript/" /* note trailing slash */

typedef unsigned char byte;
/*typedef char byte; */
typedef enum {treepen, labelpen} pentype;
typedef enum {lw,hp,tek,ibm,mac,houston,decregis,epson,oki,fig,
                citoh,toshiba,pcx,pcl,pict,ray,pov,xpreview,xbm,bmp,
                gif,idraw,vrml,winpreview,other} plottertype;
typedef enum {vertical, horizontal} growth;
typedef enum {cladogram,phenogram,curvogram,
              eurogram,swoopogram,circular} treestyle;
typedef enum {penup,pendown} pensttstype;
typedef enum {plotnow, changeparms, quitnow} winactiontype;
typedef short fonttype[fontsize];
typedef Char *striparray;
typedef striparray striptype[maxstripedepth];

struct LOC_plottext {              /* Local variables for plottext: */
  double height, compress;
  short *font;
  short coord;
  double heightfont, xfactor, yfactor, xfont, yfont, xplot, yplot, sinslope,
         cosslope, xx, yy;
  pensttstype penstatus;
} ;

typedef struct colortype {
  const char *name;
  double red, green, blue;
} colortype;

typedef struct vrmllighttype {
  double intensity, x, y, z;
} vrmllighttype;

double lengthtext(char *, long, char *, fonttype);
double heighttext(fonttype, char *);
void plotrparms(long ntips);
void   clearit(void);
void   getplotter(void);
void   getpreview(void);
const char *figfontname(int id);
boolean isfigfont(char *);
void   plot(pensttstype, double, double);
void   curvespline(double, double, double, double, boolean, long);
void swoopspline(double x1, double y1, double x2, double y2, double x3,
                 double y3, boolean sense, long segs);
void changepen(pentype pen);
void plottext(Char *pstring,long nchars,double height_,double cmpress2,
               double x,double y,double slope,short *font_,char *fontname);
void loadfont(short *font, char *application);
long allocstripe(striptype stripe, long x, long y);
boolean plotpreview(char *, double *, double *, double *, long , node *);
void initplotter(long ntips, char *fontname);
void drawit(char *fontname, double *xoffset, double *yoffset,
                        long numlines, node *root);
void   finishplotter(void);
void   write_bmp_header(FILE *, int, int);
void   turn_rows(byte *, int, int);
void   write_full_pic(byte *, int);
void translate_stripe_to_bmp(striptype *stripe, byte *full_pic,
                           int increment, int width, int div, int *total_bytes);
void   plottree(node *, node *);
void plotlabels(char *fontname);
void   pout(long);
double computeAngle(double oldx, double oldy, double newx, double newy);


/* For povray, added by Dan F. */
#define TREE_TEXTURE "T_Tree\0"
#define NAME_TEXTURE "T_Name\0"



#ifdef X
Display *display;               /* the X display */
extern Window  mainwin;         /* the main display window */
int x, y;                       /* the corner of the window */
unsigned int width, height;     /* the width and height of the window */
#define FONT "-*-new century schoolbook-medium-r-*-*-14-*"
char *fontrsc;                  /* the font resource */
XFontStruct *fontst;            /* the font strcture for the font */
XGCValues gcv;                  /* graphics context values */
GC gc1;                         /* a graphics context */
XtAppContext appcontext;
Widget toplevel;
int nargc;
char** nargv;
extern String res[16];

#define DEFGEOMETRY "600x400+20+50"
#endif
#define LARGE_BUF_LENGTH 500
extern char fontname[LARGE_BUF_LENGTH]; /* the font name to use */

#ifdef WIN32
#define DEFPLOTTER lw
#define DEFPREV winpreview
#endif
#ifdef DOS
#define DEFPLOTTER lw
#define DEFPREV ibm
#endif   
#ifdef MAC 
#ifdef OSX_CARBON
#define DEFPLOTTER lw
#define DEFPREV mac
#else
#define DEFPLOTTER pict
#define DEFPREV mac 
#endif
#endif   
#ifdef VMS
#define DEFPLOTTER lw
#define DEFPREV decregis
#endif
#ifndef DOS
#ifndef MAC 
#ifndef VMS
#ifndef WIN32
#define DEFPLOTTER lw
#ifdef X
#define DEFPREV xpreview
#endif
#ifndef X
#define DEFPREV tek
#endif
#endif
#endif   
#endif   
#endif   

/* Define SEEK_SET (needed for fseek()) for machines that haven't 
   got it already, */
#ifndef SEEK_SET
#define SEEK_SET 0
#endif

