#ifdef WIN32
#include <windows.h>
#endif

#ifdef OSX_CARBON
#include <Carbon/Carbon.h>
#include "interface.h"
#endif

#include "phylip.h"
#include "draw.h"

#ifdef QUICKC
struct videoconfig myscreen;
void   setupgraphics();
#endif

#ifdef WIN32
extern HDC hdc;
extern HPEN hPenTree;
extern HPEN hPenLabel;
extern void winplotpreview();

struct winpreviewparms_t {
  char * fn;
  double *xo, *yo, *scale;
  long nt;
  node *root;
};

struct winpreviewparms_t winpreviewparms;

#endif
#ifdef X
struct {
        char* fn;
        double *xo, *yo, *scale;
        long nt;
        node *root;
} xpreviewparms;
  
Atom wm_delete_window;
Atom wm_delete_window2;

Widget dialog;
Widget shell=NULL;
Window mainwin=0;

void init_x(void);
void  redraw(Widget w,XtPointer client, XExposeEvent *ev);
void plot_callback(Widget w,XtPointer client, XtPointer call);
void change_callback(Widget w,XtPointer client, XtPointer call);
void about_callback(Widget w,XtPointer client, XtPointer call);
void quit_callback(Widget w,XtPointer client, XtPointer call);
void close_x(void);
void do_dialog(void);
void delete_callback(Widget w, XEvent* event, String *params, int *num_params); 
void dismiss_dialog(void);
#endif


long winheight;
long winwidth;

extern winactiontype winaction;

colortype colors[7] = {
  {"White    ",1.0,1.0,1.0},
  {"Red      ",1.0,0.3,0.3},
  {"Orange   ",1.0,0.6,0.6},
  {"Yellow   ",1.0,0.9,0.4},
  {"Green    ",0.3,0.8,0.3},
  {"Blue     ",0.5,0.5,1.0},
  {"Violet   ",0.6,0.4,0.8},
};

vrmllighttype vrmllights[3] = {
  {1.0, -100.0, 100.0, 100.0},
  {0.5, 100.0, -100.0, -100.0},
  {0.3, 0.0, -100.0, 100.0},
};

long treecolor, namecolor, backcolor, bottomcolor, vrmlskycolornear, vrmlskycolorfar,
     vrmlgroundcolornear, vrmlgroundcolorfar, vrmlplotcolor;

char fontname[LARGE_BUF_LENGTH];

/* format of matrix: capheight,  length[32],length[33],..length[256]*/

byte *full_pic ;
int increment = 0 ;
int total_bytes = 0 ;

short unknown_metric[256];

static short helvetica_metric[] = { 718,
278,278,355,556,556,889,667,222,333,333,389,584,278,333,278,278,556,556,556,
556,556,556,556,556,556,556,278,278,584,584,584,556,1015,667,667,722,722,667,
611,778,722,278,500,667,556,833,722,778,667,778,722,667,611,722,667,944,667,
667,611,278,278,278,469,556,222,556,556,500,556,556,278,556,556,222,222,500,
222,833,556,556,556,556,333,500,278,556,500,722,500,500,500,334,260,334,584,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,333,556,
556,167,556,556,556,556,191,333,556,333,333,500,500,0,556,556,556,278,0,537,
350,222,333,333,556,1000,1000,0,611,0,333,333,333,333,333,333,333,333,0,333,
333,0,333,333,333,1000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1000,0,370,0,0,0,0,556,
778,1000,365,0,0,0,0,0,889,0,0,0,278,0,0,222,611,944,611,0,0,0};
static short helveticabold_metric[] = {718, /* height */ 
278,333,474,556,556,889,722,278,333,333,389,584,278,333,278,278,556,556,556,
556,556,556,556,556,556,556,333,333,584,584,584,611,975,722,722,722,722,667,
611,778,722,278,556,722,611,833,722,778,667,778,722,667,611,722,667,944,667,
667,611,333,278,333,584,556,278,556,611,556,611,556,333,611,611,278,278,556,
278,889,611,611,611,611,389,556,333,611,556,778,556,556,500,389,280,389,584,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,333,556,
556,167,556,556,556,556,238,500,556,333,333,611,611,0,556,556,556,278,0,556,
350,278,500,500,556,1000,1000,0,611,0,333,333,333,333,333,333,333,333,0,333,
333,0,333,333,333,1000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1000,0,370,0,0,0,0,611,
778,1000,365,0,0,0,0,0,889,0,0,0,278,0,0,278,611,944,611,0,0,0};
static short timesroman_metric[] = {662,
250,333,408,500,500,833,778,333,333,333,500,564,250,333,250,278,500,500,500,
500,500,500,500,500,500,500,278,278,564,564,564,444,921,722,667,667,722,611,
556,722,722,333,389,722,611,889,722,722,556,722,667,556,611,722,722,944,722,
722,611,333,278,333,469,500,333,444,500,444,500,444,333,500,500,278,278,500,
278,778,500,500,500,500,333,389,278,500,500,722,500,500,444,480,200,480,541,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,333,500,
500,167,500,500,500,500,180,444,500,333,333,556,556,0,500,500,500,250,0,453,
350,333,444,444,500,1000,1000,0,444,0,333,333,333,333,333,333,333,333,0,333,
333,0,333,333,333,1000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,889,0,276,0,0,0,0,611,
722,889,310,0,0,0,0,0,667,0,0,0,278,0,0,278,500,722,500,0,0,0};
static short timesitalic_metric[] = {660, /* height */ 
250,333,420,500,500,833,778,333,333,333,500,675,250,333,250,278,500,500,500,
500,500,500,500,500,500,500,333,333,675,675,675,500,920,611,611,667,722,611,
611,722,722,333,444,667,556,833,667,722,611,722,611,500,556,722,611,833,611,
556,556,389,278,389,422,500,333,500,500,444,500,444,278,500,500,278,278,444,
278,722,500,500,500,500,389,389,278,500,444,667,444,444,389,400,275,400,541,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,389,500,
500,167,500,500,500,500,214,556,500,333,333,500,500,0,500,500,500,250,0,523,
350,333,556,556,500,889,1000,0,500,0,333,333,333,333,333,333,333,333,0,333,
333,0,333,333,333,889,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,889,0,276,0,0,0,0,556,
722,944,310,0,0,0,0,0,667,0,0,0,278,0,0,278,500,667,500,0,0,0};
static short timesbold_metric[] = {681, /* height */ 
250,333,555,500,500,1000,833,333,333,333,500,570,250,333,250,278,500,500,500,
500,500,500,500,500,500,500,333,333,570,570,570,500,930,722,667,722,722,667,
611,778,778,389,500,778,667,944,722,778,611,778,722,556,667,722,722,1000,722,
722,667,333,278,333,581,500,333,500,556,444,556,444,333,500,556,278,333,556,
278,833,556,500,556,556,444,389,333,556,500,722,500,500,444,394,220,394,520,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,333,500,500,
167,500,500,500,500,278,500,500,333,333,556,556,0,500,500,500,250,0,540,350,
333,500,500,500,1000,1000,0,500,0,333,333,333,333,333,333,333,333,0,333,333,
0,333,333,333,1000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1000,0,300,0,0,0,0,667,778,
1000,330,0,0,0,0,0,722,0,0,0,278,0,0,278,500,722,556,0,0,0};
static short timesbolditalic_metric[] = {662, /* height */ 
250,389,555,500,500,833,778,333,333,333,500,570,250,333,250,278,500,500,500,
500,500,500,500,500,500,500,333,333,570,570,570,500,832,667,667,667,722,667,
667,722,778,389,500,667,611,889,722,722,611,722,667,556,611,722,667,889,667,
611,611,333,278,333,570,500,333,500,500,444,500,444,333,500,556,278,278,500,
278,778,556,500,500,500,389,389,278,556,444,667,500,444,389,348,220,348,570,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,389,500,
500,167,500,500,500,500,278,500,500,333,333,556,556,0,500,500,500,250,0,500,
350,333,500,500,500,1000,1000,0,500,0,333,333,333,333,333,333,333,333,0,333,
333,0,333,333,333,1000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,944,0,266,0,0,0,0,611
,722,944,300,0,0,0,0,0,722,0,0,0,278,0,0,278,500,722,500,0,0,0};
static const char *figfonts[] = {"Times-Roman","Times-Italic","Times-Bold","Times-BoldItalic",
 "AvantGarde-Book","AvantGarde-BookOblique","AvantGarde-Demi","AvantGarde-DemiOblique",
"Bookman-Light","Bookman-LightItalic","Bookman-Demi","Bookman-DemiItalic",
"Courier","Courier-Italic","Courier-Bold","Courier-BoldItalic",
"Helvetica","Helvetica-Oblique","Helvetica-Bold","Helvetica-BoldOblique",
"Helvetica-Narrow","Helvetica-Narrow-Oblique","Helvetica-Narrow-Bold","Helvetica-Narrow-BoldOblique",
"NewCenturySchlbk-Roman","NewCenturySchlbk-Italic","NewCenturySchlbk-Bold","NewCenturySchlbk-BoldItalic",
"Palatino-Roman","Palatino-Italic","Palatino-Bold","Palatino-BoldItalic",
"Symbol","ZapfChancery-MediumItalic","ZapfDingbats"};                               
                               
double oldx, oldy;

boolean didloadmetric;
long   nmoves,oldpictint,pagecount;
double labelline,linewidth,oldxhigh,oldxlow,oldyhigh,oldylow,
       vrmllinewidth, raylinewidth,treeline,oldxsize,oldysize,oldxunitspercm,
       oldyunitspercm,oldxcorner,oldycorner,oldxmargin,oldymargin,
       oldhpmargin,oldvpmargin,clipx0,clipx1,clipy0,clipy1,userxsize,userysize;
long rootmatrix[51][51];
long  HiMode,GraphDriver,GraphMode,LoMode,bytewrite;

/* externals should move to .h file later. */
extern long          strpbottom,strptop,strpwide,strpdeep,strpdiv,hpresolution;
extern boolean       dotmatrix,empty,preview,previewing,pictbold,pictitalic,
                     pictshadow,pictoutline;
extern double        expand,xcorner,xnow,xsize,xscale,xunitspercm,
                     ycorner,ynow,ysize,yscale,yunitspercm,labelrotation,
                     labelheight,xmargin,ymargin,pagex,pagey,paperx,papery,
                     hpmargin,vpmargin;
extern long          filesize;
extern growth        grows;
extern enum {yes,no} penchange,oldpenchange;
extern FILE          *plotfile;
extern plottertype   plotter,oldplotter,previewer;
extern striptype     stripe;
extern char             resopts;
pentype                     lastpen;
extern char pltfilename[FNMLNGTH];
extern char progname[FNMLNGTH];

#define NO_PLANE 666    /* To make POVRay happy */

#ifndef OLDC
/* function prototypes */
int    pointinrect(double, double, double, double, double, double);
int    rectintersects(double, double, double, double,
                double, double, double, double);
long   upbyte(long);
long   lobyte(long);

void   pictoutint(FILE *, long);
Local long SFactor(void);
long   DigitsInt(long);
Local boolean IsColumnEmpty(striparray *, long, long);
void   Skip(long Amount);
Local long FirstBlack(striparray *, long, long);
Local long FirstWhite(striparray *, long, long);
Local boolean IsBlankStrip(striparray *mystripe, long deep);

void   striprint(long, long);
long showvrmlparms(long vrmltreecolor, long vrmlnamecolor, 
                long vrmlskycolornear, long vrmlskycolorfar, 
                long vrmlgroundcolornear);
void getvrmlparms(long *vrmltreecolor, long *vrmlnamecolor, 
                long *vrmlskycolornear,long *vrmlskycolorfar, 
                long *vrmlgroundcolornear,long *vrmlgroundcolorfar, 
                long numtochange);

#ifdef QUICKC
void   setupgraphics(void);
#endif


long   showrayparms(long, long, long, long, long, long);
void   getrayparms(long *, long *, long *, long *, long *,long *, long);

int    readafmfile(char *, short *);

void   metricforfont(char *, short *);
void   plotchar(long *, struct LOC_plottext *);
void   swap_charptr(char **, char **);
void   plotpb(void);
char  *findXfont(char *, double, double *, int *);
int    macfontid(char *);
int figfontid(char *fontname);
                        
void   makebox(char *, double *, double *, double *, long);
/* function prototypes */
#endif


int pointinrect(double x,double y,double x0,double y0,double x1,double y1)
{
  double tmp;
  if (x0 > x1)
    tmp = x0,
      x0  = x1,
      x1  = tmp;
  if (y0 > y1)
    tmp = y0,
      y0  = y1,
      y1  = tmp;
  
  return ((x >= x0 && x <= x1) && (y >= y0 && y <= y1));
}  /* pointinrect */

int rectintersects(double xmin1,double ymin1,double xmax1,double ymax1,
                double xmin2,double ymin2,double xmax2,double ymax2)
{
  double temp;
  
  /* check if any of the corners of either square are contained within the  *
   * other one. This catches MOST cases, the last one (two) is two thin     *
   * bands crossing each other (like a '+' )                                */
  
  if (xmin1 > xmax1){
    temp  = xmin1;  xmin1 = xmax1;  xmax1 = temp;}
  if (xmin2 > xmax2){
    temp  = xmin2;  xmin2 = xmax2;  xmax2 = temp;}
  if (ymin1 > ymax1){
    temp  = ymin1;  ymin1 = ymax1;  ymax1 = temp;}
  if (ymin2 > ymax2){
    temp  = ymin2;  ymin2 = ymax2;  ymax2 = temp;}
  
  return (pointinrect(xmin1,ymin1,xmin2,ymin2,xmax2,ymax2) ||
          pointinrect(xmax1,ymin1,xmin2,ymin2,xmax2,ymax2) ||
          pointinrect(xmin1,ymax1,xmin2,ymin2,xmax2,ymax2) ||
          pointinrect(xmax1,ymax1,xmin2,ymin2,xmax2,ymax2) ||
          pointinrect(xmin2,ymin2,xmin1,ymin1,xmax1,ymax1) ||
          pointinrect(xmax2,ymin2,xmin1,ymin1,xmax1,ymax1) ||
          pointinrect(xmin2,ymax2,xmin1,ymin1,xmax1,ymax1) ||
          pointinrect(xmax2,ymax2,xmin1,ymin1,xmax1,ymax1) || 
          (xmin1 >= xmin2 && xmax1 <= xmax2 &&
           ymin2 >= ymin1 && ymax2 <= ymax1)                ||
          (xmin2 >= xmin1 && xmax2 <= xmax1 &&
           ymin1 >= ymin2 && ymax1 <= ymax2)); 
}  /* rectintersects */

void clearit()
{
  long i;
  
  if (previewer == tek)
    printf("%c\f", escape);
  else if (ansi || ibmpc)
#ifdef WIN32
    phyClearScreen();
#else
    printf("\033[2J\033[H");
#endif
  else {
    for (i = 1; i <= 24; i++)
      putchar('\n');
  }
#ifdef WIN32
  phyClearScreen();
#endif
}  /* clearit */


boolean isfigfont(char *fontname)
{
  int i;
  if (strcmp(fontname,"Hershey") == 0)
    return 1;
  for (i=0;i<34;++i)
    if (strcmp(fontname,figfonts[i]) == 0)
      break;
  return (i < 34);
}  /* isfigfont */


int figfontid(char *fontname)
{
  int i;
  for (i=0;i<34;++i)
    if (strcmp(fontname,figfonts[i]) == 0)
      return i;
  return -1;
}  /* figfontid */


const char *figfontname(int id)
{
  return figfonts[id];
}  /* figfontname */


void getpreview()
{
  long loopcount;
  Char ch;

  clearit();
  printf("\nWhich type of screen will it be previewed on?\n\n");
  printf("   type:       to choose one compatible with:\n\n");
  printf("        N         will not be previewed\n");
#ifdef DOS
  printf("        I         MSDOS graphics screens\n");
#endif
#ifdef MAC
  printf("        M         Macintosh screens\n");
#endif
#ifdef X
  printf("        X         X Windows display\n");
#endif
#ifdef WIN32
  printf("        W         MS Windows display\n");
#endif
  printf("        K         TeKtronix 4010 graphics terminal\n");
  printf("        D         DEC ReGIS graphics (VT240 terminal)\n");
  printf("        U         other: one you have inserted code for\n");
  loopcount = 0;
  do {
    printf(" Choose one: \n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    fflush(stdout);
    scanf("%c%*[^\n]", &ch);
    getchar();
    uppercase(&ch);
    countup(&loopcount, 10);
  }
#undef FOO
#ifdef DOS
#define FOO
  while (strchr("NIKDU",ch) == NULL);
#endif
#ifdef MAC
#define FOO
  while (strchr("NMKDU",ch) == NULL);
#endif
#ifdef X
#define FOO
  while (strchr("NXKDU",ch) == NULL);
#endif
#ifdef WIN32
#define FOO
  while (strchr("NWKDU",ch) == NULL);
#endif
#ifndef FOO
  while (strchr("NKDU",ch) == NULL);
#endif
  preview = true;
  switch (ch) {

  case 'N':
    preview = false;
    previewer = other;   /* Added by Dan F. */
    break;

  case 'I':
    previewer = ibm;
    break;

  case 'M':
    previewer = mac;
    break;

  case 'X':
    previewer = xpreview;
    break;

  case 'W':
    previewer = winpreview;
    break;

  case 'K':
    previewer = tek;
    break;

  case 'D':
    previewer = decregis;
    break;

  case 'U':
    previewer = other;
    break;
  }
  printf("\n\n\n");
}  /* getpreview */


void pout(long n)
{
#ifdef MAC
  if (previewing)
    printf("%*ld", (int)((long)(0.434295 * log((double)n) + 0.0001)), n);
  else
    fprintf(plotfile, "%*ld",
            (int)((long)(0.434295 * log((double)n) + 0.0001)), n);
#endif
#ifndef MAC
  if (previewing)
    printf("%*ld", (int)((long)(0.434295 * log((double)n) + 0.0001)), n);
  else
    fprintf(plotfile, "%*ld",
            (int)((long)(0.434295 * log((double)n) + 0.0001)), n);
#endif
}  /* pout */


long upbyte(long num)
{
  /* get upper nibble of byte */
  long Result = 0, i, j, bytenum, nibcount;
  boolean done;

  bytenum = 0;
  done = false;
  nibcount = 0;
  i = num / 16;
  i /= 16;
  j = 1;
  while (!done) {
    bytenum += (i & 15) * j;
    nibcount++;
    if (nibcount == 2) {
      Result = bytenum;
      done = true;
    } else {
      j *= 16;
      i /= 16;
    }
  }
  return Result;
}  /* upbyte */


long lobyte(long num)
{
  /* get low order nibble of byte */
  long Result = 0, i, j, bytenum, nibcount;
  boolean done;

  bytenum = 0;
  done = false;
  nibcount = 0;
  i = num;
  j = 1;
  while (!done) {
    bytenum += (i & 15) * j;
    nibcount++;
    if (nibcount == 2) {
      Result = bytenum;
      done = true;
    } else {
      j *= 16;
      i /= 16;
    }
  }
  return Result;
}  /* lobyte */


void pictoutint(FILE *file, long pictint)
{
char picthi, pictlo;

picthi = (char)(pictint / 256);
pictlo = (char)(pictint % 256);
fprintf(file, "%c%c", picthi, pictlo);
}


void initplotter(long ntips, char *fontname)
{
  long i,j, hres, vres;
  Char picthi, pictlo;
  long pictint;
  int padded_width, byte_width;
  unsigned int dummy1, dummy2;
  double viewangle;
  
  treeline = 0.18 * labelheight * yscale * expand;
  labelline = 0.06 * labelheight * yscale * expand;
  linewidth = treeline;
  if (dotmatrix ) {
    for (i = 0; i <= 50; i++) {   /* for fast circle calculations */
     for (j = 0; j <= 50; j++){
       rootmatrix[i][j] =
           (long)floor(sqrt((double)(i * i + j * j)) + 0.5);}
   }
  }
 switch (plotter) {

  case xpreview:
#ifdef X
    XGetGeometry(display,mainwin,
    &DefaultRootWindow(display),&x,&y,&width,&height,&dummy1,&dummy2);
    XClearWindow(display,mainwin);
#endif
    break;

  case tek:
    oldxhigh = -1.0;
    oldxlow = -1.0;
    oldyhigh = -1.0;
    oldylow = -1.0;
    nmoves = 0;       /* DLS/JMH -- See function  PLOT                  */
    if (previewing)   /* DLS/JMH                                        */
      printf("%c\f", escape);   /* DLS/JMH */
    else
      fprintf(plotfile, "%c\f", escape);
    break;

  case hp:
    fprintf(plotfile, "IN;SP1;VS10.0;\n");
    break;

  case ray:
    fprintf(plotfile, "report verbose\n");

    fprintf(plotfile, "screen %f %f\n", xsize, ysize);
    if (ysize >= xsize) {
      viewangle = 2 * atan(ysize / (2 * 1.21 * xsize)) * 180 / pi;
      fprintf(plotfile, "fov 45 %3.1f\n", viewangle);
      fprintf(plotfile, "light 1 point 0 %6.2f %6.2f\n",
              -xsize * 1.8, xsize * 1.5);
      fprintf(plotfile, "eyep %6.2f %6.2f %6.2f\n",
              xsize * 0.5, -xsize * 1.21, ysize * 0.55);
    } else {
      viewangle = 2 * atan(xsize / (2 * 1.21 * ysize)) * 180 / pi;
      fprintf(plotfile, "fov %3.1f 45\n", viewangle);
      fprintf(plotfile, "light 1 point 0 %6.2f %6.2f\n",
              -ysize * 1.8, ysize * 1.5);
      fprintf(plotfile, "eyep %6.2f %6.2f %6.2f\n",
              xsize * 0.5, -ysize * 1.21, ysize * 0.55);
    }

    fprintf(plotfile, "lookp %6.2f 0 %6.2f\n", xsize * 0.5, ysize * 0.5);
    fprintf(plotfile, "/* %.10s */\n", colors[treecolor - 1].name);
    fprintf(plotfile,
            "surface treecolor diffuse %5.2f%5.2f%5.2f specular 1 1 1 specpow 30\n",
            colors[treecolor - 1].red, colors[treecolor - 1].green,
            colors[treecolor - 1].blue);
    fprintf(plotfile, "/* %.10s */\n", colors[namecolor - 1].name);
    fprintf(plotfile,
            "surface namecolor diffuse %5.2f%5.2f%5.2f specular 1 1 1 specpow 30\n",
            colors[namecolor - 1].red, colors[namecolor - 1].green,
            colors[namecolor - 1].blue);
    fprintf(plotfile, "/* %.10s */\n", colors[backcolor - 1].name);
    fprintf(plotfile, "surface backcolor diffuse %5.2f%5.2f%5.2f\n\n",
            colors[backcolor - 1].red, colors[backcolor - 1].green,
            colors[backcolor - 1].blue);
    treeline = 0.27 * labelheight * yscale * expand;
    linewidth = treeline;
    raylinewidth = treeline;
    if (grows == vertical)
      fprintf(plotfile, "plane backcolor 0 0 %2.4f 0 0 1\n", ymargin);
    else
      fprintf(plotfile, "plane backcolor 0 0 %2.4f 0 0 1\n",
              ymargin - ysize / (ntips - 1));

    fprintf(plotfile, "\nname tree\n");
    fprintf(plotfile, "grid 22 22 22\n");
    break;

  case pov:
    fprintf(plotfile, "// Declare the colors\n\n");

    fprintf(plotfile, "#declare C_Tree        = color rgb<%6.2f, %6.2f, %6.2f>;\n",
            colors[treecolor-1].red,
            colors[treecolor-1].green,
            colors[treecolor-1].blue);

    fprintf(plotfile, "#declare C_Name        = color rgb<%6.2f, %6.2f, %6.2f>;\n\n",
            colors[namecolor-1].red,
            colors[namecolor-1].green,
            colors[namecolor-1].blue);

    fprintf(plotfile, "// Declare the textures\n\n");

    fprintf(plotfile, "#declare %s = texture { pigment { C_Tree }\n", TREE_TEXTURE);
    fprintf(plotfile, "\t\tfinish { phong 1 phong_size 100 }};\n");

    fprintf(plotfile, "#declare %s = texture { pigment { C_Name }\n", NAME_TEXTURE);
    fprintf(plotfile, "\t\tfinish { phong 1 phong_size 100 }};\n");

    fprintf(plotfile, "\n#global_settings { assumed_gamma 2.2 }\n\n");
    
    fprintf(plotfile, "light_source { <0, %6.2f, %6.2f> color <1,1,1> }\n\n",
            xsize * 1.8, xsize * 1.5);

    /* The camera location */
    fprintf(plotfile, "camera {\n");

    if (ysize >= xsize) {
      fprintf(plotfile, "\tlocation <%6.2f, %6.2f, %6.2f>\n",
              -xsize * 0.5, -xsize * 1.21, ysize * 0.55);
    } else {
      fprintf(plotfile, "\tlocation <%6.2f, %6.2f, %6.2f>\n",
              -xsize * 0.5, -ysize * 1.21, ysize * 0.55);
    }
      
    fprintf(plotfile, "\tlook_at <%6.2f, 0, %6.2f>\n",
            -xsize * 0.5, ysize * 0.5);

    /* Handily, we can rotate since the rayshade paradigm ain't
       exactly congruent to the povray paradigm */
    fprintf(plotfile, "\trotate z*180\n");
    fprintf(plotfile, "}\n\n");

    fprintf(plotfile, "#background { color rgb <%6.2f, %6.2f, %6.2f> }\n\n",
            colors[backcolor-1].red,
            colors[backcolor-1].green,
            colors[backcolor-1].blue);

    if (bottomcolor != NO_PLANE) {
      /* The user wants a plane on the bottom... */
      if (grows == vertical) 
        fprintf(plotfile, "plane { z, %2.4f\n", 0.0 /*ymargin*/);
      else
        fprintf(plotfile, "plane { z, %2.4f\n",
                ymargin - ysize / (ntips - 1));
      
      fprintf(plotfile, "\tpigment {color rgb <%6.2f, %6.2f, %6.2f> }}\n\n",
              colors[bottomcolor-1].red,
              colors[bottomcolor-1].green,
              colors[bottomcolor-1].blue);
    }
    treeline = 0.27 * labelheight * yscale * expand;
    linewidth = treeline;
    raylinewidth = treeline;
    fprintf(plotfile, "\n// First, the tree\n\n");
    break;

  case vrml:
    vrmllinewidth = treeline;
    break;

  case pict:
    plotfile = freopen(pltfilename,"wb",plotfile);
    for (i=0;i<512;++i)
      putc('\000',plotfile);
    pictoutint(plotfile,1000); /* size...replaced later with seek */
    pictoutint(plotfile,1);    /* bbx0   */
    pictoutint(plotfile,1);    /* bby0   */
    pictoutint(plotfile,612);  /* bbx1   */
    pictoutint(plotfile,792);  /* bby1   */
    fprintf(plotfile,"%c%c",0x11,0x01); /* version "1" (B&W) pict */
    fprintf(plotfile,"%c%c%c",0xa0,0x00,0x82);
    fprintf(plotfile,"%c",1);    /* clip rect */
    pictoutint(plotfile,10);  /* region size, bytes. */
    pictoutint(plotfile,1);   /* clip x0             */
    pictoutint(plotfile,1);   /* clip y0             */
    pictoutint(plotfile,612); /* clip x1             */
    pictoutint(plotfile,792); /* clip y1             */
   
    bytewrite+=543;

    oldpictint = 0;
    pictint = (long)(linewidth + 0.5);
    if (pictint == 0)
      pictint = 1;
    picthi = (Char)(pictint / 256);
    pictlo = (Char)(pictint % 256);
    fprintf(plotfile, "\007%c%c%c%c", picthi, pictlo, picthi, pictlo);
    /* Set pen size for drawing tree. */
    break;

  case bmp:
    plotfile = freopen(pltfilename,"wb",plotfile);
    write_bmp_header(plotfile, (int)(xsize*xunitspercm), (int)(ysize*yunitspercm));
    byte_width   = (int) ceil (xsize / 8.0);
    padded_width = ((byte_width + 3) / 4) * 4 ;
    full_pic     = (byte *) Malloc ((padded_width *2) * (int) ysize) ;
    break ;

  case xbm:  /* what a completely verbose data representation format! */

    fprintf(plotfile, "#define drawgram_width %5ld\n",
            (long)(xunitspercm * xsize));
    fprintf(plotfile, "#define drawgram_height %5ld\n",
            (long)(yunitspercm * ysize));
    fprintf(plotfile, "static char drawgram_bits[] = {\n");
    /*filesize := 53;  */
    break;

  case lw:     /* write conforming postscript */
    fprintf(plotfile,"%%!PS-Adobe-2.0\n");
    fprintf(plotfile,"%%%%Title: Phylip Tree Output\n");
    fprintf(plotfile,"%%%%DocumentFonts: (atend)\n");
    fprintf(plotfile,"%%%%Pages: %d 1\n",
           ((int)((pagex-hpmargin-0.01)/(paperx-hpmargin))+1)*
               ((int)((pagey-vpmargin-0.01)/papery-vpmargin)+1));
    fprintf(plotfile,"%%%%BoundingBox: 0 0 612 792\n");
    fprintf(plotfile,"%%%%DocumentPaperSizes: Letter\n"); /* this may not be right */
    fprintf(plotfile,"%%%%Orientation: Portrait\n"); 
    fprintf(plotfile,"%%%%EndComments\n");
    fprintf(plotfile,"/l {newpath moveto lineto stroke} def\n");
    fprintf(plotfile,"%%%%EndProlog\n%%%%\n");
    fprintf(plotfile,"%%%%Page: 1 1\n");
    fprintf(plotfile,"%%%%PageBoundingBox: 0 0 %d %d\n",
            (int)(xunitspercm*paperx),(int)(yunitspercm*papery));
    fprintf(plotfile,"%%%%PageFonts: (atend)\n%%%%BeginPageSetup\n");
    fprintf(plotfile,"%%%%PaperSize: Letter\n");
    fprintf(plotfile," 1 setlinecap \n 1 setlinejoin  \n");
    fprintf(plotfile, "%8.2f setlinewidth newpath \n", treeline);
    break;

  case idraw:
    fprintf(plotfile, "%%I Idraw 9 Grid 8 \n\n");
    fprintf(plotfile,"%%%%Page: 1 1\n\n");
    fprintf(plotfile,"Begin\n");
    fprintf(plotfile,"%%I b u\n");
    fprintf(plotfile,"%%I cfg u\n");
    fprintf(plotfile,"%%I cbg u\n");
    fprintf(plotfile,"%%I f u\n");
    fprintf(plotfile,"%%I p u\n");
    fprintf(plotfile,"%%I t\n");
    fprintf(plotfile,"[ 0.679245 0 0 0.679245 0 0 ] concat\n"); 

    fprintf(plotfile,"/originalCTM matrix currentmatrix def\n\n");
    break;

  case ibm:
#ifdef TURBOC
    initgraph(&GraphDriver,&HiMode,"");
#endif
#ifdef QUICKC
    setupgraphics();
#endif
    break;

  case mac:
#ifdef MAC
    gfxmode();
    pictint=(long)(linewidth + 0.5);
#endif
    break;

  case houston:
    break;

  case decregis:
    oldx = (double) 300;
    oldy = (double) 1;
    nmoves = 0;
    if (previewing)
          printf("%c[2J%cPpW(I3);S(A[0,0][799,479]);S(I(W))S(E);S(C0);W(I(D))\n",
                 escape,escape);
     else
       fprintf(plotfile,
               "%c[2J%cPpW(I3);S(A[0,0][799,479]);S(I(W))S(E);S(C0);W(I(D))\n",
               escape,escape);
    break;

  case epson:
    plotfile = freopen(pltfilename,"wb",plotfile);
    fprintf(plotfile, "\0333\030");
    break;

  case oki:
    plotfile = freopen(pltfilename,"wb",plotfile);
    fprintf(plotfile, "\033%%9\020");
    break;

  case citoh:
    plotfile = freopen(pltfilename,"wb",plotfile);
    fprintf(plotfile, "\033T16");
    break;

  case toshiba: /* reopen in binary since we always need \n\r on the file */
                /* and dos in text mode puts it, but unix does not        */
    plotfile = freopen(pltfilename,"wb",plotfile);
    fprintf(plotfile, "\033\032I\n\r\n\r");
    fprintf(plotfile, "\033L06\n\r");
    break;

  case pcl:
    plotfile = freopen(pltfilename,"wb",plotfile);
    if (hpresolution == 150 || hpresolution == 300)
      fprintf(plotfile, "\033*t%3ldR", hpresolution);
    else if (hpresolution == 75)
      fprintf(plotfile, "\033*t75R");
    break;

  case pcx:
    plotfile = freopen(pltfilename,"wb",plotfile);
    fprintf(plotfile,"\012\003\001\001%c%c%c%c",0,0,0,0);
  /* Manufacturer version (1 byte) version (1 byte), encoding (1 byte),
     bits per pixel (1 byte), xmin (2 bytes) ymin (2 bytes),
     Version */
    hres = strpwide;
    vres = (long)floor(yunitspercm * ysize + 0.5);
    fprintf(plotfile, "%c%c", (unsigned char)lobyte(hres - 1),
            (unsigned char)upbyte(hres - 1)); /* Xmax */
    fprintf(plotfile, "%c%c", (unsigned char)lobyte(vres - 1),
            (unsigned char)upbyte(vres - 1)); /* Ymax */
    fprintf(plotfile, "%c%c", (unsigned char)lobyte(hres),
            (unsigned char)upbyte(hres));
    /* Horizontal resolution */
    fprintf(plotfile, "%c%c", (unsigned char)lobyte(vres),
            (unsigned char)upbyte(vres));
    /* Vertical resolution */
    for (i = 1; i <= 48; i++)   /* Color Map */
      putc('\000', plotfile);
    putc('\000', plotfile);
    putc('\001', plotfile);   /* Num Planes */
    putc(hres / 8, plotfile);   /* Bytes per line */
    putc('\000',plotfile);
    for (i = 1; i <= 60; i++)   /* Filler */
      putc('\000',plotfile);
    break;

  case fig:
    fprintf(plotfile, "#FIG 2.0\n");
    fprintf(plotfile, "80 2\n");
    break;

 case gif:
 case other:
    break;
 default:        /* case vrml not handled        */
    break;
    /* initialization code for a new plotter goes here */
  }
}  /* initplotter */


void finishplotter()
{
  int padded_width, byte_width; /* For bmp code */

  switch (plotter) {

  case xpreview:
#ifdef X
         plotter=oldplotter;
         redraw(NULL,NULL,NULL);
         XtAppMainLoop(appcontext);

#endif
    break;

  case tek:
    if (previewing) {
      fflush(stdout);
      scanf("%*c%*[^\n]");
      getchar();
      printf("%c\f", escape);
    } else {
      putc('\n', plotfile);
      plot(penup, 1.0, 1.0);
    }
    break;

  case hp:
    plot(penup, 1.0, 1.0);
    fprintf(plotfile, "SP;\n");
    break;

  case ray:
    fprintf(plotfile,"end\n\nobject treecolor tree\n");
    fprintf(plotfile,"object namecolor species_names\n");
    break;

  case pov:
    break;

  case pict:
    fprintf(plotfile,"%c%c%c%c%c",0xa0,0x00,0x82,0xff,0x00);
    bytewrite+=5;
    fseek(plotfile,512L,SEEK_SET);
    pictoutint(plotfile,bytewrite);
    break;

  case lw:
    fprintf(plotfile, "stroke showpage \n\n");
    fprintf(plotfile,"%%%%PageTrailer\n");
    fprintf(plotfile,"%%%%PageFonts: %s\n",
            (strcmp(fontname,"Hershey") == 0) ? "" : fontname);
    fprintf(plotfile,"%%%%Trailer\n");
    fprintf(plotfile,"%%%%DocumentFonts: %s\n",
            (strcmp(fontname,"Hershey") == 0) ? "" : fontname);
    break;

  case idraw:
    fprintf(plotfile, "\nEnd %%I eop\n\n");
    fprintf(plotfile, "showpage\n\n");
    fprintf(plotfile, "%%%%Trailer\n\n");
    fprintf(plotfile, "end\n");
    break;

  case ibm:
#ifdef TURBOC
    getchar();
    restorecrtmode();
#endif
#ifdef QUICKC
    getchar();
    _clearscreen(_GCLEARSCREEN);
    _setvideomode(_DEFAULTMODE);
#endif
    break;

  case mac:
#ifdef MAC
    plotter=oldplotter;
    eventloop();
#endif
    break;

  case houston:
    break;

  case decregis:
    plot(penup, 1.0, 1.0);
    if (previewing)
      printf("%c\\", escape);
    else
      fprintf(plotfile, "%c\\", escape);
    if (previewing) {
      getchar();
      printf("%c[2J",escape);
    }
    break;

  case epson:
    fprintf(plotfile, "\0333$");
    break;

  case oki:
    /* blank case */
    break;

  case citoh:
    fprintf(plotfile, "\033A");
    break;

  case toshiba:
    fprintf(plotfile, "\033\032I\n\r");
    break;

  case pcl:
    fprintf(plotfile, "\033*rB");    /* Exit graphics mode */
    putc('\f', plotfile);            /* just to make sure? */
    break;

  case pcx:
    /* blank case */
    break;

  case bmp:
    byte_width = (int) ceil (xsize / 8.0);
    padded_width = ((byte_width + 3) / 4) * 4 ;
    turn_rows (full_pic, padded_width, (int) ysize);
    write_full_pic(full_pic, total_bytes);
    free (full_pic) ;
    break;

  case xbm:
    fprintf(plotfile, "}\n");
    break;

  case fig:
    /* blank case */
    break;

  case gif:
  case other:
    break;
  default:        /* case vrml not handled        */
    break;
    /* termination code for a new plotter goes here */
  }
}  /* finishplotter */


Local long SFactor()
{
  /* the dot-skip is resolution-independent. */
  /* this makes all the point-skip instructions skip the same # of dots. */
  long Result = 0;

  if (hpresolution == 150)
    Result = 2;
  if (hpresolution == 300)
    Result = 1;
  if (hpresolution == 75)
    return 4;
  return Result;
}  /* SFactor */


long DigitsInt(long x)
{
  if (x < 10)
    return 1;
  else if (x >= 10 && x < 100)
    return 2;
  else
    return 3;
}  /* DigistInt */


Local boolean IsColumnEmpty(striparray *mystripe, long pos, long deep)
{
  long j;
  boolean ok;

  ok = true;
  j = 1;
  while (ok && j <= deep) {
    ok = (ok && mystripe[j - 1][pos - 1] == null);
    j++;
  }
  return ok;
}  /* IsColumnEmpty */


void Skip(long Amount)
{
  /* assume we're not in gfx mode. */
  fprintf(plotfile, "\033&f1S");   /* Pop the graphics cursor    */
#ifdef MAC
  fprintf(plotfile, "\033*p+%*ldX",
          (int)DigitsInt(Amount * SFactor()), Amount * SFactor());
#endif
#ifndef MAC
  fprintf(plotfile, "\033*p+%*ldX",
          (int)DigitsInt(Amount * SFactor()), Amount * SFactor());
#endif
  fprintf(plotfile, "\033&f0S");   /* Push the cursor to new location */
  filesize += 15 + DigitsInt(Amount * SFactor());
}  /* Skip */


Local long FirstBlack(striparray *mystripe, long startpos, long deep)
{
  /* returns, given a strip and a position, next x with some y's nonzero */
  long i;
  boolean columnempty;

  i = startpos;
  columnempty = true;
  while (columnempty && i < strpwide / 8) {
    columnempty = (columnempty && IsColumnEmpty(mystripe, i,deep));
    if (columnempty)
      i++;
  }
  return i;
}  /* FirstBlack */


Local long FirstWhite(striparray *mystripe, long startpos, long deep)
{
  /* returns, given a strip and a position, the next x with all y's zero */
  long i;
  boolean columnempty;

  i = startpos;
  columnempty = false;
  while (!columnempty && i < strpwide / 8) {
    columnempty = IsColumnEmpty(mystripe, i,deep);
    if (!columnempty)
      i++;
  }
  return i;
}  /* FirstWhite */


Local boolean IsBlankStrip(striparray *mystripe, long deep)
{
  long i, j;
  boolean ok;

  ok = true;
  i = 1;
  while (ok && i <= strpwide / 8) {
    for (j = 0; j < (deep); j++)
      ok = (ok && mystripe[j][i - 1] == '\0');
    i++;
  }
  return ok;
}  /* IsBlankStrip */


void striprint(long div, long deep)
{
  long i, j, t, x, theend, width;
  unsigned char counter;
  boolean done;
  done = false;
  width = strpwide;
  if (plotter != pcx && plotter != pcl &&
      plotter != bmp && plotter != xbm) {
    while (!done) {
      for (i = 0; i < div; i++)
        done = done || (stripe[i] && (stripe[i][width - 1] != null));
      if (!done)
        width--;
      done = (done || width == 0);
    }
  }
  switch (plotter) {

  case epson:
    if (!empty) {
      fprintf(plotfile, "\033L%c%c", (char) width & 255, (char) width / 256);
      for (i = 0; i < width; i++)
        putc(stripe[0][i], plotfile);
      filesize += width + 4;
    }
    putc('\n', plotfile);
    putc('\r', plotfile);
    break;

  case oki:
    if (!empty) {
      fprintf(plotfile, "\033%%1%c%c", (char) width / 128, (char) width & 127);
      for (i = 0; i < width; i++)
        putc(stripe[0][i], plotfile);
      filesize += width + 5;
    }
    putc('\n', plotfile);
    putc('\r', plotfile);
    break;

  case citoh:
    if (!empty) {
      fprintf(plotfile, "\033S%04ld",width);
      for (i = 0; i < width; i++)
        putc(stripe[0][i], plotfile);
      filesize += width + 6;
    }

    putc('\n', plotfile);
    putc('\r', plotfile);
    break;

  case toshiba:
    if (!empty) {
      for (i = 0; i < width; i++) {
        for (j = 0; j <= 3; j++)
          stripe[j][i] += 64;
      }
      fprintf(plotfile, "\033;%04ld",width);

      for (i = 0; i < width; i++)
        fprintf(plotfile, "%c%c%c%c",
                stripe[0][i], stripe[1][i], stripe[2][i], stripe[3][i]);
      filesize += width * 4 + 6;
    }
    putc('\n', plotfile);
    putc('\r', plotfile);
    break;

  case pcx:
    width = strpwide / 8;
    for (j = 0; j < div; j++) {
      t = 1;
      while (1) {
        i = 0; /* i == RLE count ???? */
        while ((stripe[j][t + i - 1]) == (stripe[j][t + i])
               && t + i < width && i < 63)
          i++;
        if (i > 0) {
          counter = 192;
          counter += i;
          putc(counter, plotfile);
          putc(255 - stripe[j][t - 1], plotfile);
          t += i;
          filesize += 2;
        } else {
          if (255 - (stripe[j][t - 1] & 255) >= 192) {
            putc(193, plotfile);
            filesize++;
          }
          putc(255 - stripe[j][t - 1], plotfile);
          t++;
          filesize++;

        }
        if (t >width) break;
      }
    }
    break;

  case pcl:
    width = strpwide / 8;
    if (IsBlankStrip(stripe,deep)) {
#ifdef MAC
      fprintf(plotfile, "\033&f1S\033*p0X\033*p+%*ldY\033&f0S",
              (int)DigitsInt(deep * SFactor()), deep * SFactor());
#endif
#ifndef MAC
      fprintf(plotfile, "\033&f1S\033*p0X\033*p+%*dY\033&f0S",
              (int)DigitsInt(deep * SFactor()), (int) (deep * SFactor()));
#endif
      filesize += DEFAULT_STRIPE_HEIGHT + DigitsInt(deep * SFactor());
    } else {  /* plotting the actual strip as bitmap data */
      x = 1;
      theend = 1;
      while (x < width) {
        x = FirstBlack(stripe, x,deep);    /* all-black strip is now    */
        Skip((x - theend - 1) * 8);        /* x..theend                 */
        theend = FirstWhite(stripe, x,deep) - 1;/* like lastblack            */
        fprintf(plotfile, "\033*r1A");     /* enter gfx mode            */
        for (j = 0; j < div; j++) {
#ifdef MAC
          fprintf(plotfile, "\033*b%*ldW",
                  (int)DigitsInt(theend - x + 1), theend - x + 1);
#endif
#ifndef MAC
          fprintf(plotfile, "\033*b%*dW",
                  (int)DigitsInt(theend - x + 1), (int) (theend - x + 1));
#endif
              /* dump theend-x+1 bytes */
          for (t = x - 1; t < theend; t++)
            putc(stripe[j][t], plotfile);
          filesize += theend - x + DigitsInt(theend - x + 1) + 5;
        }
        fprintf(plotfile, "\033*rB");   /* end gfx mode */
        Skip((theend - x + 1) * 8);
        filesize += 9;
        x = theend + 1;
      }
      fprintf(plotfile, "\033&f1S");   /* Pop cursor  */
#ifdef MAC
      fprintf(plotfile, "\033*p0X\033*p+%*ldY",
              (int)DigitsInt(deep * SFactor()), deep * SFactor());
#endif
#ifndef MAC
      fprintf(plotfile, "\033*p0X\033*p+%*dY",
              (int)DigitsInt(deep * SFactor()), (int) (deep * SFactor()));
#endif
      filesize += DEFAULT_STRIPE_HEIGHT + DigitsInt(deep * SFactor());
      fprintf(plotfile, "\033&f0S");   /* Push cursor  */
    }
    break;
    /* case for hpcl code */


  case bmp:
    width = ((strpwide -1) / 8) +1;
    translate_stripe_to_bmp (&stripe, full_pic, increment++,
                             width, div, &total_bytes) ;
    break;
    /* case for bmp code */

  case xbm:
    x = 0;   /* count up # of bytes so we can put returns. */
    width = ((strpwide -1) / 8) +1;
    for (j = 0; j <  div; j++) {
      for (i = 0; i < width; i++) {
        fprintf(plotfile, "0x%02x,",(unsigned char)stripe[j][i]);
        filesize += 5;
        x++;
        if ((x % 15) == 0) {
          putc('\n', plotfile);
          filesize++;
        }
      }
    }
   putc('\n',plotfile);
   break;

  case lw:
  case hp:
  case xpreview:
  case winpreview:
  case tek:
  case ibm:
  case mac:
  case houston:
  case decregis:
  case fig:
  case pict:
  case ray:
  case pov:
  case gif:
  case idraw:
  case other:
    break;
  default:      /* case vrml not handled        */
    break;
    /* graphics print code for a new printer goes here */
  }
}  /* striprint */


#ifdef QUICKC
void setupgraphics()
{
_getvideoconfig(&myscreen);
#ifndef WATCOM
switch(myscreen.adapter){
  case _CGA:
  case _OCGA:
   _setvideomode(_HRESBW);
    break;
  case _EGA:
  case _OEGA:
    _setvideomode(_ERESNOCOLOR);
  case _VGA:
  case _OVGA:
  case _MCGA:
    _setvideomode(_VRES2COLOR);
     break;
  case _HGC:
    _setvideomode(_HERCMONO);
     break;
  default:
     printf("Your display hardware is unsupported by this program.\n");
      break;
}
#endif
#ifdef WATCOM
switch(myscreen.adapter){
  case _VGA:
  case _SVGA:
      _setvideomode(_VRES16COLOR);
      break;
  case _MCGA:
      _setvideomode(_MRES256COLOR);
      break;
  case _EGA:
     _setvideomode(_ERESNOCOLOR);
     break;
  case _CGA:
     _setvideomode(_MRES4COLOR);
     break;
  case _HERCULES:
     _setvideomode(_HERCMONO);
     break;
  default:
     printf("Your display hardware is unsupported by this program.\n");
     exxit(-1);
     break;
   }
#endif
_getvideoconfig(&myscreen);
_setlinestyle(0xffff);
xunitspercm=myscreen.numxpixels / 25;
yunitspercm=myscreen.numypixels / 17.5;
xsize = 25.0;
ysize = 17.5;
}  /* setupgraphics */
#endif


void loadfont(short *font, char *application)
{

  FILE *fontfile;
  long i, charstart = 0, dummy;
  Char ch = 'A';

  i=0;
  openfile(&fontfile,FONTFILE,"font file","r",application,NULL);

  while (!(eoff(fontfile) || ch == ' ')) {
    charstart = i + 1;
    if (fscanf(fontfile, "%c%c%ld%hd%hd", &ch, &ch, &dummy, &font[charstart + 1],
           &font[charstart + 2]) != 5) {
      printf("Error while reading fontfile\n\n");
      exxit(-1);
    }
    font[charstart] = ch;
    i = charstart + 3;
    do {
      if ((i - charstart - 3) % 10 == 0) 
        scan_eoln(fontfile);
      i++;
      if (fscanf(fontfile, "%hd", &font[i - 1]) != 1) {
        printf("Error while reading fontfile\n\n");
        exxit(-1);
      }

    } while (abs(font[i - 1]) < 10000);
    scan_eoln(fontfile);
    font[charstart - 1] = i + 1;
  }
  font[charstart - 1] = 0;
 FClose(fontfile);
}  /* loadfont */


long showrayparms(long treecolor, long namecolor, long backcolor,
                        long bottomcolor, long rx, long ry)
{
  long i, loopcount;
  Char ch,input[32];
  long numtochange;

  if (previewer == tek)
    printf("%c\f", escape);
  else {
    for (i = 1; i <= 24; i++)
      putchar('\n');
  }
  if (plotter == ray) {
    printf("Settings for Rayshade file: \n\n");
    printf(" (1)               Tree color:  %.10s\n",colors[treecolor-1].name);
    printf(" (2)      Species names color:  %.10s\n",colors[namecolor-1].name);
    printf(" (3)         Background color:  %.10s\n",colors[backcolor-1].name);
    printf(" (4)               Resolution:  %2ld X %2ld\n\n",rx,ry);
  } else if (plotter == pov) {
    printf("Settings for POVray file: \n\n");
    printf(" (1)               Tree color:  %.10s\n",colors[treecolor-1].name);
    printf(" (2)      Species names color:  %.10s\n",colors[namecolor-1].name);
    printf(" (3)         Background color:  %.10s\n",colors[backcolor-1].name);
    printf(" (4)             Bottom plane:  %.10s\n",
           bottomcolor == NO_PLANE ?
           "(none)\0" : colors[bottomcolor-1].name);
  }
    
  printf(" Do you want to accept these? (Yes or No)\n");
  loopcount = 0;
  for (;;) {
    printf(" Type Y or N or the number (1-4) of the one to change: \n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    getstryng(input);
    numtochange=atoi(input);
    uppercase(&input[0]);
    ch=input[0];
    if (ch == 'Y' || ch == 'N' || (numtochange >= 1 && numtochange <= 4))
      break;
    countup(&loopcount, 10);
  }
 return (ch == 'Y') ? -1 : numtochange;
}  /* showrayparms */


void getrayparms(long *treecolor, long *namecolor, long *backcolor,
                        long *bottomcolor, long *rx,long *ry, long numtochange)
{
  Char ch;
  long i, loopcount;

  if (numtochange == 0) {
    loopcount = 0;
    do {
      printf(" Type the number of one that you want to change (1-4):\n");
#ifdef WIN32
      phyFillScreenColor();
#endif
      fflush(stdout);
      scanf("%ld%*[^\n]", &numtochange);
      getchar();
      countup(&loopcount, 10);
    } while (numtochange < 1 || numtochange > 10);
  }
  switch (numtochange) {

  case 1:
    printf("\nWhich of these colors will the tree be?:\n");
    printf("   White, Red, Orange, Yellow, Green, Blue, or Violet\n");
    printf(" (W, R, O, Y, G, B, or V)\n");
    loopcount = 0;
    do {
      printf(" Choose one: \n");
#ifdef WIN32
      phyFillScreenColor();
#endif
      fflush(stdout);
      scanf("%c%*[^\n]", &ch);
      getchar();
      if (ch == '\n')
        ch = ' ';
      uppercase(&ch);
      (*treecolor) = 0;
      for (i = 1; i <= 7; i++) {
        if (ch == colors[i - 1].name[0]) {
          (*treecolor) = i;
          return;
        }
      }
      countup(&loopcount, 10);
    } while ((*treecolor) == 0);
    break;
    
  case 2:
    printf("\nWhich of these colors will the species names be?:\n");
    printf("   White, Red, Orange, Yellow, Green, Blue, or Violet\n");
    printf(" (W, R, O, Y, G, B, or V)\n");
    loopcount = 0;
    do {
      printf(" Choose one: \n");
#ifdef WIN32
      phyFillScreenColor();
#endif
      fflush(stdout);
      scanf("%c%*[^\n]", &ch);
      getchar();
      if (ch == '\n')
        ch = ' ';
      uppercase(&ch);
      (*namecolor) = 0;
      for (i = 1; i <= 7; i++) {
        if (ch == colors[i - 1].name[0]) {
          (*namecolor) = i;
          return;
        }
      }
      countup(&loopcount, 10);
    } while ((*namecolor) == 0);
    break;

  case 3:
    printf("\nWhich of these colors will the background be?:\n");
    printf("   White, Red, Orange, Yellow, Green, Blue, or Violet\n");
    printf(" (W, R, O, Y, G, B, or V)\n");
    loopcount = 0;
    do {
      printf(" Choose one: \n");
#ifdef WIN32
      phyFillScreenColor();
#endif
      fflush(stdout);
      scanf("%c%*[^\n]", &ch);
      getchar();
      if (ch == '\n')
        ch = ' ';
      uppercase(&ch);
      (*backcolor) = 0;
      for (i = 1; i <= 7; i++) {
        if (ch == colors[i - 1].name[0]) {
          (*backcolor) = i;
          return;
        }
      }
      countup(&loopcount, 10);
    } while ((*backcolor) == 0);
    break;

  case 4:
    /* Dimensions for rayshade, bottom plane for povray */
    if (plotter == pov) {
      printf("\nWhich of these colors will the bottom plane be?:\n");
      printf("   White, Red, Orange, Yellow, Green, Blue, Violet, or None (no plane)\n");
      printf(" (W, R, O, Y, G, B, V, or N)\n");
      loopcount = 0;
      do {
        printf(" Choose one: \n");
#ifdef WIN32
        phyFillScreenColor();
#endif
        fflush(stdout);
        scanf("%c%*[^\n]", &ch);
        getchar();
        if (ch == '\n')
          ch = ' ';
        uppercase(&ch);
        /* If the user doesn't want a bottom plane. . . */
        if (ch == 'N') {
          (*bottomcolor) = NO_PLANE;
          return;
        } else {
          (*bottomcolor) = 0;
          for (i = 1; i <= 7; i++) {
            if (ch == colors[i - 1].name[0]) {
              (*bottomcolor) = i;
              return;
            }
          }
        }
      countup(&loopcount, 10);
      } while ((*bottomcolor) == 0);

    } else if (plotter == ray) {
      printf("\nEnter the X resolution:\n");
#ifdef WIN32
      phyFillScreenColor();
#endif
      fflush(stdout);
      scanf("%ld%*[^\n]", rx);
      getchar();
      printf("Enter the Y resolution:\n");
#ifdef WIN32
      phyFillScreenColor();
#endif
      fflush(stdout);
      scanf("%ld%*[^\n]",ry);
      getchar();
    }
    break;
  }
}  /* getrayparms */


long showvrmlparms(long vrmltreecolor, long vrmlnamecolor, long vrmlskycolornear,
                        long vrmlskycolorfar, long vrmlgroundcolornear)
{
  long i, loopcount;
  Char ch,input[32];
  long numtochange;

  if (previewer == tek)
    printf("%c\f", escape);
  else {
    for (i = 1; i <= 24; i++)
      putchar('\n');
  }
  printf("Settings for VRML file: \n\n");
  printf(" (1)               Tree color:  %.10s\n",colors[vrmltreecolor-1].name);
  printf(" (2)      Species names color:  %.10s\n",colors[vrmlnamecolor-1].name);
  printf(" (3)            Horizon color:  %.10s\n",colors[vrmlskycolorfar-1].name);
  printf(" (4)             Zenith color:  %.10s\n",colors[vrmlskycolornear-1].name);
  printf(" (5)             Ground color:  %.10s\n",colors[vrmlgroundcolornear-1].name);
    
  printf(" Do you want to accept these? (Yes or No)\n");
  loopcount = 0;
  for (;;) {
    printf(" Type Y or N or the number (1-5) of the one to change: \n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    getstryng(input);
    numtochange=atoi(input);
    uppercase(&input[0]);
    ch=input[0];
    if (ch == 'Y' || ch == 'N' || (numtochange >= 1 && numtochange <= 5))
      break;
    countup(&loopcount, 10);
  }
 return (ch == 'Y') ? -1 : numtochange;
}  /* showvrmlparms */


void getvrmlparms(long *vrmltreecolor, long *vrmlnamecolor, long *vrmlskycolornear,
                        long *vrmlskycolorfar, long *vrmlgroundcolornear,
                        long *vrmlgroundcolorfar, long numtochange)
{
  Char ch;
  long i, loopcount;

  if (numtochange == 0) {
    loopcount = 0;
    do {
      printf(" Type the number of one that you want to change (1-4):\n");
#ifdef WIN32
      phyFillScreenColor();
#endif
      fflush(stdout);
      scanf("%ld%*[^\n]", &numtochange);
      getchar();
      countup(&loopcount, 10);
    } while (numtochange < 1 || numtochange > 10);
  }
  switch (numtochange) {

  case 1:
    printf("\nWhich of these colors will the tree be?:\n");
    printf("   White, Red, Orange, Yellow, Green, Blue, or Violet\n");
    printf(" (W, R, O, Y, G, B, or V)\n");
    loopcount = 0;
    do {
      printf(" Choose one: \n");
#ifdef WIN32
      phyFillScreenColor();
#endif
      fflush(stdout);
      scanf("%c%*[^\n]", &ch);
      getchar();
      if (ch == '\n')
        ch = ' ';
      uppercase(&ch);
      (*vrmltreecolor) = 0;
      for (i = 1; i <= 7; i++) {
        if (ch == colors[i - 1].name[0]) {
          (*vrmltreecolor) = i;
          return;
        }
      }
      countup(&loopcount, 10);
    } while ((*vrmltreecolor) == 0);
    break;
    
  case 2:
    printf("\nWhich of these colors will the species names be?:\n");
    printf("   White, Red, Orange, Yellow, Green, Blue, or Violet\n");
    printf(" (W, R, O, Y, G, B, or V)\n");
    loopcount = 0;
    do {
      printf(" Choose one: \n");
#ifdef WIN32
      phyFillScreenColor();
#endif
      fflush(stdout);
      scanf("%c%*[^\n]", &ch);
      getchar();
      if (ch == '\n')
        ch = ' ';
      uppercase(&ch);
      (*vrmlnamecolor) = 0;
      for (i = 1; i <= 7; i++) {
        if (ch == colors[i - 1].name[0]) {
          (*vrmlnamecolor) = i;
          return;
        }
      }
      countup(&loopcount, 10);
    } while ((*vrmlnamecolor) == 0);
    break;

  case 3:
    printf("\nWhich of these colors will the horizon be?:\n");
    printf("   White, Red, Orange, Yellow, Green, Blue, or Violet\n");
    printf(" (W, R, O, Y, G, B, or V)\n");
    loopcount = 0;
    do {
      printf(" Choose one: \n");
#ifdef WIN32
      phyFillScreenColor();
#endif
      fflush(stdout);
      scanf("%c%*[^\n]", &ch);
      getchar();
      if (ch == '\n')
        ch = ' ';
      uppercase(&ch);
      (*vrmlskycolorfar) = 0;
      for (i = 1; i <= 7; i++) {
        if (ch == colors[i - 1].name[0]) {
          (*vrmlskycolorfar) = i;
          return;
        }
      }
      countup(&loopcount, 10);
    } while ((*vrmlskycolorfar) == 0);
    break;

  case 4:
    printf("\nWhich of these colors will the zenith be?:\n");
    printf("   White, Red, Orange, Yellow, Green, Blue, or Violet\n");
    printf(" (W, R, O, Y, G, B, or V)\n");
    loopcount = 0;
    do {
      printf(" Choose one: \n");
#ifdef WIN32
      phyFillScreenColor();
#endif
      fflush(stdout);
      scanf("%c%*[^\n]", &ch);
      getchar();
      if (ch == '\n')
        ch = ' ';
      uppercase(&ch);
      (*vrmlskycolornear) = 0;
      for (i = 1; i <= 7; i++) {
        if (ch == colors[i - 1].name[0]) {
          (*vrmlskycolornear) = i;
          return;
        }
      }
      countup(&loopcount, 10);
    } while ((*vrmlskycolornear) == 0);
    break;

  case 5:
    printf("\nWhich of these colors will the ground be?:\n");
    printf("   White, Red, Orange, Yellow, Green, Blue, or Violet\n");
    printf(" (W, R, O, Y, G, B, or V)\n");
    loopcount = 0;
    do {
      printf(" Choose one: \n");
#ifdef WIN32
      phyFillScreenColor();
#endif
      fflush(stdout);
      scanf("%c%*[^\n]", &ch);
      getchar();
      if (ch == '\n')
        ch = ' ';
      uppercase(&ch);
      (*vrmlgroundcolornear) = 0;
      for (i = 1; i <= 7; i++) {
        if (ch == colors[i - 1].name[0]) {
          (*vrmlgroundcolornear) = i;
          (*vrmlgroundcolorfar) = i;
          return;
        }
      }
      countup(&loopcount, 10);
    } while ((*vrmlgroundcolornear) == 0);
    break;

  }
}  /* gevrmlparms */


void plotrparms(long ntips)
{
  /* set up initial characteristics of plotter or printer */
  long rayresx, rayresy;
  long n, loopcount;
  double xsizehold, ysizehold;

  xsizehold = xsize;
  ysizehold = ysize;
  penchange = no;
  xcorner = 0.0;
  ycorner = 0.0;
  if (dotmatrix && (!previewing))
    strpdiv = 1;
  switch (plotter) {

  case ray:
    penchange = yes;
    xunitspercm = 1.0;
    yunitspercm = 1.0;
    xsize = 10.0;
    ysize = 10.0;
    rayresx = 512;
    rayresy = 512;
    treecolor = 6;
    namecolor = 4;
    backcolor = 1;
    /* MSVC gave warning that bottomcolor was uninitialized.  Unsure what
       this should be */
    bottomcolor = 1;
    loopcount = 0;
    do {
      n=showrayparms(treecolor,namecolor,backcolor,bottomcolor,rayresx,rayresy);
      if (n != -1)
        getrayparms(&treecolor,&namecolor,&backcolor,&bottomcolor,&rayresx,&rayresy,n);
      countup(&loopcount, 10);
    } while (n != -1);
    xsize = rayresx;
    ysize = rayresy;
    break;

  case pov:
    penchange = yes;
    xunitspercm = 1.0;
    yunitspercm = 1.0;
    xsize = 10.0;
    ysize = 10.0;
    rayresx = 512;
    rayresy = 512;
    treecolor = 6;
    namecolor = 4;
    backcolor = 1;
    bottomcolor = 1;
    loopcount = 0;
    do {
      n=showrayparms(treecolor,namecolor,backcolor,bottomcolor,rayresx,rayresy);
      if (n != -1)
        getrayparms(&treecolor,&namecolor,&backcolor,&bottomcolor,&rayresx,&rayresy,n);
      countup(&loopcount, 10);
    } while (n != -1);
    xsize = rayresx;
    ysize = rayresy;
    break;

  case vrml:
#ifndef MAC
    penchange = yes;
    xunitspercm = 1.0;
    yunitspercm = 1.0;
    xsize = 10.0;
    ysize = 10.0;
    vrmlplotcolor = treecolor;
    loopcount = 0;
    do {
      n=showvrmlparms(treecolor, namecolor, vrmlskycolornear,
                        vrmlskycolorfar, vrmlgroundcolornear);
      if (n != -1)
        getvrmlparms(&treecolor, &namecolor, &vrmlskycolornear,
               &vrmlskycolorfar, &vrmlgroundcolornear, &vrmlgroundcolorfar, n);
      countup(&loopcount, 10);
    } while (n != -1);
    break;
#endif

  case pict:
    strcpy(fontname,"Times");
    penchange = yes;
    xunitspercm = 28.346456693;
    yunitspercm = 28.346456693;
    /*7.5 x 10 inch default PICT page size*/
    xsize = 19.05;
    ysize = 25.40;
    break;

  case lw:
    penchange = yes;
    xunitspercm = 28.346456693;
    yunitspercm = 28.346456693;
    xsize = pagex;
    ysize = pagey;
    break;

  case idraw:
    penchange = yes;
    xunitspercm = 28.346456693;
    yunitspercm = 28.346456693;
    xsize = 21.59;
    ysize = 27.94;
    break;

  case hp:
    penchange = no;
    xunitspercm = 400.0;
    yunitspercm = 400.0;
    xsize = 24.0;
    ysize = 18.0;
    break;

#ifdef X
  case xpreview:
    xunitspercm = 39.37;
    yunitspercm = 39.37;
    xsize = width * 0.0254;
    ysize = height * 0.0254;
    break;
#endif

#ifdef WIN32
  case winpreview:
    penchange = yes;
    xunitspercm = 28.346456693;
    yunitspercm = 28.346456693;
    xsize = winwidth / xunitspercm;
    ysize = winheight / yunitspercm;
    break;
#endif

  case tek:
    xunitspercm = 50.0;
    yunitspercm = 50.0;
    xsize = 20.46;
    ysize = 15.6;
    break;

  case ibm:
#ifdef TURBOC
    GraphDriver = 0;
    detectgraph(&GraphDriver,&GraphMode);
    getmoderange(GraphDriver,&LoMode,&HiMode);
    initgraph(&GraphDriver,&HiMode,"");
    xunitspercm = getmaxx()/25;
    yunitspercm = getmaxy() / 17.5;
    restorecrtmode();
    xsize = 25.0;
    ysize = 17.5;
#endif
#ifdef QUICKC
    setupgraphics();
#endif
    break;

  case mac:
    penchange = yes;
    penchange = yes;
    xunitspercm = 28.346456693;
    yunitspercm = 28.346456693;
    xsize = winwidth / xunitspercm;
    ysize = winheight / yunitspercm;
    break;

  case houston:
    penchange = yes;
    xunitspercm = 100.0;
    yunitspercm = 100.0;
    xsize = 24.5;
    ysize = 17.5;
    break;

  case decregis:
    xunitspercm = 30.0;
    yunitspercm = 30.0;
    xsize = 25.0;
    ysize = 15.0;
    break;

  case epson:
    penchange = yes;
    xunitspercm = 47.244;
    yunitspercm = 28.346;
    xsize = 18.70;
    ysize = 22.0;
    strpwide = 960;
    strpdeep = 8;
    strpdiv = 1;
    break;

  case oki:
    penchange = yes;
    xunitspercm = 56.692;
    yunitspercm = 28.346;
    xsize = 19.0;
    ysize = 22.0;
    strpwide = 1100;
    strpdeep = 8;
    strpdiv = 1;
    break;

  case citoh:
    penchange = yes;
    xunitspercm = 28.346;
    yunitspercm = 28.346;
    xsize = 22.3;
    ysize = 26.0;
    strpwide = 640;
    strpdeep = 8;
    strpdiv = 1;
    break;

  case toshiba:
    penchange = yes;
    xunitspercm = 70.866;
    yunitspercm = 70.866;
    xsize = 19.0;
    ysize = 25.0;
    strpwide = 1350;
    strpdeep = 24;
    strpdiv = 4;
    break;

  case pcl:
    penchange = yes;
    xsize = 21.59;
    ysize = 27.94;
    xunitspercm = 118.11023622;   /* 300 DPI = 118.1 DPC                    */
    yunitspercm = 118.11023622;
    strpwide = 2550;   /* 8.5 * 300 DPI                                     */
    strpdeep = DEFAULT_STRIPE_HEIGHT;     /* height of the strip            */
    strpdiv = DEFAULT_STRIPE_HEIGHT;      /* in this case == strpdeep       */
                       /* this is information for 300 DPI resolution        */
    switch (hpresolution) {

    case 75:
      strpwide    /= 4;
      xunitspercm /= 4.0;
      yunitspercm /= 4.0;
      break;

    case 150:
      strpwide     /= 2;
      xunitspercm /= 2.0;
      yunitspercm /= 2.0;
      break;

    case 300:
      break;
    }
    break;

  case bmp:            /* since it's resolution dependent, make 1x1 pixels  */
    penchange = yes;   /* per square cm for easier math.                    */
    xunitspercm = 1.0;
    yunitspercm = 1.0;
    xsize = userxsize / xunitspercm; 
    ysize = userysize / yunitspercm;
    strpdeep = DEFAULT_STRIPE_HEIGHT;
    strpdiv = DEFAULT_STRIPE_HEIGHT;
    strpwide = (long)(xsize * xunitspercm);
    break;

  case xbm:            /* since it's resolution dependent, make 1x1 pixels  */
    penchange = yes;   /* per square cm for easier math.                    */
    xunitspercm = 1.0;
    yunitspercm = 1.0;
    xsize = userxsize / xunitspercm;
    ysize = userysize / yunitspercm;
    strpdeep = 10;
    strpdiv = 10;
    strpwide = (long)(xsize*xunitspercm);
    break;

  case pcx:
    penchange = yes;
    xsize = 21.16;
    ysize = 15.88;
    strpdeep = 10;
    strpdiv = 10;
    xunitspercm = strpwide / xsize;

    switch (resopts) {

    case 1:
      strpwide = 640;
      yunitspercm = 350 / ysize;
      break;

    case 2:
      strpwide = 800;
      yunitspercm = 600 / ysize;
      break;
      
    case 3:
      strpwide = 1024;
      yunitspercm = 768 / ysize;
      break;
    }
    break;

  case fig:
    penchange = yes;
    xunitspercm = 31.011;
    yunitspercm = 29.78;
    xsize = 25.4;
    ysize = 20.32;
    break;

  case gif:
  case other:
    break;
  default:
    break;
    /* initial parameter settings for a new plotter go here */
  }

  if (xsizehold != 0.0 && ysizehold != 0.0) {
    xmargin = xmargin * xsize / xsizehold;
    ymargin = ymargin * ysize / ysizehold;
  }
  
  if (previewing)
    return;
}  /* plotrparms */


void getplotter()
{
  long loopcount;
  Char ch,input[100];

  clearit() ;
  printf("\nWhich plotter or printer will the tree be drawn on?\n");
  printf("(many other brands or models are compatible with these)\n\n");
  printf("   type:       to choose one compatible with:\n\n");
  printf("        L         Postscript printer file format\n");
  printf("        M         PICT format (for drawing programs)\n");
  printf("        J         HP Laserjet PCL file format\n");
  printf("        W         MS-Windows Bitmap\n");
#ifdef DOS
  printf("        I         IBM PC graphics screens\n");
#endif
  printf("        F         FIG 2.0 drawing program format          \n");
  printf("        A         Idraw drawing program format            \n");
  printf("        Z         VRML Virtual Reality Markup Language file\n");
  printf("        P         PCX file format (for drawing programs)\n");
  printf("        K         TeKtronix 4010 graphics terminal\n");
  printf("        X         X Bitmap format\n");
  printf("        V         POVRAY 3D rendering program file\n");
  printf("        R         Rayshade 3D rendering program file\n");
  printf("        H         Hewlett-Packard pen plotter (HPGL file format)\n");
  printf("        D         DEC ReGIS graphics (VT240 terminal)\n");
  printf("        E         Epson MX-80 dot-matrix printer\n");
  printf("        C         Prowriter/Imagewriter dot-matrix printer\n");
  printf("        T         Toshiba 24-pin dot-matrix printer\n");
  printf("        O         Okidata dot-matrix printer\n");
  printf("        B         Houston Instruments plotter\n");
  printf("        U         other: one you have inserted code for\n");
  loopcount = 0;
  do {
    printf(" Choose one: \n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    fflush(stdout);
    scanf("%c%*[^\n]", &ch);
    getchar();
    uppercase(&ch);
    countup(&loopcount, 10);
  }
#ifdef DOS
while (strchr("LJKHIDBECOTUAZPXRMFWV",ch) == NULL);   
#endif
#ifndef DOS
while (strchr("LJKHDBECOTAZUPXRMFWV",ch) == NULL);
#endif
  switch (ch) {

  case 'L':
    plotter = lw;
    strcpy(fontname, "Times-Roman");
    break;

  case 'A':
    plotter = idraw;
    strcpy(fontname, "Times-Bold");
    break;

  case 'M':
    plotter = pict;
    strcpy(fontname, "Times");
    break;

  case 'R':
    plotter = ray;
    strcpy(fontname, "Hershey");
    break;

  case 'V':
    plotter = pov;
    strcpy(fontname, "Hershey");
    break;

  case 'Z':
    plotter = vrml;
    strcpy(fontname, "Hershey");
    treecolor = 5;
    namecolor = 4;
    vrmlskycolornear = 6;
    vrmlskycolorfar = 6;
    vrmlgroundcolornear = 3;
    vrmlgroundcolorfar = 3;
    break;

  case 'J':
    plotter = pcl;
    strcpy(fontname, "Hershey");
    printf("Please select Laserjet resolution\n\n");
    printf("1:  75 DPI\n2:  150 DPI\n3:  300 DPI\n\n");
    loopcount = 0;
    do {
#ifdef WIN32
      phyFillScreenColor();
#endif
      getstryng(input);
      ch = atoi(input);
      countup(&loopcount, 10);
    } while (ch != 1 && ch != 2 && ch != 3);
    hpresolution = 75*(1<<(ch-1));
    /* following pcl init code copied here from plotrparms                  */
    xunitspercm = 118.11023622;   /* 300 DPI = 118.1 DPC                    */
    yunitspercm = 118.11023622;
    strpwide = 2550;   /* 8.5 * 300 DPI                                     */
    strpdeep = DEFAULT_STRIPE_HEIGHT;     /* height of the strip            */
    strpdiv = DEFAULT_STRIPE_HEIGHT;      /* in this case == strpdeep       */
                       /* this is information for 300 DPI resolution        */
    switch (hpresolution) {

    case 75:
      strpwide    /= 4;
      xunitspercm /= 4.0;
      yunitspercm /= 4.0;
      break;

    case 150:
      strpwide     /= 2;
      xunitspercm /= 2.0;
      yunitspercm /= 2.0;
      break;

    case 300:
      break;
    }

    break;

  case 'K':
    plotter = tek;
    strcpy(fontname, "Hershey");
    break;

  case 'H':
    plotter = hp;
    strcpy(fontname, "Hershey");
    break;

  case 'I':
    plotter = ibm;
    strcpy(fontname, "Hershey");
    break;

  case 'D':
    plotter = decregis;
    strcpy(fontname, "Hershey");
    break;

  case 'B':
    plotter = houston;
    strcpy(fontname, "Hershey");
    break;

  case 'E':
    plotter = epson;
    strcpy(fontname, "Hershey");
    break;

  case 'C':
    plotter = citoh;
    strcpy(fontname, "Hershey");
    break;

  case 'O':
    plotter = oki;
    strcpy(fontname, "Hershey");
    break;

  case 'T':
    plotter = toshiba;
    strcpy(fontname, "Hershey");
    break;

  case 'P':
    plotter = pcx;
    strcpy(fontname, "Hershey");
    printf("Please select the PCX file resolution\n\n");
    printf("1: EGA 640  X 350\n");
    printf("2: VGA 800  X 600\n");
    printf("3: VGA 1024 X 768\n\n");
    loopcount = 0;
    do {
#ifdef WIN32
      phyFillScreenColor();
#endif
      getstryng(input);
      ch = (char)atoi(input);
      uppercase(&ch);
      countup(&loopcount, 10);
    } while (ch != 1 && ch != 2 && ch != 3);
    switch (ch) {
      
    case 1:
      strpwide = 640;
      yunitspercm = 350 / ysize;
      resopts = 1;
      break;
      
    case 2:
      strpwide = 800;
      yunitspercm = 600 / ysize;
      resopts = 2;
      break;
      
    case 3:
      strpwide = 1024;
      yunitspercm = 768 / ysize;
      resopts = 3;
      break;
    }
    break;

  case 'W':
    plotter = bmp;
    strcpy(fontname, "Hershey");
    printf("Please select the MS-Windows bitmap file resolution\n");
    printf("X resolution?\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    fflush(stdout);
    scanf("%lf%*[^\n]", &userxsize);
    getchar();
    printf("Y resolution?\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    fflush(stdout);
    scanf("%lf%*[^\n]", &userysize);
    getchar();
    xunitspercm = 1.0;
    yunitspercm = 1.0;
    /* Assuming existing reasonable margin values, set the margins
       to be the same as those in the previous output mode/resolution.
       This corrects the problem of the tree being hard up against the border
       when large resolutions are entered. */
    xmargin = userxsize / xsize * xmargin;
    ymargin = userysize / ysize * ymargin;

    xsize = userxsize;
    ysize = userysize;

    strpdeep = DEFAULT_STRIPE_HEIGHT;
    strpdiv = DEFAULT_STRIPE_HEIGHT;
    strpwide = (long)xsize;
    break;

  case 'X':
    plotter = xbm;
    strcpy(fontname, "Hershey");
    printf("Please select the X-bitmap file resolution\n");
    printf("X resolution?\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    fflush(stdout);
    scanf("%lf%*[^\n]", &userxsize);
    getchar();
    printf("Y resolution?\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    fflush(stdout);
    scanf("%lf%*[^\n]", &userysize);
    getchar();
    xunitspercm = 1.0;
    yunitspercm = 1.0;
    /* Assuming existing reasonable margin values, set the margins
       to be the same as those in the previous output mode/resolution. 
       This corrects the problem of the tree being hard up against the border
       when large resolutions are entered. */
    xmargin = userxsize / xsize * xmargin;
    ymargin = userysize / ysize * ymargin;

    xsize = userxsize;
    ysize = userysize;
    strpdeep = DEFAULT_STRIPE_HEIGHT;
    strpdiv = DEFAULT_STRIPE_HEIGHT;
    strpwide = (long)xsize;
    break;

  case 'F':
    plotter = fig;
    strcpy(fontname, "Times-Roman");
    break;

  case 'U':
    plotter = other;
    break;
  }
  dotmatrix = (plotter == epson || plotter == oki || plotter == citoh ||
               plotter == toshiba || plotter == pcx || plotter == pcl ||
               plotter == xbm || plotter == bmp);
}  /* getplotter */


void changepen(pentype pen)
{
  Char picthi, pictlo;
  long  pictint;
  lastpen = pen;

 switch (pen) {

  case treepen:
    linewidth = treeline;
    if (plotter == hp)
      fprintf(plotfile, "SP1;\n");
    if (plotter == lw) {
      fprintf(plotfile, "stroke %8.2f setlinewidth \n", treeline);
      fprintf(plotfile, " 1 setlinecap 1 setlinejoin \n");
    }
#ifdef WIN32
    if (plotter == winpreview)
      SelectObject(hdc, hPenTree); 
#endif
    break;

  case labelpen:
    linewidth = labelline;
    if (plotter == hp)
      fprintf(plotfile, "SP2;\n");
    if (plotter == lw) {
      fprintf(plotfile, " stroke%8.2f setlinewidth \n", labelline);
      fprintf(plotfile, "1 setlinecap 1 setlinejoin \n");
    }
#ifdef WIN32
    if (plotter == winpreview)
      SelectObject(hdc, hPenLabel); 
#endif
    break;
  }
#ifdef MAC
if (plotter == mac){
      pictint = ( long)(linewidth + 0.5);
      if (pictint ==0)
           pictint = 1;
}
#endif

  if (plotter != pict)
    return;
  pictint = ( long)(linewidth + 0.5);
  if (pictint == 0)
    pictint = 1;
  picthi = (Char)(pictint / 256);
  pictlo = (Char)(pictint & 255);
  fprintf(plotfile, "\007%c%c%c%c", picthi, pictlo, picthi, pictlo);
  bytewrite += 5;
}  /* changepen */


int readafmfile(char *filename, short *metric)
{
char line[256], word1[100], word2[100];
int scanned = 1, nmetrics=0, inmetrics, charnum, charlen, i, capheight=0;
FILE *fp;

fp = fopen(filename,"r");
if (!fp)
  return 0;
inmetrics = 0;

for (i=0;i<256;i++){
  metric[i] = (short)0;
}


for (;;){
  scan_eoln(fp);
  if (scanned != 1 )
    break;
  scanned=sscanf(line,"%s %s",word1,word2);
  if (scanned == 2 && strcmp(word1,"CapHeight") == 0)
    capheight = atoi(word2);

  if (inmetrics){
    sscanf(line,"%*s %s %*s %*s %s",word1,word2);
    charnum = atoi(word1);
    charlen = atoi(word2);
    nmetrics--;
    if (nmetrics == 0)
      break;

    if (charnum != -1 && charnum >= 32)
      metric[charnum-31] = charlen;
  }
  else
    if (scanned == 2 && strcmp(word1,"StartCharMetrics") == 0)
      nmetrics = atoi(word2),
      inmetrics = 1;

  if ((strcmp(word1,"EndCharMetrics") == 0) || (feof(fp)))
    break;
}
FClose(fp);
metric[0] = capheight;
return 1;
}  /* readafmfile */


void metricforfont(char *fontname, short *fontmetric)
{
  int i;
  long loopcount;
  char afmfile[FNMLNGTH];
  
  if ((strcmp(fontname,"Helvetica") == 0) ||
      (strcmp(fontname,"Helvetica-Oblique") == 0))
    for (i=31;i<256;++i)
      fontmetric[i-31] = helvetica_metric[i-31];
  else if ((strcmp(fontname,"Helvetica-Bold") == 0) ||
      (strcmp(fontname,"Helvetica-BoldOblique") == 0))
    for (i=31;i<256;++i)
      fontmetric[i-31] = helveticabold_metric[i-31];
  else if (strcmp(fontname,"Times-Roman") == 0)
    for (i=31;i<256;++i)
      fontmetric[i-31] = timesroman_metric[i-31];
  else if (strcmp(fontname,"Times") == 0)
    for (i=31;i<256;++i)
      fontmetric[i-31] = timesroman_metric[i-31];
  else if (strcmp(fontname,"Times-Italic") == 0)
    for (i=31;i<256;++i)
      fontmetric[i-31] = timesitalic_metric[i-31];
  else if (strcmp(fontname,"Times-Bold") == 0)
    for (i=31;i<256;++i)
      fontmetric[i-31] = timesbold_metric[i-31];
  else if (strcmp(fontname,"Times-BoldItalic") == 0)
    for (i=31;i<256;++i)
      fontmetric[i-31] = timesbolditalic_metric[i-31];
  
  else if (strncmp(fontname,"Courier",7) == 0){
    fontmetric[0] = 562;  
    for (i=32;i<256;++i)
      fontmetric[i-31] = (short)600;
  }
  else {
    if (didloadmetric){
      for (i=31;i<256;++i)
        fontmetric[i-31] = unknown_metric[i-31];}
    else {
      didloadmetric = 1;
      sprintf(afmfile,"%s.afm",fontname);    /* search current dir */
      if (readafmfile(afmfile,unknown_metric)){
        for (i=31;i<256;++i)
          fontmetric[i-31] = unknown_metric[i-31];
        return;}
      sprintf(afmfile,"%s%s.afm",AFMDIR,fontname); /* search afm dir */
      if (readafmfile(afmfile,unknown_metric)){
      for (i=31;i<256;++i)
        fontmetric[i-31] = unknown_metric[i-31];
      return;}
#ifdef NeXT
      sprintf(afmfile,"%s/Library/Fonts/%s.font/%s.afm",getenv("HOME"),
              fontname,fontname);
      if (readafmfile(afmfile,unknown_metric)){
      for (i=31;i<256;++i)
        fontmetric[i-31] = unknown_metric[i-31];
      return;}
      sprintf(afmfile,"/LocalLibrary/Fonts/%s.font/%s.afm",fontname,fontname);
      if (readafmfile(afmfile,unknown_metric)){
      for (i=31;i<256;++i)
        fontmetric[i-31] = unknown_metric[i-31];
      return;}
#endif
      loopcount = 0;
      for (;;){
        printf("Enter the path of the %s.afm file, or \"none\" for best guess:",
                fontname);
        getstryng(afmfile);
        if (strcmp(afmfile,"none") == 0){
          for (i=31;i<256;++i)
            fontmetric[i-31] = timesroman_metric[i-31],
            unknown_metric[i-31] = timesroman_metric[i-31],
            didloadmetric =1;
          return;
        }
        else {
          if (readafmfile(afmfile,unknown_metric)){
            for (i=31;i<256;++i)
              fontmetric[i-31] = unknown_metric[i-31];
            return;}
          else
            printf("Can't read that file. Please re-enter.\n");
        }
        countup(&loopcount, 10);
      }
    }
  }
}  /* metricforfont */


double heighttext(fonttype font, char *fontname)
{
  short          afmetric[256];
#ifdef MAC
  FontInfo       info;
#endif

  if (strcmp(fontname,"Hershey") == 0)
    return (double)font[2];
#ifdef MAC
  else if (((plotter == pict || plotter == mac)  &&
            (((grows == vertical && labelrotation == 0.0) ||
             (grows == horizontal && labelrotation == 90.0))))){
    TextFont(macfontid(fontname)); 
    TextSize((int)(1000));
    TextFace((int)((pictbold ? 1: 0) | (pictitalic ? 2 : 0)|
        (pictoutline ? 8 : 0)|(pictshadow ? 16 : 0)));
    GetFontInfo(&info);   
    TextFont(macfontid("courier"));
    TextSize(10);
    TextFace(0);
    return (double)info.ascent;
       }
#endif
  else if (strcmp(fontname,"Hershey") == 0)
    return (double)font[2];
  else{
      metricforfont(fontname,afmetric);  
      return (double)afmetric[0];}
}  /* heighttext */


double lengthtext(char *pstring, long nchars, char *fontname,
                        fonttype font)
{  /* lengthext */
  long i, j, code;
  static double  sumlength;
  long                  sumbigunits;
  short          afmetric[256];

  sumlength = 0.0;
  if (strcmp(fontname,"Hershey") == 0) {
    for (i = 0; i < nchars; i++) {
      code = pstring[i];
      j = 1;
      while (font[j] != code && font[j - 1] != 0)
        j = font[j - 1];
      if (font[j] == code)
        sumlength += font[j + 2];
    }
    return sumlength;
  }
#ifdef MAC
  else   if  (((plotter == pict || plotter == mac) && 
          (((grows == vertical && labelrotation == 0.0) ||
          (grows == horizontal && labelrotation == 90.0))))){
  TextFont(macfontid(fontname)); 
  TextSize((int)(1000));
  TextFace((int)((pictbold ? 1: 0) | (pictitalic ? 2 : 0)|
      (pictoutline ? 8 : 0)|(pictshadow ? 16 : 0)));
   sumbigunits = 0;
    for (i = 0; i < nchars; i++) 
      sumbigunits += (long)CharWidth(pstring[i]);
      TextFace(0);
      TextSize(10);
      TextFont(macfontid("courier"));
      return (double)sumbigunits;
          }
#endif
      else {
       metricforfont(fontname,afmetric);
       sumbigunits = 0;
       for (i = 0; i < nchars; i++) 
         sumbigunits += afmetric[(int)(1+(unsigned char)pstring[i] - 32)];
    
       sumlength = (double)sumbigunits;
          }
      return sumlength;
}  /* lengthtext */


void plotchar(long *place, struct LOC_plottext *text)
{
  text->heightfont = text->font[*place + 1];
  text->yfactor = text->height / text->heightfont;
  text->xfactor = text->yfactor;
  *place += 3;
  do {
    (*place)++;
    text->coord = text->font[*place - 1];
    if (text->coord > 0)
      text->penstatus = pendown;
    else
      text->penstatus = penup;
    text->coord = abs(text->coord);
    text->coord %= 10000;
    text->xfont = (text->coord / 100 - xstart) * text->xfactor;
    text->yfont = (text->coord % 100 - ystart) * text->yfactor;
    text->xplot = text->xx + (text->xfont * text->cosslope +
                              text->yfont * text->sinslope) * text->compress;
    text->yplot = text->yy - text->xfont * text->sinslope +
      text->yfont * text->cosslope;
    plot(text->penstatus, text->xplot, text->yplot);
  } while (abs(text->font[*place - 1]) < 10000);
  text->xx = text->xplot;
  text->yy = text->yplot;
}  /* plotchar */


void swap_charptr(char **one, char **two)
{
  char *tmp = (*one);
  (*one)= (*two);
  (*two) = tmp;
}  /* swap */


void plotpb()
{
  pagecount++;
  fprintf(plotfile,"\n showpage \n%%%%PageTrailer\n");
  fprintf(plotfile,"%%%%DocumentFonts: %s\n",
              (strcmp(fontname,"Hershey") == 0) ? "" : fontname);
  fprintf(plotfile,"%%%%\n%%%%Page: %ld %ld\n",pagecount,pagecount);
  fprintf(plotfile,"%%%%PageBoundingBox: 0 0 %d %d\n",
          (int)(xunitspercm*paperx),(int)(yunitspercm*papery));
  fprintf(plotfile,"%%%%PageFonts: (atend)\n%%%%BeginPageSetup\n%%%%PaperSize: Letter\n");
  fprintf(plotfile,"0 0 moveto\n"); /* hack to make changepen work w/o errors */
  changepen(lastpen);
}  /* plotpb */


void drawit(char *fontname, double *xoffset, double *yoffset,
                        long numlines, node *root)
{
  long i, j, line, xpag, ypag;

  long test_long ;  /* To get a division out of a loop */

  (*xoffset) = 0.0;
  (*yoffset) = 0.0;

  xpag = (int)((pagex-hpmargin-0.01)/(paperx - hpmargin))+1;
  ypag = (int)((pagey-vpmargin-0.01)/(papery - vpmargin))+1;
  if (dotmatrix){
    strptop    = (long)(ysize * yunitspercm);
    strpbottom = numlines*strpdeep + 1;
  }
  else {
    pagecount = 1;
    for (j=0; j<ypag; ++j)
      for (i=0; i<xpag; ++i){
        clipx0 = (double)i*(paperx -  hpmargin);
        clipx1 = (double)(i*(paperx -  hpmargin))+(paperx - hpmargin);
        clipy0 = (double)(j*(papery -  vpmargin));
        clipy1 = (double)(j*(papery-hpmargin))+(papery+vpmargin);
        plottree(root, root);
        plotlabels(fontname);
        if (!(i == xpag - 1 && j == ypag - 1) && plotter == lw)
          plotpb(); /* page break */
      }
  }

  if (dotmatrix){
    striprint(( long)((ysize * yunitspercm)- (numlines * strpdeep)),
              ( long)((ysize * yunitspercm)- (numlines * strpdeep)));
    strptop = numlines * strpdeep;
    strpbottom = strptop - strpdeep + 1;
    printf(" writing %3ld lines ...\n", numlines);
    printf("  Line     Output file size\n");
    printf("  ----     ------ ---- ----\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    test_long = strpwide / 8 ;        /* A fix for below */

    for (line = 1; line <= numlines ; line++) {
      for (i = 0; i <= strpdeep ; i++) {
        for (j=0; j<=test_long;++j)   /* Don't do a divide every time! */
          stripe[i][j] = 0;
      }
      empty = true;
      xnow = strpwide / 2.0;
      ynow = 0.0;
      plottree(root, root);
      plotlabels(fontname);
      strptop     = strpbottom - 1;
      strpbottom -= strpdeep;

      if (strpdeep > DEFAULT_STRIPE_HEIGHT){
        /* large stripe, do in DEFAULT_STRIPE_HEIGHT (20)-line     */
        for (i=0;i<strpdeep;++i){
          swap_charptr(&stripe[i%DEFAULT_STRIPE_HEIGHT],
               &stripe[i]);
          if ((i%DEFAULT_STRIPE_HEIGHT) ==
              (DEFAULT_STRIPE_HEIGHT -1)){
            striprint(DEFAULT_STRIPE_HEIGHT,
                      DEFAULT_STRIPE_HEIGHT);}
        }
        striprint(strpdeep%DEFAULT_STRIPE_HEIGHT,
                  strpdeep%DEFAULT_STRIPE_HEIGHT);
      }
      else{                          /* small stripe, do it all now.     */
        striprint(strpdiv,strpdeep);
        if (line % 5 == 0)
          printf("%5ld%16ld\n", line, filesize);
#ifdef WIN32
          phyFillScreenColor();
#endif

      }
    }
  }
}  /* drawit */


char *findXfont(char *fontname, double pointsize, double *scale,
                        int *epointsize)
{
  static char returnval[64];

  if (strcmp(fontname,"Helvetica") == 0)
    strcpy(returnval,"*-helvetica-medium-r-*-120-75-75-*"),
    (*scale) = pointsize / 12.0,
    (*epointsize) = 12;
  else if (strcmp(fontname,"Helvetica-Oblique") == 0)
    strcpy(returnval,"*-helvetica-medium-o-*-140-75-75-*"),
    (*scale) = pointsize / 14.0,
    (*epointsize) = 14;
  else if (strcmp(fontname,"Helvetica-Bold") == 0)
    strcpy(returnval,"*-helvetica-bold-r-*-140-75-75-*"),
    (*scale) = pointsize / 14.0,
    (*epointsize) = 14;
  else if (strcmp(fontname,"Helvetica-BoldOblique") == 0)
    strcpy(returnval,"*-helvetica-medium-o-*-140-75-75-*"),
    (*scale) = pointsize / 14.0,
    (*epointsize) = 14;
  else if (strcmp(fontname,"Times-Roman") == 0)
    strcpy(returnval,"*-times-medium-r-*-140-75-75-*"),
    (*scale) = pointsize / 14.0,
    (*epointsize) = 14;
  else if (strcmp(fontname,"Times-Italic") == 0)
    strcpy(returnval,"*-times-medium-i-*-140-75-75-*"),
    (*scale) = pointsize / 14.0,
    (*epointsize) = 14;
  else if (strcmp(fontname,"Times-Bold") == 0)
    strcpy(returnval,"*-times-medium-i-*-140-75-75-*"),
    (*scale) = pointsize / 14.0,
    (*epointsize) = 14;
  else if (strcmp(fontname,"Times-BoldItalic") == 0)
    strcpy(returnval,"*-times-medium-i-*-140-75-75-*"),
    (*scale) = pointsize / 14.0,
    (*epointsize) = 14;
  else if (strcmp(fontname,"Courier") == 0)
    sprintf(returnval,"*-courier-medium-r-*-100-75-75-*"),
    (*scale) = pointsize / 12.0,
    (*epointsize) = 12;
  else if (strcmp(fontname,"Courier-Italic") == 0)
    strcpy(returnval,"*-courier-medium-r-*-120-75-75-*"),
    (*scale) = pointsize / 12.0,
    (*epointsize) = 12;
  else if (strcmp(fontname,"Courier-Bold") == 0)
    strcpy(returnval,"*-courier-bold-r-*-120-75-75-*"),
    (*scale) = pointsize / 12.0,
    (*epointsize) = 12;
  else if (strcmp(fontname,"Courier-BoldItalic") == 0)
    strcpy(returnval,"*-courier-bold-r-*-120-75-75-*"),
    (*scale) = pointsize / 12.0,
    (*epointsize) = 12;
  else
    sprintf(returnval,"*-times-medium-r-*-120-75-75-*"),
    (*scale) = pointsize / 12.0,
    (*epointsize) = 12;
  return returnval;
}  /* findXfont */


int macfontid(char *fontname)
{
char fontnam[256];
int  i;

  strcpy(fontnam,fontname);
  for (i=0;i<strlen(fontnam);++i)
    fontnam[i] = toupper(fontnam[i]);
  if (strcmp(fontnam,"NEW YORK") == 0)
    return 2;
  else if (strcmp(fontnam,"GENEVA") == 0)
    return 3;
  else if (strcmp(fontnam,"MONACO") == 0)
    return 4;
  else if (strcmp(fontnam,"VENICE") == 0)
    return 5;
  else if (strcmp(fontnam,"LONDON") == 0)
    return 6; 
  else if (strcmp(fontnam,"ATHENS") == 0)
    return 7;    
  else if (strcmp(fontnam,"SAN FRANCISCO") == 0)
    return 8; 
  else if (strcmp(fontnam,"TORONTO") == 0)
    return 9; 
  else if (strcmp(fontnam,"CAIRO") == 0)
    return 11; 
  else if (strcmp(fontnam,"LOS ANGELES") == 0)
    return 12; 
  else if (strcmp(fontnam,"TIMES") == 0)
    return 20; 
  else if (strcmp(fontnam,"TIMES-ROMAN") == 0)
    return 20; 
  else if (strcmp(fontnam,"HELVETICA") == 0)
    return 21;   
  else if (strcmp(fontnam,"COURIER") == 0)
    return 22;   
  else if (strcmp(fontnam,"SYMBOL") == 0)
    return 23; 
  else if (strcmp(fontnam,"TALIESIN") == 0)
    return 24;   
  else
    return 0;  
 }  /* macfontid */


void plottext(Char *pstring,long nchars,double height_,double cmpress2,
                        double x,double y,double slope,short *font_,char *fontname)
{
#ifdef max
#undef max
#endif
#define max(a,b) ((a > b) ? a : b)
#ifdef min
#undef min
#endif
#define min(a,b) ((a > b) ? b : a)
#define max4(a,b,c,d) (max(max(a,b),max(c,d)))
#define min4(a,b,c,d) (min(min(a,b),min(c,d)))
  struct LOC_plottext text;
  long i, j, code;
  double pointsize;
  int    epointsize;   /* effective pointsize before scale in idraw matrix */
  double iscale;
  double textlen;
  double px0,py0,px1,py1; /* square bounding box of text */
  
  text.heightfont = font_[2];
  pointsize = (((height_ / xunitspercm) / 2.54) * 72.0);
  
  if (strcmp(fontname,"Hershey") !=0)
        pointsize *= ((double)1000.0 / heighttext(font_,fontname));

  text.height = height_;
  text.compress = cmpress2;
  text.font = font_;
  text.xx = x;
  text.yy = y;
  text.sinslope = sin(pi * slope / 180.0);
  text.cosslope = cos(pi * slope / 180.0);

  if ((strcmp(fontname,"Hershey") == 0)||
      (previewing && (!(((plotter == pict) || (plotter == mac))
                        && (((grows == vertical) && (labelrotation == 0.0)) ||
                            ((grows == horizontal) && (labelrotation == 90.0))
                            ))))){
    for (i = 0; i < nchars; i++) {
      code = pstring[i];
      j = 1;
      while (text.font[j] != code && text.font[j - 1] != 0)
        j = text.font[j - 1];
      plotchar(&j,  &text);
    }
  }
 /* print native font.  idraw, PS, pict, and fig. */

  else if (plotter == fig) {
    fprintf(plotfile,"4 0 %d %d 0 -1 0 %1.5f 4 19 163 %d %d %s\001\n",
            figfontid(fontname), /* font ID    */
            (int)pointsize,      /* font size  */
            (double)0.0,         /* font rotation */
            (int)x,              /* x position    */
            (int)(606.0 - y),    /* y position    */
            pstring);
  }
  else if (plotter == lw)  {
    /* If there's NO possibility that the line intersects the square bounding
     * box of the font, leave it out. Otherwise, let postscript clip to region.
     * Compute text boundary, be REAL generous. */
    textlen = (lengthtext(pstring,nchars,fontname,font_)/1000)*pointsize;
    px0 = min4(x + (text.cosslope * pointsize),
               x - (text.cosslope * pointsize),
               x + (text.cosslope * pointsize) + (text.sinslope * textlen),
               x - (text.cosslope * pointsize) + (text.sinslope * textlen))
      /28.346;
    px1 = max4(x + (text.cosslope * pointsize),
               x - (text.cosslope * pointsize),
               x + (text.cosslope * pointsize) + (text.sinslope * textlen),
               x - (text.cosslope * pointsize) + (text.sinslope * textlen))
      /28.346;
    py0 = min4(y + (text.sinslope * pointsize),
               y - (text.sinslope * pointsize),
               y + (text.sinslope * pointsize) + (text.cosslope * textlen),
               y - (text.sinslope * pointsize) + (text.cosslope * textlen))
      /28.346;
    py1 = max4(y + (text.sinslope * pointsize),
               y - (text.sinslope * pointsize),
               y + (text.sinslope * pointsize) + (text.cosslope * textlen),
               y - (text.sinslope * pointsize) + (text.cosslope * textlen))
      /28.346;

    /* if rectangles intersect, print it. */
    if (rectintersects(px0,py0,px1,py1,clipx0,clipy0,clipx1,clipy1)) {
      fprintf(plotfile,"gsave\n");
      fprintf(plotfile,"/%s findfont %f scalefont setfont\n",fontname,
              pointsize);
      fprintf(plotfile,"%f %f translate %f rotate\n",
                x-(clipx0*xunitspercm),y-(clipy0*xunitspercm),-slope);
      fprintf(plotfile,"0 0 moveto\n");
      fprintf(plotfile,"(%s) show\n",pstring);
      fprintf(plotfile,"grestore\n");
    }
  }
else if (plotter == idraw) {
   iscale = pointsize / 12.0;
   y += text.height * text.cosslope;
   x += text.height * text.sinslope;
        fprintf(plotfile, "Begin %%I Text\n");
   fprintf(plotfile, "%%I cfg Black\n");
   fprintf(plotfile, "0 0 0 SetCFg\n");
   fprintf(plotfile, "%%I f %s\n",
             findXfont(fontname,pointsize,&iscale,&epointsize));
   fprintf(plotfile,"%s %d SetF\n",fontname,epointsize);
   fprintf(plotfile, "%%I t\n");
   fprintf(plotfile, "[ %f %f %f %f %f %f ] concat\n",
           text.cosslope*iscale, -text.sinslope*iscale,
           text.sinslope*iscale, text.cosslope*iscale,
           x+216.0 ,y+285.0);
   fprintf(plotfile, "%%I\n");
   fprintf(plotfile, "[\n(%s)\n] Text\nEnd\n\n",pstring); 
   
 }
else if (plotter == pict || plotter == mac) {
   if (previewing){
#ifdef MAC
   TextFont(macfontid(fontname)); 
   TextSize((int)(pointsize+0.5));
   TextFace((int)((pictbold ? 1: 0) | (pictitalic ? 2 : 0)|
      (pictoutline ? 8 : 0)|(pictshadow ? 16 : 0)));
      MoveTo((int)floor((double)x + 0.5),
              winheight - (long)floor((double)y + 0.5)+MAC_OFFSET); 
      putstring(pstring);
    TextFont(macfontid("courier"));
    TextSize(10);
    TextFace(0);
#endif
   }
   else {
   /* txfont: */
    fprintf(plotfile,"%c",(unsigned char)3);
    pictoutint(plotfile,macfontid(fontname));
    /* txsize: */
    fprintf(plotfile,"%c",13);
    pictoutint(plotfile,(int)(pointsize+0.5));
    /* txface: */
    fprintf(plotfile,"%c%c",4,
      (int)((pictbold ? 1: 0) | (pictitalic ? 2 : 0)|
      (pictoutline ? 8 : 0)|(pictshadow ? 16 : 0)));
    /* txfloc: */    
    fprintf(plotfile,"%c",40);
    pictoutint(plotfile,(int)floor(ysize * yunitspercm - y + 0.5));
    pictoutint(plotfile,(int)(x+0.5));
    fprintf(plotfile,"%c%s",(char)strlen(pstring),pstring);
    bytewrite+=(14+strlen(pstring));  
   }
 }
}  /* plottext */


void makebox(char *fn,double *xo,double *yo,double *scale,long ntips)
/* fn--fontname| xo,yo--x and y offsets */
{
  /* draw the box on screen which represents plotting area.        */
  char ch;
  long xpag,ypag,i,j;
  double xpagecorrection, ypagecorrection;

  if (previewer != winpreview && previewer != mac && previewer != xpreview) {
    printf("\nWe now will preview the tree.  The box that will be\n");
    printf("plotted on the screen represents the boundary of the\n");
    printf("final plotting surface.  To see the preview, press on\n");
    printf("the ENTER or RETURN key (you may need to do it twice).\n");
    printf("When finished viewing it, press on that key again.\n");
  }
  oldpenchange   = penchange;
  oldxsize       = xsize;
  oldysize       = ysize;
  oldxunitspercm = xunitspercm;
  oldyunitspercm = yunitspercm;
  oldxcorner     = xcorner;
  oldycorner     = ycorner;
  oldxmargin     = xmargin;
  oldymargin     = ymargin;
  oldhpmargin    = hpmargin;
  oldvpmargin    = vpmargin;
  oldplotter     = plotter;
  plotter        = previewer;
  if (previewer != winpreview && previewer != mac && previewer != xpreview) {
#ifdef WIN32
    phyFillScreenColor();
#endif
    fflush(stdout);
    scanf("%c%*[^\n]", &ch);
    (void)getchar();
    if (ch == '\n')
      ch = ' ';
  }
  plotrparms(ntips);
  initplotter(ntips,fn);
  xcorner += 0.05 * xsize;
  ycorner += 0.05 * ysize;
  xsize *= 0.9;
  ysize *= 0.9;
  (*scale) = ysize / oldysize;
  if (xsize / oldxsize < (*scale))
    (*scale) = xsize / oldxsize;
  xpagecorrection = oldxsize / pagex;
  ypagecorrection = oldysize / pagey;
  (*xo) = (xcorner + (xsize - oldxsize * (*scale)) / 2.0) / (*scale);
  (*yo) = (ycorner   + (ysize - oldysize * (*scale)) / 2.0) / (*scale);

  xscale = (*scale) * xunitspercm;
  yscale = (*scale) * yunitspercm;
  xmargin *= (*scale);
  ymargin *= (*scale);
  hpmargin *= (*scale);
  vpmargin *= (*scale);
  xpag = (int)((pagex-hpmargin-0.01)/(paperx - hpmargin))+1;
  ypag = (int)((pagey-vpmargin-0.01)/(papery - vpmargin))+1;
  /* draw the outer borders */

  plot(penup, xscale * (*xo), yscale * (*yo));
  plot(pendown, xscale * (*xo), yscale * ((*yo) + pagey * ypagecorrection));
  plot(pendown, xscale * ((*xo) + pagex * xpagecorrection), 
    yscale * ((*yo) + pagey * ypagecorrection));
  plot(pendown, xscale * ((*xo) + pagex * xpagecorrection), yscale * (*yo));
  plot(pendown, xscale * (*xo), yscale * (*yo));
  /* we've done the extent, now draw the dividing lines: */
  for (i=0; i<xpag; ++i){
      plot(penup,(xscale * (*xo))+xscale*i*(paperx - hpmargin)*xpagecorrection,
                  ((*yo)*yscale));
      plot(pendown,(xscale * (*xo))+xscale*i*(paperx - hpmargin)*xpagecorrection,
                  ((*yo)*yscale)+yscale*pagey*ypagecorrection);
      if (i != 0) {
        plot(penup,(xscale * (*xo))
             +xscale*i*(paperx - hpmargin)*xpagecorrection+xscale*hpmargin,((*yo)*yscale));
        plot(pendown,(xscale * (*xo))
             +xscale*i*(paperx - hpmargin)*xpagecorrection+xscale*hpmargin,
             ((*yo)*yscale)+yscale*pagey*ypagecorrection);
      }
  }
  for (j=0;j<ypag;++j){
      plot(penup,(xscale * (*xo)),
           ((*yo)*yscale)+yscale*j*(papery-vpmargin)*ypagecorrection);
      plot(pendown,(xscale * (*xo))+xscale*pagex*xpagecorrection,
           ((*yo)*yscale)+yscale*j*(papery-hpmargin)*ypagecorrection);
      if (j != 0) {
        plot(penup,(xscale * (*xo)),
             ((*yo)*yscale)+yscale*j*(papery-vpmargin)*ypagecorrection+yscale*vpmargin);
        plot(pendown,(xscale * (*xo))+xscale*pagex*xpagecorrection,
             ((*yo)*yscale)+yscale*j*(papery-hpmargin)*ypagecorrection+yscale*vpmargin);
      }
  }
}  /* makebox */


boolean plotpreview(char *fn, double *xo, double *yo, double *scale,
                        long nt, node *root)
{
  long loopcount;
  boolean canbeplotted;
  Char ch;

  previewing = true;
#ifdef WIN32
  if (previewer == winpreview) {
    winpreviewparms.fn    = fn;
    winpreviewparms.xo    = xo;
    winpreviewparms.yo    = yo;
    winpreviewparms.scale = scale;
    winpreviewparms.nt    = nt;
    winpreviewparms.root  = root;
    winplotpreview();
  }
  else {
    makebox(fn,xo,yo,scale,nt);
    plottree(root, root);
    plotlabels(fn);
    finishplotter();
  }
#endif
#ifdef MAC 
  if (previewer == mac ){

    macpreviewparms.fn    = fn;
        macpreviewparms.xo    = xo;
        macpreviewparms.yo    = yo;
        macpreviewparms.scale = scale;
        macpreviewparms.nt    = nt;
        macpreviewparms.root  = root;

        /*oldplotter=plotter;
        plotter = previewer;
        initplotter(nt,fn);*/
        gfxmode();
        paint_gfx_window();
        oldplotter=plotter;
        plotter=mac;
        finishplotter();
    } else  {
        makebox(fn,xo,yo,scale,nt);
        plottree(root, root);
        plotlabels(fn);
        finishplotter();
    }
#endif
#ifdef X
    if (previewer == xpreview ){
        xpreviewparms.fn    = fn;
        xpreviewparms.xo    = xo;
        xpreviewparms.yo    = yo;
        xpreviewparms.scale = scale;
        xpreviewparms.nt    = nt;
        xpreviewparms.root  = root;

        init_x();
        oldplotter=plotter;
        plotter=xpreview;
        finishplotter();
    } else  {
        makebox(fn,xo,yo,scale,nt);
        plottree(root, root);
        plotlabels(fn);
        finishplotter();
    }

#endif
#ifndef MAC
#ifndef WIN32
#ifndef X
  makebox(fn,xo,yo,scale,nt);
  plottree(root, root);
  plotlabels(fn);
  finishplotter();
#endif
#endif
#endif
  penchange = oldpenchange;
  xsize = oldxsize;
  ysize = oldysize;
  xunitspercm = oldxunitspercm;
  yunitspercm = oldyunitspercm;
  xscale = xunitspercm;
  yscale = yunitspercm;
  plotter = oldplotter;
  xcorner = oldxcorner;
  ycorner = oldycorner;
  xmargin = oldxmargin;
  ymargin = oldymargin;
  hpmargin = oldhpmargin;
  vpmargin = oldvpmargin;
  if (previewer == winpreview || previewer == xpreview || previewer == mac) {
    canbeplotted = (winaction == plotnow);
  } else {
    printf(" Is the tree ready to be plotted? (Answer Y or N)\n");
    loopcount = 0;
    do {
      printf("Type Y or N:\n");
#ifdef WIN32
      phyFillScreenColor();
#endif
      fflush(stdout);
      scanf("%c%*[^\n]", &ch);
      (void)getchar();
      if (ch == '\n')
        ch = ' ';
      uppercase(&ch);
      countup(&loopcount, 10);
    } while (ch != 'Y' && ch != 'N');
    canbeplotted = (ch == 'Y');
  }
  if (previewer == xpreview) {
  }
  return canbeplotted;
}  /* plotpreview */


#ifdef WIN32
winplotpreviewcore(int windowwidth, int windowheight)
{
  winheight = windowheight;
  winwidth = windowwidth;
  makebox(winpreviewparms.fn,
          winpreviewparms.xo,
          winpreviewparms.yo,
          winpreviewparms.scale,
          winpreviewparms.nt);
  plottree(winpreviewparms.root, winpreviewparms.root);
  plotlabels(winpreviewparms.fn);
  finishplotter();
  penchange = oldpenchange;
  xsize = oldxsize;
  ysize = oldysize;
  xunitspercm = oldxunitspercm;
  yunitspercm = oldyunitspercm;
  xscale = xunitspercm;
  yscale = yunitspercm;
  plotter = oldplotter;
  xcorner = oldxcorner;
  ycorner = oldycorner;
  xmargin = oldxmargin;
  ymargin = oldymargin;
  hpmargin = oldhpmargin;
  vpmargin = oldvpmargin;
}
#endif


long allocstripe(striptype stripe, long x, long y)
{
  long i;
  for (i=0 ; i<=y; ++i){
    stripe[i] = (MALLOCRETURN *)Malloc((x+1)*sizeof(Char));
    if (!stripe[i])
      break;
  }
  return i-1; 
} /* allocstripe */


#ifdef X

void  redraw(Widget w,XtPointer client, XExposeEvent *ev) {
        XGCValues values;
        values.line_width=3;
        XChangeGC(display,gc1,GCLineWidth,&values);
        makebox(xpreviewparms.fn,
             xpreviewparms.xo,
             xpreviewparms.yo,
             xpreviewparms.scale,
                xpreviewparms.nt);
        values.line_width=1;
        XChangeGC(display,gc1,GCLineWidth,&values);
        plottree(xpreviewparms.root, xpreviewparms.root);
        plotlabels(xpreviewparms.fn);

        penchange = oldpenchange;
        xsize = oldxsize;
        ysize = oldysize;
        xunitspercm = oldxunitspercm;
        yunitspercm = oldyunitspercm;
        xscale = xunitspercm;
        yscale = yunitspercm;
        plotter = oldplotter;
        xcorner = oldxcorner;
        ycorner = oldycorner;
        xmargin = oldxmargin;
        ymargin = oldymargin;
        hpmargin = oldhpmargin;
        vpmargin = oldvpmargin;


}

void plot_callback(Widget w,XtPointer client, XtPointer call) {
        winaction=plotnow;
        close_x();
}


void change_callback(Widget w,XtPointer client, XtPointer call) {
        winaction=changeparms;
        close_x();
}

void quit_callback(Widget w,XtPointer client, XtPointer call) {
        winaction=quitnow;
        close_x();
}

void about_callback(Widget w,XtPointer client, XtPointer call) {
        do_dialog();
}
 
void delete_callback(Widget w, XEvent* event, String *params, int *num_params) {
        if (event->type != ClientMessage || event->xclient.data.l[0] !=
                        wm_delete_window) 
                return;
        winaction=changeparms;
        close_x();

}

void close_x() {
        shell=NULL;
        XtAppSetExitFlag(appcontext);
        XtUnrealizeWidget(toplevel);
        XtDestroyWidget(toplevel);
        XtCloseDisplay(display);
}


void dismiss_dialog()
{
        XtDestroyWidget(shell);        
        shell=NULL;
}

void do_dialog() {
        if (shell != NULL)
                return;        
        shell=XtCreatePopupShell("About",transientShellWidgetClass,
                        toplevel,NULL,0);
        dialog=XtCreateManagedWidget("dialog",dialogWidgetClass,shell,NULL,0);
        XawDialogAddButton(dialog,"Dismiss",(XtCallbackProc)dismiss_dialog
                        ,NULL);
        XtRealizeWidget(shell);
          wm_delete_window2 = XInternAtom(XtDisplay(shell),
                  "WM_DELETE_WINDOW",0);
        XSetWMProtocols(XtDisplay(shell),XtWindow(shell),
                          &wm_delete_window2,1);
        XtMapWidget(shell);
}

static XtActionsRec draw_actions[] = {
        { "quit", (XtActionProc)delete_callback },
};

void init_x() {
  Widget paned;
  Widget menubar;
  Widget menuButton;
  Widget menu;
  Widget entry;
  Widget drawing_area; 
  XSetWindowAttributes winAttr;
  Arg wargs[7];
  unsigned int dummy1,dummy2;
  Window dummy3;
  XGCValues values;

  toplevel=XtAppInitialize(&appcontext,"phylip",NULL,0,&nargc,nargv,res,
                  NULL,0); 
  
  /* make the top level window*/
                /* this is for closing the window*/
  XtAppAddActions(appcontext,draw_actions,1);
  XtOverrideTranslations(toplevel,
          XtParseTranslationTable ("<Message>WM_PROTOCOLS: quit()"));

  /* create a form add it to toplevel */
  paned = XtCreateManagedWidget("paned",formWidgetClass,toplevel,NULL,0);
  
  /* create a menubar add it to the form*/
  menubar = XtCreateManagedWidget("menubar",boxWidgetClass,paned,NULL,0);
  
  
  /* create an area to draw in  with a size relative to the size of the screen*/
  XGetGeometry(XtDisplay(toplevel),XDefaultRootWindow(XtDisplay(toplevel)),
                  &dummy3,&x,&y,&width,&height,&dummy1,&dummy2);

  height *= 0.7;
  width = 0.75 * height;
 
  XtSetArg(wargs[0],XtNwidth,width);
  XtSetArg(wargs[1],XtNheight,height); 
                  
  drawing_area = XtCreateManagedWidget("drawing_area",coreWidgetClass,
                  paned,wargs,2);

  /* create a menubuton add it to the menubar*/
  menuButton = XtCreateManagedWidget ("File",menuButtonWidgetClass,
                   menubar,NULL,0);

  /* create a menu add it to the menubutton */
  menu = XtCreatePopupShell("menu",simpleMenuWidgetClass,menuButton,NULL,0);

  entry=XtCreateManagedWidget("Plot",smeBSBObjectClass,menu,NULL,0);
  XtAddCallback(entry,XtNcallback,plot_callback,NULL);
          
  entry=XtCreateManagedWidget("Change Parameters",smeBSBObjectClass,
                menu,NULL,0);
  XtAddCallback(entry,XtNcallback,change_callback,NULL);
          
  entry=XtCreateManagedWidget("Quit",smeBSBObjectClass,menu,NULL,0);
  XtAddCallback(entry,XtNcallback,quit_callback,NULL);
          
  menuButton = XtCreateManagedWidget("Help",menuButtonWidgetClass,
                  menubar,NULL,0);
  
  menu = XtCreatePopupShell("menu",simpleMenuWidgetClass,menuButton,NULL,0);
  
  entry=XtCreateManagedWidget("About",smeBSBObjectClass,menu,NULL,0);
  XtAddCallback(entry,XtNcallback,about_callback,NULL);
 

  /* realize the widgets */ 
  XtRealizeWidget(toplevel);
  
  wm_delete_window = XInternAtom(XtDisplay(toplevel),
                  "WM_DELETE_WINDOW",0);
  XSetWMProtocols(XtDisplay(toplevel),XtWindow(toplevel),
                  &wm_delete_window,1);
  
  values.foreground=BlackPixel(XtDisplay(toplevel),0);  
  gc1=XCreateGC (XtDisplay (toplevel), XtWindow (drawing_area),
                  GCForeground,&values);
 
  mainwin=XtWindow(drawing_area);                        
          
  XtAddEventHandler(drawing_area,ExposureMask ,FALSE,
                  (XtEventHandler)redraw,NULL);
  
  XtAddEventHandler(toplevel,StructureNotifyMask,FALSE,
                  (XtEventHandler)redraw,NULL);
  
  display=XtDisplay(toplevel);



  winAttr.backing_store = Always;
  winAttr.save_under=1;
  XChangeWindowAttributes(display,mainwin,CWBackingStore|CWSaveUnder,&winAttr);
  
  XGetGeometry(display,mainwin,&DefaultRootWindow(display),
                  &x,&y,&width,&height,&dummy1,&dummy2);

}
#endif

