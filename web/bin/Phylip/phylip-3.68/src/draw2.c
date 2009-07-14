
#include <stdio.h> /* Metrowerks for windows defines WIN32 here */
#define swap_m(x,y) temp = y,y=x,x=temp;

extern long winheight;
extern long winwidth;

#ifdef WIN32
#include <windows.h>
HDC hdc;


/******* Menu Defines *******/

#define IDM_ABOUT      1000
#define IDM_PLOT       1001
#define IDM_CHANGE     1002
#define IDM_QUIT       1003

#define XWINPERCENT    0.66
#define YWINPERCENT    0.66
#endif


#ifdef QUICKC
extern struct videoconfig myscreen;
#endif

#ifdef OSX_CARBON
#include <Carbon/Carbon.h>
#endif

#include "draw.h"
#include "phylip.h"



static long eb[]={
  0 , 1 ,2 ,3 ,55,45,46,47,22,5,37,11,12,13,14,15,16,17,18,19,60,61,50,38,
  24, 25,63,39,28,29,30,31,64,90,127,123,91,108,80,125,77,93,92,78,107,96,
  75,97,240,241,242,243,244,245,246,247,248,249,122,94,76,126,110,111, 124,
  193,194,195,196,197,198,199,200,201,209,210,211, 212,213,214,215,216,217,
  226,227,228,229,230,231,232,233,173,224,189, 95,109,121,129,130,131,132,
  133,134,135,136,137,145,146,147,148,149,150, 151, 152,153,162,163,164,165,
  166,167,168,169,192,79,208,161,7};

double oldxreal, oldyreal;
boolean didenter, didexit, curvetrue;
extern long vrmlplotcolor;

extern double oldx, oldy ;
extern node *root;
extern long    nmoves, oldpictint ;
extern long    rootmatrix[51][51];
extern long    strpbottom,strptop,strpwide,strpdeep;
extern boolean dotmatrix, empty, previewing;
extern double  ynow, ysize, xsize, yunitspercm;
extern FILE          *plotfile;
extern plottertype   plotter;
extern striptype     stripe;

extern long treecolor, namecolor, vrmlskycolorfar, vrmlskycolornear,
  vrmlgroundcolorfar, vrmlgroundcolornear;
extern colortype colors[7];
extern vrmllighttype vrmllights[3];


/* Added by Dan F. for the new previewing paradigm */
extern double labelline,linewidth,oldxhigh,oldxlow,oldyhigh,oldylow,
  vrmllinewidth, raylinewidth,treeline,oldxsize,oldysize,oldxunitspercm,
  oldyunitspercm,oldxcorner,oldycorner,clipx0,clipx1,clipy0,clipy1;

/* func. protocol added for vrml - danieyek 981111 */

extern long          strpdiv,hpresolution;
extern boolean       preview,pictbold,pictitalic,
  pictshadow,pictoutline;
extern double        expand,xcorner,xnow,xscale,xunitspercm,
  ycorner,yscale,labelrotation,
  labelheight,ymargin,pagex,pagey,paperx,papery,hpmargin,vpmargin;

extern long          filesize;
extern growth        grows;
extern enum {yes,no} penchange,oldpenchange;
extern plottertype   oldplotter,previewer;
extern char             resopts;
extern winactiontype winaction;

#ifndef OLDC
/* function prototypes */
void   plotdot(long, long);
void   circlepoints(int, int, int, int);
void   drawpen(long, long, long);
void   drawfatline(long, long, long, long, long);
void   idellipse(double, double); 
void   splyne(double,double,double,double,boolean,long,boolean,boolean);
static void   putshort(FILE *, int);

static void   putint(FILE *, int);
void   reverse_bits (byte *, int); 
void   makebox_no_interaction(char *, double *, double *, double *, long);
void   void_func(void);
/* function prototypes */
#endif


void plotdot(long ix, long iy)
{
  /* plot one dot at ix, iy */
  long ix0, iy0, iy1 = 0, iy2 = 0;

  iy0 = strptop - iy;
  if ((unsigned)iy0 > strpdeep || ix <= 0 || ix > strpwide)
    return;
  empty = false;
  ix0 = ix;
  switch (plotter) {

  case citoh:
    iy1 = 1;
    iy2 = iy0;
    break;

  case epson:
    iy1 = 1;
    iy2 = 7 - iy0;
    break;

  case oki:
    iy1 = 1;
    iy2 = 7 - iy0;
    break;

  case toshiba:
    iy1 = iy0 / 6 + 1;
    iy2 = 5 - iy0 % 6;
    break;

  case pcx:
    iy1 = iy0 + 1;
    ix0 = (ix - 1) / 8 + 1;
    iy2 = 7 - ((ix - 1) & 7);
    break;

  case pcl:
    iy1 = iy0 + 1;
    ix0 = (ix - 1) / 8 + 1;
    iy2 = 7 - ((ix - 1) & 7);
    break;

  case bmp: 
    iy1 = iy0 + 1;
    ix0 = (ix - 1) / 8 + 1;
    iy2 = 7 - ((ix - 1) & 7);
  case xbm:
  case gif:
    iy1 = iy0 + 1;
    ix0 = (ix - 1) / 8 + 1;
    iy2 = (ix - 1) & 7;
    break;

  case lw:
  case hp:
  case tek:
  case mac:
  case houston:
  case decregis:
  case fig:
  case pict:
  case ray:
  case pov:
  case idraw:
  case ibm:
  case other:
    break;
  default:        /* cases xpreview and vrml not handled        */
    break;
    /* code for making dot array for a new printer
      goes here */
  }
  stripe[iy1 - 1][ix0 - 1] |= (unsigned char)1<<iy2;
}  /* plotdot */


void circlepoints(int x, int y, int x0, int y0)
{
  /* circlepoints is consecutively passed a circle center and x,y coordinates *
   * for 1 octant of a circle. Since the circle is symmetrical, we can use    *
   * this to plot a complete circle with adjacent pixels (to avoid holes      *
   * often associated with diagonal stairstepping.                            */
  
  plotdot(x0+x,y0+y);
  plotdot(x0+x,y0+y-1);
  plotdot(x0+y,y0+x);
  plotdot(x0+y-1,y0+x);
  
  plotdot(x0-x,y0+y);
  plotdot(x0-x,y0+y-1);
  plotdot(x0-y,y0+x);
  plotdot(x0-y+1,y0+x);
  
  plotdot(x0+x,y0-y);
  plotdot(x0+x,y0-y+1);
  plotdot(x0+y,y0-x);
  plotdot(x0+y-1,y0-x);
  
  plotdot(x0-x,y0-y);
  plotdot(x0-x,y0-y+1);
  plotdot(x0-y,y0-x);
  plotdot(x0-y+1,y0-x);
}  /* circlepoints */


void drawpen(long x0, long y0, long radius)
{
 int x,y,d,deltaE,deltaSE;

  x = 0;
  y = radius;
  d = 1-radius;
  deltaE = 3;
  deltaSE = -2 * radius + 5;
  circlepoints(x,y,x0,y0);

  while (y > x){
    if (d < 0) {
      d = d + deltaE;
      deltaE += 2;
      deltaSE += 2;
      x++;
    }
    else {
      d+=deltaSE;
      deltaE += 2;
      deltaSE += 4;
      x++;
      y--;
    }
    circlepoints(x,y,x0,y0);
  }
} /* drawpen */


void drawfatline(long ixabs, long iyabs, long ixnow, long iynow,
                        long penwide)
{
  long temp, xdiff, ydiff, err, x1, y1;

  didenter = false;
  didexit = false;

  if (ixabs < ixnow) {
    temp = ixnow;
    ixnow = ixabs;
    ixabs = temp;
    temp = iynow;
    iynow = iyabs;
    iyabs = temp;
  }
  xdiff = ixabs - ixnow;
  ydiff = iyabs - iynow;
  if (ydiff >= 0) {
    if (xdiff >= ydiff) {
      err = -(xdiff / 2);
      x1 = ixnow;
      while (x1 <= ixabs && !(didenter && didexit)) {
        drawpen(x1, iynow, penwide);
        err += ydiff;
        if (err > 0) {
          iynow++;
          err -= xdiff;
        }
        x1++;
      }
      return;
    }
    err = -(ydiff / 2);
    y1 = iynow;
    while (y1 < iyabs && !(didenter && didexit)) {
      drawpen(ixnow, y1, penwide);
      err += xdiff;
      if (err > 0) {
        ixnow++;
        err -= ydiff;
      }
      y1++;
    }
    return;
  }
  if (xdiff < -ydiff) {
    err = ydiff / 2;
    y1 = iynow;
    while (y1 >= iyabs && !(didenter && didexit)) {
      drawpen(ixnow, y1, penwide);
      err += xdiff;
      if (err > 0) {
        ixnow++;
        err += ydiff;
      }
      y1--;
    }
    return;
  }
  err = -(xdiff / 2);
  x1 = ixnow;
  while (x1 <= ixabs && !(didenter && didexit)) {
    drawpen(x1, iynow, penwide);
    err -= ydiff;
    if (err > 0) {
      iynow--;
      err -= xdiff;
    }
    x1++;
  }
}  /* drawfatline */


void plot(pensttstype pen, double xabs, double yabs)
{
  long xhigh, yhigh, xlow, ylow, ixnow, iynow, ixabs, iyabs,
       cdx, cdy, temp, i;
  long pictint;
  double newx, newy, dx, dy, lscale, dxreal, dyreal;
  Char picthi, pictlo;

  /* added to give every line a name in vrml! - danieyek 981110 */
  static long lineCount = 0;
  /* Record the first node as the coordinate for viewpoint! */
  static int firstNodeP=1;
  double distance, angle;
  double episilon = 1.0e-10;

  /* For povray, added by Dan F. */
  char texture_string[7];

/* remember to respect & translate for clipping region, clip{x,y}{0,1} */

  if (!dotmatrix || previewing) {
    switch (plotter) {

    case xpreview:
      if (pen == pendown) {
#ifdef X
        XDrawLine(display,mainwin,gc1,(long)oldx,(long)(height-oldy),
                                      (long)xabs,(long)(height-yabs));
#endif
      }
      oldx = xabs;
      oldy = yabs;
      break;

    case winpreview:
#ifdef WIN32
      if (pen == pendown) {
        LineTo(hdc, (int) xabs, (int)(winheight-yabs)); 
      }
      else {
        MoveToEx(hdc, (int) xabs, (int) (winheight-yabs), (LPPOINT) NULL); 
      }
#endif
      break;

    case tek:
      if (pen == penup) {
        if (previewing)
          putchar('\035');
        else
          putc('\035', plotfile);
      }
      ixnow = (long)floor(xabs + 0.5);
      iynow = (long)floor(yabs + 0.5);
      xhigh = ixnow / 32;
      yhigh = iynow / 32;
      xlow = ixnow & 31;
      ylow = iynow & 31;
      if (!ebcdic) {
        if (yhigh != oldyhigh) {
          if (previewing)
            putchar(yhigh + 32);
          else
            putc(yhigh + 32, plotfile);
        }
        if (ylow != oldylow || xhigh != oldxhigh) {
          if (previewing)
            putchar(ylow + 96);
          else
            putc(ylow + 96, plotfile);
        }
        if (xhigh != oldxhigh) {
          if (previewing)
            putchar(xhigh + 32);
          else
            putc(xhigh + 32, plotfile);
        }
        if (previewing)
          putchar(xlow + 64);
        else
          putc(xlow + 64, plotfile);
      } else {  /* DLS/JMH -- for systems that use EBCDIC coding */
        if (yhigh != oldyhigh) {
          if (previewing)
            putchar(eb[yhigh + 32]);
          else
            putc(eb[yhigh + 32], plotfile);
        }
        if (ylow != oldylow || xhigh != oldxhigh) {
          if (previewing)
            putchar(eb[ylow + 96]);
          else
            putc(eb[ylow + 96], plotfile);
        }
        if (xhigh != oldxhigh) {
          if (previewing)
            putchar(eb[xhigh + 32]);
          else
            putc(eb[xhigh + 32], plotfile);
        }
        if (previewing)
          putchar(eb[xlow + 64]);
        else
          putc(eb[xlow + 64], plotfile);
      }

      oldxhigh = xhigh;
      oldxlow = xlow;
      oldyhigh = yhigh;
      oldylow = ylow;
      break;

    case hp:
      if (pen == pendown)
        fprintf(plotfile, "PD");
      else
        fprintf(plotfile, "PU");
      pout((long)floor(xabs + 0.5));
      putc(',', plotfile);
      pout((long)floor(yabs + 0.5));
      fprintf(plotfile, ";\n");
      break;

    case pict:
      newx = floor(xabs + 0.5);
      newy = floor(ysize * yunitspercm - yabs + 0.5);
      if (pen == pendown) {
        if (linewidth > 5) {
          dxreal = xabs - oldxreal;
          dyreal = yabs - oldyreal;
          lscale = sqrt(dxreal * dxreal + dyreal * dyreal) /
            (fabs(dxreal) + fabs(dyreal));
          pictint = (long)(lscale * linewidth + 0.5);

          if (pictint == 0)
            pictint = 1;
          if (pictint != oldpictint) {
            picthi = (Char)(pictint / 256);
            pictlo = (Char)(pictint & 255);
            fprintf(plotfile, "\007%c%c%c%c", picthi, pictlo, picthi, pictlo);
          }
          oldpictint = pictint;
        }
        fprintf(plotfile, " %c%c%c%c",
                (Char)((long) oldy / 256), (Char)((long) oldy & 255), (Char)((long) oldx / 256),
                (Char)((long) oldx & 255));
        fprintf(plotfile, "%c%c%c%c",
                (Char)((long)newy / 256), (Char)((long)newy & 255), (Char)((long)newx / 256),
                (Char)((long)newx & 255));
        }
      oldxreal = xabs;
      oldyreal = yabs;
      oldx = newx;
      oldy = newy;
      break;

    case ray:
      if (pen == pendown) {
        if (linewidth != treeline) {
          if (raylinewidth > labelline) {
            raylinewidth = labelline;
            fprintf(plotfile, "end\n\n");
            fprintf(plotfile, "name species_names\n");
            fprintf(plotfile, "grid 22 22 22\n");
          }
        }

        if (oldxreal != xabs || oldyreal != yabs) {
          raylinewidth *= 0.99999;
          fprintf(plotfile, "cylinder %8.7f %6.3f 0 %6.3f %6.3f 0 %6.3f\n",
                  raylinewidth, oldxreal, oldyreal, xabs, yabs);
          fprintf(plotfile, "sphere %8.7f %6.3f 0 %6.3f\n",
                  raylinewidth, xabs, yabs);
        }
      }
      oldxreal = xabs;
      oldyreal = yabs;
      break;

    case pov:
      /* Default to writing out tree texture... */
      strcpy (texture_string, TREE_TEXTURE);

      if (pen == pendown) {
        if (linewidth != treeline) {
          /* Change the texture to name texture */
          strcpy (texture_string, NAME_TEXTURE);

          if (raylinewidth > labelline) {
            raylinewidth = labelline;
            fprintf(plotfile, "\n// Now, the species names:\n\n");
          }
        }

        if (oldxreal != xabs || oldyreal != yabs) {
          raylinewidth *= 0.99999;
          fprintf(plotfile,
                  "cylinder { <%6.3f, 0, %6.3f,>, <%6.3f, 0, %6.3f>, %8.7f \n",
                  oldxreal, oldyreal, xabs, yabs, raylinewidth);
          fprintf(plotfile, "\ttexture { %s } }\n", texture_string);
          fprintf(plotfile, "sphere { <%6.3f, 0, %6.3f>, %8.7f \n",
                  xabs, yabs, raylinewidth);
          fprintf(plotfile, "\ttexture { %s } }\n", texture_string);

        }
      }
      oldxreal = xabs;
      oldyreal = yabs;
      break;

    case lw:
      if (pen == pendown){
        /* If there's NO possibility that the line interesects the page, 
         * leave it out. Otherwise, let postscript clip it to the page. */
          if (!((xabs > clipx1*xunitspercm &&  oldx > clipx1*xunitspercm) ||
              (xabs < clipx0*xunitspercm &&  oldx < clipx0*xunitspercm) ||
              (yabs > clipy1*yunitspercm &&  oldy > clipy1*yunitspercm) ||
              (yabs < clipy0*yunitspercm &&  oldy < clipy0*yunitspercm)))
            fprintf(plotfile, "%8.2f %8.2f %8.2f %8.2f l\n",
                    oldx-(clipx0*xunitspercm), oldy-(clipy0*yunitspercm),
                    xabs-(clipx0*xunitspercm), yabs-(clipy0*yunitspercm));
        }
        oldx     = xabs,
        oldy     = yabs;
      break;

    case idraw:
      if (pen == pendown) {
        fprintf(plotfile, "Begin %%I Line\n");
        fprintf(plotfile, "%%I b 65535\n");
        fprintf(plotfile, "%d 0 0 [] 0 SetB\n",
                ((linewidth>=1.0) ? (int)linewidth : 1));
        fprintf(plotfile, "%%I cfg Black\n");
        fprintf(plotfile, "0 0 0 SetCFg\n");
        fprintf(plotfile, "%%I cbg White\n");
        fprintf(plotfile, "1 1 1 SetCBg\n");
        fprintf(plotfile, "%%I p\n");
        fprintf(plotfile, "0 SetP\n");
        fprintf(plotfile, "%%I t\n");
        fprintf(plotfile, "[ 0.01 0 0 0.01 216 285 ] concat\n");
        fprintf(plotfile, "%%I\n");
        fprintf(plotfile, "%ld %ld %ld %ld Line\n",
                (long)(100.0 * (oldxreal+0.5)),
                (long)(100.0 * (oldyreal+0.5)),
                (long)(100.0 * (xabs+0.5)),
                (long)(100.0 * (yabs+0.5)));
        fprintf(plotfile, "End\n\n");

        if (linewidth >= 4.0) {
            fprintf(plotfile, "Begin %%I Elli\n");
            fprintf(plotfile, "%%I b 65535\n");
            fprintf(plotfile, "1 0 0 [] 0 SetB\n");
            fprintf(plotfile, "%%I cfg Black\n");
            fprintf(plotfile, "0 0 0 SetCFg\n");
            fprintf(plotfile, "%%I cbg White\n");
            fprintf(plotfile, "1 1 1 SetCBg\n");
            fprintf(plotfile, "%%I p\n");
            fprintf(plotfile, "0 SetP\n");
            fprintf(plotfile, "%%I t\n");
            fprintf(plotfile, "[ 0.01 0 0 0.01 216 285 ] concat\n");
            fprintf(plotfile, "%%I\n");
            fprintf(plotfile, "%ld %ld %ld %ld Elli\n",
                    (long)(100.0 * (oldxreal+0.5)),
                    (long)(100.0 * (oldyreal+0.5)),
                    (long)(100.0 * (linewidth/2)) - 100,
                    (long)(100.0 * (linewidth/2)) - 100);
            fprintf(plotfile, "End\n");
            
            fprintf(plotfile, "Begin %%I Elli\n");
            fprintf(plotfile, "%%I b 65535\n");
            fprintf(plotfile, "1 0 0 [] 0 SetB\n");
            fprintf(plotfile, "%%I cfg Black\n");
            fprintf(plotfile, "0 0 0 SetCFg\n");
            fprintf(plotfile, "%%I cbg White\n");
            fprintf(plotfile, "1 1 1 SetCBg\n");
            fprintf(plotfile, "%%I p\n");
            fprintf(plotfile, "0 SetP\n");
            fprintf(plotfile, "%%I t\n");
            fprintf(plotfile, "[ 0.01 0 0 0.01 216 285 ] concat\n");
            fprintf(plotfile, "%%I\n");
            fprintf(plotfile, "%ld %ld %ld %ld Elli\n",
                    (long)(100.0 * (xabs+0.5)),
                    (long)(100.0 * (yabs+0.5)),
                    (long)(100.0 * (linewidth/2)) - 100,
                    (long)(100.0 * (linewidth/2)) - 100);
            fprintf(plotfile, "End\n");
        }
      }
      oldxreal = xabs;
      oldyreal = yabs;
      break;

    case ibm:
#ifdef TURBOC
    newx = floor(xabs + 0.5);
    newy = fabs(floor(yabs) - getmaxy());
    if (pen == pendown)
      line((long)oldx,(long)oldy,(long)newx,(long)newy);
    oldx = newx;
    oldy = newy;
#endif
#ifdef QUICKC
    newx = floor(xabs + 0.5);
    newy = fabs(floor(yabs) - myscreen.numypixels);

    if (pen == pendown)
        _lineto((long)newx,(long)newy);
    else
        _moveto((long)newx,(long)newy);
    oldx = newx;
    oldy = newy;

#endif
    break;
     case mac:
#ifdef MAC
      if (pen == pendown){
        LineTo((int)floor((double)xabs + 0.5),
               winheight - (long)floor((double)yabs + 0.5)+MAC_OFFSET);}
      else{
        MoveTo((int)floor((double)xabs + 0.5),
               winheight - (long)floor((double)yabs + 0.5)+MAC_OFFSET);}
#endif

      break;

    case houston:
      if (pen == pendown)
        fprintf(plotfile, "D ");
      else
        fprintf(plotfile, "U ");
      pout((long)((long)floor(xabs + 0.5)));
      putc(',', plotfile);
      pout((long)((long)floor(yabs + 0.5)));
      putc('\n', plotfile);
      break;

    case decregis:
      newx = floor(xabs + 0.5);
      newy = fabs(floor(yabs + 0.5) - 479);
      if (pen == pendown) {
        if (previewing) {
          printf("P[");
          pout((long)oldx);
          putchar(',');
          pout((long)oldy);
          printf("]V[");
          pout((long)newx);
          putchar(',');
          pout((long)newy);
          putchar(']');
        } else {
          fprintf(plotfile, "P[");
          pout((long)oldx);
          putc(',', plotfile);
          pout((long)oldy);
          fprintf(plotfile, "]V[");
          pout((long)newx);
          putc(',', plotfile);
          pout((long)newy);
          putc(']', plotfile);
        }
        nmoves++;
        if (nmoves == 3) {
          nmoves = 0;
          if (previewing)
            putchar('\n');
          else
            putc('\n', plotfile);
        }
      }
      oldx = newx;
      oldy = newy;
      break;

    case fig:
      newx = floor(xabs + 0.5);
      newy = floor(yabs + 0.5);
      if (pen == pendown) {
        fprintf(plotfile, "2 1 0 %5ld 0 0 0 0 0.000 0 0\n",
                (long)floor(linewidth + 0.5) + 1);
        fprintf(plotfile, "%5ld%5ld%5ld%5ld 9999 9999\n",
                (long)oldx, 606 - (long) oldy, (long)newx, 606 - (long)newy);
        fprintf(plotfile,
          "1 3 0  1 0 0 0 21 0.00 1 0.0 %5ld%5ld%5ld %5ld %5ld%5ld%5ld 349\n",
                (long)oldx, 
                606 - (long) oldy, 
                (long)floor(linewidth / 2 + 0.5),
                (long)floor(linewidth / 2 + 0.5), 
                (long)oldx, 606 - (long)oldy, 
                606 - (long)oldy);
        fprintf(plotfile,
          "1 3 0  1 0 0 0 21 0.00 1 0.0 %5ld%5ld%5ld %5ld %5ld%5ld%5ld 349\n",
                (long)newx, 
                606 - (long)newy, 
                (long)floor(linewidth / 2 + 0.5),
                (long)floor(linewidth / 2 + 0.5), 
                (long)newx, 
                606 - (long)newy, 
                606 - (long)newy);
      }
      oldx = newx;
      oldy = newy;
      break;

    case vrml:

      newx = xabs;
      newy = yabs;
      /* if this is the root node, 
         use the coordinates to define the view point */
      if (firstNodeP-- == 1)
      {
        fprintf(plotfile, "#VRML V2.0 utf8\n");
        fprintf(plotfile, "    NavigationInfo {\n");
        fprintf(plotfile, "      headlight FALSE\n");
        fprintf(plotfile, "    }\n");
        fprintf(plotfile, "    Viewpoint\n");
        fprintf(plotfile, "    {\n");
        fprintf(plotfile, "      position %f %f %f\n", xsize/2, ysize/2, ysize*1.2);
        fprintf(plotfile, "      description \"Entry View\"\n");
        fprintf(plotfile, "    }\n");

        for (i=0; i<3; i++) {
          fprintf(plotfile, "    PointLight {\n");
          fprintf(plotfile, "      on TRUE\n");
          fprintf(plotfile, "      intensity %f\n",
               vrmllights[i].intensity);
          fprintf(plotfile, "      ambientIntensity 0.0\n");
          fprintf(plotfile, "      color 1.0 1.0 1.0\n");
          fprintf(plotfile, "      location %f %f %f\n",
               vrmllights[i].x,
               vrmllights[i].y,
               vrmllights[i].z);
          fprintf(plotfile, "      attenuation 0.0 0.0 0.0\n");
          fprintf(plotfile, "      radius 200.0\n");
          fprintf(plotfile, "    }\n");
        }

        fprintf(plotfile, "    Background\n");
        fprintf(plotfile, "    {\n");
        fprintf(plotfile, "      skyAngle [1.75]\n");
        fprintf(plotfile, "      skyColor [%f %f %f, %f %f %f]\n",
                  colors[vrmlskycolornear-1].red, 
                  colors[vrmlskycolornear-1].green, 
                  colors[vrmlskycolornear-1].blue,
                  colors[vrmlskycolorfar-1].red, 
                  colors[vrmlskycolorfar-1].green, 
                  colors[vrmlskycolorfar-1].blue);
        fprintf(plotfile, "      groundAngle[0 1.57 3.14]\n");
        fprintf(plotfile,
            "      groundColor [0.9 0.9 0.9, 0.7 0.7 0.7, %f %f %f]\n",
                  colors[vrmlgroundcolorfar-1].red,
                  colors[vrmlgroundcolorfar-1].green, 
                  colors[vrmlgroundcolorfar-1].blue);
        fprintf(plotfile, "    }\n");
      }

      if (pen == penup)
      {/* pen down = beginning of a new path */
      }
      else if (pen == pendown)
      {/* pen up = continue, line may not end yet. */

        if (linewidth != treeline) {
          if (vrmllinewidth > labelline) {
            vrmllinewidth = labelline;
            vrmlplotcolor = namecolor;
          }
        }

        distance = sqrt((newy - oldy)*(newy - oldy) + (newx - oldx)*(newx - oldx));
        angle = computeAngle(oldx, oldy, newx, newy);

        if (distance >= episilon)
        {
          fprintf(plotfile, "    DEF Line%ld Transform\n", lineCount++);
          fprintf(plotfile, "    {\n");
          fprintf(plotfile, "      rotation 0 0 1 %f\n", angle);
          fprintf(plotfile, "      translation %f %f 0\n", oldx, oldy);
          fprintf(plotfile, "      children\n");
          fprintf(plotfile, "      [\n");
          fprintf(plotfile, "        Shape\n");
          fprintf(plotfile, "        {\n");
          fprintf(plotfile, "          appearance Appearance\n");
          fprintf(plotfile, "          {\n");
          fprintf(plotfile, "            material Material { diffuseColor %f %f %f}\n",
                  colors[vrmlplotcolor-1].red, 
                  colors[vrmlplotcolor-1].green, 
                  colors[vrmlplotcolor-1].blue);
          fprintf(plotfile, "          }\n");
          fprintf(plotfile, "          geometry Sphere\n");
          fprintf(plotfile, "          {\n");
/*          vrmllinewidth *= 0.99999; */
          fprintf(plotfile, "            radius %f\n", vrmllinewidth);
          fprintf(plotfile, "          }\n");
          fprintf(plotfile, "        }\n");
          fprintf(plotfile, "        Transform\n");
          fprintf(plotfile, "        {\n");
          fprintf(plotfile, "          rotation 0 0 1 -1.570796327\n" );
          fprintf(plotfile, "          translation %f 0 0\n", distance/2);
          fprintf(plotfile, "          children\n");
          fprintf(plotfile, "          [\n");
          fprintf(plotfile, "            Shape\n");
          fprintf(plotfile, "            {\n");
          fprintf(plotfile, "              appearance Appearance\n");
          fprintf(plotfile, "              {\n");
          fprintf(plotfile, "                material Material { diffuseColor %f %f %f}\n",
                  colors[vrmlplotcolor-1].red, 
                  colors[vrmlplotcolor-1].green, 
                  colors[vrmlplotcolor-1].blue );
          fprintf(plotfile, "              }\n");
          fprintf(plotfile, "              geometry Cylinder\n");
          fprintf(plotfile, "              {\n");
          /* line radius affects end sphere's size */
/*          vrmllinewidth *= 0.99999; */
          fprintf(plotfile, "                radius %f\n", vrmllinewidth);
          fprintf(plotfile, "                height %f\n", distance);
          fprintf(plotfile, "              }\n");
          fprintf(plotfile, "            }\n");
          fprintf(plotfile, "          ]\n");
          fprintf(plotfile, "        }\n");
          fprintf(plotfile, "        Transform\n");
          fprintf(plotfile, "        {\n");
          fprintf(plotfile, "          translation %f 0 0\n", distance);
          fprintf(plotfile, "          children\n");
          fprintf(plotfile, "          [\n");
          fprintf(plotfile, "            Shape\n");
          fprintf(plotfile, "            {\n");
          fprintf(plotfile, "              appearance Appearance\n");
          fprintf(plotfile, "              {\n");
          fprintf(plotfile, "                material Material { diffuseColor %f %f %f}\n",
                  colors[vrmlplotcolor-1].red, 
                  colors[vrmlplotcolor-1].green, 
                  colors[vrmlplotcolor-1].blue );
          fprintf(plotfile, "              }\n");
          fprintf(plotfile, "              geometry Sphere\n");
          fprintf(plotfile, "              {\n");
          /* radius affects line size */
/*          vrmllinewidth *= 0.99999; */
          fprintf(plotfile, "                radius %f\n", vrmllinewidth);
          fprintf(plotfile, "              }\n");
          fprintf(plotfile, "            }\n");
          fprintf(plotfile, "          ]\n");
          fprintf(plotfile, "        }\n");
          fprintf(plotfile, "      ]\n");
          fprintf(plotfile, "    }\n");
        }
      }
      else
      {
        fprintf(stderr, "ERROR: Programming error in plot().");
      }

      oldx = newx;
      oldy = newy;
      break;

    case epson:
    case oki:
    case citoh:
    case toshiba:
    case pcx:
    case pcl:
    case bmp:
    case xbm:
    case gif:
    case other:
      break;
      /* code for a pen move on a new plotter goes here */
    }
    return;
  }
  if (pen == pendown) {
    ixabs = (long)floor(xabs + 0.5);
    iyabs = (long)floor(yabs + 0.5);
    ixnow = (long)floor(xnow + 0.5);
    iynow = (long)floor(ynow + 0.5);
    if (ixnow > ixabs) {
      temp = ixnow;
      ixnow = ixabs;
      ixabs = temp;
      temp = iynow;
      iynow = iyabs;
      iyabs = temp;
    }
    dx = ixabs - ixnow;
    dy = iyabs - iynow;
   /* if (dx + fabs(dy) <= 0.0)
      c = 0.0;
    else
      c = 0.5 * linewidth / sqrt(dx * dx + dy * dy); */
    cdx = (long)floor(linewidth + 0.5);
    cdy = (long)floor(linewidth + 0.5);
    if ((iyabs + cdx >= strpbottom || iynow + cdx >= strpbottom) &&
        (iyabs - cdx <= strptop || iynow - cdx <= strptop)) {
      drawfatline(ixnow,iynow,ixabs,iyabs,(long)floor(linewidth+0.5));
    }
  }

  xnow = xabs;
  ynow = yabs;

  /* Bitmap Code to plot (xnow,ynow) to (xabs,yabs)                 */
} /* plot                                                           */


void idellipse(double x, double y) 
{
  fprintf(plotfile, "Begin %%I Elli\n");
  fprintf(plotfile, "%%I b 65535\n");
  fprintf(plotfile, "1 0 0 [] 0 SetB\n");
  fprintf(plotfile, "%%I cfg Black\n");
  fprintf(plotfile, "0 0 0 SetCFg\n");
  fprintf(plotfile, "%%I cbg White\n");
  fprintf(plotfile, "1 1 1 SetCBg\n");
  fprintf(plotfile, "%%I p\n");
  fprintf(plotfile, "0 SetP\n");
  fprintf(plotfile, "%%I t\n");
  fprintf(plotfile, "[ 0.01 0 0 0.01 216 285 ] concat\n");
  fprintf(plotfile, "%%I\n");
  fprintf(plotfile, "%ld %ld %ld %ld Elli\n",
          (long)(100.0 * (x+0.5)),(long)(100.0 * (y+0.5)),
          (long)(100.0 * (linewidth/2)) - 100,
          (long)(100.0 * (linewidth/2)) - 100);
  fprintf(plotfile, "End\n");
}  /* idellipse */


void splyne(double x1, double y1, double x2, double y2, boolean sense,
                        long segs, boolean head, boolean tail)
{
/* sense is true if line departing from x1,y1 is tangential to x,
   false if tangential to y */

   
  long i,fromx,fromy,tox,toy;
  double f, g, h, x3, y3;
  long ptop, pleft, pbottom, pright, startangle, arcangle;
  double dtheta;
  double sintheta,costheta,sindtheta,cosdtheta,newsintheta,newcostheta;
  double rx,ry; /* axes of ellipse   */
  double ox,oy; /* center of ellipse */
  double prevx,prevy;
  long pictint;
  
  x1 = x1 - (clipx0 * xunitspercm);
  x2 = x2 - (clipx0 * xunitspercm);
  y1 = y1 - (clipy0 * yunitspercm);
  y2 = y2 - (clipy0 * yunitspercm); /* adjust by clipping region */

  switch (plotter) {

  case lw:
    fprintf(plotfile,"stroke %8.2f %8.2f moveto\n",x1,y1);
    if (sense)
      fprintf(plotfile,"%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f curveto\n",
              (x1+(0.55*(x2-x1))), y1, x2, (y1+(0.45*(y2-y1))),
              x2, y2);
    else
      fprintf(plotfile,"%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f curveto\n",
              x1, (y1+(0.55*(y2-y1))), (x1+(0.45*(x2-x1))), y2,
              x2, y2);
    break;

  case pict:
    {
      double dtop, dleft, dbottom, dright,temp;
      if (x1 == x2 || y1 == y2) {
        plot(penup, x1, y1);
        plot(pendown, x2, y2);
      } else {
    if (x2 > x1 && y2 < y1){ swap_m(x2,x1); swap_m(y2,y1); sense = !sense; } 
        
        y1 = (ysize * yunitspercm) - y1;
        y2 = (ysize * yunitspercm) - y2;  

        if (sense) {
          if (x2 > x1) {
            dtop = y2 - y1 + y2;
            dleft = x1 - x2 + x1;
            dbottom = y1;
            dright = x2;
            startangle = 90;
          } else {

            dtop = y2 - y1 + y2;
            dleft = x2;
            dbottom = y1;
            dright = x1 + (x1 - x2);
            startangle = 180;
          }
        }
         else {
          if (x2 > x1) {
            dtop = y1 + (y1 - y2);
            dleft = x1;
            dbottom = y2;
            dright = x2 + (x2 - x1);
            startangle = 270;
          } else {
            dtop = y2;
            dleft = x1;
            dbottom = y1 + y1 - y2;;
                dright = x2 + (x2 - x1);
            startangle = 0;
          }
        }
        arcangle = 90;
         if (dbottom < dtop) {swap_m(dbottom,dtop);}
    if (dleft> dright) {swap_m(dleft,dright);}
 
        ptop    = (long)floor((dtop - 0) + 0.5);
        pleft   = (long)floor(dleft  + 0.5);
        pbottom = (long)floor(dbottom + 0.5) + (long)floor(linewidth + 0.5);
        pright  = (long)floor(dright  + 0.5) + (long)floor(linewidth + 0.5);

    if (!sense)
        pbottom++;
     else
        if (x2 < x1)
            pright++;
         else
            pleft--;
        pictint = 1;    

        fprintf(plotfile,"\140%c%c%c%c%c%c%c%c%c%c%c%c",
                (Char)(ptop / 256), (Char)(ptop % 256),
                (Char)(pleft / 256), (Char)(pleft % 256),
                (Char)(pbottom / 256), (Char)(pbottom % 256),
                (Char)(pright / 256), (Char)(pright % 256),
                (Char)(startangle / 256), (Char)(startangle % 256),
                (Char)(arcangle / 256), (Char)(arcangle % 256));
      }
    }
    break;

  case fig:
   fromx = (long)floor(x1 + 0.5);
   fromy = (long)floor(y1 + 0.5);
   tox = (long)floor(x2 + 0.5);
   toy = (long)floor(y2 + 0.5);
    
    fprintf(plotfile, "3 0 0 %5ld 0 0 0 0 0.000 0 0\n",
            (long)floor(linewidth + 0.5) + 1);
    if (sense)
      fprintf(plotfile, "%5ld%5ld%5ld%5ld%5ld%5ld%5ld%5ld 9999 9999\n",
              fromx, 606 - fromy,
              (long)floor((x1+(0.55*(x2-x1))) + 0.5), 606 - fromy,
              tox, 606 - (long)floor((y1+(0.45*(y2-y1))) + 0.5),
              tox, 606 - toy);
    else
      fprintf(plotfile, "%5ld%5ld%5ld%5ld%5ld%5ld%5ld%5ld 9999 9999\n",
              fromx, 606 - fromy,
              fromx, 606 - (long)floor((y1+(0.55*(y2-y1))) + 0.5),
              (long)floor((x1+(0.45*(x2-x1))) + 0.5), 606 - toy,
              tox, 606 - toy);
    fprintf(plotfile, "1 3 0  1 0 0 0 21 0.00 1 0.0 ");
    fprintf(plotfile, "%5ld%5ld%5ld %5ld %5ld%5ld%5ld 349\n",
            fromx, 606 - fromy, (long)floor(linewidth / 2 + 0.5),
            (long)floor(linewidth / 2 + 0.5), fromx,
            606 - fromy, 606 - fromy);
    fprintf(plotfile, "1 3 0  1 0 0 0 21 0.00 1 0.0 ");
    fprintf(plotfile, "%5ld%5ld%5ld %5ld %5ld%5ld%5ld 349\n",
            tox, 606 - toy, (long)floor(linewidth / 2 + 0.5),
            (long)floor(linewidth / 2 + 0.5), tox,
            606 - toy, 606 - toy);
    break;

  case idraw:
    
    if (head){
      fprintf(plotfile,"Begin %%I Pict\n%%I b u\n%%I cfg u\n%%I cbg u\n");
      fprintf(plotfile,"%%I f u\n%%I p u \n%%I t u\n\n");
      idellipse(x1,y1);
      fprintf(plotfile, "Begin %%I BSpl\n");
      fprintf(plotfile, "%%I b 65535\n");
      fprintf(plotfile, "%ld 0 0 [] 0 SetB\n",
              ((linewidth>=1.0) ? (long)linewidth : 1));
      fprintf(plotfile, "%%I cfg Black\n");
      fprintf(plotfile, "0 0 0 SetCFg\n");
      fprintf(plotfile, "%%I cbg White\n");
      fprintf(plotfile, "1 1 1 SetCBg\n");
      fprintf(plotfile, "none SetP %%I p n\n");
      fprintf(plotfile, "%%I t\n");
      fprintf(plotfile, "[ 0.01 0 0 0.01 216 285 ] concat\n");
      if (tail)
        fprintf(plotfile,"%%I %ld\n",segs+1);
      else
        fprintf(plotfile,"%%I %ld\n",(segs*2)+1);
      fprintf(plotfile, "%ld %ld\n", (long)(100.0 * (x1+0.5)), 
              (long)(100.0 * (y1+0.5))); 
    }
    rx = (fabs(x2 - x1));
    ry = (fabs(y2 - y1));

    if (!sense){
      if (x2 < x1)
        sintheta  = 0.0,
        costheta  = 1.0,
        dtheta = 90.0 / ((double)segs),
        ox = x2,
        oy = y1;
      else
        sintheta  = 0.0,
        costheta  = -1.0,
        dtheta = -90.0 / ((double)segs),
        ox = x2,
        oy = y1;
    }
    else{
      if (x2 < x1)
        sintheta  = -1.0,
        costheta  = 0.0,
        dtheta = -90.0 / ((double)segs),
        ox = x1,
        oy = y2;
      else
        sintheta  = -1.0,
        costheta  = 0.0,
        dtheta = 90.0 / ((double)segs),
        ox = x1,
        oy = y2;
        }
    x3        = x1;
    y3        = y1;
    sindtheta = sin(dtheta * (3.1415926535897932384626433 / 180.0));
    cosdtheta = cos(dtheta * (3.1415926535897932384626433 / 180.0));
    
    for (i = 1; i <= segs; i++) {
      prevx = x3;
      prevy = y3;
      newsintheta = (sintheta * cosdtheta) + (costheta * sindtheta);
      newcostheta = (costheta * cosdtheta) - (sintheta * sindtheta);
      sintheta = newsintheta;
      costheta = newcostheta;
      x3 = ox + (costheta * rx);
      y3 = oy + (sintheta * ry);

      /* adjust spline for better aesthetics: */
      if (i == 1){
        if (sense)
          y3 = (y3 + prevy)  / 2.0;
        else
          x3 = (x3 + prevx) / 2.0;}
      else if (i == segs - 1){
        if (sense)
          x3 = (x3 + x2) / 2.0;
        else
          y3 = (y2 + y3) / 2.0;
      }
      fprintf(plotfile, "%ld %ld\n", (long)(100.0 * (x3+0.5)), 
              (long)(100.0 * (y3+0.5))); 
  }
    if (head && tail) 
      fprintf(plotfile," BSpl\nEnd\n\n"); /* changed for gcc */
      /*fprintf(plotfile,"%ld BSpl\nEnd\n\n"); This is the original */
    else if (tail) 
             fprintf(plotfile," BSpl \nEnd\n\n");  /* changed for gcc */
        /*fprintf(plotfile,"%ld BSpl\nEnd\n\n"); This is the original */
    if (tail)
      idellipse(x2,y2),
      fprintf(plotfile,"\nEnd %%I eop\n\n");
    break;

  case hp:
    plot(penup,x1,y1);    
    if (sense){
      if (x2 > x1)
        fprintf(plotfile,"PD;AA%ld,%ld,90,1;\n",(long)x1,(long)y2);
      else
        fprintf(plotfile,"PD;AA%ld,%ld,-90,1;\n",(long)x1,(long)y2);
    }
    else {
      if (x2 > x1)
        fprintf(plotfile,"PD;AA%ld,%ld,-90,1;\n",(long)x2,(long)y1);
      else
        fprintf(plotfile,"PD;AA%ld,%ld,90,1;\n",(long)x2,(long)y1);
    }
    plot(penup,x2,y2); fprintf(plotfile,"PD;PU;");

/*    else
      fprintf(plotfile,"PD;AA%ld,%ld,90,1;\n",(long)x2,(int)y1); */
    plot(penup,x2,y2);
    break;
  default:
    for (i = 1; i <= 2*segs; i++) {
      f = (double)i / (2*segs);
      g = (double)i / (2*segs);
      h = 1.0 - sqrt(1.0 - g * g);
      if (sense) {
        x3 = x1 * (1.0 - f) + x2 * f;
        y3 = y1 + (y2 - y1) * h;
      } else {
        x3 = x1 + (x2 - x1) * h;
        y3 = y1 * (1.0 - f) + y2 * f;
      }
      plot(pendown, x3, y3);
    }
    break;
  }
}  /* splyne */


void swoopspline(double x1, double y1, double x2, double y2, double x3,
                        double y3, boolean sense, long segs)
{
  splyne(x1,y1,x2,y2,sense,segs/4,true,false); 
  splyne(x2,y2,x3,y3,(boolean)(!sense),segs/4,false,true); 
}  /* swoopspline */


void curvespline(double x1, double y1, double x2, double y2,
                        boolean sense, long segs)
{
  splyne(x1,y1,x2,y2,sense,segs/2,true,true); 
}  /* curvespline */


/*******************************************/
static void putshort(FILE *fp, int i)
{
  int c, c1;

  c = ((unsigned int ) i) & 0xff;  c1 = (((unsigned int) i)>>8) & 0xff;
  putc(c, fp);   putc(c1,fp);
}  /* putshort */
/*******************************************/


static void putint(FILE *fp, int i)
{
  int c, c1, c2, c3;
  c  = ((unsigned int ) i)      & 0xff;  
  c1 = (((unsigned int) i)>>8)  & 0xff;
  c2 = (((unsigned int) i)>>16) & 0xff;
  c3 = (((unsigned int) i)>>24) & 0xff;

  putc(c, fp);   putc(c1,fp);  putc(c2,fp);  putc(c3,fp);
}  /* ptint */


void write_bmp_header (FILE *plotfile,int width,int height)
{
  /*
   *  write a 1-bit image header
   *
   */

  byte r1[2],g1[2],b1[2] ;

  int i, bperlin;

  r1[0] = (long) 255;   /* Black */
  g1[0] = (long) 255;
  b1[0] = (long) 255;

  r1[1] = 0;
  g1[1] = 0;
  b1[1] = 0;

  bperlin = ((width + 31) / 32) * 4;   /* # bytes written per line */

  putc('B', plotfile);  
  putc('M', plotfile);        /* BMP file magic number */

  /* compute filesize and write it */
  i = 14 +                    /* size of bitmap file header */
      40 +                    /* size of bitmap info header */
      8 +                     /* size of colormap */
      bperlin * height;       /* size of image data */

  putint(plotfile, i);
  putshort(plotfile, 0);          /* reserved1 */
  putshort(plotfile, 0);          /* reserved2 */
  putint(plotfile, 
         14 + 40 + 8);            /* offset from BOfile to BObitmap */
  putint(plotfile, 40);           /* biSize: size of bitmap info header */
  putint(plotfile, width);        /* Width */
  putint(plotfile, height);       /* Height */
  putshort(plotfile, 1);          /* Planes:  must be '1' */
  putshort(plotfile, 1);          /* BitCount: 1 */
  putint(plotfile, 0);            /* Compression:  BI_RGB = 0 */
  putint(plotfile, bperlin*height);/* SizeImage:  size of raw image data */
  putint(plotfile, 75 * 39);      /* XPelsPerMeter: (75dpi * 39 in. per meter) */
  putint(plotfile, 75 * 39);      /* YPelsPerMeter: (75dpi * 39 in. per meter) */
  putint(plotfile, 2);            /* ClrUsed: # of colors used in cmap */
  putint(plotfile, 2);            /* ClrImportant: same as above */

  /* write out the colormap */
  for (i = 0 ; i < 2 ; i++) {
    putc(b1[i],plotfile);
    putc(g1[i],plotfile);
    putc(r1[i],plotfile);
    putc(0,    plotfile);
  }
}  /* write_bmp_header */


void reverse_bits (byte *full_pic, int location) 
{
  /* Reverse all the bits at location */
  int i, loop_end ; 
  byte orig, reversed;

  /* initialize...*/
  orig = full_pic[location] ;
  reversed = (byte) '\0'; 
  loop_end = sizeof (byte) * 8 ;

  if (orig == (byte) '\0') {
    /* No need to do anything for 0 bytes, */
    return ;
  } else {
    for (i = 0 ; i < loop_end ; i++) {
      reversed = (reversed << 1) | (orig & 1) ;
      orig   >>= 1 ;
    }
    full_pic[location] = reversed ;
  }
}  /* reverse_bits */


void turn_rows (byte *full_pic, int padded_width, int height)
{
  int i, j;
  int midpoint = padded_width / 2 ;
  byte temp ; /* For the swap call */

  for (j = 0 ; j < height ; j++) {
    for (i = 0 ; i < midpoint ; i++) {

      reverse_bits (full_pic, (j * padded_width) + i);
      reverse_bits (full_pic, (j * padded_width) + (padded_width - i));
      swap_m (full_pic[(j * padded_width) + i],
            full_pic[(j * padded_width) + (padded_width - i)]) ;
    }
    /* Then do the midpoint */
    reverse_bits (full_pic, (j * padded_width) + midpoint);
  }
}  /* turn_rows */


void translate_stripe_to_bmp(striptype *stripe, byte *full_pic,
                        int increment, int width, int div, int *total_bytes) 
{
  int padded_width, i, j, offset, pad_size,
    total_stripes, last_stripe_offset, truncated_stripe_height ;
  if (div == 0)
    /* For some reason this is called once without valid data */
    return ;
  else if (div == DEFAULT_STRIPE_HEIGHT) {
    /* For a non-last-stripe, figure out if the last stripe is going
       to be shorter than the others, to know how far from the bottom
       things should be offset. */

    truncated_stripe_height = (int) ysize % DEFAULT_STRIPE_HEIGHT;
    
    if (truncated_stripe_height != 0)
      /* The last stripe isn't default height */
      last_stripe_offset = DEFAULT_STRIPE_HEIGHT - ((int) ysize %
                                                    DEFAULT_STRIPE_HEIGHT) ;
    else
      /* Stripes are all default height */
      last_stripe_offset = 0 ;

  } else {
    /* For the last stripe, */
    last_stripe_offset = 0 ; 
  }

  total_stripes        = (int) ceil (ysize / (double) DEFAULT_STRIPE_HEIGHT);

  /* width, padded to be a multiple of 32 bits, or 4 bytes */
  padded_width = ((width + 3)/4) * 4;  
  pad_size     = padded_width - width;

  /* Include pad_size here, as it'll be turned horizontally later */
  offset       = ((total_stripes - increment) *
                  (padded_width * DEFAULT_STRIPE_HEIGHT))
    - (padded_width * last_stripe_offset)
    + pad_size ;

  for (j = div; j >= 0; j--) {
    for (i = 0; i < width; i++) {
      full_pic[offset +        
              (((div-j) * padded_width) 
               + (width-i))] = (byte) (*stripe)[j][i];
      (*total_bytes)++ ;
    }

    /* Take into account the padding */
    (*total_bytes) += pad_size ;
  }
}  /* translate_stripe_to_bmp */


void write_full_pic(byte *full_pic, int total_bytes)
{
  int i ;
  for (i = 0; i < total_bytes; i++) {
    putc (full_pic[i], plotfile);
  }
}  /* write_full_pic */


void makebox_no_interaction(char *fn, double *xo, double *yo,
                        double *scale, long ntips)
/* fn--fontname        xo,yo--x and y offsets */
{
  /* draw the box on screen which represents plotting area.        */

  long xpag,ypag,i,j;

  oldpenchange   = penchange;
  oldxsize       = xsize;
  oldysize       = ysize;
  oldxunitspercm = xunitspercm;
  oldyunitspercm = yunitspercm;
  oldxcorner     = xcorner;
  oldycorner     = ycorner;
  oldplotter     = plotter;

  plotrparms(ntips);
  xcorner += 0.05 * xsize;
  ycorner += 0.05 * ysize;
  xsize *= 0.9;
  ysize *= 0.9;
  (*scale) = ysize / oldysize;
  if (xsize / oldxsize < (*scale))
    (*scale) = xsize / oldxsize;
  (*xo) = (xcorner + (xsize - oldxsize * (*scale)) / 2.0) / (*scale);
  (*yo) = (ycorner   + (ysize - oldysize * (*scale)) / 2.0) / (*scale);

  xscale = (*scale) * xunitspercm;
  yscale = (*scale) * yunitspercm;
  initplotter(ntips,fn);
  plot(penup, xscale * (*xo), yscale * (*yo));
  plot(pendown, xscale * (*xo), yscale * ((*yo) + oldysize));
  plot(pendown, xscale * ((*xo) + oldxsize), yscale * ((*yo) + oldysize));
  plot(pendown, xscale * ((*xo) + oldxsize), yscale * (*yo));
  plot(pendown, xscale * (*xo), yscale * (*yo));
  /* we've done the extent, now draw the dividing lines: */
  xpag = (int)((pagex-hpmargin-0.01)/(paperx - hpmargin))+1;
  ypag = (int)((pagey-vpmargin-0.01)/(papery - vpmargin))+1;
  for (i=0;i<xpag;++i){
    plot(penup,(xscale * (*xo))+xscale*i*(paperx - hpmargin),((*yo)*yscale)+0);
    plot(pendown,(xscale * (*xo))+xscale*i*(paperx - hpmargin),((*yo)*yscale)+yscale*pagey);
    }
  for (j=0;j<ypag;++j){
    plot(penup,(xscale * (*xo)),((*yo)*yscale)+yscale*j*(papery-vpmargin));
    plot(pendown,(xscale * (*xo))+xscale*pagex,((*yo)*yscale)+yscale*j*(papery-hpmargin));
    }
}  /* makebox_no_interaction */


void void_func()
{
    fprintf(plotfile, "// Declare the colors\n\n");
    fprintf(plotfile, "#declare C_White       = color rgb<1, 1, 1>\n");

    fprintf(plotfile, "#declare C_White_trans = color rgbt<1, 1, 1, 0.7>\n");

    fprintf(plotfile, "#declare C_Red         = color rgb<1, 0, 0>\n");
    
    fprintf(plotfile, "#declare C_Yellow      = color rgb<1, 1, 0>\n");
    
    fprintf(plotfile, "#declare C_Green       = color rgb<0, 1, 0>\n");
    
    fprintf(plotfile, "#declare C_Black       = color rgb<0, 0, 0>\n");
    
    fprintf(plotfile, "#declare C_Blue        = color rgb<0, 0, 1>\n");
    
    fprintf(plotfile, "\n// Declare the textures\n\n");
    fprintf(plotfile, "#declare T_White = texture { pigment { C_White }}\n");
    fprintf(plotfile, "#declare T_White_trans = texture { pigment { C_White_trans }}\n");
    fprintf(plotfile, "#declare T_Red = texture { pigment { C_Red }\n");
    fprintf(plotfile, "\tfinish { phong 1 phong_size 100 }}\n");
    fprintf(plotfile, "#declare T_Red_trans = texture { pigment { C_Red filter 0.7 }\n");
    fprintf(plotfile, "\tfinish { phong 1 phong_size 100 }}\n");
    fprintf(plotfile, "#declare T_Green = texture { pigment { C_Green }\n");
    fprintf(plotfile, "\tfinish { phong 1 phong_size 100 }}\n");

    fprintf(plotfile, "#declare T_Green_trans = texture { \n");
    fprintf(plotfile, "\tpigment { C_Green filter 0.7 }\n");
    fprintf(plotfile, "\tfinish { phong 1 phong_size 100 }}\n");

    fprintf(plotfile, "#declare T_Blue = texture { pigment { C_Blue }\n");
    fprintf(plotfile, "\tfinish { phong 1 phong_size 100 }}\n");

    fprintf(plotfile, "#background { color rgb<1, 1, 1> }\n");
}  /* void_func */


/* added for vrml - danieyek 981111 */
/* Returned angle in radian */
/* A related function is "double angleBetVectors(Xu, Yu, Xv, Yv)"
   in drawtree.c */
double computeAngle(double oldx, double oldy, double newx, double newy)
{
  double angle;

  if ((newx-oldx) == 0 )
  {
    /* pi/2 or -pi/2! */
    if (newy > oldy) angle = pi/2;
    else if (newy < oldy) angle = -pi/2;
    else 
    {
      /* added - danieyek 990130 */
      /* newx = oldx; newy = oldy; one point on top of the other! 
         If new and old correspond to 2 points, changes are that the 2 coordinates
         are not identical under double precision value. */
      fprintf(stderr, 
      "ERROR: Angle can't be computed, 2 points on top of each other in computeAngle()!\n");
      angle = 0;
    }
  }
  else
  {
    angle = atan( (newy-oldy)/(newx-oldx) );

    if (newy >= oldy && newx >= oldx)
    {
      /* First quardrant - no adjustment */
    }
    else if (newx <= oldx)
    {
      /* Second (angle = negative) and
         third (angle = positive) quardrant */
      angle = pi + angle;
    }
    else if (newy <= oldy && newx >= oldx)
    {
      /* Fourth quardrant; "angle" is negative! */
      angle = 2*pi + angle;
    }
    else
    {
      /* Should never get here. */
      fprintf(stderr, "ERROR: Programming error in computeAngle()!\n");
    }
  }
  return angle;
}  /* computeAngle */

#ifdef WIN32
#include <windows.h>

/*********************  Prototypes  ***********************/

LRESULT WINAPI MainWndProc( HWND, UINT, WPARAM, LPARAM );
LRESULT WINAPI AboutDlgProc( HWND, UINT, WPARAM, LPARAM );

/*******************  Global Variables ********************/
extern void winplotpreviewcore();
HANDLE ghInstance;
HPEN hPenTree, hPenLabel, hPenBackground, hPenOld; 

/********************************************************************\
*  Comments: Register window class, create and display the main      *
*            window, and enter message loop.                         *
\********************************************************************/

winplotpreview()
{
   WNDCLASS wc;
   MSG msg;
   HWND hWnd;
   int screenXres, screenYres, winXres, winYres;

   winaction = quitnow;

   wc.lpszClassName = "GenericAppClass";
   wc.lpfnWndProc = MainWndProc;
   wc.style = CS_OWNDC | CS_VREDRAW | CS_HREDRAW;
   wc.hInstance = NULL;
   wc.hIcon = LoadIcon( NULL, IDI_APPLICATION );
   wc.hCursor = LoadCursor( NULL, IDC_ARROW );
   wc.hbrBackground = (HBRUSH)( COLOR_WINDOW+1 );
   wc.lpszMenuName = "GenericAppMenu";
   wc.cbClsExtra = 0;
   wc.cbWndExtra = 0;

   RegisterClass( &wc );

   ghInstance = NULL;

   screenXres = GetSystemMetrics(SM_CXSCREEN);
   winXres = (int)((float)(screenXres)*XWINPERCENT);
   screenYres = GetSystemMetrics(SM_CYSCREEN);
   winYres = (int)((float)(screenYres)*YWINPERCENT);

   hWnd = CreateWindow( "GenericAppClass",
      "Tree Preview",
      WS_OVERLAPPEDWINDOW,
      0,
      0,
      winXres,
      winYres,
      NULL,
      NULL,
      NULL,
      NULL
   );

   ShowWindow( hWnd, SW_SHOWNORMAL );

   while( GetMessage( &msg, NULL, 0, 0 ) ) {
      TranslateMessage( &msg );
      DispatchMessage( &msg );
   }

   return msg.wParam;
}

/*********************     *
*                                                                    *
* Comments: The following messages are processed                     *
*                                                                    *
*           WM_PAINT                                                 *
*           WM_COMMAND                                               *
*           WM_DESTROY                                               *
*                                                                    *
*                                                                    *
\********************************************************************/

LRESULT CALLBACK MainWndProc( HWND hWnd, UINT msg, WPARAM wParam,
   LPARAM lParam )
{
   PAINTSTRUCT ps;
   LOGBRUSH lb;
   HBRUSH bgbrush; 
   RECT lpRect;
   int windowwidth, windowheight;

   switch( msg ) {
/**************************************************************\
*     WM_ACTIVATE:                                             *
\**************************************************************/

      case WM_ACTIVATE:
         if (wParam != WA_INACTIVE)
            BringWindowToTop(hWnd);
      break;

/**************************************************************\
*     WM_PAINT:                                                *
\**************************************************************/

      case WM_PAINT:
         hdc = BeginPaint( hWnd, &ps );
         /* Initialize the pen's brush. */
         lb.lbStyle = BS_SOLID; 
         lb.lbColor = RGB(0,0,0); 
         lb.lbHatch = 0;
         /* 2 pixel pen for the tree */
         hPenTree = ExtCreatePen(PS_GEOMETRIC | PS_SOLID | PS_ENDCAP_ROUND,
                            (DWORD)2, &lb, 0, NULL); 
         /* 1 pixel pen for labels */
         hPenLabel = ExtCreatePen(PS_GEOMETRIC | PS_SOLID | PS_ENDCAP_ROUND,
                            (DWORD)1, &lb, 0, NULL); 
         /* light blue pen for outline of background rectangle */
         lb.lbColor = RGB(204,255,255);
         hPenBackground = ExtCreatePen(PS_GEOMETRIC | PS_SOLID,
                            (DWORD)1, &lb, 0, NULL); 
         /* light blue brush for interior of background rectangle */
         bgbrush = CreateSolidBrush(RGB(204,255,255));
         /* GetClientRect returns the size of that part of the window
            that is actually ours to draw in.
          */
         GetClientRect(hWnd, &lpRect);
         windowwidth = lpRect.right;
         windowheight = lpRect.bottom;
         /* select background pen and brush */ 
         SelectObject(hdc, hPenBackground);
         SelectObject(hdc, bgbrush);
         /* fill background */
         Rectangle(hdc, 0, 0, windowwidth, windowheight);
         /* select tree pen */
         hPenOld = SelectObject(hdc, hPenTree); 
         /* winplotpreviewcore calls makebox, plottree, plotlabels and
            finishplotter
          */
         winplotpreviewcore(windowwidth, windowheight);
         /* delete pens to recover memory */
         DeleteObject(hPenTree); 
         DeleteObject(hPenLabel);
         DeleteObject(hPenBackground);
         DeleteObject(bgbrush);
         EndPaint( hWnd, &ps );
         break;

/**************************************************************\
*     WM_COMMAND:                                              *
\**************************************************************/

      case WM_COMMAND:
         switch( wParam ) {
            case IDM_ABOUT:
               DialogBox( ghInstance, "AboutDlg", hWnd, (DLGPROC)
                          AboutDlgProc );
            break;
            case IDM_PLOT: // "Plot" menu item
               winaction = plotnow;
               DestroyWindow(hWnd);
            break;
            case IDM_CHANGE: // "Change Parameters" menu item
               winaction = changeparms;
               DestroyWindow(hWnd);
            break;
            case IDM_QUIT: // "Quit" menu item
               winaction = quitnow;
               DestroyWindow(hWnd);
            break;
         }
      break;

/**************************************************************\
*     WM_DESTROY: PostQuitMessage() is called                  *
\**************************************************************/

      case WM_DESTROY:
         PostQuitMessage( 0 );
         break;

/**************************************************************\
*     Let the default window proc handle all other messages    *
\**************************************************************/

      default:
         return( DefWindowProc( hWnd, msg, wParam, lParam ));
   }

   return 0;
}

/********************************************************************\
* Function: LRESULT CALLBACK AboutDlgProc(HWND, UINT, WPARAM, LPARAM)*
*                                                                    *
*  Purpose: Processes "About" Dialog Box Messages                    *
*                                                                    *
* Comments: The About dialog box is displayed when the user clicks   *
*           About from the Help menu.                                *
*                                                                    *
\********************************************************************/

LRESULT CALLBACK AboutDlgProc( HWND hDlg, UINT uMsg, WPARAM wParam, LPARAM lParam )
{
   switch( uMsg ) {
      case WM_INITDIALOG:
         return TRUE;
      case WM_COMMAND:
         switch( wParam ) {
            case IDOK:
               EndDialog( hDlg, TRUE );
               return TRUE;
         }
      break;
   }

   return FALSE;
}


#endif
