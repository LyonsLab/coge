/* stderr */
/* Interface
   version 3.5c. (c) Copyright 1992-2004 by the University of Washington.
   Written by Sean T. Lamont and Michal Palczewski.
   For use with the Macintosh version of the Phylogeny Inference Package,
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed.

   Functions you need to know how to use:
   macsetup(char *name):  initializes the interface, brings up a window of
                          the name of the argument, for I/O.
   textmode(); hides the graphics window
   gfxmode(); shows the graphics window.
 */
#ifdef __MWERKS__
#include <sioux.h>
#endif
#include <ctype.h>
#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#ifdef OSX_CARBON
#include <Carbon/Carbon.h>
void resize_gfx_window(EventRecord ev);
#ifndef __MWERKS__
OSErr CPSEnableForegroundOperation(ProcessSerialNumber *);
#endif
void makebox(char *, double *, double *, double *, long);
#endif
#include "draw.h"
#include "interface.h"
#define MAX(a,b) (a) > (b) ? (a) : (b)

#define TEXT 0
#define GFX 1

void do_about (void);
void handle_dlg_event(WindowPtr dlg);
void paint_dialog(WindowPtr about_dialog);

extern winactiontype winaction;

extern long winheight;
extern long winwidth;
Rect rect = { 0, 0, 16000, 16000 }; /* a nice big rect very convenient */
extern char* about_message;

/* These are all external variables from other files that this program needs
to access in order to draw and resize */
#define boolean char
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
 
extern double oldx, oldy;

extern boolean didloadmetric;
extern long   nmoves,oldpictint,pagecount;
extern double labelline,linewidth,oldxhigh,oldxlow,oldyhigh,oldylow,
       raylinewidth,treeline,oldxsize,oldysize,oldxunitspercm,
       oldyunitspercm,oldxcorner,oldycorner,oldxmargin,oldymargin,
       oldhpmargin,oldvpmargin,clipx0,clipx1,clipy0,clipy1,userxsize,userysize;
extern long rootmatrix[51][51];
extern long  HiMode,GraphDriver,GraphMode,LoMode,bytewrite;


void handlemouse(WindowPtr win,EventRecord ev);

/* Global variables used in many functions*/

#ifdef OSX_CARBON
/* Process Serial Number in text mode. Typically the OSX Terminal application */
ProcessSerialNumber textModePSN;
Boolean textModePSN_valid = false;
#endif 

int mode = TEXT;
int quitmac;
WindowPtr gfx_window;
ControlHandle plot_button;
ControlHandle change_button;
ControlHandle quit_button;
ControlHandle about_button;
Rect gfxBounds = { 50, 10, 400, 260 };        /* position and size of gfx_window */
Rect resizeBounds = {100+MAC_OFFSET,340,16000,16000}; /* how much the window can be resized*/

RGBColor background; 
RGBColor foreground;

/* saved parameters needed to make a call to makebox*/
mpreviewparams macpreviewparms;

/* initialize general stuff*/
void
macsetup (char *tname, char *gname)
{
  Str255 buf1, buf2;                        
  /* Str255 title = ""; */
  Rect plot_rec={10,10,30,110}; 
  Rect change_rec={10,120,60,220};
  Rect quit_rec={40,10,60,110};
  Rect about_rec={10,230,30,330};
  unsigned char* plot=(unsigned char*)"\4Plot";
  unsigned char* change=(unsigned char*)"\7Change\rParameters";
  unsigned char* quit=(unsigned char*)"\4Quit";
  unsigned char* about=(unsigned char*)"\5About";

#ifndef __MWERKS__
  ProcessSerialNumber myProcess;
#endif

#ifdef OSX_CARBON
  OSStatus err;
  BitMap screenBits;

  /* Save the frontmost process, which should be the process containing the console. */
  err = GetFrontProcess(&textModePSN);
  if (err != noErr)
    printf("Cannot get PSN of front process: %d\n", err);

  GetQDGlobalsScreenBits(&screenBits);
#endif

#ifndef OSX_CARBON
  change[0]=0x11;
#endif

  background.red=0xcc00; /* #ccffff phylip color */
  background.green=0xffff;
  background.blue=0xffff;
  
#undef fontsize

#ifndef OSX_CARBON  /* this looks bad in OS X */
  SIOUXSettings.fontsize= log(qd.screenBits.bounds.right);
#endif

#ifdef __MWERKS__
  SIOUXSettings.autocloseonquit = true;
#endif

  putchar('\n'); /* initialize sioux and let sioux initialize the toolbox and menus*/
  
  strcpy ((char *) buf1 + 1, tname);
  strcpy ((char *) buf2 + 1, gname);
  buf1[0] = strlen (tname);
  buf2[0] = strlen (gname);
  
#ifdef OSX_CARBON
  gfxBounds.bottom=screenBits.bounds.bottom*.7;
#else
  gfxBounds.bottom=qd.screenBits.bounds.bottom*.7;
#endif
  gfxBounds.right=MAX((gfxBounds.bottom-MAC_OFFSET)*.7,340);
  winheight=gfxBounds.bottom-gfxBounds.top-MAC_OFFSET;
  winwidth=gfxBounds.right-gfxBounds.left;

#ifndef __MWERKS__
  /* BUGBUG: This is an undocumented hack so that our windows will work 
     correctly as an unbunbled console app -db */
  GetCurrentProcess(&myProcess);
  CPSEnableForegroundOperation(&myProcess);
#endif

#ifdef OSX_CARBON
  CreateNewWindow(kDocumentWindowClass,
                  kWindowCloseBoxAttribute | kWindowFullZoomAttribute | 
                  kWindowCollapseBoxAttribute | kWindowResizableAttribute
                  ,&gfxBounds, &gfx_window);
#else
  gfx_window = NewCWindow (0L, &gfxBounds, buf2, false, documentProc,
                           (WindowPtr) - 1L, true, 0);
#endif

  plot_button = NewControl(gfx_window,&plot_rec,plot,1,0,0,1,pushButProc,0);
  change_button = NewControl(gfx_window,&change_rec,change,1,0,0,1,pushButProc,0);
  quit_button = NewControl(gfx_window,&quit_rec,quit,1,0,0,1,pushButProc,0);
  about_button = NewControl(gfx_window,&about_rec,about,1,0,0,1,pushButProc,0);
  
  foreground.red=0x0000;  /* black foreground */
  foreground.green=0x0000;
  foreground.blue=0x0000;
}

/* event loop for the preview window */
void
eventloop ()
{
  int status=1;
  quitmac=0;

  while (status > 0 && quitmac == 0)
        {
          status = handleevent ();
          if (status <= 0 || quitmac)
                textmode();

        }
}


/* event handler */
int
handleevent ()
{
  EventRecord ev;
  WindowPtr win;
#ifdef __MWERKS__
  int SIOUXDidEvent;
#endif
int where, ok;

#ifdef OSX_CARBON
  ok = WaitNextEvent (everyEvent, &ev, 0x7FFFFFFF, NULL);
#else
  ok = GetNextEvent (everyEvent, &ev);
#endif
  if (!ok) return 1;
  where = FindWindow (ev.where, &win);
#ifdef __MWERKS__
  if (win != gfx_window || where  == inMenuBar) {
          SIOUXDidEvent = SIOUXHandleOneEvent(&ev);
          if (SIOUXDidEvent) return 1;
  }
#endif
  if (win != gfx_window) return 1;
  
  if ((ev.what == keyDown) &&
          (ev.modifiers & cmdKey) && 
          (toupper ((char) (ev.message & charCodeMask)) == 'W'))
        return 0;
  else if (ev.what == activateEvt) {
#ifdef OSX_CARBON
    InvalWindowRect(gfx_window, &rect);
#else
    InvalRect(&rect);
#endif
  }
  else if (ev.what == updateEvt && win == gfx_window )
    paint_gfx_window();
  else if (ev.what == mouseDown && where == inContent) 
    handlemouse(win,ev);
#ifndef OSX_CARBON
  else if (ev.what == mouseDown && where == inSysWindow)
    SystemClick (&ev, win);
#endif
  else if (ev.what == mouseDown && where == inDrag)
    DragWindow (win, ev.where, &rect);
  else if (ev.what == mouseDown && where == inGrow) 
    resize_gfx_window(ev);
  else if (ev.what == mouseDown && where == inGoAway)
    if ( TrackGoAway( win, ev.where ) ) {
      winaction = changeparms;
      return 0;        
    }
  if (ev.what == mouseDown ) {
    SelectWindow(win);
  }

  return 1;
}

/*Handle mouse down event */
void handlemouse(WindowPtr win,EventRecord ev) {
        Point mouse;
        ControlHandle control;
        int part;
        
        
        mouse=ev.where;
        GlobalToLocal(&mouse);
        part=FindControl(mouse,win,&control);
        if (part != 10) return;
        TrackControl(control,mouse,NULL);
        if (control == plot_button) {
                quitmac=1;
                winaction=plotnow;
        } else if (control == quit_button){
                quitmac=1;
            winaction=quitnow;
        } else if (control == change_button) {
                winaction = changeparms;
                quitmac=1;
        } else if (control == about_button) {
                do_about();
        }
}

/* Draw a string to the graphics window */
void
putstring (string)
         char *string;
{
  unsigned char buf[256];
  strncpy ((char *) buf + 1, string, 253);        
  buf[0] = strlen (string);
  DrawString (buf);
}

/* go into text mode */
void
textmode ()
{
#ifdef OSX_CARBON
  ProcessSerialNumber psn;
  const Boolean visible = false;
  OSErr err;
  ProcessInfoRec info;
  unsigned char processName[32];

  memset(processName, '\0', 32);

#  ifdef DEBUG
  /* Print the process name which we have saved. (Should be "Terminal" or somesuch.) */
  if (textModePSN_valid == false)
    assert(0);

  info.processInfoLength = sizeof(ProcessInfoRec);
  info.processName = processName;
  info.processAppSpec = NULL;

  err = GetProcessInformation(&textModePSN, &info);
  if (err == noErr) {
    printf("Text mode process: \"%s\" %d\n", info.processName + 1, info.processName[0]);
  }
#  endif

  /* Bring the terminal process to the front again */
  
  err = SetFrontProcess(&textModePSN);
  if (err != noErr)
    printf("Warning: cannot bring text-mode process to front: %d\n", err);
#endif

  /* Hide the graphics window */
  HideWindow (gfx_window);

  mode = TEXT;
}

/* go into graphics mode */
void gfxmode() {
  OSErr err;
  ProcessSerialNumber psn;

  InitCursor();

#ifdef OSX_CARBON
  SetPort(GetWindowPort(gfx_window));
#else
  SetPort(gfx_window);
#endif

  /* Display the graphics window and bring it to the front of the application */
  ShowWindow(gfx_window);
  SelectWindow(gfx_window);

  /* Bring our process to the foreground, in case we were operating from Terminal.app */
  err = GetCurrentProcess(&psn);
  if (err == noErr) {
    ShowHideProcess(&psn, true);
    SetFrontProcess(&psn);
  }

  mode = GFX;
}


/*call this function to paint the graphics window*/
void paint_gfx_window () {

  BeginUpdate(gfx_window);
  RGBBackColor(&background);
  RGBForeColor(&foreground);
  EraseRect(&rect);
  
  PenSize(1,1);        
  makebox(macpreviewparms.fn,
          macpreviewparms.xo,
          macpreviewparms.yo,
          macpreviewparms.scale,
          macpreviewparms.nt);
  if (linewidth < 5) linewidth = 5;
  PenSize(linewidth/5,linewidth/5);
  plottree(macpreviewparms.root, macpreviewparms.root);
  plotlabels(macpreviewparms.fn);
  /* changed to support OSX.  Seems to work fine w/ OS 9 -db
     UpdateControls(gfx_window,gfx_window->visRgn); */
  DrawControls(gfx_window);
  EndUpdate(gfx_window);
  
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

/* resize the graphics window*/
void resize_gfx_window(EventRecord ev) {
  long windowsize;
#ifdef OSX_CARBON
  Rect newBox;
  windowsize=ResizeWindow(gfx_window,ev.where,&resizeBounds, &newBox);
  if (windowsize != 0 ) {
    winheight=newBox.bottom-newBox.top - MAC_OFFSET;
    winwidth=newBox.right-newBox.left;
    InvalWindowRect(gfx_window, &rect);
  }
#else
  windowsize=GrowWindow(gfx_window,ev.where,&resizeBounds);
  if (windowsize != 0 ) {
    SizeWindow(gfx_window,LoWord(windowsize),HiWord(windowsize),TRUE);
    winheight=HiWord(windowsize)-MAC_OFFSET;
    winwidth=LoWord(windowsize);
    InvalRect(&rect);
  }
#endif
}




void do_about () {
  WindowPtr about_dialog;
  ControlHandle okbutton;
  short txt;
  
  Rect gfxBounds = {50,50,300,600};
  /* Rect labelRect = {25,15,150,400}; */
  Rect buttonRect = {200,350,220,450};
#ifdef OSX_CARBON
  CreateNewWindow (kMovableModalWindowClass,0, &gfxBounds, &about_dialog);
#else
    about_dialog = NewCWindow (0L, &gfxBounds, (unsigned char*)"\5About", false, dBoxProc,
                           (WindowPtr) - 1L, true, 0);
#endif
  SetWTitle(about_dialog,(unsigned char*)"\5About");
  okbutton = NewControl(about_dialog, &buttonRect,(unsigned char*)"\2Ok",1,0,0,1,kControlPushButtonProc,101);
  ShowWindow(about_dialog);
#ifdef OSX_CARBON
  SetPort(GetWindowPort(about_dialog));
#else
  SetPort(about_dialog);
#endif
  paint_dialog(about_dialog);
  handle_dlg_event(about_dialog);
  HideWindow(about_dialog);
  DisposeWindow(about_dialog);
#ifdef OSX_CARBON
  SetPort(GetWindowPort(gfx_window));
#else
  SetPort(gfx_window);
#endif
}


void paint_dialog(WindowPtr about_dialog) {
  /* int i; */
  char *cur;
  char *what,*tmpstr,*orig;
  int go=1;
  int curx = 10,cury = 20;

  /* No strdup on mac! need a editable string */
  cur = malloc(sizeof(char) * (strlen(about_message) + 1));
  strcpy(cur,about_message);
  orig = cur;


  BeginUpdate(about_dialog);
  MoveTo(curx,cury);
  while (go) {
    what = strchr(cur,'\r');
    if ( what == NULL) {
      what = cur;
      go = 0;
    }
    else {
      *what = '\0';
      tmpstr = cur;
      cur = what;
      what = tmpstr;
      cur += 1;
    }
    putstring(what);
    cury += 20;
    MoveTo(curx,cury);
  }
  free(orig);
  DrawControls(about_dialog);
  EndUpdate(about_dialog);
}


void handle_dlg_event(WindowPtr dlg)
{
  EventRecord ev;
  int ok,where;
  WindowPtr win;
  char key;
  char keycode;
  ControlHandle h;
  int go=1;
  short w;
  /* int i; */
  Point mouse;

  while (go) {
    ok = GetNextEvent(everyEvent,&ev);
    where = FindWindow(ev.where,&win);
    if (ev.what == mouseDown && (win == dlg)) {
      mouse = ev.where;
#ifdef OSX_CARBON
      SetPort(GetWindowPort(dlg));
#else
          SetPort(dlg);
#endif
      GlobalToLocal(&mouse);
      if (where == inSysWindow) {
      /*  SystemClick(&ev,win);*/
      }
      else if (where == inDrag)
            DragWindow(win, ev.where, &rect);
      else {
        FindControl(mouse,dlg,&h);
        w = TrackControl(h,mouse,NULL);
        if ( h == NULL ) continue;
        
        if ( w == 0) continue;
        /*i = h->contrlRfCon ;
        if (i == 101) */
          go = 0;
      }
    }
    else if (ev.what == updateEvt && win == dlg)
      paint_dialog(dlg);
    else if (ev.what == updateEvt && win == gfx_window)
      paint_gfx_window();
    else if (ev.what == keyDown) {
      key = (char)(ev.message & charCodeMask);
      keycode = (short)((ev.message &keyCodeMask) >> 8);
      if ( key == '\r' ) {
        go = 0;
      }
    }
  }
}
