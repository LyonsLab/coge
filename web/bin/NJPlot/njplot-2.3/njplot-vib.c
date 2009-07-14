/* 

one of WIN_MOTIF  WIN_MAC  WIN_MSWIN must be defined for normal compilation

under WIN_MOTIF, one can also define NO_PDF so the PDFLibLite library is not used,
and one can also define HELPFILENAME to the full pathname (with enclosing "") of the help file

Special compilations:

compile with -DTTY to have it output text-only tree without GUI

compile with -DADDROOT to have it root the unrooted input tree and write tree on stdout (no GUI)

compile with -DNO_GUI to have it work without GUI at all and convert tree to PDF/PostScript file 


WIN_MOTIF resources:
Vibrant.systemfont : helvetica,14,b
Njplot.systemfont : times,14,b
Vibrant.programfont : fixed,14
placed in file
$HOME/.Xdefaults
$XAPPLRESDIR/Vibrant
/usr/X11R6/lib/X11/app-defaults/Vibrant    (Linux)
/usr/openwin/lib/X11/app-defaults/Vibrant   (Solaris)
*/
#define NJPLOTVERSION "2.3"

#include <math.h>
#include <time.h>
#include <ctype.h>


#if defined(WIN_MSWIN)  
#define WITH_PDF 1
#elif defined(WIN_MOTIF)  && ! defined(NO_PDF)
#define WITH_PDF 1
#elif defined(NO_GUI)  && ! defined(NO_PDF)
#define WITH_PDF 1
#endif


#ifdef WITH_PDF
#include "pdflib.h"
PDF *pdf = NULL;
int pdf_font;
#endif
#define MAX_FRAC 0.95

#if defined(TTY) || defined (ADDROOT)
#define NO_GUI
#endif

#ifdef NO_GUI

#define FALSE 0
#define TRUE (!FALSE)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#else

#ifdef WIN_MAC
#include <Carbon.h> /* must be before <vib*.h> or it does not compile */
#endif
#include <vibtypes.h>
#include <vibprocs.h>

#if defined(WIN_MSWIN) || defined (WIN_MAC)
#define  MICRO
#endif

#ifdef WIN_MAC
#define myWaitAndProcessNextEvent Nlm_ProcessAnEvent
extern int Nlm_textScrapFull;

#elif defined(WIN_MSWIN)
#define HDC void *
extern HDC Nlm_currentHDC;
extern HDC Nlm_GetPicWinHDC ( void );
char *get_prog_dir(void);
#define myWaitAndProcessNextEvent Nlm_ProcessAnEvent
#define tempnam _tempnam

#else  /* unix */

#undef Boolean
#include <X11/Intrinsic.h>

extern XtAppContext Nlm_appContext;
extern Nlm_CharPtr Nlm_XrmGetResource(const Nlm_Char PNTR _resource);
void set_systemfont(char *systemfont);

void myWaitAndProcessNextEvent (void)
{
XEvent       event;
if( XtAppPending(Nlm_appContext) & XtIMXEvent ) {
	Nlm_ProcessAnEvent();
	}
else	{
/* wait until next event and catch it */
    	XtAppNextEvent (Nlm_appContext, &event);
    	XtDispatchEvent (&event); 
	}
}

#endif


#include <document.h> /* need this include after the other mac includes */

#ifdef WIN_MAC /* need the mac version of each of these functions */
#undef MoveTo
#undef LineTo
#undef ClipRect
#undef TextWidth
#undef DrawString
#endif

struct win_help_extra {
	int tag;
	DoC document;
	FonT police;
	};

struct my_RecT4 {
	Nlm_Int4 left, top, right, bottom;
	}  print_rect;
#endif /* NO_GUI */

/* typedefs */

typedef enum { times = 1, helvetica, courier} font_name;
enum { font_tiny = 1, font_small,
	font_normal, font_medium, font_large, font_bold, font_italic };
typedef enum { A4 = 1, LETTER } 
	paper_item;

struct noeud {
	double l1,l2,l3;
	struct noeud *v1,*v2,*v3;
	};
	
typedef struct branche { /* une branche definie par ses deux extremites */
	struct noeud *bouta;
	struct noeud *boutb;
	char *br_label;
	} branche;

#define s_noeud sizeof(struct noeud)
#define myrint(d) ((int)(d + 0.5)) /* not correct if d < 0 but irrelevant for plotting */
struct nom {
	double x,y;
	char *nom;
	char disp_option;
	};
struct mon_point {
	double x,y;
	int number;
	};
struct trait {
	double xd,yd,xf,yf;
	};

typedef struct _moments {
	int N; /* nbre de fils */
	double somme, carres; /* somme des longueurs jusqu'aux feuilles, et somme des carres */
	} moments;
enum operations {show_tree = 1, depl_racine, permutation, subtree};

typedef struct {
#ifndef NO_GUI
	int tag;
	Nlm_WindoW nj_plot;
	Nlm_IteM open_button;
	Nlm_PrompT click_message;
	Nlm_MenU tree_font_menu;
	Nlm_MenU paper_menu;
	Nlm_IteM win_menu_item;
	Nlm_IteM save_plot_button;
	Nlm_IteM save_tree_button;
	Nlm_IteM save_unrooted_button;
	Nlm_IteM SearchItem;
	Nlm_IteM AgainItem;
	Nlm_ChoicE choix_font, choix_taille;
	Nlm_IteM bold_item, italic_item;
	Nlm_ButtoN branch_length_button;
	Nlm_ButtoN bootstrap_button;
	Nlm_GrouP choix_group;
	Nlm_ButtoN reroot_button;
	Nlm_ButtoN subtree_up_button;
	Nlm_IteM exit_button;
	Nlm_PaneL tree_plot;
	float zoomvalue;
	Nlm_PrompT zoomprompt;
	Nlm_RecT old_rect;
	Nlm_FonT current_font;
#endif /* NO_GUI */
	char *tree_name;
	int notu;
	int totbranches;
	int has_br_length;
	int has_internal;
	int root_num;
	int rooted;
	struct noeud **tabtax;
	struct noeud *racine;
	char **tabnames;
	char **labels;
	struct branche *branches;
	struct nom *noms;
	struct trait *traits;
	struct mon_point *points;
	char *br_length_txt;
	int totnoms;
	int totpoints;
	int tottraits;
	double deltay;
	int show_bootstrap;
	int need_runtree;
	enum operations choix;
	int *widnames;
	double *profs;
	int long_arbre_parenth; /* long de la chaine decrivant l'arbre lu */
	struct noeud *subtree_center, *subtree_ascend;
	int subtree_notu;
	int plot_br_l;
	double root_br_l;
	int char_height;
	int ascent;
	int font_size_rank;
	int font_family_rank;
	char font_bold_italic[3];
	} FD_nj_plot;


/* prototypes of included functions */
#ifndef NO_GUI
FD_nj_plot *create_win_nj_plot(void);
void cre_help_win(IteM item);
FILE *open_path(char *fname);
void topic_callback(LisT list);
void help_ok_callback(ButtoN bouton);
#ifdef WIN_MAC
void show_apropos_njplot(Nlm_IteM item);
void mac_timer(void);
void window_callback(IteM item);
#endif
void load_help_topics(FILE *fich, LisT topic_list);
void subtree_up_callback(ButtoN bouton);
 void open_callback(Nlm_IteM );
void new_callback(IteM ob);
 void font_callback(Nlm_ChoicE);
 void paper_callback(Nlm_ChoicE);
 void save_plot_callback(Nlm_IteM);
 void save_tree_callback(Nlm_IteM);
 void toggle_branch_callback(Nlm_ButtoN);
 void toggle_bootstrap_callback(Nlm_ButtoN);
 void operation_callback(Nlm_GrouP);
 void subtree_callback(void);
 void swap_callback(void);
 void new_outgroup_callback(void);
 void show_tree_callback(void);
 void exit_callback(Nlm_IteM);
 void tree_draw_proc(Nlm_PaneL);
 void tree_click_proc(Nlm_PaneL, Nlm_PoinT);
static void closeproc(WindoW);
void closefrommenu(Nlm_IteM item);
void win_resize_proc(WindoW);
void change_panel_size(PaneL panel, int change_x, int change_y);
void string_callback(ButtoN item);
void search_callback(IteM item);
void process_keys(Nlm_Char key);
void change_page_count(ButtoN item);
void page_count_callback(IteM item);
void clear_tree(IteM);
void scrollnotu(int count);
void paste_tree(IteM);
void scrollcallback(Nlm_BaR b, Nlm_SlatE g, Nlm_Int2 before, Nlm_Int2 after);
void zoomcallback(Nlm_SwitcH s, Nlm_Int2 before, Nlm_Int2 after);
#ifdef MICRO
void copy_plot(IteM);
void print_title(int x, int y, char *text, Nlm_FonT title_font, int p,int totp);
void print_plot(IteM);
#endif
#ifdef WIN_MAC
extern void add_apropos(char *progname, Nlm_ActnProc my_show_apropos);
void create_win_if_needed(char *fname);
int is_macosx(void);
int crefpict(char *, PicHandle );
void example_tree(IteM unused);
void pict_plot(void);
extern int MG_GetInputFName(char *fname, int maxl);
extern int MG_GetOutputFName(char *fname, int maxl, char *dfault);
extern char *mac_fname_to_roman(char *);
extern QDPictRef MyPictToQDPict(PicHandle mypicture);
extern int MyQDPictToPDFfile (QDPictRef picture, char *fname);
extern char *MG_GetBundleResourcesDir(void);
#elif defined(WIN_MSWIN)
extern int MG_GetInputFName(char *fname, int maxl);
extern int MG_GetOutputFName(char *fname, int maxl, char *dfault); 
#else 
#define MG_GetOutputFName Nlm_GetOutputFileName
#define MG_GetInputFName(a,b) Nlm_GetInputFileName(a, b, NULL, NULL)
#endif
extern Nlm_WindoW Nlm_GetNext(Nlm_WindoW);
#endif /* NO_GUI */

void process_args(int *argc, char *argv[], char **pred, int *pboot_arg, int *p_plot_br_l, int *p_font_size_rank);
void prepare_fonts(void);
void prepare_pdf_font(int fontnum, int use_bold, int use_italic);
int direct_pdf_plot(char *fname);
int calc_otu(struct noeud *pere, struct noeud *centre);
void do_plot(void);
int calc_text_size(char *text, int *pheight, int *pascent);
void plotstring(char *nom);
void dir_moveto(double x,double y);
void dir_lineto(double x,double y);
void scale_window(double lxmin, double lxmax, double lymin, double lymax,
 	double pxmin, double pxmax, double pymin, double pymax);
void ch_echelle(double lx, double ly, double *px, double *py);
void convert_mem_point(void);
void init_tree(char *fname, char *displayname);
char *preptree(char *fname);
char *check_alloc(int nbrelt, int sizelt);
void err_message(const char *text);
void loadphylip(char *arbre, char *last_bootstrap);
struct noeud *unrootedset(char *deb, char *fin, branche **p_int_br);
char *extract_filename(char *name);
char *nextpar(char *pospar);
void make_binary_or_unrooted(char *arbre);
int make_binary(char *arbre, char *debut, char *fin, int go_down);
void mydrawstring(double x, double y, char *nom, char option);
void moveto(double x,double y);
void lineto(double x,double y);
int calc_brl_for_lengthless(struct noeud *centre, struct noeud *pere);
void add_value_downstream(struct noeud *centre, int value);
void place_midpoint_root(void);
void parcourir_branches(struct noeud *centre, struct noeud *origine);
void process_branche(struct noeud *cote1, struct noeud *cote2, double length);
double get_length_down(struct noeud *, struct noeud *);
moments stat_from_node(struct noeud *pere, struct noeud *racine);
void runtree(void);
void mem_plot(struct noeud *pere, struct noeud *centre, double currx,
	double *curry);
void mem_point(double x, double y, int number);
void mem_nom(double x, double y, char *nom, char option);
void mem_trait(double xd, double yd, double xf, double yf);
char *get_br_label(struct noeud *a, struct noeud *b);
int get_br_from_bouts(struct noeud *a, struct noeud *b);
void free_tree(void);
char *ecrit_arbre_parenth_unrooted(struct noeud *root);
char *ecrit_arbre_parenth(struct noeud *root);
char *recur_ecrit_arbre(struct noeud *centre, char *arbre, char *finarbre);
void removeroot(void);
void remove_arg(int target, int *argc, char *argv[]);
double length_log_phys(double p);
double length_phys_log(double p);
double arrondi_echelle(double x);
double calc_echelle(double larg);
void draw_scale(void);
int calc_n_desc(struct noeud *pere);
void bad_format(char *);
void majuscules(char *p);
void postscript_plot(void); 
#ifdef WITH_PDF
void plot_to_pdf(void);
#endif
void tty_plot(void);




/* globals */
#ifndef NO_GUI
extern Nlm_WindoW Nlm_desktopWindow;
int tag_njplot, tag_help;
#endif
FD_nj_plot *fd_nj_plot;
int page_count = 1;
char *end_br_length;
double	physx,physy,physx_min,physy_min,physx_corr,physy_corr;
/* window scaling variables */
double tek_xmin,tek_xmax,tek_ymin,tek_ymax,tek_dx,tek_dy;
double tek_minx,tek_maxx,tek_miny,tek_maxy;
double maxx, maxy, nexty;
/* max length of taxa names */
int pdf_plot = 0, file_plot,  
	num_noeud, nextotu, no_title = FALSE,
	swap = 0;
char current_ps_font[40];
int list_font_size[] = { 0, 6, 8, 10, 12, 14, 18 };
#ifdef WIN_MSWIN
char list_font_name[4][40] = {"", "Times New Roman", "Lucida Sans", "Courier New"};
int default_font_size_rank = font_normal;
#else
char list_font_name[4][15] = {"", "Times", "Helvetica", "Courier"};
int default_font_size_rank = font_medium;
#endif
char list_ps_font_name[4][15] = {"", "Times", "Helvetica", "Courier"};
paper_item paper_choice = A4;
char plotfilename[200], rooted_fname[200];
FILE *plotfile;
FILE *help_file = NULL;
int converted;
int ps_width = 500, ps_height = 0 /* 0 means not defined by this variable */;
/* variables globales pour memoriser la meilleure branche:
celle qui partage l'arbre en 2 parties les plus egales possibles en prof moyenne
*/
double current_best_diff, current_br_length, current_balance;
struct noeud *current_cote1=NULL, *current_cote2=NULL;
#define VERY_BIG 9.e99

int doing_copy = FALSE;
int doing_print = FALSE;
float postscript_ratio = 0.60;
char *tty_prefix = NULL;

char **tty_page;
int tty_x, tty_y;
char *outotu = NULL; /* value of optional argument -outotu name */


#ifdef WIN_MAC
Rect		myrect;
PicHandle 	mypicture;
int 		margin;
#endif

#ifdef NO_GUI
int pdf_plot_only = TRUE;
int main(int argc, char **argv)
{
char *red, *fname; 
int boot_arg, plot_br_l;
fd_nj_plot = (FD_nj_plot *)calloc(1, sizeof(FD_nj_plot));

#else
int pdf_plot_only = FALSE;

#ifdef hpux
extern int argc;
extern char **argv;
#else
int argc;
char **argv;
#endif

Nlm_Int2 Nlm_Main(void)
{
char *red, *fname; 
int boot_arg, plot_br_l;
#ifndef hpux
argc = Nlm_GetArgc();
argv = Nlm_GetArgv();
#endif
#endif /* NO_GUI */

fname = NULL;
#ifdef WIN_MAC
if( ! is_macosx() ) argc = 1;
else if(argc > 1 && strncmp(argv[1], "-psn_", 5) == 0) remove_arg(1, &argc, argv); 
#endif
process_args(&argc, argv, &red, &boot_arg, &plot_br_l, &default_font_size_rank);
if(argc>=2) {
	fname = argv[1];
	}
else	{
	fname = NULL;
	}
#ifdef unix
putenv("LC_NUMERIC=C");
#endif
   if(pdf_plot_only) {
 	int status = 0;
	fd_nj_plot->tree_name = fname;
	fd_nj_plot->show_bootstrap = boot_arg;
	fd_nj_plot->plot_br_l = plot_br_l;
	fd_nj_plot->font_size_rank = default_font_size_rank;
	fd_nj_plot->font_family_rank = times;
 	if( fname != NULL) status = direct_pdf_plot(fname);
	return status;
	}
#ifndef NO_GUI
   Nlm_WatchCursor();
#ifdef WIN_MOTIF
	{char *font;
	font = Nlm_XrmGetResource("systemfont");
	if(font != NULL) set_systemfont(font);
	font = Nlm_XrmGetResource("programfont");
	if(font != NULL) Nlm_programFont = Nlm_ParseFont(font); 
	}
#endif
#ifdef WIN_MAC
add_apropos("njplot", show_apropos_njplot);
#endif
#ifdef WIN_MAC
	{
	char *p;
	p = MG_GetBundleResourcesDir();
	strcat(p, "/njplot.help");
	help_file = fopen(p, "r");
	}
#elif defined(HELPFILENAME)
	help_file = fopen(HELPFILENAME, "r");
#else
	help_file = open_path("njplot.help");
#endif
Nlm_KeyboardView(process_keys);
memcpy(&tag_njplot, "NJPL", sizeof(int));
memcpy(&tag_help, "HELP", sizeof(int));
   fd_nj_plot = create_win_nj_plot();
   SetStatus(fd_nj_plot->bootstrap_button, boot_arg);
   SetStatus(fd_nj_plot->branch_length_button, plot_br_l);
   if(fname != NULL) init_tree( fname, NULL );
   else SetTitle(fd_nj_plot->nj_plot, "njplot");
   Show(fd_nj_plot->nj_plot);
   Nlm_ObjectRect(fd_nj_plot->nj_plot, &(fd_nj_plot->old_rect));
   if(red != NULL) search_callback( (Nlm_IteM)red );
   Nlm_ArrowCursor();
#ifdef WIN_MAC
   Nlm_Metronome(mac_timer);
#endif
   ProcessEvents();
#endif
   return 0;
}


#ifndef NO_GUI
FD_nj_plot *create_win_nj_plot(void)
{
  WindoW win;
  MenU menu;
  ChoicE choix_paper;
  GrouP group;
  FD_nj_plot *fdui = (FD_nj_plot *) calloc(1, sizeof(FD_nj_plot));
  PoinT position;
  Nlm_Handle obj; int left, right;
  static int countwin = 0;
  
  
fdui->tag = tag_njplot;
fdui->nj_plot = win = DocumentWindow(-50 - countwin, -33 - countwin, -3, -3, "Njplot", 
		closeproc, win_resize_proc);
countwin += 3; countwin %= 15;
#ifdef WIN_MAC
static FD_nj_plot *mac_menu_bar_extra = NULL;
if(mac_menu_bar_extra != NULL) {
	fdui->open_button = mac_menu_bar_extra->open_button;
	fdui->save_plot_button = mac_menu_bar_extra->save_plot_button;
	fdui->save_tree_button = mac_menu_bar_extra->save_tree_button;
	fdui->save_unrooted_button = mac_menu_bar_extra->save_unrooted_button;
	fdui->SearchItem = mac_menu_bar_extra->SearchItem;
	fdui->AgainItem = mac_menu_bar_extra->AgainItem;
	fdui->tree_font_menu = mac_menu_bar_extra->tree_font_menu;
	fdui->choix_font = mac_menu_bar_extra->choix_font;
	fdui->choix_taille = mac_menu_bar_extra->choix_taille;
	fdui->bold_item = mac_menu_bar_extra->bold_item;
	fdui->italic_item = mac_menu_bar_extra->italic_item;
	fdui->paper_menu = mac_menu_bar_extra->paper_menu;
	goto after_menus;
	}
#define OPEN_LABEL      "Open        O"
#define NEW_LABEL       "New          N"
#define CLOSE_LABEL     "Close        W"
#define SAVE_PLOT_LABEL "Save as PDF   S"
#define PRINT_LABEL     "Print         P"
#define QUIT_LABEL      "Quit          Q"
#define COPY_LABEL      "Copy   C"
#define PASTE_LABEL     "Paste  V"
#define FIND_LABEL      "Find    F"
#define AGAIN_LABEL     "Again  A"
#define PLACE_OF_MENUS NULL
#else
#define OPEN_LABEL      "Open        ^O"
#define NEW_LABEL       "New         ^N"
#define CLOSE_LABEL     "Close        ^W"
#define PRINT_LABEL     "Print       ^P"
#define SAVE_PLOT_LABEL "Save as PDF   ^S"
#define QUIT_LABEL      "Quit          ^Q"
#define COPY_LABEL      "Copy     ^C"
#define PASTE_LABEL     "Paste     ^V"
#define FIND_LABEL      "Find       ^F"
#define AGAIN_LABEL     "Again     ^A"
#define PLACE_OF_MENUS fdui->nj_plot
#endif
  menu = PulldownMenu(PLACE_OF_MENUS, "File");
  fdui->open_button = CommandItem(menu, OPEN_LABEL, open_callback);
  CommandItem(menu, NEW_LABEL, new_callback);
fdui->save_plot_button = 
#if ( ! defined(WIN_MOTIF) ) || defined(WITH_PDF) 
  	CommandItem(menu, SAVE_PLOT_LABEL, save_plot_callback);
#endif
#ifdef WIN_MOTIF
  	CommandItem(menu, "Save as PostScript", save_plot_callback);
#endif
fdui->save_tree_button = CommandItem(menu, "Save Rooted Tree", save_tree_callback);
fdui->save_unrooted_button = CommandItem(menu, "Save Unrooted Tree", save_tree_callback);
#ifdef MICRO
  CommandItem(menu, PRINT_LABEL, print_plot);
#endif
  CommandItem(menu, CLOSE_LABEL, closefrommenu);
#ifdef WIN_MAC
if( ! is_macosx() )
#endif
  CommandItem(menu, QUIT_LABEL, exit_callback);
  Advance(PLACE_OF_MENUS);
  menu = PulldownMenu(PLACE_OF_MENUS, "Edit");
  CommandItem(menu, "Clear", clear_tree);
#ifdef MICRO
  CommandItem(menu, COPY_LABEL, copy_plot);
#endif
  CommandItem(menu, PASTE_LABEL, paste_tree);
#ifdef WIN_MAC
  SeparatorItem(menu);
  CommandItem(menu, "Example", example_tree);
#endif
  SeparatorItem(menu);
  fdui->SearchItem = CommandItem(menu, FIND_LABEL, search_callback);
  fdui->AgainItem = CommandItem(menu, AGAIN_LABEL, search_callback);
/*  Disable(fdui->AgainItem); */
  Advance(PLACE_OF_MENUS);
  fdui->tree_font_menu = PulldownMenu(PLACE_OF_MENUS, "Font"); 
  fdui->choix_font = ChoiceGroup(fdui->tree_font_menu, font_callback);
  Nlm_ChoiceItem(fdui->choix_font, list_font_name[1]);
  Nlm_ChoiceItem(fdui->choix_font, list_font_name[2]);
  Nlm_ChoiceItem(fdui->choix_font, list_font_name[3]);
  SeparatorItem(fdui->tree_font_menu);
  fdui->choix_taille = ChoiceGroup(fdui->tree_font_menu, font_callback);
  ChoiceItem(fdui->choix_taille, "6");
  ChoiceItem(fdui->choix_taille, "8");
  ChoiceItem(fdui->choix_taille, "10");
  ChoiceItem(fdui->choix_taille, "12");
  ChoiceItem(fdui->choix_taille, "14");
  ChoiceItem(fdui->choix_taille, "18");
  SeparatorItem(fdui->tree_font_menu);
  fdui->bold_item = StatusItem(fdui->tree_font_menu, "Bold", 
  	(Nlm_ItmActnProc) font_callback);
  fdui->italic_item = StatusItem(fdui->tree_font_menu, "Italic", 
  	(Nlm_ItmActnProc) font_callback);
  SetStatus(fdui->bold_item, FALSE);
  SetStatus(fdui->italic_item, FALSE);
  
  Advance(PLACE_OF_MENUS);
  fdui->paper_menu = PulldownMenu(PLACE_OF_MENUS, "Paper"); 
  choix_paper = ChoiceGroup(fdui->paper_menu, paper_callback);
  ChoiceItem(choix_paper, "A4");
  ChoiceItem(choix_paper, "US letter");
  SetValue(choix_paper, paper_choice);
  SeparatorItem(fdui->paper_menu);
#ifdef WIN_MSWIN /* SetTitle ne marche pas sur PC! */
  	CommandItem(fdui->paper_menu, "Page count", page_count_callback);
#elif defined(WIN_MAC) /* en 2 temps car () ds menu ne marche pas sur Mac ! */
  {	Nlm_IteM item; char nom[25];
  	item = CommandItem(fdui->paper_menu, "tmp", page_count_callback);
	sprintf(nom, "Page count (%d)", page_count);
  	SetTitle(item, nom);
  }
#else /* unix */
  {	char nom[25];
	sprintf(nom, "Page count (%d)", page_count);
  	CommandItem(fdui->paper_menu, nom, page_count_callback);
  }
#endif  
#ifdef WIN_MAC
  Advance(NULL);
  static Nlm_MenU window_menu;
  window_menu = PulldownMenu(NULL, "Window");
#endif
  Advance(PLACE_OF_MENUS);
  obj = PulldownMenu(PLACE_OF_MENUS, "Help");
  CommandItem(obj, "Help", cre_help_win);

  Break(PLACE_OF_MENUS);
#ifdef WIN_MAC
mac_menu_bar_extra = (FD_nj_plot *) calloc(1, sizeof(FD_nj_plot));
*mac_menu_bar_extra = *fdui;
after_menus:
	{static Nlm_IteM mac_win_item[50];
	static int mac_win_item_count = 0; int i;
	for(i=0; i<mac_win_item_count; i++) {
		if(!Nlm_Enabled(mac_win_item[i])) {
			Nlm_Enable(mac_win_item[i]);
			fdui->win_menu_item = mac_win_item[i];
			break;
			}
		}
	if(i>=mac_win_item_count) {
		fdui->win_menu_item = mac_win_item[mac_win_item_count] = CommandItem(window_menu, "njplot", window_callback);
		mac_win_item_count++;
		}
	}
#endif
  fdui->choix_group = NormalGroup(win, 4, 1, "Operation", programFont, operation_callback);
  RadioButton(fdui->choix_group, "Full tree");
  fdui->reroot_button = RadioButton(fdui->choix_group, "New outgroup");
  RadioButton(fdui->choix_group, "Swap nodes");
  RadioButton(fdui->choix_group, "Subtree");
  SetValue(fdui->choix_group, show_tree);

  Break(win);
  group = NormalGroup(win, 2, 1, "Display", programFont, NULL);
  fdui->branch_length_button = CheckBox(group, "Branch lengths", 
		toggle_branch_callback);
  SetStatus(fdui->branch_length_button, fdui->plot_br_l);
  fdui->bootstrap_button = CheckBox(group, "Bootstrap values", 
			toggle_bootstrap_callback);

  Advance(win);
  group = NormalGroup(win, 2, 1, "Zoom", programFont, NULL);
obj = Nlm_UpDownSwitch(group, 0, zoomcallback);
fdui->zoomvalue = 1.;
Nlm_SetSwitchParams( (Nlm_SwitcH)obj, 1, 100);
fdui->zoomprompt = Nlm_StaticPrompt(group, "100%", 35, Nlm_stdLineHeight, Nlm_programFont, 'l');

  Advance(win);
  GetNextPosition(win, &position);
  position.x += 5;
  position.y += 10;
  SetNextPosition(win, position);
  obj = fdui->subtree_up_button = PushButton(win, "Subtree Up", subtree_up_callback);
  Disable(fdui->subtree_up_button);
  

Break(win);
Nlm_GetPosition(obj, &(fdui->old_rect)); 
right = fdui->old_rect.right;
Nlm_GetNextPosition(win, &position); 
left = position.x;

#ifdef WIN_MSWIN
#define PLOT_HEIGHT 430
#elif defined(WIN_MAC)
#define PLOT_HEIGHT 300
#else
#define PLOT_HEIGHT 350
#endif

fdui->tree_plot = Nlm_AutonomousPanel(win, right - left + 1, PLOT_HEIGHT, 
		tree_draw_proc, scrollcallback, NULL, sizeof(int), NULL, NULL); 
	{int hidebar = TRUE; /* put TRUE iff scrollbar should be hidden */
	Nlm_SetPanelExtra(fdui->tree_plot, &hidebar);
	}
SetPanelClick(fdui->tree_plot, tree_click_proc, NULL, NULL, NULL);
obj = Nlm_GetSlateVScrollBar((SlatE)fdui->tree_plot);
Nlm_Hide(obj);
Nlm_SetBarValue(obj, 0);
Nlm_SetBarMax(obj, 1);

{int space = 0; /* some space at the bottom of the window */
#if defined(__APPLE__) || defined(WIN_MAC)
space = Nlm_hScrollBarHeight;
#elif defined(WIN_MOTIF)
space = 3;
#endif
if(space > 0) {
	Nlm_Break(win);
	Nlm_StaticPrompt(win, "", 4, space, Nlm_programFont, 'l');
	}
}

fdui->old_rect.left = -1;
SetWindowExtra(win, fdui, NULL);
fdui->font_size_rank = default_font_size_rank;
fdui->font_family_rank = times;
Nlm_SetValue(fdui->choix_taille, default_font_size_rank);
Nlm_SetValue(fdui->choix_font, times);
Nlm_SetStatus(fdui->bold_item, FALSE);
Nlm_SetStatus(fdui->italic_item, FALSE);
#ifdef WIN_MSWIN
	{ /* contre bug fenetre mange un peu du panel */
	Int2 right;
  	Advance(win);
  	GetNextPosition(win, &position);
  	right = position.x + 5;
  	SetNextPosition(win, position);
  	Break(win);
  	GetNextPosition(win, &position);
  	position.x = right;
  	position.y += 3;
  	SetNextPosition(win, position);
  	StaticPrompt(win, "", 1, 1, programFont, 'l');
  	}
#endif
  return fdui;
}


int update_tree_w_data(Nlm_Handle h)
{
FD_nj_plot *f;
Nlm_WindoW w;
#if defined(WIN_MAC) 
w = Nlm_FrontWindow();
#elif defined(WIN_MSWIN)
if(h == NULL) w = Nlm_FrontWindow();
else w = Nlm_ParentWindow(h);
#else
if(h == NULL) w = fd_nj_plot->nj_plot;
else w = Nlm_ParentWindow(h);
#endif
f = (FD_nj_plot *)Nlm_GetWindowExtra(w);
if(f== NULL || f->tag != tag_njplot) return FALSE;
fd_nj_plot = f;
return TRUE;
}

FILE *open_path(char *fname)  /* to open in read-only file fname searching for 
				 it through all path directories */
{
#define Mxdir 200
        static char dir[Mxdir+1];
        char *path, *deb, *fin;
        FILE *fich;
        int lf, ltot;
#ifdef __VMS
	static char vmspath[] = "njplot";
	path = vmspath;
#elif defined(unix)
        path = getenv("PATH"); 	/* get the list of path directories, 
					separated by :
    				*/
#elif defined(WIN_MSWIN)
/* try first dir where program was launched */
	path = get_prog_dir();
#else
	path = NULL;
#endif
        if (path == NULL ) {
        	fich = fopen(fname,"r");
        	return fich;
        	}
        lf=strlen(fname);
        deb=path;
        do
                {
#ifdef WIN_MSWIN
		fin = NULL;
#else
                fin = strchr(deb,':');
#endif
                if(fin!=NULL)
                        { strncpy(dir,deb,fin-deb); ltot=fin-deb; }
                else
                        { strcpy(dir,deb); ltot=strlen(dir); }
                /* now one directory is in string dir */
                if( ltot + lf + 1 <= Mxdir)
                        {
#ifdef __VMS
                        dir[ltot]=':';
#elif defined(WIN_MSWIN)
			if(dir[ltot-1] == '\\') ltot--;
                        else 	dir[ltot] = '\\';
#else 
                        dir[ltot]='/';
#endif
                        strcpy(dir+ltot+1,fname); /* now dir is appended with filename */
                        fich = fopen(dir,"r");
                        if( fich  != NULL) break;
                        }
                else fich = NULL;
                deb=fin+1;
                }
        while (fin != NULL);
	if(fich == NULL) {
			fich=fopen(fname,"r"); /* try also current directory */
			}
        return fich;
#undef Mxdir
}


void cre_help_win(IteM bouton)
{
static WindoW win;
LisT topic_list;
static struct win_help_extra win_data;
static int first = TRUE;
int larg;
GrouP group;

if(first) {
	if(help_file == NULL) {
		Disable(bouton);
		err_message(
#ifdef WIN_MAC
	"Sorry, no help resource available."
#else
	"Sorry, help file njplot.help is not found."
#endif
		);
		return;
		}
	first = FALSE;
	win = DocumentWindow( -50, -33, -3, -3, "Njplot Help", (Nlm_WndActnProc)help_ok_callback, NULL);
	group = HiddenGroup(win, -1, -2, NULL);
	SetGroupSpacing(group, 0, 5);
	topic_list = SingleList(group, 7, 8, topic_callback);
	load_help_topics(help_file, topic_list);
	DefaultButton(group, "Ok", help_ok_callback);
	Advance(win);
#ifdef WIN_MSWIN
	win_data.police = ParseFont("Courier New,10");
#else
	win_data.police = programFont;
#endif
	SelectFont(win_data.police);
	larg = 81 * Nlm_TextWidth("Q", 1);
	win_data.document = DocumentPanel(win, larg, 20 * Nlm_LineHeight() );
	win_data.tag = tag_help;
	SetWindowExtra(win, &win_data, NULL);
	}
Show(win);
Select(win);
}


void topic_callback(LisT list)
{
int num_topic;
char ligne[100], *p;
static int last_topic = 0;
Int2 total;
static char txt[5000];
WindoW win;
DoC help_doc;
struct win_help_extra *p_win_data;

num_topic = GetValue(list);
if(num_topic == 0 || num_topic == last_topic) return;
last_topic = num_topic;
win = ParentWindow(list);
p_win_data = (struct win_help_extra *)GetWindowExtra(win);
help_doc = p_win_data->document;
Hide(help_doc);
GetDocParams(help_doc, &total, NULL);
while(total > 0) DeleteItem(help_doc, total--);
rewind(help_file);
while(fgets(ligne, sizeof(ligne), help_file) != NULL) {
	if(strncmp(ligne, ">>>", 3) != 0) continue;
	if( --num_topic == 0) break;
	}
if(num_topic != 0) return;
p = txt;
while(fgets(ligne, sizeof(ligne), help_file) != NULL) {
	if(strncmp(ligne, ">>>", 3) == 0) break;
	if(strncmp(ligne, "Version:", 8) == 0) sprintf(ligne, "Version: %s\n", NJPLOTVERSION);
	strcpy(p, ligne); p += strlen(ligne);
	}
AppendText(help_doc, txt, NULL, NULL, p_win_data->police);
Show(help_doc);
}


void help_ok_callback(ButtoN bouton)
{
Nlm_WindoW h;

h = Nlm_ParentWindow(bouton);
Hide(h);
}

void load_help_topics(FILE *fich, LisT topic_list)
{
char ligne[100], topic[50];
int l;

rewind(fich);
while(fgets(ligne, sizeof(ligne), fich) != NULL) {
	if(strncmp(ligne, ">>>", 3) != 0) continue;
	strcpy(topic, ligne + 3);
	l = strlen(topic) - 1;
	topic[l] = 0;
	ListItem(topic_list, topic);
	}
}


void subtree_up_callback(ButtoN bouton)
{
fd_nj_plot = (FD_nj_plot *)Nlm_GetWindowExtra( Nlm_ParentWindow(bouton) );
if(fd_nj_plot->subtree_center == NULL) return;
fd_nj_plot->subtree_center = fd_nj_plot->subtree_ascend;
if(fd_nj_plot->subtree_center != NULL) {
	fd_nj_plot->subtree_ascend = fd_nj_plot->subtree_center->v3;
	if(fd_nj_plot->subtree_ascend == NULL) {
		SetValue(fd_nj_plot->choix_group, show_tree);
		Disable(fd_nj_plot->subtree_up_button);
		fd_nj_plot->choix = show_tree;
		show_tree_callback();
		return;
		}
	fd_nj_plot->subtree_notu = calc_n_desc(fd_nj_plot->subtree_center) - 1;
	scrollnotu(fd_nj_plot->subtree_notu);
	}
fd_nj_plot->need_runtree = TRUE;
tree_draw_proc(fd_nj_plot->tree_plot);
}


void process_keys(Nlm_Char key)
{
#ifdef WIN_MAC
	/* combinaison command-lettre */
#define test_key(a) (Nlm_cmmdKey && Nlm_currentKey == a)
#else
	/* combinaison control-lettre */
#define test_key(a)  key == a - 96
#endif
#ifdef MICRO
if( (!(test_key('o') || test_key('n'))) && ! update_tree_w_data(NULL) ) return;
#endif
if( test_key('q') ) 
	exit(0);
else if( test_key('o') ) {
	if(fd_nj_plot == NULL) fd_nj_plot = create_win_nj_plot();
	open_callback(fd_nj_plot->open_button);
	}
else if( test_key('n') )
	new_callback(NULL);
else if( test_key('w') )
	closefrommenu(fd_nj_plot->open_button);
else if( test_key('s') )
	save_plot_callback(fd_nj_plot->save_plot_button);
else if( test_key('f') )
	search_callback(fd_nj_plot->SearchItem);
else if( test_key('a') ) {
	if(Enabled(fd_nj_plot->AgainItem) ) 
		search_callback(fd_nj_plot->AgainItem);
	}
#ifdef MICRO
else if( test_key('c') )
	copy_plot(NULL);
else if( test_key('p') )
	print_plot(NULL);
#endif
else if( test_key('v') )
	paste_tree(NULL);
return;
}


static void closeproc(WindoW i)
{
Nlm_WindoW w;
FD_nj_plot *tofree;
fd_nj_plot = (FD_nj_plot *)Nlm_GetWindowExtra( i );
tofree = fd_nj_plot; /* necessary because with MSWIN fd_nj_plot changes after Remove(i) */
free_tree();
#ifdef WIN_MAC
Nlm_Disable(tofree->win_menu_item);
#endif
Nlm_Remove(i);
free(tofree);
w = Nlm_desktopWindow;
while(w != NULL) { /* is there another tree window ? */
	fd_nj_plot = (FD_nj_plot *)Nlm_GetWindowExtra( w );
	if( fd_nj_plot != NULL && fd_nj_plot->tag == tag_njplot) {
		Nlm_Select(fd_nj_plot->nj_plot);
		return;
		}
	w = Nlm_GetNext(w);
	}
#ifdef WIN_MAC
fd_nj_plot = NULL;
#else
exit(0);
#endif
}


void closefrommenu(Nlm_IteM item)
{
#ifdef WIN_MAC
if( update_tree_w_data(NULL) ) closeproc(fd_nj_plot->nj_plot);
#else
if( update_tree_w_data(item) ) closeproc(Nlm_ParentWindow(item));
#endif
}

void win_resize_proc(WindoW win)
{
RecT *old_rect, rect;
int change_x, change_y, old_height, new_height;

fd_nj_plot = (FD_nj_plot *)GetWindowExtra(win);
old_rect = &(fd_nj_plot->old_rect);
if(old_rect->left == -1) return;
ObjectRect(win, &rect);
old_height = old_rect->bottom - old_rect->top;
new_height = rect.bottom - rect.top;
change_y = (new_height - old_height);
old_height = old_rect->right - old_rect->left;
new_height = rect.right - rect.left;
change_x = (new_height - old_height);
if(change_x != 0 || change_y != 0) {
	change_panel_size(fd_nj_plot->tree_plot, change_x, change_y);
	*old_rect = rect;
	}
}


void change_panel_size(PaneL panel, int change_x, int change_y)
{
Nlm_RecT r;
#ifdef WIN_MOTIF
Nlm_RecT r_group;
GrouP group;
Nlm_Boolean in_group;
Nlm_BaR bar;
int visible;
static int recursive = FALSE;
#endif

GetPosition(panel, &r);
#ifdef WIN_MOTIF
if(recursive) return;
recursive = TRUE;
bar = Nlm_GetSlateVScrollBar( (Nlm_SlatE)panel);
visible = Nlm_Visible(bar);
group = Parent(panel);
in_group = (group != ParentWindow(panel));
if(in_group) {
	GetPosition(group, &r_group);
	Hide(group);
	}
else
	Hide(panel); 
#endif
r.bottom += change_y;
r.right += change_x;
SetPosition(panel, &r);
#ifdef WIN_MOTIF
if(in_group) {
	r_group.bottom += change_y;
	r_group.right += change_x;
	SetPosition( group, &r_group);
	Show(group);
	}
else	Show(panel);
if( !visible) Nlm_Hide(bar); 
recursive = FALSE;
#endif
}


void font_callback(ChoicE choix)
{
if( ! update_tree_w_data(choix) ) return;
prepare_fonts();
fd_nj_plot->need_runtree = TRUE;
tree_draw_proc(fd_nj_plot->tree_plot);
}


void prepare_fonts(void)
{
Nlm_Boolean use_bold, use_italic;
char font_full_name[40], *p;

/* prepare vibrant font */
use_bold = GetStatus(fd_nj_plot->bold_item);
use_italic = GetStatus(fd_nj_plot->italic_item);
fd_nj_plot->font_family_rank = GetValue(fd_nj_plot->choix_font);
fd_nj_plot->font_size_rank = GetValue(fd_nj_plot->choix_taille);
sprintf(font_full_name, "%s,%d", list_font_name[fd_nj_plot->font_family_rank], 
	list_font_size[fd_nj_plot->font_size_rank]);
if(use_bold || use_italic) strcat(font_full_name, ",");
if(use_bold) strcat(font_full_name, "b");
if(use_italic) strcat(font_full_name, "i");
fd_nj_plot->current_font = ParseFont(font_full_name);
if((p = strchr(font_full_name, ',')) != NULL) strcpy(fd_nj_plot->font_bold_italic, p + 1);
else fd_nj_plot->font_bold_italic[0] = 0;
prepare_pdf_font(fd_nj_plot->font_family_rank, use_bold, use_italic);
}


/* callbacks for form nj_plot */
#if defined(WIN_MOTIF) && defined(__APPLE__)
/* bug sur mac/darwin : Nlm_GetOutputFileName crashes if Nlm_GetInputFileName was not used before ! */
static int can_run_getoutputfname = FALSE;
#else
static int can_run_getoutputfname = TRUE;
#endif


void open_callback(IteM ob)
{
char fname[200];

can_run_getoutputfname = TRUE;
if(!MG_GetInputFName(fname, sizeof(fname))) {
	if(fd_nj_plot != NULL) tree_draw_proc(fd_nj_plot->tree_plot);
	return;
	}
if(fd_nj_plot == NULL || fd_nj_plot->notu != 0) fd_nj_plot = create_win_nj_plot();
fd_nj_plot->notu = 0;
init_tree( fname, NULL );
if(fd_nj_plot->notu == 0) {
	closeproc(fd_nj_plot->nj_plot);
	return;
	}
Nlm_Show(fd_nj_plot->nj_plot);
Nlm_ObjectRect(fd_nj_plot->nj_plot, &(fd_nj_plot->old_rect));
}


void new_callback(IteM ob)
{
fd_nj_plot = create_win_nj_plot();
fd_nj_plot->notu = 0;
Nlm_Show(fd_nj_plot->nj_plot);
Nlm_ObjectRect(fd_nj_plot->nj_plot, &(fd_nj_plot->old_rect));
}


void paper_callback(ChoicE ob)
{
paper_choice = (paper_item) GetValue(ob);
}


void save_plot_callback(IteM ob)
{
char *p;
if(!update_tree_w_data(ob)) return;
strcpy(plotfilename, fd_nj_plot->tree_name);
p = strchr(extract_filename(plotfilename), '.'); 
if(p == NULL) p = plotfilename + strlen(plotfilename);
#ifdef WIN_MAC
strcpy(p,".pdf");
#else
#ifdef WITH_PDF
pdf_plot = (ob == fd_nj_plot->save_plot_button);
if(pdf_plot) 
	strcpy(p,".pdf"); 
else 
#endif
strcpy(p,".ps");
#endif
if(ob != NULL) {
if(can_run_getoutputfname) {
	if( ! MG_GetOutputFName(plotfilename, sizeof(plotfilename), 
		extract_filename(plotfilename) ) ) return;
	}
	}

#ifdef WIN_MAC
pict_plot();
#else
#ifdef WITH_PDF
if(pdf_plot)
	plot_to_pdf();
else
#endif
    postscript_plot(); 
#endif
}


void save_tree_callback(IteM ob)
{
char *p, *arbre, ajout[10], *name_part;
FILE *out;
int unrooted;

if(!update_tree_w_data(ob)) return;
if(fd_nj_plot->notu==0)return;
unrooted = (ob == fd_nj_plot->save_unrooted_button);
/* build out file name */
strcpy(rooted_fname, fd_nj_plot->tree_name);
name_part = extract_filename(rooted_fname);
/* ajouter _root. dans le nom */
strcpy(ajout, unrooted ? "_noroot." : "_root.");
if(strstr(name_part, ajout) == NULL) {
	p = strchr(name_part, '.'); 
	if(p == NULL) p = name_part + strlen(name_part);
	strcpy(p, ajout);
	}
name_part = extract_filename(fd_nj_plot->tree_name);
p = strchr(name_part, '.');
if(p == NULL)
	strcat(rooted_fname, "ph");
else
	strcpy( strchr(extract_filename(rooted_fname), '.'), p);
if( !MG_GetOutputFName(rooted_fname, sizeof(rooted_fname), extract_filename(rooted_fname) ) ) 
		return;
if(unrooted) arbre = ecrit_arbre_parenth_unrooted(fd_nj_plot->racine);
else arbre = ecrit_arbre_parenth(fd_nj_plot->racine);
if(arbre == NULL) {
	err_message("Sorry, not enough memory");
	return;
	}
out = fopen(rooted_fname,"w");
if(out != NULL) {
	fputs(arbre, out); putc('\n', out);
	fclose(out);
	}
else	{
	char mess[250];
	sprintf(mess, "Sorry, cannot write to file %s", 
		extract_filename(rooted_fname) );
	err_message(mess);
	}
free(arbre);
}


void toggle_branch_callback(ButtoN ob)
{
fd_nj_plot = (FD_nj_plot *)Nlm_GetWindowExtra( Nlm_ParentWindow(ob) );
	fd_nj_plot->plot_br_l = GetStatus(ob);
	fd_nj_plot->need_runtree = TRUE;
	tree_draw_proc(fd_nj_plot->tree_plot);
}


void toggle_bootstrap_callback(ButtoN ob)
{
fd_nj_plot = (FD_nj_plot *)Nlm_GetWindowExtra( Nlm_ParentWindow(ob) );
fd_nj_plot->show_bootstrap = !fd_nj_plot->show_bootstrap;
fd_nj_plot->need_runtree = TRUE;
tree_draw_proc(fd_nj_plot->tree_plot);
}


void operation_callback(GrouP group)
{
fd_nj_plot = (FD_nj_plot *)Nlm_GetWindowExtra( Nlm_ParentWindow(group) );
fd_nj_plot->choix = (enum operations) GetValue(group);
if(fd_nj_plot->choix == permutation)
	swap_callback();
else if(fd_nj_plot->choix == depl_racine)
	new_outgroup_callback();
else if(fd_nj_plot->choix == subtree)
	subtree_callback();
else {
	show_tree_callback();
	}
}

void subtree_callback(void)
{
	fd_nj_plot->choix = subtree;
/*	Disable(fd_nj_plot->save_plot_button); */
	Disable(fd_nj_plot->save_tree_button);
	Disable(fd_nj_plot->save_unrooted_button);
	Enable(fd_nj_plot->subtree_up_button);
	if(fd_nj_plot->has_br_length) {
		SetStatus(fd_nj_plot->branch_length_button, FALSE);
		Disable(fd_nj_plot->branch_length_button);
		fd_nj_plot->plot_br_l=0;
		}
	SetStatus(fd_nj_plot->bootstrap_button, FALSE);
	fd_nj_plot->show_bootstrap = FALSE;
	Disable(fd_nj_plot->bootstrap_button);
	Disable(fd_nj_plot->reroot_button);
fd_nj_plot->need_runtree = TRUE;
tree_draw_proc(fd_nj_plot->tree_plot);
}

void swap_callback(void)
{
	fd_nj_plot->choix = permutation;
	Disable(fd_nj_plot->save_plot_button);
	Disable(fd_nj_plot->save_tree_button);
	Disable(fd_nj_plot->save_unrooted_button);
	if(fd_nj_plot->has_br_length) {
		SetStatus(fd_nj_plot->branch_length_button, FALSE);
		Disable(fd_nj_plot->branch_length_button);
		fd_nj_plot->plot_br_l=0;
		}
	SetStatus(fd_nj_plot->bootstrap_button, FALSE);
	fd_nj_plot->show_bootstrap = FALSE;
	Disable(fd_nj_plot->bootstrap_button);
fd_nj_plot->need_runtree = TRUE;
tree_draw_proc(fd_nj_plot->tree_plot);
}

void new_outgroup_callback(void)
{
	fd_nj_plot->choix = depl_racine;
	Disable(fd_nj_plot->save_plot_button);
	Disable(fd_nj_plot->save_tree_button);
	Disable(fd_nj_plot->save_unrooted_button);
	if(fd_nj_plot->has_br_length) {
		SetStatus(fd_nj_plot->branch_length_button, FALSE);
		Disable(fd_nj_plot->branch_length_button);
		fd_nj_plot->plot_br_l=0;
		}
	SetStatus(fd_nj_plot->bootstrap_button, FALSE);
	fd_nj_plot->show_bootstrap = FALSE;
	Disable(fd_nj_plot->bootstrap_button);
fd_nj_plot->need_runtree = TRUE;
tree_draw_proc(fd_nj_plot->tree_plot);
}


void show_tree_callback(void)
{
	fd_nj_plot->choix = show_tree;
	Enable(fd_nj_plot->save_plot_button);
	Enable(fd_nj_plot->save_tree_button);
	Enable(fd_nj_plot->save_unrooted_button);
	Disable(fd_nj_plot->subtree_up_button);
	if(fd_nj_plot->has_br_length) {
		Enable(fd_nj_plot->branch_length_button);
		}
	if(fd_nj_plot->has_internal) {
		Enable(fd_nj_plot->bootstrap_button);
		}
	Enable(fd_nj_plot->reroot_button);
fd_nj_plot->need_runtree = TRUE;
fd_nj_plot->subtree_notu = fd_nj_plot->notu;
if(fd_nj_plot->subtree_center != NULL) scrollnotu(fd_nj_plot->notu);
fd_nj_plot->subtree_center = NULL;
tree_draw_proc(fd_nj_plot->tree_plot);
}


void exit_callback(IteM ob)
{
exit(0);
}


void tree_draw_proc(PaneL panel)
{
RecT ob_rect;
Int2 width, height;
static Int2 previous_h = 0, previous_w = 0;
Nlm_BaR bar;
int hidebar;

fd_nj_plot = (FD_nj_plot *)Nlm_GetWindowExtra( Nlm_ParentWindow(panel) );
WatchCursor();
bar = Nlm_GetSlateVScrollBar((Nlm_SlatE)panel);
Nlm_GetPanelExtra(panel, &hidebar);
if(hidebar && Nlm_Visible(bar)) Nlm_Hide(bar);
else if( !hidebar && !Nlm_Visible(bar)) Nlm_Show(bar);
Select (panel);
ObjectRect(panel, &ob_rect);
EraseRect(&ob_rect);
FrameRect(&ob_rect);
Nlm_ClipRect(&ob_rect);
	if(fd_nj_plot->notu == 0) {
		ArrowCursor();
		return;
		}
	width = ob_rect.right - ob_rect.left;
	height = ob_rect.bottom - ob_rect.top;
	SelectFont(fd_nj_plot->current_font);
	Black();
	if( fd_nj_plot->need_runtree || /* recalcul des positions des # si window re-size */
		(fd_nj_plot->choix != show_tree && (height != previous_h || width != previous_w))) {
		runtree();
			}
	fd_nj_plot->need_runtree = FALSE;
	swap = 0;
	do_plot();
	previous_w = width;
	previous_h = height;
Nlm_ResetClip();
ArrowCursor();
}


void tree_click_proc(PaneL panel, PoinT click)
{
	double eps, x, y;
	int i, found=0, node_num;
	RecT ob_rect;
	Int2 panel_x, panel_y, panel_h;
	fd_nj_plot = (FD_nj_plot *)Nlm_GetWindowExtra( Nlm_ParentWindow(panel) );
	if(fd_nj_plot->notu == 0) return;
	if(fd_nj_plot->choix != permutation && fd_nj_plot->choix != depl_racine && fd_nj_plot->choix != subtree) 
		return;
	ObjectRect(panel, &ob_rect);
	panel_x = ob_rect.left; panel_y = ob_rect.top; panel_h = ob_rect.bottom - ob_rect.top;
	eps = fd_nj_plot->char_height/2;
	x = click.x - panel_x;
	click.y -= panel_y;
	y = panel_h - click.y;
	for(i=0; i<=fd_nj_plot->totpoints; i++ ) {
		if(fabs(x - (fd_nj_plot->points+i)->x) > eps)
			 continue;
		if(fabs(y - (fd_nj_plot->points+i)->y) > eps)
			 continue;
		found=1;
		break;
		}
	if(!found)return;
	node_num=(fd_nj_plot->points+i)->number;
	if(fd_nj_plot->choix==depl_racine) {
		if(node_num >= 1 && node_num <= 2*fd_nj_plot->notu && 
						node_num-1 != fd_nj_plot->root_num) {
			fd_nj_plot->root_num = node_num-1;
			removeroot();
			}
		}
	else if(fd_nj_plot->choix==permutation) {
		if(node_num >= fd_nj_plot->notu+2 && node_num <= 2*fd_nj_plot->notu+1) 
			swap= node_num-1;
		else
			swap=0;
		}
	else if(fd_nj_plot->choix==subtree) {
		if(node_num < fd_nj_plot->notu+2 || node_num > 2*fd_nj_plot->notu+1) return;
		fd_nj_plot->subtree_center = fd_nj_plot->tabtax[node_num - 1];
		fd_nj_plot->subtree_ascend = fd_nj_plot->subtree_center->v3;
		fd_nj_plot->subtree_notu = calc_n_desc(fd_nj_plot->subtree_center) - 1;
		scrollnotu(fd_nj_plot->subtree_notu);
		fd_nj_plot->choix = show_tree;
		SetValue(fd_nj_plot->choix_group, 0);
		Nlm_InvalObject(fd_nj_plot->choix_group);
		if(fd_nj_plot->has_br_length) {
			Enable(fd_nj_plot->branch_length_button);
			Nlm_InvalObject(fd_nj_plot->branch_length_button);
			}
		if(fd_nj_plot->has_internal) {
			Enable(fd_nj_plot->bootstrap_button);
			Nlm_InvalObject(fd_nj_plot->bootstrap_button);
			}
		}
	fd_nj_plot->need_runtree = TRUE;
	tree_draw_proc(panel);
}


void paste_tree(IteM ob)
{
FILE *tmp;
char *buff, *tmpfname;
if(!update_tree_w_data(ob)) return;
if(fd_nj_plot->notu != 0) {
	err_message("Do \"Edit:clear\" before pasting tree data");
	return;
	}

#ifdef WIN_MAC
Nlm_textScrapFull = FALSE; /* arrange un bug de Nlm_ClipboardHasString */
#endif
if( ! Nlm_ClipboardHasString() ) {
	err_message("No text describing tree present in clipboard");
	return;
	}
buff = Nlm_ClipboardToString();
tmpfname = tempnam(NULL, "njplottemp_"); 
tmp=fopen(tmpfname,"w");
fwrite(buff, 1, strlen(buff), tmp);
fclose(tmp);
init_tree(tmpfname, "pasted tree");
remove(tmpfname);
}



void clear_tree(IteM i)
{
if(!update_tree_w_data(i)) return;
if(fd_nj_plot->notu != 0) free_tree();
fd_nj_plot->notu = 0;
SetTitle(fd_nj_plot->nj_plot, "njplot");
#ifdef WIN_MAC
SetTitle(fd_nj_plot->win_menu_item, "njplot");
#endif
tree_draw_proc(fd_nj_plot->tree_plot);
}


void change_page_count(ButtoN item)
{
Nlm_Boolean *pdone;

pdone = (Nlm_Boolean *)GetWindowExtra( Nlm_ParentWindow(item) );
*pdone = TRUE;
}


void page_count_callback(IteM item)
{
static WindoW win;
static Nlm_Boolean done;
static TexT select_box;
char select[100];
static int first = TRUE;
	if(first) {
		win = FixedWindow(-50,-50,-5,-5,"Page count:",NULL);
		sprintf(select, "%d", page_count);
		select_box=DialogText(win,select,8,NULL);
		Advance(win);
		DefaultButton(win,"ok",change_page_count);
		SetWindowExtra(win, &done, NULL);
		first = FALSE;
		}
	Show(win);
	Select(win); Select(select_box);
	done=FALSE;
	while(! done) {
		myWaitAndProcessNextEvent();
		}
	GetTitle(select_box, select, sizeof(select));
	sscanf(select, "%d", &page_count);
	Hide(win);
#ifndef WIN_MSWIN /* changement titre menu item impossible sur PC ! */
	sprintf(select,"Page count (%d)", page_count);
	Nlm_SetTitle(item, select);
#endif	
	Select(fd_nj_plot->nj_plot);
}




void string_callback(ButtoN item)
{
Nlm_Boolean *pdone;

pdone = (Nlm_Boolean *)GetWindowExtra( Nlm_ParentWindow(item) );
*pdone = TRUE;
}


void search_callback(IteM item)
{
static WindoW win;
static Nlm_Boolean search_done;
static TexT select_box;
static char select[500];
char aux[500];
int num, trouve;
static int first = TRUE;

if(!update_tree_w_data(item)) return;
if(item == fd_nj_plot->SearchItem) {
	if(first) {
		win = FixedWindow(-50,-50,-5,-5,"Name searched:",NULL);
		select_box=DialogText(win,"",15,NULL);
		Advance(win);
		DefaultButton(win,"ok",string_callback);
		SetWindowExtra(win, &search_done, NULL);
		first = FALSE;
		}
	Show(win);
	Select(win); Select(select_box);
	search_done=FALSE;
	while(! search_done) {
		myWaitAndProcessNextEvent();
		}
	GetTitle(select_box, select, sizeof(select));
	majuscules(select);
	Hide(win);
/*	Enable(fd_nj_plot->AgainItem); */
	Select(fd_nj_plot->nj_plot);
	}
else if(item != fd_nj_plot->AgainItem) { /* happens only when called after start  */
	strcpy(select, (char *)item);
	majuscules(select);
	}
if(strlen(select) == 0) return;
trouve = FALSE;
for(num = 0; num <= fd_nj_plot->totnoms; num++) {
	strcpy(aux, (fd_nj_plot->noms+num)->nom);
	majuscules(aux);
	if(strstr( aux, select) != NULL) {
		(fd_nj_plot->noms+num)->disp_option = 'r';
		trouve = TRUE;
		}
	}
if(trouve) {
	tree_draw_proc(fd_nj_plot->tree_plot);
	}
}


void init_tree(char *fname, char *displayname)
{
char *pname, *p;
char titre_l[300];

if(displayname == NULL) displayname = extract_filename(fname);
fd_nj_plot->choix = show_tree;
WatchCursor();
/* read tree file */
if( (pname=preptree(fname)) != NULL ) {
	char *message;
	size_t taille;
	fd_nj_plot->tree_name = NULL;
	fd_nj_plot->notu = 0;
	SetTitle(fd_nj_plot->nj_plot, "njplot");
	tree_draw_proc(fd_nj_plot->tree_plot);
	taille = strlen(pname) + strlen(fname) + 3;
	message = check_alloc(1, taille);
	strcpy(message, pname);
	strcat(message, " ");
	strcat(message, extract_filename(fname) );
	err_message(message);
	free(message);
	return;
	}
maxy = 1000.;
prepare_fonts();
fd_nj_plot->deltay = maxy/fd_nj_plot->notu;
if(fd_nj_plot->notu > 0) fd_nj_plot->tree_name = strdup(fname);
fd_nj_plot->need_runtree = TRUE;
fd_nj_plot->subtree_center = NULL;

strcpy(titre_l, displayname );
if(fd_nj_plot->notu > 0) sprintf(titre_l + strlen(titre_l), " (%d tips)", fd_nj_plot->notu + 1);
p = titre_l;
#ifdef WIN_MAC
p = mac_fname_to_roman(titre_l); 
SetTitle(fd_nj_plot->win_menu_item, displayname);
#endif
SetTitle(fd_nj_plot->nj_plot, p);

tree_draw_proc(fd_nj_plot->tree_plot);
if(fd_nj_plot->has_br_length) {
	Enable(fd_nj_plot->branch_length_button);
	}
else
	Disable(fd_nj_plot->branch_length_button);
if(fd_nj_plot->has_internal) {
	Enable(fd_nj_plot->bootstrap_button);
	}
else
	Disable(fd_nj_plot->bootstrap_button);
}


/* fenetre message d'alerte.
Operation Select dans la boucle des evenements qui ramene toujours la fenetre
au premier plan. Lent, mais interessant.
*/
static void alert_ok_action(ButtoN bouton)
{
Nlm_Boolean *alert_done;
alert_done = (Nlm_Boolean *)GetWindowExtra(Parent(bouton));
*alert_done = TRUE;
}


void err_message(const char *in)
{
Nlm_WindoW alwin;
Nlm_Boolean alert_done;
Nlm_TexT t;
int nl;
char *p, *q, *r;
static char texte[1000];
int w, width;

#ifdef WIN_MAC
alwin = Nlm_FixedWindow(-50,-50,-5,-5,"MESSAGE",NULL); // Bug: ModalWindow crashes when return key hit
#else
alwin = Nlm_ModalWindow(-50,-50,-5,-5,NULL);
Nlm_SetTitle(alwin, "Message");
#endif
Nlm_SelectFont (Nlm_systemFont);
nl = 0; p = (char *)in; r = texte; width = 0;
while(TRUE)	{
	q = strchr(p,'\n');
	if(q == NULL) q = p + strlen(p);
	nl++; 
	if( (r - texte) + (q - p) + 2 >= sizeof(texte) ) break;
	if(q > p) {
		memcpy(r, p, q - p); 
		w = Nlm_TextWidth(p, q - p);
		if(w > width) width = w;
		}
	r += q - p;
#ifdef WIN_MSWIN
	*(r++) = '\r';
#endif
#ifdef WIN_MAC
	*(r++) = '\r';
#else
	*(r++) = '\n';
#endif
	*r = 0;
	q = strchr(p,'\n');
	if(q == NULL) break;
	p = q + 1; 
	}

width = ( width + 20 ) / Nlm_stdCharWidth;
t = Nlm_ScrollText(alwin, width, nl , Nlm_systemFont, TRUE, NULL);

Nlm_SetTitle(t, texte);
Nlm_Break(alwin);
Nlm_DefaultButton(alwin,"ok",alert_ok_action);
Nlm_Show(alwin);
alert_done=FALSE;
Nlm_SetWindowExtra(alwin, &alert_done, NULL);
Nlm_ArrowCursor();
while(! alert_done) {
#ifdef WIN_MAC
	if(! Nlm_InFront(alwin) ) 
		Nlm_Select(alwin); // to emulate modal window on mac
#endif
	myWaitAndProcessNextEvent();
	}
Nlm_Remove(alwin);
}

#else /* NO_GUI */

void err_message(const char *in)
{
fprintf(stderr, "%s\n", in);
}

#endif /* NO_GUI */


void prepare_pdf_font(int font_num, int use_bold, int use_italic)
{
/* prepare PDF name of font */
strcpy(current_ps_font, list_ps_font_name[font_num]);
if( use_bold && use_italic )
	if(font_num == times) 
		strcat(current_ps_font,"-BoldItalic");
	else
		strcat(current_ps_font,"-BoldOblique");
else if( use_bold )
	strcat(current_ps_font,"-Bold");
else if( use_italic )  {
	if(font_num == times) 
		strcat(current_ps_font,"-Italic");
	else
		strcat(current_ps_font,"-Oblique");
	}
else if(font_num == times) 
	strcat(current_ps_font,"-Roman");
}


int direct_pdf_plot(char *fname)
{
char *p, *q;

fd_nj_plot->choix = show_tree;
/* read tree file */
if( (p = preptree(fname) ) != NULL ) {
	fprintf(stderr, "%s\n", p);
	return 1;
	}
strcpy(plotfilename, fname);
q = extract_filename(plotfilename);
p = strchr(q, '.'); 
if(p == NULL) p = plotfilename + strlen(plotfilename);
file_plot = TRUE;
#ifdef TTY
strcpy(p, ".txt");
#else
#ifdef WITH_PDF
strcpy(p, ".pdf" );
pdf_plot = TRUE;
#else
strcpy(p, ".ps");
#endif
prepare_pdf_font(times, FALSE, FALSE);
#endif
maxy = 1000.;
fd_nj_plot->deltay = maxy/fd_nj_plot->notu;
runtree();
#ifdef ADDROOT
/* calcul du nombre de feuilles (tres obligatoire !!!) */
calc_otu(NULL, fd_nj_plot->racine);
/*  calcul de son format phylip */
p = ecrit_arbre_parenth(fd_nj_plot->racine);
if(p == NULL) {
	fprintf(stderr, "Sorry, not enough memory\n");
	exit(1);
	}
/* ecriture de son format phylip */
printf("%s\n", p);
#elif defined(TTY)
tty_plot();
#elif !defined(WIN_MAC)
#ifdef WITH_PDF
plot_to_pdf();
#else
postscript_plot();
#endif
#endif
return 0;
}

int calc_otu(struct noeud *pere, struct noeud *centre)
{
struct noeud *gauche, *droite;
double bg, bd, bpere;

if(centre == NULL) return 0;

/* orienter le noeud centre de maniere standard: centre->v3=pere */
if( centre->v1 == pere ) {
	gauche =centre->v2; droite = centre->v3;
	bg = centre->l2; bd = centre->l3; bpere = centre->l1;
	}
else if( centre->v2 == pere ) {
	gauche =centre->v1; droite = centre->v3;
	bg = centre->l1; bd = centre->l3; bpere = centre->l2;
	}
else	{
	gauche =centre->v1; droite = centre->v2;
	bg = centre->l1; bd = centre->l2; bpere = centre->l3;
	}
centre->v3=pere; centre->v1=gauche; centre->v2=droite;
centre->l3=bpere; centre->l1=bg; centre->l2=bd;

return calc_otu(centre, centre->v1) + calc_otu(centre, centre->v2);
}


void remove_arg(int target, int *argc, char *argv[])
{
int num;
for(num = target; num < *argc - 1; num++) 
	argv[num] = argv[num+1];
(*argc)--;
}

#ifdef NO_PDF
#define PDFONLY "-psonly"
#else
#define PDFONLY "-pdfonly"
#endif

void process_args(int *argc, char *argv[], char **pred, int *pboot_arg, int *p_plot_br_l, int *p_font_size_rank)
{
int num, taille;

for(num = 1; num < *argc; num++) {
	if( strncmp(argv[num], "-h", 2) == 0 ){
		char message[1000];
		sprintf(message,
		"Usage: %s [-options] [tree_file_name]\n"
"where options are:\n"
"-h             print out this message\n"
#ifdef ADDROOT
"-outotu name   name of species to use as outgroup for rooting\n"
#elif defined(TTY)
"-txtw n        width of output text (default = 120)\n"
"-prefix string string to prefix to all leaf names in output tree\n"
#else
"-us            PDF or PostScript tree file prepared for US Letter paper size\n"
#ifndef NO_GUI
PDFONLY"       no window interface, just write the PDF/PostScript tree plot\n"
"               to file named as input tree file with .pdf/.ps extension\n"
"-red xxxx      highlight in red all taxon names containing string xxxx\n"
#endif
"-pc n          number of pages for PDF/PostScript output\n"
"-size n        font size n used for taxon names\n"
"-lengths       show branch lengths if they appear in tree file\n"
"-boot          show bootstrap values if they appear in tree file\n"
"-psize         size of page for PDF/PostScript expressed as WIDTHxHEIGHT\n"
"-notitle       don't include title in PDF/PostScript output\n"
#endif
"\n"
"and where tree_file_name is the name of a tree file in the Newick format\n" 
#ifdef TTY
"the output goes to file <tree_file_name>.txt\n"
#elif defined(ADDROOT)
"the output goes to stdout\n"
#elif defined(NO_GUI)
"the output goes to file <tree_file_name>.pdf or <tree_file_name>.ps\n"
#endif
, extract_filename(argv[0]) );

#ifdef MICRO
		err_message(
#else
		fprintf(stderr,
#endif
		message);
		exit(0);
		}
	}

*pboot_arg = FALSE; *p_plot_br_l = FALSE;
for(num = 1; num < *argc; num++) {
	if( strncmp(argv[num], "-us", 3) == 0) {
		paper_choice = LETTER;
		remove_arg(num, argc, argv);
		break;
		}
	}

for(num = 1; num < *argc; num++) {
	if( strncmp(argv[num], "-size", 5) == 0) {
		const int maxi = sizeof(list_font_size) / sizeof(int);
		taille = 12;
		if(num + 1 < *argc) {
			sscanf(argv[num + 1], "%d", &taille);
			remove_arg(num + 1, argc, argv);
			}
		remove_arg(num, argc, argv);
		for(num = 0; num < maxi; num++) 
			if(taille <= list_font_size[num]) break;
		if(num >= maxi) 
			num = maxi - 1;
		*p_font_size_rank = num;
		break;
		}
	}
for(num = 1; num < *argc; num++) {
	if( strncmp(argv[num], "-pc", 3) == 0) {
		if(num + 1 < *argc) {
			sscanf(argv[num + 1], "%d", &page_count);
			remove_arg(num + 1, argc, argv);
			}
		remove_arg(num, argc, argv);
		break;
		}
	}
for(num = 1; num < *argc; num++) {
	if( strncmp(argv[num], PDFONLY, 9) == 0) {
		pdf_plot_only = TRUE;
		remove_arg(num, argc, argv);
		break;
		}
	}
for(num = 1; num < *argc; num++) {
	if( strncmp(argv[num], "-lengths", 8) == 0) {
		*p_plot_br_l = TRUE;
		remove_arg(num, argc, argv);
		break;
		}
	}
for(num = 1; num < *argc; num++) {
	if( strncmp(argv[num], "-boot", 5) == 0) {
		*pboot_arg = TRUE;
		remove_arg(num, argc, argv);
		break;
		}
	}
for(num = 1; num < *argc; num++) {
	if( strncmp(argv[num], "-psize", 6) == 0) {
		if(num + 1 < *argc) {
			sscanf(argv[num + 1], "%dx%d", &ps_width, &ps_height);
			remove_arg(num + 1, argc, argv);
			}
		remove_arg(num, argc, argv);
		break;
		}
	}
for(num = 1; num < *argc; num++) {
	if( strncmp(argv[num], "-notitle", 8) == 0) {
		no_title = TRUE;
		remove_arg(num, argc, argv);
		break;
		}
	}
*pred = NULL;
#ifndef NO_GUI
for(num = 1; num < *argc; num++) {
	if( strncmp(argv[num], "-red", 4) == 0) {
		if(num + 1 < *argc) {
			*pred = argv[num+1];
			remove_arg(num + 1, argc, argv);
			}
		remove_arg(num, argc, argv);
		break;
		}
	}
#elif defined(TTY)
physx = 120; /* default width of text output */
for(num = 1; num < *argc; num++) {
	if( strncmp(argv[num], "-txtw", 5) == 0) {
		if(num + 1 < *argc) {
			int val;
			sscanf(argv[num+1], "%d", &val);
			physx = (double)val;
			remove_arg(num + 1, argc, argv);
			}
		remove_arg(num, argc, argv);
		break;
		}
	}
for(num = 1; num < *argc; num++) {
	if( strncmp(argv[num], "-prefix", 7) == 0) {
		if(num + 1 < *argc) {
			tty_prefix = argv[num+1];
			remove_arg(num + 1, argc, argv);
			}
		remove_arg(num, argc, argv);
		break;
		}
	}
#elif defined(ADDROOT)
for(num = 1; num < *argc; num++) {
	if(strcmp(argv[num], "-outotu") == 0) {
		if(num + 1 < *argc) {
			outotu = argv[num+1];
			remove_arg(num + 1, argc, argv);
			}
		remove_arg(num, argc, argv);
		break;
		}
	}
#endif
}


double length_log_phys(double p)
{
return p * tek_dx;
}

double length_phys_log(double p)
{
return  p / tek_dx;
}


double arrondi_echelle(double x)
{ /* arrondi x a une valeur 1, 2, 5 pour echelle */
double l, n;
int r;
static int corresp[] = {1,1,2,2,5,5,5,10,10,10,10,10};
                     /* 0,1,2,3,4,5,6, 7, 8, 9,10,11 */
l = log10(x);
n = floor(l);
l = x * pow(10., -n);  /* 10. plutot que 10 necessaire pour alpha! */
r = myrint(l); r = corresp[r];
return r * pow(10., n);  /* 10. plutot que 10 necessaire pour alpha! */
}


void scale_window(double lxmin, double lxmax, double lymin, double lymax,
 double pxmin, double pxmax, double pymin, double pymax)
/* to scale the plot window for logical coords l... in physical coords p... */
{
tek_xmin=lxmin; tek_xmax=lxmax; tek_ymin=lymin; tek_ymax=lymax;
tek_minx=pxmin; tek_maxx=pxmax; tek_miny=pymin; tek_maxy=pymax;
tek_dx = (tek_maxx-tek_minx)/(tek_xmax-tek_xmin < 1e-4 ? 1e-4 : tek_xmax-tek_xmin);
tek_dy = (tek_maxy-tek_miny)/(tek_ymax-tek_ymin);
}


void ch_echelle(double lx, double ly, double *px, double *py)
/* conversion coord logique lx,ly en coord physique px,py */
{
*px = (lx-tek_xmin)*tek_dx + tek_minx;
*py = (ly-tek_ymin)*tek_dy + tek_miny;
}


void convert_mem_point(void)
{
int i;
if(converted) return;
for (i=0; i<=fd_nj_plot->totpoints; i++) {
   	ch_echelle( (fd_nj_plot->points+i)->x, (fd_nj_plot->points+i)->y, &((fd_nj_plot->points+i)->x), 
			&((fd_nj_plot->points+i)->y) );
	(fd_nj_plot->points+i)->x += fd_nj_plot->char_height / 6.;
	}
converted = TRUE;
}


double calc_echelle(double larg)
{ /* rend taille logique pour echelle optimale */
double log_val, phys_val;
phys_val = larg/10;
log_val = length_phys_log(phys_val);
log_val = arrondi_echelle(log_val);
return log_val;
}


void draw_scale(void)
{
char ech_name[20];
double phys_w, y, xd, xf, lc;
double log_val;

log_val = calc_echelle(physx);
phys_w = myrint(length_log_phys(log_val));
y = physy - 1.35 * fd_nj_plot->char_height;
xf = physx * 0.95;
xd = xf - phys_w;
sprintf(ech_name, "%.1g", log_val);
lc = calc_text_size(ech_name, NULL, NULL);
dir_moveto(xd, y);
dir_lineto(xf, y);
#ifndef TTY
dir_moveto(xd, y - fd_nj_plot->char_height/3.);
dir_lineto(xd, y + fd_nj_plot->char_height/3.);
dir_moveto(xf, y - fd_nj_plot->char_height/3.);
dir_lineto(xf, y + fd_nj_plot->char_height/3.);
#endif
dir_moveto( (xd + xf)/2 - lc/2, 
	y + fd_nj_plot->char_height/6. ); 
plotstring(ech_name);
}


#ifndef NO_GUI
void scrollcallback(Nlm_BaR b, Nlm_SlatE g, Nlm_Int2 after, Nlm_Int2 before)
{
fd_nj_plot = (FD_nj_plot *)Nlm_GetWindowExtra( Nlm_ParentWindow(b) );
tree_draw_proc(fd_nj_plot->tree_plot);
}


void zoomcallback(Nlm_SwitcH s, Nlm_Int2 after, Nlm_Int2 before)
{
char aux[50];
int hide;

fd_nj_plot = (FD_nj_plot *)Nlm_GetWindowExtra( Nlm_ParentWindow(s) );
if(after == 1) fd_nj_plot->zoomvalue = 1.;
else if(after > before) fd_nj_plot->zoomvalue *= 1.4;
else if(after < before) fd_nj_plot->zoomvalue /= 1.4;
sprintf(aux, "%d%%", (int)(100. * fd_nj_plot->zoomvalue + 0.5) );
Nlm_SetTitle(fd_nj_plot->zoomprompt, aux);
Nlm_InvalObject(fd_nj_plot->zoomprompt);
scrollnotu(fd_nj_plot->subtree_notu);
hide = (after == 1);
Nlm_SetPanelExtra(fd_nj_plot->tree_plot, &hide);
tree_draw_proc(fd_nj_plot->tree_plot);
}


void scrollnotu(int count)
{
Nlm_BaR obj; int lines;
float old;

	obj = Nlm_GetSlateVScrollBar( (SlatE)fd_nj_plot->tree_plot );
	old = Nlm_GetBarValue(obj) / (float)Nlm_GetBarMax(obj);
	lines = (count + 1) * (1 - 1/fd_nj_plot->zoomvalue) + 0.5; 
	if(lines < 3) lines = 3;
	Nlm_CorrectBarMax(obj, lines);
	Nlm_CorrectBarValue(obj, old * lines);
	Nlm_CorrectBarPage(obj, lines/fd_nj_plot->zoomvalue, lines/fd_nj_plot->zoomvalue);
}

#endif

void do_plot(void)
{
struct trait *p;
int num;
double factor, value, m = 1, v = 0, zoom = 1;

#ifndef NO_GUI
/* size of the tree plot window  */
if( ! (file_plot || doing_print) ) {
		Nlm_RecT ob_rect; Nlm_BaR bar;
		Nlm_ObjectRect(fd_nj_plot->tree_plot, &ob_rect);
		bar = Nlm_GetSlateVScrollBar((Nlm_SlatE)fd_nj_plot->tree_plot);
		v = (double)Nlm_GetBarValue(bar);
		m = (double)Nlm_GetBarMax(bar);
		physx = ob_rect.right - ob_rect.left;
		physy = (ob_rect.bottom - ob_rect.top);
		zoom = fd_nj_plot->zoomvalue;
	}
#endif

factor = 1e99;
for(num=0; num<=fd_nj_plot->notu; num++) {
	if(fd_nj_plot->profs[num] == 0) continue;
	if(fd_nj_plot->profs[num] <= 0) continue; /* new */
	value = (physx - fd_nj_plot->widnames[num] - 2 * fd_nj_plot->ascent) / fd_nj_plot->profs[num];
	if(value < factor) factor = value;
	}

physx_min  = fd_nj_plot->ascent;
physx_corr = maxx * factor + physx_min;
if(physx_corr <= physx_min + 1) {
	physx_corr = physx_min + 1;
	}
physy_corr = physy - 2*fd_nj_plot->char_height;
physy_min  = fd_nj_plot->char_height;
#ifdef TTY
physy_corr = physy - 3;
physy_min  = 0;
#endif
/* scale_window(0, maxx, 0, maxy, 
	physx_min, physx_corr, physy_min, physy_corr); */
scale_window(0, maxx, 
	((m-v)/m) * (maxy - maxy/zoom), 
	((m-v)/m) * (maxy - maxy/zoom) + maxy/zoom, 
	physx_min, physx_corr, physy_min, physy_corr);

for(num=0; num<=fd_nj_plot->totnoms; num++) {
	mydrawstring( (fd_nj_plot->noms+num)->x, (fd_nj_plot->noms+num)->y , (fd_nj_plot->noms+num)->nom,
			(fd_nj_plot->noms+num)->disp_option );
	}
for(num=0; num<=fd_nj_plot->tottraits; num++) {
	p= fd_nj_plot->traits+num;
	moveto(p->xd,p->yd);
	lineto(p->xf,p->yf);
	}

if(fd_nj_plot->has_br_length){
	/* echelle */
	draw_scale();
	}
convert_mem_point();
} /* end of do_plot */



int calc_text_size(char *text, int *pheight, int *pascent)
{
int ascent, descent, width;

if(file_plot) {
	int font_size;
	font_size = list_font_size[fd_nj_plot->font_size_rank];
#ifdef WIN_MAC
	{
	static FontInfo macfontinfo;
	width = TextWidth(text, 0, strlen(text));
	GetFontInfo(&macfontinfo);
	ascent = macfontinfo.ascent;
	descent = macfontinfo.descent;
	}
#elif defined(TTY)
	width = strlen(text);
	ascent = 1; descent = 0;
#elif defined(NO_GUI)
	width = strlen(text) * font_size * postscript_ratio + 0.5;
	ascent = font_size; descent = 0;
#else

#ifdef WITH_PDF
	if(pdf_plot && pdf != NULL) width = PDF_stringwidth(pdf, text, pdf_font, font_size);
	else
#endif /* faire calculer taille en mode ecran et appliquer ratio ecran->postscript */
		width = Nlm_TextWidth(text, strlen(text) ) * postscript_ratio + 0.5;
	ascent = font_size; descent = 0;
#endif
	}
#ifndef NO_GUI
else	{
#ifdef WIN_MSWIN
	HDC picHDC = Nlm_currentHDC;
	if(doing_copy) {
		Nlm_currentHDC = Nlm_GetPicWinHDC();
		Nlm_SelectFont(fd_nj_plot->current_font);
		}
#endif
    ascent = Nlm_Ascent(); descent = Nlm_Descent();
	width = Nlm_TextWidth(text, strlen(text) );
#ifdef WIN_MSWIN
	if(doing_copy) Nlm_currentHDC = picHDC;
#endif
	}
#endif
if(pheight != NULL) *pheight = ascent + descent;
if(pascent != NULL) *pascent = ascent;
return width;
}


static int lastmovex, lastmovey;
#define bound_int2(i) ( i > INT2_MAX ? INT2_MAX : (i < INT2_MIN ? INT2_MIN : i) )

void plotstring(char *nom)
{
if(file_plot) { 
#ifdef WIN_MAC
	static char copy[255];
	int l;
	l = strlen(nom);
	*copy = l;
	memcpy(copy + 1, nom, l + 1);
	DrawString( (ConstStr255Param) copy);
#elif defined(TTY)
	memcpy(tty_page[tty_y]+tty_x+1, nom, strlen(nom));
#else
#ifdef WITH_PDF
if(pdf_plot) {
	PDF_stroke(pdf);
	PDF_show_xy(pdf, nom, lastmovex, lastmovey);
	}
else
#endif
	fprintf(plotfile,"(%s) show\n",nom);
#endif
	}
#ifndef NO_GUI
else 	{
#ifdef WIN_MSWIN
	int y; y = lastmovey - (doing_copy ? fd_nj_plot->ascent : 0) ;
	y = bound_int2(y);
  	Nlm_PaintStringEx(nom, lastmovex, y );
#else
        Nlm_PaintString(nom);
#endif
	}
#endif
}


void dir_moveto(double x,double y)
/* move to physical coord x,y */
{
if(file_plot) {
	int xi,yi;
	xi=x;
#ifdef WIN_MAC
	xi += margin;
	yi = physy - y + margin;
	MoveTo(xi, yi);
#elif defined(TTY)
	tty_x = x + 0.5; tty_y = y + 0.5;
#else
	yi=y;
#ifdef WITH_PDF
if(pdf_plot) {
	PDF_moveto(pdf, xi, yi);
	lastmovex = xi; lastmovey = yi;
	}
else
#endif
	fprintf(plotfile,"%d %d moveto\n",xi,yi);
#endif
}
#ifndef NO_GUI
else	{
	if(doing_print) {
		lastmovex = myrint(x) + print_rect.left;
		lastmovey = print_rect.bottom - myrint(y);
		Nlm_MoveTo(lastmovex, bound_int2(lastmovey) );	
		}
	else {
		Nlm_RecT obrect;
		ObjectRect(fd_nj_plot->tree_plot, &obrect);
		lastmovex = myrint(x) + obrect.left;
		lastmovey = obrect.bottom - myrint(y);
		if(lastmovey >= INT2_MAX) 
			lastmovey = INT2_MAX - 1;
		else if(lastmovey < 0) 
			lastmovey = 0;
		Nlm_MoveTo((Nlm_Int2)lastmovex, (Nlm_Int2)lastmovey );	
		}
	}
#endif
}


void dir_lineto(double x,double y)
/* draw line from current pos to physical coord x,y */
{
int xi,yi;
char *p;
if(file_plot) {
	xi=x;
#ifdef WIN_MAC
	xi += margin;
	yi = physy - y + margin;
	LineTo(xi, yi);
#elif defined(TTY)
	if(tty_y == (int)(y + 0.5)) {
		xi = x + 0.5;
		while(tty_x != xi) { 
			p = tty_page[tty_y] + tty_x; 
			*p = (*p == '|' ? '+' : '-'); 
			if(tty_x < xi) ++tty_x; else --tty_x;
			}
		p = tty_page[tty_y] + tty_x; 
		*p = (*p == '|' ? '+' : '-'); 
		}
	else {
		yi = y + 0.5;
		while(tty_y != yi) { 
			p = tty_page[tty_y] + tty_x; 
			*p = (*p == '-' ? '+' : '|'); 
			if(tty_y < yi) ++tty_y; else --tty_y;
			}
		p = tty_page[tty_y] + tty_x; 
		*p = (*p == '-' ? '+' : '|'); 
		}
#else
	yi=y;
#ifdef WITH_PDF
if(pdf_plot) {
	PDF_lineto(pdf, xi, yi);
	PDF_stroke(pdf);
	}
else
#endif
	fprintf(plotfile,"%d %d lineto stroke\n",xi,yi);
#endif
	}
#ifndef NO_GUI
else	{
	if(doing_print) {
		xi =	myrint(x) + print_rect.left; yi = print_rect.bottom - myrint(y) ;
		yi = bound_int2(yi);
		Nlm_LineTo(	xi, yi );
		}
	else {
		Nlm_RecT obrect;
		ObjectRect(fd_nj_plot->tree_plot, &obrect);
			{int j;
			j = obrect.bottom - myrint(y);
			if(j >= INT2_MAX)
				j = INT2_MAX - 1;
			else if( j < 0)
				j = 0;
			Nlm_LineTo(myrint(x) + obrect.left, (Nlm_Int2)j);
			}
		}
	}
#endif
}



#ifdef WIN_MAC

static void closewinproc(Nlm_WindoW w)
{
Nlm_Remove(Nlm_ParentWindow(w));
}

void show_apropos_njplot(Nlm_IteM item)
{
Nlm_WindoW win; int larg; Nlm_DoC doc;
char ligne[100], *p;
static char txt[5000];

if(help_file == NULL) return;
win = Nlm_FixedWindow( -50, -33, -3, -3, "about Njplot", closewinproc);
	Nlm_SelectFont(Nlm_programFont);
	larg = 81 * Nlm_TextWidth("Q", 1);
	doc = Nlm_DocumentPanel(win, larg, 25 * Nlm_LineHeight() );
rewind(help_file);
fgets(ligne, sizeof(ligne), help_file);
p = txt;
while(fgets(ligne, sizeof(ligne), help_file) != NULL) {
	if(strncmp(ligne, ">>>", 3) == 0) break;
	if(strncmp(ligne, "Version:", 8) == 0) sprintf(ligne, "Version: %s\n", NJPLOTVERSION);
	strcpy(p, ligne); p += strlen(ligne);
	}
Nlm_AppendText(doc, txt, NULL, NULL, Nlm_programFont);
Nlm_Break(win);
Nlm_DefaultButton(win, "ok", (Nlm_BtnActnProc)closewinproc);
Nlm_Show(doc);
Nlm_Show(win);
}


void mac_timer(void)
{
static FD_nj_plot *current, *previous = NULL;
current = (FD_nj_plot *)Nlm_GetWindowExtra(Nlm_FrontWindow());
if(current != previous && current != NULL && current->tag == tag_njplot ) {
		previous = current;
		Nlm_SetValue(current->choix_taille, current->font_size_rank);
		Nlm_SetValue(current->choix_font, current->font_family_rank);
		Nlm_SetStatus(current->bold_item, strchr(current->font_bold_italic, 'b') != NULL);
		Nlm_SetStatus(current->italic_item, strchr(current->font_bold_italic, 'i') != NULL);
		}
}


void window_callback(IteM item)
{
Nlm_WindoW w = Nlm_desktopWindow;
while(w != NULL) {
	FD_nj_plot *data = Nlm_GetWindowExtra(w);
	if(data != NULL && data->tag == tag_njplot && data->win_menu_item == item) {
		fd_nj_plot = data;
		Nlm_Select(fd_nj_plot->nj_plot);
		break;
		}
	w = Nlm_GetNext(w);
	}
}


void create_win_if_needed(char *fname)
{
if( (!update_tree_w_data(NULL)) || fd_nj_plot->notu != 0) fd_nj_plot =  create_win_nj_plot();
init_tree( fname, NULL);
Nlm_Show(fd_nj_plot->nj_plot);
Nlm_ObjectRect(fd_nj_plot->nj_plot, &(fd_nj_plot->old_rect));
}

int is_macosx(void)
{
long val;

Gestalt (gestaltSystemVersion, &val);
return (val >= 0x1000);
}


int crefpict(char *fname, PicHandle picture)
{
QDPictRef myqdpict;
int err;

myqdpict = MyPictToQDPict(picture);
if(myqdpict != NULL) err = MyQDPictToPDFfile(myqdpict, fname);
else err = TRUE;
KillPicture(picture);
return err;
}


void example_tree(IteM unused)
{
char *p;
FILE *in;
if(!update_tree_w_data(NULL)) return;
if(fd_nj_plot->notu > 0) free_tree();
fd_nj_plot->notu = 0;
p = MG_GetBundleResourcesDir();
strcat(p, "/example.phb");
in = fopen(p, "r");
if(in == NULL) {
	err_message("Sorry, no example available");
	return;
	}
fclose(in);
init_tree(p, NULL);
}


void copy_plot(IteM unused)
{
extern void MyCopyPictToClipboard (PicHandle mypicture);

if(!update_tree_w_data(NULL)) return;
if(fd_nj_plot->notu == 0) return;
doing_copy = TRUE;
pict_plot();
	
MyCopyPictToClipboard (mypicture );

KillPicture(mypicture);
file_plot = FALSE;
fd_nj_plot->need_runtree = TRUE;
doing_copy = FALSE;
ArrowCursor();
}


void pict_plot(void)
{
char mess[250];
double currx, curry;
int erreur;
	int macfont;
	Style mystyle=normal;
	margin = 30;
	myrect.top=10;
	myrect.left=10;
	myrect.bottom=760 * page_count;
	myrect.right=500;
	ClipRect(&myrect);
	mypicture = OpenPicture(&myrect);
	PenNormal();
	macfont = GetValue(fd_nj_plot->choix_font);
	if(macfont == courier)
		macfont=22;
	else if(macfont == helvetica)
		macfont=21;
	else if(macfont == times)
		macfont=20;
	else
		macfont=22;
	TextFont(macfont);
	TextSize( list_font_size[fd_nj_plot->font_size_rank] );
    if( GetStatus(fd_nj_plot->bold_item) )
    	mystyle += bold;
    if( GetStatus(fd_nj_plot->italic_item) )
    	mystyle += italic;
    TextFace(mystyle);
	physx = myrect.right-myrect.left - 2 * margin;
	physy = myrect.bottom-myrect.top - 2 * margin;

WatchCursor();
file_plot = TRUE;
fd_nj_plot->totnoms = fd_nj_plot->tottraits = fd_nj_plot->totpoints = -1;
end_br_length = fd_nj_plot->br_length_txt;
currx = 0.;
maxx = 0.;
nexty = -fd_nj_plot->deltay;
if(fd_nj_plot->subtree_center == NULL)
	mem_plot(NULL, fd_nj_plot->racine, currx, &curry);
else
	mem_plot(fd_nj_plot->subtree_ascend, fd_nj_plot->subtree_center, currx, &curry);
do_plot();
file_plot = FALSE; /* pour eviter boucle infernale! */
	ClosePicture();
if(doing_copy) return;
	erreur = crefpict(plotfilename, mypicture);
	if( !erreur ) sprintf(mess, "Tree plot is now in file %s in PDF format", 
			extract_filename(plotfilename) );
	else sprintf(mess, "Error while writing to file %s", 
			extract_filename(plotfilename) );
	err_message(mess);

fd_nj_plot->need_runtree = TRUE;
ArrowCursor();
}

#endif


#ifdef MICRO

#ifdef WIN_MSWIN
#define TITLEFONT "Times New Roman,10"
#else
#define TITLEFONT "Times,10"
#endif

void print_title(int x, int y, char *text, 
		Nlm_FonT title_font, int p, int totp)
{
static char ligne[200];
time_t heure;
int h;

time (&heure);
#ifdef WIN_MAC
text = mac_fname_to_roman(text);
#endif
sprintf(ligne, "njplot    %s    %s", text, ctime(&heure) );
h = strlen(ligne) - 1; ligne[h] = 0;
if(totp > 1) sprintf(ligne + h, " Page %d of %d", p, totp);
Nlm_Black();
Nlm_SelectFont(title_font);
Nlm_PaintStringEx(ligne, x, y);
}


void print_plot(IteM item)
{
double currx, curry;
Nlm_WindoW w;
int h, page, superpos;
Nlm_RecT true_print_rect, title_rect;
Nlm_FonT title_font;

if(!update_tree_w_data(item)) return;
if(fd_nj_plot->notu == 0) return;
w = Nlm_StartPrinting();
if(w == NULL) return;
Nlm_WatchCursor();
doing_print = TRUE;
Nlm_PrintingRect(&true_print_rect);
title_rect = true_print_rect;
title_font = Nlm_ParseFont(TITLEFONT);
Nlm_SelectFont(title_font);
calc_text_size("Mq", &fd_nj_plot->char_height, &fd_nj_plot->ascent);
title_rect.bottom = title_rect.top + fd_nj_plot->ascent - 1;
true_print_rect.top += fd_nj_plot->char_height;

print_rect.left = true_print_rect.left; print_rect.right = true_print_rect.right;
print_rect.top = true_print_rect.top; print_rect.bottom = true_print_rect.bottom;
h = print_rect.bottom - print_rect.top + 1;
superpos = h / 30;
print_rect.bottom = print_rect.top + (page_count - 1) * (h - superpos) + h - 1;

	fd_nj_plot->totnoms = fd_nj_plot->tottraits = fd_nj_plot->totpoints = -1;
	end_br_length = fd_nj_plot->br_length_txt;
	currx = 0.;
	maxx = 0.;
	nexty = -fd_nj_plot->deltay;
	if(fd_nj_plot->subtree_center == NULL)
		mem_plot(NULL, fd_nj_plot->racine, currx, &curry);
	else
		mem_plot(fd_nj_plot->subtree_ascend, fd_nj_plot->subtree_center, currx, &curry);

for(page = 0; page < page_count; page++) {
	Nlm_StartPage();
	Nlm_ClipPrintingRect( &title_rect);
	print_title(title_rect.left, title_rect.bottom, fd_nj_plot->tree_name, 
		title_font, page + 1, page_count);
	Nlm_ClipPrintingRect( &true_print_rect);
	Nlm_SelectFont(fd_nj_plot->current_font);
	Nlm_Gray();
	Nlm_FrameRect(&true_print_rect);
	Nlm_Black();
	physx = print_rect.right - print_rect.left;
	physy = print_rect.bottom - print_rect.top;
	do_plot();
	Nlm_EndPage();
	print_rect.top -= (h - superpos); print_rect.bottom -= (h - superpos);
	}
Nlm_EndPrinting(w);
fd_nj_plot->need_runtree = TRUE;
doing_print = FALSE;
ArrowCursor();
}

#endif


#ifdef WIN_MSWIN
void copy_plot(IteM item)
{
double currx, curry;
Nlm_WindoW w;
Nlm_RecT r;

if(!update_tree_w_data(item)) return;
if(fd_nj_plot->notu == 0) return;
/* dans la picture, les calculs de taille de caractere ne marchent pas!
il faut les faire dans le DC de la fenetre
*/
WatchCursor();
doing_copy = TRUE;
Nlm_GetPosition( fd_nj_plot->tree_plot , &r);
w = Nlm_StartPicture(&r);
fd_nj_plot->totnoms = fd_nj_plot->tottraits = fd_nj_plot->totpoints = -1;
end_br_length = fd_nj_plot->br_length_txt;
currx = 0.;
maxx = 0.;
nexty = -fd_nj_plot->deltay;
Nlm_Black();
Nlm_SelectFont(fd_nj_plot->current_font);
if(fd_nj_plot->subtree_center == NULL)
	mem_plot(NULL, fd_nj_plot->racine, currx, &curry);
else
	mem_plot(fd_nj_plot->subtree_ascend, fd_nj_plot->subtree_center, currx, &curry);
do_plot();
Nlm_EndPicture(w);
fd_nj_plot->need_runtree = TRUE;
doing_copy = FALSE;
ArrowCursor();
}
#endif


char *preptree(char *fname)
{
int i, c, steparbre, maxarbre, maxlname, v;
FILE *njfile;
char *arbre, *der_arbre, *finarbre, *tmp;
static char message[200];
char *last_bootstrap, *p;

if( (njfile=fopen(fname,"r")) == NULL ){
	sprintf(message, "Tree file %s not found.", fname);
	return (message);
	}
/* recherche du debut de la description de l'arbre */
do c=fgetc(njfile);
while ( isspace(c) );
/* for fastDNAml format, skip initial [comment] */
if(c == '[') {
	do c=fgetc(njfile); while (c != ']');
	do c=fgetc(njfile); while ( isspace(c) );
	}
if ( c != '(') {
	fclose(njfile);
	goto erreur;
	}
/* lecture de l'arbre par paquets de steparbre caracteres*/
steparbre=1000;
maxarbre=steparbre;
arbre=check_alloc(maxarbre,1);
der_arbre = arbre+maxarbre;
*arbre=c;
fd_nj_plot->notu=2; i=3; v = 0;
finarbre=arbre;
while( (c=fgetc(njfile)) != EOF && c != ';') {
	if( c=='\n' || c=='\r') continue;
	if(++finarbre >= der_arbre-2) {
		maxarbre += steparbre;
		tmp=check_alloc(maxarbre,1);
		memcpy(tmp,arbre,finarbre-arbre);
		finarbre=tmp+(finarbre-arbre);
		free(arbre);
		arbre=tmp;
		der_arbre=arbre+maxarbre;
		}
	*(finarbre)=c;
	if(c == ')') fd_nj_plot->notu++;
	if(c == '(') i++;
	if(c == ',') v++;
	}
*(finarbre+1)='\0';
fclose(njfile);
if(i != fd_nj_plot->notu)goto erreur;
finarbre = nextpar(arbre) + 1;
/* memorize bootstrap value after last ) */
while(isspace(*finarbre)) finarbre++;
if(*finarbre != ';' && *finarbre != 0) {
	last_bootstrap = strdup(finarbre);
	p = strchr(last_bootstrap, ';');
	if(p != NULL) *p = 0;
	}
else last_bootstrap = NULL;

arbre = (char *)realloc(arbre, strlen(arbre) + 4 * v + 5 ); /* worst case add 4 chars for each , */
make_binary_or_unrooted(arbre);
fd_nj_plot->long_arbre_parenth = strlen(arbre);
fd_nj_plot->notu = v ; /* after this fd_nj_plot->notu = number of OTUs - 1  */
fd_nj_plot->totbranches= -1;
/* allocate all memory */
fd_nj_plot->tabtax = (struct noeud **)check_alloc(2*fd_nj_plot->notu +1,sizeof(struct noeud *));
fd_nj_plot->branches = (branche *)check_alloc(fd_nj_plot->notu, sizeof(branche));
fd_nj_plot->tabnames = (char **)check_alloc(2*fd_nj_plot->notu+1, sizeof(char *));
for(i=0; i<2*fd_nj_plot->notu+1; i++) *(fd_nj_plot->tabtax+i)=
			(struct noeud *)check_alloc(1,s_noeud);
fd_nj_plot->noms = (struct nom *)check_alloc(5*fd_nj_plot->notu+1,sizeof(struct nom));
fd_nj_plot->points = (struct mon_point *)check_alloc(4*fd_nj_plot->notu+1,sizeof(struct mon_point));
fd_nj_plot->traits = (struct trait *)check_alloc(3*fd_nj_plot->notu,sizeof(struct trait));
fd_nj_plot->labels = (char **)check_alloc(fd_nj_plot->notu+1, sizeof(char *));
fd_nj_plot->widnames = (int *)check_alloc(fd_nj_plot->notu+1, sizeof(int));
fd_nj_plot->profs = (double *)check_alloc(fd_nj_plot->notu+1, sizeof(double));
fd_nj_plot->br_length_txt = (char *)check_alloc(2*fd_nj_plot->notu,10);
loadphylip(arbre, last_bootstrap);
free(arbre);
maxlname = 0; /* largeur max des noms des feuilles */

if(nextotu != fd_nj_plot->notu) bad_format("Error: incorrect tree file");

#ifdef TTY
if(tty_prefix != NULL) {
	int l, lp, i;
	char *p;
	lp = strlen(tty_prefix);
	for(i=0; i <= fd_nj_plot->notu; i++) {
		l = strlen(fd_nj_plot->labels[i]) + lp;
		p = (char *)check_alloc(l + 1, 1);
		strcpy(p, tty_prefix);
		strcat(p, fd_nj_plot->labels[i]);
		free(fd_nj_plot->labels[i]);
		fd_nj_plot->labels[i] = p;
		}
	}
#endif

for(i = 0; i <= nextotu; i++) {
	c = strlen(fd_nj_plot->labels[i]);
	if(c > maxlname) maxlname = c;
	}
for(i = 0; i < 2*fd_nj_plot->notu+1; i++) {
	fd_nj_plot->tabnames[i] = (char *)check_alloc(maxlname + 2 + 1, 1); /*2=place pour # */
	}
	
if(!fd_nj_plot->rooted) {
	fd_nj_plot->racine = *(fd_nj_plot->tabtax+(++num_noeud));
	if(num_noeud >= 2*fd_nj_plot->notu + 1) bad_format("Error: incorrect tree file");
	if(fd_nj_plot->has_br_length) {
		place_midpoint_root();
		fd_nj_plot->root_num = -1;
		}
	else	{
/* ancienne version: derniere espece est groupe externe
*/
		fd_nj_plot->racine->v3 = NULL;
		fd_nj_plot->root_num = fd_nj_plot->notu; 
		}
	}
else	{
	fd_nj_plot->racine = *(fd_nj_plot->tabtax+num_noeud);
	fd_nj_plot->root_br_l= fd_nj_plot->racine->l1 + fd_nj_plot->racine->l2;
	fd_nj_plot->root_num = num_noeud;
	if(!fd_nj_plot->has_br_length) calc_brl_for_lengthless(fd_nj_plot->racine, NULL);
/* y a-t-il un bootstrap sur l'une des branches racine ? */
	i = get_br_from_bouts(fd_nj_plot->racine, fd_nj_plot->racine->v1); 
	if(i == -1) i = get_br_from_bouts(fd_nj_plot->racine, fd_nj_plot->racine->v2);
	if(i != -1 && get_br_from_bouts(fd_nj_plot->racine, NULL) == -1) {
		fd_nj_plot->branches[i].bouta = fd_nj_plot->racine->v1;
		fd_nj_plot->branches[i].boutb = fd_nj_plot->racine->v2;
		}
	} 
if(fd_nj_plot->notu+1<3) return ("Tree should contain at least 3 elements.");
fd_nj_plot->subtree_notu = fd_nj_plot->notu;
#ifndef NO_GUI
if(!pdf_plot_only) scrollnotu(fd_nj_plot->notu);
#endif
return NULL;

erreur:
#ifdef WIN_MAC
return ("File or pasted data does not contain correct tree.");
#else
return ("File does not contain correct tree data.");
#endif
} /* end of preptree */


char *check_alloc(int nbrelt, int sizelt)
{
char *retval;
if( (retval=(char *)calloc(nbrelt,sizelt)) != NULL ) return retval;
err_message("ERROR: Not enough memory.");
exit(1);
}


void loadphylip(char *arbre, char *last_bootstrap)
{
char *deba,*debb,*debc, *finarbre;
struct noeud *p1, *p2, *p3, *p;
branche *int_br_g, *int_br_d;

fd_nj_plot->has_br_length = 2;
fd_nj_plot->has_internal = FALSE;
/* ignore all stuff after last closing parenthesis 
(needed for fastDNAml output)
*/
finarbre= nextpar(arbre);
fd_nj_plot->rooted=0;
deba=arbre+1;
debb=deba;
while(*debb != ',') {
	if(*debb == 0) bad_format("Incorrect tree file");
	if(*debb == '(')debb=nextpar(debb);
	debb++;
	}
debb++;
debc=debb;
while(*debc != ',' && debc<finarbre) {
	if(*debc == '(')debc=nextpar(debc);
	debc++;
	}
if(*debc==',') {
/* the tree is unrooted <==> it has 3 subtrees at its bottommost level */
	debc++;
	}
else	{
/* the tree is rooted */
	debc=finarbre+1;
	fd_nj_plot->rooted=1;
/*
	fd_nj_plot->notu--;  //useless now
*/
	}

num_noeud=fd_nj_plot->notu;
nextotu= -1;
p1=unrootedset(deba,debb-2,&int_br_g);
p2=unrootedset(debb,debc-2,&int_br_d);
p = *(fd_nj_plot->tabtax+(++num_noeud));
if(num_noeud >= 2*fd_nj_plot->notu + 1) bad_format("Error: incorrect tree file");
p->v1=p1; p1->v3=p; p->l1=p1->l3;
if(int_br_g!=NULL) { int_br_g->bouta=p; int_br_g->boutb=p1; }
p->v2=p2; p2->v3=p; p->l2=p2->l3;
if(int_br_d!=NULL) { int_br_d->bouta=p; int_br_d->boutb=p2; }
if(!fd_nj_plot->rooted) {
	p3=unrootedset(debc,finarbre-1,&int_br_g);
	if(int_br_g!=NULL) { int_br_g->bouta=p; int_br_g->boutb=p3; }
	p->v3=p3; p3->v3=p; p->l3=p3->l3;
	}
else	{
	p->v3=NULL;
/* recherche d'un dernier label interne */
	debc=finarbre+1;
	while(*debc!=0 && *debc!=':' && *debc!='[') debc++;
	if(debc-finarbre>1) {
		int l=debc-finarbre-1;
		fd_nj_plot->has_internal = TRUE;
		fd_nj_plot->totbranches++;
		fd_nj_plot->branches[fd_nj_plot->totbranches].br_label=check_alloc(l+1,1);
		memcpy(fd_nj_plot->branches[fd_nj_plot->totbranches].br_label,finarbre+1,l);
		fd_nj_plot->branches[fd_nj_plot->totbranches].br_label[l]=0;
		fd_nj_plot->branches[fd_nj_plot->totbranches].bouta=p1;
		fd_nj_plot->branches[fd_nj_plot->totbranches].boutb=p2;
		}
	}
if(fd_nj_plot->rooted && last_bootstrap != NULL) {
/* attach last_bootstrap to branch racine <--> NULL */
		fd_nj_plot->totbranches++;
		fd_nj_plot->branches[fd_nj_plot->totbranches].br_label = strdup(last_bootstrap);
		fd_nj_plot->branches[fd_nj_plot->totbranches].bouta=p;
		fd_nj_plot->branches[fd_nj_plot->totbranches].boutb=NULL;
	}
if(last_bootstrap != NULL) free(last_bootstrap);
}


struct noeud *unrootedset(char *deb, char *fin, branche **p_int_br)
{
struct noeud *p, *pp;
char *virg, *ferme;
branche *int_br;
static int l;
static double brlength;

*p_int_br=NULL;
while(*deb==' ')deb++;
while(*fin==' ')fin--;
if(*deb != '(') { /* une feuille */
	virg = strchr(deb, ':');
	if(virg != NULL && virg < fin) {
//		if(fd_nj_plot->has_br_length == 0) goto problem;
		sscanf(virg+1, "%le", &brlength);
		fd_nj_plot->has_br_length=1;
		}
	else	{
//		if(fd_nj_plot->has_br_length == 1) goto problem;
		brlength = 1;
		fd_nj_plot->has_br_length=0;
		virg = fin + 1;
		}
	virg--;
	while(*deb==' ' || *deb=='\'')deb++;
	if( virg-1 >= deb && *virg == '\'' ) virg--;
	l = virg-deb+1;
	fd_nj_plot->labels[ ++nextotu] = (char *)check_alloc(l + 1, 1);
	memcpy(fd_nj_plot->labels[nextotu], deb, l);
	fd_nj_plot->labels[nextotu][l] = 0;
/*
if(strchr(fd_nj_plot->labels[nextotu], '(') != NULL || strchr(fd_nj_plot->labels[nextotu], ')') != NULL) {
		char *tmp = (char *)malloc(strlen(fd_nj_plot->labels[nextotu]) + 60);
		sprintf(tmp, "Error: parentheses in name: %s", fd_nj_plot->labels[nextotu]);
		bad_format(tmp);
		}
*/
	p = *(fd_nj_plot->tabtax + nextotu);
	p->l3 = brlength;
	p->v1 = p->v2 = p->v3 = NULL;
	return p;
	}
/* un noeud */
num_noeud++;
if(num_noeud >= 2*fd_nj_plot->notu + 1) bad_format("Error: incorrect tree file");
p = *(fd_nj_plot->tabtax + num_noeud);
ferme =  nextpar(deb);
virg=deb + 1;
while(*virg != ',' && virg < fin) {
	if(*virg == '(') virg=nextpar(virg);
	virg++;
	}
if(virg>=ferme) bad_format("Error: incorrect tree file");
pp = unrootedset(deb + 1, virg - 1, &int_br);
p->v1 = pp; pp->v3 = p; p->l1 = pp->l3;
if(int_br != NULL) { int_br->bouta = p; int_br->boutb = pp; }
pp = unrootedset(virg + 1, ferme - 1, &int_br);
p->v2 = pp; pp->v3 = p; p->l2 = pp->l3;
if(int_br != NULL) { int_br->bouta = p; int_br->boutb = pp; }
virg = strchr(ferme, ':');
if(virg != NULL && virg < fin) { /* traitement longueur */
//	if(fd_nj_plot->has_br_length == 0) goto problem;
	sscanf(virg+1, "%le", &brlength);
	fd_nj_plot->has_br_length=1;
	if(*fin == ']') { /* bootstrap entre [] apres longueurs */
		static char *q;
		q = fin - 1;
		while(q > virg && *q != '[') q--;
		if(*q == '[' && fin - q >= 2) {
			fd_nj_plot->has_internal = TRUE;
			fd_nj_plot->totbranches++;
			l = fin - q - 1;
			fd_nj_plot->branches[fd_nj_plot->totbranches].br_label =
				check_alloc(l+1,1);
			memcpy(fd_nj_plot->branches[fd_nj_plot->totbranches].br_label,q+1,l);
			fd_nj_plot->branches[fd_nj_plot->totbranches].br_label[l]=0;
			*p_int_br= &fd_nj_plot->branches[fd_nj_plot->totbranches];
			}
		}
	}
else	{
//	if(fd_nj_plot->has_br_length == 1) goto problem;
	brlength = 1;
	fd_nj_plot->has_br_length=0;
	virg = fin + 1;
	}
/* recherche bootstrap (internal label) */
l=virg-ferme-1;
if(l>0) {
	fd_nj_plot->has_internal = TRUE;
	fd_nj_plot->totbranches++;
	fd_nj_plot->branches[fd_nj_plot->totbranches].br_label=
		check_alloc(l+1,1);
	memcpy(fd_nj_plot->branches[fd_nj_plot->totbranches].br_label,ferme+1,l);
	fd_nj_plot->branches[fd_nj_plot->totbranches].br_label[l]=0;
	*p_int_br= &fd_nj_plot->branches[fd_nj_plot->totbranches];
	}
p->l3 = brlength;
return p;

problem:
err_message("Error: Inconsistent tree file for branch lengths.");
exit(1);
}


char *nextpar(char *pospar)
{
char *pos;
pos=pospar+1;
while(*pos != ')') {
	if(*pos == 0) bad_format("Error: unbalanced parentheses");
	if(*pos == '(') pos=nextpar(pos);
	pos++;
	}
return pos;
}


void bad_format(char *mess)
{
err_message(mess);
exit(1);
}



void make_binary_or_unrooted(char *arbre)
{
char *finarbre, *deba, *debb, *debc;

finarbre= nextpar(arbre);
*(finarbre + 1) = 0;
deba=arbre+1;
debb=deba;
while(*debb != ',') {
	if(*debb == 0) bad_format("Incorrect tree file");
	if(*debb == '(')debb=nextpar(debb);
	debb++;
	}
debb++;
debc=debb;
while(*debc != ',' && debc<finarbre) {
	if(*debc == '(')debc=nextpar(debc);
	debc++;
	}
if(*debc == ',') {
/* the tree is unrooted <==> it has 3 subtrees or more at its bottommost level */
	make_binary(arbre , debb, finarbre - 1, TRUE);
	make_binary(arbre , deba, debb - 2, TRUE);
	}
else make_binary(arbre , deba, finarbre - 1, TRUE);
return;
}


int make_binary(char *arbre, char *debut, char *fin, int go_down)
{
int virg, l, retval;
char *p, *q;

p = debut; virg = 0; retval = 0;
while(p < fin) {
	if(*p == ',') virg++;
	else if(*p == '(')	{
		q = nextpar(p);
		if(go_down) {
			l = make_binary(arbre, p+1, q-1, TRUE);
			fin += l;
			retval += l;
			p = q + l;
			}
		else p = q;
		}
	p++;
	}
if(virg > 1) { /* multifurcation */
	/* recherche de la 2eme virgule */
	p = debut; l = 0;
	while(TRUE) {
		if(*p == ',') {
			l++;
			if(l == 2) break;
			}
		else if(*p == '(')	{
			p = nextpar(p);
			}
		p++;
		}
	l = strlen(p);
	memmove(p + 4, p, l + 1);
	memmove(debut + 1, debut, p - debut);
	*debut = '(';
	memcpy(p + 1, "):0", 3);
	fin += 4;
	retval += 4;
	if(virg > 2) retval += make_binary(arbre, debut, fin, FALSE);  
	}
return retval;
}


void mydrawstring(double x, double y, char *nom, char option)
{
static double px,py;

ch_echelle(x,y,&px,&py);
#ifndef TTY
if(option == '1') {
/* ecrire une chaine en la montant d'un chouia */
	py += fd_nj_plot->char_height/6.;
	}
else if(option == 'c' || option == 'r')	{
/* ecrire une chaine en la centrant vert. sur cette position*/
	py -= fd_nj_plot->char_height/3.;
	px += fd_nj_plot->char_height/6.;
	}
else if(option == 't')	{
	py -=  (5./6.) * fd_nj_plot->char_height;
	}
else if(option == 'b')	{
	py += fd_nj_plot->char_height/6.;
	}
#endif
dir_moveto(px,py);
#ifndef NO_GUI
if(option == 'r') 
	Nlm_Red();
#endif
plotstring(nom);
#ifndef NO_GUI
if(option == 'r') 
	Nlm_Black();
#endif
}


void moveto(double x,double y)
/* move from current pos to logical coord x,y */
{
static double px,py;
ch_echelle(x,y,&px,&py);
dir_moveto(px,py);
}


void lineto(double x,double y)
/* draw line from current pos to logical coord x,y */
{
static double px,py;
ch_echelle(x,y,&px,&py);
dir_lineto(px,py);
}

					
int calc_brl_for_lengthless(struct noeud *centre, struct noeud *pere)
/* Recursively computes branch lengths of a lengthless tree having some branches fixed to 0 to allow
multifurcations so that all tips align to the right of the plot.
*/
{
int n1 = 0, n2, depth;
volatile double l;
if(centre->v1 == NULL && centre->v2 == NULL) return 0;//a leaf
//rearrange with centre->v3 towards root
if(centre->v1 == pere) {
centre->v1 = centre->v3;
centre->v3 = pere;
l = centre->l3;
centre->l3 = centre->l1;
centre->l1 = l;
}
else if(centre->v2 == pere) {
centre->v2 = centre->v3;
centre->v3 = pere;
l = centre->l3;
centre->l3 = centre->l2;
centre->l2 = l;
}
n1 = calc_brl_for_lengthless(centre->v1, centre);
n2 = calc_brl_for_lengthless(centre->v2, centre);
depth = n1;
if(centre->l1 != 0) depth++;
if(depth < n2) depth = n2;
if(centre->l2 != 0 && n2 + 1 > depth) depth++;
if(centre->l1 != 0) {
centre->l1 = depth - n1; 
centre->v1->l3 = depth - n1;
}
else if(depth - n1 > 0) add_value_downstream(centre->v1, depth - n1);
if(centre->l2 != 0) {
centre->l2 = depth - n2; 
centre->v2->l3 = depth - n2;
}
else if(depth - n2 > 0) add_value_downstream(centre->v2, depth - n2);
return depth;
}

void add_value_downstream(struct noeud *centre, int value)
{
if(centre->l1 != 0) {
	centre->l1 += value;
	centre->v1->l3 = centre->l1;
	}
else add_value_downstream(centre->v1, value);
if(centre->l2 != 0) {
	centre->l2 += value;
	centre->v2->l3 = centre->l1;
	}
else add_value_downstream(centre->v2, value);
}
					

void place_midpoint_root(void)
/* enraciner l'arbre sans racine en cherchant son centre 
*/
{
struct noeud *aux;
double laux;

current_best_diff = VERY_BIG;
if(outotu != NULL) {
	int i; 
	for(i=0; i <= fd_nj_plot->notu; i++) {
		if(fd_nj_plot->labels[i] != NULL && strcmp(fd_nj_plot->labels[i], outotu) == 0) break;	
		}
	if(i > fd_nj_plot->notu) {
		fprintf(stderr, "Error: outotu %s not found in tree\n", outotu);
		exit(1);
		}
	process_branche(fd_nj_plot->tabtax[i], fd_nj_plot->tabtax[i]->v3, fd_nj_plot->tabtax[i]->l3);
	}
else {
	current_cote1=current_cote2=NULL;
	parcourir_branches(*fd_nj_plot->tabtax, NULL);
	}
fd_nj_plot->rooted = TRUE;
fd_nj_plot->root_br_l = current_br_length;
/* il faut toujours que la racine soit telle que racine->v1->v3=racine */
if (current_cote1->v1 == current_cote2 ) {
/* echanger les voisins v1 et v3 de cote1 */
	aux=current_cote1->v1;
	current_cote1->v1=current_cote1->v3;
	current_cote1->v3=aux;
	laux=current_cote1->l1;
	current_cote1->l1=current_cote1->l3;
	current_cote1->l3=laux;
	}
else if (current_cote1->v2 == current_cote2) {
/* echanger les voisins v2 et v3 de cote1 */
	aux=current_cote1->v2;
	current_cote1->v2=current_cote1->v3;
	current_cote1->v3=aux;
	laux=current_cote1->l2;
	current_cote1->l2=current_cote1->l3;
	current_cote1->l3=laux;
	}
current_cote1->v3 = fd_nj_plot->racine;

if (current_cote2->v1 == current_cote1 )
	current_cote2->v1 = fd_nj_plot->racine;
else if (current_cote2->v2 == current_cote1)
	current_cote2->v2 = fd_nj_plot->racine;
else
	current_cote2->v3 = fd_nj_plot->racine;
fd_nj_plot->racine->v1=current_cote1;
fd_nj_plot->racine->v2=current_cote2;
fd_nj_plot->racine->v3=NULL;
fd_nj_plot->racine->l3=0;
/* avoid very unbalanced division of root branch */
if(current_balance > MAX_FRAC) current_balance = MAX_FRAC;
else if(current_balance < 1 - MAX_FRAC) current_balance = 1 - MAX_FRAC;
fd_nj_plot->racine->l1 = current_br_length * current_balance;
fd_nj_plot->racine->l2 = current_br_length - fd_nj_plot->racine->l1;
}


void parcourir_branches(struct noeud *centre, struct noeud *origine)
/* parcourir recursivement toutes les branches de l'arbre sans racine
a partir de centre et dans la direction opposee a son voisin origine
*/
{
if(centre==NULL) return;
if(centre->v1!=origine) {
	process_branche(centre,centre->v1,centre->l1);
	parcourir_branches(centre->v1,centre);
	}
if(centre->v2!=origine) {
	process_branche(centre,centre->v2,centre->l2);
	parcourir_branches(centre->v2,centre);
	}
if(centre->v3!=origine) {
	process_branche(centre,centre->v3,centre->l3);
	parcourir_branches(centre->v3,centre);
	}
}


void process_branche(struct noeud *cote1, struct noeud *cote2, double length)
/* calculer la variance des distances racine - feuilles si racine est sur branche cote1 - cote2
et memoriser la meilleure branche dans les variables globales
*/
{
double x ; /* partage de branche par facteur x */
moments m1, m2;
double A, B, C; /* variance(x) = A x x + B x + C */
double mini_val; /* meilleure variance pour toutes valeurs de x entre 0 et 1 */

if(cote1 == NULL || cote2 == NULL) return;
if( length < 0 ) length = 0;

m1 = stat_from_node(cote2, cote1);
m2 = stat_from_node(cote1, cote2);
A = 4 * m1.N * ( m2.N * length ) * length;
B = 4 * length * ( m2.N * m1.somme - m1.N * m2.somme - length * m1.N * m2.N);
C = (m1.N + m2.N) * (m1.carres + m2.carres) + m1.N * length * m2.N * length +
	2 * m1.N * length * m2.somme - 2 * m2.N * length * m1.somme -
	(m1. somme + m2.somme) * (m1. somme + m2.somme);
if(A < 1e-20) {
	x = 0.5;
	mini_val = VERY_BIG * 0.99 ;
	}
else {
	x = - B / (2 * A);
	mini_val = C - (B * B) / 4 / A;
	}
if( x < 0 ) {
	x = 0; mini_val = C;
	}
else if( x > 1) {
	x = 1; mini_val = A + B + C;
	}

if(mini_val < current_best_diff ) {
	current_best_diff = mini_val;
	current_cote1 = cote1;
	current_cote2 = cote2;
	current_br_length = length;
	current_balance = x;
	}
}


/* no longer useful */
double get_length_down(struct noeud *pere, struct noeud *racine)
/* compute the average length of the tree down a node */
{
if(racine == NULL) return 0.0;
else	{
	struct noeud *gauche, *droite;
	double bg,bd,lg,ld;
	if( racine->v1 == pere ) {
		gauche =racine->v2; droite = racine->v3;
		bg = (racine->l2); bd = (racine->l3);
		}
	else if( racine->v2 == pere ) {
		gauche =racine->v1; droite = racine->v3;
		bg = (racine->l1); bd = (racine->l3);
		}
	else	{
		gauche =racine->v1; droite = racine->v2;
		bg = (racine->l1); bd = (racine->l2);
		}
 	/* conserver cette ecriture sinon plante sur PC */
 	lg = ld = 0;
 	if(gauche != NULL) lg = get_length_down(racine,gauche);
 	lg += bg;
 	if(droite != NULL) ld = get_length_down(racine,droite);
 	ld += bd;
	return ((lg+ld)/2);
	}
}


moments stat_from_node(struct noeud *pere, struct noeud *racine)
/* compute the moments down a node to all descending tips */
{
struct noeud *gauche, *droite;
moments m;
 static  moments mtmp;
double bg, bd;

	if( racine->v1 == pere ) {
		gauche =racine->v2; droite = racine->v3;
		bg = (racine->l2); bd = (racine->l3);
		}
	else if( racine->v2 == pere ) {
		gauche =racine->v1; droite = racine->v3;
		bg = (racine->l1); bd = (racine->l3);
		}
	else	{
		gauche =racine->v1; droite = racine->v2;
		bg = (racine->l1); bd = (racine->l2);
		}
 	if(gauche != NULL) {
		mtmp = stat_from_node(racine, gauche);
		m.N = mtmp.N; /* nbre de fils */
		m.somme = mtmp.somme + bg * mtmp.N;
		m.carres = mtmp.carres + 2 * bg * mtmp.somme + mtmp.N * bg * bg;

		mtmp = stat_from_node(racine, droite);
		m.N += mtmp.N; /* nbre de fils */
		m.somme += mtmp.somme + bd * mtmp.N;
		m.carres += mtmp.carres + 2 * bd * mtmp.somme + mtmp.N * bd * bd;
 		}
 	else {
 		m.N = 1;
 		m.somme = 0;
 		m.carres = 0;
 		}

return m;
}


void runtree(void)
{
double currx, curry;
int i;
struct noeud *p1, *p2;
double b1,b2,frac_gauche;
if(!fd_nj_plot->rooted) {
	/* place root at user-chosen place: node # fd_nj_plot->root_num */
	fd_nj_plot->rooted = !fd_nj_plot->rooted;
	p1 = *(fd_nj_plot->tabtax+fd_nj_plot->root_num);
	p2 = p1->v3;
	fd_nj_plot->root_br_l = p1->l3;
	if(fd_nj_plot->has_br_length) {
		current_best_diff = VERY_BIG;
		process_branche(p1, p2, fd_nj_plot->root_br_l);
		}
	p1->v3 = fd_nj_plot->racine;
	if (p2->v1 == p1 )
		p2->v1 = fd_nj_plot->racine;
	else if (p2->v2 == p1)
		p2->v2 = fd_nj_plot->racine;
	else
		p2->v3 = fd_nj_plot->racine;
	fd_nj_plot->racine->v1=p1;
	fd_nj_plot->racine->v2=p2;
	fd_nj_plot->racine->v3=NULL;
	if(fd_nj_plot->has_br_length) {
		frac_gauche = current_balance;
		if(frac_gauche>MAX_FRAC) frac_gauche=MAX_FRAC;
		else if(frac_gauche<1-MAX_FRAC) frac_gauche=1-MAX_FRAC;
		b1 = frac_gauche*fd_nj_plot->root_br_l;  b2 = fd_nj_plot->root_br_l - b1;
		fd_nj_plot->racine->l1=b1; fd_nj_plot->racine->l2=b2;
		}
	else	{
		fd_nj_plot->racine->l1 = p1->l3;
		if (p2->v1 == fd_nj_plot->racine )
					fd_nj_plot->racine->l2 = p2->l1;
		else if (p2->v2 == fd_nj_plot->racine)
					fd_nj_plot->racine->l2 = p2->l2;
		else
					fd_nj_plot->racine->l2 = p2->l3;
		calc_brl_for_lengthless(fd_nj_plot->racine, NULL);
		}
	}
/* initialize leave and node names */
for(i = 0; i <= 2*fd_nj_plot->notu; i++) fd_nj_plot->tabnames[i][0] = 0;
if (fd_nj_plot->choix == depl_racine) {
	for(i=0; i <= 2*fd_nj_plot->notu - 1; i++) {
//		if(i == fd_nj_plot->root_num) continue;//skip current root
		if( !fd_nj_plot->has_br_length) {
			if(fd_nj_plot->tabtax[i]->l3 == 0) continue;//skip internal multifurcation nodes that can't be broken by root
			}
		sprintf(fd_nj_plot->tabnames[i], "# ");
		}
	}
else if(fd_nj_plot->choix == permutation)
	for(i=fd_nj_plot->notu+1; i<=2*fd_nj_plot->notu; i++) sprintf(fd_nj_plot->tabnames[i],"#");
else if(fd_nj_plot->choix == subtree)
	for(i=fd_nj_plot->notu+1; i<=2*fd_nj_plot->notu-1; i++) sprintf(fd_nj_plot->tabnames[i],"#");

for(i=0; i<=fd_nj_plot->notu; i++) {
	int j;
	j = strlen(fd_nj_plot->tabnames[i]);
	strcpy(fd_nj_plot->tabnames[i] + j, fd_nj_plot->labels[i]);
	}
#ifndef ADDROOT
fd_nj_plot->totnoms = fd_nj_plot->tottraits = fd_nj_plot->totpoints = -1;
end_br_length = fd_nj_plot->br_length_txt;
currx = 0.;
maxx = 0.;
memset(fd_nj_plot->profs, 0, (fd_nj_plot->notu+1)*sizeof(double));
if(fd_nj_plot->subtree_center != NULL)
	fd_nj_plot->deltay = maxy / fd_nj_plot->subtree_notu;
else {
	fd_nj_plot->deltay = maxy / fd_nj_plot->notu;
	fd_nj_plot->subtree_notu = fd_nj_plot->notu;
	}
nexty = -fd_nj_plot->deltay;
calc_text_size("Mq", &fd_nj_plot->char_height, &fd_nj_plot->ascent);
if(fd_nj_plot->subtree_center == NULL)
	mem_plot(NULL, fd_nj_plot->racine, currx, &curry);
else
	mem_plot(fd_nj_plot->subtree_ascend, fd_nj_plot->subtree_center, currx, &curry);
converted = FALSE;
#endif
} /* end of runtree */


void mem_plot(struct noeud *pere, struct noeud *centre, double currx,
double *curry)
/* prepare tree by memorizing all graphic requests */
{
int num=0;
char *p;
struct noeud *gauche, *droite;
double bg, bd, bpere;

while( *(fd_nj_plot->tabtax+num) != centre) num++;
if(pere != NULL && 
	(centre->v1 == NULL || centre->v2 == NULL || centre->v3 == NULL) ) {
	/* orienter le noeud centre de maniere standard: centre->v3=pere */
	if( centre->v1 == pere ) {
		gauche =centre->v2; droite = centre->v3;
		bg = centre->l2; bd = centre->l3; bpere = centre->l1;
		}
	else if( centre->v2 == pere ) {
		gauche =centre->v1; droite = centre->v3;
		bg = centre->l1; bd = centre->l3; bpere = centre->l2;
		}
	else	{
		gauche =centre->v1; droite = centre->v2;
		bg = centre->l1; bd = centre->l2; bpere = centre->l3;
		}
	centre->v3=pere; centre->v1=gauche; centre->v2=droite;
	centre->l3=bpere; centre->l1=bg; centre->l2=bd;
	/* write taxon name */
	nexty += fd_nj_plot->deltay;
	mem_nom(currx,nexty,fd_nj_plot->tabnames[num],'c');
	fd_nj_plot->profs[num] = currx;
	fd_nj_plot->widnames[num] = calc_text_size(fd_nj_plot->tabnames[num], NULL, NULL);
	if(fd_nj_plot->choix==depl_racine) mem_point(currx,nexty,num+1);
	*curry = nexty;
	}
else 	{
	double yg, yd, xg, xd;

	static int doswap;
	static struct noeud *tmp;
	/* doswap vrai ssi permutation de 2 descendants necessaire ici */
	doswap = (swap != 0 && *(fd_nj_plot->tabtax+swap) == centre);
	if( centre->v1 == pere ) {
		if(doswap) {
			tmp= centre->v2;centre->v2= centre->v3;centre->v3= tmp;
			bg= centre->l2; centre->l2= centre->l3; centre->l3= bg;
			}
		gauche =centre->v2; droite = centre->v3;
		bg = centre->l2; bd = centre->l3; bpere = centre->l1;
		}
	else if( centre->v2 == pere ) {
		if(doswap) {
			tmp= centre->v1;centre->v1= centre->v3;centre->v3= tmp;
			bg= centre->l1; centre->l1= centre->l3; centre->l3= bg;
			}
		gauche =centre->v1; droite = centre->v3;
		bg = centre->l1; bd = centre->l3; bpere = centre->l2;
		}
	else	{
		if(doswap) {
			tmp= centre->v1;centre->v1= centre->v2;centre->v2= tmp;
			bg= centre->l1; centre->l1= centre->l2; centre->l2= bg;
			}
		gauche =centre->v1; droite = centre->v2;
		bg = centre->l1; bd = centre->l2; bpere = centre->l3;
		}
	/* orienter le noeud centre de maniere standard: centre->v3=pere */
	centre->v3=pere; centre->v1=gauche; centre->v2=droite;
	centre->l3=bpere; centre->l1=bg; centre->l2=bd;
	xg=currx+bg; xd=currx+bd;

	mem_plot(centre,gauche,xg,&yg);
	mem_plot(centre,droite,xd,&yd);

	mem_trait(currx,yg,xg,yg);
	mem_trait(currx,yd,xd,yd);
	mem_trait(currx,yg,currx,yd);
	*curry = (yg+yd)/2;

/* write internal fd_nj_plot->labels */
	if(fd_nj_plot->choix==show_tree && fd_nj_plot->show_bootstrap) {
		if(pere != NULL) { /* for all but root node */
			if( (p=get_br_label(centre, gauche)) != NULL ) mem_nom(currx,yg,p,'t');
			if( (p=get_br_label(centre, droite)) != NULL ) mem_nom(currx,yd,p,'b');
			}
		else {/* for root node */
			if( (p=get_br_label(centre, NULL)) != NULL ) { /* if root bootstrap exists */
				mem_nom(currx,*curry,p,'c');
				if( (p=get_br_label(centre, gauche)) != NULL ) mem_nom(currx,yg,p,'t');
				if( (p=get_br_label(centre, droite)) != NULL ) mem_nom(currx,yd,p,'b');
				}
			else {/* if no root bootstrap any bootstrap is to be centered */
				if( (p=get_br_label(centre, gauche)) != NULL ) mem_nom(currx,*curry,p,'c');
				if( (p=get_br_label(centre, droite)) != NULL ) mem_nom(currx,*curry,p,'c');
				if( (p=get_br_label(gauche, droite)) != NULL ) mem_nom(currx,*curry,p,'c');
				}
			}
		}

/* write node number */
	mem_nom(currx,*curry,fd_nj_plot->tabnames[num],'c');
	mem_point(currx,*curry,num+1);
	if(fd_nj_plot->plot_br_l && fd_nj_plot->has_br_length) {
/* write branch length */
		double min_br_l = 0.008; /* minimum branch length displayed */
		if(bg>min_br_l) {
			sprintf(end_br_length,"%.3f",bg);
			mem_nom(currx+bg/10,yg,end_br_length,'1');
			end_br_length += (strlen(end_br_length)+1);
			}
		if(bd>min_br_l) {
			sprintf(end_br_length,"%.3f",bd);
			mem_nom(currx+bd/10,yd,end_br_length,'1');
			end_br_length += (strlen(end_br_length)+1);
			}
		}
	}
}  /* end of mem_plot */

void mem_point(double x, double y, int number)
{
++fd_nj_plot->totpoints;
(fd_nj_plot->points+fd_nj_plot->totpoints)->x = x;
(fd_nj_plot->points+fd_nj_plot->totpoints)->y = y;
(fd_nj_plot->points+fd_nj_plot->totpoints)->number = number;
}

void mem_nom(double x, double y, char *nom, char option)
/* x,y: logical coordinates of beginning of string 
   nom: address of string to be displayed later
   option: 'c' use for position of center of character height
           '1' put bottom of characters at y+1 pixel
	   't' use for position of top of character
	   'b' use for position of bottom of character
*/
{
if(strlen(nom) != 0) {
	fd_nj_plot->totnoms++;
	if(x > maxx) maxx = x;
	(fd_nj_plot->noms+fd_nj_plot->totnoms)->x = x;
	(fd_nj_plot->noms+fd_nj_plot->totnoms)->y = y;
	(fd_nj_plot->noms+fd_nj_plot->totnoms)->nom = nom;
	(fd_nj_plot->noms+fd_nj_plot->totnoms)->disp_option = option;
	}
}


void mem_trait(double xd, double yd, double xf, double yf)
{
fd_nj_plot->tottraits++;
if(xd>maxx) maxx=xd;
if(xf>maxx) maxx=xf;
(fd_nj_plot->traits+fd_nj_plot->tottraits)->xd = xd;
(fd_nj_plot->traits+fd_nj_plot->tottraits)->yd = yd;
(fd_nj_plot->traits+fd_nj_plot->tottraits)->xf = xf;
(fd_nj_plot->traits+fd_nj_plot->tottraits)->yf = yf;
}


char *get_br_label(struct noeud *a, struct noeud *b)
{
int i;
for(i=0; i<=fd_nj_plot->totbranches; i++) {
	if(fd_nj_plot->branches[i].bouta==a && fd_nj_plot->branches[i].boutb==b) 
		return fd_nj_plot->branches[i].br_label;
	if(fd_nj_plot->branches[i].boutb==a && fd_nj_plot->branches[i].bouta==b) 
		return fd_nj_plot->branches[i].br_label;
	}
return NULL;
}


int get_br_from_bouts(struct noeud *a, struct noeud *b)
{
int i;
for(i=0; i<=fd_nj_plot->totbranches; i++) {
	if(fd_nj_plot->branches[i].bouta==a && fd_nj_plot->branches[i].boutb==b) 
		return i;
	if(fd_nj_plot->branches[i].boutb==a && fd_nj_plot->branches[i].bouta==b) 
		return i;
	}
return -1;
}


void free_tree(void)
{
int i;
if(fd_nj_plot->notu == 0) return;
/* de-allocate all memory */
for(i=0; i<2*fd_nj_plot->notu+1; i++)
	free(fd_nj_plot->tabtax[i]);
free(fd_nj_plot->tabtax);
for(i=0; i<fd_nj_plot->notu; i++)
	if(fd_nj_plot->branches[i].br_label != NULL) free(fd_nj_plot->branches[i].br_label);
free(fd_nj_plot->branches);
for(i = 0; i < 2 * fd_nj_plot->notu + 1; i++) free(fd_nj_plot->tabnames[i]);
free(fd_nj_plot->tabnames);
free(fd_nj_plot->widnames);
free(fd_nj_plot->noms);
free(fd_nj_plot->points);
free(fd_nj_plot->traits);
free(fd_nj_plot->profs);
for(i = 0; i < fd_nj_plot->notu + 1; i++) free(fd_nj_plot->labels[i]);
free(fd_nj_plot->labels);
free(fd_nj_plot->br_length_txt);
if(fd_nj_plot->tree_name != NULL) free(fd_nj_plot->tree_name);
}


/* ecriture d'un arbre non racine au format phylip, multifurcations allowed */
char *ecrit_arbre_parenth_unrooted(struct noeud *root)
{
struct noeud *t1, *t2, *t3;
char *p1, *p2, *p3, *p, *bootstrap;
int l;
float l1, l2, l3;	

t1 = root->v1->v1; t2 = root->v1->v2; t3 = root->v2;
l1 = root->v1->l1; l2 = root->v1->l2; 
if(t1 == NULL) {
	t1 = root->v2->v1; t2 = root->v2->v2; t3 = root->v1;
	l1 = root->v2->l1; l2 = root->v2->l2; 
	}
l3 = root->l1 + root->l2;
bootstrap = get_br_label(root->v1, root->v2);
if(bootstrap == NULL) bootstrap = "";

p1 = ecrit_arbre_parenth(t1);  *(p1 + strlen(p1) - 1) = 0;
p2 = ecrit_arbre_parenth(t2);  *(p2 + strlen(p2) - 1) = 0;
p3 = ecrit_arbre_parenth(t3);  *(p3 + strlen(p3) - 1) = 0;
l = strlen(p1) + strlen(p2)+ strlen(p3) + 150;
p = (char *)check_alloc(l, 1);
if(fd_nj_plot->has_br_length) sprintf(p, "(%s:%.5f,%s:%.5f,%s%s:%.5f);", p1, l1, p2, l2, p3, bootstrap, l3);
else sprintf(p, "(%s,%s,%s%s);", p1, p2, p3, bootstrap);
free(p1); free(p2); free(p3);
return p;
}


/* ecriture d'un arbre racine au format phylip, multifurcations allowed */
char *ecrit_arbre_parenth(struct noeud *root)
{
char *arbre, *fin, *p;
int l, maxarbre = 2 * fd_nj_plot->long_arbre_parenth + 1000;

arbre=check_alloc( maxarbre+20,1);
fin=recur_ecrit_arbre(root,arbre,arbre+maxarbre-1);
/* ecriture du dernier label interne */
if( fin != NULL && (p=get_br_label(root->v1,root->v2)) != NULL ) {
	l=strlen(p);
	if(fin+l>=arbre+maxarbre) fin= NULL;
	else	{
		memcpy(fin+1,p,l);
		fin+=l;
		}
	}
if(fin == NULL) {
	free(arbre);
	return NULL;
	}
strcpy(fin+1,";");
return arbre;
}


char *recur_ecrit_arbre(struct noeud *centre, char *arbre, char *finarbre)
{
int num, l;
char *p, *q;

if(centre->v1==NULL && centre->v2==NULL) {
	num=0;
	while( *(fd_nj_plot->tabtax+num) != centre) num++;
	l=strlen(fd_nj_plot->tabnames[num]);
	if(arbre+l>=finarbre) return NULL;	
	memcpy(arbre,fd_nj_plot->tabnames[num],l);
	arbre += l-1;
	}
else	{
	*arbre='(';
	p = arbre;
	arbre=recur_ecrit_arbre(centre->v1,arbre+1,finarbre);
	if(arbre==NULL) return NULL;
	if(fd_nj_plot->has_br_length) {
		if(arbre+10>=finarbre) return NULL;
		sprintf(++arbre,":%.5f",centre->l1);
		while(*arbre!=0) arbre++;
		}
	else arbre++;
	*arbre=',';
	arbre=recur_ecrit_arbre(centre->v2,arbre+1,finarbre);
	if(arbre==NULL) return NULL;
	if(fd_nj_plot->has_br_length) {
		if(arbre+10>=finarbre) return NULL;
		sprintf(++arbre,":%.5f",centre->l2);
		while(*arbre!=0) arbre++;
		}
	else arbre++;
	*arbre=')';
	/* ecriture des fd_nj_plot->labels internes */
	if( (q=get_br_label(centre,centre->v3)) != NULL  && (fd_nj_plot->has_br_length || (centre->l3 != 0) ) ) {
		l=strlen(q);
		if(arbre+l>=finarbre) return NULL;
		memcpy(arbre+1,q,l);
		arbre+=l;
		}
	else if(centre->v3 != NULL && (!fd_nj_plot->has_br_length) && (centre->l3 == 0) ) {//multibranches processed here
		memmove(p, p + 1, arbre - p);
		arbre -= 2;
		}
	}
return arbre;
}


void removeroot(void)
{
struct noeud *p1, *p2;
p1=fd_nj_plot->racine->v1;
p2=fd_nj_plot->racine->v2;
if(p1->v1 == fd_nj_plot->racine )
	{p1->v1 = p2; p1->l1 = fd_nj_plot->root_br_l;}
else if (p1->v2 == fd_nj_plot->racine)
	{p1->v2 = p2; p1->l2 = fd_nj_plot->root_br_l;}
else
	{p1->v3 = p2; p1->l3 = fd_nj_plot->root_br_l;}
if(p2->v1 == fd_nj_plot->racine )
	{p2->v1 = p1; p2->l1 = fd_nj_plot->root_br_l;}
else if (p2->v2 == fd_nj_plot->racine)
	{p2->v2 = p1; p2->l2 = fd_nj_plot->root_br_l;}
else
	{p2->v3 = p1; p2->l3 = fd_nj_plot->root_br_l;}
fd_nj_plot->rooted = 0;
}


int calc_n_desc(struct noeud *pere)
{
if(pere->v1 == NULL) return 1;
else return calc_n_desc(pere->v1) + calc_n_desc(pere->v2);
}


void postscript_plot(void)
{
double currx, curry;
int page, erreur = FALSE, offset;
int paper_h, sheets_h, recouvre, font_size;
const int a4_h = 750, letter_h = 700;
time_t heure;
float screen_w, ps_w;
static int first = TRUE;
const char alphabet[]="ABCDEFGHIJKLMNOPQRSTUVWXYZ";
/* partie dessinable 500(ps_width) x a4_h decalee de 50 x 50 par rapport a coin inferieur gauche du papier
cadre exterieur de 10 sur les 4 cotes, correspond aussi a zone de clip
*/

if(fd_nj_plot->notu == 0) return;
if(paper_choice == A4) {
	paper_h = a4_h;
	}
else {
	paper_h = letter_h;
	}
if(ps_height != 0) paper_h = ps_height;

#ifndef NO_GUI
Nlm_WatchCursor();
#endif
file_plot = TRUE;
font_size = list_font_size[fd_nj_plot->font_size_rank];

if(first) {
	first = FALSE;
#ifndef NO_GUI
#ifdef WIN_MSWIN
	Nlm_SelectFont( Nlm_ParseFont( "Courier New,12" ) );
#else
	Nlm_SelectFont( Nlm_ParseFont( "Courier,12" ) );
#endif
	screen_w = Nlm_TextWidth(alphabet, 26);/*taille sur ecran*/
	ps_w = 12 * 0.60 * 26; /*taille en postscript pour courier,12*/
	postscript_ratio = ps_w / screen_w ;
	Nlm_SelectFont(fd_nj_plot->current_font);/*retour a police courante pour ecran */
#endif
	}

/* recouvre = paper_h / 30; */
recouvre = 0;
sheets_h = (paper_h - recouvre) * page_count + recouvre;
offset = paper_h - recouvre;
physx = ps_width;
physy = sheets_h;

plotfile=fopen(plotfilename,"w");
if(plotfile == NULL) {
	char tmp[300];
	sprintf(tmp, "Error opening %s for writing\n", plotfilename);
#ifndef NO_GUI
	err_message(tmp);
#else
	fputs(tmp, stderr);
#endif
	exit(1);
	}
fprintf(plotfile,"%%!PS-Adobe-1.0\n"
	"%%%%DocumentFonts: Times-Roman %s\n"
	"%%%%Creator: njplot\n"
	"%%%%Pages: %d\n"
	"%%%%Title: %s\n"
	"%%%%BeginFeature: *PageSize \n"
 	"%s\n"
	"%%%%EndFeature\n"
	"%%%%EndComments\n",
	current_ps_font, page_count, fd_nj_plot->tree_name, paper_choice == A4 ? "a4" : "letter" );
fprintf(plotfile,
	"/setpacking where {true setpacking} if\n"
	"1 setlinecap 1 setlinejoin 1 setlinewidth 0 setgray\n"
	"/basefont /%s findfont %d scalefont def\n",
	current_ps_font, font_size);
fprintf(plotfile,"/titlefont /Times-Roman findfont 12 scalefont def\n");
fprintf(plotfile, "/setclip {40 40 moveto %d 40 lineto %d %d lineto 40 %d "
	"lineto closepath clip newpath} def\n", 
	ps_width + 60, ps_width + 60, paper_h+60, paper_h+60);
time(&heure);
fprintf(plotfile, "/title {titlefont setfont\n"
	"40 %d moveto (%s   %s) show (  Page ) show show ( of %d) show\n"
	"} def\n", 
	paper_h+65, extract_filename(fd_nj_plot->tree_name), ctime(&heure), page_count);
fprintf(plotfile,"%%%%EndProlog\n");

fd_nj_plot->totnoms = fd_nj_plot->tottraits = fd_nj_plot->totpoints = -1;
end_br_length = fd_nj_plot->br_length_txt;
currx = 0.;
maxx = 0.;
nexty = -fd_nj_plot->deltay;
if(fd_nj_plot->subtree_center == NULL)
	mem_plot(NULL, fd_nj_plot->racine, currx, &curry);
else
	mem_plot(fd_nj_plot->subtree_ascend, fd_nj_plot->subtree_center, currx, &curry);
for(page = 0; page < page_count ; page++) {
	fprintf(plotfile,"%%%%Page: ? %d\n", page+1);
	if( ! no_title) fprintf(plotfile,"(%d) title ", page+1);
	fprintf(plotfile,"setclip\n");
	fprintf(plotfile,"0 %d translate\n", - (page_count - page - 1) * offset);
	fprintf(plotfile,"basefont setfont\n");
	fprintf(plotfile,"50 50 translate\n");
	fprintf(plotfile,"0.7 setgray -10 -10 moveto %d -10 lineto %d %d lineto "
		 "-10 %d lineto closepath stroke 0 setgray \n", 
		 ps_width + 10, ps_width + 10, sheets_h+10, sheets_h+10);
	do_plot();
	fprintf(plotfile,"showpage\n");
	}
fprintf(plotfile, "%%%%Trailer\n");
if(ferror(plotfile)) erreur = TRUE;
if(fclose(plotfile) != 0) erreur = TRUE;

fd_nj_plot->need_runtree = TRUE;
file_plot = FALSE;
#ifndef NO_GUI
ArrowCursor();
if(!pdf_plot_only ) {		
	char mess[250];
	if(!erreur  ) {		
		sprintf(mess, 
			"Tree plot is now in file %s in Postscript format", 
			extract_filename(plotfilename) );
		}
	else {
		sprintf(mess, "Error while writing to file %s", 
			extract_filename(plotfilename) );
		}
	err_message(mess);
	}
#endif
}


#ifdef WITH_PDF
void plot_to_pdf(void)
{
int t_font, font_size;
char *encoding;

double currx, curry;
int page, offset;
double paper_h; 
int sheets_h, recouvre;
const int a4_h = 750, letter_h = 700;
time_t heure;
/* partie dessinable 500(ps_width) x a4_h decalee de 50 x 50 par rapport a coin inferieur gauche du papier
cadre exterieur de 10 sur les 4 cotes, correspond aussi a zone de clip
*/
char bigline[500];
#ifdef WIN_MAC
    encoding	= "macroman";
#else
    encoding	= "iso8859-1";
#endif

if(fd_nj_plot->notu == 0) return;

pdf = PDF_new();
if(pdf == NULL) return;

if( PDF_begin_document(pdf, plotfilename, 0, "compatibility=1.3") == -1) {
	PDF_delete(pdf);
	char tmp[300];
	sprintf(tmp, "Error opening %s for writing\n", plotfilename);
#ifndef NO_GUI
	err_message(tmp);
#else
	fputs(tmp, stderr);
#endif
	return;
	}

PDF_TRY(pdf) {

PDF_set_info(pdf, "Title", fd_nj_plot->tree_name );
PDF_set_info(pdf, "Creator", "njplot");
t_font = PDF_load_font(pdf, "Times-Roman", 0, encoding, "");
pdf_font = PDF_load_font(pdf, current_ps_font, 0, encoding, "");
font_size = list_font_size[fd_nj_plot->font_size_rank];

if(paper_choice == A4) { paper_h = a4_h; }
else /* LETTER */ { paper_h = letter_h; }
if(ps_height != 0) paper_h = ps_height;

#ifndef NO_GUI
Nlm_WatchCursor();
#endif
pdf_plot = TRUE; file_plot = TRUE;


/* recouvre = paper_h / 30; */
recouvre = 0;
sheets_h = (paper_h - recouvre) * page_count + recouvre;
offset = paper_h - recouvre;
physx = ps_width;
physy = sheets_h;

time(&heure);

fd_nj_plot->totnoms = fd_nj_plot->tottraits = fd_nj_plot->totpoints = -1;
end_br_length = fd_nj_plot->br_length_txt;
currx = 0.;
maxx = 0.;
nexty = -fd_nj_plot->deltay;
if(fd_nj_plot->subtree_center == NULL)
	mem_plot(NULL, fd_nj_plot->racine, currx, &curry);
else
	mem_plot(fd_nj_plot->subtree_ascend, fd_nj_plot->subtree_center, currx, &curry);
for(page = 0; page < page_count ; page++) {
	PDF_begin_page_ext(pdf, 
		paper_choice == A4 ? a4_width : letter_width,
		paper_choice == A4 ? a4_height : letter_height,
		"");
	if( ! no_title) {
		PDF_setfont(pdf, t_font, 12);
		sprintf(bigline, "%s  %s  Page %d of %d", extract_filename(fd_nj_plot->tree_name), ctime(&heure), page+1, page_count);
		PDF_show_xy(pdf, bigline, 40, paper_h+65);
		}

	PDF_setlinewidth(pdf, 1);
	PDF_moveto(pdf, 40, 40); PDF_lineto(pdf, ps_width + 60, 40); PDF_lineto(pdf, ps_width + 60, paper_h+60);
	PDF_lineto(pdf, 40, paper_h+60); 
	PDF_closepath(pdf); PDF_clip(pdf);
	PDF_translate(pdf, 0, - (page_count - page - 1) * offset);
	PDF_setfont(pdf, pdf_font, font_size);
	PDF_translate(pdf, 50, 50);
	PDF_setcolor(pdf, "stroke", "gray", 0.7, 0,0,0);
	PDF_moveto(pdf, -10, -10); 
	PDF_lineto(pdf, ps_width + 10, -10); 
	PDF_lineto(pdf, ps_width + 10, sheets_h+10); 
	PDF_lineto(pdf, -10, sheets_h+10); 
	PDF_closepath_stroke(pdf); 
	PDF_setcolor(pdf, "stroke", "gray", 0, 0,0,0);
	do_plot();
	PDF_end_page_ext(pdf, "");
	}
PDF_end_document(pdf, "");
PDF_delete(pdf);
} /* end of PDF_TRY */
PDF_CATCH(pdf) {
        sprintf(bigline, "Error while writing pdf file: [%d] %s: %s",
	  	  PDF_get_errnum(pdf), PDF_get_apiname(pdf), PDF_get_errmsg(pdf) );
		  err_message(bigline);
        PDF_delete(pdf);
        return; 
}

fd_nj_plot->need_runtree = TRUE;
file_plot = FALSE;
#ifndef NO_GUI
Nlm_ArrowCursor();
if(!pdf_plot_only ) {		
	char mess[250];
	sprintf(mess, 
			"Tree plot is now in file %s in PDF format", 
			extract_filename(plotfilename) );
	err_message(mess);
	}
#endif
}
#endif

void tty_plot(void)
{
double currx, curry;
int h, erreur = FALSE;
char *p;

if(fd_nj_plot->notu == 0) return;
file_plot = TRUE;

/* physx was set in process_args */
physy = 2*(fd_nj_plot->notu + 1) + 4;
tty_page = (char **)malloc(physy*sizeof(char *));
for(h=0; h < physy; h++) {
	tty_page[h] = (char *)malloc(physx + 1);
	memset(tty_page[h], ' ', physx);
	tty_page[h][(int)physx] = 0;
	}

fd_nj_plot->totnoms = fd_nj_plot->tottraits = fd_nj_plot->totpoints = -1;
end_br_length = fd_nj_plot->br_length_txt;
currx = 0.;
maxx = 0.;
nexty = -fd_nj_plot->deltay;
if(fd_nj_plot->subtree_center == NULL)
	mem_plot(NULL, fd_nj_plot->racine, currx, &curry);
else
	mem_plot(fd_nj_plot->subtree_ascend, fd_nj_plot->subtree_center, currx, &curry);
do_plot();
for(h=0; h < physy; h++) {
	p = tty_page[h] + strlen(tty_page[h]) - 1;
	while(p > tty_page[h] && *p == ' ') *(p--) = 0;
	}

plotfile=fopen(plotfilename, "w");
for(h=1; h <= physy; h++) {
	fprintf(plotfile, "%s\n", tty_page[(int)physy - h]);
	}
if(ferror(plotfile)) erreur = TRUE;
if(fclose(plotfile) != 0) erreur = TRUE;

fd_nj_plot->need_runtree = TRUE;
file_plot = FALSE;
}


void majuscules(char *p)
{
while(*p != 0) {*p = toupper(*p); p++;}
}


#ifdef WIN_MAC
#if TARGET_RT_MAC_MACHO
#define DIR_DELIM '/'
#else
#define DIR_DELIM ':'
#endif
#elif defined(WIN_MSWIN)
#define DIR_DELIM '\\'
#else
#define DIR_DELIM '/'
#endif

char *extract_filename(char *name)
{
char *p;
while( ( p = strchr(name, DIR_DELIM) ) != NULL) name = p + 1;
return name;
}


#ifdef WIN_MOTIF
#include <ncbiport.h>


static void Nlm_GetFontData (Nlm_FonT f, Nlm_FontData * fdata)

{
  Nlm_FntPtr  fp;

  if (f != NULL && fdata != NULL) {
    fp = (Nlm_FntPtr) Nlm_HandLock (f);
    *fdata = *fp;
    Nlm_HandUnlock (f);
  }
}


void set_systemfont(char *fname)
{
		XFontStruct   *font;
		XmFontListEntry entry;
		Nlm_FontData fontdata;
		extern XmFontList Nlm_XfontList;
		
		Nlm_systemFont = Nlm_ParseFont(fname);
        Nlm_GetFontData(Nlm_systemFont, &fontdata);
        font= (XFontStruct *)fontdata.handle;
		entry = XmFontListEntryCreate(XmFONTLIST_DEFAULT_TAG, 
					XmFONT_IS_FONT, (XtPointer)font);
        Nlm_XfontList = XmFontListAppendEntry(NULL, entry);
                
        Nlm_SelectFont(Nlm_systemFont);
        Nlm_stdLineHeight=Nlm_LineHeight();
        Nlm_stdFontHeight=Nlm_FontHeight();
        Nlm_stdAscent=Nlm_Ascent();
        Nlm_stdDescent=Nlm_Descent();
        Nlm_stdLeading=Nlm_Leading();
        Nlm_stdCharWidth = Nlm_MaxCharWidth ();
}
#endif
