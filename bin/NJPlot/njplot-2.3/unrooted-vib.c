#include <vibtypes.h>
#include <vibprocs.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#if defined(WIN_MSWIN) || defined (WIN_MAC)
#define  MICRO
#endif

#ifdef WIN_MAC
#include <Quickdraw.h>
#include <Files.h>
#include <Scrap.h>
#include <Types.h>
#ifndef M_PI
#define M_PI           3.14159265358979323846  /* pi */
#endif
#if TARGET_RT_MAC_MACHO
#define DIR_DELIM '/'
#else
#define DIR_DELIM ':'
#endif
#define myWaitAndProcessNextEvent Nlm_ProcessAnEvent
extern int Nlm_textScrapFull;


#elif defined(WIN_MSWIN)
#define HDC void *
extern HDC Nlm_currentHDC;
extern HDC Nlm_GetPicWinHDC ( void );
char *get_prog_dir(void);
#define M_PI pi
#define DIR_DELIM '\\'
extern const double pi;
#define myWaitAndProcessNextEvent Nlm_ProcessAnEvent
#define tempnam _tempnam

#include <conio.h> //temporaire passage CW8
int __getche(void)
{
return getche();
}

#else  /* unix */


#include <time.h>
#define DIR_DELIM '/'
#undef Boolean
#include <X11/Intrinsic.h>
void set_systemfont(char *systemfont);
extern Nlm_CharPtr Nlm_XrmGetResource(const Nlm_Char PNTR _resource);


extern XtAppContext Nlm_appContext;
void myWaitAndProcessNextEvent (void)
{
XEvent       event;
if( XtAppPending(Nlm_appContext) & XtIMXEvent ) {
        Nlm_ProcessAnEvent();
        }
else    {
/* wait until next event and catch it */
        XtAppNextEvent (Nlm_appContext, &event);
        XtDispatchEvent (&event); 
        }
}


#endif


#ifndef M_PI
#define M_PI           3.14159265358979323846  /* pi */
#endif



/* typedefs */
typedef enum { NO_HIT, LEFT_HIT, RIGHT_HIT, TOP_HIT, BOTTOM_HIT } hit;


enum {times = 1, helvetica, courier};


enum { font_tiny = 1, font_small,
        font_normal, font_medium, font_large, font_bold, font_italic };
        
typedef enum { A4 = 1, LETTER } 
        paper_item;


typedef struct {
        int x, y;
        } point;


typedef struct _cp_point {
        double x, y, r, angle;
        } cp_point;


typedef struct _bignoeud {
        struct _bignoeud *v1, *v2, *v3;
        double l1, l2, l3;
        cp_point position;
        char *nom;
        } bignoeud;


struct noeud {
        double l1,l2,l3;
        struct noeud *v1,*v2,*v3;
        char *nom;
        };


typedef struct _branche {
        cp_point debut, fin;
        int color;
        char *nom;
        } branche;


typedef struct _plot_data {
        int x, y, w, h;
        int ps_plot;
        int lstyle, lcol, lsize;
        int comp_phys_bounds;
        } plot_data;
        
typedef struct {
        Nlm_WindoW unrooted;
        Nlm_IteM open_button;
        Nlm_MenU tree_font_menu;
        Nlm_MenU paper_menu;
        Nlm_IteM save_plot_button;
		Nlm_IteM AgainItem;
        Nlm_ChoicE choix_font, choix_taille;
        Nlm_IteM bold_item, italic_item;
        Nlm_IteM exit_button;
        Nlm_PaneL tree_plot;
        char *tree_name;
} FD_unrooted;



/* prototypes of included functions */
int direct_ps_plot(void);
FD_unrooted *create_win_unrooted(void);
void open_callback(Nlm_IteM ob);
void win_resize_proc(Nlm_WindoW);
void prepare_fonts(void);
void change_panel_size(Nlm_PaneL panel, int change_x, int change_y);
void err_message(char *text);
void place_midpoint_root(void);
char *extract_filename(char *fname);
void save_plot_callback(Nlm_IteM ob);
void font_callback(Nlm_ChoicE choix);
void paper_callback(Nlm_ChoicE ob);
void exit_callback(Nlm_IteM ob);
void tree_draw_proc(Nlm_PaneL panel);
void set_tree_w_title(Nlm_WindoW w, char *file, int tips);
void init_tree(char *titre);
void debut_arbre(void);
void my_draw_line(int dx, int dy, int fx, int fy, int col, int ps, int h);
void my_draw_text(int align, int x, int y, int col, int size, int style, 
        char *chaine, int ps, int h, int ascent);
hit draw_line(branche *br, plot_data *ob, int doit, int char_width,
        int char_height, int descend, int ascent);
hit draw_name_angle(cp_point *phys_pos, char *nom, plot_data *ob, double angle,
        int doit, int char_width, int char_height, int descend, int color,
        int ascent);
hit draw_letter(point *position, char lettre, plot_data *ob, int doit, 
        int char_width, int color, int ascent);
void draw_scale(plot_data *ob, int char_height, int ascent);
int my_get_char_height(int *ascend, int *descend);
int my_get_string_width(char *nom);
char *check_alloc(int nbrelt, int sizelt);
void majuscules(char *p);
void search_callback(Nlm_IteM item);
void reset_callback(Nlm_IteM item);
void string_callback(Nlm_ButtoN item);
void process_keys(Nlm_Char key);
void clear_tree(Nlm_IteM unused);
void paste_tree(Nlm_IteM unused);
#ifdef MICRO
void copy_callback(Nlm_IteM);
int print_title(int x, int y, char *text);
void print_plot(Nlm_IteM);
#endif
#ifdef WIN_MAC
int crefpict(char *, PicHandle );
extern void add_apropos(char *progname, Nlm_ActnProc my_show_apropos);
extern void show_apropos_unrooted(Nlm_IteM item);
extern char *mac_fname_to_roman(char *);
#else
void process_args(int *argc, char *argv[]);
void remove_arg(int target, int *argc, char *argv[]);
#endif

#ifdef WIN_MAC
extern int MG_GetInputFName(char *fname, int maxl);
extern int MG_GetOutputFName(char *fname, int maxl, char *dfault);
extern QDPictRef MyPictToQDPict(PicHandle mypicture);
extern int MyQDPictToPDFfile (QDPictRef picture, char *fname);
#elif defined(WIN_MSWIN)
extern int MG_GetInputFName(char *fname, int maxl);
extern int MG_GetOutputFName(char *fname, int maxl, char *dfault);
#else 
#define MG_GetOutputFName Nlm_GetOutputFileName
#define MG_GetInputFName(a,b) Nlm_GetInputFileName(a, b, NULL, NULL)
#endif



/* external functions */
double arrondi_echelle(double x);
double calc_echelle(int larg);
bignoeud *cre_new_tree(struct noeud *debut, struct noeud *parent, 
                        bignoeud *bigparent);
void remove_big_root(bignoeud *bigracine);
double calc_dist_centre_feuilles(bignoeud *debut, bignoeud *parent);
double proc_null_neg_branches(bignoeud *debut, bignoeud *parent);
int set_angles_noeuds(bignoeud *debut, bignoeud *parent, double delta,
                                double *p_current_angle, double rayon);
void calc_cartesienne( cp_point *p);
void calc_polaire( cp_point *p);
cp_point calc_point_direction( cp_point *depart, cp_point *direction, 
        double longueur);
branche *calc_position_noeuds(bignoeud *debut, bignoeud *parent,
        branche *curr_branche);
void mem_line(cp_point *debut, cp_point *fin, bignoeud *noeud_term, 
        branche *br);
void draw_tree(plot_data *ob, branche *branches, int comp_phys_bounds);
void log_to_phys(cp_point *log_pos, cp_point *phys_pos);
double length_log_phys(double p);
double length_phys_log(double p);
char *preptree(char *fname);
void free_tree(void);
int myrint(double val);



/* globals */
FD_unrooted *fd_unrooted;
Nlm_FonT current_font;
char current_ps_font[40];
#ifdef WIN_MSWIN
int current_font_size_rank = font_small;
#else
int current_font_size_rank = font_normal;
#endif
int list_font_size[] = { 0, 8, 10, 12, 14, 18 };
#ifdef WIN_MSWIN
char list_font_name[4][40] = {"", "Times New Roman", "Lucida Sans", 
	"Courier New"};
#else
char list_font_name[4][15] = {"", "Times", "Helvetica", "Courier"};
#endif
paper_item paper_choice = A4;
char fname[200];
plot_data plot_aux_data;
branche *branches;
double log_min_x, log_min_y, log_max_x, log_max_y, mini_br_length;
int phys_min_x, phys_min_y, phys_max_x, phys_max_y;
int doing_copy = FALSE;
int doing_print = FALSE; Nlm_RecT print_rect;
int ps_plot_only = FALSE;
FILE *plotfile;
#ifdef WIN_MAC
PicHandle       mypicture;
int             margin;
#endif



/* external variables */
extern int notu, has_br_length;
extern struct noeud *racine;


int argc;
char **argv;




Nlm_Int2 Nlm_Main(void)
{
Nlm_RecT *p_rect;


*fname = 0;
#ifndef WIN_MAC
argc = Nlm_GetArgc();
argv = Nlm_GetArgv();
process_args(&argc, argv);
if(argc>=2) {
        strcpy(fname,argv[1]);
        }
else    {
        *fname=0;
        }
#endif
   if(ps_plot_only) {
        if ( *fname != 0 ) exit( direct_ps_plot() );
        else exit(1);
        }
   Nlm_WatchCursor();
#ifdef WIN_MOTIF
	{char *font;
	font = Nlm_XrmGetResource("systemfont");
	if(font != NULL) set_systemfont(font);
	}
#endif
   fd_unrooted = create_win_unrooted();
   prepare_fonts();
   fd_unrooted->tree_name = fname;
   /* show the first form */
   if(*fname != 0) debut_arbre();
   set_tree_w_title(fd_unrooted->unrooted, (*fname == 0 ? "unrooted" : fname), notu );
   Nlm_Show(fd_unrooted->unrooted);
   p_rect = (Nlm_RecT *)Nlm_GetWindowExtra(fd_unrooted->unrooted);
   Nlm_ObjectRect(fd_unrooted->unrooted, p_rect);
   Nlm_ArrowCursor();
   Nlm_ProcessEvents();
   return 0;
}



int direct_ps_plot(void)
{
        /* valeurs initiales des extr physiques des lignes du graphique */
        phys_max_x = 500;
        phys_max_y = 750;
        phys_min_x = 0;
        phys_min_y = 0;
debut_arbre();
strcpy(current_ps_font, "Times-Roman");
plot_aux_data.comp_phys_bounds = TRUE;
save_plot_callback(NULL);
return 0;
}



FD_unrooted *create_win_unrooted(void)
{
  Nlm_WindoW win;
  Nlm_MenU menu;
  Nlm_ChoicE choix_paper;
  Nlm_WindoW place_menu;
  FD_unrooted *fdui = (FD_unrooted *) calloc(1, sizeof(FD_unrooted));
  static Nlm_RecT winrect;
  int width, height;
  
#ifdef WIN_MAC
add_apropos("unrooted", show_apropos_unrooted);
#endif
Nlm_KeyboardView(process_keys);
fdui->unrooted = win = Nlm_DocumentWindow(-50, -33, -3, -3, "Unrooted", 
(Nlm_WndActnProc) exit_callback, win_resize_proc);
#ifdef WIN_MAC
place_menu = NULL;
#else
place_menu = win;
#endif
#ifdef WIN_MAC
#define OPEN_LABEL      "Open        O"
#define SAVE_PLOT_LABEL "Save Plot S"
#define PRINT_LABEL     "Print        P"
#define FIND_LABEL      "Find     F"
#define AGAIN_LABEL     "Again  A"
#define RESET_LABEL     "Reset  R"
#define COPY_LABEL      "Copy   C"
#define PASTE_LABEL     "Paste  V"
#else
#define OPEN_LABEL      "Open       ^O"
#define SAVE_PLOT_LABEL "Save Plot ^S"
#define PRINT_LABEL     "Print     ^P"
#define QUIT_LABEL      "Quit       ^Q"
#define FIND_LABEL      "Find       ^F"
#define AGAIN_LABEL     "Again     ^A"
#define RESET_LABEL     "Reset     ^R"
#define COPY_LABEL      "Copy      ^C"
#define PASTE_LABEL     "Paste     ^V"
#endif
  menu = Nlm_PulldownMenu(place_menu, "File");
  fdui->open_button = Nlm_CommandItem(menu, OPEN_LABEL, open_callback);
  fdui->save_plot_button = Nlm_CommandItem(menu, SAVE_PLOT_LABEL, save_plot_callback);
#ifdef MICRO
  Nlm_CommandItem(menu, PRINT_LABEL, print_plot);
#endif
#ifndef WIN_MAC
  Nlm_CommandItem(menu, QUIT_LABEL, exit_callback);
#endif
  Nlm_Advance(place_menu);
  menu = Nlm_PulldownMenu(place_menu, "Edit"); 
#if defined(WIN_MAC) || defined(WIN_MSWIN)
#endif
  Nlm_CommandItem(menu, "Clear", clear_tree);
#ifdef MICRO
  Nlm_CommandItem(menu, COPY_LABEL, copy_callback);
#endif
  Nlm_CommandItem(menu, PASTE_LABEL, paste_tree);
  Nlm_SeparatorItem(menu);
  Nlm_CommandItem(menu, FIND_LABEL, search_callback);
  fdui->AgainItem = Nlm_CommandItem(menu, AGAIN_LABEL, search_callback);
  Nlm_CommandItem(menu, RESET_LABEL, reset_callback);
  Nlm_Disable(fdui->AgainItem);
  Nlm_Advance(place_menu);
  fdui->tree_font_menu = Nlm_PulldownMenu(place_menu, "Font"); 
  fdui->choix_font = Nlm_ChoiceGroup(fdui->tree_font_menu, font_callback);
  Nlm_ChoiceItem(fdui->choix_font, list_font_name[1]);
  Nlm_ChoiceItem(fdui->choix_font, list_font_name[2]);
  Nlm_ChoiceItem(fdui->choix_font, list_font_name[3]);
  Nlm_SetValue(fdui->choix_font, times);
  Nlm_SeparatorItem(fdui->tree_font_menu);
  fdui->choix_taille = Nlm_ChoiceGroup(fdui->tree_font_menu, font_callback);
  Nlm_ChoiceItem(fdui->choix_taille, "8");
  Nlm_ChoiceItem(fdui->choix_taille, "10");
  Nlm_ChoiceItem(fdui->choix_taille, "12");
  Nlm_ChoiceItem(fdui->choix_taille, "14");
  Nlm_ChoiceItem(fdui->choix_taille, "18");
  Nlm_SetValue(fdui->choix_taille, current_font_size_rank);
  Nlm_SeparatorItem(fdui->tree_font_menu);
  fdui->bold_item = Nlm_StatusItem(fdui->tree_font_menu, "Bold", 
        (Nlm_ItmActnProc) font_callback);
  fdui->italic_item = Nlm_StatusItem(fdui->tree_font_menu, "Italic", 
        (Nlm_ItmActnProc) font_callback);
  Nlm_SetStatus(fdui->bold_item, FALSE);
  Nlm_SetStatus(fdui->italic_item, FALSE);
  
  Nlm_Advance(place_menu);
  fdui->paper_menu = Nlm_PulldownMenu(place_menu, "Paper"); 
  choix_paper = Nlm_ChoiceGroup(fdui->paper_menu, paper_callback);
  Nlm_ChoiceItem(choix_paper, "A4");
  Nlm_ChoiceItem(choix_paper, "US letter");
  if(paper_choice == A4)
        Nlm_SetValue(choix_paper, 1);
  else
        Nlm_SetValue(choix_paper, 2);


  Nlm_Break(win);

#ifdef WIN_MOTIF 
/* bug in LessTif : resize does not work if nothing except panel */
{ Nlm_Handle obj; Nlm_PoinT position;
Nlm_GetNextPosition(win, &position);
obj = Nlm_PushButton(win,"void",NULL); Nlm_Hide(obj);
Nlm_SetNextPosition(win, position);
}
#endif

#if defined(unix) || defined(WIN_MSWIN)
        width = height = 510;
#else
        width = height = 450;
#endif
fdui->tree_plot = Nlm_AutonomousPanel(win, width, height, tree_draw_proc, NULL, 
                NULL, 0, NULL, NULL);
        winrect.left = -1;
        Nlm_SetWindowExtra(win, &winrect, NULL);
  return fdui;
}




void win_resize_proc(Nlm_WindoW win)
{
Nlm_RecT *old_rect, rect;
int change_x, change_y, old_height, new_height;


old_rect = (Nlm_RecT *)Nlm_GetWindowExtra(win);
if(old_rect->left == -1) return;
Nlm_ObjectRect(win, &rect);
old_height = old_rect->bottom - old_rect->top;
new_height = rect.bottom - rect.top;
change_y = (new_height - old_height);
old_height = old_rect->right - old_rect->left;
new_height = rect.right - rect.left;
change_x = (new_height - old_height);
if(change_x != 0 || change_y != 0) {
        change_panel_size(fd_unrooted->tree_plot, change_x, change_y);
        *old_rect = rect;
        }
}



void change_panel_size(Nlm_PaneL panel, int change_x, int change_y)
{
Nlm_RecT r, r_group;
Nlm_GrouP group;
Nlm_Boolean in_group;
static int recursive = FALSE;


Nlm_GetPosition(panel, &r);
#ifdef WIN_MOTIF
if(recursive) return;
recursive = TRUE;
group = Nlm_Parent(panel);
in_group = (group != Nlm_ParentWindow(panel));
if(in_group) {
        Nlm_GetPosition(group, &r_group);
        Nlm_Hide(group);
        }
else
        Nlm_Hide(panel);
#endif
r.bottom += change_y;
r.right += change_x;
Nlm_SetPosition(panel, &r);
#ifdef WIN_MOTIF
if(in_group) {
        r_group.bottom += change_y;
        r_group.right += change_x;
        Nlm_SetPosition( group, &r_group);
        Nlm_Show(group);
        }
else    Nlm_Show(panel);
recursive = FALSE;
#endif
}



void font_callback(Nlm_ChoicE choix)
{
prepare_fonts();
plot_aux_data.comp_phys_bounds = TRUE;
tree_draw_proc(fd_unrooted->tree_plot);
}



void prepare_fonts(void)
{
Nlm_Boolean use_bold, use_italic;
Nlm_Int2 font_num;
char font_full_name[40];


/* prepare vibrant font */
use_bold = Nlm_GetStatus(fd_unrooted->bold_item);
use_italic = Nlm_GetStatus(fd_unrooted->italic_item);
font_num = Nlm_GetValue(fd_unrooted->choix_font);
current_font_size_rank = Nlm_GetValue(fd_unrooted->choix_taille);
sprintf(font_full_name, "%s,%d", list_font_name[font_num], 
        list_font_size[current_font_size_rank]);
if(use_bold || use_italic) strcat(font_full_name, ",");
if(use_bold) strcat(font_full_name, "b");
if(use_italic) strcat(font_full_name, "i");
current_font = Nlm_ParseFont(font_full_name);


/* prepare postscript name of font */
strcpy(current_ps_font, list_font_name[font_num]);
current_ps_font[0] = toupper(current_ps_font[0]);
if( use_bold )
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



#ifndef WIN_MAC
void remove_arg(int target, int *argc, char *argv[])
{
int num;
for(num = target; num < *argc - 1; num++) 
        argv[num] = argv[num+1];
(*argc)--;
}


void process_args(int *argc, char *argv[])
{
int num, taille;


for(num = 1; num < *argc; num++) {
        if( strncmp(argv[num], "-h", 2) ==0 ){
                fprintf(stderr,"Usage: unrooted [options] [tree_file_name]\n"
"where options are:\n"
"-h             print out this message\n"
"-psonly        no window interface, just write the postscript tree file\n"
"-us            postcript tree file prepared for US Letter paper size\n"
"-size n        font size n used for taxon names\n"
"\n"
"and where tree_file_name is the name of a tree file using nested parentheses\n");
                exit(0);
                }
        }


for(num = 1; num < *argc; num++) {
        if( strncmp(argv[num], "-us", 3) == 0) {
                paper_choice = LETTER;
                remove_arg(num, argc, argv);
                break;
                }
        }


for(num = 1; num < *argc; num++) {
        if( strncmp(argv[num], "-psonly", 7) == 0) {
                ps_plot_only = TRUE;
                remove_arg(num, argc, argv);
                break;
                }
        }


for(num = 1; num < *argc; num++) {
        if( strncmp(argv[num], "-size", 5) == 0) {
                taille = 12;
                if(num + 1 < *argc) {
                        sscanf(argv[num + 1], "%d", &taille);
                        remove_arg(num + 1, argc, argv);
                        }
                remove_arg(num, argc, argv);
                for(num = 0; num < 5; num++) 
                        if(taille <= list_font_size[num]) break;
                current_font_size_rank = num;
                break;
                }
        }


}
#endif



void open_callback(Nlm_IteM ob)
{
Nlm_Boolean reponse;
char *p;

reponse = MG_GetInputFName(fname, sizeof(fname) );
if(!reponse) {
	tree_draw_proc(fd_unrooted->tree_plot);
        return;
        }
Nlm_WatchCursor();
if(notu != 0) free_tree();
notu = 0;
debut_arbre();
plot_aux_data.comp_phys_bounds = TRUE;
set_tree_w_title(fd_unrooted->unrooted, fname, notu);
tree_draw_proc(fd_unrooted->tree_plot);
}


void set_tree_w_title(Nlm_WindoW w, char *file, int tips)
{
char tmp[200];

file = extract_filename(file);
#ifdef WIN_MAC
file = mac_fname_to_roman(file);
#endif
if(tips != 0) sprintf(tmp, "%s (%d tips)", file, tips + 1);
Nlm_SetTitle(w, (tips == 0 ? "unrooted" : tmp ) );
}


void paper_callback(Nlm_ChoicE ob)
{
paper_choice = (paper_item) Nlm_GetValue(ob);
}



void exit_callback(Nlm_IteM ob)
{
exit(0);
}



void tree_draw_proc(Nlm_PaneL panel)
{
Nlm_RecT ob_rect;
Nlm_Int2 width, height;
static Nlm_Int2 previous_h = 0, previous_w = 0;

if(doing_print) 
	ob_rect = print_rect;
else	{
	Nlm_Select(panel);
	Nlm_ObjectRect(panel, &ob_rect);
	}
#ifdef WIN_MOTIF
ob_rect.bottom -= 3;
#endif
Nlm_Black();
Nlm_EraseRect(&ob_rect);
Nlm_FrameRect(&ob_rect);
if(notu == 0) return;
width = ob_rect.right - ob_rect.left + 1;
height = ob_rect.bottom - ob_rect.top + 1;
if( width != previous_w || height != previous_h ) {
/* si on a redimensionne la fenetre ou la 1ere fois */
                previous_w = width; previous_h = height;
                plot_aux_data.comp_phys_bounds = TRUE;
                plot_aux_data.x = ob_rect.left;
                plot_aux_data.y = ob_rect.top;
                plot_aux_data.w = width;
                plot_aux_data.h = height;
                tree_draw_proc(panel);
                return;
                }
Nlm_SelectFont(current_font);
if(!doing_print) Nlm_ClipRect(&ob_rect);
if( plot_aux_data.comp_phys_bounds ) {
        /* valeurs initiales des extr physiques des lignes du graphique */
        phys_max_x = ob_rect.right;
        phys_max_y = ob_rect.bottom;
        phys_min_x = ob_rect.left;
        phys_min_y = ob_rect.top;
        }
plot_aux_data.x = ob_rect.left;
plot_aux_data.y = ob_rect.top;
plot_aux_data.w = width;
plot_aux_data.h = height;
Nlm_WatchCursor();
draw_tree(&plot_aux_data, branches, plot_aux_data.comp_phys_bounds);
Nlm_ArrowCursor();
if(!doing_print) Nlm_ResetClip();
plot_aux_data.comp_phys_bounds = FALSE;
}



/* fenetre message d'alerte.
Operation Select dans la boucle des evenements qui ramene toujours la fenetre
au premier plan. Lent, mais interessant.
*/
static void alert_ok_action(Nlm_ButtoN bouton)
{
Nlm_Boolean *alert_done;
alert_done = (Nlm_Boolean *)Nlm_GetWindowExtra(Nlm_Parent(bouton));
*alert_done = TRUE;
}



void err_message(char *texte)
{
Nlm_WindoW alwin;
Nlm_Boolean alert_done;


#ifdef WIN_MAC
alwin = Nlm_FixedWindow(-50,-50,-5,-5,"MESSAGE",NULL); // Bug: ModalWindow crashes when return key hit
#else
alwin = Nlm_ModalWindow(-50,-50,-5,-5,NULL);
#endif
Nlm_StaticPrompt(alwin,texte,0,0,Nlm_programFont,'l');
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



char *extract_filename(char *name)
{
char *p;
while( ( p = strchr(name, DIR_DELIM) ) != NULL) name = p + 1;
return name;
}



#ifdef WIN_MAC
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


void copy_callback(Nlm_IteM unused)
{
extern void MyCopyPictToClipboard (PicHandle mypicture);

if(notu == 0) return;
doing_copy = TRUE;
save_plot_callback(NULL);

MyCopyPictToClipboard (mypicture );
KillPicture(mypicture);
doing_copy = FALSE;
}
#endif

#ifdef WIN_MSWIN

void copy_callback(Nlm_IteM unused)
{
Nlm_RecT r;
Nlm_WindoW w;

if(notu == 0) return;
/* dans la picture, les calculs de taille de caractere ne marchent pas!
il faut les faire dans le DC de la fenetre
*/
doing_copy = TRUE;
Nlm_GetPosition( fd_unrooted->tree_plot , &r);
w = Nlm_StartPicture(&r);
tree_draw_proc(fd_unrooted->tree_plot);
Nlm_EndPicture(w);
doing_copy = FALSE;
}
#endif


void paste_tree(Nlm_IteM unused)
{
char *buff;
FILE *tmp;
char *tmpfname;
char pasted_title[] = "pasted tree";

if(notu != 0) {
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
tmpfname = tempnam(NULL, "unrootedtemp_");
tmp = fopen(tmpfname,"w");
fwrite(buff, 1, strlen(buff), tmp);
fclose(tmp);
strcpy(fname, tmpfname);
debut_arbre();
remove(tmpfname);
if(notu > 0) {
	set_tree_w_title(fd_unrooted->unrooted, pasted_title, notu);
	strcpy(fname, pasted_title);
	plot_aux_data.comp_phys_bounds = TRUE;
	tree_draw_proc(fd_unrooted->tree_plot);
	}
else	fname[0] = 0;
}




void save_plot_callback(Nlm_IteM ob)
{
const int a4_h = 750, letter_h = 700;
int paper_h, font_size;
char *sp, *p1, message[100];
char plotfilename[200];
int myw, myh, myx, myy;
#ifdef WIN_MAC
int macfont;
Style mystyle=normal;
Rect myrect;
#else
time_t heure;
#endif


if(notu==0)return;
myw = plot_aux_data.w; myh = plot_aux_data.h; 
myx = plot_aux_data.x; myy = plot_aux_data.y;
if(!doing_copy) {
        strcpy(plotfilename,fname);
        sp = strchr(plotfilename,'.'); 
        if(sp == NULL) sp = plotfilename + strlen(plotfilename);
#ifdef WIN_MAC
        strcpy(sp,".pdf");
#else
        strcpy(sp,".ps");
#endif
        p1 = extract_filename(plotfilename);
        if(!ps_plot_only) {
                if( !MG_GetOutputFName(plotfilename, 
                        sizeof(plotfilename), p1) ) return;
                }
        }
#ifdef WIN_MAC
        margin = 30;
        myrect.top=10;
        myrect.left=10;
        myrect.bottom=760;
        myrect.right=500;
        ClipRect(&myrect);
        mypicture = OpenPicture(&myrect);
        PenNormal();
        macfont = Nlm_GetValue(fd_unrooted->choix_font);
        if(macfont == courier)
                macfont=22;
        else if(macfont == helvetica)
                macfont=21;
        else if(macfont == times)
                macfont=20;
        else
                macfont=22;
        TextFont(macfont);
        TextSize( list_font_size[current_font_size_rank] );
    if( Nlm_GetStatus(fd_unrooted->bold_item) )
        mystyle += bold;
    if( Nlm_GetStatus(fd_unrooted->italic_item) )
        mystyle += italic;
    TextFace(mystyle);
        plot_aux_data.w = myrect.right-myrect.left - 2 * margin;
        plot_aux_data.h = myrect.bottom-myrect.top - 2 * margin;
#else
        font_size = list_font_size[current_font_size_rank];
        plotfile = fopen(plotfilename,"w");
        if(plotfile == NULL) {
                return;
                }
        if(paper_choice == A4)
                paper_h = a4_h;
        else
                paper_h = letter_h;
plot_aux_data.h = paper_h;
plot_aux_data.w = 500; 
plot_aux_data.x = plot_aux_data.y = 0;
fprintf(plotfile,"%%!\n1 setlinecap 1 setlinejoin 1 setlinewidth 0 setgray\n");
fprintf(plotfile,"/basefont /%s findfont %d scalefont def\n",
        current_ps_font, font_size);
fprintf(plotfile,"/titlefont /Times-Roman findfont 12 scalefont def\n");
fprintf(plotfile,"50 50 translate\n");
fprintf(plotfile,"-10 -10 moveto 510 -10 lineto 510 %d lineto -10 %d lineto \
-10 -10 lineto stroke\n", paper_h+10, paper_h+10);
time(&heure);
fprintf(plotfile,"titlefont setfont\n");
fprintf(plotfile, "0 %d moveto (%s     %s) show\n", paper_h+15, fname, ctime(&heure));
fprintf(plotfile,"basefont setfont\n");
#endif
plot_aux_data.ps_plot = TRUE;
Nlm_WatchCursor();
draw_tree(&plot_aux_data, branches, TRUE);
Nlm_ArrowCursor();
plot_aux_data.w = myw; plot_aux_data.h = myh;
plot_aux_data.x = myx; plot_aux_data.y = myy;
plot_aux_data.ps_plot = FALSE;


#ifdef WIN_MAC
        ClosePicture();
        if( !doing_copy ) 
                if( crefpict(plotfilename, mypicture) ) {
                        sprintf(message, "Problem writing PDF file %s", plotfilename);
                        err_message(message);
                        }
#else
fprintf(plotfile,"showpage\n");
fclose(plotfile);
if(ps_plot_only) return;
sprintf(message, "The tree is now in file %s "
                "in Postscript format", plotfilename);
err_message(message);
#endif
plot_aux_data.comp_phys_bounds = TRUE;
if( !doing_copy ) tree_draw_proc(fd_unrooted->tree_plot);
}


void create_win_if_needed(char *titre)
{
strcpy(fname, titre);
debut_arbre();
set_tree_w_title(fd_unrooted->unrooted, titre, notu);
tree_draw_proc(fd_unrooted->tree_plot);
}


void debut_arbre(void)
{
bignoeud *centre, *bigracine;
char *pname;
int tot;
double delta, current_angle, width, height;
branche *fin_branche, *br;
double radius;


/* read tree file */
if( (pname=preptree(fname)) != NULL ) {
        char mess[150];
        strcpy(mess, pname); strcat(mess, fname);
        err_message(mess);
        notu = 0;
		Nlm_UseWindow(fd_unrooted->unrooted);
        return;
        }
bigracine = cre_new_tree(racine, NULL, NULL);
if(!has_br_length) 
        bigracine->l2 = 0;
centre = bigracine->v1;
/* changer de centre si c'est une feuille */
tot=0;
if(centre->v1 == NULL) tot++;
if(centre->v2 == NULL) tot++;
if(centre->v3 == NULL) tot++;
if(tot>= 2) centre = bigracine->v2;
remove_big_root(bigracine);
radius = calc_dist_centre_feuilles(centre, NULL);
/* on va remplacer les branches nulles par des branches
100 fois plus courtes que la plus courte de l'arbre
*/
mini_br_length = proc_null_neg_branches(centre, NULL);
mini_br_length /= 100; 
delta = 2*M_PI / (notu+1);
current_angle = 0;
set_angles_noeuds(centre, NULL, delta, &current_angle, radius);
centre->position.r = 0;
centre->position.angle = 0;
calc_cartesienne(&(centre->position));
branches = (branche *)check_alloc(2*(notu+1)-3, sizeof(branche));
fin_branche = calc_position_noeuds(centre, NULL, branches);
/* calcul extremites du graphique des branches */
log_min_x = log_max_x = branches->debut.x;
log_min_y = log_max_y = branches->debut.y;
for(br = branches; br < fin_branche; br++) {
        if(br->debut.x < log_min_x) log_min_x = br->debut.x;
        if(br->debut.x > log_max_x) log_max_x = br->debut.x;
        if(br->fin.x < log_min_x) log_min_x = br->fin.x;
        if(br->fin.x > log_max_x) log_max_x = br->fin.x;
        if(br->debut.y < log_min_y) log_min_y = br->debut.y;
        if(br->debut.y > log_max_y) log_max_y = br->debut.y;
        if(br->fin.y < log_min_y) log_min_y = br->fin.y;
        if(br->fin.y > log_max_y) log_max_y = br->fin.y;
        }
width = log_max_x - log_min_x;
height = log_max_y - log_min_y;
if( width > height ) {
        log_min_y -= (width - height)/2;
        log_max_y += (width - height)/2;
        }
else if( height > width ) {
        log_min_x -= (height - width)/2;
        log_max_x += (height - width)/2;
        }
}



void my_draw_line(int dx, int dy, int fx, int fy, int col, int ps, int h)
{
if(ps) {
#ifdef WIN_MAC
        MoveTo(dx + margin, dy + margin);
        LineTo(fx + margin, fy + margin);
#else
        fprintf(plotfile, "%d %d moveto ", dx, h - dy);
        fprintf(plotfile, "%d %d lineto stroke\n", fx, h - fy);
#endif
        }
else {
        Nlm_MoveTo(dx, dy);
        Nlm_LineTo(fx, fy);
        }
}



void my_draw_text(int align, int x, int y, int color, int size, int style, 
        char *chaine, int ps, int h, int ascent)
{
int w;
if( align == ALIGN_CENTER) 
        w = my_get_string_width(chaine)/2;
else    w = 0;
if( ps ) {
#ifdef WIN_MAC
        static char copy[255];
        int l;
        MoveTo(x + margin - w, y + margin);
        l = strlen(chaine);
        *copy = l;
        memcpy(copy + 1, chaine, l + 1);
        DrawString( (ConstStr255Param) copy);
#else
        fprintf(plotfile, "%d %d moveto ", x - w, h - y);
        fprintf(plotfile,"(%s) show\n",chaine);
#endif
        }
else 	{
        if(color) Nlm_Red();
#ifdef WIN_MSWIN
        Nlm_PaintStringEx(chaine, x-w, (doing_copy ? y-ascent : y) );
#else
        Nlm_MoveTo(x - w, y);
        Nlm_PaintString(chaine);
#endif
        if(color) Nlm_Black();
        }
}



hit draw_line(branche *br, plot_data *ob, int doit, int char_width, 
        int char_height, int descend, int ascent)
{
double h, w;
double angle;
cp_point phys_debut, phys_fin;
hit result;


log_to_phys(&br->debut, &phys_debut);
log_to_phys(&br->fin, &phys_fin);
w = phys_fin.x - phys_debut.x;
h = phys_fin.y - phys_debut.y;
if( doit ) {
        my_draw_line(phys_debut.x + ob->x, phys_debut.y + ob->y, 
                phys_fin.x + ob->x, phys_fin.y + ob->y, ob->lcol, ob->ps_plot,
                ob->h);
        }
if(br->nom != NULL) {
/* calcule de l'angle tel que vu sur le dessin: il faut utiliser -h
   car les coord en y sont calculees avec 0 en haut
   aussi mettre angle dans [0 , 2*pi[
*/
        angle = atan( (-h) / w );
        if( w < 0 ) angle = M_PI + angle;
        if(angle < 0) angle += 2*M_PI;
        result = draw_name_angle(&phys_fin, br->nom, ob, angle, doit, 
                char_width, char_height, descend, br->color, ascent);
        if( result != NO_HIT ) 
                return result;
        }
return NO_HIT;
}



hit draw_name_angle(cp_point *phys_pos, char *nom, plot_data *ob, double angle,
        int doit, int char_width, int char_height, int descend, int color,
        int ascent)
{
double x, y, x_delta, y_delta;
point new_position;
int stepx;
hit result;
static char lettre[2] = "A";


#ifdef WIN_MAC
if(doit && ob->ps_plot) PicComment(140, 0, nil); /* pour grouper les caract d'un nom */
#endif


x = phys_pos->x;
y = phys_pos->y;
if((angle>= 0 && angle <= M_PI/4) || (angle >= 7*M_PI/4 && angle <= 2*M_PI)) {
        stepx = TRUE;
        }
else if((angle > M_PI/4 && angle < M_PI/2) || 
                (angle > 3*M_PI/2 && angle < 7*M_PI/4)) {
        stepx = FALSE;
        }
else if((angle >= M_PI/2 && angle <= 3*M_PI/4) || 
                (angle >= 5*M_PI/4 && angle <= 3*M_PI/2)) {
        stepx = FALSE;
        y_delta = char_height;
        y_delta -= descend;
        y_delta *= ( (int) strlen(nom) );
        if(angle >= M_PI/2 && angle <= 3*M_PI/4) y -= y_delta;
        else    y += y_delta;
        x_delta = y_delta / fabs( tan(angle) );
        x -= x_delta;
        angle -= M_PI;
        if(angle < 0) angle += 2*M_PI;
        }
else    {
        stepx = TRUE;
        x_delta = my_get_string_width(nom);
        y_delta = x_delta * tan(angle);
        x -= x_delta;
        y += y_delta;
        angle -= M_PI;
        if(angle < 0) angle += 2*M_PI;
        }
do      {
        new_position.x = x + 0.5;
        new_position.y = y + 0.5;
        result = draw_letter(&new_position, *nom, ob, doit, char_width, 
        	color, ascent);
        if( result != NO_HIT ) 
                return result;
        if(stepx) {
                lettre[0] = nom[0];
                x_delta = my_get_string_width(lettre);
                x += x_delta;
                y_delta = x_delta * tan(angle);
                y -= y_delta;
                }
        else    {
                y_delta = char_height;
                y_delta -= descend;
                x_delta = y_delta / fabs( tan(angle) );
                x += x_delta;
                if(angle <= M_PI/2) y -= y_delta;
                else    y += y_delta;
                }
        nom++;
        }
while(*nom != 0);
#ifdef WIN_MAC
if(doit && ob->ps_plot) PicComment(141, 0, nil); /* fin du groupe */
#endif
return NO_HIT;
}



hit draw_letter(point *position, char lettre, plot_data *ob, int doit, 
        int char_width, int color, int ascent)
{
static char chaine[2] = " ";
double maxx;
if(ob->ps_plot) {
        maxx = ob->w;
        }
else    {
        maxx = ob->w - char_width;
        }
/* clipping */
if( position->x <= 0 ) return LEFT_HIT;
if( position->x >= maxx ) return RIGHT_HIT;
if( position->y - ascent < 0 )  return TOP_HIT;
if( position->y >= ob->h ) return BOTTOM_HIT;
chaine[0] = lettre;
if(doit && lettre != ' ') 
        my_draw_text(ALIGN_LEFT, position->x + ob->x, position->y + ob->y,
                color, ob->lsize, ob->lstyle, chaine, ob->ps_plot, 
                ob->h, ascent);
return NO_HIT;
}


void draw_scale(plot_data *ob, int charheight, int ascent)
{
char ech_name[20];
int phys_w, y, xd, xf;
double log_val;


log_val = calc_echelle(ob->w);
phys_w = myrint(length_log_phys(log_val));
y = 1.5 * charheight;
xf = ob->w * 0.97;
xd = xf - phys_w;
my_draw_line(xd + ob->x, y + ob->y, xf + ob->x, y + ob->y, ob->lcol, 
        ob->ps_plot, ob->h);
sprintf(ech_name, "%.1g", log_val);
my_draw_text(ALIGN_CENTER, (xd + xf)/2 + ob->x,
        y - charheight/3 + ob->y, 
        ob->lcol, ob->lsize, ob->lstyle, ech_name, ob->ps_plot, ob->h,
        ascent);
my_draw_line(xd + ob->x, y - charheight/3 + ob->y, 
        xd + ob->x, y + charheight/3 + ob->y, ob->lcol, ob->ps_plot, ob->h);
my_draw_line(xf + ob->x, y - charheight/3 + ob->y, 
        xf + ob->x, y + charheight/3 + ob->y, ob->lcol, ob->ps_plot, ob->h);
}



int my_get_char_height(int *ascend, int *descend)
{
if(plot_aux_data.ps_plot) {
#ifdef WIN_MAC
        static FontInfo myinfo;
        GetFontInfo(&myinfo);
        *ascend = myinfo.ascent; *descend = myinfo.descent;
#else
        int font_size;
        font_size = list_font_size[current_font_size_rank];
        *ascend = 0.80 * font_size;
        *descend = font_size - *ascend;
#endif
        }
else	{
#ifdef WIN_MSWIN
	HDC picHDC = Nlm_currentHDC;
	if(doing_copy) {
		Nlm_currentHDC = Nlm_GetPicWinHDC();
		Nlm_SelectFont(current_font);
		}
#endif
        *ascend = Nlm_Ascent(); *descend = Nlm_Descent();
#ifdef WIN_MSWIN
	if(doing_copy) Nlm_currentHDC = picHDC;
#endif
	}
return *ascend + *descend ;
}


int my_get_string_width(char *nom)
{
int w;

if(plot_aux_data.ps_plot) {
#ifdef WIN_MAC
        static char nom255[256];
        w = strlen(nom); memcpy(nom255+1, nom, w); nom255[0] = w;
        w = StringWidth( (ConstStr255Param) nom255);
#else   
        w = list_font_size[current_font_size_rank];
        w = w * 0.61 * (int)strlen(nom);
#endif
        }
else	{
#ifdef WIN_MSWIN
	HDC picHDC = Nlm_currentHDC;
	if(doing_copy) {
		Nlm_currentHDC = Nlm_GetPicWinHDC();
		Nlm_SelectFont(current_font);
		}
#endif
        w = Nlm_TextWidth(nom, strlen(nom)) ;
#ifdef WIN_MSWIN
	if(doing_copy) Nlm_currentHDC = picHDC;
#endif
        }
return w;
}



char *check_alloc(int nbrelt, int sizelt)
{
char *retval;
if( (retval = calloc(nbrelt,sizelt)) != NULL ) return retval;
err_message("ERROR: not enough memory.");
}




void process_keys(Nlm_Char key)
{
#ifdef WIN_MAC
        /* combinaison command-lettre */
        #define test_key(a) Nlm_cmmdKey && Nlm_currentKey == a
#else
        /* combinaison control-lettre */
        #define test_key(a)  key == a - 96
#endif
if( test_key('q') ) 
        exit(0);
else if( test_key('o') )
        open_callback(fd_unrooted->open_button);
else if( test_key('s') )
        save_plot_callback(fd_unrooted->save_plot_button);
#ifdef MICRO
else if( test_key('c') )
        copy_callback(NULL);
else if( test_key('p') )
	print_plot(NULL);
#endif
else if( test_key('v') )
		paste_tree(NULL);
else if( test_key('f') )
        search_callback(NULL);
else if( test_key('a') ) {
        if(Nlm_Enabled(fd_unrooted->AgainItem) ) 
                search_callback(fd_unrooted->AgainItem);
        }
else if( test_key('r') )
        reset_callback(NULL);
return;
} 


void string_callback(Nlm_ButtoN item)
{
Nlm_Boolean *pdone;

pdone = (Nlm_Boolean *)Nlm_GetWindowExtra( Nlm_ParentWindow(item) );
*pdone = TRUE;
}


void majuscules(char *p)
{
while(*p != 0) {*p = toupper(*p); p++;}
}


void search_callback(Nlm_IteM item)
{
static Nlm_WindoW win;
static Nlm_Boolean search_done;
static Nlm_TexT select_box;
static char select[100];
char aux[100];
int num, trouve, totnoms;
static int first = TRUE;

if(item != fd_unrooted->AgainItem) {
	if(first) {
		win = Nlm_FixedWindow(-50,-50,-5,-5,"Name searched:",NULL);
		select_box=Nlm_DialogText(win,"",15,NULL);
		Nlm_Advance(win);
		Nlm_DefaultButton(win,"ok",string_callback);
		Nlm_SetWindowExtra(win, &search_done, NULL);
		first = FALSE;
		}
	Nlm_Show(win);
	Nlm_Select(win); Nlm_Select(select_box);
	search_done=FALSE;
	while(! search_done) {
		myWaitAndProcessNextEvent();
		}
	Nlm_GetTitle(select_box, select, sizeof(select));
	majuscules(select);
	Nlm_Hide(win);
	Nlm_Enable(fd_unrooted->AgainItem);
	Nlm_Select(fd_unrooted->unrooted);
	Nlm_FlushEvents();
	}
if(strlen(select) == 0) return;
trouve = FALSE; totnoms = 2 * (notu + 1) - 3;
for(num = 0; num < totnoms; num++) {
	if(branches[num].nom == NULL) continue;
	strcpy(aux, branches[num].nom);
	majuscules(aux);
	if(strstr( aux, select) != NULL) {
		trouve = TRUE;
		branches[num].color = TRUE;
		}
	}
if(trouve) {
	tree_draw_proc(fd_unrooted->tree_plot);
	}
}


void reset_callback(Nlm_IteM item)
{
int num;
for(num = 0; num < 2 * (notu + 1) - 3; num++) 
		branches[num].color = FALSE;
tree_draw_proc(fd_unrooted->tree_plot);
}


void clear_tree(Nlm_IteM unused)
{
if(notu != 0) free_tree();
notu = 0;
Nlm_SetTitle(fd_unrooted->unrooted, "unrooted");
tree_draw_proc(fd_unrooted->tree_plot);
}


#ifdef MICRO

int print_title(int x, int y, char *text)
{
static char ligne[200];
time_t heure;
int h, ascent, descent;

time (&heure);
#ifdef WIN_MAC
text = mac_fname_to_roman(text);
#endif
sprintf(ligne, "unrooted    %s    %s", text, ctime(&heure) );
h = strlen(ligne); ligne[h - 1] = 0;
h = my_get_char_height(&ascent, &descent);
Nlm_PaintStringEx(ligne, x, y + ascent);
return h;
}


void print_plot(Nlm_IteM unused)
{
Nlm_WindoW w;
int h;

if(notu == 0) return;
w = Nlm_StartPrinting();
if(w == NULL) return;
doing_print = TRUE;
Nlm_PrintingRect(&print_rect);
Nlm_StartPage();
Nlm_ClipPrintingRect( &print_rect);
Nlm_SelectFont(current_font);
h = print_title(print_rect.left, print_rect.top, fname);
print_rect.top += h;
plot_aux_data.comp_phys_bounds = TRUE;
tree_draw_proc(fd_unrooted->tree_plot);
Nlm_EndPage();
Nlm_EndPrinting(w);
doing_print = FALSE;
}

#endif

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

