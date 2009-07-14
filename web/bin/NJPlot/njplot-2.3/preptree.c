#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef FALSE
#define FALSE 0
#define TRUE (!FALSE)
#endif

#ifndef M_PI
#define M_PI           3.14159265358979323846  /* pi */
#endif
#ifndef M_PI_2
#define	M_PI_2		(M_PI/2)
#endif


/* typedefs */

typedef struct {
	unsigned int left : 1;
	unsigned int right : 1;
	unsigned int top : 1;
	unsigned int bottom : 1;
	} hits_field;

struct noeud {
	double l1,l2,l3;
	struct noeud *v1,*v2,*v3;
	char *nom;
	} **tabtax;
	
typedef struct { /* une branche definie par ses deux extremites */
	struct noeud *bouta;
	struct noeud *boutb;
	char *br_label;
	} branche_noeud;

typedef struct _cp_point {
	double x, y, r, angle;
	} cp_point;

typedef struct _branche {
	cp_point debut, fin;
        int color;
	char *nom;
	} branche;

typedef struct _bignoeud {
	struct _bignoeud *v1, *v2, *v3;
	double l1, l2, l3;
	cp_point position;
	char *nom;
	} bignoeud;
	
typedef struct _plot_data {
	int x, y, w, h;
	int ps_plot;
	int lstyle, lcol, lsize;
	int comp_phys_bounds;
	} plot_data;
	
typedef enum { NO_HIT, LEFT_HIT, RIGHT_HIT, TOP_HIT, BOTTOM_HIT } hit;

#define s_noeud sizeof(struct noeud)

char *preptree(char *fname);
void make_binary_or_unrooted(char *arbre);
int make_binary(char *arbre, char *debut, char *fin, int go_down);
void loadphylip(char *arbre);
struct noeud *unrootedset(char *arbre, char *deb, char *fin, branche_noeud **p_int_br);
char *nextpar(char *text, char *pospar);
void free_tree(void);
double tree_bal(struct noeud *centre, struct noeud *p1, struct noeud *p2);
void get_next_br(struct noeud *pere, struct noeud *fils, 
		double **abdown, double **abup, struct noeud **f1, struct noeud **f2);
void place_midpoint_root(void);
void parcourir_branches(struct noeud *centre, struct noeud *origine);
void process_branche(struct noeud *cote1, struct noeud *cote2, double length);
double get_length_down(struct noeud *pere, struct noeud *racine);
double arrondi_echelle(double x);
double calc_echelle(int larg);
int myrint(double val);
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
int chg_phys(plot_data *ob, hits_field hits, int char_height);
void log_to_phys(cp_point *log_pos, cp_point *phys_pos);
double length_log_phys(double p);
double length_phys_log(double p);

/* external functions */
char *check_alloc(int nbrelt, int sizelt);
void err_message(char *text);
hit draw_line(branche *br, plot_data *ob, int doit, int char_width,
	int char_height, int descend, int ascent);
int my_get_char_height(int *ascend, int *descend);
int my_get_string_width(char *nom);
void draw_scale(plot_data *ob, int char_height, int ascent);

/* globals */
struct noeud *racine;
double root_br_l;
int has_br_length = 0, notu, totbranches, 
	num_noeud, rooted, nextotu;
static branche_noeud *branches;
/* variables globales pour memoriser la meilleure branche:
celle qui partage l'arbre en 2 parties les plus egales possibles en prof moyenne
*/
static double current_best_diff, current_br_length, current_balance;
static struct noeud *current_cote1=NULL, *current_cote2=NULL;

/* externals */
extern int phys_min_x, phys_min_y, phys_max_x, phys_max_y;
extern double log_min_x, log_min_y, log_max_x, log_max_y, mini_br_length;



char *preptree(char *fname)
{
int i, c, v, steparbre, maxarbre;
FILE *njfile;
char *arbre, *der_arbre, *finarbre, *tmp;

if( (njfile=fopen(fname,"r")) == NULL )return ("Tree file not found.");
/* recherche du debut de la description de l'arbre */
do c=fgetc(njfile);
while (c == '\n'  || c == '\r' || c == ' ' );
/* for fastDNAml format, skip initial [comment] */
if(c == '[') {
	do c=fgetc(njfile); while (c != ']');
	do c=fgetc(njfile); while (c == '\n' || c == '\r' || c == ' ' );
	}
if ( c != '(') goto erreur;
	/* lecture de l'arbre par paquets de steparbre caracteres*/
	steparbre=1000;
	maxarbre=steparbre;
	arbre=check_alloc(maxarbre,1);
	der_arbre = arbre+maxarbre;
	*arbre=c;
	notu=2; i=3; v = 0;
	finarbre=arbre;
	while( (c=fgetc(njfile)) != EOF && c != ';') {
		if( /* c==' ' || */ c=='\n' || c=='\r') continue;
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
		if(c==')')notu++;
		if(c=='(')i++;
		if(c == ',') v++;
		}
	*(finarbre+1)='\0';
fclose(njfile);
if(i!=notu)goto erreur;

arbre = (char *)realloc(arbre, strlen(arbre) + 4 * v + 5 ); /* worst case add 4 chars for each , */
make_binary_or_unrooted(arbre);
/* notu--; */
notu = v ; /* after this notu = number of OTUs - 1  */
totbranches= -1;

/* allocate all memory */
tabtax = (struct noeud **)check_alloc(2*notu +1,sizeof(struct noeud *));
if(notu > 2) branches = 
	(branche_noeud *)check_alloc(notu-2, sizeof(branche_noeud));
for(i=0; i<2*notu+1; i++) *(tabtax+i)=
			(struct noeud *)check_alloc(1,s_noeud);
loadphylip(arbre);
free(arbre);
if(notu == 0) goto erreur;
if(!rooted) {
	racine = *(tabtax+(++num_noeud));
	if(num_noeud >= 2*notu + 1) return("Error: incorrect tree file");
	if(has_br_length) {
		place_midpoint_root();
		}
	else	{
		struct noeud *centre = *(tabtax + num_noeud - 1);
		racine->v1 = centre->v1;
		racine->v2 = centre;
		racine->v3 = NULL;
		centre->v1 = racine;
		if(racine->v1->v1 == centre)
			racine->v1->v1 = racine;
		else if(racine->v1->v2 == centre)
			racine->v1->v2 = racine;
		else
			racine->v1->v3 = racine;
		}
	}
else	{
	racine = *(tabtax+num_noeud);
	root_br_l= racine->l1 + racine->l2;
	if(!has_br_length) tree_bal(racine,racine->v1,racine->v2);
	}
if(notu+1<3) return ("Tree should contain at least 3 elements.");
return (NULL);

erreur:
return ("Not a valid tree file: ");
} /* end of preptree */


void make_binary_or_unrooted(char *arbre)
{
char *finarbre, *deba, *debb, *debc;

finarbre= nextpar(arbre,arbre);
*(finarbre + 1) = 0;
deba=arbre+1;
debb=deba;
while(*debb != ',') {
	if(*debb == 0) { err_message("Incorrect tree file"); exit(1); }
	if(*debb == '(')debb=nextpar(debb,debb);
	debb++;
	}
debb++;
debc=debb;
while(*debc != ',' && debc<finarbre) {
	if(*debc == '(')debc=nextpar(debc,debc);
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
		q = nextpar(p,p);
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
	/* make sure tree has branch lengths */
	if( memchr(debut, ':', fin - debut + 1) == NULL) {
		err_message("Cannot process multibranched tree without branch lengths");
		exit(1);
		}
	/* recherche de la 2eme virgule */
	p = debut; l = 0;
	while(TRUE) {
		if(*p == ',') {
			l++;
			if(l == 2) break;
			}
		else if(*p == '(')	{
			p = nextpar(p,p);
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


void loadphylip(char *arbre)
{
char *deba,*debb,*debc, *finarbre;
struct noeud *p1, *p2, *p3, *p;
branche_noeud *int_br_g, *int_br_d;
int i;

has_br_length=2;
/* ignore all stuff after last closing parenthesis 
(needed for fastDNAml output)
*/
finarbre= nextpar(arbre,arbre);
rooted=0;
deba=arbre+1;
debb=deba;
while(*debb != ',') {
	if(*debb == '(')debb=nextpar(arbre,debb);
	debb++;
	}
debb++;
debc=debb;
while(*debc != ',' && debc<finarbre) {
	if(*debc == '(')debc=nextpar(arbre,debc);
	debc++;
	}
if(*debc==',') {
/* the tree is unrooted <==> it has 3 subtrees at its bottommost level */
	debc++;
	}
else	{
/* the tree is rooted */
	debc=finarbre+1;
	rooted=1;
/*
	notu--;  //useless now
*/
	}
num_noeud=notu;
nextotu= -1;
p1=unrootedset(arbre,deba,debb-2,&int_br_g);
if(p1 != NULL) p2=unrootedset(arbre,debb,debc-2,&int_br_d);
if(p1 == NULL || p2 == NULL || num_noeud + 1 >= 2*notu + 1) {
	notu = 0;
	return;
	}
p = *(tabtax+(++num_noeud));
if(!has_br_length) {
	p1->l3 = 0.5*p1->l3;
	p2->l3 = 0.5*p2->l3;
	}
p->v1=p1; p1->v3=p; p->l1=p1->l3;
if(int_br_g!=NULL) { int_br_g->bouta=p; int_br_g->boutb=p1; }
p->v2=p2; p2->v3=p; p->l2=p2->l3;
if(int_br_d!=NULL) { int_br_d->bouta=p; int_br_d->boutb=p2; }
if(!rooted) {
	p3=unrootedset(arbre,debc,finarbre-1,&int_br_g);
	if(p3 == NULL ) {
		notu = 0;
		return;
		}
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
		totbranches++;
		branches[totbranches].br_label=check_alloc(l+1,1);
		memcpy(branches[totbranches].br_label,finarbre+1,l);
		branches[totbranches].br_label[l]=0;
		branches[totbranches].bouta=p1;
		branches[totbranches].boutb=p2;
		}
	}
}


struct noeud *unrootedset(char *arbre, char *deb, char *fin, branche_noeud **p_int_br)
{
struct noeud *p;
char *virg;
int l;
branche_noeud *int_br_g, *int_br_d;

*p_int_br=NULL;
while(*deb==' ' || *deb=='\'')deb++;
while(*fin==' ')fin--;
virg=deb;
while(*virg != ',' && virg < fin) {
	if(*virg == '(') virg=nextpar(arbre,virg);
	virg++;
	}
if(virg>fin) virg=deb;
if(*virg == ',') {
	struct noeud *p1,*p2;
	p1=unrootedset(arbre,deb,virg-1,&int_br_g);
if(p1 == NULL) return NULL;
	p2=unrootedset(arbre,virg+1,fin,&int_br_d);
if(p2 == NULL) return NULL;
if(num_noeud + 1 >= 2*notu + 1) 
			return NULL;
	p = *(tabtax+(++num_noeud));
	p->v1=p1; p1->v3=p; p->l1=p1->l3;
	if(int_br_g!=NULL) { int_br_g->bouta=p; int_br_g->boutb=p1; }
	p->v2=p2; p2->v3=p; p->l2=p2->l3;
	if(int_br_d!=NULL) { int_br_d->bouta=p; int_br_d->boutb=p2; }
	}
else	{
	double brlength;
	virg=deb;
	while(*virg != ':' && virg < fin) {
		if(*virg=='(')virg=nextpar(arbre,virg);
		virg++;
		}
	if(virg>fin) virg=deb;
	if(*virg == ':') {
		if(has_br_length == 0) goto problem;
		sscanf(virg+1,"%le",&brlength);
		virg--;
		has_br_length=1;
		}
	else	{
		if(has_br_length == 1) goto problem;
		brlength=1;
		has_br_length=0;
		}
	if(*deb == '(') {
		char *fpar;
		branche_noeud *prov;
		fpar=nextpar(arbre,deb)-1;
		p=unrootedset(arbre,deb+1,fpar,&prov);
if(p == NULL) return NULL;
/* recherche internal label */
		l=virg-fpar-1;
		if(l>0) {
			totbranches++;
			branches[totbranches].br_label=
				check_alloc(l+1,1);
			memcpy(branches[totbranches].br_label,fpar+2,l);
			branches[totbranches].br_label[l]=0;
			*p_int_br= &branches[totbranches];
			}
		}
	else	{
		size_t n;
		if( virg-1>=deb && *virg=='\'' )virg--;
		n=virg-deb+1;
		++nextotu;
		p= *(tabtax+nextotu);
		p->nom = (char *)check_alloc(n+1,1);
		memcpy(p->nom, deb, n); p->nom[n] = 0;
		p->v1=p->v2=p->v3=NULL;
		}
	p->l3=brlength;
	}
return p;
problem:
err_message("Error: Inconsistent tree file for branch lengths.");
}


char *nextpar(char *text, char *pospar)
{
char *pos;
pos=pospar+1;
while(*pos != ')') {
	if(*pos == '(') pos=nextpar(text,pos);
	pos++;
	}
return pos;
}


void free_tree(void)
{
int i;
if(notu == 0) return;
/* de-allocate all memory */
for(i=0; i<2*notu+1; i++) {
	if(tabtax[i]->nom != NULL) free(tabtax[i]->nom);
	free(tabtax[i]);
	}
free(tabtax);
for(i=0; i<notu-2; i++)
	if(branches[i].br_label != NULL) free(branches[i].br_label);
free(branches);
}


/* pour un arbre sans longueur de branche, les calculer de facon a ce
que toutes les feuilles arrivent a la meme profondeur */
double tree_bal(struct noeud *centre, struct noeud *p1, struct noeud *p2)
{
double ld, lg, pg, pd, *abup1, *abup2, *abdown1, *abdown2;
struct noeud *f1, *f2;
get_next_br(centre,p1,&abup1,&abdown1,&f1,&f2);
if(f1!=NULL) lg=tree_bal(p1,f1,f2);
else lg=0.0;
get_next_br(centre,p2,&abup2,&abdown2,&f1,&f2);
if(f1!=NULL) ld=tree_bal(p2,f1,f2);
else ld=0.0;
pg=lg+1; pd=ld+1;
if(pg>pd) pd=pg;
else pg=pd;
*abup1 = *abdown1 = pg-lg;
*abup2 = *abdown2 = pd-ld;
return (pg);
}


void get_next_br(struct noeud *pere, struct noeud *fils, 
double **abdown, double **abup, struct noeud **f1, struct noeud **f2)
/* pour une branche pere->fils donnee, calculer dans *f1, *f2 les deux autres
voisins de fils et dans *abdown et *abup les adresses des longeurs de la
branche pere->fils dans les deux sens */
{
if(pere->v1==fils) *abdown= &(pere->l1);
else if(pere->v2==fils) *abdown=&(pere->l2);
else *abdown=&(pere->l3);
if(fils->v1==pere) {
	*abup=&(fils->l1);
	*f1=fils->v2;
	*f2=fils->v3;
	}
else if(fils->v2==pere) {
	*abup=&(fils->l2);
	*f1=fils->v1;
	*f2=fils->v3;
	}
else 	{
	*abup=&(fils->l3);
	*f1=fils->v1;
	*f2=fils->v2;
	}
}


void place_midpoint_root(void)
/* enraciner l'arbre sans racine en cherchant son centre 
*/
{
struct noeud *aux;
double laux;

current_best_diff= 9.e99;
current_cote1=current_cote2=NULL;
parcourir_branches(*tabtax,NULL);
rooted = TRUE;
root_br_l = current_br_length;
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
current_cote1->v3 = racine;

if (current_cote2->v1 == current_cote1 )
	current_cote2->v1 = racine;
else if (current_cote2->v2 == current_cote1)
	current_cote2->v2 = racine;
else
	current_cote2->v3 = racine;
racine->v1=current_cote1;
racine->v2=current_cote2;
racine->v3=NULL;
racine->l3=0;
racine->l1=current_br_length*current_balance;
racine->l2=current_br_length - racine->l1;
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
/* calculer la prof moyenne des 2 cotes de la branche cote1<-->cote2
de longueur length
et memoriser la meilleure branche dans les variables globales
*/
{
double b1, b2, x, dist_root_side1, dist_root_side2, diff_betw_sides;

if(cote1==NULL || cote2==NULL) return;
b1=get_length_down(cote2,cote1);
b2=get_length_down(cote1,cote2);
if( fabs(length) > 1.e-5 )
	x=(b2-b1+length)/(2*length);
else
	x=0;
if(x<0) x=0;
if(x>1) x=1;
dist_root_side1=length*x+b1;
dist_root_side2=length*(1-x)+b2;
diff_betw_sides=fabs(dist_root_side1-dist_root_side2);
if(diff_betw_sides < current_best_diff ) {
	current_best_diff=diff_betw_sides;
	current_cote1=cote1;
	current_cote2=cote2;
	current_br_length=length;
	current_balance=x;
	}
}


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


bignoeud *cre_new_tree(struct noeud *debut, struct noeud *parent, 
			bignoeud *bigparent)
{
bignoeud *nouveau;
size_t l;
if(debut == NULL) return NULL;
nouveau = (bignoeud *)check_alloc(1, sizeof(bignoeud) );
if(debut->v1 == parent) {
	nouveau->v1 = bigparent;
	nouveau->v2 = cre_new_tree(debut->v2, debut, nouveau);
	nouveau->v3 = cre_new_tree(debut->v3, debut, nouveau);
	}
else if(debut->v2 == parent) {
	nouveau->v2 = bigparent;
	nouveau->v1 = cre_new_tree(debut->v1, debut, nouveau);
	nouveau->v3 = cre_new_tree(debut->v3, debut, nouveau);
	}
else	{
	nouveau->v3 = bigparent;
	nouveau->v1 = cre_new_tree(debut->v1, debut, nouveau);
	nouveau->v2 = cre_new_tree(debut->v2, debut, nouveau);
	}
if(has_br_length) {
	nouveau->l1 = debut->l1;
	nouveau->l2 = debut->l2;
	nouveau->l3 = debut->l3;
	}
else	{
	if(nouveau->v1 != NULL) nouveau->l1 = 1;
	if(nouveau->v2 != NULL) nouveau->l2 = 1;
	if(nouveau->v3 != NULL) nouveau->l3 = 1;
	}
if(debut->nom != NULL) { /*ajouter espace en tete du nom pour plus joli dessin*/
	l = strlen(debut->nom);
	nouveau->nom = (char *)check_alloc(l+2, 1);
	nouveau->nom[0] = ' ';
	memcpy(nouveau->nom + 1, debut->nom, l+1);
	free(debut->nom);
	debut->nom = nouveau->nom;
	}
return nouveau;
}


void remove_big_root(bignoeud *bigracine)
{
bignoeud *p1, *p2;
double root_br_l;

p1=bigracine->v1;
p2=bigracine->v2;
root_br_l = bigracine->l1 + bigracine->l2;
if(p1->v1 == bigracine )
	{p1->v1 = p2; p1->l1 = root_br_l;}
else if (p1->v2 == bigracine)
	{p1->v2 = p2; p1->l2 = root_br_l;}
else
	{p1->v3 = p2; p1->l3 = root_br_l;}
if(p2->v1 == bigracine )
	{p2->v1 = p1; p2->l1 = root_br_l;}
else if (p2->v2 == bigracine)
	{p2->v2 = p1; p2->l2 = root_br_l;}
else
	{p2->v3 = p1; p2->l3 = root_br_l;}
}


double calc_dist_centre_feuilles(bignoeud *debut, bignoeud *parent)
{
double valeur, current;
if(debut == NULL) return 0;
valeur = 0;  /*  !!!!! conserver cette ecriture sinon plante sur PC !!!!!! */
if(debut->v1 != parent) {
	current = calc_dist_centre_feuilles(debut->v1, debut);
	current += debut->l1;
	if(current > valeur) valeur = current;
	}
if(debut->v2 != parent) {
	current = calc_dist_centre_feuilles(debut->v2, debut);
	current += debut->l2;
	if(current > valeur) valeur = current;
	}
if(debut->v3 != parent) {
	current = calc_dist_centre_feuilles(debut->v3, debut);
	current += debut->l3;
	if(current > valeur) valeur = current;
	}
return valeur;
}


double proc_null_neg_branches(bignoeud *debut, bignoeud *parent)
/* mettre branches negatives a 0 !attention dans un seul sens!
retourner la + petite branche non nulle
*/
{
double valeur, current;
if(debut == NULL) return 0;
valeur = 1e50;
if(debut->v1 != parent) {
	if(debut->l1 < 0) debut->l1 = 0;
	current = debut->l1;
	if(current < valeur && current > 0) valeur = current;
	current = proc_null_neg_branches(debut->v1, debut);
	if(current < valeur && current > 0) valeur = current;
	}
if(debut->v2 != parent) {
	if(debut->l2 < 0) debut->l2 = 0;
	current = debut->l2;
	if(current < valeur && current > 0) valeur = current;
	current = proc_null_neg_branches(debut->v2, debut);
	if(current < valeur && current > 0) valeur = current;
	}
if(debut->v3 != parent) {
	if(debut->l3 < 0) debut->l3 = 0;
	current = debut->l3;
	if(current < valeur && current > 0) valeur = current;
	current = proc_null_neg_branches(debut->v3, debut);
	if(current < valeur && current > 0) valeur = current;
	}
return valeur;
}


int set_angles_noeuds(bignoeud *debut, bignoeud *parent, double delta,
				double *p_current_angle, double rayon)
{
int feuille = FALSE;
int poids, poids1, poids2;
double angle1 = -1, angle2 = -1;
static int totfeuilles=0;
char *nom;

if(debut == NULL) return 0;
if(debut->v1 != parent) {
	poids = set_angles_noeuds(debut->v1, debut, delta, p_current_angle, rayon);
	if(debut->v1 == NULL) feuille = TRUE;
	else if (parent != NULL) {
		angle1 = debut->v1->position.angle;
		poids1 = poids;
		}
	}
if(debut->v2 != parent) {
	poids = set_angles_noeuds(debut->v2, debut, delta, p_current_angle, rayon);
	if(debut->v2 == NULL) feuille = TRUE;
	else if (parent != NULL) {
		if ( angle1 == -1 ) {
			angle1 = debut->v2->position.angle;
			poids1 = poids;
			}
		else	{
			angle2 = debut->v2->position.angle;
			poids2 = poids;
			}
		}
	}
if(debut->v3 != parent) {
	poids = set_angles_noeuds(debut->v3, debut, delta, p_current_angle, rayon);
	if(debut->v3 == NULL) feuille = TRUE;
	else if (parent != NULL) {
		angle2 = debut->v3->position.angle;
		poids2 = poids;
		}
	}
if( feuille ) { totfeuilles++;
	debut->position.angle = *p_current_angle;
	*p_current_angle += delta;
	poids = 1;
	}
else if(parent != NULL)	{ /* faire angle moyen modulo 2.pi */
	debut->position.angle = (poids1*angle1 + poids2*angle2)/(poids1+poids2);
	poids = poids1 + poids2;
	if( angle1 > angle2 )
		debut->position.angle -= M_PI;
	}
debut->position.r = rayon;
calc_cartesienne(&(debut->position));
return poids;
}


void calc_cartesienne( cp_point *p)
{
p->x = p->r * cos(p->angle);
p->y = p->r * sin(p->angle);
return;
}


void calc_polaire( cp_point *p)
{
p->r = sqrt(p->x * p->x + p->y * p->y);
if( p->x == 0 ) {
	if( p->y == 0) p->angle = 0;
	else if (p->y > 0) p->angle = M_PI_2;
	else 	p->angle = 3*M_PI_2;
	}
else if (p->x > 0) {
	p->angle = atan( p->y / p->x );
	if( p->angle < 0 ) p->angle += 2*M_PI;
	}
else	{
	p->angle = M_PI - atan( p->y / p->x );
	}
return;
}


cp_point calc_point_direction( cp_point *depart, cp_point *direction, 
	double longueur)
{
static cp_point retour;
double lac, eps, tmp1, tmp2;

tmp1 = direction->x - depart->x;
tmp2 = direction->y - depart->y;
lac = sqrt( tmp1*tmp1 + tmp2*tmp2 );
/* on remplace les branches nulles par des branches tres courtes pour que le
calcul de l'angle en double soit bon mais que le dessin en entier soit le
meme
*/
if(longueur == 0) longueur = mini_br_length;
eps = longueur / lac;
retour.x = depart->x + eps * tmp1;
retour.y = depart->y + eps * tmp2;
calc_polaire(&retour);
return retour;
}


branche *calc_position_noeuds(bignoeud *debut, bignoeud *parent,
	branche *curr_branche)
{
if(debut->v1 != parent && debut->v1 != NULL ) {
	debut->v1->position = calc_point_direction( &(debut->position), &(debut->v1->position), debut->l1);
	curr_branche = calc_position_noeuds(debut->v1, debut, curr_branche);
	mem_line(&debut->position, &debut->v1->position, debut->v1, 
		curr_branche);
	curr_branche++;
	}
if(debut->v2 != parent && debut->v2 != NULL ) {
	debut->v2->position = calc_point_direction( &(debut->position), &(debut->v2->position), debut->l2);
	curr_branche = calc_position_noeuds(debut->v2, debut, curr_branche);
	mem_line(&debut->position, &debut->v2->position, debut->v2, 
		curr_branche);
	curr_branche++;
	}
if(debut->v3 != parent && debut->v3 != NULL ) {
	debut->v3->position = calc_point_direction( &(debut->position), &(debut->v3->position), debut->l3);
	curr_branche = calc_position_noeuds(debut->v3, debut, curr_branche);
	mem_line(&debut->position, &debut->v3->position, debut->v3, 
		curr_branche);
	curr_branche++;
	}
return curr_branche;
}


void mem_line(cp_point *debut, cp_point *fin, bignoeud *noeud_term, branche *br)
{
br->debut = *debut;
br->fin = *fin;
br->nom = noeud_term->nom;
}


void draw_tree(plot_data *ob, branche *branches, int comp_phys_bounds)
{
int char_width, char_height, num, dernier, encore = FALSE, ascend, descend;
hit result;
hits_field hits;
char_width = my_get_string_width("M");
char_height = my_get_char_height(&ascend, &descend);
dernier = 2*(notu+1)-3;
if( comp_phys_bounds ) {
/* rendre dimensions carrees */
	if (phys_max_x - phys_min_x > phys_max_y - phys_min_y)
		phys_max_x = phys_min_x + (phys_max_y - phys_min_y);
	else if(phys_max_y - phys_min_y > phys_max_x - phys_min_x)
		phys_max_y = phys_min_y + (phys_max_x - phys_min_x);

	do	{
		hits.left = hits.right = hits.top = hits.bottom = 0;
		for(num = 0; num < dernier; num++) {
			result = draw_line(branches+num, ob, FALSE, 
				char_width, char_height, descend, ascend);
			if(result == LEFT_HIT) hits.left = 1;
			else if(result == RIGHT_HIT) hits.right = 1;
			else if(result == TOP_HIT) hits.top = 1;
			else if(result == BOTTOM_HIT) hits.bottom = 1;
			}
		encore = chg_phys(ob, hits, char_height);
		}
	while (encore);
	for(num = 0; num < dernier; num++) branches[num].color = FALSE;
	}
for(num=0; num< dernier; num++) {
	draw_line(branches+num, ob, TRUE, char_width, char_height, 
		descend, ascend);
	}
if(has_br_length) draw_scale(ob, char_height, ascend);
}



int chg_phys(plot_data *ob, hits_field hits, int char_height)
{
int encore = FALSE;
const float delta = 0.05;
int value;
static int lr = 0, tb = 0, oldvaluex = 0, oldvaluey = 0;

if(hits.left == 1 && hits.right == 0) {
	if(oldvaluex == 0) value = (phys_max_x - phys_min_x) * delta; 
	else value = oldvaluex;
	if(value < 2) value = 2;
	if(lr == 2) { /* alternance left/right */
		if(oldvaluex <= 2) {
			hits.right = 1;
			return chg_phys(ob, hits, char_height);
			}
		value = oldvaluex/2; 
		}
	lr = 1;
	oldvaluex = value;
	phys_min_x += value;
	phys_max_x += value;
	if(phys_max_x <= phys_min_x) {
		phys_max_x = phys_min_x + 1;
		return FALSE;
		}
	encore = TRUE;
	}
else if(hits.left == 0 && hits.right == 1) {
	if(oldvaluex == 0) value = (phys_max_x - phys_min_x) * delta; 
	else value = oldvaluex;
	if(value < 2) value = 2;
	if(lr == 1) { /* alternance left/right */
		if(oldvaluex <= 2) {
			hits.left = 1;
			return chg_phys(ob, hits, char_height);
			}
		value = oldvaluex/2; 
		}
	lr = 2;
	oldvaluex = value;
	phys_max_x -= value;
	phys_min_x -= value;
	if(phys_max_x <= phys_min_x) {
		phys_max_x = phys_min_x + 1;
		return FALSE;
		}
	encore = TRUE;
	}
else if(hits.left == 1 && hits.right == 1) {
	lr = tb = 0;
	oldvaluex = 0;
	value = (phys_max_x - phys_min_x) * delta / 2; if(value < 2) value = 2;
	phys_min_x += value;
	phys_max_x -= value;
	if(phys_max_x <= phys_min_x) {
		phys_max_x = phys_min_x + 1;
		return FALSE;
		}
	encore = TRUE;
	phys_max_y = phys_min_y + (phys_max_x - phys_min_x);
	}

if(hits.top == 0 && hits.bottom == 1) {
	if(oldvaluey == 0) value = (phys_max_y - phys_min_y) * delta; 
	else value = oldvaluey;
	if(value < 2) value = 2;
	if(tb == 2) {
		if(oldvaluey <= 2) {
			hits.top = 1;
			return chg_phys(ob, hits, char_height);
			}
		value = oldvaluey/2; 
		}
	tb = 1;
	oldvaluey = value;
	phys_max_y -= value;
	phys_min_y -= value;
	if(phys_max_y <= phys_min_y) {
		phys_max_y = phys_min_y + 1;
		return FALSE;
		}
	encore = TRUE;
	}
else if(hits.top == 1 && hits.bottom == 0) {
	if(oldvaluey == 0) value = (phys_max_y - phys_min_y) * delta; 
	else value = oldvaluey;
	if(value < 2) value = 2;
	if(tb == 1) {
		if(oldvaluey <= 2) {
			hits.bottom = 1;
			return chg_phys(ob, hits, char_height);
			}
		value = oldvaluey/2; 
		}
	tb = 2;
	oldvaluey = value;
	phys_min_y += value;
	phys_max_y += value;
	if(phys_max_y <= phys_min_y) {
		phys_max_y = phys_min_y + 1;
		return FALSE;
		}
	encore = TRUE;
	}
else if(hits.top == 1 && hits.bottom == 1) {
	value = (phys_max_y - phys_min_y) * delta/2; if(value < 2) value = 2;
	lr = tb = 0;
	oldvaluey = 0;
	phys_min_y += value;
	phys_max_y -= value;
	if(phys_max_y <= phys_min_y) {
		phys_max_y = phys_min_y + 1;
		return FALSE;
		}
	encore = TRUE;
	phys_max_x = phys_min_x + (phys_max_y - phys_min_y);
	}

return encore;
}


void log_to_phys(cp_point *log_pos, cp_point *phys_pos)
{
double factorx, factory;

factorx = (phys_max_x - phys_min_x) / (log_max_x - log_min_x);
factory = (phys_max_y - phys_min_y) / (log_max_y - log_min_y);
phys_pos->x = factorx * ( log_pos->x - log_min_x ) + phys_min_x;
phys_pos->y = factory * ( log_pos->y - log_min_y ) + phys_min_y;
}


double length_log_phys(double p)
{
double factor;
factor = (phys_max_x - phys_min_x) / (log_max_x - log_min_x);
return p * factor;
}

double length_phys_log(double p)
{
double factor;
factor = (phys_max_x - phys_min_x) / (log_max_x - log_min_x);
return  p / factor;
}


double arrondi_echelle(double x)
{ /* arrondi x a une valeur 1, 2, 5 pour echelle */
double l, n;
int r;
static int corresp[] = {1,1,2,2,5,5,5,10,10,10,10,10};
                     /* 0,1,2,3,4,5,6, 7, 8, 9,10,11 */
l = log10(x);
n = floor(l);
l = x * pow(10, -n);
r = myrint(l); r = corresp[r];
return r * pow(10, n);
}


double calc_echelle(int larg)
{ /* rend taille logique pour echelle optimale */
double log_val, phys_val;
phys_val = larg/10;
log_val = length_phys_log(phys_val);
log_val = arrondi_echelle(log_val);
return log_val;
}


int myrint(double val)
{
if(val >= 0)
	return (int)(val + 0.5);
else
	return (int)(val - 0.5);
}



