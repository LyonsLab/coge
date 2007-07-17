#include <stdio.h>
#include <ctype.h>
#include "pratique.h"

/* ------------------------------------------------------------*/

void erreur(char *message)
{
printf("%s\n", message);
exit(1);
}

/* ------------------------------------------------------------*/

void *allouer(size_t taille)
{
void *pointeur;

if (taille == 0)	taille = 1;

pointeur = (void *) malloc(taille);

if (pointeur==NULL)
 	erreur("out of memory");
 else	return(pointeur);
}

/* ------------------------------------------------------------*/

void *reallouer(void *pointeur, size_t taille)
{
void *p;

p = pointeur;

pointeur = (void *) realloc(pointeur, taille);

if (pointeur==NULL)
  	erreur("out of memory");
else	return(pointeur);
}

/* ------------------------------------------------------------*/

void liberer(void *pointeur)
{

free(pointeur);
}

/* ------------------------------------------------------------*/

void **callouer_mat(size_t t_elt, size_t nb_lig, size_t nb_col)
{
void **pointeur;
int i;

pointeur = (void **) allouer(nb_lig * sizeof(void *));

for (i=0; i < nb_lig; i++)
	pointeur[i] = (void *) allouer(nb_col * t_elt);

return(pointeur);
}

/* ------------------------------------------------------------*/

void **recallouer_mat(void **pointeur, size_t t_elt, size_t anc_nb_lig, 
		size_t nb_lig, size_t nb_col)
{
int i;

if (anc_nb_lig == nb_lig)	return(pointeur);

for (i=nb_lig; i < anc_nb_lig; i++)
	liberer(pointeur[i]);

pointeur = (void **) reallouer(pointeur, nb_lig * sizeof(void *));

for (i=anc_nb_lig; i < nb_lig; i++)
	pointeur[i] = (void *) allouer(nb_col * t_elt);

return(pointeur);
}

/* ------------------------------------------------------------*/

void **recallouer_mat2(void **pointeur, size_t t_elt, 
		size_t anc_nb_lig, size_t nb_lig, size_t nb_col)
{
int i;
for (i=nb_lig; i < anc_nb_lig; i++)
	liberer(pointeur[i]);

pointeur = (void **) reallouer(pointeur, nb_lig * sizeof(void *));

for (i=0; i < min(anc_nb_lig, nb_lig); i++)
	pointeur[i] = (void *) reallouer(pointeur[i], nb_col * t_elt);

for (i=anc_nb_lig; i < nb_lig; i++)
	pointeur[i] = (void *) allouer(nb_col * t_elt);

return(pointeur);
}

/* ------------------------------------------------------------*/

void liberer_mat(void **pointeur, size_t nb_lig)
{
int i;

for (i=0; i < nb_lig; i++)
	liberer(pointeur[i]);

liberer(pointeur);
}

/* ------------------------------------------------------------*/

FILE *ouvrir(char *nomfich, char *mode)
{
FILE *f;

if ((f = fopen(nomfich, mode)) == NULL) 
   {
 	printf("fopen(\"%s\",\"%s\"): ", nomfich, mode);
	erreur("enable to open file");
   }
else	return f;
}

/* ------------------------------------------------------------*/

void fermer(FILE *f)
{
if (fclose(f) == EOF)
  	erreur("enable to close file");
}

/* ------------------------------------------------------------*/

void fcopie(FILE *fdestination, FILE *fsource)
{
char line[TAILLE_MAX_LIGNE_FICHIER];

while (fgets(line, TAILLE_MAX_LIGNE_FICHIER, fsource) != NULL)
 	fputs(line, fdestination);
 
fflush(fdestination);
}

/* ------------------------------------------------------------*/

void strmin(char *p)
{
char c;

for (; (c=*p); p++)
	*p = tolower(c);
}

/* ------------------------------------------------------------*/

void strmaj(char *p)
{
char c;

for (; (c=*p); p++)
	*p = toupper(c);
}

/* ------------------------------------------------------------*/

