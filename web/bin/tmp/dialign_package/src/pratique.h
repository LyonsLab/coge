#ifndef _PRATIQUE_H
#define _PRATIQUE_H

#include <stdio.h>
#include <stdlib.h>

#define	true		1
#define	false		0

#define	min(a,b)	((a)<(b)?(a):(b))
#define	max(a,b)	((a)>(b)?(a):(b))

#define	TAILLE_MAX_LIGNE_FICHIER	10000

void erreur(char *message);

void *allouer(size_t taille);
void *reallouer(void *pointeur, size_t taille);
void liberer(void *pointeur);
void **callouer_mat(size_t t_elt, size_t nb_lig, size_t nb_col);
void **recallouer_mat(void **pointeur, size_t t_elt, size_t anc_nb_lig, size_t nb_lig, size_t nb_col);
void **recallouer_mat2(void **pointeur, size_t t_elt, size_t anc_nb_lig, size_t nb_lig, size_t nb_col);
void liberer_mat(void **pointeur, size_t nb_lig);

FILE *ouvrir(char *nomfich, char *mode);
void fermer(FILE *f);
void fcopie(FILE *fdestination, FILE *fsource);

void strmin(char *p);
void strmaj(char *p);

#endif /* _PRATIQUE_H */
