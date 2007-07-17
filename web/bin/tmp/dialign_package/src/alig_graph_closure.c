
/************************************************************/
/************************************************************/
/**             GABIOS-LIB 1.0 (1999)                      **/
/** A library for Greedy Alignment of BIOlogical Sequences **/
/**             Developed by Said Abdeddaim                **/
/**          Said.Abdeddaim@dir.univ-rouen.fr              **/
/************************************************************/
/************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
 
#include "pratique.h"
#include "alig_graph_closure.h"

void alloc_closure(CLOSURE *clos);
void free_closure(CLOSURE *clos);
void realloc_closure(CLOSURE *clos);
void computeClosure(CLOSURE *clos);
int path(CLOSURE *clos, int x, int i, int y, int j);
void moveAligSet(CLOSURE *clos, int n1, int n2);
void read_closure(CLOSURE *clos, int nbreancr, int **ancrages);
void init_closure(CLOSURE *clos, int nbreancr, int **ancrages);
void init_seq(CLOSURE *clos, int nbreseq, int *longseq);
void desinit_seq(CLOSURE *clos);

int print_aligSets(CLOSURE *clos, int nseq, int i);

char DEBUG=0;

void computeClosure(CLOSURE *clos)
{
int **Succ, **Pred, *NSucc, *NPred, *npred;
int nsucc, ni, nj, s, top, bottom, n0, p, n, i, k, pos_n;
int x;

Succ = (int **) callouer_mat(sizeof(int), clos->nbrAligSets+2, clos->seqNbr);
Pred = (int **) callouer_mat(sizeof(int), clos->nbrAligSets+2, clos->seqNbr);
NSucc = (int *) allouer((clos->nbrAligSets+2) * sizeof(int));
NPred = (int *) allouer((clos->nbrAligSets+2) * sizeof(int));
npred = (int *) allouer((clos->nbrAligSets+2) * sizeof(int));
clos->topolog = (int *) allouer(sizeof(int));

/* C A L C U L des Succ[n][x] et NPred[n] */

for (n=1; n <= clos->nbrAligSets; n++)	NPred[n] = 0;

for (n=1; n <= clos->nbrAligSets; n++)
   {
	nsucc = 0;
	for (x=0; x < clos->seqNbr; x++)
	  if (clos->aligSet[n].pos[x] > 0)
	     {
		pos_n = clos->aligSet[n].pos[x];
		for (i=pos_n+1; i <= clos->seq[x].longueur && clos->seq[x].aligSetNbr[i] == 0; i++)
			clos->seq[x].predAligSetPos[i] = pos_n;
		if (i <= clos->seq[x].longueur)
		    {
			clos->seq[x].predAligSetPos[i] = pos_n;
		 	if (clos->aligSet[clos->seq[x].aligSetNbr[i]].nbr > 0)
			    {
				n0 = Succ[n][nsucc] = clos->seq[x].aligSetNbr[i];
				clos->aligSet[n0].nbr = - clos->aligSet[n0].nbr;
				nsucc++;
			    }
		    }
		for (i=pos_n-1; i > 0 && clos->seq[x].aligSetNbr[i] == 0; i--)
			clos->seq[x].succAligSetPos[i] = pos_n;
		if (i > 0) 	
			clos->seq[x].succAligSetPos[i] = pos_n;
	      }
	for (p=0; p < nsucc; p++)
	     {
		n0 = Succ[n][p];
		Pred[n0][NPred[n0]] = n;
		NPred[n0]++;

		clos->aligSet[n0].nbr = - clos->aligSet[n0].nbr;
	     }
	NSucc[n] = nsucc;
   }

/* C A L C U L de clos->topolog */

clos->topolog = (int *) reallouer(clos->topolog, (clos->nbrAligSets+2) * sizeof(int));

bottom = top = 0;

for (n=1; n <= clos->nbrAligSets; n++)
    {
	npred[n] = NPred[n];
	if (npred[n] == 0)	{ top++; clos->topolog[top] = n; }
    }


while( bottom != top)
   {
	bottom++;
	ni = clos->topolog[bottom];
	for (s=0; s < NSucc[ni]; s++)
	   {
		nj = Succ[ni][s];
		npred[nj]--;
		if (npred[nj] == 0)
		   {
			top++;
			clos->topolog[top] = nj;
		   }
	   }
   }


for (x=0; x < clos->seqNbr; x++)
    {
	clos->predFrontier[0][x] = 0;
	clos->succFrontier[clos->nbrAligSets+1][x] = clos->seq[x].longueur+1;
    }


for (k=1; k <= clos->nbrAligSets; k++)
    {
	n0 = clos->topolog[k];
	for (x=0; x < clos->seqNbr; x++)
	    {
		if (clos->aligSet[n0].pos[x] > 0)
			clos->predFrontier[n0][x] = clos->aligSet[n0].pos[x];
	    	else
		  for (p=0, clos->predFrontier[n0][x]=0; p < NPred[n0]; p++)
		    {
			n = Pred[n0][p];
			if (clos->predFrontier[n][x] > clos->predFrontier[n0][x]) 
				clos->predFrontier[n0][x] = clos->predFrontier[n][x];
		    }
	    }
    }


for (k=clos->nbrAligSets; k > 0; k--)
    {
	n0 = clos->topolog[k];
	for (x=0; x < clos->seqNbr; x++)
	    {
		if (clos->aligSet[n0].pos[x] > 0)
			clos->succFrontier[n0][x] = clos->aligSet[n0].pos[x];
		else
		  for (p=0, clos->succFrontier[n0][x]=clos->seq[x].longueur+1; 
					p < NSucc[n0]; p++)
		    {
			n = Succ[n0][p];
			if (clos->succFrontier[n][x] < clos->succFrontier[n0][x]) 
				clos->succFrontier[n0][x] = clos->succFrontier[n][x];
		    }
	    }
    }

liberer(npred); liberer(NPred); liberer(NSucc); 
liberer_mat((void **) Pred, clos->nbrAligSets+2); 
liberer_mat((void **) Succ, clos->nbrAligSets+2); 
liberer(clos->topolog);
}

void moveAligSet(CLOSURE *clos, int n1, int n2)
{
int x;
int k;

for (x=0; x < clos->seqNbr; x++)
  {
	k = clos->aligSet[n1].pos[x] = clos->aligSet[n2].pos[x];
	if (k > 0)	clos->seq[x].aligSetNbr[k] = n1;

	clos->predFrontier[n1][x] = clos->predFrontier[n2][x];
	clos->succFrontier[n1][x] = clos->succFrontier[n2][x];
  }

clos->aligSet[n1].nbr = clos->aligSet[n2].nbr;
}

void read_closure(CLOSURE *clos,  int nbreancr, int **ancrages)
{
FILE *f;
int x;
int i, ind, k, n;
int **Succ, **Pred, *NSucc, *NPred, *npred;

for (n=0; n < nbreancr; n++)
   {
	clos->nbrAligSets++;
	realloc_closure(clos);

	clos->aligSet[clos->nbrAligSets].nbr = 0;
	for (x=0; x < clos->seqNbr; x++)
	    {
		ind = clos->aligSet[clos->nbrAligSets].pos[x] = ancrages[n][x];
		if (ind > 0)
		   {
			clos->aligSet[clos->nbrAligSets].nbr++;
			clos->seq[x].aligSetNbr[ind] = clos->nbrAligSets;
		   }
	   }
   }

computeClosure(clos);
}

void init_closure(CLOSURE *clos, int nbreancr, int **ancrages)
{
int x;
int i, *longsequ;

longsequ = (int *) allouer(clos->seqNbr * sizeof(int));

for (x=0; x < clos->seqNbr; x++)
   {
	longsequ[x] = clos->seq[x].longueur;
	for (i=1; i <= clos->seq[x].longueur; i++)
		clos->seq[x].aligSetNbr[i] = clos->seq[x].succAligSetPos[i] 
			= clos->seq[x].predAligSetPos[i] = 0;
   }

clos->nbrAligSets = 0;

if (nbreancr > 0)
	read_closure(clos, nbreancr, ancrages);

for (x=0; x < clos->seqNbr; x++)
	clos->seq[x].longueur = longsequ[x];

liberer(longsequ);
}


void alloc_closure(CLOSURE *clos)
{
long nmax, na;
int x;

clos->predFrontier = (int **) callouer_mat(sizeof(int), clos->maxLong+2, clos->seqNbr+1); 	/* sera re'alloue' */
clos->succFrontier = (int **) callouer_mat(sizeof(int), clos->maxLong+2, clos->seqNbr+1); 	/* sera re'alloue' */

clos->aligSet = (positionSet *) allouer((clos->maxLong+2) * sizeof(positionSet));		/* sera re'alloue' */
for (na=0; na <= clos->maxLong+1; na++)
    {
	clos->aligSet[na].pos = (int *) allouer(clos->seqNbr * sizeof(int));
    }
clos->oldNbrAligSets = clos->maxLong;

for (x=0; x < clos->seqNbr; x++)
   {
	clos->seq[x].aligSetNbr = (int *) allouer((clos->seq[x].longueur+2)*sizeof(int));
	clos->seq[x].predAligSetPos = (int *) 
				allouer((clos->seq[x].longueur+2)*sizeof(int));
	clos->seq[x].succAligSetPos = (int *) 
				allouer((clos->seq[x].longueur+2)*sizeof(int));
   }

clos->gauche1 = (int *) allouer(clos->seqNbr * sizeof(int));
clos->gauche2 = (int *) allouer(clos->seqNbr * sizeof(int));
clos->droite1 = (int *) allouer(clos->seqNbr * sizeof(int));
clos->droite2 = (int *) allouer(clos->seqNbr * sizeof(int));
clos->pos_ = (int **) callouer_mat(sizeof(int), clos->seqNbr, clos->seqNbr);
}

void free_closure(CLOSURE *clos)
{
long nmax, na;
int x;

liberer(clos->gauche1); liberer(clos->gauche2); liberer(clos->droite1); liberer(clos->droite2); 
liberer_mat((void **) clos->pos_, clos->seqNbr); 

liberer_mat((void **) clos->succFrontier, clos->oldNbrAligSets+2); 
liberer_mat((void **) clos->predFrontier, clos->oldNbrAligSets+2);

for (x=0; x < clos->seqNbr; x++)
   {
	liberer(clos->seq[x].aligSetNbr); 
	liberer(clos->seq[x].predAligSetPos); 
	liberer(clos->seq[x].succAligSetPos); 
   }

for (na=0; na <= clos->oldNbrAligSets+1; na++)
   {
	liberer(clos->aligSet[na].pos); 
    }
liberer(clos->aligSet);
}

void realloc_closure(CLOSURE *clos)
{
int na;

if (clos->nbrAligSets > clos->oldNbrAligSets)
   {
	clos->predFrontier = (int **) recallouer_mat((void **) clos->predFrontier, sizeof(int), 
					clos->oldNbrAligSets+2, clos->nbrAligSets+2, clos->seqNbr+1);
	clos->succFrontier = (int **) recallouer_mat((void **) clos->succFrontier, sizeof(int), 
					clos->oldNbrAligSets+2, clos->nbrAligSets+2, clos->seqNbr+1);
	clos->aligSet = (positionSet *) reallouer(clos->aligSet, (clos->nbrAligSets+2) * sizeof(positionSet));
	for (na=clos->oldNbrAligSets+2; na <= clos->nbrAligSets+1; na++)
	    {
		clos->aligSet[na].pos = (int *) allouer(clos->seqNbr * sizeof(int));
	    }
	clos->oldNbrAligSets = clos->nbrAligSets;
   }
}


int print_aligSets(CLOSURE *clos, int nseq, int i)
{
char nouveau_, terminer;
int n, ng, nd, nn, k;
int x, y;

n = ng = nd = clos->seq[nseq].aligSetNbr[i];

if (ng == 0)
	   {
		k = clos->seq[nseq].predAligSetPos[i];
		if (k > 0) ng = clos->seq[nseq].aligSetNbr[k];
		k = clos->seq[nseq].succAligSetPos[i];
		if (k > 0) nd = clos->seq[nseq].aligSetNbr[k];
	   }

printf("echelle %d: ", n);
if (n != 0)	
   for (x=0; x < clos->seqNbr; x++) 
	printf("%d ", clos->aligSet[n].pos[x]);

printf("\nfrontiere clos->gauche %d: ", ng);
if (ng != 0)
   for (x=0; x < clos->seqNbr; x++) 
	printf("%d ", clos->predFrontier[ng][x]);

printf("\nfrontiere clos->droite %d: ", nd);
if (nd != 0)
   for (x=0; x < clos->seqNbr; x++) 
	printf("%d ", clos->succFrontier[nd][x]);

printf("\n");

}

void init_seq(CLOSURE *clos, int nbreseq, int *longseq)
{
int x;

clos->seqNbr = nbreseq;

clos->seq = (sequence *) allouer(clos->seqNbr * sizeof(sequence));

for (x=clos->maxLong=0; x < clos->seqNbr; x++)
   {
	clos->seq[x].longueur = longseq[x];
	if (clos->maxLong < longseq[x])
		clos->maxLong = longseq[x];
   }
}

void desinit_seq(CLOSURE *clos)
{
int x;

liberer(clos->seq);
}

/*********************************************************/
/************** EXTERN FONCTIONS *************************/
/*********************************************************/

CLOSURE *newAligGraphClosure(int nbreseq, int *longseq, 
				int nbreancr, int **ancrages)
{

CLOSURE *clos = (CLOSURE *) allouer(sizeof(CLOSURE));

init_seq(clos, nbreseq, longseq);

alloc_closure(clos); /* utilise clos->maxLong */

init_closure(clos, nbreancr, ancrages); 

return	clos;
}

void freeAligGraphClosure(CLOSURE *clos)
{
free_closure(clos);

desinit_seq(clos);

liberer(clos);
}

int addAlignedPositions(CLOSURE *clos, int seq1, int i, int seq2, int j)
{
char nouveau_, terminer;
int n, n1, n2, ng1, ng2, nd1, nd2, nn, k;
int x, y;

n1 = ng1 = nd1 = clos->seq[seq1].aligSetNbr[i];	n2 = ng2 = nd2 = clos->seq[seq2].aligSetNbr[j];

if (n1 == 0 || n2 == 0 || n1 != n2)
     {
	if (ng1 == 0)
	   {
		k = clos->seq[seq1].predAligSetPos[i];
		if (k > 0) ng1 = clos->seq[seq1].aligSetNbr[k];
		k = clos->seq[seq1].succAligSetPos[i];
		if (k > 0) nd1 = clos->seq[seq1].aligSetNbr[k];
	   }
	if (ng2 == 0)
	   {
		k = clos->seq[seq2].predAligSetPos[j];
		if (k > 0) ng2 = clos->seq[seq2].aligSetNbr[k];
		k = clos->seq[seq2].succAligSetPos[j];
		if (k > 0) nd2 = clos->seq[seq2].aligSetNbr[k];
	   }

	if (ng1 == 0)	for (x=0; x < clos->seqNbr; x++) clos->gauche1[x] = 0;
	else		for (x=0; x < clos->seqNbr; x++) clos->gauche1[x] = clos->predFrontier[ng1][x];
	if (nd1 == 0)	for (x=0; x < clos->seqNbr; x++) clos->droite1[x] = clos->seq[x].longueur + 1;
	else		for (x=0; x < clos->seqNbr; x++) clos->droite1[x] = clos->succFrontier[nd1][x];
	if (ng2 == 0)	for (x=0; x < clos->seqNbr; x++) clos->gauche2[x] = 0;
	else		for (x=0; x < clos->seqNbr; x++) clos->gauche2[x] = clos->predFrontier[ng2][x];
	if (nd2 == 0)	for (x=0; x < clos->seqNbr; x++) clos->droite2[x] = clos->seq[x].longueur + 1;
	else		for (x=0; x < clos->seqNbr; x++) clos->droite2[x] = clos->succFrontier[nd2][x];

	clos->gauche1[seq1] = clos->droite1[seq1] = i;
	clos->gauche2[seq2] = clos->droite2[seq2] = j;

	nn = clos->nbrAligSets + 1;

	for (x=0; x < clos->seqNbr; x++)	
	   {
		clos->aligSet[nn].pos[x] = 0;
		if (n1 > 0 && clos->aligSet[n1].pos[x] > 0)
			clos->aligSet[nn].pos[x] = clos->aligSet[n1].pos[x];
		else { if (n2 > 0 && clos->aligSet[n2].pos[x] > 0)
			clos->aligSet[nn].pos[x] = clos->aligSet[n2].pos[x];}

		if (clos->aligSet[nn].pos[x] == 0)
		    {
			clos->predFrontier[nn][x] = max(clos->gauche1[x], clos->gauche2[x]);
			clos->succFrontier[nn][x] = min(clos->droite1[x], clos->droite2[x]);
		    }
		else	clos->predFrontier[nn][x] = clos->succFrontier[nn][x] = clos->aligSet[nn].pos[x];
	   }
	clos->predFrontier[nn][seq1] = clos->succFrontier[nn][seq1] = clos->aligSet[nn].pos[seq1] = i;
	clos->predFrontier[nn][seq2] = clos->succFrontier[nn][seq2] = clos->aligSet[nn].pos[seq2] = j;


	for (x=clos->aligSet[nn].nbr=0; x < clos->seqNbr; x++)
	  if (clos->aligSet[nn].pos[x] > 0)
	    {
		k = clos->aligSet[nn].pos[x]; 
		clos->seq[x].aligSetNbr[k] = nn;
		clos->aligSet[nn].nbr++;
	    }

	for (x=0; x < clos->seqNbr; x++)
          if (clos->droite1[x] != clos->droite2[x])  /* => la front. clos->gauche peut changer */
	    for (y=0; y < clos->seqNbr; y++)
	      {
		clos->pos_[x][y] = 0;
		k = clos->succFrontier[nn][x];
		if (k == clos->aligSet[nn].pos[x]) 	
			k = clos->seq[x].succAligSetPos[k];
		if (k <= clos->seq[x].longueur)
		  while (k > 0)
		    {
			n = clos->seq[x].aligSetNbr[k];
			if (clos->predFrontier[n][y] < clos->predFrontier[nn][y])
			    {
				clos->pos_[x][y] = k;
				k = clos->seq[x].succAligSetPos[k];
			    }
			else	k = 0;
		    }
	      }

	for (x=0; x < clos->seqNbr; x++)
          if (clos->droite1[x] != clos->droite2[x])  
		/* => la front. gauche peut changer */
	    for (y=0; y < clos->seqNbr; y++)
	      {
		k = clos->succFrontier[nn][x];
		if (k == clos->aligSet[nn].pos[x]) 
			k = clos->seq[x].succAligSetPos[k];
		if (clos->pos_[x][y] > 0)
		  while (k > 0 && k <= clos->pos_[x][y])
		    {
			n = clos->seq[x].aligSetNbr[k];
			clos->predFrontier[n][y] = clos->predFrontier[nn][y];
			k = clos->seq[x].succAligSetPos[k];
		    }
	      }

	for (x=0; x < clos->seqNbr; x++)
          if (clos->gauche1[x] != clos->gauche2[x])  
		/* => la front. droite peut changer */
	    for (y=0; y < clos->seqNbr; y++)
	      {
		clos->pos_[x][y] = 0;
		k = clos->predFrontier[nn][x];
		if (k > 0 && k == clos->aligSet[nn].pos[x])
			k = clos->seq[x].predAligSetPos[k];
		while (k > 0)
		    {
			n = clos->seq[x].aligSetNbr[k];
			if (clos->succFrontier[n][y] > clos->succFrontier[nn][y])
			    {
				clos->pos_[x][y] = k;
				k = clos->seq[x].predAligSetPos[k];
			    }
			else	k = 0;
		    }
	      }

	for (x=0; x < clos->seqNbr; x++)
          if (clos->gauche1[x] != clos->gauche2[x])  /* => la front. clos->droite peut changer */
	    for (y=0; y < clos->seqNbr; y++)
	      {
		k = clos->predFrontier[nn][x];
		if (k > 0 && k == clos->aligSet[nn].pos[x])
			k = clos->seq[x].predAligSetPos[k];
		if (clos->pos_[x][y] > 0)
		  while (k >= clos->pos_[x][y])
		    {
			n = clos->seq[x].aligSetNbr[k];
			clos->succFrontier[n][y] = clos->succFrontier[nn][y];
			k = clos->seq[x].predAligSetPos[k];
		    }
	      }

	if (n1 == 0)
	   {
		for (k=i-1; k > 0 && clos->seq[seq1].aligSetNbr[k] == 0; k--)
			clos->seq[seq1].succAligSetPos[k] = i;
		if (k > 0) 	
			clos->seq[seq1].succAligSetPos[k] = i;
		for (k=i+1; k <= clos->seq[seq1].longueur 
					&& clos->seq[seq1].aligSetNbr[k] == 0; k++)
			clos->seq[seq1].predAligSetPos[k] = i;
		if (k <= clos->seq[seq1].longueur) 	
			clos->seq[seq1].predAligSetPos[k] = i;
	   }

	if (n2 == 0)
	   {
		for (k=j-1; k > 0 && clos->seq[seq2].aligSetNbr[k] == 0; k--)
			clos->seq[seq2].succAligSetPos[k] = j;
		if (k > 0) 	
			clos->seq[seq2].succAligSetPos[k] = j;
		for (k=j+1; k <= clos->seq[seq2].longueur 
					&& clos->seq[seq2].aligSetNbr[k] == 0; k++)
			clos->seq[seq2].predAligSetPos[k] = j;
		if (k <= clos->seq[seq2].longueur) 	
			clos->seq[seq2].predAligSetPos[k] = j;
	   }


	if (n1 > n2)  {	n = n1; n1 = n2; n2 = n; }

	if (n2 == 0)	
	     {
		clos->nbrAligSets++;

		realloc_closure(clos);
	     }
	else { 
		if (n1 == 0)
		     {
			moveAligSet(clos, n2, nn);
		     }
		else  
		     {
			moveAligSet(clos, n1, nn);

			if (n2 < clos->nbrAligSets)	moveAligSet(clos, n2, clos->nbrAligSets);
			clos->nbrAligSets--;

			realloc_closure(clos);
		     }
	     }
      }
}

int path(CLOSURE *clos, int x, int i, int y, int j)
{
int n2, k;

if (x == y) return(i <= j);

n2 = clos->seq[y].aligSetNbr[j];

if (n2 == 0) 
    {
	k = clos->seq[y].predAligSetPos[j];
	if (k > 0) n2 = clos->seq[y].aligSetNbr[k];
    }

if (n2 == 0)	return(false);
else		return(i <= clos->predFrontier[n2][x]);
}

int alignedPositions(CLOSURE *clos, int x, int i, int y, int j)
{

return (x == y && i == j) || (clos->seq[x].aligSetNbr[i] != 0 &&
	clos->seq[x].aligSetNbr[i] == clos->seq[y].aligSetNbr[j]);
}

int alignablePositions(CLOSURE *clos, int x, int i, int y, int j)
{

if (path(clos, x, i, y, j))
	return(path(clos, y, j, x, i));
else		
	return(!path(clos, y, j, x, i));
}

int addAlignedSegments(CLOSURE *clos, int x, int i, int y, int j, int l)
{
int k;

for (k=0; k < l; i++, j++, k++)
	addAlignedPositions(clos, x, i, y, j);

}

int alignableSegments(CLOSURE *clos, int x, int i, int y, int j, int l)
{
int k;

for (k=0; k < l && alignablePositions(clos, x, i, y, j); i++, j++, k++);

return(k==l);
}

int alignedSegments(CLOSURE *clos, int x, int i, int y, int j, int l)
{
int k;

for (k=0; k < l && alignedPositions(clos, x, i, y, j); i++, j++, k++);

return(k==l);
}

int predFrontier(CLOSURE *clos, int x, int i, int y) /* on suppose que x!=y */
{
int n, k;

n = clos->seq[x].aligSetNbr[i];

if (n == 0)
   {
	k = clos->seq[x].predAligSetPos[i];
	if (k > 0) n = clos->seq[x].aligSetNbr[k];
   }

if (n > 0)	return(clos->predFrontier[n][y]);
else		return(0);
}

int succFrontier(CLOSURE *clos, int x, int i, int y) /* on suppose que x!=y */
{
int n, k;

n = clos->seq[x].aligSetNbr[i];

if (n == 0)
   {
	k = clos->seq[x].succAligSetPos[i];
	if (k > 0) n = clos->seq[x].aligSetNbr[k];
   }

if (n > 0)	return(clos->succFrontier[n][y]);
else		return(clos->seq[y].longueur+1);
}
