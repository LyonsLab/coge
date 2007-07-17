/************************************************************/
/************************************************************/
/**             GABIOS-LIB 1.0 (1999)                      **/
/** A library for Greedy Alignment of BIOlogical Sequences **/
/**             Developed by Said Abdeddaim                **/
/**          Said.Abdeddaim@dir.univ-rouen.fr              **/
/************************************************************/
/************************************************************/                         

#ifndef _ALIG_GRAPH_CLOSURE_H
#define _ALIG_GRAPH_CLOSURE_H


typedef struct {
	int *pos;
	int nbr;
	} positionSet;

typedef struct {
		int longueur;

		int *aligSetNbr, *predAligSetPos, *succAligSetPos;
	      } sequence;

typedef struct {
		int seqNbr;
		sequence *seq;
		int maxLong;

		positionSet *aligSet;
		int nbrAligSets, oldNbrAligSets;

		int **predFrontier, **succFrontier;

		int *topolog;
		int *gauche1, *gauche2, *droite1, *droite2, **pos_;

		} CLOSURE;
		

CLOSURE *newAligGraphClosure(int nbreseq, int *longseq, 
				int nbreancr, int **ancrages);

void freeAligGraphClosure(CLOSURE *clos);

int addAlignedPositions(CLOSURE *clos, int x, int i, int y, int j);

int alignablePositions(CLOSURE *clos, int x, int i, int y, int j);

int alignedPositions(CLOSURE *clos, int x, int i, int y, int j);

int addAlignedSegments(CLOSURE *clos, int x, int i, int y, int j, int l);

int alignableSegments(CLOSURE *clos, int x, int i, int y, int j, int l);

int alignedSegments(CLOSURE *clos, int x, int i, int y, int j, int l);

int predFrontier(CLOSURE *clos, int x, int i, int y);

int succFrontier(CLOSURE *clos, int x, int i, int y);


#endif /* _ALIG_GRAPH_CLOSURE_H */
