/* version 3.6. (c) Copyright 1993-2005 by the University of Washington.
   Written by Dan Fineman, Joseph Felsenstein, Mike Palczewski, Hisashi Horino,
   Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee
   is charged for it and provided that this copyright notice is not removed. */

#include "phylip.h"
#include "cons.h"


typedef enum { SYMMETRIC, BSD } distance_type;

/* The following extern's refer to things declared in cons.c */
extern int tree_pairing;
extern Char outfilename[FNMLNGTH], intreename[FNMLNGTH], intree2name[FNMLNGTH], outtreename[FNMLNGTH];
extern node *root;

extern long numopts, outgrno, col;
extern long maxgrp;    /* max. no. of groups in all trees found  */

extern boolean trout, firsttree, noroot, outgropt, didreroot, prntsets,
               progress, treeprint, goteof;
extern pointarray treenode, nodep;
extern group_type **grouping, **grping2, **group2;/* to store groups found  */
extern long **order, **order2, lasti;
extern group_type *fullset;
extern node *grbg;
extern long tipy;

extern double **timesseen, **tmseen2, **times2;
extern double trweight, ntrees;

static distance_type  dtype;
static long output_scheme;

#ifndef OLDC
/* function prototpes */
void    assign_tree(group_type **, pattern_elm ***, long, long *); 
boolean group_is_null(group_type **, long);
void    compute_distances(pattern_elm ***, long, long);
void    free_patterns(pattern_elm ***, long); 
void    produce_square_matrix(long, long *);
void    produce_full_matrix(long, long, long *);
void    output_submenu(void);
void    pairing_submenu(void);
void    read_second_file(pattern_elm ***, long, long);
void    getoptions(void);
void    assign_lengths(double **lengths, pattern_elm ***pattern_array, 
             long tree_index);
void    print_header(long trees_in_1, long trees_in_2);
void    output_distances(long trees_in_1, long trees_in_2);
void    output_long_distance(long diffl, long tree1, long tree2,
                             long trees_in_1, long trees_in_2);
void    output_matrix_long(long diffl, long tree1, long tree2, long trees_in_1,
                long trees_in_2);
void    output_matrix_double(double diffl, long tree1, long tree2, long trees_in_1,
                long trees_in_2);
void    output_double_distance(double diffd, long tree1, long tree2,
                long trees_in_1, long trees_in_2);
long    symetric_diff(group_type **tree1, group_type **tree2,
                long ntree1, long ntree2, long patternsz1, long patternsz2);
double  bsd_tree_diff(group_type **tree1, group_type **tree2, long ntree1,
                long ntree2, double* lengths1, double *lengths2,
                long patternsz1, long patternsz2);
void    tree_diff(group_type **tree1, group_type **tree2, double *lengths1,
                double* lengths2, long patternsz1, long patternsz2, long ntree1,
                long ntree2, long trees_in_1, long trees_in_2);
void    print_line_heading(long tree);
int     get_num_columns(void);
void    print_matrix_heading(long tree, long maxtree);
/* function prototpes */
#endif


void assign_lengths(double **lengths, pattern_elm ***pattern_array, long
                        tree_index)
{
  *lengths = pattern_array[0][tree_index]->length;
}


void assign_tree(group_type **treeN, pattern_elm ***pattern_array,
                long tree_index, long *pattern_size) 
{ /* set treeN to be the tree_index-th tree in pattern_elm */
  long i;

  for ( i = 0 ; i < setsz ; i++ ) {
    treeN[i] = pattern_array[i][tree_index]->apattern;
  }
  *pattern_size = *pattern_array[0][tree_index]->patternsize;
}  /* assign_tree */


boolean group_is_null(group_type **treeN, long index)
{
  /* Check to see if a given index to a tree array points to an empty
     group */
  long i;

  for ( i = 0 ; i < setsz ; i++ )
    if (treeN[i][index] != (group_type) 0)
      return false;

  /* If we've gotten this far, then the index is to an empty group in
     the tree. */
  return true;
}  /* group_is_null */


double bsd_tree_diff(group_type **tree1, group_type **tree2,
                     long ntree1, long ntree2, double *lengths1,
                     double* lengths2, long patternsz1, long patternsz2)
{
  /* Compute the difference between 2 given trees. Return
     that value as a double. */

  long index1, index2;
  double return_value = 0;
  boolean match_found;
  long i;

  if ( group_is_null(tree1, 0) || group_is_null(tree2, 0) ) {
    printf ("Error computing tree difference between tree %ld and tree %ld\n",
             ntree1, ntree2);
    exxit(-1);
  }

  for ( index1 = 0; index1 < patternsz1; index1++ ) {
    if ( !group_is_null(tree1, index1) ) {
      if ( lengths1[index1] == -1 ) {
        printf(
          "Error: tree %ld is missing a length from at least one branch\n",
          ntree1
        );

        exxit(-1);
      }
    }
  }

  for ( index2 = 0; index2 < patternsz2; index2++ ) {
    if ( !group_is_null(tree2, index2) ) {
      if ( lengths2[index2] == -1 ) {
        printf(
          "Error: tree %ld is missing a length from at least one branch\n",
          ntree2
        );
        exxit(-1);
      }
    }
  }

  for ( index1 = 0 ; index1 < patternsz1; index1++ ) {
    /* For every element in the first tree, see if there's
        a match to it in the second tree. */
    match_found = false;
    
    if ( group_is_null(tree1, index1) ) {
      /* When we've gone over all the elements in tree1, greater
          number of elements in tree2 will constitute that much more
          of a difference... */
      while ( !group_is_null(tree2, index1) ) {
        return_value += pow(lengths1[index1], 2);
        index1++;
      }
      break;
    }

    for ( index2 = 0 ; index2 < patternsz2 ; index2++ ) {
      /* For every element in the second tree, see if any match
          the current element in the first tree. */
      if ( group_is_null(tree2, index2) ) {
        /* When we've gone over all the elements in tree2 */
        match_found = false;
        break;
      }
      else {
        /* Tentatively set match_found; will be changed later if
            neccessary. . . */
        match_found = true;  

        for ( i = 0 ; i < setsz ; i++ ) {
            /* See if we've got a match, */ 
            if ( tree1[i][index1] != tree2[i][index2] )
              match_found = false;
        }

        if ( match_found == true ) {
          break;
        }
      }
    }

    if ( match_found == false ) {
        return_value += pow(lengths1[index1], 2);
    }
  }

  for ( index2 = 0 ; index2 < patternsz2 ; index2++ ) {
    /* For every element in the second tree, see if there's
        a match to it in the first tree. */
    match_found = false;
    if ( group_is_null(tree2, index2) ) {
      /* When we've gone over all the elements in tree2, greater
          number of elements in tree1 will constitute that much more
          of a difference... */

      while ( !group_is_null(tree1, index2) ) {
        return_value += pow(lengths2[index2], 2);
        index2++;
      }
      break;
    }

    for ( index1 = 0 ; index1 < patternsz1 ; index1++ ) {
      /* For every element in the first tree, see if any match
          the current element in the second tree. */
      if ( group_is_null (tree1, index1) ) {
        /* When we've gone over all the elements in tree2 */
        match_found = false;
        break;
      }
      else {
        /* Tentatively set match_found; will be changed later if
            neccessary. . . */
        match_found = true;  

        for ( i = 0 ; i < setsz ; i++ ) {
            /* See if we've got a match, */ 
            if ( tree2[i][index2] != tree1[i][index1] )
              match_found = false;
        }

        if ( match_found == true ) {
          return_value += pow(lengths1[index1] - lengths2[index2], 2);
          break;
        }
      }
    }

    if ( match_found == false ) {
        return_value += pow(lengths2[index2], 2);
    }
  }
  if ( return_value > 0.0 )
    return_value = sqrt(return_value);
  else
    return_value = 0.0;
  return return_value;
}


long symetric_diff(group_type **tree1, group_type **tree2,
                long ntree1, long ntree2, long patternsz1, long patternsz2)
{
  /* Compute the symmetric difference between 2 given trees. Return
     that value as a long. */

  long index1, index2, return_value = 0;
  boolean match_found;
  long i;

  if ( group_is_null(tree1, 0) || group_is_null(tree2, 0) ) {
    printf ("Error computing tree difference.\n");
    return 0;
  }

  for ( index1 = 0 ; index1 < patternsz1 ; index1++ ) {
    /* For every element in the first tree, see if there's
        a match to it in the second tree. */
    match_found = false;
    if ( group_is_null (tree1, index1) ) {
      /* When we've gone over all the elements in tree1, greater
          number of elements in tree2 will constitute that much more
          of a difference... */

      while ( !group_is_null(tree2, index1) ) {
        return_value++;
        index1++;
      }
      break;
    }

    for ( index2 = 0 ; index2 < patternsz2 ; index2++ ) {
      /* For every element in the second tree, see if any match
          the current element in the first tree. */
      if ( group_is_null(tree2, index2) ) {
        /* When we've gone over all the elements in tree2 */
        match_found = false;
        break;
      }
      else {
        /* Tentatively set match_found; will be changed later if
            neccessary. . . */
        match_found = true;  

        for ( i = 0 ; i < setsz ; i++ ) {
            /* See if we've got a match, */ 
            if ( tree1[i][index1] != tree2[i][index2] )
              match_found = false;
        }

        if ( match_found == true ) {
          /* If the previous loop ran from 0 to setsz without setting
              match_found to false, */
          break;
        }
      }
    }

    if ( match_found == false ) {
      return_value++;
    }
  }
  return return_value;
}  /* symetric_diff */


void output_double_distance(double diffd, long tree1, long tree2,
    long trees_in_1, long trees_in_2)
{
  switch ( tree_pairing ) {
    case ADJACENT_PAIRS: 
      if ( output_scheme == VERBOSE ) {
        fprintf (outfile, "Trees %ld and %ld:    %e\n", tree1, tree2, diffd);
      }
      else if (output_scheme == SPARSE) {
        fprintf (outfile, "%ld %ld %e\n", tree1, tree2, diffd);
      }
      break;

    case ALL_IN_FIRST: 
      if ( output_scheme == VERBOSE ) {
        fprintf (outfile, "Trees %ld and %ld:    %e\n", tree1, tree2, diffd);
      }
      else if ( output_scheme == SPARSE ) {
        fprintf (outfile, "%ld %ld %e\n", tree1, tree2, diffd );
      }
      else if ( output_scheme == FULL_MATRIX ) {
        output_matrix_double(diffd, tree1, tree2, trees_in_1, trees_in_2);
      }
      break;

    case CORR_IN_1_AND_2:
      if ( output_scheme == VERBOSE ) {
        fprintf(outfile, "Tree pair %ld:    %e\n", tree1, diffd);
      }
      else if ( output_scheme == SPARSE ) {
        fprintf(outfile, "%ld %e\n", tree1, diffd);
      }
      break; 

    case ALL_IN_1_AND_2:
      if ( output_scheme == VERBOSE )    
        fprintf(outfile, "Trees %ld and %ld:    %e\n", tree1, tree2, diffd);
      else if ( output_scheme == SPARSE )    
        fprintf(outfile, "%ld %ld %e\n", tree1, tree2, diffd);
      else if ( output_scheme == FULL_MATRIX ) {
        output_matrix_double(diffd, tree1, tree2, trees_in_1, trees_in_2);
      }
      break; 
  }
} /* output_double_distance */


void print_matrix_heading(long tree, long maxtree)
{
  long i;

  if ( tree_pairing == ALL_IN_1_AND_2 ) {
    fprintf(outfile, "\n\nFirst\\  Second tree file:\n");
    fprintf(outfile, "tree  \\\n");
    fprintf(outfile, "file:  \\");
  }
  else
    fprintf(outfile, "\n\n      ");
      
  for ( i = tree ; i <= maxtree ; i++ ) {
    if ( dtype == SYMMETRIC ) 
      fprintf(outfile, "%5ld ", i);
    else 
      fprintf(outfile, "    %7ld ", i);
  }
  fprintf(outfile, "\n");
  if ( tree_pairing == ALL_IN_1_AND_2 ) 
    fprintf(outfile, "        \\");
  else
    fprintf(outfile, "      \\");
  for ( i = tree ;  i <= maxtree ; i++ ) {
    if ( dtype == SYMMETRIC )
      fprintf(outfile, "------");
    else
      fprintf(outfile, "------------");
  }
}


void print_line_heading(long tree) 
{
  if ( tree_pairing == ALL_IN_1_AND_2 ) 
    fprintf(outfile, "\n%4ld    |", tree);
  else
    fprintf(outfile, "\n%5ld |", tree);
}


void output_matrix_double(double diffl, long tree1, long tree2, long trees_in_1,
    long trees_in_2)
{
  if ( tree1 == 1 && ((tree2 - 1) % get_num_columns() == 0 || tree2 == 1 ) ) {
    if ( (tree_pairing == ALL_IN_FIRST
           && tree2 + get_num_columns() - 1 < trees_in_1
         )
         || (tree_pairing == ALL_IN_1_AND_2
           && tree2 + get_num_columns() - 1 < trees_in_2
         )
       ) {
      print_matrix_heading(tree2, tree2 + get_num_columns() - 1);
    }
    else {
      if ( tree_pairing == ALL_IN_FIRST )
        print_matrix_heading(tree2, trees_in_1);
      else
        print_matrix_heading(tree2, trees_in_2);
    }
  }
  if ( (tree2 - 1) % get_num_columns() == 0 || tree2 == 1 ) {
    print_line_heading(tree1);
  }
  fprintf(outfile, " %9g  ", diffl);
  if ((tree_pairing == ALL_IN_FIRST && tree1 == trees_in_1 && tree2 == 
      trees_in_1) || (tree_pairing == ALL_IN_1_AND_2 && tree1 == trees_in_1 &&
        tree2 == trees_in_2))
    fprintf(outfile, "\n\n\n");
} /* output_matrix_double */


void output_matrix_long(long diffl, long tree1, long tree2, long trees_in_1,
    long trees_in_2)
{
  if ( tree1 == 1 && ((tree2 - 1) % get_num_columns() == 0 || tree2 == 1 )) {
    if ( (tree_pairing == ALL_IN_FIRST && tree2 + get_num_columns() - 1
          < trees_in_1) ||
         (tree_pairing == ALL_IN_1_AND_2 && tree2 + get_num_columns() - 1
          < trees_in_2)) {
      print_matrix_heading(tree2, tree2 + get_num_columns() - 1);
    } else {
      if ( tree_pairing == ALL_IN_FIRST)
        print_matrix_heading(tree2, trees_in_1);
      else
        print_matrix_heading(tree2, trees_in_2);
    }
  }
  if ( (tree2 - 1) % get_num_columns() == 0 || tree2 == 1) {
    print_line_heading(tree1);
  }
  fprintf(outfile, "%4ld  ", diffl);
  if ((tree_pairing == ALL_IN_FIRST && tree1 == trees_in_1 && tree2 == 
      trees_in_1) || (tree_pairing == ALL_IN_1_AND_2 && tree1 == trees_in_1 &&
        tree2 == trees_in_2))
    fprintf(outfile, "\n\n\n");
} /* output_matrix_long */


void output_long_distance(long diffl, long tree1, long tree2, long trees_in_1,
    long trees_in_2)
{
  switch (tree_pairing) {
    case ADJACENT_PAIRS: 
      if (output_scheme == VERBOSE ) {
        fprintf (outfile, "Trees %ld and %ld:    %ld\n", tree1, tree2, diffl);
      } else if (output_scheme == SPARSE) {
        fprintf (outfile, "%ld %ld %ld\n", tree1, tree2, diffl);
      }
      break;

    case ALL_IN_FIRST: 
      if (output_scheme == VERBOSE) {
        fprintf (outfile, "Trees %ld and %ld:    %ld\n", tree1, tree2, diffl);
      } else if (output_scheme == SPARSE) {
        fprintf (outfile, "%ld %ld %ld\n", tree1, tree2, diffl );
      } else if (output_scheme == FULL_MATRIX) {
        output_matrix_long(diffl, tree1, tree2, trees_in_1, trees_in_2); 
      }
      break;

    case CORR_IN_1_AND_2:

      if (output_scheme == VERBOSE) {
        fprintf (outfile, "Tree pair %ld:    %ld\n", tree1, diffl);
      } else if (output_scheme == SPARSE) {
        fprintf (outfile, "%ld %ld\n", tree1, diffl);
      }
      break; 

    case ALL_IN_1_AND_2:
      if (output_scheme == VERBOSE)
        fprintf (outfile, "Trees %ld and %ld:    %ld\n", tree1, tree2, diffl);
      else if (output_scheme == SPARSE)
        fprintf (outfile, "%ld %ld %ld\n", tree1, tree2, diffl);
      else if (output_scheme == FULL_MATRIX ) {
        output_matrix_long(diffl, tree1, tree2, trees_in_1, trees_in_2); 
      }
      break;
  }
}


void tree_diff(group_type **tree1, group_type **tree2, double *lengths1,
                double* lengths2, long patternsz1, long patternsz2,
                long ntree1, long ntree2, long trees_in_1, long trees_in_2)
{
  long diffl;
  double diffd;

  switch (dtype) {
    case SYMMETRIC:
      diffl = symetric_diff (tree1, tree2, ntree1, ntree2,
                              patternsz1, patternsz2);
      diffl += symetric_diff (tree2, tree1, ntree1, ntree2,
                               patternsz2, patternsz1);
      output_long_distance(diffl, ntree1, ntree2, trees_in_1, trees_in_2);
      break;
    case BSD:
      diffd = bsd_tree_diff(tree1, tree2, ntree1, ntree2,
                             lengths1, lengths2, patternsz1, patternsz2);
      output_double_distance(diffd, ntree1, ntree2, trees_in_1, trees_in_2);
      break;
  }
} /* tree_diff */


int get_num_columns(void) 
{
  if ( dtype == SYMMETRIC )
    return 10;
  else return 7;
} /* get_num_columns */


void compute_distances(pattern_elm ***pattern_array, long trees_in_1,
                long trees_in_2)
{
  /* Compute symmetric distances between arrays of trees */

  long  tree_index, end_tree, index1, index2, diff_index, index3;
  group_type **treeA, **treeB;
  long patternsz1, patternsz2;
  double *length1 = NULL, *length2 = NULL; 
  int num_columns = get_num_columns();

  diff_index = 0;
  index1 = 0;

  /* Put together space for treeA and treeB */
  treeA = (group_type **) Malloc (setsz * sizeof (group_type *));
  treeB = (group_type **) Malloc (setsz * sizeof (group_type *));

  print_header(trees_in_1, trees_in_2);

  switch (tree_pairing) {

    case ADJACENT_PAIRS: 
      /* For every tree, compute the distance between it and the tree
          at the next location; do this in both directions */
      end_tree = trees_in_1 - 1;
      for (tree_index = 0 ; tree_index < end_tree ; tree_index += 2) {

        assign_tree (treeA, pattern_array, tree_index, &patternsz1);
        assign_tree (treeB, pattern_array, tree_index + 1, &patternsz2);
        assign_lengths(&length1, pattern_array, tree_index);
        assign_lengths(&length2, pattern_array, tree_index + 1);

        tree_diff (treeA, treeB, length1, length2, patternsz1, patternsz2,
            tree_index+1, tree_index+2, trees_in_1, trees_in_2);

        if (tree_index + 2 == end_tree)
          printf("\nWARNING: extra tree at the end of input tree file.\n");
      }
      break;
  
    case ALL_IN_FIRST: 
      /* For every tree, compute the distance between it and every
          other tree in that file. */
      end_tree   = trees_in_1;
      if ( output_scheme != FULL_MATRIX ) {
        /* verbose or sparse output */
        for (index1 = 0 ; index1 < end_tree ; index1++) {
          assign_tree (treeA, pattern_array, index1, &patternsz1);
          assign_lengths(&length1, pattern_array, index1);
  
          for (index2 = 0 ; index2 < end_tree ; index2++) {
            assign_tree (treeB, pattern_array, index2, &patternsz2);
            assign_lengths(&length2, pattern_array, index2);
            tree_diff (treeA, treeB, length1, length2, patternsz1, patternsz2,
                index1 + 1, index2 + 1, trees_in_1, trees_in_2);
          }
        }
      }
      else {
        /* full matrix output */
        for ( index3 = 0 ; index3 < trees_in_1 ; index3 += num_columns) {
          for ( index1 = 0 ; index1 < trees_in_1 ; index1++) {
          assign_tree (treeA, pattern_array, index1, &patternsz1);
          assign_lengths(&length1, pattern_array, index1);
            for ( index2 = index3 ; 
                  index2 < index3 + num_columns && index2 < trees_in_1 ;
                  index2++
                ) {
              assign_tree (treeB, pattern_array, index2, &patternsz2);
              assign_lengths(&length2, pattern_array, index2);
              tree_diff (treeA, treeB, length1, length2, patternsz1, patternsz2,
                  index1 + 1, index2 + 1, trees_in_1, trees_in_2);
            }
          }
        }
      }
      break;

    case CORR_IN_1_AND_2:
      if (trees_in_1 != trees_in_2) {
        /* Set end tree to the smaller of the two totals. */
        end_tree = trees_in_1 > trees_in_2 ? trees_in_2 : trees_in_1;
	
        /* Print something out to the outfile and to the terminal. */
        fprintf(outfile,
            "\n\n"
            "*** Warning: differing number of trees in first and second\n"
            "*** tree files.  Only computing %ld pairs.\n"
            "\n",
            end_tree
        );
        printf(
            "\n"
            " *** Warning: differing number of trees in first and second\n"
            " *** tree files.  Only computing %ld pairs.\n"
            "\n",
            end_tree
        );
  
      }
      else
        end_tree = trees_in_1;
  
      for (tree_index = 0 ; tree_index < end_tree ; tree_index++) {
      /* For every tree, compute the distance between it and the
            tree at the parallel location in the other file; do this in
            both directions */
  
        assign_tree(treeA, pattern_array, tree_index, &patternsz1);
        assign_lengths(&length1, pattern_array, tree_index);
        /* (tree_index + trees_in_1) will be the corresponding tree in
            the second file. */
        assign_tree(treeB, pattern_array, tree_index + trees_in_1, &patternsz2);
        assign_lengths(&length2, pattern_array, tree_index + trees_in_1);
        tree_diff(
          treeA, treeB, length1, length2, patternsz1, patternsz2,
          tree_index + 1, 0, trees_in_1, trees_in_2
        );
      }
      break; 

    case ALL_IN_1_AND_2:
      end_tree = trees_in_1 + trees_in_2;
      
      if ( output_scheme != FULL_MATRIX ) {
        for (tree_index = 0 ; tree_index < trees_in_1 ; tree_index++) {
          /* For every tree in the first file, compute the distance
              between it and every tree in the second file. */
          assign_tree (treeA, pattern_array, tree_index, &patternsz1);
          assign_lengths(&length1, pattern_array, tree_index);
          for (index2 = trees_in_1 ; index2 < end_tree ; index2++) {
            assign_tree (treeB, pattern_array, index2, &patternsz2);
            assign_lengths(&length2, pattern_array, index2);
            tree_diff(treeA, treeB, length1, length2, patternsz1, patternsz2,
                tree_index + 1 , index2 + 1, trees_in_1, trees_in_2);
          }
        }
        for ( ; tree_index < end_tree ; tree_index++) {
          /* For every tree in the second file, compute the distance
              between it and every tree in the first file. */
  
          assign_tree (treeA, pattern_array, tree_index, &patternsz1);
          assign_lengths(&length1, pattern_array, tree_index);
  
          for (index2 = 0 ; index2 < trees_in_1 ; index2++) {
            assign_tree (treeB, pattern_array, index2, &patternsz2);
            assign_lengths(&length2, pattern_array, index2);
            tree_diff (treeA, treeB, length1, length2 , patternsz1, patternsz2,
                tree_index + 1, index2 + 1, trees_in_1, trees_in_2);
          }
        }
      } else {
        for ( index3 = trees_in_1 ; index3 < end_tree ; index3 += num_columns) {
          for ( index1 = 0 ; index1 < trees_in_1 ; index1++) {
          assign_tree (treeA, pattern_array, index1, &patternsz1);
          assign_lengths(&length1, pattern_array, index1);
            for ( index2 = index3 ; 
                index2 < index3 + num_columns && index2 < end_tree ;
                index2++) {
              assign_tree (treeB, pattern_array, index2, &patternsz2);
              assign_lengths(&length2, pattern_array, index2);
              tree_diff (treeA, treeB, length1, length2, patternsz1, patternsz2,
                  index1 + 1, index2 - trees_in_1 + 1, trees_in_1, trees_in_2);
            }
          }
        }
      }
      break; 
  }
  /* Free up treeA and treeB */
  free (treeA);
  free (treeB);        
}  /* compute_distances */


void free_patterns(pattern_elm ***pattern_array, long total_trees) 
{
  long i, j;

  /* Free each pattern array, */
  for (i=0 ; i < setsz ; i++) {
    for (j = 0 ; j < total_trees ; j++) {
      free (pattern_array[i][j]->apattern);
      free (pattern_array[i][j]->patternsize);
      free (pattern_array[i][j]->length);
      free (pattern_array[i][j]);
    }
    free (pattern_array[i]);
  }
  free (pattern_array);
}  /* free_patterns */


void print_header(long trees_in_1, long trees_in_2) 
{
  long end_tree;

  switch (tree_pairing) {
    case ADJACENT_PAIRS: 
      end_tree = trees_in_1 - 1;

      if (output_scheme == VERBOSE) {
        fprintf(outfile,
            "\n"
            "Tree distance program, version %s\n\n", VERSION
               );
        if (dtype == BSD)
          fprintf(outfile, 
            "Branch score distances between adjacent pairs of trees:\n"
            "\n"
                 );
        else
          fprintf (outfile, 
            "Symmetric differences between adjacent pairs of trees:\n\n");
      }
      else if ( output_scheme != SPARSE)
        printf ("Error -- cannot output adjacent pairs into a full matrix.\n");
      break;

    case ALL_IN_FIRST: 
      end_tree = trees_in_1;

      if (output_scheme == VERBOSE) {
        fprintf(outfile, "\nTree distance program, version %s\n\n", VERSION);
        if (dtype == BSD)
          fprintf (outfile, 
        "Branch score distances between all pairs of trees in tree file\n\n"
                  );
        else
          fprintf (outfile, 
        "Symmetric differences between all pairs of trees in tree file:\n\n");
      }
      else if (output_scheme == FULL_MATRIX) {
        fprintf(outfile, "\nTree distance program, version %s\n\n", VERSION);
        if (dtype == BSD)
          fprintf (outfile, 
        "Branch score distances between all pairs of trees in tree file:\n\n");
        else
          fprintf (outfile, 
          "Symmetric differences between all pairs of trees in tree file:\n\n");
      }
      break;

    case CORR_IN_1_AND_2:

      if (output_scheme == VERBOSE) {
        fprintf(outfile, "\nTree distance program, version %s\n\n", VERSION);
      if (dtype == BSD) {
          fprintf (outfile, 
            "Branch score distances between corresponding pairs of trees\n");
          fprintf (outfile, "   from first and second tree files:\n\n");
        }
        else {
          fprintf (outfile, 
              "Symmetric differences between corresponding pairs of trees\n");
          fprintf (outfile, "   from first and second tree files:\n\n");
        }
      }
      else if (output_scheme != SPARSE)
        printf (
          "Error -- cannot output corresponding pairs into a full matrix.\n");
      break;

    case (ALL_IN_1_AND_2) :
      if ( output_scheme == VERBOSE) {
        fprintf(outfile, "\nTree distance program, version %s\n\n", VERSION);
        if (dtype == BSD) {
          fprintf (outfile, 
            "Branch score distances between all pairs of trees\n");
          fprintf (outfile, "   from first and second tree files:\n\n");
        }
        else {
         fprintf(outfile,"Symmetric differences between all pairs of trees\n");
          fprintf(outfile,"   from first and second tree files:\n\n");
        }
      } else if ( output_scheme == FULL_MATRIX) {
        fprintf(outfile, "\nTree distance program, version %s\n\n", VERSION);
      }
      break;
  }
} /* print_header */


void output_submenu()
{
  /* this allows the user to select a different output of distances scheme. */
  long loopcount;
  boolean done = false;
  Char    ch;

  if (tree_pairing == NO_PAIRING)
    return;

  loopcount = 0;
  while (!done) {
    printf ("\nDistances output options:\n");

    if ((tree_pairing == ALL_IN_1_AND_2) ||
          (tree_pairing == ALL_IN_FIRST))
      printf (" F     Full matrix.\n");
    
    printf (" V     One pair per line, verbose.\n");
    printf (" S     One pair per line, sparse.\n");
    
    if ((tree_pairing == ALL_IN_1_AND_2) ||
	(tree_pairing == ALL_IN_FIRST))
      printf ("\n Choose one: (F,V,S)\n");
    else
      printf ("\n Choose one: (V,S)\n");

    fflush(stdout);
    scanf("%c%*[^\n]", &ch);
    getchar();
    uppercase(&ch);

    if (strchr("FVS", ch) != NULL) {
      switch (ch) {
        case 'F':
          if ((tree_pairing == ALL_IN_1_AND_2) ||
              (tree_pairing == ALL_IN_FIRST))
            output_scheme = FULL_MATRIX;
          else
            /* If this can't be a full matrix... */
            continue;
          break;

        case 'V':
          output_scheme = VERBOSE;
          break;

        case 'S':
          output_scheme = SPARSE;
          break;
      }
      done = true;
    }
    countup(&loopcount, 10);
  }
}  /* output_submenu */


void pairing_submenu()
{
  /* this allows the user to select a different tree pairing scheme. */
  long    loopcount;
  boolean done = false;
  Char    ch;

  loopcount = 0;
  while (!done) {
    cleerhome();
    printf(
      "Tree Pairing Submenu:\n"
      " A     Distances between adjacent pairs in tree file.\n"
      " P     Distances between all possible pairs in tree file.\n"
      " C     Distances between corresponding pairs in one tree file and another.\n"
      " L     Distances between all pairs in one tree file and another.\n"
      "\n"
      " Choose one: (A,P,C,L)\n"
    );
      
    fflush(stdout);
    scanf("%c%*[^\n]", &ch);
    getchar();
    uppercase(&ch);

    if ( strchr("APCL", ch) != NULL ) {
      switch ( ch ) {
        case 'A':
          tree_pairing = ADJACENT_PAIRS;
          break;
              
        case 'P':
          tree_pairing = ALL_IN_FIRST;
          break;

        case 'C':
          tree_pairing = CORR_IN_1_AND_2;
          break;
              
        case 'L':
          tree_pairing = ALL_IN_1_AND_2;
          break;
      }
      output_submenu();
      done = true;
    }
    countup(&loopcount, 10);
  }
}  /* pairing_submenu */


void read_second_file(
    pattern_elm ***pattern_array,
    long trees_in_1,
    long trees_in_2
)
{
  boolean firsttree2, haslengths;
  long nextnode, trees_read=0;
  long j;

  firsttree2 = false;
  while ( !eoff(intree2) ) {
    goteof = false;
    nextnode = 0;
    haslengths = false;
    allocate_nodep(&nodep, &intree2, &spp);
    treeread(
        intree2, &root, treenode, &goteof, &firsttree2, nodep, &nextnode,
        &haslengths, &grbg, initconsnode, false, -1
    );
    missingname(root);
    reordertips();
    if (goteof)
      continue;
    ntrees += trweight;
    if (noroot) {
      reroot(nodep[outgrno - 1], &nextnode);
      didreroot = outgropt;
    }
    accumulate(root);
    gdispose(root);
    for (j = 0; j < 2*(1 + spp); j++)
      nodep[j] = NULL;
    free(nodep);

    store_pattern(pattern_array, trees_in_1 + trees_read);
    trees_read++;
  }
}  /* read_second_file */


void getoptions()
{
  /* interactively set options */
  long loopcount;
  Char ch;
  boolean done;

  /* Initial settings */
  dtype          = BSD;
  tree_pairing   = ADJACENT_PAIRS;
  output_scheme  = VERBOSE;
  ibmpc          = IBMCRT;
  ansi           = ANSICRT;
  didreroot      = false;
  spp            = 0;
  grbg           = NULL;
  col            = 0;

  putchar('\n');
  noroot = true;
  numopts = 0;
  outgrno = 1;
  outgropt = false;
  progress = true;

  /* The following are not used by treedist, but may be used
     in functions in cons.c, so we set them here. */
  treeprint = false;
  trout = false;
  prntsets = false;

  loopcount = 0;
  do {
    cleerhome();
    printf("\nTree distance program, version %s\n\n", VERSION);
    printf("Settings for this run:\n");
   
    printf(" D                         Distance Type:  ");
    switch (dtype) {
      case SYMMETRIC:
        printf("Symmetric Difference\n");
        break;
      case BSD:
        printf("Branch Score Distance\n");
        break;
    }
    printf(" R         Trees to be treated as Rooted:");
    if (noroot)
      printf("  No\n");
    else
      printf("  Yes\n");
    printf(" T    Terminal type (IBM PC, ANSI, none):");
    if (ibmpc)
      printf("  IBM PC\n");
    if (ansi)
      printf("  ANSI\n");
    if (!(ibmpc || ansi))
      printf("  (none)\n");
    printf(" 1  Print indications of progress of run:  %s\n",
           (progress ? "Yes" : "No"));

    printf(" 2                 Tree distance submenu:");
    switch (tree_pairing) {
      case NO_PAIRING: 
        printf("\n\nERROR: Unallowable option!\n\n");
        exxit(-1);
        break;

      case ADJACENT_PAIRS:
        printf("  Distance between adjacent pairs\n");
        break;

      case CORR_IN_1_AND_2: 
        printf("  Distances between corresponding \n");
        printf("                                             pairs in first and second tree files\n");
        break;

      case ALL_IN_FIRST: 
        printf("  Distances between all possible\n");
        printf("                                             pairs in tree file.\n");
        break;

      case ALL_IN_1_AND_2: 
        printf("  Distances between all pairs in\n");
        printf("                                              first and second tree files\n");
        break;
    }

    printf("\nAre these settings correct? (type Y or the letter for one to change)\n");
    fflush(stdout);
    scanf("%c%*[^\n]", &ch);
    getchar();
    uppercase(&ch);
    done = (ch == 'Y');
    if (!done) {
      if ((noroot && (ch == 'O')) || strchr("RTD12",ch) != NULL) {
        switch (ch) {

        case 'D':
          if ( dtype == SYMMETRIC )
            dtype = BSD;
          else if ( dtype == BSD )
            dtype = SYMMETRIC;
          break;

        case 'R':
          noroot = !noroot;
          break;

        case 'T':
          initterminal(&ibmpc, &ansi);
          break;

        case '1':
          progress = !progress;
          break;
        
        case '2':
          pairing_submenu();
          break;
        }
      } else
        printf("Not a possible option!\n");
    }
    countup(&loopcount, 100);
  } while (!done);
}  /* getoptions */


int main(int argc, Char *argv[])
{  
  pattern_elm  ***pattern_array;
  char *s; /* for getenv */
  long trees_in_1 = 0, trees_in_2 = 0;
  long tip_count = 0;
  long alt_maxgrp = -1;
  node * p;

#ifdef MAC
  argc = 1;                /* macsetup("Treedist", "");        */
  argv[0] = "Treedist";
#endif
  init(argc, argv);
  /* Open in binary: ftell() is broken for UNIX line-endings under WIN32 */
  openfile(&intree, INTREE, "input tree file", "rb", argv[0], intreename);
  openfile(&outfile, OUTFILE, "output file", "w", argv[0], outfilename);

  /* Initialize option-based variables, then ask for changes regarding
     their values. */
  getoptions();

  ntrees = 0.0;
  lasti  = -1;

  /* read files to determine size of structures we'll be working with */
  trees_in_1 = countsemic(&intree);
  countcomma(&intree,&tip_count);
  tip_count++; /* countcomma does a raw comma count, tips is one greater */


  if ( (tree_pairing == ALL_IN_1_AND_2)
       || (tree_pairing == CORR_IN_1_AND_2) ) {
    /* If another intree file should exist, */
    /* Open in binary: ftell() is broken for UNIX line-endings under WIN32 */
    openfile(&intree2, INTREE2, "input tree file 2", "rb", 
             argv[0], intree2name);
    trees_in_2 = countsemic(&intree2);
  }

  /* 
   * EWFIX.BUG.756
   * 
   * inside cons.c there are several arrays which are allocated
   * to size "maxgrp", the maximum number of groups (sets of
   * tips more closely connected than the rest of the tree) we
   * can see as the code executes.
   *
   * We have two measures we use to determine how much space to
   * allot:
   *  (1) based on the tip count of the trees in the infile
   *  (2) based on total number of trees in infile, and 
   *
   * (1) -- Tip Count Method
   * Since each group is a subset of the set of tips we must
   * represent at most pow(2,tips) different groups. (Technically
   * two fewer since we don't store the empty or complete subsets,
   * but let's keep this simple.
   *
   * (2) -- Total Tree Size Method
   * Each tree we read results in 
   *      singleton groups for each tip, plus
   *      a group for each interior node except the root
   * Since the singleton tips are identical for each tree, this gives
   * a bound of #tips + ( #trees * (# tips - 2 ) )
   *
   *
   * Ignoring small terms where expedient, either of the following should
   * result in an adequate allocation:
   *       pow(2,#tips)
   *       (#trees + 1) * #tips 
   *
   * Since "maxgrp" is a limit on the number of items we'll need to put
   * in a hash, we double it to make space for quick hashing
   */
  maxgrp = pow(2,tip_count); 
  alt_maxgrp = (trees_in_1 + trees_in_2 + 1) * tip_count;
  if(alt_maxgrp < maxgrp) maxgrp = alt_maxgrp;
  maxgrp *= 2;

  /* EWFIX.BUG.756 -- stealing an initial value from consense.
   * This should really be defined in a header, and a warning
   * should be given if the other values are bigger -- presumably
   * that means that we're likely to run out of space.
   *
   * However, since we have the hashing machinery and we don't know
   * how big the user's system is, this should be fine.
   */
  if(maxgrp > 32767) maxgrp = 32767;

  /* Read the (first) tree file and put together grouping, order, and
   * timesseen */
  read_groups (&pattern_array, trees_in_1 + trees_in_2, tip_count, intree);

  if ( (tree_pairing == ADJACENT_PAIRS)
       || (tree_pairing == ALL_IN_FIRST) ) {

    /* Here deal with the adjacent or all-in-first pairing
       difference computation */
    compute_distances (pattern_array, trees_in_1, 0);

  }
  else if ( (tree_pairing == CORR_IN_1_AND_2)
            || (tree_pairing == ALL_IN_1_AND_2) ) {
    /* Here, open the other tree file, parse it, and then put
         together the difference array */
    read_second_file(pattern_array, trees_in_1, trees_in_2);

    compute_distances (pattern_array, trees_in_1, trees_in_2);

  } else if (tree_pairing == NO_PAIRING) {
    /* Compute the consensus tree. */
    putc('\n', outfile);
    /* consensus();         Reserved for future development */
  }

  if (progress)
    printf("\nOutput written to file \"%s\"\n\n", outfilename);

  FClose(outtree);
  FClose(intree);
  FClose(outfile);

  if ((tree_pairing == ALL_IN_1_AND_2) || 
      (tree_pairing == CORR_IN_1_AND_2))
    FClose(intree2);



#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(outtreename);
#endif


  free_patterns (pattern_array, trees_in_1 + trees_in_2);
  clean_up_final();
  /* clean up grbg */
  p = grbg;
  while (p != NULL) {
     node * r = p;
     p = p->next;
     free(r->nodeset);
     free(r->view);
     free(r);
  }


  printf("Done.\n\n");
#ifdef WIN32
  phyRestoreConsoleAttributes();
#endif

  return 0;
}  /* main */

