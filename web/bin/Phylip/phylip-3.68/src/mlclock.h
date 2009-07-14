#include "phylip.h"

/* Uncomment this line to dump details to dnamlk.log */
/* #define DEBUG */

#ifdef DEBUG
extern double    dump_likelihood_graph(FILE *fp, node *p, double min, double max, double step);
extern double dump_likelihood_graph_twonode(FILE *fp, node *p1, double npoints);
double dump_likelihood_graph_2d(FILE *fp, node *p1, double npoints);
extern void    get_limits(node **nodea, double *min, double *max);
extern double  set_tyme_evaluate(node *, double);
#endif /* DEBUG */

typedef double (*evaluator_t)(node *);

extern const double MIN_BRANCH_LENGTH;
extern const double MIN_ROOT_TYME;

/* module initialization */
extern void mlclock_init(tree *t, evaluator_t f);

/* check or fix node tymes and branch lengths */
extern boolean all_tymes_valid(node *, double, boolean);

/* change node tymes */
extern void setnodetymes(node* p, double newtyme);

/* limits of node movement */
extern double min_child_tyme(node *p);
extern double parent_tyme(node *p);
extern boolean valid_tyme(node *p, double tyme);

/* save/restore tymes */
extern void save_tymes(tree* save_tree, double tymes[]);
extern void restore_tymes(tree *load_tree, double tymes[]);

/* optimize a node tyme */
double maximize(double min_tyme, double cur, double max_tyme, double(*f)(double), double eps, boolean *success);
extern boolean makenewv(node *p);
