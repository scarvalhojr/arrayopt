#ifndef _QAPPR_H
#define _QAPPR_H

#include "qapgrasp.h"

/* --------------------------------------------------------------
 * Path Relinking implementation
 *    Used by GRASP in the intensification phase.
 * ------------------------------------------------------------- */

/* Elit set structure:
 *   represent the pool of ``good'' solutions 
 * */
typedef struct {
   int size, cur_size;
   int worst, best;
   qap_sol **sol;
	int  *diff;   /* probabilities given to each element */
	grasp *g;
} elite;
elite *pr_new(grasp *g, int size);
void pr_delete(elite *e);
int pr_copy(elite *e, elite *e2);
void pr_local_search(grasp * g);
void extra_ls(grasp *g);
int execute_pr(qap_sol * s1, qap_sol * s2, grasp * g);
int pr_get_guiding_sol(elite *e, qap_sol *sol);
void pr_rev_run(elite *e, grasp *g);
void pr_insert_sol(elite *e, grasp *g);
void pr_run(elite *e, grasp *g);
void pr_update(elite * e, grasp * g);
int pr_post_optimization(elite * e, grasp * g);
int pr_diversity(elite *e);

#endif
