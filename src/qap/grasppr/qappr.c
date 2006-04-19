#include "qappr.h"
#include "common.h"
#include "timer.h"
#include "randgen.h"
#include <stdlib.h>

// ******************************************************************
// Commented out: no use of following libraries here
//
// #include <stdio.h>
// #include <limits.h>
//
// Code writing to the standard output was deliberately commented out
// since, when integrated into arrayopt, no output will be seen.
//
// -Sergio A. de Carvalho Jr.
// ******************************************************************

/*
 * maximum number of iterations without improvement to remove
 * some elements of elite set
 */
#define MAX_ITER_NO_IMPROV 20
/*
 *  minimum difference required from solutions on the elite set
 *  needed for a solution to be included.
 */
#define MIN_DIFF	3


/* creates a new elite set.
 *   returns zero if an error ocurrs. */
elite *pr_new(grasp *g, int size)
{
   int i;
   elite *e = New(elite);
   if (!e) return 0;
   e->size = size;
   e->sol = NewVec(qap_sol*, size);
   if (!e->sol) return 0;
	e->diff = NewVec(int, size);
	if (!e->diff) goto error1;
   for (i=0; i<size; i++) {
      e->sol[i] = qs_new(g->q->n);
		if (!e->sol[i]) goto error2;
   }
   e->cur_size = 0;
   e->worst = e->best = 0;
	e->g = g;
   return e;
error2:
	for (i--; i>=0; i--)
		qs_delete(e->sol[i]);
	Delete(e->diff);
error1:
	Delete(e->sol);
	return 0;
}

/* delete the elite data structure and all
 * associated data. */
void pr_delete(elite *e)
{
   int i;
   for (i=0; i<e->size; i++)
      qs_delete(e->sol[i]);
	Delete(e->diff);
   Delete(e);
}

/* copies elite set e2 to e1 */
int pr_copy(elite *e, elite *e2)
{
	int i;
	for (i=0; i<e2->cur_size; i++) {
		qs_copy(e->sol[i], e2->sol[i]);
	}
   e->size = e2->size;
	e->cur_size = e2->cur_size;
   e->worst = e2->worst;
	e->best = e2->best;
	return 1;
}

/* the local search algorithm. */
void pr_local_search(grasp * g)
{
	int n = g->q->n;
	/* NB.: values of p,q just for the compiler */
	int i, j, p =-1, q =-1;
	int nochange = 0, cost, best;
	while (nochange < 20) {
		if (g->max_time > 0 && get_time() >= g->max_time)
			break;
		i = genrandint() % n;
		best = 0;
		for (j = 0; j < n; j++) {
			if (j != i) {
				cost = g_improvement(g->s, g->q, i, j);
				if (cost > best) {
					best = cost;
					p = i;
					q = j;
				}
			}
		}
		if (best > 0) {
			qs_swap(g->s, p, q);
			g->s->cost -= best;
			nochange = 0;
		} else
			nochange++;
	}
}

/* do an extra local search by using varying neigborhoods.
 * To implement this we need just to flip assinments two times. */
void extra_ls(grasp *g)
{
   int i;
   for (i=0; i<10; i++) {
      pr_local_search(g);
      update(g);
      /* do 2 local search steps and continue */
      ls_step(g);
      ls_step(g);
   }
}

/*
 *  execute the path relinking.
 *  s1 is the starting solution
 *  s2 is the ending solution
 *  The resulting solution is given by the following rule:
 *  	if it has a solution better than the current or the guiding,
 *		then this is the solution
 *		else
 *			if there is solution s which is better then the previous and
 *				worse than the next
 *			then return s
 *			else return the best of the current and guiding solutions
 **/
int execute_pr(qap_sol * s1, qap_sol * s2, grasp * g)
{
	qap_sol *sbest, *es, *s, *prev, *lopt;
	/* cost of the previous and prev-prev solution */
	int cpp = -1, cp = -1;
	int improv, i, j, n = g->s->n;
	/* create the starting solutions  */
	prev = qs_new(s1->n);
	if (!prev) return 0;
	lopt = qs_new(s1->n);
	if (!lopt) goto error1;
	lopt->cost = -1;  /* mark this solution as not used */
	s = qs_new(s1->n);
	if (!s) goto error2;
	qs_copy(s, s1);
	/* s2 is the guiding solution */
	es = s2;
	/* save the best solution */
	sbest = qs_new(s1->n);
	if (!sbest) goto error3;
	qs_copy(sbest, s1);
	for (i = 0; i < n; i++)
		if (s->p[i] != es->p[i]) {
			/* copy previous values */
			cpp = cp;
			cp = s->cost;
			qs_copy(prev, s);
			/* find the facility that must be changed */
			j = s->rev[es->p[i]];
			improv = g_improvement(s, g->q, i, j);
			qs_swap(s, i, j);
			s->cost -= improv;
			if (s->cost < sbest->cost) {
				qs_copy(g->s, s);
				pr_local_search(g);
				qs_copy(sbest, g->s);
			}
			if (cpp!= -1 && cp < cpp && cp < s->cost
				&& (lopt->cost == -1 || cp < lopt->cost)) {
				/* the previous solution is a local optimum */
				qs_copy(lopt, prev);
				/* printf("found a local sol %d, %d, %d\n", cpp, cp, s->cost); */
			}
		}
	if (sbest->cost < s2->cost) {
		qs_copy(g->s, sbest);
	} else if (lopt->cost > -1) {
		qs_copy(g->s, lopt);
	} else {
		/* no option other then return this */
		qs_copy(g->s, sbest);
	}
	qs_delete(sbest);
	qs_delete(s);
	qs_delete(lopt);
	qs_delete(prev);
	return 1;
error3:
	qs_delete(s);
error2:
	qs_delete(lopt);
error1:
	qs_delete(prev);
	return 0;
}

/*
 * returns the number of a  solution chosen from the
 * elite pool. This propability of a solution $s$ be chosen is
 * proportional to the difference between $s$ and
 * the current solution.
 */
int pr_get_guiding_sol(elite *e, qap_sol *sol)
{
	int i, r, total=0;
	for (i=0; i<e->size; i++) {
		e->diff[i] = qs_similarity(sol, e->sol[i]);
		total += e->diff[i];
	}
	for (i=1; i<e->size; i++)
		e->diff[i] += e->diff[i-1];
   r = genrandint()%total;
	for (i=0; i<e->size && r > e->diff[i]; i++)
		;
	return i;
}

/*  do the real job of reverse path relinking */
void pr_rev_run(elite *e, grasp *g)
{
	int guide;
   /* check if we have enough solutions */
   if (e->cur_size == e->size) {
		guide = pr_get_guiding_sol(e, g->s);
      execute_pr(e->sol[guide], g->s, g);
	}
}

/* inset a new solution into the elite set */
void pr_insert_sol(elite * e, grasp * g)
{
	/* just add the current solution here */
	qs_copy(e->sol[e->cur_size], g->s);
	if (e->cur_size > 0) {
		if (e->sol[e->cur_size]->cost > e->sol[e->worst]->cost)
			e->worst = e->cur_size;
		if (e->sol[e->cur_size]->cost < e->sol[e->best]->cost)
			e->best = e->cur_size;
	}
	e->cur_size++;
}


/* do the path relinking */
void pr_run(elite *e, grasp *g)
{
	int guide, i, reject=0;
   /* check if there are enough solutions in the elite set */
   if (e->cur_size < e->size) {
      /* add the current solution here if it is different of the
			existing solutions */
		for (i=0; i<e->cur_size; i++) {
			if (qs_similarity(g->s, e->sol[i]) < MIN_DIFF)
				reject = 1;
		}
		if (!reject)
      	pr_insert_sol(e, g);
      return;
   } else {
		/* find a guiding solution */
		guide = pr_get_guiding_sol(e, g->s);
      execute_pr(g->s, e->sol[guide], g);
   }
}

/*
 * test if solution s is in the pool of elite solutions
 * (in fact, this will return 1 if the solution is equal of
 * very close (< MIN_DIFF) positions )
 */
int pr_sol_in_e(elite * e, qap_sol * s)
{
	int i;
	for (i = 0; i < e->size; i++)
		if (e->sol[i]->cost == s->cost
			|| qs_similarity(e->sol[i], s) < MIN_DIFF)
			return 1;
	return 0;
}

/* used by quicksort */
int pr_compar(const void *a, const void *b)
{
	const qap_sol *sa = *((qap_sol**)a), *sb = *((qap_sol**)b);
	if (sa->cost > sb->cost) return -1;
	if (sa->cost < sb->cost) return 1;
	return 0;
}

/*
 * mark half of the elit set with infinite cost
 */
static int mark_half_set_infinite(elite * e, grasp * g)
{
	qap_sol **scp, *sol;
	int i, j, best, worst;
	scp = NewVec(qap_sol *, e->size);
	if (!scp) goto error0;
	for (i = 0; i < e->size; i++) {
		scp[i] = qs_new(g->q->n);
		if (!scp[i]) goto error1;
		qs_copy(scp[i], e->sol[i]);
		/* use cur_size to mark the original position */
		scp[i]->n = i;
	}
	qsort(scp, e->size, sizeof(qap_sol *), pr_compar);
	for (i = 0; i < e->size / 2; i++) {
		sol = e->sol[scp[i]->n];
		/* mark this position to deletion */
		e->sol[scp[i]->n]->cost = -1;
		e->cur_size--;
	}
	/* copy remaining solutions to empty slots, to avoid holes */
	for (i=0, j=e->size-1; i < j; i++) {
		if (e->sol[i]->cost == -1) {
			while (e->sol[j]->cost == -1 && i<j)
				j--;
			if (i>=j) continue;
			/* exchange these two positions */
			sol = e->sol[i];
			e->sol[i] = e->sol[j];
			e->sol[j] = sol;
			j--;
		}
	}
	/* update best and worst */
	best = 999999;
	worst = 0;
	for (i = 0; i < e->cur_size; i++) {
		if (e->sol[i]->cost < best) {
			e->best = i;
			best = e->sol[i]->cost;
		}
		if (e->sol[i]->cost > worst) {
			e->worst = i;
			worst = e->sol[i]->cost;
		}
	}
	/* reset our counter */
	g->last_improv_iter = g->curr_iter;
	/* delete allocated memory */
	for (i = 0; i < e->size; i++)
		qs_delete(scp[i]);
	Delete(scp);
	return 1;
error1:
	for (i--; i >= 0; i--)
		qs_delete(scp[i]);
	Delete(scp);
error0:

	// *****************************************************
	// Commented out: nothing written to stdout will be seen
	//
	// printf("error allocating mem in pr_update\n");
	//
	// -Sergio A. de Carvalho Jr.
	// *****************************************************

	return 0;
}

/* update the elite set with the current solution */
void pr_update(elite * e, grasp * g)
{
	int i, diff, sim, position = 0;
	if (e->cur_size != e->size)
		return;
	if (g->s->cost < e->sol[e->best]->cost
	    || (g->s->cost < e->sol[e->worst]->cost &&
		!pr_sol_in_e(e, g->s))) {
		/*
		 * find a place to insert the solution: the solution most
		 * similar to the current one.
		 */
		diff = g->s->n + 1;  /* assume the biggest  difference */
		for (i = 0; i < e->size; i++) {
			/*
			 * check only for solutions with cost greater than or
			 * equal to the current solution
			 */
			if (e->sol[i]->cost >= g->s->cost) {
				sim = qs_similarity(e->sol[i], g->s);
				if (diff > sim) {
					diff = sim;
					position = i;
				}
				/*
				 * if there is no difference, I just flip a
				 * coin to decide
				 */
				if (diff == sim && (genrandint() % 2 == 0)) {
					position = i;
				}
			}
		}
		/* insert the new solution */
		qs_copy(e->sol[position], g->s);
		/* update field best */
		if (g->s->cost < e->sol[e->best]->cost)
			e->best = position;
		/* update field worst */
		if (g->s->cost > e->sol[e->worst]->cost)
			e->worst = position;
	} else if (g->curr_iter - g->last_improv_iter >=  MAX_ITER_NO_IMPROV) {
		mark_half_set_infinite(e, g);
	}
}

/*
 * used by post optimization to run the path relinking
 */
static void call_pr(qap_sol * s1, qap_sol * s2, elite * e, grasp *g)
{
	execute_pr(s1, s2, g);
	if (e->cur_size < e->size) {
		pr_insert_sol(e, g);
		return;
	} else {
		pr_update(e, g);
	}
	update(g);
}

/* this is executed at the end of the program, to garantee
 * local optimality among the elements of the elite set.
 * */
int pr_post_optimization(elite * e, grasp * g)
{
	int i, j, change, cost;
	elite *ecp;
	ecp = pr_new(g, e->size);
	if (!ecp)
		return 0;
	if (e->cur_size < e->size)
		return 0;
	change = 1;
	while (change) {
		change = 0;
		pr_copy(ecp, e);
		/* this makes the elite set empty */
		e->cur_size = 0;
		cost = g->best->cost;
		for (i = 0; i < e->size; i++) {
			for (j = 0; j < e->size; j++) {
				if (i == j) continue;
				call_pr(ecp->sol[i], ecp->sol[j], e, g);
				call_pr(ecp->sol[j], ecp->sol[i], e, g);
			}
		}
		if (g->best->cost < cost)
			{
				change = 1;

				// *****************************************************
				// Commented out: nothing written to stdout will be seen
				//
				// printf("Post-opt improv: %d\n",g->best->cost);
				//
				// -Sergio A. de Carvalho Jr.
				// *****************************************************
			}
	}
	pr_delete(ecp);
	return 1;
}

/* returns a number representing the diversity of the
 * solutions in the elite set. Is the sum of the differences
 * between each pair of solutions
 * */
int pr_diversity(elite * e)
{
	int i, j, d = 0;
	for (i = 0; i < e->cur_size; i++) {
		for (j = 0; j < e->cur_size; j++) {
			d += qs_similarity(e->sol[i], e->sol[j]);
		}
	}
	return d;
}
