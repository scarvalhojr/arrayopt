#include "qapgrasp.h"
#include "common.h"
#include "sort.h"
#include "randgen.h"
#include "timer.h"

// ******************************************************************
// Commented out: no use of following libraries here
//
// #include <stdio.h>
// #include <stdlib.h>
//
// Code writing to the standard output was deliberately commented out
// since, when integrated into arrayopt, no output will be seen.
//
// Added: use of memset function in string library was implicit
//
#include <string.h>
//
// -Sergio A. de Carvalho Jr.
// ******************************************************************

/*
 * creates a new grasp data structure.
 * returns zero on failure,
 *   the pointer to the instance if successful
 */
grasp *g_new(qap_inst *q)
{
   int n = q->n;
   grasp *g = New(grasp);
   if (!g) return 0;
   g->q = q;
   g->s = qs_new(n);
   g->best = qs_new(n);
   if (!g->best || !g->s) return 0;
   g->f = NewVec(pair, n*n-n+1);
   g->d = NewVec(pair, n*n-n+1);
   if (!g->f || !g->d) return 0;
   g->ralpha = 0;    /* this means that alpha is *not* random */
   g->alpha = .10;   /* default values for alpha */
   g->beta = .5;     /* and beta */
   g->assigned = NewVec(int, n);
   if (!g->assigned) return 0;
   g->done = NewVec(int, n);
   if (!g->done) return 0;
   g->ldone = NewVec(int, n);
   if (!g->ldone) return 0;
   g->elite_size = 10;
   g->max_time = 0;
   g->l4 = -1;        /* default: don't look for any value */
   g->print = 0;
   g->print2 = 0;
   g->curr_iter = 0;
   g->no_pr = 0;
	g->last_improv_iter = 0;
   return g;
}

void g_delete(grasp *g)
{
   if (g->s)
      qs_delete(g->s);
   if (g->best)
      qs_delete(g->best);
   if (g->f) Delete(g->f);
   if (g->d) Delete(g->d);
   Delete(g);
}

/*
 * set of macros to define different heap sort functions.
 * */

#define pairinc pair
#define heap_compare(p,q) (a[p].cost > a[q].cost)
#define heap_swap(p,q) \
  do { pairinc x = a[p]; a[p] = a[q]; a[q] = x; } while (0)
define_heap_operation(pairinc)
#undef heap_compare
#undef heap_swap
#undef pairinc


#define pairdec pair
#define heap_compare(p,q) (a[p].cost < a[q].cost)
#define heap_swap(p,q) \
  do { pairdec x = a[p]; a[p] = a[q]; a[q] = x; } while (0)
define_heap_operation(pairdec)
#undef heap_compare
#undef heap_swap
#undef pairdec

#define doublepair pair
#define S ((doublepair*)s)
#define heap_compare(p,q) (a[p].cost*S[p].cost < a[q].cost*S[q].cost)
#define heap_swap(p,q) \
  do { doublepair x = a[p], y = S[p]; a[p] = a[q]; a[q] = x; \
                                      S[p] = S[q]; S[q] = y;  } while (0)
/* I am using another definition bec. we need to sort just
 * some values, not all of them */
define_extended_heap_operation(doublepair)
#undef heap_compare
#undef S
#undef heap_swap
#undef doublepair

/* helper function to sort the vectors in the constructor phase 1 */
void sort_vector(pair *d, int **s, int n, int inc)
{
   int i, j, k;
   for (k=0, i=0; i<n; i++)
      for (j=0; j<n; j++)
	 if (i!=j) {
	    d[k].cost = s[i][j];
	    d[k].i = i;
	    d[k].j = j;
	    k++;
	 }
   /* sort the resulting vector
      OBS.:note that we sort a zero based vector  */
   d[k] = d[0];
   if (inc) heap_sort_pairinc(d, 0, k);
   else     heap_sort_pairdec(d, 0, k);
   for (i=0; i<k; i++) d[i] = d[i+1];
}

/* assign q to positon p */
void g_assign(grasp *g, int p, int q)
{
   qs_assign(g->s, p, q);
   g->assigned[g->n_assigned] = p;
   g->n_assigned++;
   g->done[p] = g->ldone[q] = 1;
}

/* In this phase we find the two initial assignments. */
void phase1(grasp * g)
{
	int pos, n = g->q->n;
	int last = (int) ((n * n - n) * g->beta);
	int rcl_size = (int) (last * g->alpha);
	/* avoid this to be 0, bec. of divisions */
	if (!rcl_size)
		rcl_size++;
	sort_vector(g->d, g->q->D, n, 1);
	sort_vector(g->f, g->q->F, n, 0);
	/* choose the position in RCL */
	pos = genrandint() % rcl_size;
	pos++;  /* we avoid chosing 0 because of the sorting bellow */
	/* now construct a heap and remove the pos-th smallest */
	g->f[n] = g->f[0];
	g->d[n] = g->d[0];
	heap_sortn_doublepair(g->f, g->d, n * n - n, n * n - n);
	/*
	 * NOTE: heapsort work by putting the last elements at the end of the
	 * vector, take care below...
	 */
	pos = n * n - n - pos + 1;
	/* now we have the initial two assignments */
	memset(g->done, 0, sizeof(int) * n);
	memset(g->ldone, 0, sizeof(int) * n);
	g->n_assigned = 0;
	g_assign(g, g->f[pos].i, g->d[pos].i);
	g_assign(g, g->f[pos].j, g->d[pos].j);
}

/* compute the cost of assign q to position p relative
 * to the previous assignments. */
int g_compute_cost(grasp *g, int p, int q)
{
   int i, cost = 0, k = g->n_assigned;
   int a, b;
   for (i=0; i<k; i++) {
      a = g->assigned[i];
      b = g->s->p[a];
      cost += g->q->F[a][p] * g->q->D[b][q];
   }
   return cost;
}

/* this is used only in sparse instances: it can speedup the
 * process by computing assignment with zero cost relative to
 * the previous assignments. */
void sparse_phase2(grasp * g)
{
	pair *delta, *ndelta, *tmp, last;
	int i, j, k, l, r, nd, n = g->q->n;
	for (l = 0, i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			if (!g->done[i] && !g->ldone[j] &&	/* not assigned yet */
			    g_compute_cost(g, i, j) == 0) {
				/* add (i,j) to delta */
				g->d[l].i = i;
				g->d[l].j = j;
				l++;
			}
	if (l == 0) return;
	r = genrandint() % l;
	last = g->d[r];
	g_assign(g, last.i, last.j);
	delta = g->d;
	ndelta = g->f;
	for (k = 0; k < n - 3; k++) {
		nd = l;
		for (l = 0, i = 0; i < nd; i++)
			/* we insert the elements in the previous delta
			 * which are not assigned yet */
			if (!g->done[delta[i].i] && !g->ldone[delta[i].j] &&
			/* and give zero cost relative to last element
			 * included */
			    (g->q->F[delta[i].i][last.i] *
			     g->q->D[delta[i].j][last.j]  ) == 0) {
				ndelta[l].i = delta[i].i;
				ndelta[l].j = delta[i].j;
				l++;
			}
		if (l == 0)
			return;
		r = genrandint() % l;
		last = ndelta[r];
		g_assign(g, last.i, last.j);
		tmp = delta;
		delta = ndelta;
		ndelta = tmp;
	}
}

/* in this phase we do the remaining assignments */
void phase2(grasp *g)
{
   int pos, i, j, k, l, n = g->q->n, p, max;
   i = g->assigned[0];
   j = g->s->p[i];
   k = g->assigned[1];
   l = g->s->p[k];
   if (g->q->F[i][k] * g->q->D[j][l] == 0) { /* there is some sparsity */
      sparse_phase2(g);
   }
   max = n-g->n_assigned-1;
   for (k=0; k<max; k++) {
      for (l=0, i=0; i<n; i++)
	 for (j=0; j<n; j++)
	    if (!g->done[i] && !g->ldone[j]) {
	       g->f[l].cost = g_compute_cost(g, i, j);
	       g->f[l].i = i;
	       g->f[l].j = j;
	       l++;
	    }
      /* sort the costs */
      g->f[l] = g->f[0];
      /* FIXME: can be improved by sorting only first pos elements */
      heap_sort_pairinc(g->f, 0, l);
      p = (int)(l * g->alpha);
      if (!p) p++;    /* this is to avoid zero the next operation. */
      pos = genrandint()% p;
      pos++;  /* NOTE: we avoid selecting 0 here, bec. of sorting! */
      g_assign(g, g->f[pos].i, g->f[pos].j);
   }
}

/* this is the greed randomized adptive constructor. */
void constructor(grasp *g)
{
   phase1(g);
   phase2(g); /* O(n^3) */
   qs_objective(g->s, g->q);
}

/* compute the improvement of interchanging the assignments in
 * positions i and j. */
int g_improvement(qap_sol *s, qap_inst *q, int i, int j)
{
   /* a is permutation */
   int *a = s->p;
   int **d = q->D, **f = q->F;
   int k, n = q->n, xgain = 0;
   for (k = 0; k<n; k++) {
      if (k!=i  && k!=j) {
	 xgain +=
	    (d[k][i]-d[k][j])*(f[a[k]][a[i]]-f[a[k]][a[j]]) +
	    (d[i][k]-d[j][k])*(f[a[i]][a[k]]-f[a[j]][a[k]]);
      }
   }
   xgain += (d[i][j]-d[j][i])*(f[a[i]][a[j]]-f[a[j]][a[i]]);
   return xgain;
}

void local_search(grasp * g)
{
	int change = 1, i, j, cost, n = g->q->n;
	while (change) {
		change = 0;
		/* try all possible changes */
		for (i = 0; i < n - 1; i++) {
			if (g->max_time > 0 && get_time() >= g->max_time)
				return;
			for (j = i + 1; j < n; j++) {
				cost = g_improvement(g->s, g->q, i, j);
				if (cost > 0) {
					qs_swap(g->s, i, j);
					g->s->cost -= cost;
					change = 1;
				}
			}
		}
	}
}

/* update the best solution found by GRASP */
void update(grasp * g)
{
	if (g->s->cost < g->best->cost) {
		 qs_copy(g->best, g->s);
		/* update the improvement counter */
		g->last_improv_iter = g->curr_iter;

		// *****************************************************************
		// Commented out: nothing written to stdout will be seen
		//
		// if (g->print) {
		// 	printf("\nn: %d iter: %d time: %.2f best: %d\n",
		// 		g->curr_prog_iter, g->curr_iter, get_time(), g->best->cost);
		// 	qs_print(g->best);
		// 	fflush(stdout);
		// }
		//
		// -Sergio A. de Carvalho Jr.
		// *****************************************************************
	}
}

/* executes a single step of local search */
void ls_step(grasp *g)
{
   int p, q, n = g->q->n;
   p = genrandint()%n;
   do {
		q = genrandint()%n;
	} while (p == q);
   g->s->cost -= g_improvement(g->s, g->q, p, q);
   qs_swap(g->s, p, q);
}
