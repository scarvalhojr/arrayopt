#include "qapsol.h"
#include "common.h"
#include "randgen.h"

// ******************************************************************
// Commented out: no use of following libraries here
//
// #include <stdio.h>
// #include <stdlib.h>
//
// Code writing to the standard output was deliberately commented out
// since, when integrated into arrayopt, no output will be seen.
//
// -Sergio A. de Carvalho Jr.
// ******************************************************************

qap_sol *qs_new(int n)
{
   int i;
   qap_sol *s = New(qap_sol);
   if (n < 1) {

		// *****************************************************************
		// Commented out: nothing written to stdout will be seen
		//
		// fprintf(stderr, "error: wrong argumet for qs_new: %d\n", n);
		//
		// -Sergio A. de Carvalho Jr.
		// *****************************************************************

      return 0;
   }
   if (!s) return 0;
   s->p = NewVec(int, n);
   if (!s->p) return 0;
   s->rev = NewVec(int, n);
   if (!s->rev) return 0;
   /* init the permutation */
   for (i=0; i<n; i++) {
      s->p[i] = s->rev[i] = i;
   }
   s->n = n;
   return s;
}

void qs_delete(qap_sol *s)
{
   if (s->p)
      Delete(s->p);

	// *****************************************************************
	// Added: should also delete the rev array
	//
	if (s->rev)
		Delete (s->rev);
	//
	// -Sergio A. de Carvalho Jr.
	// *****************************************************************

   Delete(s);
}

/* swap the two positions a and b in the permutation */
void qs_swap(qap_sol *s, int a, int b)
{
   int t;
   s->rev[s->p[a]] = b;
   s->rev[s->p[b]] = a;
   /* exchange */
   t = s->p[a];
   s->p[a] = s->p[b];
   s->p[b] = t;
}

/* make an assignment of val to position pos */
void qs_assign(qap_sol *s, int pos, int val)
{
   qs_swap(s, pos, s->rev[val]);
}

/* create a random solution */
void qs_random(qap_sol *s)
{
   int i, r, t;
   for (i=0; i<s->n; i++)
      s->p[i] = i;
   for (i=0; i<s->n; i++) {
      r = genrandint()%(s->n-i);
      /* exchange it */
      t = s->p[i];
      s->p[i] = s->p[r];
      s->p[r] = t;
   }

	// *****************************************************************
	// Added: need to update the reverse array
	//
	for (i = 0; i < s->n; i++)
		s->rev[s->p[i]] = i;
	//
	// -Sergio A. de Carvalho Jr.
	// *****************************************************************
}

// ***********************************************
// Commented out: no need of input/output routines
//
// /* show the contents of the solution */
// void qs_print(qap_sol *s)
// {
//    int i;
//    for (i=0; i<s->n; i++) {
//       printf("%d ", s->p[i]+1);
//    }
//    printf("\n");
// }
//
// -Sergio A. de Carvalho Jr.
// ***********************************************

/* copy the contents of solution ns to s */
void qs_copy(qap_sol *s, qap_sol *ns)
{
   int i;
   for (i=0; i<s->n; i++) {
      s->p[i] = ns->p[i];
      s->rev[i] = ns->rev[i];
   }
   s->cost = ns->cost;
   s->n = ns->n;
}

/* 1 if the two solutions are equal
 * 0 otherwise */
int qs_equal(qap_sol *sa, qap_sol *sb)
{
   int i;
   if (sa->cost != sb->cost || sa->n != sb->n)
      return 0;
   for (i=0; i<sa->n; i++)
      if (sa->p[i] != sb->p[i])
	 return 0;
   return 1;
}

/* returns a value representing the similarity of solutions
 * s1 and s2. This is measured as the number of different assignments
 * in the two solutions.
 * If the result is zero, the solutions are equal.
 * Higher values represent more different solutions. */
int qs_similarity(qap_sol * s1, qap_sol * s2)
{
	int i, d = 0;
	for (i = 0; i < s1->n; i++) {
		if (s1->p[i] != s2->p[i]) {
			d++;
		}
	}
	return d;
}

/* compute the objective function for the permutation in
 * solution s */
void qs_objective(qap_sol *s, qap_inst *qi)
{
   int i, j, n = qi->n, cost=0;
   for (i=0; i<n; i++)
      for (j=0; j<n; j++)
	 cost += qi->D[i][j] * qi->F[s->p[i]][s->p[j]];
   s->cost = cost;
}
