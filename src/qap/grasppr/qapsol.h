#ifndef _QAPSOL_H
#define _QAPSOL_H

#include "qapinst.h"

/* This represents the solution: a permutation (and reverse
 * link) and its cost */
typedef struct {
   int n;
   int *p, *rev;
   int cost;
} qap_sol;

qap_sol *qs_new(int n);
void qs_delete(qap_sol *s);
void qs_swap(qap_sol *s, int a, int b);
void qs_assign(qap_sol *s, int pos, int val);
void qs_random(qap_sol *s);

// ***********************************************
// Commented out: no need of input/output routines
//
// void qs_print(qap_sol *s);
//
// -Sergio A. de Carvalho Jr.
// ***********************************************

void qs_copy(qap_sol *s, qap_sol *ns);
int qs_equal(qap_sol *sa, qap_sol *sb);
int qs_similarity(qap_sol *s1, qap_sol *s2);
void qs_objective(qap_sol *s, qap_inst *qi);

#endif
