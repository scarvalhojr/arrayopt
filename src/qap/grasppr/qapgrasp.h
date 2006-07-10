#ifndef _QAPGRASP_H
#define _QAPGRASP_H

#include "qapsol.h"

/* --------------------------------------------------------------
 * GRASP structure and functions
 * -------------------------------------------------------------- */

/* pairs are used in the constructor to sort assignments */
typedef struct {
   int cost;
   int i, j;
} pair;

/*
 * This is the GRASP structure, holding data needed by
 * GRASP.
 * */
typedef struct {
   int ralpha;
   float alpha, beta;
   qap_sol *best;
   qap_sol *s;
   qap_inst *q;
   pair *f;
   pair *d;
   int n_assigned; /* no. of assigned pairs */
   int *assigned;  /* list of facilities assigned */
   int *done;      /* done[i] == 1 iff fac. i assigned */
   int *ldone;     /* done[i] == 1 iff loc. i assigned */
   int max_time;   /* maximum time */
   int elite_size;
   int l4;         /* value to look for */
   int print; /* true if we print each improvement */
   int print2; /* true if we print each n^i iterations */
   int curr_iter; /* current iteration */
   int curr_prog_iter; /* current iteration of the program */
   int no_pr;  /* if true, pr won't run */
	int last_improv_iter; /* records the last improving interations */
} grasp;

grasp *g_new(qap_inst *q);
void g_delete(grasp *g);
void phase2(grasp *g);
void constructor(grasp *g);
int g_improvement(qap_sol *s, qap_inst *q, int i, int j);
void local_search(grasp *g);
void update(grasp *g);
void ls_step(grasp *g);
void extra_ls(grasp *g);

#endif
