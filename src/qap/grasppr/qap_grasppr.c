/*
 * qap_grasppr.c
 *
 * $Revision$
 *
 * $Date$
 *
 * Copyright 2005 Sergio Anibal de Carvalho Junior
 *
 * This file is part of ArrayOpt.
 *
 * --- License ----------------------------------------------------------------
 * ArrayOpt is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option)
 * any later version.
 *
 * ArrayOpt is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * ArrayOpt; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307, USA.
 * ----------------------------------------------------------------------------
 *
 * This is the result of a PhD work developed at the Universitaet Bielefeld
 * under the supervision of Dr. Sven Rahmann. Proper attribution of the author
 * as the source of the software is appreciated.
 *
 * Sergio Anibal de Carvalho Jr.  http://www.cebitec.uni-bielefeld.de/~scarvalh
 * AG Genominformatik             http://gi.cebitec.uni-bielefeld.de
 * Universitaet Bielefeld         http://www.uni-bielefeld.de
 *
 */

/*
 * This program provides an interface from ArrayOpt's API to an implementation
 * of the following QAP (Quadratic Assignment Problem) solver:
 *
 * Method:       GRASP for dense QAP
 * Language:     C
 * Authors:      Carlos A. S. Oliveira (AT&T Research, New Jersey)
 *               Panos M. Pardalos (ISE dept., University of Florida)
 *               Mauricio G. C. Resende (AT&T Research, New Jersey)
 * Available at: http://www.research.att.com/~mgcr/exp/gqapspr/
 *
 */

#include <arrayopt_qap_GraspPathRelinking.h>

// GRASP libraries
#define NDEBUG 1
#include "qapinst.h"
#include "qapsol.h"
#include "qapgrasp.h"
#include "qappr.h"
#include "randgen.h"
#include "timer.h"

#define ERROR_CODE		-1
#define MIN_DIMENSION	2

// ****************************************************************************
// Function prototypes
int run_grasp_pr (grasp*, int);

// ****************************************************************************
// Java to C interface

JNIEXPORT jlong JNICALL Java_arrayopt_qap_GraspPathRelinking_qap_1grasppr
	(JNIEnv *env, jclass obj, jint n, jintArray flow, jintArray dist,
	jintArray sol, jfloat alpha, jfloat beta, jint runs, jint max_itr,
	jint look4, jint elite_size, jint max_time, jintArray in_out)
{
	qap_inst	*qap_data;
	grasp		*grasp_param;
	int			i, n2, seed, best_cost;

	// check min dimension
	if (n < MIN_DIMENSION) return ERROR_CODE;

	// length of linear arrays
	// used to store matrices
	n2 = n * n;

	// get reference to Java arrays (input)
	jint *f  = (*env)->GetIntArrayElements(env, flow, 0);
	jint *d  = (*env)->GetIntArrayElements(env, dist, 0);
	jint *s	 = (*env)->GetIntArrayElements(env, sol, 0);
	jint *io = (*env)->GetIntArrayElements(env, in_out, 0);

	// workaround: use an integer array to pass in/out arguments
	// (passing a pointer to an integer doesn't work)

	// TODO add new parameter for seed instead of using this in/out array
	// since this method will not return an updated seed
	seed = io[0];

	// create QAP instance
	qap_data = qi_new (n, f, d);
	if (!qap_data) return -1;

	// create and set GRASP parameters
	grasp_param = g_new(qap_data);
	if (!grasp_param) return -1;
	grasp_param->alpha = alpha;
	grasp_param->beta = beta;
	grasp_param->elite_size = elite_size;
	grasp_param->max_time = max_time;
	grasp_param->l4 = look4;
	grasp_param->curr_prog_iter = 1;

	// no printing of any kind
	grasp_param->print = 0;
	grasp_param->print2 = 0;

	// no random alpha (parameter not used
	// here: alway use same alpha
	grasp_param->ralpha = 0;

	// use path-relink (parameter not used
	// here: always run with path-relink)
	grasp_param->no_pr = 0;

	// set random generator's seed
	sgenrand (seed);

	// run GRASP with path-relinking
	best_cost = run_grasp_pr (grasp_param, max_itr);

	// copy best permutation
	for (i = 0; i < n; i++)
		s[i] = grasp_param->best->p[i];

	// workaround: copy in/out parameters back to the array
	// (passing a pointer doesn't work)
	io[1] = grasp_param->curr_iter;

	g_delete(grasp_param);
	qi_delete(qap_data);

	// release Java arrays' memory
	(*env)->ReleaseIntArrayElements(env, dist, d, 0);
	(*env)->ReleaseIntArrayElements(env, flow, f, 0);
	(*env)->ReleaseIntArrayElements(env, sol, s, 0);
	(*env)->ReleaseIntArrayElements(env, in_out, io, 0);

	// return best cost
	return best_cost;
}

int run_grasp_pr (grasp *param, int max_iter)
{
	int		i;
	elite	*elite_set;

	// create elite set
	elite_set = pr_new (param, param->elite_size);
	if (!elite_set) return -1;

	set_initial_time();

	param->curr_iter = 0;

	// instead of generating a random permutation
	// we start with the natural permutation (0, 1, ...n-1);
	// the implications of this decision for arrayopt depend
	// on how the QAP solver is being used:
	// - in case a new layout is being generated, the natural
	//   permutation is likely to be random anyway
	// - if an existing layout is being optimized,
	//   the current layout is a good starting point, and it
	//   guarantees that the layout will never be worsened

	// compute the cost of such permutation
	qs_objective (param->best, param->q);

	for (i = 0; i < max_iter; i++)
	{
		param->curr_iter++;

		// stop if current cost is equal or better then what we are looking for
		if (param->l4 >= 0 && param->l4 < param->best->cost)
			break;

		// stop if exceeded maximum time
		if (param->max_time > 0 && get_time() >= param->max_time)
			break;

		// main GRASP routines
		constructor (param);
		extra_ls (param);
		pr_run (elite_set, param);
		pr_rev_run (elite_set, param);
		pr_update (elite_set, param);
		update (param);
	}

	pr_post_optimization (elite_set, param);
	update (param);
	pr_delete (elite_set);

	return param->best->cost;
}
