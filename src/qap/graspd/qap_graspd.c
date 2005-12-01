/*
 * qap_graspd.c
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
 * Method: GRASP for dense QAP
 * Implementation: Fortran
 * Provided by: Mauricio G. C. Resende
 * Available at: http://www.research.att.com/~mgcr/
 *
 */

#include <cfortran.h>
#include <arrayopt_qap_GraspDense.h>

#define ERROR_CODE		-1
#define MIN_DIMENSION	2

// ****************************************************************************
// C to Fortran prototype

	PROTOCCALLSFSUB24 (GQAPD,gqapd, INT,INT,INT,FLOAT,FLOAT,INT,PINT,INTV,	\
							INTV,INTV,INTV,INTV,INTV,INTV,INTV,INTV,INTV,	\
							INTV,INTV,INTV,INTV,INTV,PINT,PINT)

#define QAP_GRASP_DENSE(N,N2,NITER,ALPHA,BETA,LOOK4,SEED,F,					\
							D,A,B,SRTF,SRTIF,SRTD,SRTID,SRTC,SRTIC,			\
							INDEXD,INDEXF,COST,FDIND,OPTA,BESTV,ITER)		\
		CCALLSFSUB24 (GQAPD,gqapd, INT,INT,INT,FLOAT,FLOAT,INT,PINT,INTV,	\
							INTV,INTV,INTV,INTV,INTV,INTV,INTV,INTV,INTV,	\
							INTV,INTV,INTV,INTV,INTV,PINT,PINT,				\
							N,N2,NITER,ALPHA,BETA,LOOK4,SEED,F,				\
							D,A,B,SRTF,SRTIF,SRTD,SRTID,SRTC,SRTIC,			\
							INDEXD,INDEXF,COST,FDIND,OPTA,BESTV,ITER)

// ****************************************************************************
// Java to C interface

JNIEXPORT jlong JNICALL Java_arrayopt_qap_GraspDense_qap_1graspd
	(JNIEnv *env, jclass obj, jint n, jint niter, jfloat alpha, jfloat beta,
		jint look4, jintArray dist, jintArray flow, jintArray sol,
		jintArray in_out)
{
	int	n2, seed, bestv, iter;
	int	*a, *b, *srtf, *srtif, *srtd, *srtid, *srtc, *srtic, *idxd, *idxf,
		*cost, *fdind;

	// check min dimension
	if (n < MIN_DIMENSION) return ERROR_CODE;

	// length of linear arrays
	// used to store matrices
	n2 = n * n;

	// allocate temporary working space
	if (!(a 	= malloc (n  * sizeof(int)))) return ERROR_CODE;
	if (!(b 	= malloc (n  * sizeof(int)))) return ERROR_CODE;
	if (!(srtf	= malloc (n2 * sizeof(int)))) return ERROR_CODE;
	if (!(srtif	= malloc (n2 * sizeof(int)))) return ERROR_CODE;
	if (!(srtd	= malloc (n2 * sizeof(int)))) return ERROR_CODE;
	if (!(srtid	= malloc (n2 * sizeof(int)))) return ERROR_CODE;
	if (!(srtc	= malloc (n2 * sizeof(int)))) return ERROR_CODE;
	if (!(srtic	= malloc (n2 * sizeof(int)))) return ERROR_CODE;
	if (!(idxd	= malloc (n2 * sizeof(int)))) return ERROR_CODE;
	if (!(idxf	= malloc (n2 * sizeof(int)))) return ERROR_CODE;
	if (!(cost	= malloc (n2 * sizeof(int)))) return ERROR_CODE;
	if (!(fdind	= malloc (n2 * sizeof(int)))) return ERROR_CODE;

	// get reference to Java arrays
	jint *d    = (*env)->GetIntArrayElements(env, dist, 0);
	jint *f    = (*env)->GetIntArrayElements(env, flow, 0);
	jint *opta = (*env)->GetIntArrayElements(env, sol, 0);
	jint *io   = (*env)->GetIntArrayElements(env, in_out, 0);

	// workaround (passing a pointer doesn't work)
	seed = io[0];

	// call fortran subroutine
	QAP_GRASP_DENSE (n,n2,niter,alpha,beta,look4,seed,f,d,a,b,
						srtf,srtif,srtd,srtid,srtc,srtic,idxd,
						idxf,cost,fdind,opta,bestv,iter);

	// workaround (passing a pointer doesn't work)
	io[0] = seed;
	io[1] = iter;

	// release Java arrays memory
	(*env)->ReleaseIntArrayElements(env, dist, d, 0);
	(*env)->ReleaseIntArrayElements(env, flow, f, 0);
	(*env)->ReleaseIntArrayElements(env, sol, opta, 0);
	(*env)->ReleaseIntArrayElements(env, in_out, io, 0);

	// release allocated memory
	free(a);
	free(b);
	free(srtf);
	free(srtif);
	free(srtd);
	free(srtid);
	free(srtc);
	free(srtic);
	free(idxd);
	free(idxf);
	free(cost);
	free(fdind);

	return bestv;
}
