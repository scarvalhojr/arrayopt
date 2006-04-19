/*
 * GraspDense.java
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
 * ÂrrayOpt is distributed in the hope that it will be useful, but WITHOUT ANY
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

package arrayopt.qap;

/**
 *
 */
public class GraspDense extends QAPSolverAlgorithm
{
	static
	{
		String lib_name = "qap_graspd";

		try
		{
			System.loadLibrary(lib_name);
		}
		catch (UnsatisfiedLinkError e)
		{
			System.err.println ("Unable to load external library '" +
				System.mapLibraryName(lib_name) + "'.\nLibrary path is: " +
				System.getProperty("java.library.path"));

			throw e;
		}
	}

	protected int max_iter;

	protected float alpha;

	protected float beta;

	protected int seed;

	protected int last_num_iter;

	protected int in_out[];

	public static final int DEFAULT_MAX_ITERACTIONS = 32;

	public static final float DEFAULT_ALPHA = .1f;

	public static final float DEFAULT_BETA = .4f;

	public static final int DEFAULT_SEED = 270001;

	public GraspDense ()
	{
		this (DEFAULT_MAX_ITERACTIONS, DEFAULT_ALPHA, DEFAULT_BETA, DEFAULT_SEED);
	}

	public GraspDense (int max_iter)
	{
		this (max_iter, DEFAULT_ALPHA, DEFAULT_BETA, DEFAULT_SEED);
	}

	public GraspDense (int max_iter, float alpha, float beta, int seed)
	{
		this.max_iter = max_iter;
		this.alpha = alpha;
		this.beta = beta;
		this.seed = seed;
		this.in_out = new int [2];
	}

	/**
	 * Native method (C->Fortran).
	 */
	private native long qap_graspd (int dim, int niter, float alpha_,
		float beta_, int look4, int dist[], int flow[], int sol[],
		int in_out_[]);
	// TODO remove runs parameter (not used)

	@Override
	long solveQAP (int dim, int dist[], int flow[], int sol[])
	{
		long cost;
		int look4 = -1;

		this.in_out[0] = seed;

		cost = qap_graspd (dim, max_iter, alpha, beta, look4,
							dist, flow, sol, in_out);
		
		this.seed = in_out[0];
		this.last_num_iter = in_out[1];

		return cost;
	}

	public void setAlpha (float alpha)
	{
		this.alpha = alpha;
	}

	public float getAlpha ()
	{
		return this.alpha;
	}

	public void setBeta (float beta)
	{
		this.beta = beta;
	}

	public float getBeta ()
	{
		return this.beta;
	}

	public void setSeed (int seed)
	{
		this.seed = seed;
	}

	public int getCurrentSeed ()
	{
		return this.seed;
	}

	public int getLastNumberOfIteractions ()
	{
		return this.last_num_iter;
	}
}
