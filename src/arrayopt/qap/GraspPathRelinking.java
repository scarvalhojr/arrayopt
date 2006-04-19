/*
 * GraspPathRelinking.java
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
 * ÃrrayOpt is distributed in the hope that it will be useful, but WITHOUT ANY
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
public class GraspPathRelinking extends QAPSolverAlgorithm
{
	static
	{
		String lib_name = "qap_grasppr";

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
	
	protected int runs;
	
	protected int max_iter;

	protected float alpha;

	protected float beta;

	protected int seed;
	
	protected int elite_size;

	protected int last_num_iter;

	protected int in_out[];
	
	public static final int DEFAULT_RUNS = 1;

	public static final int DEFAULT_MAX_ITERACTIONS = 32;

	public static final float DEFAULT_ALPHA = .1f;

	public static final float DEFAULT_BETA = .4f;
	
	public static final int DEFAULT_ELITE_SIZE = 10;

	public static final int DEFAULT_SEED = 270001;

	public GraspPathRelinking ()
	{
		this (DEFAULT_RUNS, DEFAULT_MAX_ITERACTIONS);
	}

	public GraspPathRelinking (int runs, int max_iter)
	{
		this (runs, max_iter, DEFAULT_ALPHA, DEFAULT_BETA, DEFAULT_ELITE_SIZE);
	}

	public GraspPathRelinking (int runs, int max_iter, float alpha, float beta,
		int elite_size)
	{
		this.runs = runs;
		this.max_iter = max_iter;
		this.alpha = alpha;
		this.beta = beta;
		
		// TODO use random seed?
		// this.seed = (int) (Short.MAX_VALUE * Math.random());
		this.seed = DEFAULT_SEED;
		
		this.elite_size = elite_size;
		this.in_out = new int [2];
	}

	/**
	 * Native method (in C).
	 */
	private native long qap_grasppr (int dim, int flow[], int dist[],
		int sol[], float alpha_, float beta_, int runs_, int max_itr_,
		int look4, int elite_size_, int max_time_, int in_out_[]);

	/**
	 *
	 */
	@Override
	long solveQAP (int dim, int dist[], int flow[], int sol[])
	{
		long cost;
		int look4 = -1;
		int max_time = 0;

		// TODO remove this as the C code does not return an updated seed
		this.in_out[0] = seed;

		cost = qap_grasppr (dim, flow, dist, sol, alpha, beta, runs, max_iter,
							look4, elite_size, max_time, in_out);

		// TODO generate a new seed for the next call?
		// this.seed = (int) (Short.MAX_VALUE * Math.random());
		
		this.last_num_iter = in_out[1];
		
		return cost;
	}

	public void setPhase1Parameters (float alpha, float beta)
	{
		this.alpha = alpha;
		this.beta = beta;
	}

	public float getAlpha ()
	{
		return this.alpha;
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
