/*
 * QAPSolverAlgorithm.java
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
 * This class defines an algorithm able to solve a quadratic assignment problem
 * (QAP). It serves mainly to define the minimum functionality that any such
 * implementation must provide.
 * 
 * <P>A classic example of a QAP is the facility location problem, which is the
 * problem of assigning facilities to locations. The facility location problem
 * has two components, a flow matrix (which stores the flow of materials between
 * the facilities) and a distance matrix (which contains the distance between
 * the locations). The flow of materials has a cost, and thus the problem is to
 * find an assignment of facilities to locations with minimum cost.</P>
 * 
 * <P>Every concrete sub-class must be able to solve a QAP instance that
 * consists of the flow and distance matrices passed as integer linear arrays.
 * This functionality is defined by the {@link #solveQAP} method. The result
 * must be a permutation of the integer numbers that corresponds to an
 * assignment of facilities to locations.</P>
 * 
 * <P>Sub-classes can implement their own solver or call an external method,
 * perhaps implemented in another language. Solvers can be exact (return the
 * optimal solution) or heuristic (return a good but not necessarily optimal
 * solution). However, QAP is notoriously hard and exact solutions are unlikely
 * to be found in a resonable time for large instances.</P>
 * 
 * <P>This class also provides two simple functionalites: computing the cost of
 * a given permutation ({@link #computeCost}) and swapping the flow and distance
 * matrices before calling the {@link #solveQAP} method. The latter may be
 * useful in several situations depending on the underlying solver and
 * characteristics of the matrices such as their sparsity (number of zero
 * entries).</P> 
 * 
 * @author Sergio A. de Carvalho Jr.
 */
public abstract class QAPSolverAlgorithm
{
	/**
	 * Swap the distance and flow matrices before calling the internal
	 * {@link #solveQAP} method.
	 */
	private boolean swap = false;

	/**
	 * Set or unset the swap of flow and distance matrices before the QAP is
	 * solved.
	 * 
	 * @param swap true if matrices should be swapped, false otherwise
	 */
	public void setSwapMatrices (boolean swap)
	{
		this.swap = swap;
	}

	/**
	 * Checks whether the swap of flow and distance matrices is enabled or not.
	 * 
	 * @return true is swap is enabled, false otherwise
	 */
	public boolean isSwapMatrices ()
	{
		return swap;
	}
	
	/**
	 * Solve an instance of a QAP. This method calls the (abstract) method
	 * {@link #solveQAP} that sub-classes must provide. If swap is enabled,
	 * the flow and distance matrices are swapped before the call. In this case,
	 * the permutation returned is also inverted to preserve its correct
	 * meaning.   
	 * 
	 * @param dim dimension of QAP
	 * @param dist distance matrix as an integer array
	 * @param flow flow matrix an integer array
	 * @param sol permutation with the QAP solution (returned)
	 * @return cost of (best) solution found
	 */
	public long solve (int dim, int dist[], int flow[], int sol[])
	{
		if (swap)
		{
			int tmp[] = new int [sol.length];
			
			// use tmp permutation, and swap flow and distance matrices
			long cost = solveQAP (dim, flow, dist, tmp);
			
			// invert permutation
			for (int i = 0; i < dim; i++)
				sol[tmp[i]] = i;
			
			return cost;
		}
		
		// else: no swap
		return solveQAP (dim, dist, flow, sol);
	}

	/**
	 * Internal method that sub-classes must implement. This is the method that
	 * actually must solve the QAP instance.
	 * 
	 * @param dim dimension of QAP
	 * @param dist distance matrix as an integer array
	 * @param flow flow matrix an integer array
	 * @param sol permutation with the QAP solution (returned)
	 * @return cost of (best) solution found
	 */
	abstract long solveQAP (int dim, int dist[], int flow[], int sol[]);

	/**
	 * Computes the cost of a solution for the given QAP instance.
	 *  
	 * @param dim dimension of QAP
	 * @param dist distance matrix as an integer array
	 * @param flow flow matrix an integer array
	 * @param sol permutation with the proposed QAP solution
	 * @return cost of the given solution
	 */
	public long computeCost (int dim, int dist[], int flow[], int sol[])
	{
		long cost = 0;

		for (int i = 0; i < dim; i++)
			for (int j = 0; j < dim; j++)
				cost += dist[i * dim + j] * flow[sol[i] * dim + sol[j]];

		return cost;
	}
}
