/*
 * SlidingWindowOptimization.java
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

package arrayopt.layout;

/**
 *
 */
public class SlidingWindowOptimization implements PostPlacementAlgorithm
{
	/**
	 * document this
	 */
	protected IteractiveOptimizationAlgorithm alg;

	/**
	 * document this
	 */
	protected RectangularRegion region;

	/**
	 * document this
	 *
	 * how many rows (or columns) the window is shifted
	 */
	protected int shift;

	/**
	 * document this
	 *
	 * maximum number of iteractions that must be performed
	 */
	protected int max_iter;

	/**
	 * document this
	 *
	 * default maximum number of iteractions (an iteraction is a full pass over
	 * the chip surface)
	 */
	protected static final int DEFAULT_MAX_ITER = 10;

	/**
	 * document this
	 *
	 * can be specified as to stop the iteractions if improvement falls below it
	 */
	protected int threshold;

	/**
	 * document this
	 *
	 * default threshold
	 */
	protected static final float DEFAULT_THRESHOLD = .03f;

	/**
	 * document this
	 *
	 * dimension of sliding window and shift are required; no default values
	 * because of different QAP solvers' caracteristics)
	 */
	public SlidingWindowOptimization (IteractiveOptimizationAlgorithm alg,
		int window_dim, int shift)
	{
		// to do
	}

	/**
	 * document this
	 */
	public SlidingWindowOptimization (IteractiveOptimizationAlgorithm alg,
		int window_dim, int shift, int max_iter)
	{
		// to do
	}

	/**
	 * document this
	 */
	public SlidingWindowOptimization (IteractiveOptimizationAlgorithm alg,
		 int window_dim, int shift, int max_iter, float threshold)
	{
		// to do

		// validate arguments

		// instantiate arrays to appropriate dimension
		// and initialize instance variables

		// instantiate rectangular region object, initially pointing to
		// coordinates (0,0) with the requested dimension (this object
		// will be reused as the window slides over the chip surface)
	}

	/**
	 * document this
	 */
	public synchronized void optimizeLayout (Chip chip)
	{
		// loop for desired number of iteraction

			// slide window from left to right, top to bottom

				// update window coordinates

				// run iteractive algorithm

				// check improvement

			// quit if maximum improvement of last iteraction
			// (full pass over the chip) falls below threshold
	}
}
