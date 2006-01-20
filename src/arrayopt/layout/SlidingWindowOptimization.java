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
	protected IteractiveOptimizationAlgorithm optimizer;

	/**
	 * document this
	 */
	protected RectangularRegion window;

	/**
	 * document this
	 */
	protected int window_dim;

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
	protected float threshold;

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
		this (alg, window_dim, shift, DEFAULT_MAX_ITER);
	}

	/**
	 * document this
	 */
	public SlidingWindowOptimization (IteractiveOptimizationAlgorithm alg,
		int window_dim, int shift, int max_iter)
	{
		this (alg, window_dim, shift, max_iter, DEFAULT_THRESHOLD);
	}

	/**
	 * document this
	 */
	public SlidingWindowOptimization (IteractiveOptimizationAlgorithm alg,
		 int window_dim, int shift, int max_iter, float threshold)
	{
		if (window_dim < 2)
			throw new IllegalArgumentException
				("Window dimension cannot be smaller than 2.");

		if (shift < 1)
			throw new IllegalArgumentException
				("Shift cannot be smaller than 1.");

		if (max_iter < 1)
			throw new IllegalArgumentException
				("Maximum number of iteractions cannot be smaller than 1.");

		this.optimizer = alg;
		this.window_dim = window_dim;
		this.shift = shift;
		this.max_iter = max_iter;
		this.threshold = threshold;
		
		// instantiate the window as a rectangular region object
		// (appropriate coordinates will be set when needed)
		window = new RectangularRegion (0, 0, 0, 0);
	}

	/**
	 * document this
	 */
	public synchronized void optimizeLayout (Chip chip)
	{
		RectangularRegion	region;
		float				impr, total_impr = 0;
		int					ncols, nrows, ncalls;
		
		region = chip.getChipRegion();

		// window dimension must not be greater than the chip's area
		if (window_dim > region.last_row - region.first_row + 1 ||
			window_dim > region.last_col - region.first_col + 1)
			throw new IllegalStateException
				("Dimension of sliding window is larger than the chip's region.");
				
		// compute how many times the algorithm must run
		// to fully cover the chip's surface
		ncols = region.last_col - region.first_col + 1;
		nrows = region.last_row - region.first_row + 1;
		ncalls  = 1 + (int) Math.ceil((ncols - window_dim) / (float) shift);
		ncalls *= 1 + (int) Math.ceil((nrows - window_dim) / (float) shift);
		
		System.err.println("ncalls: " + ncalls);

		// set window's position to top left corner
		window.first_col = region.first_col;
		window.last_col = region.first_col + window_dim - 1;
		window.first_row = region.first_row;
		window.last_row = region.first_row + window_dim - 1;					

		// loop for desired number of iteraction
		for (int i = 0; i < max_iter;)
		{
			// call opt algorithm
			impr = optimizer.optimizeLayout (chip, window);
			
			// update total improvement
			total_impr += impr / ncalls;

			// check if can shift window to the right
			if (window.last_col != region.last_col)
			{
				// right shift is possible
				window.last_col += shift;
				
				// correct if went outside the chip's surface
				if (window.last_col > region.last_col)
				{
					window.last_col = region.last_col;
					window.first_col = region.last_col - window_dim + 1;
				}
				else
					window.first_col += shift;
			}
			else
			{
				// no, need to move window back to the left border
				window.first_col = region.first_col;
				window.last_col = region.first_col + window_dim - 1;
				
				// and now shift it downwards
				
				// but first check if not reached bottom already
				if (window.last_row == region.last_row)
				{
					// yes... then move window to the top
					window.first_row = region.first_row;
					window.last_row = region.first_row + window_dim - 1;					

					// this means another iteraction
					// has been completed
					i++;
					
					System.err.println ("Total improvement: " + total_impr);
					
					// quit if total improvement is insufficient
					if (total_impr < threshold) return;
					
					// otherwise just reset it
					total_impr = 0;
				}
				else
				{
					// ok, just shift it down
					window.last_row += shift;

					// correct if went outside the chip's surface
					if (window.last_row > region.last_row)
					{
						window.last_row = region.last_row;
						window.first_row = region.last_row - window_dim + 1;
					}
					else
						window.first_row += shift;
				}
			}
		}
	}
}
