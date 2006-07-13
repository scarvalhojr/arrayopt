/*
 * GreedyEmbeddingsPlacer.java
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
 * �rrayOpt is distributed in the hope that it will be useful, but WITHOUT ANY
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

import arrayopt.util.ArrayIndexedCollection;
import arrayopt.util.QuickSort;

/**
 *
 */
public class GreedyEmbeddingsPlacer implements PlacementAlgorithm, FillingAlgorithm
{
	private int window_size;
	
	private int mode;
	
	public static final int BORDER_LENGTH_MIN =
								OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN;
	
	public static final int CONFLICT_INDEX_MIN =
								OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN;
		
	private boolean randomize_input;
	
	private OptimumSingleProbeEmbedding ospe;
	
	private ProbeSorting probe_sort;
	
	private long probe_rank[];

	public GreedyEmbeddingsPlacer ()
	{
		this(BORDER_LENGTH_MIN);
	}

	public GreedyEmbeddingsPlacer (int mode)
	{
		this(mode, 0);
	}

	public GreedyEmbeddingsPlacer (int mode, int window_size)
	{
		this(mode, window_size, false);
	}

	public GreedyEmbeddingsPlacer (int mode, int window_size, boolean randomize)
	{
		switch (mode)
		{
			case BORDER_LENGTH_MIN:
			case CONFLICT_INDEX_MIN:
				this.mode = mode;
				break;
				
			default:
				throw new IllegalArgumentException
					("Illegal value for argument 'mode'.");
		}

		this.randomize_input = randomize;
		this.window_size = window_size;
	}

	/**
	 *
	 */
	public int makeLayout (Chip chip)
	{
		int		id[];

		// reset current layout (if any)
		chip.resetLayout();

		// get list of movable probes
		id = chip.getMovableProbes ();

		return fillRegion (chip, chip.getChipRegion(), id, 0, id.length - 1);
	}

	/**
	 *
	 */
	public int fillRegion (Chip chip, Region region, int probe_id[])
	{
		return fillRegion (chip, region, probe_id, 0, probe_id.length - 1);
	}

	/**
	 *
	 */
	public int fillRegion (Chip chip, Region region, int probe_id[], int start,
		int end)
	{
		RectangularRegion r;

		if (!(region instanceof RectangularRegion))
			throw new IllegalArgumentException
				("Only rectangular regions are supported.");

		r = (RectangularRegion) region;
		
		if (randomize_input)
		{
			randomizeInput (probe_id, start, end);
		}

		// create OSPE object
		ospe = OptimumSingleProbeEmbedding.createEmbedder(chip, mode);
		
		// sort probes lexicographically
		probe_rank = chip.computeProbeRanks(probe_id, 0, probe_id.length - 1);
		probe_sort = new ProbeSorting (probe_id, probe_rank);
		QuickSort.sort(probe_sort, start, end - start + 1);
		
		if (chip instanceof SimpleChip)
		{
			return fillRegion ((SimpleChip) chip, r, probe_id, start, end);
		}
		else if (chip instanceof AffymetrixChip)
		{
			// TODO this case needs to be tested!
			
			return fillRegion ((AffymetrixChip) chip, r, probe_id, start, end);
		}
		else
			throw new IllegalArgumentException ("Unsupported chip type.");
	}

	/**
	 *
	 */
	private int fillRegion (SimpleChip chip, RectangularRegion region,
		int probe_id[], int start, int end)
	{
		int opt;
		
		// TODO place first a pivot (probe with min number of embeddings)
		
		for (int r = region.first_row; r <= region.last_row; r ++)
		{
			for (int c = region.first_col; c <= region.last_col; c++)
			{
				// skip spot if not empty
				if (chip.spot[r][c] != Chip.EMPTY_SPOT)
					continue;

				// skip fixed spot
				if (chip.isFixedSpot(r, c))
					continue;

				// find probe whose embedding produce minimum conflicts
				opt = findOptimalEmbedding (r, c, probe_id, start, end);
				
				// swap selected probe with the first element of the list
				probe_sort.swap(start, opt);
				
				// place selected probe
				chip.spot[r][c] = probe_id[start];
				
				start ++;

				if (start > end)
					// all probes were placed
					return 0;
			}
		}

		// some probe could not be placed
		return (end - start + 1);
	}

	/**
	 *
	 */
	private int fillRegion (AffymetrixChip chip, RectangularRegion region,
		int probe_id[], int start, int end)
	{
		int opt;
		
		// TODO place first a pivot (probe with min number of embeddings)
		
		for (int r = region.first_row; r < region.last_row; r ++)
		{
			for (int c = region.first_col; c <= region.last_col; c++)
			{
				// skip spot pair if not empty
				if (chip.spot[r][c] != Chip.EMPTY_SPOT ||
					chip.spot[r+1][c] != Chip.EMPTY_SPOT)
					continue;

				// skip fixed spot pair
				if (chip.isFixedSpot(r, c) || chip.isFixedSpot(r+1, c))
					continue;

				// find probe pair whose embeddings produce minimum conflicts
				opt = findOptimalEmbedding (r, c, probe_id, start, end);

				// swap selected probe with the first element of the list
				probe_sort.swap(start, opt);
				
				// place selected probe pair
				chip.spot[r][c] = probe_id[start];
				chip.spot[r+1][c] = probe_id[start] + 1;
				
				start ++;

				if (start > end)
					// all probes were placed
					return 0;
			}
		}

		// some probe pairs could not be placed
		return 2 * (end - start + 1);
	}

	private int findOptimalEmbedding (int row, int col, int probe_id[],
			int start, int end)
	{
		double	cost, min;
		int		best;
		
		if (window_size > 0)
		{
			// limit search space if probe list is larger than window
			end = end - start + 1 > window_size ? start + window_size - 1 : end;
		}
		
		// compute cost of placing first element in the spot
		min = ospe.minDistanceSpot(row, col, probe_id[start]);
		
		// if min conflict is zero, probably the spot has only empty neighbors,
		// and all probes will have zero cost 
		if (min <= 0) return start;
		
		best = start;
		for (int i = start + 1; i <= end; i++)
		{
			cost = ospe.minDistanceProbe(probe_id[i]);

			if (cost < min)
			{
				min = cost;
				best = i;
			}
		}
		
		// re-embed best probe optimally
		ospe.reembedOptimally(probe_id[best]);
		
		return best;
	}
	
	
	private void randomizeInput (int probe_id[], int start, int end)
	{
		int tmp, rand;
		
		for (int i = start; i < end; i++)
		{
			// select a random element of the list
			rand = i + (int) ((end - i + 1) * Math.random());

			// swap the selected element with
			// the one at the current position
			tmp = probe_id[i];
			probe_id[i] = probe_id[rand];
			probe_id[rand] = tmp;
		}		
	}
	
	private class ProbeSorting implements ArrayIndexedCollection
	{
		private int probe_id[];
		
		private long rank[];
		
		private long pivot;
		
		ProbeSorting (int probe_id[], long probe_rank[])
		{
			this.probe_id = probe_id;
			this.rank = probe_rank;
		}
		
		public int compare (int i, int j)
		{
			return rank[i] < rank[j] ? -1 :
					rank[i] == rank[j] ? 0 : +1;
		}
		
		public void swap (int i, int j)
		{
			int tmp1;
			tmp1 = probe_id[i];
			probe_id[i] = probe_id[j];
			probe_id[j] = tmp1;
			
			long tmp2;
			tmp2 = rank[i];
			rank[i] = rank[j];
			rank[j] = tmp2;
		}
		
		public void setPivot (int i)
		{
			this.pivot = rank[i];
		}
		
		public int compareToPivot (int i)
		{
			return rank[i] < pivot ? -1 :
					rank[i] == pivot ? 0 : +1;
		}
		
		public int medianOfThree (int i, int j, int k)
		{
			long rank_i = rank[i];
			long rank_j = rank[j];
			long rank_k = rank[k];

			return rank_i < rank_j ?
					(rank_j < rank_k ? j : (rank_i < rank_k ? k : i)) :
					(rank_j > rank_k ? j : (rank_i > rank_k ? k : i));
		}
	}
}