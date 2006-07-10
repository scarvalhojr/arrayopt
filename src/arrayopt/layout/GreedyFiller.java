/*
 * GreedyFiller.java
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

import arrayopt.util.ArrayIndexedCollection;
import arrayopt.util.QuickSort;

/**
 *
 */
public class GreedyFiller implements PlacementAlgorithm, FillingAlgorithm
{
	private int window_size;
	
	private int mode;
	
	public static final int MODE_BORDER_LENGTH = 0;
	
	public static final int MODE_CONFLICT_INDEX = 1;
	
	private static final int NO_PREPROCESSING = 0;
	
	public static final int SORT_EMBEDDINGS = 1;

	public static final int RANDOMIZE_INPUT = 2;
	
	private boolean sort_embeddings;
	
	private boolean randomize_input;
	
	private RectangularRegion chip_region;
	
	private int embed_len;
	
	private int probe_len;
	
	private int ci_dim;
	
	private double pos_weight[];
	
	private double conflict[];

	public GreedyFiller ()
	{
		this(MODE_BORDER_LENGTH);
	}

	public GreedyFiller (int mode)
	{
		this(mode, 0);
	}

	public GreedyFiller (int mode, int window_size)
	{
		this(mode, window_size, NO_PREPROCESSING);
	}

	public GreedyFiller (int mode, int window_size, int options)
	{
		switch (mode)
		{
			case MODE_BORDER_LENGTH:
			case MODE_CONFLICT_INDEX:
				this.mode = mode;
				break;
				
			default:
				throw new IllegalArgumentException
					("Illegal value for argument 'mode'.");
		}

		// set pre-processing options
		sort_embeddings = (options == SORT_EMBEDDINGS ? true : false);
		randomize_input = (options == RANDOMIZE_INPUT ? true : false);
				
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
		
		if (sort_embeddings)
		{
			// sort embeddings lexicographically (as binary strings)
			QuickSort.sort(new EmbeddingSort(chip, probe_id),
							start, end - start +1);
		}
		else if (randomize_input)
		{
			randomizeInput (probe_id, start, end);
		}

		if (chip instanceof SimpleChip)
		{
			if (mode == MODE_CONFLICT_INDEX)
			{
				// prepare for faster conflict index calculations
				ci_dim = ConflictIndex.dimConflictRegion();
				chip_region = chip.getChipRegion();
				embed_len = chip.getEmbeddingLength();
				probe_len = chip.getProbeLength();
				
				if (conflict == null)
				{
					conflict = new double [embed_len];
					pos_weight = new double [probe_len + 1];
				}
				else if (conflict.length != embed_len)
				{
					conflict = new double [embed_len];
					pos_weight = new double [probe_len + 1];
				}
				
				// store position weights locally
				for (int b = 0; b <= probe_len; b++)
					pos_weight[b] = ConflictIndex.positionWeight(b, probe_len);
			}
			
			return fillRegion ((SimpleChip) chip, r, probe_id, start, end);
		}
		else if (chip instanceof AffymetrixChip)
		{
			// TODO implement fast conflict index calculation for Affymetrix
			
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
		int s, tmp;
		
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

				// find probe with minimum cost
				if (mode == MODE_BORDER_LENGTH)
					s = findMinBorderLength (chip, r, c, probe_id, start, end);
				else // (mode == MODE_CONFLICT_INDEX)
					s = findMinConflictIndex (chip, r, c, probe_id, start, end);

				// move selected probe to the front of the list
				tmp = probe_id[start];
				probe_id[start] = probe_id[s];
				probe_id[s] = tmp;
				
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
		int s, tmp;
		
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

				// find probe pair with minimum cost
				if (mode == MODE_BORDER_LENGTH)
					s = findMinBorderLength (chip, r, c, probe_id, start, end);
				else // (mode == MODE_CONFLICT_INDEX)
					s = findMinConflictIndex (chip, r, c, probe_id, start, end);

				// move selected probe pair to the front of the list
				tmp = probe_id[start];
				probe_id[start] = probe_id[s];
				probe_id[s] = tmp;
				
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

	private int findMinBorderLength (SimpleChip chip, int row, int col,
		int probe_id[], int start, int end)
	{
		long	cost, min;
		int		best;

		// if window size is limited
		if (window_size > 0)
		{
			// limit search space if probe list is larger than window
			end = end - start + 1 > window_size ? start + window_size - 1 : end;  
		}

		min = LayoutEvaluation.borderLength (chip, row, col, probe_id[start]);
		if (min == 0) return start;
		
		best = start;

		for (int i = start + 1; i <= end; i++)
		{
			cost = LayoutEvaluation.borderLength (chip, row, col, probe_id[i]);
			
			if (cost < min)
			{
				min = cost;
				best = i;
			}
		}
		
		return best;
	}

	private int findMinConflictIndex (SimpleChip chip, int row, int col,
			int probe_id[], int start, int end)
	{
		boolean	empty;
		double	cost, min;
		int		best;
		
		// prepare conflict index costs
		empty = examineNeighbors (chip, row, col);
		
		// if region around the spot is empty,
		// place any probe and return
		if (empty) return start;

		// if window size is limited
		if (window_size > 0)
		{
			// limit search space if probe list is larger than window
			end = end - start + 1 > window_size ? start + window_size - 1 : end;  
		}

		min = conflictIndex (chip, probe_id[start], Double.POSITIVE_INFINITY);
		
		best = start;

		for (int i = start + 1; i <= end; i++)
		{
			cost = conflictIndex(chip, probe_id[i], min);

			if (cost < min)
			{
				min = cost;
				best = i;
			}
		}
		
		return best;
	}
	
	private boolean examineNeighbors (SimpleChip chip, int row, int col)
	{
		boolean empty = true;
		double delta;
		int r, c, id, step, word, bitmask = 0;
		
		// reset costs
		for (step = 0; step < embed_len; step++)
			conflict[step] = 0;
		
		// define region around the spot that needs to be examined
		int min_row = Math.max(row - ci_dim, chip_region.first_row);
		int max_row = Math.min(row + ci_dim, chip_region.last_row);
		int min_col = Math.max(col - ci_dim, chip_region.first_col);
		int max_col = Math.min(col + ci_dim, chip_region.last_col);
		
		for (r = min_row; r <= max_row; r++)
		{
			for (c = min_col; c <= max_col; c++)
			{
				// skip if neighbor is empty
				if ((id = chip.spot[r][c]) == Chip.EMPTY_SPOT)
					continue;
				
				// get distance-dependent weight
				delta = ConflictIndex.distanceWeight(r, c, row, col);
				
				if (delta == 0)
					continue;
				
				empty = false;
				
				for (step = 0, word = - 1; step < embed_len; step++)
				{
					if (step % Integer.SIZE == 0)
					{
						word++;
						bitmask = 0x01 << (Integer.SIZE - 1);
					}
					else
						bitmask >>>= 1;

					// check state of embedding at current step
					if ((chip.embed[id][word] & bitmask) != 0)
					{
						// spot is in an unmasked step
						conflict[step] += delta;
					}
				}
			}
		}

		return empty;
	}
	
	private double conflictIndex (SimpleChip chip, int id, double max)
	{
		int base, step, word, bitmask = 0;
		double ci = 0;
		
		for (base = 0, step = 0, word = - 1; step < embed_len; step++)
		{
			if (step % Integer.SIZE == 0)
			{
				word++;
				bitmask = 0x01 << (Integer.SIZE - 1);
			}
			else
				bitmask >>>= 1;
			
			// check embedding state at current step 
			if ((chip.embed[id][word] & bitmask) != 0)
			{
				// spot is in an unmasked step

				// increment the number of synthesized bases
				base++;
			}
			else
			{
				// masked step
				ci += pos_weight[base] * conflict[step];
				
				// stop if CI exceeds limit
				if (ci > max) break;
			}
		}

		return ci;
	}

	private int findMinBorderLength (AffymetrixChip chip, int row, int col,
			int probe_id[], int start, int end)
	{
		long	cost, min;
		int		best;

		// if window size is limited
		if (window_size > 0)
		{
			// limit search space if probe list is larger than window
			end = end - start + 1 > window_size ? start + window_size - 1 : end;  
		}

		min = LayoutEvaluation.borderLength (chip, row, col, probe_id[start]);
		if (min == 0) return start;
		
		best = start;

		for (int i = start + 1; i <= end; i++)
		{
			cost = LayoutEvaluation.borderLength (chip, row, col, probe_id[i]);
			
			if (cost < min)
			{
				min = cost;
				best = i;
			}
		}
		
		return best;
	}

	private int findMinConflictIndex (AffymetrixChip chip, int row, int col,
			int probe_id[], int start, int end)
	{
		double	cost, min;
		int		best;

		// if window size is limited
		if (window_size > 0)
		{
			// limit search space if probe list is larger than window
			end = end - start + 1 > window_size ? start + window_size - 1 : end;  
		}

		min = LayoutEvaluation.conflictIndex (chip, row, col, probe_id[start]);
		if (min == 0) return start;
		
		best = start;

		for (int i = start + 1; i <= end; i++)
		{
			cost = LayoutEvaluation.conflictIndex(chip, row, col, probe_id[i]);

			if (cost < min)
			{
				min = cost;
				best = i;
			}
		}
		
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
	
	private class EmbeddingSort implements ArrayIndexedCollection
	{
		private Chip chip;
		
		private int words;
		
		private int probe_id[];
		
		private int pivot;
		
		EmbeddingSort (Chip chip, int probe_id[])
		{
			this.chip = chip;
			this.words = chip.embed[0].length;
			this.probe_id = probe_id;
		}
		
		public int compare (int i, int j)
		{
			int id_i = probe_id[i];
			int id_j = probe_id[j];
			
			for (int w = 0; w < words; w++)
			{
				// compare first bit (signal)
				// first bit is 1 => negative number
				// first bit is 0 => non-negative
				if (chip.embed[id_i][w] < 0)
				{
					if (chip.embed[id_j][w] >= 0)
						return 1;
				}
				else
				{
					if (chip.embed[id_j][w] < 0)
						return -1;
				}
				
				// compare remaining bits if both have same signal
				if (chip.embed[id_i][w] < chip.embed[id_j][w])
					return -1;
				else if (chip.embed[id_i][w] > chip.embed[id_j][w])
					return 1;
			}
			
			return 0;
		}
		
		public void swap (int i, int j)
		{
			int tmp;
			tmp = probe_id[i];
			probe_id[i] = probe_id[j];
			probe_id[j] = tmp;
		}
		
		public void setPivot (int i)
		{
			this.pivot = probe_id[i];
		}
		
		public int compareToPivot (int i)
		{
			int id_i = probe_id[i];

			for (int w = 0; w < words; w++)
			{
				// compare first bit (signal)
				// first bit is 1 => negative number
				// first bit is 0 => non-negative
				if (chip.embed[id_i][w] < 0)
				{
					if (chip.embed[pivot][w] >= 0)
						return 1;
				}
				else
				{
					if (chip.embed[pivot][w] < 0)
						return -1;
				}
				
				// compare remaining bits if both have same signal
				if (chip.embed[id_i][w] < chip.embed[pivot][w])
					return -1;
				else if (chip.embed[id_i][w] > chip.embed[pivot][w])
					return 1;
			}
			
			return 0;
		}
		
		public int medianOfThree (int i, int j, int k)
		{
			int median;
			
			if (compare(i,j) < 0)
			{
				if (compare(j,k) < 0)
					median = j;
				else if (compare(i,k) < 0)
					median = k;
				else
					median = i;
			}
			else
			{
				if (compare(j,k) > 0)
					median = j;
				else if (compare(i,k) > 0)
					median = k;
				else
					median = i;
			}
			
			return median;
		}
	}
}
