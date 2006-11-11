/*
 * KThreadingPlacer.java
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
public class KThreadingPlacer implements PlacementAlgorithm, FillingAlgorithm
{
	/**
	 * TODO document this
	 */
	private int kvalue;
	
	/**
	 * TODO document this
	 */
	private boolean sort_embeddings;

	/**
	 * TODO document this
	 */
	public KThreadingPlacer (int kvalue)
	{
		this(kvalue, false);
	}

	/**
	 * TODO document this
	 */
	public KThreadingPlacer (int kvalue, boolean sort_embeddings)
	{
		this.kvalue = kvalue;
		this.sort_embeddings = sort_embeddings;
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

		if (chip instanceof SimpleChip)
			return fillRegion ((SimpleChip) chip, r, probe_id, start, end);

		else if (chip instanceof AffymetrixChip)
			return fillRegion ((AffymetrixChip) chip, r, probe_id, start, end);

		else
			throw new IllegalArgumentException ("Unsupported chip type.");
	}

	/**
	 *
	 */
	protected int fillRegion (SimpleChip chip, RectangularRegion region,
		int probe_id[], int start, int end)
	{
		int row, r, c, dir = -1;
		int delta, move, UP = 0, DOWN = 1;
		
		for (row = region.first_row; row <= region.last_row; row += kvalue + 1)
		{
			// alternates filling direction
			if ((dir = -dir) == +1)
				// left to right
				c = region.first_col;
			else
				// right to left
				c = region.last_row;
			
			// k-threading
			move = DOWN;
			delta = -1;
			
			while (true)
			{
				// k-threading
				if (move == DOWN)
				{
					if (delta == kvalue)
					{
						c += dir;
						move = UP;
					}
					else
						delta++;
				}
				else // (move == UP)
				{
					if (delta == 0)
					{
						c += dir;
						move = DOWN;
					}
					else
						delta--;
				}
				
				// stop when column gets out of region
				if (c < region.first_col || c > region.last_col)
					break;
				
				// compute row with k-threading
				r = row + delta;
				
				// skip if k-threading row gets out of the region
				if (r > region.last_row)
					continue;
				
				// skip non-empty spot
				if (chip.spot[r][c] != Chip.EMPTY_SPOT)
					continue;
				
				// skip fixed spot
				if (chip.isFixedSpot(r, c))
					continue;
				
				// place first probe of the list 
				chip.spot[r][c] = probe_id[start];

				// advance list pointer
				if (++start > end)
					// return if reached end of the list
					return 0;					
			}
		}

		// some probes could not be placed
		return end - start + 1;
	}

	/**
	 *
	 */
	protected int fillRegion (AffymetrixChip chip, RectangularRegion region,
		int probe_id[], int start, int end)
	{
		// TODO implement this
		// look at this.fillRegion (SimpleChip...)
		//     and Sequential.fillRegion (AffymetrixChip...)
		
		return 2 * (end - start + 1);
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
