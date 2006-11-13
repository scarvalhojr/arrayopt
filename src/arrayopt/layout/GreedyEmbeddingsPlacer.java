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
public class GreedyEmbeddingsPlacer implements PlacementAlgorithm, FillingAlgorithm
{	
	private int window_size;
	
	private int mode;
	
	public static final int BORDER_LENGTH_MIN =
								OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN;
	
	public static final int CONFLICT_INDEX_MIN =
								OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN;
		
	private static final int NO_PREPROCESSING = 0;
	
	public static final int SORT_PROBES = 1;

	public static final int RANDOMIZE_INPUT = 2;
	
	private boolean sort_probes;
	
	private boolean randomize_input;
	
	private OptimumSingleProbeEmbedding ospe;
	
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
		this(mode, window_size, NO_PREPROCESSING);
	}

	public GreedyEmbeddingsPlacer (int mode, int window_size, int options)
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

		// set pre-processing options
		this.sort_probes = (options == SORT_PROBES ? true : false);
		this.randomize_input = (options == RANDOMIZE_INPUT ? true : false);

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
		MyLinkedList head, prev, curr;
		RectangularRegion r;
		
		if (end < start) return 0;

		if (!(region instanceof RectangularRegion))
			throw new IllegalArgumentException
				("Only rectangular regions are supported.");

		r = (RectangularRegion) region;

		// create OSPE object
		ospe = OptimumSingleProbeEmbedding.createEmbedder(chip, mode);

		if (sort_probes)
		{
			// sort probes lexicographically
			long probe_rank[] = chip.computeProbeRanks(probe_id, 0,
					probe_id.length - 1);
			QuickSort.sort(new ProbeSorting (probe_id, probe_rank), start,
					end - start + 1);
		}
		else if (randomize_input)
		{
			randomizeInput (probe_id, start, end);
		}

		// create a linked list with the probe IDs
		head = prev = new MyLinkedList (probe_id[start], null);
		for (int i = start + 1; i <= end; i++)
		{
			curr = new MyLinkedList (probe_id[i], prev);
			prev.next = curr;
			prev = curr;
		}

		if (chip instanceof SimpleChip)
		{
			return fillRegion ((SimpleChip) chip, r, head);
		}		
		else if (chip instanceof AffymetrixChip)
		{
			// TODO this case needs to be tested!
			return fillRegion ((AffymetrixChip) chip, r, head);
		}
		else
			throw new IllegalArgumentException ("Unsupported chip type.");
	}

	/**
	 *
	 */
	private int fillRegion (SimpleChip chip, RectangularRegion region,
			MyLinkedList list)
	{
		MyLinkedList n;
		int count;
		
		// TODO place first a pivot (probe with min number of embeddings)
		
		// TODO fill spots using k-threading 

		for (int r = region.first_row; r <= region.last_row; r ++)
		{
			for (int c = region.first_col; c <= region.last_col; c++)
			{
				// if spot is not empty, re-embed its probe optimally
				if (chip.spot[r][c] != Chip.EMPTY_SPOT)
				{
					ospe.reembedSpot(r, c);
					continue;
				}

				// skip (empty) fixed spots
				if (chip.isFixedSpot(r, c))
					continue;

				// place probe whose embedding produce minimum conflicts
				list = placeOptimalEmbedding (chip, r, c, list);
				
				if (list == null)
					// all probes were placed
					return 0;
			}
		}

		// count how many elements are left in the list
		count = 1;
		for (n = list.prev; n != null; n = n.prev) count++;
		for (n = list.next; n != null; n = n.next) count++;
		
		return count;
	}

	/**
	 *
	 */
	private int fillRegion (AffymetrixChip chip, RectangularRegion region,
			MyLinkedList list)
	{
		int count;
		
		// TODO re-write with new algorithm (using double-linked list)
		return 0;
		
		/*
		for (int r = region.first_row; r < region.last_row; r ++)
		{
			for (int c = region.first_col; c <= region.last_col; c++)
			{
				// if spot pair is not empty, re-embed its probes optimally
				if (chip.spot[r][c] != Chip.EMPTY_SPOT ||
					chip.spot[r+1][c] != Chip.EMPTY_SPOT)
				{
					// TODO this needs to be tested!
					// ospe.reembedSpot(r,c);
					continue;
				}

				// skip (empty) fixed spot pair
				if (chip.isFixedSpot(r, c) || chip.isFixedSpot(r+1, c))
					continue;

				// place probe pair whose embedding produce minimum conflicts
				list = placeOptimalEmbedding (chip, r, c, list);
				
				if (list == null)
					// all probes were placed
					return 0;
			}
		}

		// count how many elements are left in the list
		for (count = 0; list != null; list = list.next)
			count++;
		
		return count;
		*/
	}

	private MyLinkedList placeOptimalEmbedding (SimpleChip chip, int row,
			int col, MyLinkedList list)
	{
		MyLinkedList best, node;
		double cost, min;
		int count;

		// first check if there are (window_size/2) elements
		// to the right of the last placed probe
		count = (int) Math.floor(window_size / (double) 2);
		for (node = list; node.next != null && count > 0; count--)
			node = node.next;

		// we want to start the search (window_size/2 - 1) elements
		// to the left of the last placed probe
		
		// but if count > 0, there are only (window_size/2 - count)
		// elements to the right of the last placed probe, so
		// we compensate by starting the search (window_size/2 - 1 + count)
		// to the left of the last placed probe
		count = count + (int) Math.ceil(window_size / (double) 2) - 1; 
		for (node = list; node.prev != null && count > 0; count--)
			node = node.prev;
		
		// compute cost of placing first element in the spot
		min = ospe.minDistanceSpot(row, col, node.info);
		best = node;
		
		if (min > 0)
		{
			node = node.next;
			
			for (count = 1; node != null && count < window_size; count++)
			{
				if ((cost = ospe.minDistanceProbe(node.info, min)) < min)
				{
					min = cost;
					best = node;
				}

				node = node.next;
			}

			// re-embed best probe optimally
			ospe.reembedProbe(best.info);
		}
		
		// place best probe on the spot
		chip.spot[row][col] = best.info;
		
		node = null;
		
		// and delete it from the list
		if (best.next != null)
		{
			best.next.prev = best.prev;
			node = best.next;
		}

		if (best.prev != null)
		{
			best.prev.next = best.next;
			node = best.prev;
		}
		
		return node;
	}
	
	private MyLinkedList placeOptimalEmbedding (AffymetrixChip chip, int row,
			int col, MyLinkedList head)
	{
		
		// TODO re-write with new algorithm (using double-linked list)
		
		return head;
		
		/*
		MyLinkedList best, last, curr;
		double cost, min;
		int best_id, count = 1;
		
		// compute cost of placing first element in the spot
		min = ospe.minDistanceSpot(row, col, head.info);
		
		best = null;
		last = head;
		curr = head.next;
		
		if (min > 0)
		{
			while (curr != null)
			{
				cost = ospe.minDistanceProbe(curr.info);
				
				if (cost < min)
				{
					min = cost;
					best = last;
				}
				
				if (++count > window_size) break;
				
				last = curr;
				curr = last.next;
			}
		}
		
		if (best == null)
		{
			best_id = head.info;
			head = head.next;
		}
		else
		{
			best_id = best.next.info;
			best.next = best.next.next;
		}
		
		// place best probe on the spot
		chip.spot[row][col] = best_id;
		chip.spot[row + 1][col] = best_id + 1;
		
		// and re-embed them optimally
		ospe.reembedProbe(best_id);
		
		return head;
		*/
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
	
	private class MyLinkedList
	{
		MyLinkedList prev;
		MyLinkedList next;
		
		int info;
		
		MyLinkedList (int info, MyLinkedList prev)
		{
			this.info = info;
			this.prev = prev;
			this.next = null;
		}
	}
}
