/*
 * GreedyPlusPlacer.java
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
 * This class implements the Greedy Plus placement algorithm.
 * 
 * <P>The Greedy Plus algorithm creates a new layout by filling spots
 * sequentially using a k-threading pattern identical to the one implemented in
 * {@link KThreadingPlacer}. The amplitude of the k-threading is set with the
 * {@link #kvalue} variable (at instantiation time). For every spot
 * <CODE>s</CODE>, it finds a probe <CODE>p</CODE> with minimum cost to fill
 * <CODE>s</CODE>.</P>
 * 
 * <P>The difference between the Greedy Plus and the {@link GreedyPlacer}
 * algorithm is that the former considers all possible embeddings of a probe
 * candidate for filling a given spot. This is done by using the OSPE algorithm
 * ({@link OptimumSingleProbeEmbedding}). Spots are then filled in an attempt
 * to minimize the sum of border conflicts
 * ({@link OptimumSingleProbeEmbedding#BORDER_LENGTH_MIN}) or conflict indices
 * around the spot ({@link OptimumSingleProbeEmbedding#CONFLICT_INDEX_MIN}).</P>
 * 
 * <P>The algorithm looks at {@link #window_size} probe candidate for filling
 * each spot. All probes are initially lexicographically sorted and a
 * doubly-linked list is created, where each node of the list contains a single
 * probe ID. Once a probe is placed, it is removed from the list, and the next
 * search of probe candidates starts at the next or previous element in the
 * list. The search will examine at most {@link #window_size} elements,
 * preferably <CODE>window_size/2</CODE> to the left and
 * <CODE>window_size/2</CODE> to the right of the last placed probe.</P>
 * 
 * <P>The sorted list allows more candidates to be examined because many rows of
 * the dynamic programming matrix can be skipped if two candidates have a common
 * prefix. The sorted list also improves the chances of finding probes similar
 * to its neighbors to fill a given spot. A doubly-linked list is used to
 * maintain the sorting during the whole algorith.</P>
 * 
 * @author Sergio A. de Carvalho Jr.
 */
public class GreedyPlusPlacer implements LayoutAlgorithm, FillingAlgorithm
{	
	/**
	 * This variable stores the current minimization mode used by the algorithm.
	 * Possible values are {@link OptimumSingleProbeEmbedding#BORDER_LENGTH_MIN}
	 * and {@link OptimumSingleProbeEmbedding#CONFLICT_INDEX_MIN}. 
	 */
	private int mode;
	
	/**
	 * Maximum number of probe candidades considered for each spot.
	 */
	private int window_size;
	
	/**
	 * The amplitude of the k-threading. This gives the number of upward and
	 * downward movements between the right-to-left and left-to-right paths over
	 * the chip.
	 */
	private int kvalue;
	
	private OptimumSingleProbeEmbedding ospe;
	
	private ProbeOrderingAlgorithm ordering;
	
	/**
	 * Creates an instance of the Greedy Plus placement algorithm with the
	 * desired minimization mode, window size and k value.
	 * 
	 * @param mode minimization mode, either border length or conflict index
	 * @param window_size maximum number of candidades examined for each spot
	 * @param kvalue amplitude of k-threading
	 */
	public GreedyPlusPlacer (int mode, int window_size, int kvalue)
	{
		switch (mode)
		{
			case OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN:
			case OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN:
				this.mode = mode;
				break;
				
			default:
				throw new IllegalArgumentException
					("Illegal value for argument 'mode'.");
		}
		
		if (kvalue < 0)
			throw new IllegalArgumentException ("Invalid k value: " + kvalue);

		this.window_size = window_size;
		this.kvalue = kvalue;
		this.ordering = new SortedSequencesOrdering();
	}

	/**
	 * Creates a new layout of a microarray chip using the Greedy Plus placement
	 * algorithm.
	 * 
	 * @param chip instance of a microarray chip
	 */
	public void changeLayout (Chip chip)
	{
		int		id[];

		// reset current layout (if any)
		chip.resetLayout();

		// get list of movable probes
		id = chip.getMovableProbes ();

		fillRegion (chip, chip.getChipRegion(), id, 0, id.length - 1);
	}

	/**
	 * Fills the spots of a region with the given list of probes.
	 * 
	 * @param chip instance of a microarray chip
	 * @param region region to be filled
	 * @param probe_id list of probe IDs
	 */
	public int fillRegion (Chip chip, Region region, int probe_id[])
	{
		return fillRegion (chip, region, probe_id, 0, probe_id.length - 1);
	}

	/**
	 * Fills the spots of a region with the given list of probes delimited by
	 * a starting and ending positions.
	 * 
	 * @param chip instance of a microarray chip
	 * @param region region to be filled
	 * @param probe_id list of probe IDs
	 * @param start first element of the list
	 * @param end last element of the list
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

		// region to be filled
		r = (RectangularRegion) region;

		// create OSPE object
		ospe = OptimumSingleProbeEmbedding.createEmbedder(chip, mode);
		
		// sort probes lexicographically
		ordering.orderProbes(chip, probe_id, start, end);

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
		// else
			throw new IllegalArgumentException ("Unsupported chip type.");
	}

	private int fillRegion (SimpleChip chip, RectangularRegion region,
			MyLinkedList list)
	{
		int row, r, c, dir = -1;
		int delta, move, UP = 0, DOWN = 1;
		
		// start placement with a pivot (a probe with min number of embeddings)
		list = placePivot (chip, region.first_row, region.first_col, list);
		
		for (row = region.first_row; row <= region.last_row; row += kvalue + 1)
		{
			// alternates filling direction
			if ((dir = -dir) == +1)
				// left to right
				c = region.first_col;
			else
				// right to left
				c = region.last_col;
			
			// k-threading
			move = DOWN; delta = -1;
			
			while (true)
			{
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
				if (c < region.first_col || c > region.last_col) break;
				
				// skip if k-threading row gets out of the region
				if ((r = row + delta) > region.last_row) continue;
				
				// if spot is not empty, re-embed its probe optimally
				if (chip.spot[r][c] != Chip.EMPTY_SPOT)
				{
					ospe.reembedSpot(r, c);
					continue;
				}

				// skip (empty) fixed spots
				if (chip.isFixedSpot(r, c)) continue;

				// place probe with embedding resulting in minimum conflicts
				list = placeOptimalEmbedding (chip, r, c, list);
				
				if (list == null)
					// all probes were placed
					return 0;
			}
		}

		// count how many elements are left in the list
		MyLinkedList n;
		int count = 1;
		for (n = list.prev; n != null; n = n.prev) count++;
		for (n = list.next; n != null; n = n.next) count++;
		
		return count;
	}

	private MyLinkedList placeOptimalEmbedding (SimpleChip chip, int row,
			int col, MyLinkedList node)
	{
		MyLinkedList best;
		double	cost, min;
		int		count;
		
		// find node so that the search will examine 
		// window_size elements around the last placed probe 
		node = findStartingNode (node);
		
		// compute cost of placing first element in the spot
		min = ospe.minDistanceSpot(row, col, node.info);
		best = node;
		
		// TODO in case of ties, choose probe "closer" to last placed probe
		// (see code of GreedyPlacer.minConflictIndex) 
		
		if (min > 0)
		{
			node = node.next;
			
			for (count = 1; node != null && count < window_size; count++)
			{
				cost = ospe.minDistanceProbe(node.info, min);
				if (cost < min)
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

	private MyLinkedList placePivot (SimpleChip chip, int row, int col,
			MyLinkedList node)
	{
		MyLinkedList best;
		long num_embed, min;
		
		best = node;
		
		for (min = Long.MAX_VALUE; node != null; node = node.next)
		{
			num_embed = ospe.numberOfEmbeddings(node.info);
			
			if (num_embed == 1)
			{
				best = node;
				break;
			}
			
			if (num_embed < min)
			{
				best = node;
				min = num_embed;
			}
		}
		
		// place pivot on the spot
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
	
	private MyLinkedList findStartingNode (MyLinkedList list)
	{
		MyLinkedList node;
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

		return node;
	}
	
	/**
	 * Returns the algorithm's name together with current options.
	 * 
	 * @return algorithm's name and configurable options
	 */
	@Override
	public String toString ()
	{
		String m;
		
		switch (this.mode)
		{
			case OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN:
				m = "-BL-";
				break;
				
			case OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN:
				m = "-CI-";
				break;
			
			default:
				m = "-?";
				break;
		}
		
		return this.getClass().getSimpleName() + m + window_size + "-" + kvalue;
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
