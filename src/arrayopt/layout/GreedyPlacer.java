/*
 * GreedyPlacer.java
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
 * This class implements the Greedy placement algorithm.
 * 
 * <P>The Greedy algorithm creates a new layout by filling spots sequentially
 * using a k-threading pattern identical to the one implemented in
 * {@link KThreadingPlacer}. The amplitude of the k-threading is set with the
 * {@link #kvalue} variable (at instantiation time). For every spot
 * <CODE>s</CODE>, it finds a probe <CODE>p</CODE> with minimum cost to fill
 * <CODE>s</CODE>. The cost is either the sum of border conflicts
 * ({@link #BORDER_LENGTH_MIN}) or sum of conflict indices
 * ({@link #CONFLICT_INDEX_MIN}) with filled neighbors.</P>
 * 
 * <P>The algorithm looks at {@link #window_size} probe candidate for filling
 * each spot. All probes are initially ordered according to the specified order
 * and a doubly-linked list is created, where each node of the list contains a
 * single probe ID. Once a probe is placed, it is removed from the list, and the
 * next search of probe candidates starts at the next or previous element in the
 * list. The search will examine at most {@link #window_size} elements,
 * preferably <CODE>window_size/2</CODE> to the left and
 * <CODE>window_size/2</CODE> to the right of the last placed probe. This
 * effectively improves the chances of finding probes similar to its neighbors
 * to fill a given spot.</P>
 * 
 * The doubly-linked list is used to maintain the desired order, which must be
 * specified at instantiation time with the following constants:
 * {@link #KEEP_ORDER}, {@link #RANDOM_ORDER}, {@link #SORT_EMBEDDINGS},
 * {@link #SORT_SEQUENCES} and {@link #TSP_ORDER}.</P>
 * 
 * @author Sergio A. de Carvalho Jr.
 */
public class GreedyPlacer implements LayoutAlgorithm, FillingAlgorithm
{
	/**
	 * Constant used to indicate that the algorithm should try to minimize the
	 * total border length of the microarray layout.
	 */
	public static final int BORDER_LENGTH_MIN = 0;
	
	/**
	 * Constant used to indicate that the algorithm should try to minimize the
	 * sum of conflict indices of the microarray layout.
	 */
	public static final int CONFLICT_INDEX_MIN = 1;

	/**
	 * Constant to indicate that no ordering of probes should be performed
	 * before placement/filling. This actually means that spots are filled with
	 * probes ordered as they were received in the list.
	 */
	public static final int KEEP_ORDER = 0;
	
	/**
	 * Constant to indicate that the order of the probes should be randomized
	 * before placement/filling. The randomization is performed by the
	 * {@link RandomOrdering} class.
	 */
	public static final int RANDOM_ORDER = 1;
	
	/**
	 * Constant to indicate that probes should be ordered lexicographically by
	 * their sequences before placement/filling. The ordering is performed by
	 * the {@link SortedSequencesOrdering} class.
	 */
	public static final int SORT_SEQUENCES = 2;
	
	/**
	 * Constant to indicate that probes should be ordered lexicographically by
	 * their binary embeddings before placement/filling. The ordering is
	 * performed by the {@link SortedEmbeddingsOrdering} class.
	 */
	public static final int SORT_EMBEDDINGS = 3;
	
	/**
	 * Constant to indicate that a TSP-tour of the probes must be computed
	 * before placement/filling, in order to minimize border conflicts between
	 * neighboring probes. The TSP-tour is computed on a graph where nodes
	 * represent the probes, and edges between two probes contain the number of
	 * border conflicts between their embeddings. The ordering is performed by
	 * the {@link TSPOrdering} class.
	 */
	public static final int TSP_ORDER = 4;
	
	/**
	 * This variable stores the current minimization mode used by the algorithm.
	 * Possible values are {@link #BORDER_LENGTH_MIN} and
	 * {@link #CONFLICT_INDEX_MIN}. 
	 */
	private int mode;

	/**
	 * Indicate the ordering of the probes used for placement/filling.
	 */
	private int order;
	
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
	
	private ProbeOrderingAlgorithm ordering;
	
	private RectangularRegion chip_region;
	
	private int embed_len;
	
	private int probe_len;
	
	private int ci_dim;
	
	private double pos_weight[];
	
	private double m_cost[];
	
	private double u_cost[];
	
	/**
	 * Creates an instance of the Greedy placement algorithm with the desired
	 * minimization mode, probe ordering and window size.
	 * 
	 * @see #BORDER_LENGTH_MIN
	 * @see #CONFLICT_INDEX_MIN
	 * @see #KEEP_ORDER
	 * @see #RANDOM_ORDER
	 * @see #SORT_EMBEDDINGS
	 * @see #SORT_SEQUENCES
	 * @see #TSP_ORDER
	 * @see #window_size
	 * @param mode minimization mode, either border length or conflict index
	 * @param order order of probes during placement/filling
	 * @param window_size maximum number of candidades examined for each spot
	 */
	public GreedyPlacer (int mode, int window_size, int kvalue, int order)
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
		
		switch(order)
		{
			case KEEP_ORDER:
				ordering = null;
				break;
				
			case RANDOM_ORDER:
				ordering = new RandomOrdering();
				break;
				
			case SORT_SEQUENCES:
				ordering = new SortedSequencesOrdering();
				break;
				
			case SORT_EMBEDDINGS:
				ordering = new SortedEmbeddingsOrdering();
				break;
				
			case TSP_ORDER: 
				ordering = new TSPOrdering();
				break;
				
			default:
				throw new IllegalArgumentException ("Unknown probe ordering.");
		}
		
		if (kvalue < 0)
			throw new IllegalArgumentException ("Invalid k value: " + kvalue);
		
		this.window_size = window_size;
		this.kvalue = kvalue;
		this.order = order;
	}

	/**
	 * Creates a new layout of a microarray chip using the Greedy placement
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
	 * a starting and ending position.
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
		
		// saves a reference to the full chip region
		this.chip_region = chip.getChipRegion();

		// apply probe ordering
		if (ordering != null)
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

		if (mode == CONFLICT_INDEX_MIN)
			// prepare for local conflict index calculations
			conflictIndexSetup (chip);
		
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
				
				// skip spot if not empty
				if (chip.spot[r][c] != Chip.EMPTY_SPOT) continue;
				
				// skip fixed spot
				if (chip.isFixedSpot(r, c)) continue;

				// find probe with minimum cost
				if (mode == BORDER_LENGTH_MIN)
					list = minBorderLength (chip, r, c, list);
				else // (mode == CONFLICT_INDEX_MIN)
					list = minConflictIndex (chip, r, c, list);
				
				// place selected probe
				list = fillSpot (chip, r, c, list);
				
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
	
	private MyLinkedList fillSpot (SimpleChip chip, int r, int c,
			MyLinkedList list)
	{
		MyLinkedList n = null;
		
		// place selected probe
		chip.spot[r][c] = list.info;
		
		// delete element from the list
		if (list.next != null)
		{
			list.next.prev = list.prev;
			n = list.next;
		}
		if (list.prev != null)
		{
			list.prev.next = list.next;
			n = list.prev;
		}

		return n;
	}

	private MyLinkedList minBorderLength (SimpleChip chip, int row, int col,
			MyLinkedList node)
	{
		MyLinkedList best;
		long	cost, min;
		int		count;
		
		// find node so that the search will examine 
		// window_size elements around the last placed probe 
		node = findStartingNode (node);

		// compute cost of first element
		min = LayoutEvaluation.borderLength (chip, row, col, node.info);
		if (min == 0) return node;
		
		best = node;
		node = node.next;
		
		for (count = 1; node != null && count < window_size; count++)
		{
			cost = LayoutEvaluation.borderLength (chip, row, col, node.info);
			if (cost < min)
			{
				min = cost;
				best = node;
			}

			node = node.next;
		}
		
		return best;
	}

	private MyLinkedList minConflictIndex (SimpleChip chip, int row, int col,
			MyLinkedList node)
	{
		MyLinkedList best;
		double	cost, min;
		int		count;
		boolean	empty;
		
		// find node so that the search will examine 
		// window_size elements around the last placed probe 
		node = findStartingNode (node);
		
		// prepare conflict index costs
		empty = examineNeighbors (chip, row, col);
		
		// if region around the spot is empty,
		// place any probe and return
		if (empty) return node;

		// compute cost of first element
		min = conflictIndex (chip, node.info, Double.POSITIVE_INFINITY);
		if (min == 0) return node;
		
		best = node;
		node = node.next;
		
		for (count = 1; node != null && count < window_size; count++)
		{
			cost = conflictIndex(chip, node.info, min);
			if (cost < min)
			{
				min = cost;
				best = node;
			}

			node = node.next;
		}
		
		return best;
	}
	
	private void conflictIndexSetup (SimpleChip chip)
	{
		// prepare for faster conflict index calculations
		ci_dim = ConflictIndex.dimConflictRegion();
		embed_len = chip.getEmbeddingLength();
		probe_len = chip.getProbeLength();
		
		// instantiate local arrays (if necessary)
		if (m_cost == null)
		{
			m_cost = new double [embed_len];
			u_cost = new double [embed_len];
			pos_weight = new double [probe_len + 1];
		}
		else if (m_cost.length != embed_len)
		{
			m_cost = new double [embed_len];
			u_cost = new double [embed_len];
			pos_weight = new double [probe_len + 1];
		}
		
		// store position weights locally
		for (int b = 0; b <= probe_len; b++)
			pos_weight[b] = ConflictIndex.positionWeight(b, probe_len);
	}
	
	private boolean examineNeighbors (SimpleChip chip, int row, int col)
	{
		boolean empty = true;
		double delta;
		int r, c, id, base, step, word, bitmask = 0;
		
		// reset costs
		for (step = 0; step < embed_len; step++)
			m_cost[step] = u_cost[step] = 0;
		
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
				
				for (base = 0, step = 0, word = - 1; step < embed_len; step++)
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
						m_cost[step] += delta;
						
						base++;
					}
					else
					{
						// spot is in a masked step
						u_cost[step] += pos_weight[base] * delta;
					}
				}
			}
		}

		return empty;
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
				ci += u_cost[step];

				// increment the number of synthesized bases
				base++;
			}
			else
			{
				// spot is in a masked step
				ci += pos_weight[base] * m_cost[step];
				
				// stop if CI exceeds limit
				if (ci > max) break;
			}
		}

		return ci;
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
		String ord;
		
		switch (this.mode)
		{
			case BORDER_LENGTH_MIN:
				m = "-BL-";
				break;
				
			case CONFLICT_INDEX_MIN:
				m = "-CI-";
				break;
			
			default:
				m = "-?";
				break;
		}
		
		switch (this.order)
		{
			case RANDOM_ORDER:
				ord = "-Random";
				break;
				
			case SORT_SEQUENCES:
				ord = "-SortSequences";
				break;
				
			case SORT_EMBEDDINGS:
				ord = "-SortEmbeddings";
				break;
				
			case TSP_ORDER: 
				ord = "-TSP";
				break;
			
			default:
				ord = "-KeepOrder";
				break;
		}

		return this.getClass().getSimpleName() + m + window_size + ord;
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
