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

/**
 * This class implements the k-threading placement algorithm. The algorithm is
 * described in the paper:</P>
 * 
 * <P>Hannenhalli, S.; Hubell, E.; Lipshutz, R. & Pevzner, P.A.: Combinatorial
 * algorithms for design of DNA arrays, Adv Biochem Eng Biotechnol, 2002, 77,
 * 1-19.</P>
 * 
 * <P>A k-threading is a variation of the standard row-by-row threading (see
 * {@link SequentialPlacer}, which fills the spots in an alternating
 * right-to-left and left-to-right paths interspaced with upward and downward
 * movements over k sites. The parameter k (@link #kvalue} is called the
 * amplitude of the threading.</P>
 * 
 * <P>The order of probes placed at each spot can be configured at instantiation
 * time with the following constants: {@link #KEEP_ORDER},
 * {@link #RANDOM_ORDER}, {@link #SORT_EMBEDDINGS}, {@link #SORT_SEQUENCES} and
 * {@link #TSP_ORDER}.</P> 
 * 
 * @author Sergio A. de Carvalho Jr.
 */
public class KThreadingPlacer implements LayoutAlgorithm, FillingAlgorithm
{
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
	 * Indicate the ordering of the probes used for placement/filling.
	 */
	private int order;
	
	/**
	 * The amplitude of the k-threading. This gives the number of upward and
	 * downward movements between the right-to-left and left-to-right paths over
	 * the chip.
	 */
	private int kvalue;
	
	private ProbeOrderingAlgorithm ordering;
	
	/**
	 * Creates an instance of the KThreadingPlacer algorithm with the given
	 * value of k and the default probe ordering set to {@link #KEEP_ORDER}.
	 */
	public KThreadingPlacer (int kvalue)
	{
		this(kvalue, KEEP_ORDER);
	}

	/**
	 * Creates an instance of the KThreadingPlacer algorithm with the given
	 * value of k and probe ordering.
	 * 
	 * @param kvalue amplitude of k-threading
	 * @param order order of probes during placement/filling
	 * @see #KEEP_ORDER
	 * @see #RANDOM_ORDER
	 * @see #SORT_EMBEDDINGS
	 * @see #SORT_SEQUENCES
	 * @see #TSP_ORDER
	 */
	public KThreadingPlacer (int kvalue, int order)
	{
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
		
		this.order = order;
		this.kvalue = kvalue;
	}
	
	/**
	 * Uses the k-threading placement algorithm to re-place the probes. All
	 * probes are initially removed from the spots, and a list with the desired
	 * order (see {@link #order}) is passed to the 
	 * {@link #fillRegion(Chip, Region, int[]) method.
	 * 
	 * @param chip chip instance
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
	 * Fills the spots with a given list a probes.
	 * 
	 * @param chip chip instance
	 * @param region region of the chip to be filled
	 * @param probe_id list of probe IDs 
	 */
	public int fillRegion (Chip chip, Region region, int probe_id[])
	{
		return fillRegion (chip, region, probe_id, 0, probe_id.length - 1);
	}

	/**
	 * Fills the spots with a given list a probes. The elements of the list as
	 * bounded by the given parameters. 
	 * 
	 * @param chip chip instance
	 * @param region region of the chip to be filled
	 * @param probe_id list of probe IDs 
	 * @param start first element of the list
	 * @param end last element of the list
	 */
	public int fillRegion (Chip chip, Region region, int probe_id[], int start,
		int end)
	{
		RectangularRegion r;

		if (!(region instanceof RectangularRegion))
			throw new IllegalArgumentException
				("Only rectangular regions are supported.");

		r = (RectangularRegion) region;

		// apply probe ordering
		if (ordering != null)
			ordering.orderProbes(chip, probe_id, start, end);

		if (chip instanceof SimpleChip)
			return fillRegion ((SimpleChip) chip, r, probe_id, start, end);
		// else
			throw new IllegalArgumentException ("Unsupported chip type.");
	}

	private int fillRegion (SimpleChip chip, RectangularRegion region,
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
				
				// skip spot if not empty
				if (chip.spot[r][c] != Chip.EMPTY_SPOT) continue;
				
				// skip fixed spot
				if (chip.isFixedSpot(r, c)) continue;
				
				// place next element of the list
				chip.spot[r][c] = probe_id[start];

				// advance list pointer
				if (++start > end)
					// all probes were placed
					return 0;
			}
		}

		// some probes could not be placed
		return end - start + 1;
	}

	/**
	 * Returns the algorithm's name together with current options.
	 * 
	 * @return algorithm's name and configurable options
	 */
	@Override
	public String toString ()
	{
		String ord;
		
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
		
		return kvalue + "-threading" + ord;
	}
}
