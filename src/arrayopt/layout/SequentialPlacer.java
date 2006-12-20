/*
 * SequentialPlacer.java
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
 * Class that implements the Sequential placement algorithm. The algorithm fills
 * spots sequentially, row-by-row, from top to bottom, left-to-right.
 * 
 * <P>The order of probes placed at each spot can be configured at instantiation
 * time with the following constants: {@link #KEEP_ORDER},
 * {@link #RANDOM_ORDER}, {@link #SORT_EMBEDDINGS}, {@link #SORT_SEQUENCES} and
 * {@link #TSP_ORDER}.</P> 
 * 
 * @author Sergio A. de Carvalho Jr.
 */
public class SequentialPlacer implements LayoutAlgorithm, FillingAlgorithm
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
	 * Creates an instance of the SequentialPlacer algorithm with the default
	 * probe ordering set to {@link #KEEP_ORDER}.
	 */
	public SequentialPlacer ()
	{
		this(KEEP_ORDER);
	}

	/**
	 * Creates an instance of the SequentialPlacer algorithm with the specified
	 * probe ordering.
	 * 
	 * @param order order of probes during placement
	 * @see #KEEP_ORDER
	 * @see #RANDOM_ORDER
	 * @see #SORT_EMBEDDINGS
	 * @see #SORT_SEQUENCES
	 * @see #TSP_ORDER
	 */
	public SequentialPlacer (int order)
	{
		switch (order)
		{
			case KEEP_ORDER:
			case RANDOM_ORDER:
			case SORT_SEQUENCES:
			case SORT_EMBEDDINGS:
			case TSP_ORDER:
				this.order = order;
				break;

			default:
				throw new IllegalArgumentException ("Unknown probe ordering.");
		}
	}
	
	/**
	 * Uses the Sequential placer algorithm to re-place the probes. All probes
	 * are initially removed from the spots, and a list with the desired order
	 * (see {@link #order}) is passed to the
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
	 * Fills the spots sequentially with a given list a probes.
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
	 * Fills the spots sequentially with a given list a probes. The elements of
	 * the list as bounded by the given parameters. 
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
		ProbeOrderingAlgorithm ordering = null;
		RectangularRegion r;

		if (!(region instanceof RectangularRegion))
			throw new IllegalArgumentException
				("Only rectangular regions are supported.");

		r = (RectangularRegion) region;
		
		switch(this.order)
		{
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
		}
		
		if (ordering != null)
			ordering.orderProbes(chip, probe_id, start, end);
		
		if (chip instanceof SimpleChip)
			return fillRegion ((SimpleChip) chip, r, probe_id, start, end);

		else if (chip instanceof AffymetrixChip)
			return fillRegion ((AffymetrixChip) chip, r, probe_id, start, end);

		else
			throw new IllegalArgumentException ("Unsupported chip type.");
	}

	protected int fillRegion (SimpleChip chip, RectangularRegion region,
		int probe_id[], int start, int end)
	{
		for (int r = region.first_row; r <= region.last_row; r ++)
		{
			for (int c = region.first_col; c <= region.last_col; c++)
			{
				// skip spot if not empty
				if (chip.spot[r][c] != Chip.EMPTY_SPOT)
					continue;
				
				// skip fixed spots
				if (chip.isFixedSpot(r, c))
					continue;
				
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

	protected int fillRegion (AffymetrixChip chip, RectangularRegion region,
		int probe_id[], int start, int end)
	{
		for (int r = region.first_row; r < region.last_row; r ++)
		{
			for (int c = region.first_col; c <= region.last_col; c++)
			{
				// skip spot pair if not empty
				if (chip.spot[r][c] != Chip.EMPTY_SPOT || chip.spot[r+1][c] != Chip.EMPTY_SPOT)
					continue;

				// skip fixed spot pairs
				if (chip.isFixedSpot(r, c) || chip.isFixedSpot(r+1, c))
					continue;
				
				// place next element of the list
				chip.spot[r][c] = probe_id[start];
				chip.spot[r+1][c] = probe_id[start] + 1;

				// advance list pointer
				if (++start > end)
					// all probes were placed
					return 0;
			}
		}

		// some probe pairs could not be placed
		return 2 * (end - start + 1);
	}
	
	/**
	 * Returns the algorithm's name together with current options.
	 * 
	 * @return algorithm's name and configurable options
	 */
	@Override
	public String toString ()
	{
		String ordering;
		
		switch (this.order)
		{
			case RANDOM_ORDER:
				ordering = "-Random";
				break;
				
			case SORT_SEQUENCES:
				ordering = "-SortSequences";
				break;
				
			case SORT_EMBEDDINGS:
				ordering = "-SortEmbeddings";
				break;
				
			case TSP_ORDER: 
				ordering = "-TSP";
				break;
			
			default:
				ordering = "-KeepOrder";
				break;
		}
		
		return this.getClass().getSimpleName() + ordering;
	}
}
