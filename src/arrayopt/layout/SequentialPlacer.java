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
 *
 */
public class SequentialPlacer implements LayoutAlgorithm, FillingAlgorithm
{
	/**
	 * Constant to indicate that no ordering of probes should be performed
	 * before placement/filling.
	 */
	public static final int KEEP_ORDER = 0;
	
	/**
	 * Constant to indicate that the order of the probes should be randomized
	 * before placement/filling.
	 */
	public static final int RANDOM = 1;
	
	/**
	 * Constant to indicate that probes should be ordered lexicographically by
	 * their sequences before placement/filling.
	 */
	public static final int SORT_SEQUENCES = 2;
	
	/**
	 * Constant to indicate that probes should be ordered lexicographically by
	 * their binary embeddings before placement/filling.
	 */
	public static final int SORT_EMBEDDINGS = 3;
	
	/**
	 * Constant to indicate that a TSP-tour of the probes must be computed
	 * before placement/filling, in order to minimize border conflicts between
	 * neighboring probes. The TSP-tour is computed on a graph where nodes
	 * represent the probes, and edges between two probes contain the number of
	 * border conflicts between their embeddings.
	 */
	public static final int TSP_ORDER = 4;
	
	/**
	 * Indicate the ordering of the probes used for placement/filling.
	 */
	private int order;
	
	/**
	 * TODO document this
	 */
	public SequentialPlacer ()
	{
		this(KEEP_ORDER);
	}

	/**
	 * TODO document this
	 */
	public SequentialPlacer (int order)
	{
		switch (order)
		{
			case KEEP_ORDER:
			case RANDOM:
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
	 *
	 */
	public void changeLayout (Chip chip)
	{
		int		id[];

		// reset current layout (if any)
		chip.resetLayout();

		// get list of movable probes
		id = chip.getMovableProbes ();

		fillRegion (chip, chip.getChipRegion(), id, 0, id.length - 1);
		
		return;
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
		ProbeOrderingAlgorithm ordering = null;
		RectangularRegion r;

		if (!(region instanceof RectangularRegion))
			throw new IllegalArgumentException
				("Only rectangular regions are supported.");

		r = (RectangularRegion) region;
		
		switch(this.order)
		{
			case RANDOM:
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

	/**
	 *
	 */
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
				if (chip.isFixedSpot(r, c)) continue;

				chip.spot[r][c] = probe_id[start];

				start ++;

				if (start > end)
					// all probes were placed
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
	
	/**
	 * Returns the algorithm's name together with current options.
	 */
	@Override
	public String toString ()
	{
		String ordering;
		
		switch (this.order)
		{
			case RANDOM:
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
