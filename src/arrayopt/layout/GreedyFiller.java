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

/**
 *
 */
public class GreedyFiller implements PlacementAlgorithm, FillingAlgorithm
{
	private int window_size;

	private boolean randomize_first;

	public static boolean DEFAULT_RANDOMIZE_FIRST = false;

	public GreedyFiller ()
	{
		this(0, DEFAULT_RANDOMIZE_FIRST);
	}

	public GreedyFiller (int window_size)
	{
		this(window_size, DEFAULT_RANDOMIZE_FIRST);
	}

	public GreedyFiller (int window_size, boolean randomize_first)
	{
		this.window_size = window_size;
		this.randomize_first = randomize_first;
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
		// to do

		return 0;
	}

	/**
	 *
	 */
	protected int fillRegion (AffymetrixChip chip, RectangularRegion region,
		int probe_id[], int start, int end)
	{
		int rand, tmp;

		if (randomize_first)
		{
			// randomize the list first
			for (int i = start; i < end; i++)
			{
				// select a random element of the list
				rand = start + (int) ((end - start + 1) * Math.random());

				// swap the selected element with
				// the one at the current position
				tmp = probe_id[rand];
				probe_id[i] = probe_id[rand];
				probe_id[rand] = tmp;
			}
		}

		for (int r = region.first_row; r < region.last_row; r ++)
		{
			for (int c = region.first_col; c <= region.last_col; c++)
			{
				// skip spot pair if not empty
				if (chip.spot[r][c] != chip.EMPTY_SPOT ||
					chip.spot[r+1][c] != chip.EMPTY_SPOT)
					continue;

				// skip fixed spot pair
				if (chip.isFixedSpot(r, c) || chip.isFixedSpot(r+1, c))
					continue;

				// select probe pair which minimizes conflict
				// (and move it to position 'start' in the list)
				selectNextPair (chip, r, c, probe_id, start, end);

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

	protected void selectNextPair (AffymetrixChip chip, int row, int col,
		int probe_id[], int start, int end)
	{
		double	conflict, min_conflict;
		int		best_probe, tmp;

		// if window size is limited
		if (window_size > 0)
		{
			// check if probe list is larger than window
			if (end - start + 1 > window_size)
				// yes: limit search space
				end = start + window_size - 1;
		}

		min_conflict = LayoutEvaluation.spotConflict(chip, row, col,
						probe_id[start]);

		best_probe = start;

		for (int i = start + 1; i <= end; i++)
		{
			conflict = LayoutEvaluation.spotConflict(chip, row, col,
						probe_id[i]);

			if (conflict < min_conflict)
			{
				min_conflict = conflict;
				best_probe = i;
			}
		}

		// move best probe ID to the front of the list
		tmp = probe_id[start];
		probe_id[start] = probe_id[best_probe];
		probe_id[best_probe] = tmp;
	}
}
