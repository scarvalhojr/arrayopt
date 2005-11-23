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
public class GreedyFiller extends FillingAlgorithm
{
	private int window_size;

	private boolean randomize_first;

	public static boolean DEFAULT_RANDOMIZE_FIRST = false;

	private boolean DEBUG = false;

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
	 * Note that this method has package access only. To use it, classes should
	 * call the {@link Chip#placeProbes placeProbes} method on the
	 * {@linkplain Chip} class passing an instance of this class.
	 */
	int placeProbes (Chip chip, int first_probe, int last_probe)
	{
		return fillRegion (chip, chip.getChipRegion(), first_probe, last_probe);
	}

	/**
	 * Note that this method has package access only. It should only be used
	 * internally or by an instance of another {@linkplain PlacementAlgorithm}.
	 */
	int fillRegion (Chip chip, Region region, int first_probe, int last_probe)
	{
		if (!(chip instanceof AffymetrixChip))
			throw new IllegalArgumentException ("Unsupported chip type.");

		if (!(region instanceof RectangularRegion))
			throw new IllegalArgumentException ("Only rectangular regions are supported.");

		return fillRegion ((AffymetrixChip) chip, (RectangularRegion) region, first_probe, last_probe);
	}

	protected int fillRegion (AffymetrixChip chip, RectangularRegion region, int first_probe, int last_probe)
	{
		int rand, tmp;

		if (DEBUG)
		{
			System.out.println ("Filling region " + region);
			System.out.println ("\tProbe IDs: " + first_probe + " - " + last_probe);
		}

		if (randomize_first)
		{
			// randomize the list first
			for (int i = first_probe; i < last_probe; i++)
			{
				// select a random element of the list
				rand = first_probe + (int) ((last_probe - first_probe + 1) * Math.random());

				// swap the selected element with
				// the one at the current position
				tmp = chip.probe_list[rand];
				chip.probe_list[i] = chip.probe_list[rand];
				chip.probe_list[rand] = tmp;
			}
		}

		for (int r = region.first_row; r < region.last_row; r ++)
		{
			for (int c = region.first_col; c <= region.last_col; c++)
			{
				// skip spot pair if not empty
				if (chip.spot[r][c] != chip.EMPTY_SPOT || chip.spot[r+1][c] != chip.EMPTY_SPOT)
					continue;

				// skip fixed spot pair
				if (chip.isFixedSpot(r, c) || chip.isFixedSpot(r+1, c))
					continue;

				// select probe pair which minimizes conflict
				// (and move it to position 'first_probe' in the list)
				selectNextPair (chip, r, c, first_probe, last_probe);

				chip.spot[r][c] = chip.probe_list[first_probe];
				chip.spot[r+1][c] = chip.probe_list[first_probe] + 1;

				first_probe ++;

				if (first_probe > last_probe)
					// all probes were placed
					return 0;
			}
		}

		// some probes could not be placed
		return 2 * (last_probe - first_probe + 1);
	}

	protected void selectNextPair (AffymetrixChip chip, int row, int col, int first_probe, int last_probe)
	{
		double	conflict, min_conflict;
		int		best_probe, tmp;

		// if window size is limited
		if (window_size > 0)
		{
			// check if probe list is larger than window
			if (last_probe - first_probe + 1 > window_size)
				// yes: limit search space
				last_probe = first_probe + window_size - 1;
		}

		min_conflict = LayoutEvaluation.spotConflict(chip, row, col, chip.probe_list[first_probe]);
		best_probe = first_probe;

		for (int i = first_probe + 1; i <= last_probe; i++)
		{
			conflict = LayoutEvaluation.spotConflict(chip, row, col, chip.probe_list[i]);
			if (conflict < min_conflict)
			{
				min_conflict = conflict;
				best_probe = i;
			}
		}

		// move best probe ID to the front of the list
		tmp = chip.probe_list[first_probe];
		chip.probe_list[first_probe] = chip.probe_list[best_probe];
		chip.probe_list[best_probe] = tmp;
	}
}
