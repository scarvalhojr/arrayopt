/*
 * RandomFiller.java
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
public class RandomFiller implements PlacementAlgorithm, FillingAlgorithm
{
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
		int rand, tmp;

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

				// select probe randomly
				// (and move it to position 'start' in the list)
				rand = start + (int) ((end - start + 1) * Math.random());

				// move probe ID to the front of the list
				tmp = probe_id[start];
				probe_id[start] = probe_id[rand];
				probe_id[rand] = tmp;

				chip.spot[r][c] = probe_id[start];

				start ++;

				if (start > end)
					// all probes were placed
					return 0;
			}
		}

		// some probes could not be placed
		return 2 * (end - start + 1);
	}

	/**
	 *
	 */
	protected int fillRegion (AffymetrixChip chip, RectangularRegion region,
		int probe_id[], int start, int end)
	{
		int rand, tmp;

		for (int r = region.first_row; r < region.last_row; r ++)
		{
			for (int c = region.first_col; c <= region.last_col; c++)
			{
				// skip spot pair if not empty
				if (chip.spot[r][c] != Chip.EMPTY_SPOT || chip.spot[r+1][c] != Chip.EMPTY_SPOT)
					continue;

				// skip fixed spot pair
				if (chip.isFixedSpot(r, c) || chip.isFixedSpot(r+1, c))
					continue;

				// select probe pair randomly
				// (and move it to position 'start' in the list)
				rand = start + (int) ((end - start + 1) * Math.random());

				// move probe ID to the front of the list
				tmp = probe_id[start];
				probe_id[start] = probe_id[rand];
				probe_id[rand] = tmp;

				chip.spot[r][c] = probe_id[start];
				chip.spot[r+1][c] = probe_id[start] + 1;

				start ++;

				if (start > end)
					// all probes were placed
					return 0;
			}
		}

		// some probes could not be placed
		return 2 * (end - start + 1);
	}
}
