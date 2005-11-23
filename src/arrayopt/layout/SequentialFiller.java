/*
 * SequentialFiller.java
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
public class SequentialFiller extends FillingAlgorithm
{
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
		if (!(region instanceof RectangularRegion))
			throw new IllegalArgumentException ("Only rectangular regions are supported.");

		if (chip instanceof SimpleChip)
			return fillRegion ((SimpleChip) chip, (RectangularRegion) region, first_probe, last_probe);

		else if (chip instanceof AffymetrixChip)
			return fillRegion ((AffymetrixChip) chip, (RectangularRegion) region, first_probe, last_probe);

		else
			throw new IllegalArgumentException ("Unsupported chip type.");
	}

	protected int fillRegion (SimpleChip chip, RectangularRegion region, int first_probe, int last_probe)
	{
		for (int r = region.first_row; r <= region.last_row; r ++)
		{
			for (int c = region.first_col; c <= region.last_col; c++)
			{
				// skip spot if not empty
				if (chip.spot[r][c] != chip.EMPTY_SPOT)
					continue;

				// skip fixed spots
				if (chip.isFixedSpot(r, c)) continue;

				chip.spot[r][c] = chip.probe_list[first_probe];

				first_probe ++;

				if (first_probe > last_probe)
					// all probes were placed
					return 0;
			}
		}

		// some probes could not be placed
		return last_probe - first_probe + 1;
	}

	protected int fillRegion (AffymetrixChip chip, RectangularRegion region, int first_probe, int last_probe)
	{
		for (int r = region.first_row; r < region.last_row; r ++)
		{
			for (int c = region.first_col; c <= region.last_col; c++)
			{
				// skip spot pair if not empty
				if (chip.spot[r][c] != chip.EMPTY_SPOT || chip.spot[r+1][c] != chip.EMPTY_SPOT)
					continue;

				// skip fixed spot pairss
				if (chip.isFixedSpot(r, c) || chip.isFixedSpot(r+1, c))
					continue;

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
}
