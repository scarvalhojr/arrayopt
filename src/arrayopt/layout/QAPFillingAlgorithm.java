/*
 * QAPFillingAlgorithm.java
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

import arrayopt.qap.*;

/**
 *
 */
public class QAPFillingAlgorithm implements PlacementAlgorithm, FillingAlgorithm
{
	protected int[] spot_dist;

	protected int[] probe_dist;

	protected int[] perm;

	protected QAPSolverAlgorithm solver;

	public QAPFillingAlgorithm (QAPSolverAlgorithm solver)
	{
		this.solver = solver;
	}

	/**
	 * TEMPORARY: FOR TESTING ONLY
	 *
	 * REMOVE THIS AND THE 'implements PlacementAlgorithm' above!
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
	public int fillRegion (Chip chip, Region region, int probe_id[], int start,
		int end)
	{
		RectangularRegion	rect;
		int					dim, num_probes, unplaced;
		long				curr_cost, sol_cost, comp_cost;

		if (!(region instanceof RectangularRegion))
			throw new IllegalArgumentException
				("Only rectangular regions are supported.");

		rect = (RectangularRegion) region;

		// create arrays
		dim = createArrays (chip, rect);

		// compute spot distance matrix
		if (chip instanceof SimpleChip)
			QAPHelper.computeSpotDistance ((SimpleChip) chip, dim, rect,
				spot_dist);
		else if (chip instanceof AffymetrixChip)
			QAPHelper.computeSpotDistance ((AffymetrixChip) chip, dim, rect,
				spot_dist);
		else
			throw new IllegalArgumentException
				("Unsupported chip type.");

		// check if number of probes fit in the region
		num_probes = end - start + 1;

		if (num_probes > dim)
		{
			// exceeding probes have to be ignored
			unplaced = num_probes - dim;
			end = start + dim - 1;
		}
		else
		{
			// all probes can be placed
			unplaced = 0;
		}

		// compute probe distance matrix
		if (chip instanceof SimpleChip)
			QAPHelper.computeProbeDistance ((SimpleChip) chip, dim, probe_id,
											start, end, probe_dist);
		else
			QAPHelper.computeProbeDistance ((AffymetrixChip) chip, dim, probe_id,
											start, end, probe_dist);
		// REMOVE THIS!
		curr_cost = 0;
		for (int i = 0; i < dim; i++)
			for (int j = 0; j < dim; j++)
				curr_cost += spot_dist[i * dim + j] * probe_dist[i * dim + j];

		// solve QAP
		sol_cost = solver.solve (dim, spot_dist, probe_dist, perm);

		// place probes according to the optimal permutation
		if (chip instanceof SimpleChip)
			QAPHelper.applyPermutation ((SimpleChip) chip, rect, dim, probe_id,
										start, end, perm);
		else
			QAPHelper.applyPermutation ((AffymetrixChip) chip, rect, dim, probe_id,
										start, end, perm);

		return unplaced;
	}

	protected int createArrays (Chip chip, RectangularRegion r)
	{
		int		dim, num_rows, num_cols, rows_per_probe;

		num_rows = r.last_row - r.first_row + 1;
		num_cols = r.last_col - r.first_col + 1;

		if (chip instanceof AffymetrixChip)
			rows_per_probe = 2;
		else
			rows_per_probe = 1;

		dim = ((int) (num_rows / rows_per_probe)) * num_cols;

		// create spot distance matrix (if necessary)
		if (spot_dist == null)
		{
			spot_dist = new int [dim * dim];
		}
		else
		{
			if (spot_dist.length != dim * dim)
			{
				spot_dist = new int [dim * dim];
			}
		}

		// create probe distance matrix (if necessary)
		if (probe_dist == null)
		{
			probe_dist = new int [dim * dim];
		}
		else
		{
			if (probe_dist.length != dim * dim)
			{
				probe_dist = new int [dim * dim];
			}
		}

		// create permutation array (if necessary)
		if (perm == null)
		{
			perm = new int [dim];
		}
		else
		{
			if (perm.length != dim)
			{
				perm = new int [dim];
			}
		}

		return dim;
	}
}
