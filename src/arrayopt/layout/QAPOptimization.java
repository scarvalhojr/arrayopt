/*
 * QAPOptimization.java
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
 * in case of concurrent access callers should synchronize on this object
 */
public class QAPOptimization implements IteractiveOptimizationAlgorithm,
	FillingAlgorithm, PlacementAlgorithm
{
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
	 * document this
	 */
	protected QAPSolverAlgorithm solver;

	/**
	 * document this
	 */
	protected int[] probe_id;

	/**
	 * document this
	 */
	protected int[] spot_dist;

	/**
	 * document this
	 */
	protected int[] probe_dist;

	/**
	 * document this
	 */
	protected int[] perm;

	/**
	 * document this
	 *
	 * dimension of current/last instance
	 */
	protected int dim;

	/**
	 * document this
	 *
	 * height of last rectangular region
	 */
	protected int last_height;

	/**
	 * document this
	 *
	 * width of last rectangular region
	 */
	protected int last_width;

	/**
	 * document this
	 */
	public QAPOptimization (QAPSolverAlgorithm solver)
	{
		this.solver = solver;
		this.dim = -1;
		this.last_height = -1;
		this.last_width = -1;
	}

	/**
	 * document this
	 */
	public float optimizeLayout (Chip chip, Region region)
	{
		// to do

		return 0;
	}

	/**
	 * document this
	 */
	public int fillRegion (Chip chip, Region region, int probe_id[])
	{
		return fillRegion (chip, region, probe_id, 0, probe_id.length - 1);
	}

	/**
	 * document this
	 */
	public int fillRegion (Chip chip, Region region, int probe_id[], int start,
		int end)
	{
		RectangularRegion	rect;
		int					rows_per_probe, dist_empty, num_probes, unplaced;
		long				curr_cost, sol_cost, comp_cost;

		if (!(region instanceof RectangularRegion))
			throw new IllegalArgumentException
				("Only rectangular regions are supported.");

		// chip type configuration
		if (chip instanceof SimpleChip)
		{
			// simple chips use 1 row per probe
			rows_per_probe = 1;

			// distance from string to any probe
			// equals probe length
			dist_empty = chip.getProbeLength();
		}
		else if (chip instanceof AffymetrixChip)
		{
			// affymetrix chips use 2 rows per probe (pair)
			rows_per_probe = 2;

			// distance from string to any probe
			// equals probe length + 1 since probes appear
			// in pairs which differ in the middle base
			dist_empty = chip.getProbeLength() + 1;
		}
		else
			throw new IllegalArgumentException
				("Unsupported chip type.");

		rect = (RectangularRegion) region;

		allocateArrays (chip, rect);

		computeSpotDistance (rect, rows_per_probe);

		// check if number of probes fit in the region
		if ((num_probes = end - start + 1) > dim)
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
		computeProbeDistance (chip, probe_id, start, end, dist_empty);

		// REMOVE THIS!
		curr_cost = 0;
		for (int i = 0; i < dim; i++)
			for (int j = 0; j < dim; j++)
				curr_cost += spot_dist[i * dim + j] * probe_dist[i * dim + j];

		// solve QAP
		sol_cost = solver.solve (dim, spot_dist, probe_dist, perm);

		// place probes according to the optimal permutation
		if (chip instanceof SimpleChip)
			unplaced += applyPermutation ((SimpleChip) chip, rect, probe_id,
											start, end);
		else
			unplaced += applyPermutation ((AffymetrixChip) chip, rect, probe_id,
											start, end);

		return unplaced;
	}

	/**
	 * document this
	 */
	protected void allocateArrays (Chip chip, RectangularRegion r)
	{
		int		new_dim, num_rows, num_cols, rows_per_probe;

		num_rows = r.last_row - r.first_row + 1;
		num_cols = r.last_col - r.first_col + 1;

		if (chip instanceof AffymetrixChip)
			rows_per_probe = 2;
		else
			rows_per_probe = 1;

		new_dim = ((int) (num_rows / rows_per_probe)) * num_cols;

		if (new_dim != dim)
		{
			this.dim = new_dim;
			spot_dist = new int [dim * dim];
			probe_dist = new int [dim * dim];
			perm = new int [dim];
		}
	}

	/**
	 * document this
	 */
	protected void computeSpotDistance (RectangularRegion r, int rows_per_probe)
	{
		int	s1, s1_row, s1_col, s2, s2_row, s2_col, v_dist, h_dist, weight;

		// check if needs to recompute matrix
		if (r.last_row - r.first_row + 1 == last_height &&
			r.last_col - r.first_col + 1 == last_width)
		{
			// no (same dimension of last computed region)
			return;
		}

		s1 = 0;
		s1_row = r.first_row;
		s1_col = r.first_col;
		while (true)
		{
			s2 = s1;
			s2_row = s1_row;
			s2_col = s1_col;

			while (true)
			{
				v_dist = s2_row - s1_row;
				h_dist = s2_col - s1_col;

				weight = 0;

				if (Math.abs(h_dist) <= 3 && Math.abs(v_dist) <= 3)
				{
					// we need integer values
					weight = (int) (1000 * LayoutEvaluation.WEIGHT_DIST[3 + v_dist][3 + h_dist]);
				}

				spot_dist[s1 * dim + s2] = weight;
				spot_dist[s2 * dim + s1] = weight;

				// s2: next spot
				if (++s2_col > r.last_col)
				{
					if (++s2_row * rows_per_probe > r.last_row)
						break;

					s2_col = r.first_col;
				}
				s2++;
			}

			// s1: next spot
			if (++s1_col > r.last_col)
			{
				if (++s1_row * rows_per_probe > r.last_row)
					break;

				s1_col = r.first_col;
			}
			s1++;
		}
	}

	/**
	 * document this
	 */
	protected void computeProbeDistance (Chip chip, int probe_id[],
		int start, int end, int dist_empty)
	{
		int i, j, dist, id_i;

		for (i = start; i <= end; i++)
		{
			probe_dist[i * dim + i] = 0;

			id_i = probe_id[i];

			for (j = i + 1; j <= end; j++)
			{
				// to do: DISTANCE MUST TAKE INTO ACCOUNT POSITION-DEPENDENT WEIGHTS
				// dist = LayoutEvaluation.weightedDistance (chip, id_i, probe_id[j]);

				dist = LayoutEvaluation.hammingDistance (chip, id_i, probe_id[j]);

				probe_dist[i * dim + j] = dist;
				probe_dist[j * dim + i] = dist;
			}

			for (j = end + 1; j < start + dim; j++)
			{
				probe_dist[i * dim + j] = dist_empty;
				probe_dist[j * dim + i] = dist_empty;
			}
		}

		for (i = end + 1; i < start + dim; i++)
		{
			probe_dist[i * dim + i] = 0;

			for (j = i + 1; j < start + dim; j++)
			{
				probe_dist[i * dim + j] = 0;
				probe_dist[j * dim + i] = 0;
			}
		}
	}

	/**
	 * package access
	 */
	protected int applyPermutation (SimpleChip chip, RectangularRegion region,
		int probe_id[], int start, int end)
	{
		// to do

		return 0;
	}

	/**
	 * package access
	 */
	protected int applyPermutation (AffymetrixChip chip, RectangularRegion region,
		int probe_id[], int start, int end)
	{
		int	r, c, num_probes, unplaced = 0;

		num_probes = end - start + 1;
		r = region.first_row;
		c = region.first_col;

		for (int i = 0; i < dim; i++)
		{
			if (chip.isFixedSpot(r, c) ||  chip.isFixedSpot(r + 1, c))
			{
				unplaced++;
			}
			else
			{
				if (perm[i] > num_probes)
				{
					// empty probe
					chip.spot[r][c] = chip.EMPTY_SPOT;
					chip.spot[r + 1][c] = chip.EMPTY_SPOT;
				}
				else
				{
					chip.spot[r][c] = probe_id[start + perm[i] - 1];
					chip.spot[r + 1][c] = probe_id[start + perm[i] - 1] + 1;
				}
			}

			// should be row-by-row, left to right
			r += 2;

			if (r > region.last_row)
			{
				c ++;
				r = region.first_row;
			}
		}

		return unplaced;
	}
}
