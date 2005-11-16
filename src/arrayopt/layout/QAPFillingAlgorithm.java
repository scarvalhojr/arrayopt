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

/**
 *
 */
public abstract class QAPFillingAlgorithm extends FillingAlgorithm
{
	protected int[] spot_dist;

	protected int[] probe_dist;

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
		RectangularRegion	rect;
		int					dim, num_probes, unplaced;
		int[]				perm;
		long				curr_cost, sol_cost, comp_cost;

		if (!(region instanceof RectangularRegion))
			throw new IllegalArgumentException ("Only rectangular regions are supported.");

		rect = (RectangularRegion) region;

		// create and compute spot distance matrix
		dim = computeSpotDistanceMatrix (chip, rect);

		// create probe distance matrix
		createProbeDistanceMatrix (dim);

		// check if number of probes fit in the region
		num_probes = last_probe - first_probe + 1;

		if (num_probes > dim)
		{
			// exceeding probes have to be ignored
			unplaced = num_probes - dim;
			last_probe = first_probe + dim - 1;
		}
		else
		{
			// all probes can be placed
			unplaced = 0;
		}

		if (chip instanceof SimpleChip)
			computeProbeDistanceMatrix ((SimpleChip) chip, dim, first_probe, last_probe);

		else if (chip instanceof AffymetrixChip)
			computeProbeDistanceMatrix ((AffymetrixChip) chip, dim, first_probe, last_probe);

		else
			throw new IllegalArgumentException ("Unsupported chip type.");

		curr_cost = 0;
		for (int i = 0; i < dim; i++)
			for (int j = 0; j < dim; j++)
				curr_cost += spot_dist[i * dim + j] * probe_dist[i * dim + j];

		// run sub-classes' specific method
		perm = new int[dim];
		sol_cost = solveQAP (dim, spot_dist, probe_dist, perm);

		// place probes according to the optimal permutation
		if (chip instanceof SimpleChip)
			applyPermutation ((SimpleChip) chip, rect, dim, first_probe, last_probe, perm);
		else
			applyPermutation ((AffymetrixChip) chip, rect, dim, first_probe, last_probe, perm);

		return unplaced;
	}

	protected int computeSpotDistanceMatrix (Chip chip, RectangularRegion region)
	{
		int		dim, num_rows, num_cols, rows_per_probe;
		int		s1, s1_row, s1_col, s2, s2_row, s2_col, v_dist, h_dist, weight;

		if (chip instanceof AffymetrixChip)
			rows_per_probe = 2;
		else
			rows_per_probe = 1;

		num_rows = region.last_row - region.first_row + 1;
		num_cols = region.last_col - region.first_col + 1;

		dim = ((int) (num_rows / rows_per_probe)) * num_cols;

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

		s1 = 0;
		s1_row = region.first_row;
		s1_col = region.first_col;
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
				if (++s2_col > region.last_col)
				{
					if (++s2_row * rows_per_probe > region.last_row)
						break;

					s2_col = region.first_col;
				}
				s2++;
			}

			// s1: next spot
			if (++s1_col > region.last_col)
			{
				if (++s1_row * rows_per_probe > region.last_row)
					break;

				s1_col = region.first_col;
			}
			s1++;
		}

		return dim;
	}

	protected void createProbeDistanceMatrix (int dim)
	{
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
	}

	protected void computeProbeDistanceMatrix (SimpleChip chip, int dim, int first_probe, int last_probe)
	{
		// to do
	}

	protected void computeProbeDistanceMatrix (AffymetrixChip chip, int dim, int first_probe, int last_probe)
	{
		int i, j, dist, dist_empty, id_i;

		dist_empty = chip.getProbeLength() + 1;

		for (i = first_probe; i <= last_probe; i++)
		{
			probe_dist[i * dim + i] = 0;

			id_i = chip.probe_list[i];

			for (j = i + 1; j <= last_probe; j++)
			{
				// to do: DISTANCE MUST TAKE INTO ACCOUNT POSITION-DEPENDENT WEIGHTS
				// dist = LayoutEvaluation.weightedDistance (chip, id_i, chip.probe_list[j]);

				dist = LayoutEvaluation.hammingDistance (chip, id_i, chip.probe_list[j]);

				probe_dist[i * dim + j] = dist;
				probe_dist[j * dim + i] = dist;
			}

			for (j = last_probe + 1; j < first_probe + dim; j++)
			{
				probe_dist[i * dim + j] = dist_empty;
				probe_dist[j * dim + i] = dist_empty;
			}
		}

		for (i = last_probe + 1; i < first_probe + dim; i++)
		{
			probe_dist[i * dim + i] = 0;

			for (j = i + 1; j < first_probe + dim; j++)
			{
				probe_dist[i * dim + j] = 0;
				probe_dist[j * dim + i] = 0;
			}
		}
	}

	public abstract long solveQAP (int dim, int dist[], int flow[], int sol[]);

	protected void applyPermutation (SimpleChip chip, RectangularRegion region, int dim, int first_probe, int last_probe, int[] perm)
	{
	}

	protected void applyPermutation (AffymetrixChip chip, RectangularRegion region, int dim, int first_probe, int last_probe, int[] perm)
	{
		int	r, c, num_probes;

		num_probes = last_probe - first_probe + 1;
		r = region.first_row;
		c = region.first_col;

		for (int i = 0; i < dim; i++)
		{
			if (perm[i] > num_probes)
			{
				// empty probe
				chip.spot[r][c] = chip.EMPTY_SPOT;
				chip.spot[r + 1][c] = chip.EMPTY_SPOT;
			}
			else
			{
				chip.spot[r][c] = chip.probe_list[first_probe + perm[i] - 1];
				chip.spot[r + 1][c] = chip.probe_list[first_probe + perm[i] - 1] + 1;
			}

			// should be row-by-row, left to right
			r += 2;

			if (r > region.last_row)
			{
				c ++;
				r = region.first_row;
			}
		}
	}
}
