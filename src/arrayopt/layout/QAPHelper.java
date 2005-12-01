/*
 * QAPHelper.java
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
public class QAPHelper
{
	/**
	 * package access
	 */
	static void computeSpotDistance (SimpleChip c, int dim,
		RectangularRegion r, int m[])
	{
		// SimpleChip objects use 1 row per probe
		computeSpotDistance (dim, r, 1, m);
	}

	/**
	 * package access
	 */
	static void computeSpotDistance (AffymetrixChip c, int dim,
		RectangularRegion r, int m[])
	{
		// AffymetrixChip objects use 2 row per probe (pair)
		computeSpotDistance (dim, r, 2, m);
	}

	/**
	 *
	 */
	protected static void computeSpotDistance (int dim, RectangularRegion r,
		int rows_per_probe, int m[])
	{
		int	s1, s1_row, s1_col, s2, s2_row, s2_col, v_dist, h_dist, weight;

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

				m[s1 * dim + s2] = weight;
				m[s2 * dim + s1] = weight;

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
	 * package access
	 */
	static void computeProbeDistance (SimpleChip chip, int dim,
		int probe_id[], int start, int end, int m[])
	{
		computeProbeDistance (chip, dim, probe_id, start, end, m,
			chip.getProbeLength());
	}

	/**
	 * package access
	 */
	static void computeProbeDistance (AffymetrixChip chip, int dim,
		int probe_id[], int start, int end, int m[])
	{
		computeProbeDistance (chip, dim, probe_id, start, end, m,
			chip.getProbeLength() + 1);
	}

	/**
	 *
	 */
	protected static void computeProbeDistance (Chip chip, int dim,
		int probe_id[], int start, int end, int m[], int dist_empty)
	{
		int i, j, dist, id_i;

		for (i = start; i <= end; i++)
		{
			m[i * dim + i] = 0;

			id_i = probe_id[i];

			for (j = i + 1; j <= end; j++)
			{
				// to do: DISTANCE MUST TAKE INTO ACCOUNT POSITION-DEPENDENT WEIGHTS
				// dist = LayoutEvaluation.weightedDistance (chip, id_i, probe_id[j]);

				dist = LayoutEvaluation.hammingDistance (chip, id_i, probe_id[j]);

				m[i * dim + j] = dist;
				m[j * dim + i] = dist;
			}

			for (j = end + 1; j < start + dim; j++)
			{
				m[i * dim + j] = dist_empty;
				m[j * dim + i] = dist_empty;
			}
		}

		for (i = end + 1; i < start + dim; i++)
		{
			m[i * dim + i] = 0;

			for (j = i + 1; j < start + dim; j++)
			{
				m[i * dim + j] = 0;
				m[j * dim + i] = 0;
			}
		}
	}

	/**
	 * package access
	 */
	static void applyPermutation (SimpleChip chip, RectangularRegion region,
		int dim, int probe_id[], int start, int end, int[] perm)
	{
		// to do
	}

	/**
	 * package access
	 */
	static void applyPermutation (AffymetrixChip chip,	RectangularRegion region,
		int dim, int probe_id[], int start, int end, int[] perm)
	{
		int	r, c, num_probes;

		num_probes = end - start + 1;
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
				chip.spot[r][c] = probe_id[start + perm[i] - 1];
				chip.spot[r + 1][c] = probe_id[start + perm[i] - 1] + 1;
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
