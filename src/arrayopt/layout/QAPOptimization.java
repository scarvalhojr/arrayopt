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
	protected int mode;

	/**
	 * document this
	 */
	public static final int MODE_BORDER_LENGTH = 0;

	/**
	 * document this
	 */
	public static final int MODE_CONFLICT_INDEX = 1;

	/**
	 * document this
	 */
	public static final int DEFAULT_MODE = MODE_CONFLICT_INDEX;

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
	 */
	protected int[] p_id;

	/**
	 * document this
	 *
	 * dimension of current/last instance
	 */
	protected int dim;

	/**
	 * document this
	 */
	protected RectangularRegion region;

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
	protected int rows_per_probe;

	/**
	 * document this
	 */
	protected int empty_dist;

	/**
	 * document this
	 */
	public QAPOptimization (QAPSolverAlgorithm solver)
	{
		this (solver, DEFAULT_MODE);
	}

	/**
	 * document this
	 */
	public QAPOptimization (QAPSolverAlgorithm solver, int mode)
	{
		this.solver = solver;
		this.dim = -1;
		this.last_height = -1;
		this.last_width = -1;
		
		switch (mode)
		{
			case MODE_BORDER_LENGTH:
			case MODE_CONFLICT_INDEX:
				this.mode = mode;
				break;
			
			default:
				throw new IllegalArgumentException
					("Illegal value for argument 'mode'.");
		}
	}

	/**
	 * document this
	 *
	 * should be synchronized
	 */
	public float optimizeLayout (Chip chip, Region region)
	{
		int		num_probes, unplaced;
		long	curr_cost, sol_cost, new_cost;
		
		// prepare object to handle the
		// current problem's dimension
		configure (chip, region);
		
		// collect the probe IDs currently
		// placed on the region
		num_probes = collectProbeIDs (chip);
		
		computeProbeDistance (chip, this.p_id, 0, num_probes - 1);

		// compute current cost
		curr_cost = 0;
		for (int i = 0; i < dim; i++)
			for (int j = 0; j < dim; j++)
				curr_cost += spot_dist[i * dim + j] * probe_dist[i * dim + j];

		System.err.println ("current cost: " + curr_cost);
		
		// solve QAP
		sol_cost = solver.solve (dim, spot_dist, probe_dist, perm);
		
		System.err.println ("solution cost: " + sol_cost);

		// compute new cost
		new_cost = 0;
		for (int i = 0; i < dim; i++)
			for (int j = 0; j < dim; j++)
				new_cost += spot_dist[i * dim + j] * probe_dist[(perm[i] - 1) * dim + (perm[j] - 1)];
		
		if (new_cost != sol_cost)
			System.err.println ("CAREFUL: new cost is different: " + new_cost);
		
		if (sol_cost > curr_cost)
		{
			// current layout could not be improved
			return 0;
		}

		// place probes according to the optimal permutation
		unplaced = 0;
		if (chip instanceof SimpleChip)
			unplaced += applyPermutation ((SimpleChip) chip, this.p_id,
											0, num_probes - 1);
		else
			unplaced += applyPermutation ((AffymetrixChip) chip, this.p_id,
											0, num_probes - 1);
		if (unplaced > 0)
			System.err.println (unplaced + " probes lost in the relocation");
			
		// return the improvement rate
		return (curr_cost - sol_cost) / (float) curr_cost;
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
	 *
	 * should be synchronized
	 */
	public int fillRegion (Chip chip, Region region, int probe_id[], int start,
		int end)
	{
		int		num_probes, unplaced;
		long	curr_cost, sol_cost, comp_cost;
		
		System.err.println("fillRegion; region = " + region);
		
		// prepare object to handle the
		// current problem's dimension
		configure (chip, region);
		
		// check if number of probes fit in the region
		if ((num_probes = end - start + 1) > dim)
		{
			// exceeding probes have to be ignored
			unplaced = num_probes - dim;
			end -= unplaced;
		}
		else
		{
			// all probes can be placed
			unplaced = 0;
		}
		
		// compute probe distance matrix
		computeProbeDistance (chip, probe_id, start, end);

		// REMOVE THIS!
		curr_cost = 0;
		for (int i = 0; i < dim; i++)
			for (int j = 0; j < dim; j++)
				curr_cost += spot_dist[i * dim + j] * probe_dist[i * dim + j];

		// solve QAP
		sol_cost = solver.solve (dim, spot_dist, probe_dist, perm);

		// place probes according to the optimal permutation
		if (chip instanceof SimpleChip)
			unplaced += applyPermutation ((SimpleChip) chip, probe_id,
											start, end);
		else
			unplaced += applyPermutation ((AffymetrixChip) chip, probe_id,
											start, end);

		return unplaced;
	}

	/**
	 * document this
	 */
	protected void configure (Chip chip, Region r)
	{
		int	curr_height, curr_width, curr_dim;
		int					num_probes, unplaced;
		int					start, end;
		long				curr_cost, sol_cost, new_cost;
		
		if (!(r instanceof RectangularRegion))
			throw new IllegalArgumentException
				("Only rectangular regions are supported.");
		
		// save references to the problem's data
		this.region = (RectangularRegion) r;

		if (chip instanceof SimpleChip)
		{
			// simple chips use 1 row per probe
			this.rows_per_probe = 1;

			// distance from empty string to any probe
			// equals probe length
			this.empty_dist = chip.getProbeLength();
		}
		else if (chip instanceof AffymetrixChip)
		{
			// affymetrix chips use 2 rows per probe (pair)
			this.rows_per_probe = 2;

			// distance from empty string to any probe
			// equals probe length + 1 since probes appear
			// in pairs which differ in the middle base
			this.empty_dist = chip.getProbeLength() + 1;
		}
		else
			throw new IllegalArgumentException
				("Unsupported chip type.");

		// problem's dimension (number of spots in region)
		curr_height = (region.last_row - region.first_row + 1) / rows_per_probe;
		curr_width = region.last_col - region.first_col + 1;
		curr_dim = curr_height * curr_width;

		// alocate new arrays if current dimension differs
		// from the dimension of the last solved problem
		if (curr_dim != this.dim)
		{
			this.dim = curr_dim;
			this.spot_dist = new int [dim * dim];
			this.probe_dist = new int [dim * dim];
			this.perm = new int [dim];
			this.p_id = new int [dim];
		}
		
		// check if need to recompute
		// the spot distance matrix
		if (this.last_height != curr_height || this.last_width != curr_width)
		{
			computeSpotDistance ();
			
			// save dimensions to allow
			// future reuse of the matrix
			this.last_height = curr_height;
			this.last_width = curr_width;
		}
	}
	
	protected int collectProbeIDs (Chip chip)
	{
		int r, c, i;
		
		i = 0;
		r = region.first_row;
		c = region.first_col;

		while (c <= region.last_col)
		{
			// skip fixed and empty spots
			if (!chip.isFixedSpot(r, c) && chip.spot[r][c] != chip.EMPTY_SPOT)
			{
				this.p_id[i++] = chip.spot[r][c];
			}
				
			r += rows_per_probe;

			if (r > region.last_row)
			{
				c ++;
				r = region.first_row;
			}
		}
		
		return i;
	}

	/**
	 * document this
	 */
	protected void computeSpotDistance ()
	{
		int		s1, s1_row, s1_col, s2, s2_row, s2_col, v_dist, h_dist, d;
		
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
				d = 0;
				
				v_dist = Math.abs(s2_row - s1_row);
				h_dist = Math.abs(s2_col - s1_col);
				
				if (mode == MODE_BORDER_LENGTH)
				{
					if ((h_dist == 1 && v_dist == 0) ||
						(h_dist == 0 && v_dist == 1))
							d = 1;
				}
				else // mode == MODE_CONFLICT_INDEX
				{
					if (h_dist <= 3 && v_dist <= 3)
						// we need integer values
						d = (int) (1000 * LayoutEvaluation.WEIGHT_DIST[3 + v_dist][3 + h_dist]);
				}

				spot_dist[s1 * dim + s2] = d;
				spot_dist[s2 * dim + s1] = d;

				// s2: next spot
				if ((s2_row += rows_per_probe) > region.last_row)
				{
					if (++s2_col > region.last_col)
						break;

					s2_row = region.first_row;
				}
				s2++;
			}

			// s1: next spot
			if ((s1_row += rows_per_probe) > region.last_row)
			{
				if (++s1_col > region.last_col)
					break;

				s1_row = region.first_row;
			}
			s1++;
		}
	}

	/**
	 * document this
	 */
	protected void computeProbeDistance (Chip chip, int id[],
		int start, int end)
	{
		int i, j, num_probes, dist, id_i;
		
		num_probes = end - start + 1;

		for (i = 0; i < num_probes; i++)
		{
			probe_dist[i * dim + i] = 0;

			id_i = id[start + i];

			for (j = i + 1; j < num_probes; j++)
			{
				// to do: DISTANCE MUST TAKE INTO ACCOUNT POSITION-DEPENDENT WEIGHTS
				// dist = LayoutEvaluation.weightedDistance (chip, id_i, id[start + j]);

				dist = LayoutEvaluation.hammingDistance (chip, id_i, id[start + j]);

				probe_dist[i * dim + j] = dist;
				probe_dist[j * dim + i] = dist;
			}

			for (j = num_probes; j < dim; j++)
			{
				probe_dist[i * dim + j] = empty_dist;
				probe_dist[j * dim + i] = empty_dist;
			}
		}

		for (i = num_probes; i < dim; i++)
		{
			probe_dist[i * dim + i] = 0;

			for (j = i + 1; j < dim; j++)
			{
				probe_dist[i * dim + j] = 0;
				probe_dist[j * dim + i] = 0;
			}
		}
	}

	/**
	 * package access
	 */
	protected int applyPermutation (SimpleChip chip, int id[], int start,
		int end)
	{
		int	r, c, num_probes, unplaced = 0;

		num_probes = end - start + 1;
		r = region.first_row;
		c = region.first_col;

		for (int i = 0; i < dim; i++)
		{
			if (chip.isFixedSpot(r, c))
			{
				unplaced++;
			}
			else
			{
				if (perm[i] > num_probes)
				{
					// empty probe
					chip.spot[r][c] = chip.EMPTY_SPOT;
				}
				else
				{
					chip.spot[r][c] = id[start + perm[i] - 1];
				}
			}

			// should be row-by-row, left to right
			r ++;

			if (r > region.last_row)
			{
				c ++;
				r = region.first_row;
			}
		}

		return unplaced;
	}

	/**
	 * package access
	 */
	protected int applyPermutation (AffymetrixChip chip, int id[], int start,
		int end)
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
					chip.spot[r][c] = id[start + perm[i] - 1];
					chip.spot[r + 1][c] = id[start + perm[i] - 1] + 1;
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
