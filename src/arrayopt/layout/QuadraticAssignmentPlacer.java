/*
 * QuadraticAssignmentPlacer.java
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
 * This class models and solves the placement problem as an instance of a
 * quadratic assignment problem (QAP). This class makes a parallel between the
 * probe placement problem and a classic example of a QAP known as the facility
 * location problem, which is the problem of assigning facilities to locations.
 * The facility location problem has two components, a flow matrix (which stores
 * the flow of materials between the facilities) and a distance matrix (which
 * contains the distance between the locations). The flow of materials has a
 * cost, and thus the problem is to find an assignment of facilities to
 * locations with minimum cost.</P>
 * 
 * <P>This class can be combined with several QAP solvers available in the
 * {@link arrayopt.qap} package to work as a {@link PlacementAlgorithm} or as a
 * {@link FillingAlgorithm}, for instance. The placement problem is modelled as
 * decribed in the paper:<BR>
 * 
 * TODO add reference to the paper with QAP modelling</P>
 * 
 * <P>Two matrices are created: the flow matrix and the distance matrix
 * (mimicking the facility location). These matrices are then passed to the
 * {@link arrayopt.qap.QAPSolverAlgorithm} that was given to the constructor of
 * this class. The solution is a permutation that is mapped back to an
 * assignment of probes to spots.</P>
 * 
 * <P>The matrices can be constructed to model the placement problem as to
 * minimize the sum of border lengths ({@link #MODE_BORDER_LENGTH}) or conflict
 * indices ({@link #MODE_CONFLICT_INDEX}). The public constants are used to
 * configure the class behavior at instantiation time.</P> 
 * 
 * <P>Note that the QAP solvers, at the moment, only work with integer values.
 * This is not a problem when the goal is to minimize the total border length.
 * For conflict index minimization, however, it means that it is necessary to
 * map the decimal numbers returned by {@link ConflictIndex} to integer numbers.
 * This is done by multiplying the values by a constant and rounding to the
 * nearest integer, which may impair the results.</P>
 *
 */
public class QuadraticAssignmentPlacer implements LayoutAlgorithm,
	FillingAlgorithm, IteractiveOptimizationAlgorithm
{
	// TODO handle fixed spots
	
	// TODO swap the flow and distance matrix (to be consistent with the
	// approach described in the paper)
	
	/**
	 * Constant to indicate that the goal is to reduce the sum of border
	 * lengths. 
	 */
	public static final int MODE_BORDER_LENGTH = 0;

	/**
	 * Constant to indicate that the goal is to reduce the sum of conflict
	 * indices.
	 */
	public static final int MODE_CONFLICT_INDEX = 1;

	/**
	 * An instance of a {@link QAPSolverAlgorithm} used to solve the resulting
	 * QAP instances. 
	 */
	private QAPSolverAlgorithm solver;

	/**
	 * Indicates whether the goal is to minimize the sum of border lengths or
	 * conflict indices.
	 */
	private int mode;

	/**
	 * Constant used to map double values to integers.
	 */
	private static final int DOUBLE2INT_MULT = 100;

	/**
	 * Distance between the spots (flow matrix of QAP).
	 */
	private int[] spot_dist;

	/**
	 * Distance between the probes (distance matrix of QAP).
	 */
	private int[] probe_dist;

	/**
	 * Permutation representing the (optimal or approximate) solution to the
	 * resulting QAP. 
	 */
	private int[] perm;

	/**
	 * IDs of probes that are to be assigned to spots.
	 */
	private int[] p_id;

	/**
	 * Dimension of current or QAP instance being solved (or last solved
	 * instance).
	 */
	private int dim;

	/**
	 * Region of the chip where probes must be placed.
	 */
	private RectangularRegion region;

	/**
	 * Height (number of rows) of the (rectangular) region related to the
	 * current or last QAP instance.
	 */
	private int last_height;

	/**
	 * Width (number of columns) of the (rectangular) region related to the
	 * current or last QAP instance.
	 */
	private int last_width;

	/**
	 * Number of rows taken by each probe of the chip, i.e. 1 in case of a
	 * {@link SimpleChip} or 2 in case of a {@link AffymetrixChip}.
	 */
	private int rows_per_probe;

	/**
	 * Creates an QAP placer using a given {@link QAPSolverAlgorithm} with the
	 * desired minimization goal.
	 *   
	 * @param solver an instance of QAP solver algorithm
	 * @param mode border length or conflict index minimization
	 */
	public QuadraticAssignmentPlacer (QAPSolverAlgorithm solver, int mode)
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
	 * Returns the current minimization mode, i.e. border length
	 * ({@link #MODE_BORDER_LENGTH}) or conflict index
	 * ({@link #MODE_CONFLICT_INDEX}).
	 * 
	 * @return current minimization mode
	 */
	public int getMode()
	{
		return mode;
	}

	/**
	 * Sets the minimization mode, i.e. border length
	 * ({@link #MODE_BORDER_LENGTH}) or conflict index
	 * ({@link #MODE_CONFLICT_INDEX}).
	 * 
	 * @param mode minimization mode
	 */
	public void setMode (int mode)
	{
		this.mode = mode;
	}

	/**
	 * Creates the layout of a chip using the {@link QAPSolverAlgorithm} with
	 * the aim of minimizing the sum of border lengths or conflict indices.
	 * 
	 * @param chip chip to be designed
	 * @return number of unplaced probes
	 */
	public void changeLayout (Chip chip)
	{
		int		id[];

		// reset current layout (if any)
		chip.resetLayout();

		// get list of all probes
		id = chip.getAllProbes ();
		
		fillRegion (chip, chip.getChipRegion(), id, 0, id.length - 1);
	}
	
	/**
	 * Optimizes the chip's current layout using the {@link QAPSolverAlgorithm}
	 * with the aim of minimizing the sum of border lengths or conflict indices.
	 * 
	 * @param chip chip to be optimized 
	 */
	public void optimizeLayout (Chip chip)
	{
		optimizeLayout (chip, chip.getChipRegion());
	}
	
	/**
	 * Optimizes the layout of a region of the chip using the
	 * {@link QAPSolverAlgorithm} with the aim of minimizing the sum of border
	 * lengths or conflict indices.
	 * 
	 * @param chip a chip 
	 * @param r region to be optimized
	 * @return reduction in border length or conflict index 
	 */
	public float optimizeLayout (Chip chip, Region r)
	{
		int		num_probes, unplaced = 0;
		long	curr_cost, sol_cost;
		
		// prepare object to handle the current problem's dimension
		configure (chip, r);
		
		// collect the probe IDs currently placed on the region
		num_probes = collectProbeIDs (chip);
		
		computeProbeDistance (chip, this.p_id, 0, num_probes - 1);
		
		curr_cost = solver.computeCost (dim, spot_dist, probe_dist, perm);

		sol_cost = solver.solve (dim, spot_dist, probe_dist, perm);
		
		// place probes according to the optimal permutation
		if (chip instanceof SimpleChip)
		{
			unplaced += applyPermutation ((SimpleChip) chip, this.p_id,
											0, num_probes - 1);
		}
		else
		{
			unplaced += applyPermutation ((AffymetrixChip) chip, this.p_id,
											0, num_probes - 1);
		}
		
		if (unplaced > 0)
			throw new IllegalStateException (unplaced +
					" probes lost in the optimization.");
			
		// return the improvement rate
		return (curr_cost - sol_cost) / (float) curr_cost;
	}

	/**
	 * Fills a given region of a chip with a set of probes using the
	 * {@link QAPSolverAlgorithm} with the aim of minimizing the sum of border
	 * lengths or conflict indices.
	 * 
	 * @param chip a chip
	 * @param r region where probes must be placed
	 * @param probe_id set of probe IDs
	 * @return number of unplaced probes
	 */
	public int fillRegion (Chip chip, Region r, int probe_id[])
	{
		return fillRegion (chip, r, probe_id, 0, probe_id.length - 1);
	}

	/**
	 * Fills a given region of a chip with a set of probes using the
	 * {@link QAPSolverAlgorithm} with the aim of minimizing the sum of border
	 * lengths or conflict indices.
	 * 
	 * @param chip a chip
	 * @param r region where probes must be placed
	 * @param probe_id array containing the probe IDs
	 * @param start first element of the array to be placed
	 * @param end last element of the array to be placed
	 * @return number of unplaced probes
	 */
	public int fillRegion (Chip chip, Region r, int probe_id[], int start,
		int end)
	{
		int num_probes, unplaced;
		
		// prepare internal structures
		configure (chip, r);
		
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
		
		// create initial permutation
		for (int i = 0; i < dim; i++)
			this.perm[i] = i;
		
		// compute probe distance matrix
		computeProbeDistance (chip, probe_id, start, end);

		// solve QAP
		solver.solve (dim, spot_dist, probe_dist, perm);

		// place probes according to the optimal permutation
		if (chip instanceof SimpleChip)
		{
			unplaced += applyPermutation ((SimpleChip) chip, probe_id,
											start, end);
		}
		else
		{
			unplaced += applyPermutation ((AffymetrixChip) chip, probe_id,
											start, end);
		}

		return unplaced;
	}

	/**
	 * Prepare the internal data structures to handle the QAP instance to be
	 * solved.
	 */
	private void configure (Chip chip, Region r)
	{
		int		curr_height, curr_width, curr_dim;
		
		if (!(r instanceof RectangularRegion))
			throw new IllegalArgumentException
				("Only rectangular regions are supported.");
		
		this.region = (RectangularRegion) r;

		if (chip instanceof SimpleChip)
		{
			// simple chips use 1 row per probe
			this.rows_per_probe = 1;
		}
		else if (chip instanceof AffymetrixChip)
		{
			// affymetrix chips use 2 rows per probe (pair)
			this.rows_per_probe = 2;
		}
		else
			throw new IllegalArgumentException
				("Unsupported chip type.");

		// problem's dimension (number of spots in region)
		curr_height = (region.last_row - region.first_row + 1) / rows_per_probe;
		curr_width = region.last_col - region.first_col + 1;

		if (curr_height == this.last_height && curr_width == this.last_width)
			// region dimensions are the same as last time
			return;
		
		// alocate new arrays if current problem's dimension
		// is not equal to the last solved problem
		if ((curr_dim = curr_height * curr_width) != this.dim)
		{
			this.dim = curr_dim;
			this.spot_dist = new int [dim * dim];
			this.probe_dist = new int [dim * dim];
			this.perm = new int [dim];
			this.p_id = new int [dim];
		}
		
		// recompute spot distance matrix
		computeSpotDistance ();
			
		// save dimensions to allow future
		// reuse of the spot distance matrix
		this.last_height = curr_height;
		this.last_width = curr_width;
	}
	
	/**
	 * Collects the probe IDs currently placed on the chip.
	 */
	private int collectProbeIDs (Chip chip)
	{
		int r, c, n, last;
		
		n = 0;
		r = region.first_row;
		c = region.first_col;
		last = dim - 1;

		for (int s = 0; c <= region.last_col; s++)
		{
			// cannot handle fixed spots at the moment
			if (chip.isFixedSpot(r, c))
				throw new IllegalArgumentException
					("Current implementation cannot handle fixed spots.");

			if (chip.spot[r][c] != Chip.EMPTY_SPOT)
			{
				this.perm[s] = n;
				this.p_id[n] = chip.spot[r][c];
				n++;
			}
			else
			{
				this.perm[s] = last--;
			}
			
			r += rows_per_probe;

			if (r > region.last_row)
			{
				c ++;
				r = region.first_row;
			}
		}
		
		// return number of probes found
		return n;
	}

	/**
	 * Compute the distance between the spots (the flow matrix).
	 */
	private void computeSpotDistance ()
	{
		int	ci_dim, s1, s1_row, s1_col, s2, s2_row, s2_col, v_dist, h_dist, d;
		
		ci_dim = ConflictIndex.dimConflictRegion();
		
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
					if (h_dist <= ci_dim && v_dist <= ci_dim)
						// have to use integer values 
						d = (int) Math.round(DOUBLE2INT_MULT *
									ConflictIndex.distanceWeight(s1_row, s1_col,
										s2_row, s2_col));
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
	 * Compute the distance between the probe embeddings (the distance matrix).
	 */
	private void computeProbeDistance (Chip chip, int id[], int start,
			int end)
	{
		if (mode == MODE_BORDER_LENGTH)
			computeBorderLengthDistance (chip, id, start, end);
		else // if (mode == MODE_CONFLICT_INDEX)
			computeConflictIndexDistance (chip, id, start, end);
	}

	/**
	 * Compute the distance between the probe embeddings (the distance matrix)
	 * according to the border length minimization model.
	 */
	private void computeBorderLengthDistance (Chip chip, int id[], int start,
			int end)
	{
		int i, j, num_probes, dist, id_i;
		
		num_probes = end - start + 1;

		for (i = 0; i < num_probes; i++)
		{
			probe_dist[i * dim + i] = 0;

			id_i = id[start + i];

			for (j = i + 1; j < num_probes; j++)
			{
				dist = LayoutEvaluation.hammingDistance (chip, id_i,
						id[start + j]);

				probe_dist[i * dim + j] = dist;
				probe_dist[j * dim + i] = dist;
			}

			for (j = num_probes; j < dim; j++)
			{
				// distance to "empty" probes is always zero 
				probe_dist[i * dim + j] = probe_dist[j * dim + i] = 0;
			}
		}

		// distance between "empty" probes is always zero
		for (i = num_probes; i < dim; i++)
		{
			probe_dist[i * dim + i] = 0;

			for (j = i + 1; j < dim; j++)
				probe_dist[i * dim + j] = probe_dist[j * dim + i] = 0;
		}
	}

	/**
	 * Compute the distance between the probe embeddings (the distance matrix)
	 * according to the conflict index minimization model.
	 */
	private void computeConflictIndexDistance (Chip chip, int id[], int start,
			int end)
	{
		int i, j, num_probes, id_i;
		double d;
		
		num_probes = end - start + 1;

		for (i = 0; i < num_probes; i++)
		{
			probe_dist[i * dim + i] = 0;

			id_i = id[start + i];

			for (j = i + 1; j < num_probes; j++)
			{
				d = LayoutEvaluation.weightedDistance (chip, id_i,
						id[start + j]);
				
				// convert to integer
				probe_dist[i * dim + j] = (int) Math.round(d * DOUBLE2INT_MULT);

				d = LayoutEvaluation.weightedDistance (chip, id[start + j],
						id_i);

				// convert to integer
				probe_dist[j * dim + i] = (int) Math.round(d * DOUBLE2INT_MULT);
			}

			for (j = num_probes; j < dim; j++)
			{
				// distance to "empty" probes is always zero
				probe_dist[i * dim + j] = probe_dist[j * dim + i] = 0;
			}
		}

		// distance between "empty" probes is always zero
		for (i = num_probes; i < dim; i++)
		{
			probe_dist[i * dim + i] = 0;

			for (j = i + 1; j < dim; j++)
				probe_dist[i * dim + j] = probe_dist[j * dim + i] = 0;
		}
	}

	/**
	 * Use the QAP solution to assign probes to spots.
	 */
	private int applyPermutation (SimpleChip chip, int id[], int start,
		int end)
	{
		int	r, c, num_probes, unplaced = 0;

		num_probes = end - start + 1;
		r = region.first_row;
		c = region.first_col;

		for (int i = 0; i < dim; i++)
		{
			if (perm[i] >= num_probes)
			{
				// empty probe
				chip.spot[r][c] = Chip.EMPTY_SPOT;
			}
			else
			{
				chip.spot[r][c] = id[start + perm[i]];
			}

			if (++r > region.last_row)
			{
				c ++;
				r = region.first_row;
			}
		}

		return unplaced;
	}

	/**
	 * Use the QAP solution to assign probes to spots.
	 */
	private int applyPermutation (AffymetrixChip chip, int id[], int start,
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
				if (perm[i] >= num_probes)
				{
					// empty probe
					chip.spot[r][c] = Chip.EMPTY_SPOT;
					chip.spot[r + 1][c] = Chip.EMPTY_SPOT;
				}
				else
				{
					chip.spot[r][c] = id[start + perm[i]];
					chip.spot[r + 1][c] = id[start + perm[i]] + 1;
				}
			}

			if ((r += 2) > region.last_row)
			{
				c ++;
				r = region.first_row;
			}
		}

		return unplaced;
	}
}
