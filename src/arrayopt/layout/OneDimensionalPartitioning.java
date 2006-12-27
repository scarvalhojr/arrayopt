/*
 * OneDimensionalPartitioning.java
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
 * This class implements the 1-dimensional Partitioning. The algorithm
 * partitions the chip recursively based on the stated of the embeddings at a
 * particular synthesis step.
 * 
 * <P>The 1-D partitioning uses the first synthesis steps to divide the probe
 * set into sub-sets according to the state of their embeddings (masked or
 * unmasked). The chip is divided into sub-regions (vertical layers)
 * proportionally to the number of masked and unmasked embeddings.</P>
 * 
 * <P>The partitioning of a region stops when one of the following conditions is
 * true: 1) the number of probes in the region is less than two; 2) all
 * synthesis steps have already been used for dividing the probe set; 3) the
 * region width (number of columns) is less than {@link #stop_dim} (the stopping
 * dimension can be specified at instantiation time; default value is given by
 * the {@link #DEFAULT_STOP_DIMENSION} constant).</P>
 * 
 * <P>In the end, each sub-region is filled with its assigned probes using the
 * {@link FillingAlgorithm} specified at instantiation time.</P>
 * 
 * <P>The 1-D Partitioning produces a layout in which the first masks are highly
 * optimized, but the last masks are left with high levels of conflicts. In
 * order to reduce conflicts in the last masks, the stopping dimension can be
 * increased so that the filling algorithm will have more freedom on where to
 * place the probes.</P>
 * 
 * <P>{@link TwoDimensionalPartitioning} extends this algorithm to two
 * dimensions.</P>
 * 
 * <P>TODO take into account fixed spots when checking if probes fit into a
 * partition.</P>
 * 
 * @author Sergio A. de Carvalho Jr.
 */
public class OneDimensionalPartitioning implements LayoutAlgorithm
{
	/**
	 * Default stopping dimension. The stopping dimension is the minimum
	 * dimension of a region (number of rows and column) necessary to perform
	 * a partitioning.
	 */
	private static final int DEFAULT_STOP_DIMENSION = 4;
	
	/**
	 * This constant is used to avoid partitionings of the set of probes that
	 * results in a one very small sub-set. 
	 */
	private static final double MIN_DIV_RATE = .1;
	
	/**
	 * Filling algorithm used to place the probes in each final sub-region.
	 */
	private FillingAlgorithm filler;

	private int stop_dim;

	private Chip chip;
	
	private int probe_id[];
	
	private int rows_per_probe;
	
	/**
	 * Creates an instance of the 1-D Partitioning algorithm with the specified
	 * filling algorithm.
	 * 
	 * @param filler filling algorithm used in the final sub-regions
	 */
	public OneDimensionalPartitioning (FillingAlgorithm filler)
	{
		this(filler, DEFAULT_STOP_DIMENSION);
	}

	/**
	 * Creates an instance of the 1-D Partitioning algorithm with the specified
	 * filling algorithm and stopping dimension.
	 * 
	 * @param filler filling algorithm used in the final sub-regions
	 * @param stop_dim stopping dimension
	 */
	public OneDimensionalPartitioning (FillingAlgorithm filler, int stop_dim)
	{
		this.filler = filler;
		this.stop_dim = (stop_dim < 1) ? 1 : stop_dim;
	}

	/**
	 * Creates a new layout of a microarray chip using the 1-D Partitioning
	 * algorithm.
	 * 
	 * @param chip instance of a microarray chip
	 */
	public void changeLayout (Chip c)
	{
		RectangularRegion region;
		
		// save reference to chip
		this.chip = c;
		
		if (chip instanceof AffymetrixChip)
			rows_per_probe = 2;
		else if (chip instanceof SimpleChip)
			rows_per_probe = 1;
		else
			throw new IllegalArgumentException ("Unsupported chip type.");
		
		// reset current layout (if any)
		chip.resetLayout();

		// get movable probes and chip region
		this.probe_id = chip.getMovableProbes ();
		region = chip.getChipRegion();
		
		divide (region, 0, 0, 0, probe_id.length - 1);
	}

	private void divide (RectangularRegion r, int step, int par,
			int start, int end)
	{
		RectangularRegion l_region, r_region;
		int		probe_div, overflow, m_spots, u_spots;
		int		col_div, m_cols, u_cols, qt_rows;
		double	div_rate;
		
		if (end - start + 1 < 2)
		{
			// insufficient number of probes for partitioning
			fillRegion (r, start, end);
			return;
		}
		
		if (step >= chip.getEmbeddingLength())
		{
			// no more synthesis steps to partition the probe set
			fillRegion (r, start, end);
			return;
		}
		
		if (r.last_col - r.first_col + 1 <= stop_dim)
		{
			// region too small to be partitioned
			fillRegion (r, start, end);
			return;
		}
		
		// split the probes into two groups, M and U, according to
		// whether they are Masked or Unmasked at the current step
		probe_div = divideProbes (step, start, end);
		
		step++;
		
		div_rate = ((double)(probe_div - start)) / (end - start + 1);
		
		if (div_rate < MIN_DIV_RATE)
		{
			divide (r, step, 1 - par, start, end);
			return;
		}
		else if (1 - div_rate < MIN_DIV_RATE)
		{
			divide (r, step, par, start, end);
			return;
		}
		
		// compute how many columns are needed to fit the masked probes
		m_cols = (int) Math.round(div_rate * (r.last_col - r.first_col + 1));

		// make sure partitions get a minimum number of columns
		if (m_cols == 0)
			m_cols++;
		else if (m_cols == r.last_col - r.first_col + 1)
			m_cols--;

		u_cols = (r.last_col - r.first_col + 1) - m_cols;
		
		// make sure probe sets fit into sub-regions
		qt_rows = r.last_row - r.first_row + 1;
		m_spots = m_cols * qt_rows / rows_per_probe;
		u_spots = u_cols * qt_rows / rows_per_probe;
		if ((overflow = (probe_div - start) - m_spots) > 0)
			probe_div -= overflow;
		else if ((overflow = (end - probe_div + 1) - u_spots) > 0)
			probe_div += overflow;

		// region will be divided into two sub-regions: left and right
		
		if (par == 0)
			// masked probes will be assigned to left sub-region
			col_div = r.first_col + m_cols;
		else
			// unmasked probes will be assigned to left sub-region
			col_div = r.first_col + u_cols;

		l_region = new RectangularRegion (r.first_row, r.last_row,
											r.first_col, col_div - 1);
		r_region = new RectangularRegion (r.first_row, r.last_row,
											col_div, r.last_col);

		if (par == 0)
		{
			// assign masked probes to left sub-region
			divide (l_region, step, 0, start, probe_div - 1);

			// assign unmasked probes to right sub-region
			divide (r_region, step, 1, probe_div, end);
		}
		else
		{
			// assign masked probes to right sub-region
			divide (r_region, step, 1, start, probe_div - 1);

			// assign unmasked probes to left sub-region
			divide (l_region, step, 0, probe_div, end);
		}
	}

	private int divideProbes (int step, int start, int end)
	{
		int	p_id, w, mask;

		// which 4-byte word will be interrogated?
		w = (int) Math.floor((double) step / Integer.SIZE);

		// prepare mask to interrogate corresponding bit
		mask = 0x01 << (Integer.SIZE - 1 - (step - w * Integer.SIZE));

		while (start <= end)
		{
			p_id = probe_id[start];

			if ((chip.embed[p_id][w] & mask) == 0)
			{
				// probe is masked at this step:
				// proceed to next embedding on the list
				start ++;
			}
			else
			{
				// probe is unmasked at this step:
				// exchange it with the last element on the list
				probe_id[start] = probe_id[end];
				probe_id[end] = p_id;
				end --;
			}
		}

		return start;
	}
	
	private void fillRegion (RectangularRegion r, int start, int end)
	{
		filler.fillRegion(this.chip, r, this.probe_id, start, end);
	}

	/**
	 * Returns the algorithm's name together with current options.
	 * 
	 * @return algorithm's name and configurable options
	 */
	@Override
	public String toString ()
	{		
		return "2DPartitioning-" + stop_dim + "-" + filler;
	}
}
