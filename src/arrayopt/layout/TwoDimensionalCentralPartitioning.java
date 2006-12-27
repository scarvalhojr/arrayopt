/*
 * TwoDimensionalCentralPartitioning.java
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
 * This class implements a variation of the 2-dimensional Partitioning that
 * optimizes the central masks. Like the {@link TwoDimensionalPartitioning},
 * this algorithm partitions the chip recursively based on the stated of the
 * embeddings at a particular synthesis step. The difference is that here the
 * middle steps are used to partition the probe sets (instead of the first
 * synthesis steps). This results in highly optimized central masks (those that
 * are more likely to synthesize the middle bases of the probes). 
 * 
 * <P>For more information about the basic algorithm, see
 * {@link TwoDimensionalPartitioning}.</P>
 * 
 * <P>TODO take into account fixed spots when checking if probes fit into a
 * partition.</P>
 * 
 * @author Sergio A. de Carvalho Jr.
 */
public class TwoDimensionalCentralPartitioning implements LayoutAlgorithm
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
	 * Creates an instance of the 2-D Partitioning algorithm with central mask
	 * optimization using the the specified filling algorithm.
	 * 
	 * @param filler filling algorithm used in the final sub-regions
	 */
	public TwoDimensionalCentralPartitioning (FillingAlgorithm filler)
	{
		this(filler, DEFAULT_STOP_DIMENSION);
	}

	/**
	 * Creates an instance of the 2-D Partitioning algorithm with central mask
	 * optimization using the specified filling algorithm and stopping
	 * dimension.
	 * 
	 * @param filler filling algorithm used in the final sub-regions
	 * @param stop_dim stopping dimension
	 */
	public TwoDimensionalCentralPartitioning (FillingAlgorithm filler,
			int stop_dim)
	{
		this.filler = filler;
		this.stop_dim = (stop_dim < 1) ? 1 : stop_dim;
	}

	/**
	 * Creates a new layout of a microarray chip using the 2-D Partitioning
	 * with central mask optimization.
	 * 
	 * @param chip instance of a microarray chip
	 */
	public void changeLayout (Chip c)
	{
		RectangularRegion region;
		int m;
		
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
		
		// get middle step
		m = (int) Math.ceil(chip.getEmbeddingLength() / (double) 2);
		
		horizontalDivide (region, m - 1, m, 0, 0, 0, probe_id.length - 1);
	}

	private void horizontalDivide (RectangularRegion r, int hstep, int vstep,
			int hpar, int vpar, int start, int end)
	{
		RectangularRegion t_region, b_region;
		int		probe_div, overflow, m_spots, u_spots;
		int		row_div, m_rows, u_rows, qt_cols;
		double	div_rate;
		
		if (end - start + 1 < 2)
		{
			// insufficient number of probes for partitioning
			fillRegion (r, start, end);
			return;
		}
		
		if (hstep < 0 || hstep >= chip.getEmbeddingLength())
		{
			// no more synthesis steps to partition the probe set
			fillRegion (r, start, end);
			return;
		}
		
		if (r.last_row - r.first_row + 1 <= stop_dim * rows_per_probe)
		{
			if (r.last_col - r.first_col + 1 <= stop_dim)
			{
				// region too small to be partitioned
				fillRegion (r, start, end);
				return;
			}

			// region can still be vertically partitioned
			verticalDivide (r, hstep, vstep, hpar, vpar, start, end);
			return;
		}
		
		// split the probes into two groups, M and U, according to
		// whether they are Masked or Unmasked at the current step
		probe_div = divideProbes (hstep, start, end);
		
		hstep--;
		
		div_rate = ((double)(probe_div - start)) / (end - start + 1);
		
		if (div_rate < MIN_DIV_RATE)
		{
			verticalDivide (r, hstep, vstep, 1 - hpar, vpar, start, end);
			return;
		}
		else if (1 - div_rate < MIN_DIV_RATE)
		{
			verticalDivide (r, hstep, vstep, hpar, vpar, start, end);
			return;
		}
		
		// compute how many rows are needed to fit the masked probes
		m_rows = (int) Math.round(div_rate * (r.last_row - r.first_row + 1));

		// make sure partitions get a minimum number of rows
		if (m_rows == 0)
			m_rows += rows_per_probe;
		else if (m_rows == r.last_row - r.first_row + 1)
			m_rows -= rows_per_probe;

		// make sure partitions get an appropriate number of rows
		if (m_rows % rows_per_probe != 0)
			if (div_rate < .5)
			{
				m_rows += rows_per_probe - (m_rows % rows_per_probe);
			}
			else
			{
				m_rows -= m_rows % rows_per_probe;
			}

		u_rows = (r.last_row - r.first_row + 1) - m_rows;

		// make sure probe sets fit into sub-regions
		qt_cols = r.last_col - r.first_col + 1;
		m_spots = qt_cols * m_rows / rows_per_probe;
		u_spots = qt_cols * u_rows / rows_per_probe;
		if ((overflow = (probe_div - start) - m_spots) > 0)
			probe_div -= overflow;
		else if ((overflow = (end - probe_div + 1) - u_spots) > 0)
			probe_div += overflow;

		// region will be divided into two sub-regions: top and bottom
		
		if (hpar == 0)
			// masked probes will be assigned to top sub-region
			row_div = r.first_row + m_rows;
		else
			// unmasked probes will be assigned to top sub-region
			row_div = r.first_row + u_rows;

		t_region = new RectangularRegion (r.first_row, row_div - 1,
											r.first_col, r.last_col);
		b_region = new RectangularRegion (row_div, r.last_row,
											r.first_col, r.last_col);

		if (hpar == 0)
		{
			// assign masked probes to top sub-region
			verticalDivide (t_region, hstep, vstep, 0, vpar, start,
							probe_div - 1);

			// assign unmasked probes to bottom sub-region
			verticalDivide (b_region, hstep, vstep, 1, vpar, probe_div, end);
		}
		else
		{
			// assign masked probes to bottom sub-region
			verticalDivide (b_region, hstep, vstep, 1, vpar, start,
							probe_div - 1);

			// assign unmasked probes to top sub-region
			verticalDivide (t_region, hstep, vstep, 0, vpar, probe_div, end);
		}
	}

	private void verticalDivide (RectangularRegion r, int hstep, int vstep,
			int hpar, int vpar, int start, int end)
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
		
		if (vstep < 0 || vstep >= chip.getEmbeddingLength())
		{
			// no more synthesis steps to partition the probe set
			fillRegion (r, start, end);
			return;
		}
		
		if (r.last_col - r.first_col + 1 <= stop_dim)
		{
			if (r.last_row - r.first_row + 1 <= stop_dim * rows_per_probe)
			{
				// region too small to be partitioned
				fillRegion (r, start, end);
				return;
			}

			// region can still be horizontally partitioned
			horizontalDivide (r, hstep, vstep, hpar, vpar, start, end);
			return;
		}
		
		// split the probes into two groups, M and U, according to
		// whether they are Masked or Unmasked at the current step
		probe_div = divideProbes (vstep, start, end);
		
		vstep++;
		
		div_rate = ((double)(probe_div - start)) / (end - start + 1);
		
		if (div_rate < MIN_DIV_RATE)
		{
			horizontalDivide (r, hstep, vstep, hpar, 1 - vpar, start, end);
			return;
		}
		else if (1 - div_rate < MIN_DIV_RATE)
		{
			horizontalDivide (r, hstep, vstep, hpar, vpar, start, end);
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
		
		if (vpar == 0)
			// masked probes will be assigned to left sub-region
			col_div = r.first_col + m_cols;
		else
			// unmasked probes will be assigned to left sub-region
			col_div = r.first_col + u_cols;

		l_region = new RectangularRegion (r.first_row, r.last_row,
											r.first_col, col_div - 1);
		r_region = new RectangularRegion (r.first_row, r.last_row,
											col_div, r.last_col);

		if (vpar == 0)
		{
			// assign masked probes to left sub-region
			horizontalDivide (l_region, hstep, vstep, hpar, 0, start,
								probe_div - 1);

			// assign unmasked probes to right sub-region
			horizontalDivide (r_region, hstep, vstep, hpar, 1, probe_div, end);
		}
		else
		{
			// assign masked probes to right sub-region
			horizontalDivide (r_region, hstep, vstep, hpar, 1, start,
								probe_div - 1);

			// assign unmasked probes to left sub-region
			horizontalDivide (l_region, hstep, vstep, hpar, 0, probe_div, end);
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
		return "2DCentralPartitioning-" + stop_dim + "-" + filler;
	}
}
