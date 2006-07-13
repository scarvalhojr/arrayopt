/*
 * TwoDimensionalCenteredPartitioning.java
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
 * �rrayOpt is distributed in the hope that it will be useful, but WITHOUT ANY
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
public class TwoDimensionalCenteredPartitioning implements PlacementAlgorithm
{
	public static final int DEFAULT_STOP_DIMENSION = 4;

	private FillingAlgorithm filler;

	private int stop_dimension;

	private double MIN_DIV_RATE = .1;

	public TwoDimensionalCenteredPartitioning (FillingAlgorithm filler)
	{
		this(filler, DEFAULT_STOP_DIMENSION);
	}

	public TwoDimensionalCenteredPartitioning (FillingAlgorithm filler, int stop_dimension)
	{
		this.filler = filler;
		this.stop_dimension = (stop_dimension < 1) ? 1 : stop_dimension;
	}

	/**
	 *
	 */
	public int makeLayout (Chip chip)
	{
		int		id[], m;
		int		rows_per_probe;

		if (chip instanceof AffymetrixChip)
			rows_per_probe = 2;
		else if (chip instanceof SimpleChip)
			rows_per_probe = 1;
		else
			throw new IllegalArgumentException ("Unsupported chip type.");

		// reset current layout (if any)
		chip.resetLayout();

		// get list of movable probes
		id = chip.getMovableProbes ();
		
		m = chip.getEmbeddingLength() / 2;

		return horizontalDivide (chip, m -1, m, chip.getChipRegion(),
				rows_per_probe,	0, 0, id, 0, id.length - 1);
	}

	protected int horizontalDivide (Chip chip, int vstep, int hstep,
			RectangularRegion r, int rows_per_probe, int hor_par, int ver_par,
			int probe_id[],	int start, int end)
	{
		RectangularRegion	t_region, b_region;
		int					row_div, probe_div, m_rows, u_rows, qt_cols, overflow, unplaced;
		double				div_rate;

		if (end - start + 1 < 2)
		{
			// insufficient number of probes for partitioning:
			// place probes on the specified region and return
			return fillRegion (chip, r, probe_id, start, end);
		}

		if (hstep < 0)
		{
			// cannot partition anymore:
			// place probes on the specified region and return
			return fillRegion (chip, r, probe_id, start, end);
		}

		if (r.last_row - r.first_row + 1 <= stop_dimension * rows_per_probe)
		{
			if (r.last_col - r.first_col + 1 <= stop_dimension)
			{
				// region cannot be partitioned anymore:
				// place probes on the specified region and return
				return fillRegion (chip, r, probe_id, start, end);
			}

			// region can still be vertically partitined
			return verticalDivide (chip, hstep, vstep, r, rows_per_probe,
					hor_par, ver_par, probe_id, start, end);
		}

		// split the probes into two groups, M and U, according to
		// whether they are masked or unmasked at this step
		probe_div = divideProbes (chip, hstep, probe_id, start, end);

		div_rate = ((double)(probe_div - start)) / (end - start + 1);

		if (div_rate < MIN_DIV_RATE || 1 - div_rate < MIN_DIV_RATE)
		{
			return horizontalDivide (chip, hstep - 1, vstep, r, rows_per_probe,
									hor_par, ver_par, probe_id, start, end);
		}

		qt_cols = r.last_col - r.first_col + 1;

		// compute how many rows will be assigned to masked and unmasked probes
		m_rows = (int) Math.round(div_rate * (r.last_row - r.first_row + 1));

		// make sure that both partitions get a minimum number of rows
		// (this might be useless depending on the minimum division rate)
		if (m_rows == 0)
			m_rows += rows_per_probe;
		else if (m_rows == r.last_row - r.first_row + 1)
			m_rows -= rows_per_probe;

		// make sure the partitions will get an appropriate number of rows
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

		// check if each region will be able to handle the number of probes
		if ((overflow = (probe_div - start) - (m_rows * qt_cols) / rows_per_probe) > 0)
		{
			probe_div -= overflow;
		}
		else if ((overflow = (end - probe_div + 1) - (u_rows * qt_cols) / rows_per_probe) > 0)
		{
			probe_div += overflow;
		}

		// divide the current region into two _horizontal_ areas, T (top) and B (bottom),
		// proportionally to the number of masked and unmasked probes
		if (hor_par == 0)
		{
			// masked probes will be assigned to region T
			row_div = r.first_row + m_rows;
		}
		else
		{
			// unmasked probes will be assigned to region T
			row_div = r.first_row + u_rows;
		}

		// create T (top) and B (bottom) regions
		t_region = new RectangularRegion (r.first_row, row_div - 1, r.first_col, r.last_col);
		b_region = new RectangularRegion (row_div, r.last_row, r.first_col, r.last_col);

		hstep--;

		// check horizontal parity
		if (hor_par == 0)
		{
			// assign masked probes to region T
			unplaced = verticalDivide (chip, hstep, vstep, t_region,
					rows_per_probe,	0, ver_par, probe_id, start, probe_div - 1);

			// assign unmasked probes to region B
			unplaced += verticalDivide (chip, hstep, vstep, b_region,
					rows_per_probe, 1, ver_par, probe_id, probe_div, end);
		}
		else
		{
			// assign masked probes to region B
			unplaced = verticalDivide (chip, hstep, vstep, b_region,
					rows_per_probe, 1, ver_par, probe_id, start, probe_div - 1);

			// assign unmasked probes to region T
			unplaced += verticalDivide (chip, hstep, vstep, t_region,
					rows_per_probe, 0, ver_par, probe_id, probe_div, end);
		}

		// return number of unplaced probes
		return unplaced;
	}

	protected int verticalDivide (Chip chip, int hstep, int vstep,
			RectangularRegion r, int rows_per_probe, int hor_par, int ver_par,
			int probe_id[], int start, int end)
	{
		RectangularRegion	l_region, r_region;
		int					col_div, probe_div, m_cols, u_cols, qt_rows, overflow, unplaced;
		double				div_rate;

		if (end - start + 1 < 2)
		{
			// insufficient number of probes for partitioning:
			// place probes on the specified region and return
			return fillRegion (chip, r, probe_id, start, end);
		}

		if (vstep >= chip.getEmbeddingLength())
		{
			// cannot partition anymore:
			// place probes on the specified region and return
			return fillRegion (chip, r, probe_id, start, end);
		}

		if (r.last_col - r.first_col + 1 <= stop_dimension)
		{
			if (r.last_row - r.first_row + 1 <= stop_dimension * rows_per_probe)
			{
				// region cannot be partitioned anymore:
				// place probes on the specified region and return
				return fillRegion (chip, r, probe_id, start, end);
			}

			// region can still be horizontally partitined
			return horizontalDivide (chip, hstep, vstep, r, rows_per_probe,
								hor_par, ver_par, probe_id, start, end);
		}

		// split the probes into two groups, M and U, according to
		// whether they are masked or unmasked at this step
		probe_div = divideProbes (chip, vstep, probe_id, start, end);

		div_rate = ((double)(probe_div - start)) / (end - start + 1);

		if (div_rate < MIN_DIV_RATE || 1 - div_rate < MIN_DIV_RATE)
		{
			return verticalDivide (chip, hstep, vstep + 1, r, rows_per_probe,
					hor_par, ver_par, probe_id, start, end);
		}

		qt_rows = r.last_row - r.first_row + 1;

		// compute how many columns are needed to the masked probes
		m_cols = (int) Math.round(div_rate * (r.last_col - r.first_col + 1));

		// make sure that both partitions get a minimum number of columns
		// (this might be useless depending on the minimum division rate)
		if (m_cols == 0)
			m_cols += rows_per_probe;
		else if (m_cols == r.last_col - r.first_col + 1)
			m_cols -= rows_per_probe;

		u_cols = (r.last_col - r.first_col + 1) - m_cols;

		// check if each region will be able to handle the number of probes
		if ((overflow = (probe_div - start) - (m_cols * qt_rows) / rows_per_probe) > 0)
		{
			probe_div -= overflow;
		}
		else if ((overflow = (end - probe_div + 1) - (u_cols * qt_rows) / rows_per_probe) > 0)
		{
			probe_div += overflow;
		}

		// divide the current region into two _vertical_ areas, L (left) and R (right),
		// proportionally to the number of masked and unmasked probes
		if (ver_par == 0)
		{
			// masked probes will be assigned to region L
			col_div = r.first_col + m_cols;
		}
		else
		{
			// unmasked probes will be assigned to region L
			col_div = r.first_col + u_cols;
		}

		// create L (left) and R (right) regions
		l_region = new RectangularRegion (r.first_row, r.last_row, r.first_col, col_div - 1);
		r_region = new RectangularRegion (r.first_row, r.last_row, col_div, r.last_col);

		vstep++;

		// check vertical parity
		if (ver_par == 0)
		{
			// assign masked probes to region L
			unplaced = horizontalDivide (chip, hstep, vstep, l_region,
					rows_per_probe, hor_par, 0, probe_id, start, probe_div - 1);

			// assign unmasked probes to region R
			unplaced += horizontalDivide (chip, hstep, vstep, r_region,
					rows_per_probe, hor_par, 1, probe_id, probe_div, end);
		}
		else
		{
			// assign masked probes to region R
			unplaced = horizontalDivide (chip, hstep, vstep, r_region,
					rows_per_probe, hor_par, 1, probe_id, start, probe_div - 1);

			// assign unmasked probes to region L
			unplaced += horizontalDivide (chip, hstep, vstep, l_region,
					rows_per_probe, hor_par, 0, probe_id, probe_div, end);
		}

		// return number of unplaced probes
		return unplaced;
	}

	protected int divideProbes (Chip chip, int step, int probe_id[],
		int start, int end)
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
	
	private OptimumSingleProbeEmbedding ospe;
	
	private int fillRegion (Chip chip, RectangularRegion r, int probe_id[],
			int start, int end)
	{
		/*
		// TODO do this nicely!
		if (ospe == null)
			// TODO make mode configurable!
			ospe = OptimumSingleProbeEmbedding.createEmbedder(chip,
					OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN);
		
		if (end > start)
		{
			ospe.reembedProbe(probe_id[start+1], probe_id[start]);
			for (int i = start + 2; i <= end; i++)
				ospe.reembedProbe(probe_id[i]);
		}
		//*/
		
		return filler.fillRegion (chip, r, probe_id, start, end);
	}
}