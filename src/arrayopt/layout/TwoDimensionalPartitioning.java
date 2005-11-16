/*
 * TwoDimensionalPartitioning.java
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
public class TwoDimensionalPartitioning extends PartitioningAlgorithm
{
	public static final int DEFAULT_STOP_DIMENSION = 4;

	private FillingAlgorithm filler;

	private int stop_dimension;

	private double MIN_DIV_RATE = .1;

	public TwoDimensionalPartitioning (FillingAlgorithm filler)
	{
		this(filler, DEFAULT_STOP_DIMENSION);
	}

	public TwoDimensionalPartitioning (FillingAlgorithm filler, int stop_dimension)
	{
		this.filler = filler;
		this.stop_dimension = (stop_dimension < 1) ? 1 : stop_dimension;
	}

	public int placeProbes (Chip chip, int first_probe, int last_probe)
	{
		int rows_per_probe;

		if (chip instanceof AffymetrixChip)
			rows_per_probe = 2;
		else
			rows_per_probe = 1;

		return horizontalDivide (chip, 0, chip.getChipRegion(), rows_per_probe, 0, 0, first_probe, last_probe);
	}

	protected int horizontalDivide (Chip chip, int step, RectangularRegion r, int rows_per_probe, int hor_par, int ver_par, int first_probe, int last_probe)
	{
		RectangularRegion	t_region, b_region;
		int					row_div, probe_div, m_rows, u_rows, qt_cols, overflow, unplaced;
		double				div_rate;

		if (last_probe - first_probe + 1 < 2)
		{
			// insufficient number of probes for partitioning:
			// place probes on the specified region and return
			return filler.fillRegion (chip, r, first_probe, last_probe);
		}

		if (step >= chip.getEmbeddingLength())
		{
			// cannot partition anymore:
			// place probes on the specified region and return
			return filler.fillRegion (chip, r, first_probe, last_probe);
		}

		if (r.last_row - r.first_row + 1 <= stop_dimension * rows_per_probe)
		{
			if (r.last_col - r.first_col + 1 <= stop_dimension)
			{
				// region cannot be partitioned anymore:
				// place probes on the specified region and return
				return filler.fillRegion (chip, r, first_probe, last_probe);
			}
			else
			{
				// region can still be vertically partitined
				return verticalDivide (chip, step, r, rows_per_probe, hor_par, ver_par, first_probe, last_probe);
			}
		}

		// split the probes into two groups, M and U, according to
		// whether they are masked or unmasked at this step
		probe_div = divideProbes (chip, step, first_probe, last_probe);

		div_rate = ((double)(probe_div - first_probe)) / (last_probe - first_probe + 1);

		if (div_rate < MIN_DIV_RATE || 1 - div_rate < MIN_DIV_RATE)
		{
			return horizontalDivide (chip, step + 1, r, rows_per_probe, hor_par, ver_par, first_probe, last_probe);
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
		if ((overflow = (probe_div - first_probe) - (m_rows * qt_cols) / rows_per_probe) > 0)
		{
			probe_div -= overflow;
		}
		else if ((overflow = (last_probe - probe_div + 1) - (u_rows * qt_cols) / rows_per_probe) > 0)
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

		// check horizontal parity
		if (hor_par == 0)
		{
			// assign masked probes to region T
			unplaced = verticalDivide (chip, step + 1, t_region, rows_per_probe, 0, ver_par, first_probe, probe_div - 1);

			// assign unmasked probes to region B
			unplaced += verticalDivide (chip, step + 1, b_region, rows_per_probe, 1, ver_par, probe_div, last_probe);
		}
		else
		{
			// assign masked probes to region B
			unplaced = verticalDivide (chip, step + 1, b_region, rows_per_probe, 1, ver_par, first_probe, probe_div - 1);

			// assign unmasked probes to region T
			unplaced += verticalDivide (chip, step + 1, t_region, rows_per_probe, 0, ver_par, probe_div, last_probe);
		}

		// return number of unplaced probes
		return unplaced;
	}

	protected int verticalDivide (Chip chip, int step, RectangularRegion r, int rows_per_probe, int hor_par, int ver_par, int first_probe, int last_probe)
	{
		RectangularRegion	l_region, r_region;
		int					col_div, probe_div, m_cols, u_cols, qt_rows, overflow, unplaced;
		double				div_rate;

		if (last_probe - first_probe + 1 < 2)
		{
			// insufficient number of probes for partitioning:
			// place probes on the specified region and return
			return filler.fillRegion (chip, r, first_probe, last_probe);
		}

		if (step >= chip.getEmbeddingLength())
		{
			// cannot partition anymore:
			// place probes on the specified region and return
			return filler.fillRegion (chip, r, first_probe, last_probe);
		}

		if (r.last_col - r.first_col + 1 <= stop_dimension)
		{
			if (r.last_row - r.first_row + 1 <= stop_dimension * rows_per_probe)
			{
				// region cannot be partitioned anymore:
				// place probes on the specified region and return
				return filler.fillRegion (chip, r, first_probe, last_probe);
			}
			else
			{
				// region can still be horizontally partitined
				return horizontalDivide (chip, step, r, rows_per_probe, hor_par, ver_par, first_probe, last_probe);
			}
		}

		// split the probes into two groups, M and U, according to
		// whether they are masked or unmasked at this step
		probe_div = divideProbes (chip, step, first_probe, last_probe);

		div_rate = ((double)(probe_div - first_probe)) / (last_probe - first_probe + 1);

		if (div_rate < MIN_DIV_RATE || 1 - div_rate < MIN_DIV_RATE)
		{
			return verticalDivide (chip, step + 1, r, rows_per_probe, hor_par, ver_par, first_probe, last_probe);
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
		if ((overflow = (probe_div - first_probe) - (m_cols * qt_rows) / rows_per_probe) > 0)
		{
			probe_div -= overflow;
		}
		else if ((overflow = (last_probe - probe_div + 1) - (u_cols * qt_rows) / rows_per_probe) > 0)
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

		// check vertical parity
		if (ver_par == 0)
		{
			// assign masked probes to region L
			unplaced = horizontalDivide (chip, step + 1, l_region, rows_per_probe, hor_par, 0, first_probe, probe_div - 1);

			// assign unmasked probes to region R
			unplaced = horizontalDivide (chip, step + 1, r_region, rows_per_probe, hor_par, 1, probe_div, last_probe);
		}
		else
		{
			// assign masked probes to region R
			unplaced = horizontalDivide (chip, step + 1, r_region, rows_per_probe, hor_par, 1, first_probe, probe_div - 1);

			// assign unmasked probes to region L
			unplaced = horizontalDivide (chip, step + 1, l_region, rows_per_probe, hor_par, 0, probe_div, last_probe);
		}

		// return number of unplaced probes
		return unplaced;
	}

	protected int divideProbes (Chip chip, int step, int first_probe, int last_probe)
	{
		int	probe_id, w, mask, backup_min, backup_max;

		backup_min = first_probe;
		backup_max = last_probe;

		// which 4-byte word will be interrogated?
		w = (int) Math.floor((double) step / Integer.SIZE);

		// prepare mask to interrogate corresponding bit
		mask = 0x01 << (Integer.SIZE - 1 - (step - w * Integer.SIZE));

		while (first_probe <= last_probe)
		{
			probe_id = chip.probe_list[first_probe];

			if ((chip.embed[probe_id][w] & mask) == 0)
			{
				// probe is masked at this step:
				// proceed to next embedding on the list
				first_probe ++;
			}
			else
			{
				// probe is unmasked at this step:
				// exchange it with the last element on the list
				chip.probe_list[first_probe] = chip.probe_list[last_probe];
				chip.probe_list[last_probe] = probe_id;
				last_probe --;
			}
		}

		return first_probe;
	}
}
