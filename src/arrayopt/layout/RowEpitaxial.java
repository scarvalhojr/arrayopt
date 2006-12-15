/*
 * RowEpitaxial.java
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
 * Row-epitaxial algorithm due to Kahng et al. This is a free implementation of
 * the algorithm described in the paper:
 * 
 * <P>"Engineering a Scalable Placement Heuristic for DNA Probe Arrays",
 * A. Kahng, I. Mandoiu, P. Pevzer, S. Reda, and A. Zelikovsky. In Proc. 7th
 * Annual Int. Conference on Research in Computational Molecular Biology
 * (RECOMB), 2003, pp. 148–156.</P>
 * 
 * <P>This algorithm is similar to the {@linkplain GreedyPlacer} implementation.
 * The main difference is that the row-epitaxial is restricted to border length
 * minimization, while the {@linkplain GreedyPlacer} is also able produce a
 * layout reducing the overall conflict index (see {@link ConflictIndex}).</P>
 * 
 * <P>The row-epitaxial produces an initial layout placing the probes
 * sequentially on the spots in whatever order they appear. (Alternatively, it
 * is possible to force a randomization of the input; see
 * {@link #RowEpitaxial(int, boolean)}). Then, it scans the chip
 * top-to-bottom, left-to-right, and, for every non-empty spot <CODE>s</CODE>
 * with a probe <CODE>p</CODE>, it finds a probe <CODE>q</CODE> with minimum
 * cost. The cost is defined as the sum of the Hamming distances to the probes
 * placed on the top and left neighbors of <CODE>s</CODE>. Probe <CODE>q</CODE>
 * is searched on the next <CODE>n</CODE> spots of <CODE>s</CODE>
 * (<CODE>n</CODE> is the look-ahead parameter that can be specified at
 * instantiation time; default value is defined by
 * {@link #DEFAULT_LOOK_AHEAD}).</P>  
 * 
 * <P>Note that, if the region to be filled is completely empty, the layout
 * produced by the row-epitaxial algorithm will be identical to the one produced
 * by the {@linkplain GreedyPlacer} algorithm with border length minimization.
 * If some spots are non-empty, the results are likely to be different. This is
 * because the row-epitaxial only considers the top and left neighbors of a
 * spot, while {@linkplain GreedyPlacer} also considers the bottom and right
 * neighbors (if they are not empty, that is). As a result, the
 * {@linkplain GreedyPlacer} is slightly slower but it is also likely to produce
 * marginally better results in such cases.</P>
 * 
 * @author Sergio A. de Carvalho Jr.
 */
public class RowEpitaxial implements LayoutAlgorithm, FillingAlgorithm
{
	public static final int DEFAULT_LOOK_AHEAD = 20000;
	
	private int look_ahead;
	
	private boolean randomize_first;
	
	private RectangularRegion chip_region;

	public RowEpitaxial ()
	{
		this(DEFAULT_LOOK_AHEAD);
	}

	public RowEpitaxial (int look_ahead)
	{
		this(look_ahead, false);
	}

	public RowEpitaxial (int look_ahead, boolean randomize_first)
	{
		this.look_ahead = look_ahead > 0 ? look_ahead : DEFAULT_LOOK_AHEAD;
		this.randomize_first = randomize_first;
	}

	/**
	 *
	 */
	public void changeLayout (Chip chip)
	{
		int		id[];

		// reset current layout (if any)
		chip.resetLayout();

		// get list of movable probes
		id = chip.getMovableProbes ();

		fillRegion (chip, chip.getChipRegion(), id, 0, id.length - 1);
	}

	/**
	 *
	 */
	public int fillRegion (Chip chip, Region region, int probe_id[])
	{
		return fillRegion (chip, region, probe_id, 0, probe_id.length - 1);
	}

	/**
	 *
	 */
	public int fillRegion (Chip chip, Region region, int probe_id[], int start,
		int end)
	{
		RectangularRegion r;

		if (!(region instanceof RectangularRegion))
			throw new IllegalArgumentException
				("Only rectangular regions are supported.");

		r = (RectangularRegion) region;
		
		this.chip_region = chip.getChipRegion();

		if (chip instanceof SimpleChip)
			return fillRegion ((SimpleChip) chip, r, probe_id, start, end);

		else if (chip instanceof AffymetrixChip)
			return fillRegion ((AffymetrixChip) chip, r, probe_id, start, end);

		else
			throw new IllegalArgumentException ("Unsupported chip type.");
	}

	/**
	 *
	 */
	private int fillRegion (SimpleChip chip, RectangularRegion region,
		int probe_id[], int start, int end)
	{
		int id, unplaced;
		
		unplaced = initialPlacement (chip, region, probe_id, start, end);

		for (int r = region.first_row; r <= region.last_row; r ++)
		{
			for (int c = region.first_col; c <= region.last_col; c++)
			{
				// skip spot is empty
				if ((id = chip.spot[r][c]) == Chip.EMPTY_SPOT)
					continue;

				// skip fixed spot
				if (chip.isFixedSpot(r, c))
					continue;

				// replaced current probe with probe resulting in minimum cost
				replace (chip, region, r, c, id);
			}
		}

		return unplaced;
	}

	/**
	 *
	 */
	private int fillRegion (AffymetrixChip chip, RectangularRegion region,
		int probe_id[], int start, int end)
	{
		int unplaced;
		
		unplaced = initialPlacement (chip, region, probe_id, start, end);

		// TODO implement this

		return unplaced;
	}

	private int initialPlacement (SimpleChip chip, RectangularRegion region,
			int probe_id[], int start, int end)
	{
		if (randomize_first) randomizeProbes (probe_id, start, end);
		
		for (int r = region.first_row; r <= region.last_row; r++)
			for (int c = region.first_col; c <= region.last_col; c++)
			{
				// skip fixed spots
				if (chip.isFixedSpot(r, c))
					continue;
				
				// spot must be empty 
				if (chip.spot[r][c] != Chip.EMPTY_SPOT)
					continue;

				// place probe
				chip.spot[r][c] = probe_id[start++];
				
				// check if all probes have been placed
				if (start > end)
					// yes: zero unplaced probes
					return 0;
			}
		
		return (end - start + 1);
	}

	private int initialPlacement (AffymetrixChip chip,
			RectangularRegion region, int probe_id[], int start, int end)
	{
		if (randomize_first) randomizeProbes (probe_id, start, end);
		
		// TODO implement this
		
		return 0;
	}
	
	private void replace (SimpleChip chip, RectangularRegion region, int row,
			int col, int curr_id)
	{
		long	cost, min;
		int		top, left, id, best_row, best_col;
		
		// get ID of probe placed on the top spot
		if (row > chip_region.first_row)
			top = chip.spot[row - 1][col];
		else
			top = Chip.EMPTY_SPOT;

		// get ID of probe placed on the left spot
		if (col > chip_region.first_col)
			left = chip.spot[row][col - 1];
		else
			left = Chip.EMPTY_SPOT;

		best_row = row;
		best_col = col;
		
		// compute distance of current probe to top and left spots 
		min = 0;
		if (top != Chip.EMPTY_SPOT)
			min += LayoutEvaluation.hammingDistance(chip, curr_id, top);
		if (left != Chip.EMPTY_SPOT)
			min += LayoutEvaluation.hammingDistance(chip, curr_id, left);
		
		int r = row;
		int c = col;
		
		// check probes placed on the next 'look_ahead' spots
		for (int count = 0; count <= look_ahead; count++)
		{
			// next spot
			if (++c > region.last_col)
			{
				if (++r > region.last_row)
					break;
				
				c = region.first_col;
			}
			
			// skip fixed spots
			if (chip.isFixedSpot(r, c))
				continue;
			
			// get ID of candidate probe (skip if spot is empty) 
			if ((id = chip.spot[r][c]) == Chip.EMPTY_SPOT)
				continue;
						
			// compute cost of candidate probe 
			if (top != Chip.EMPTY_SPOT)
				cost = LayoutEvaluation.hammingDistance(chip, id, top);
			else
				cost = 0;
			if (left != Chip.EMPTY_SPOT)
				cost += LayoutEvaluation.hammingDistance(chip, id, left);
			
			// check if found better option
			if (cost < min)
			{
				min = cost;
				best_row = r;
				best_col = c;				
			}
		}
		
		// swap current probe with best candidate
		chip.spot[row][col] = chip.spot[best_row][best_col];
		chip.spot[best_row][best_col] = curr_id;
	}

	private void randomizeProbes (int probe_id[], int start, int end)
	{
		int tmp, rand;
		
		for (int i = start; i < end; i++)
		{
			// select a random element of the list
			rand = i + (int) ((end - i + 1) * Math.random());

			// swap the selected element with
			// the one at the current position
			tmp = probe_id[i];
			probe_id[i] = probe_id[rand];
			probe_id[rand] = tmp;
		}
	}
}
