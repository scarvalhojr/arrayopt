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
 * This class implements the Row-epitaxial algorithm. This is a free
 * implementation of the algorithm described in the paper:
 * 
 * <P>Kahng, A.B.; Mandoiu, I.; Pevzner, P.; Reda, S. & Zelikovsky, A.:
 * Engineering a scalable placement heuristic for DNA probe arrays, Proceedings
 * of the seventh annual international conference on research in computational
 * molecular biology (RECOMB), ACM Press, 2003, 148-156.</P>
 *
 * <P>The algorithm was initially developed for border length minimization, but
 * this implementaion can work with two types of conflict minimization:
 * border length ({@link #BORDER_LENGTH_MIN}) or conflict index
 * ({@link #CONFLICT_INDEX_MIN}) minimization.</P>
 * 
 * <P>The row-epitaxial improves an initial layout by re-locating probes. Spots
 * are re-filled sequentially, from top to bottom, left to right. For every
 * spot <CODE>s</CODE> with a probe <CODE>p</CODE>, it finds a probe
 * <CODE>q</CODE> with minimum cost to fill <CODE>s</CODE>. The cost is either
 * the sum of border conflicts or sum of conflict indices with the re-filled
 * neighbors. Probe <CODE>q</CODE> is searched on the next <CODE>n</CODE> spots
 * of <CODE>s</CODE>, where <CODE>n</CODE> is the look-ahead parameter that must
 * be specified at instantiation time.</P>
 * 
 * @author Sergio A. de Carvalho Jr.
 */
public class RowEpitaxial implements LayoutAlgorithm
{
	/**
	 * Constant used to indicate that the algorithm should try to minimize the
	 * total border length of the microarray layout.
	 */
	public static final int BORDER_LENGTH_MIN = 0;
	
	/**
	 * Constant used to indicate that the algorithm should try to minimize the
	 * sum of conflict indices of the microarray layout.
	 */
	public static final int CONFLICT_INDEX_MIN = 1;

	/**
	 * Maximum number of probe candidades considered for each spot.
	 */
	private int look_ahead;
	
	/**
	 * This variable stores the current minimization mode used by the algorithm.
	 * Possible values are {@link #BORDER_LENGTH_MIN} and
	 * {@link #CONFLICT_INDEX_MIN}. 
	 */
	private int mode;
	
	private RectangularRegion chip_region;
	
	private int embed_len;
	
	private int probe_len;
	
	private int ci_dim;
	
	private double pos_weight[];
	
	private double m_cost[];
	
	private double u_cost[];
	
	/**
	 * Creates an instance of the Row-epitaxial algorithm with the desired
	 * minimization mode and look-ahead value.
	 * 
	 * @see #BORDER_LENGTH_MIN
	 * @see #CONFLICT_INDEX_MIN
	 * @see #look_ahead
	 * @param mode minimization mode, either border length or conflict index
	 * @param look_ahead maximum number of candidades examined for each spot
	 */
	public RowEpitaxial (int mode, int look_ahead)
	{
		switch (mode)
		{
			case BORDER_LENGTH_MIN:
			case CONFLICT_INDEX_MIN:
				this.mode = mode;
				break;
				
			default:
				throw new IllegalArgumentException
					("Invalid minimization mode");
		}

		if (look_ahead < 0)
			throw new IllegalArgumentException ("Invalid look-ahead: " +
					look_ahead);
		
		this.look_ahead = look_ahead;
	}

	/**
	 * Attempts to improve the layout of a microarray chip using the
	 * Row-epitaxial algorithm.
	 * 
	 * @param chip instance of a microarray chip
	 */
	public void changeLayout (Chip chip)
	{
		this.chip_region = chip.getChipRegion();

		if (chip instanceof SimpleChip)
		{
			changeLayout ((SimpleChip) chip);
		}
		else
			throw new IllegalArgumentException ("Unsupported chip type.");
	}

	private void changeLayout (SimpleChip chip)
	{
		int id;
		
		if (mode == CONFLICT_INDEX_MIN)
		{
			// prepare for faster conflict index calculations
			ci_dim = ConflictIndex.dimConflictRegion();
			embed_len = chip.getEmbeddingLength();
			probe_len = chip.getProbeLength();
			
			if (m_cost == null)
			{
				m_cost = new double [embed_len];
				u_cost = new double [embed_len];
				pos_weight = new double [probe_len + 1];
			}
			else if (m_cost.length != embed_len)
			{
				m_cost = new double [embed_len];
				u_cost = new double [embed_len];
				pos_weight = new double [probe_len + 1];
			}
			
			// store position weights locally
			for (int b = 0; b <= probe_len; b++)
				pos_weight[b] = ConflictIndex.positionWeight(b, probe_len);
		}
		
		for (int r = chip_region.first_row; r <= chip_region.last_row; r ++)
		{
			for (int c = chip_region.first_col; c <= chip_region.last_col; c++)
			{
				// skip fixed spot
				if (chip.isFixedSpot(r, c))
					continue;

				// TODO should empty spots be really skipped?
				
				// skip if spot is empty
				if ((id = chip.spot[r][c]) == Chip.EMPTY_SPOT)
					continue;

				// replaced current probe with probe resulting in minimum cost
				if (mode == BORDER_LENGTH_MIN)
					minBorderLength (chip, r, c, id);
				else // if (mode == CONFLICT_INDEX_MIN)
					minConflictIndex (chip, r, c, id);
			}
		}
	}

	private void minBorderLength (SimpleChip chip, int row, int col, int curr)
	{
		long	cost, min;
		int		r, c, id, best_row, best_col, top, left;
		
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
		
		// if region around the spot is empty,
		// there is no criteria to choose a replacement
		if (top == Chip.EMPTY_SPOT && left == Chip.EMPTY_SPOT) return;

		// compute border conflicts of current probe with top and left neighbors 
		min = 0;
		if (top != Chip.EMPTY_SPOT)
			min += LayoutEvaluation.hammingDistance(chip, curr, top);
		if (left != Chip.EMPTY_SPOT)
			min += LayoutEvaluation.hammingDistance(chip, curr, left);

		best_row = r = row;
		best_col = c = col;
		
		// check probes placed on the next 'look_ahead' spots
		for (int count = 0; count <= look_ahead;)
		{
			// next spot
			if (++c > chip_region.last_col)
			{
				if (++r > chip_region.last_row)
					break;
				
				c = chip_region.first_col;
			}
			
			// skip fixed spots
			if (chip.isFixedSpot(r, c))
				continue;
			
			// get ID of candidate probe (skip if spot is empty) 
			if ((id = chip.spot[r][c]) == Chip.EMPTY_SPOT)
				continue;
			
			count++;
						
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
		chip.spot[best_row][best_col] = curr;
	}

	private void minConflictIndex (SimpleChip chip, int row, int col, int curr)
	{
		boolean	empty;
		double	cost, min;
		int		r, c, id, best_row, best_col;
		
		// prepare conflict index costs
		empty = examineNeighbors (chip, row, col);
		
		// if region around the spot is empty,
		// there is no criteria to choose a replacement
		if (empty) return;

		// compute current conflict index
		min = conflictIndex (chip, curr, Double.POSITIVE_INFINITY);
		best_row = r = row;
		best_col = c = col;
		
		// check probes placed on the next 'look_ahead' spots
		for (int count = 0; count <= look_ahead;)
		{
			// next spot
			if (++c > chip_region.last_col)
			{
				if (++r > chip_region.last_row)
					break;
				
				c = chip_region.first_col;
			}
			
			// skip fixed spots
			if (chip.isFixedSpot(r, c))
				continue;
			
			// get ID of candidate probe (skip if spot is empty) 
			if ((id = chip.spot[r][c]) == Chip.EMPTY_SPOT)
				continue;
			
			count++;
						
			// compute cost of candidate probe
			cost = conflictIndex (chip, id, min);
			
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
		chip.spot[best_row][best_col] = curr;
	}
	
	private boolean examineNeighbors (SimpleChip chip, int row, int col)
	{
		boolean empty = true;
		double delta;
		int r, c, id, base, step, word, bitmask = 0;
		
		// reset costs
		for (step = 0; step < embed_len; step++)
			m_cost[step] = u_cost[step] = 0;
		
		// define region around the spot that needs to be examined
		// (only look at spots to the left or above)
		int min_row = Math.max(row - ci_dim, chip_region.first_row);
		int max_row = row;
		int min_col = Math.max(col - ci_dim, chip_region.first_col);
		int max_col = Math.min(col + ci_dim, chip_region.last_col);
		
		for (r = min_row; r <= max_row; r++)
		{
			if (r == max_row) max_col = col - 1;
			
			for (c = min_col; c <= max_col; c++)
			{
				// skip if neighbor is empty
				if ((id = chip.spot[r][c]) == Chip.EMPTY_SPOT)
					continue;
				
				// get distance-dependent weight
				delta = ConflictIndex.distanceWeight(r, c, row, col);
				
				if (delta == 0)
					continue;
				
				empty = false;
				
				for (base = 0, step = 0, word = - 1; step < embed_len; step++)
				{
					if (step % Integer.SIZE == 0)
					{
						word++;
						bitmask = 0x01 << (Integer.SIZE - 1);
					}
					else
						bitmask >>>= 1;

					// check state of embedding at current step
					if ((chip.embed[id][word] & bitmask) != 0)
					{
						// spot is in an unmasked step
						m_cost[step] += delta;
						
						base++;
					}
					else
					{
						// spot is in a masked step
						u_cost[step] += pos_weight[base] * delta;
					}
				}
			}
		}

		return empty;
	}
	
	private double conflictIndex (SimpleChip chip, int id, double max)
	{
		int base, step, word, bitmask = 0;
		double ci = 0;
		
		for (base = 0, step = 0, word = - 1; step < embed_len; step++)
		{
			if (step % Integer.SIZE == 0)
			{
				word++;
				bitmask = 0x01 << (Integer.SIZE - 1);
			}
			else
				bitmask >>>= 1;
			
			// check embedding state at current step 
			if ((chip.embed[id][word] & bitmask) != 0)
			{
				// spot is in an unmasked step
				ci += u_cost[step];

				// increment the number of synthesized bases
				base++;
			}
			else
			{
				// spot is in a masked step
				ci += pos_weight[base] * m_cost[step];
				
				// stop if CI exceeds limit
				if (ci > max) break;
			}
		}

		return ci;
	}

	/**
	 * Returns the algorithm's name together with current options.
	 * 
	 * @return algorithm's name and configurable options
	 */
	@Override
	public String toString ()
	{
		String m;
		
		switch (this.mode)
		{
			case BORDER_LENGTH_MIN:
				m = "-BL-";
				break;
				
			case CONFLICT_INDEX_MIN:
				m = "-CI-";
				break;
			
			default:
				m = "-?";
				break;
		}
		
		return this.getClass().getSimpleName() + m + look_ahead;
	}
}
