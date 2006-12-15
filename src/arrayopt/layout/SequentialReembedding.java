/*
 * SequentialReembedding.java
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
 * This class implements the Sequential re-embedding post-placement
 * optimization. The alorithm is described in:
 * 
 * <P>"Evaluation of placement techniques for DNA probe array", A. Kahng,
 * I. Mandoiu, S. Reda, X. Xu, and A. Zelikovsky. In Proc. IEEE/ACM
 * Int. Conference on Computer-Aided Design (ICCAD 2003), pp. 262–269.</P>
 * 
 * <P>The algorithm scans the chip top to bottom, left to right, and, for each
 * spot, optimally re-embeds its probe in regards to its neighbors using the
 * Optimum Single-Probe Embedding (OSPE) algorithm (implemented by
 * {@linkplain OptimumSingleProbeEmbedding}). The desired conflict minimization
 * function (whether border length or conflict index minimization) must be
 * specified at instantiation time (to the constructor).<P>
 * 
 * <P>The algorithm actually scans the chip several times and, since the OSPE
 * never increases the amount of conflicts (sum of border lengths or conflict
 * indices), it is guaranteed to never produce a worse embedding. However, this
 * also means that the improvements decrease with each iteration until no more
 * improvements are possible. In other words, the algorithm finds a local
 * optimum solution.</P>
 * 
 * <P>Since the improvements always decrease with each iteration, the algorithm
 * can be configured to stop once the reduction of conflicts drops below a given
 * threshold. The threshold is set to 0.1% by default.</P>
 * 
 * <P>When computing an optimal embedding of a probe, the algorithm considers
 * all of its neighbors. However, in the first pass, it is possible to configure
 * it to only consider those probes that have already been re-embedded. This
 * feature is called "reset first" and my lead to better solutions depending on
 * the current embeddings of the probes.<P>
 * 
 * @author Sergio A. de Carvalho Jr.
 */
public class SequentialReembedding implements LayoutAlgorithm
{
	public static final double DEFAULT_THRESHOLD = 0.001d;
	
	private OptimumSingleProbeEmbedding embedder;
	
	private int mode;
	
	private boolean reset_first;
	
	private double threshold;
	
	private int spot_copy[][];
	
	private int num_rows;
	
	private int num_cols;
	
	/**
	 * Creates a new instance of the Sequential Re-embedding algorithm with the
	 * default threshold and without the "reset first" feature. 
	 * 
	 * <P>Mode of operation can be border length
	 * (@link OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN) or conflict index
	 * (@link OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN) minimization.</P>
	 * 
	 * @param mode conflict minimization mode
	 */
	public SequentialReembedding (int mode)
	{
		this(mode, false);
	}

	/**
	 * Creates a new instance of the Sequential Re-embedding algorithm with the
	 * default threshold.
	 * 
	 * <P>Mode of operation can be border length
	 * (@link OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN) or conflict index
	 * (@link OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN) minimization. The
	 * "reset first" feature allows the first pass to re-embed the probes
	 * considering only those that have already been re-embedded.</P>
	 * 
	 * @param mode conflict minimization mode
	 * @param reset_first whether to consider only re-embedded probes on the
	 * first pass over the chip
	 */
	public SequentialReembedding (int mode, boolean reset_first)
	{
		this (mode, reset_first, DEFAULT_THRESHOLD);
	}

	/**
	 * Creates a new instance of the Sequential Re-embedding algorithm with the
	 * specified configuration.
	 * 
	 * <P>Mode of operation can be border length
	 * (@link OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN) or conflict index
	 * (@link OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN) minimization. The
	 * "reset first" feature allows the first pass to re-embed the probes
	 * considering only those that have already been re-embedded. The threshold
	 * sets a limit on the minimum reduction of conflicts that must be achieved
	 * in one pass in order to continue. In other words, when the improvement
	 * drops below the threshold, the algorithm stops.</P>
	 * 
	 * @param mode conflict minimization mode
	 * @param reset_first whether to consider only re-embedded probes on the
	 * first pass over the chip
	 * @param threshold the minimum improvement that must be gained in on pass
	 * in order to continue the algorithm
	 */
	public SequentialReembedding (int mode, boolean reset_first,
			double threshold)
	{
		if (mode != OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN &&
			mode != OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN)
				throw new IllegalArgumentException
					("Unknown distance mode: " + mode);
			
		this.mode = mode;
		this.reset_first = reset_first;
		this.threshold = threshold;
	}

	/**
	 * Optimizes the layout of the chip by running the Sequential Re-embedding
	 * algorithm several times until no further improvement is possible (or
	 * until the improvement drops below the threshold.
	 * 
	 * @param chip chip instance to be optimized
	 */
	public void changeLayout (Chip chip)
	{
		this.embedder = OptimumSingleProbeEmbedding.createEmbedder(chip, mode);
		this.num_rows = chip.getNumberOfRows();
		this.num_cols = chip.getNumberOfColumns();
		
		if (chip instanceof SimpleChip)
		{
			optimize ((SimpleChip) chip);
		}
		else if (chip instanceof AffymetrixChip)
		{
			optimize ((AffymetrixChip) chip);
		}
		else
		{
			throw new IllegalArgumentException ("Unsupported chip type.");
		}
	}
	
	private void optimize (SimpleChip chip)
	{
		double	last_conf, curr_conf;
		boolean	reset;
		
		if (mode == OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN)
			curr_conf = LayoutEvaluation.borderLength(chip);
		else // CONFLICT_INDEX_MIN
			curr_conf = LayoutEvaluation.averageConflictIndex(chip);
		
		reset = this.reset_first;
		
		do
		{
			if (reset)
			{
				// reset only in the first execution
				copyAndResetSpots(chip);
				incrementalOptimization (chip);
				reset = false;
			}
			else
				totalOptimization (chip);
			
			last_conf = curr_conf;
			
			if (mode == OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN)
				curr_conf = LayoutEvaluation.borderLength(chip);
			else // CONFLICT_INDEX_MIN
				curr_conf = LayoutEvaluation.averageConflictIndex(chip);
			
		} while ((last_conf - curr_conf) / last_conf > threshold);
	}

	private void optimize (AffymetrixChip chip)
	{
		double	last_conf, curr_conf;
		boolean	reset;
		
		if (mode == OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN)
			curr_conf = LayoutEvaluation.borderLength(chip);
		else // CONFLICT_INDEX_MIN
			curr_conf = LayoutEvaluation.averageConflictIndex(chip);
		
		reset = this.reset_first;
		
		do
		{
			if (reset)
			{
				// reset only in the first execution
				copyAndResetSpots(chip);
				incrementalOptimization (chip);
				reset = false;
			}
			else
				totalOptimization (chip);
			
			last_conf = curr_conf;
			
			if (mode == OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN)
				curr_conf = LayoutEvaluation.borderLength(chip);
			else // CONFLICT_INDEX_MIN
				curr_conf = LayoutEvaluation.averageConflictIndex(chip);
			
		} while ((last_conf - curr_conf) / last_conf > threshold);
	}

	private void copyAndResetSpots (Chip chip)
	{
		this.spot_copy = new int [num_rows][num_cols];
		
		for (int r = 0; r < num_rows; r++)
			for (int c = 0; c < num_cols; c++)
			{
				this.spot_copy[r][c] = chip.spot[r][c];
				chip.spot[r][c] = Chip.EMPTY_SPOT;
			}
	}
	
	private void incrementalOptimization (SimpleChip chip)
	{
		int id;
		
		for (int r = 0; r < num_rows; r++)
			for (int c = 0; c < num_cols; c++)
			{
				if ((id = spot_copy[r][c]) == Chip.EMPTY_SPOT)
					continue;
				
				chip.spot[r][c] = id;
				
				embedder.reembedSpot(r, c, id);
			}
	}

	private void incrementalOptimization (AffymetrixChip chip)
	{
		int id;
		
		for (int r = 0; r < num_rows; r++)
			for (int c = 0; c < num_cols; c++)
			{
				if ((id = spot_copy[r][c]) == Chip.EMPTY_SPOT)
					continue;
				
				if (!chip.isPMProbe(id))
					continue;
				
				chip.spot[r][c] = id;
				chip.spot[r + 1][c] = spot_copy[r + 1][c];
				
				embedder.reembedSpot(r, c, id);
			}
	}
	
	private void totalOptimization (SimpleChip chip)
	{
		int id;
		
		for (int r = 0; r < num_rows; r++)
			for (int c = 0; c < num_cols; c++)
			{
				if ((id = chip.spot[r][c]) != Chip.EMPTY_SPOT)
					embedder.reembedSpot(r, c, id);
			}
	}

	private void totalOptimization (AffymetrixChip chip)
	{
		int id;
		
		for (int r = 0; r < num_rows; r++)
			for (int c = 0; c < num_cols; c++)
			{
				if ((id = chip.spot[r][c]) == Chip.EMPTY_SPOT)
					continue;
				
				if (chip.isPMProbe(id))
					embedder.reembedSpot(r, c, id);
			}
	}
	
	/**
	 * Returns the algorithm's name together with current options.
	 */
	@Override
	public String toString ()
	{
		String m, r;
		
		switch (this.mode)
		{
			case OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN:
				m = "-BL";
				break;
				
			case OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN:
				m = "-CI";
				break;
			
			default:
				m = "-?";
		}
		
		if (reset_first)
			r = "-Reset";
		else
			r = "-NoReset";
		
		return this.getClass().getSimpleName() + m + r;
	}
}
