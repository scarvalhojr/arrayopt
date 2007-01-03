/*
 * PriorityReembedding.java
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

import java.util.*;

/**
 * This class implements the new Priority Re-embedding post-placement
 * optimization.
 * 
 * <P>The algorithm is in many ways similiar to the
 * {@linkplain SequentialReembedding}: every probe of the chip is examined and
 * optimally re-embedded in regards to its neighbors using the Optimum
 * Single-Probe Embedding (OSPE) algorithm (implemented by
 * {@linkplain OptimumSingleProbeEmbedding}). The desired conflict minimization
 * function (whether border length or conflict index minimization) must also be
 * specified at instantiation time (to the constructor). See
 * {@link ConflictIndex} for the exact definition of conflict index).</P>  
 * 
 * <P>However, instead of scanning the chip top to bottom, left to right, this
 * algorithm can examine the spots according to one of three pre-defined orders
 * (the priorities). The idea of priorities derives from the following
 * observation: while some probe can be embedded into the deposition sequence in
 * several ways, others may have only a few valid embeddings. Those with a
 * lesser number of embeddings are less flexible and should be re-embedded
 * first. Another similar observation is that the number of re-embedded
 * neighbors also limit the re-embedding of a probe. Thus, probes with a greater
 * number of re-embedded probes should be processed first.</P>
 * 
 * <P>Three pre-defined priorities are available:
 *
 * <UL>
 * <LI>{@link #PRIORITY_NUM_OF_EMBEDDINGS}: probes with a lesser number of
 * embeddings are re-embedded first; in case of ties, probes with a greater
 * number of re-embedded neighbors are re-embedded first 
 * <LI>{@link #PRIORITY_NUM_OF_EMBEDDINGS}: probes with a greater number of
 * re-embedded neighbors are re-embedded first; in case of ties, probes with a
 * lesser number of embeddings are re-embedded first;
 * <LI>{@link #PRIORITY_BALANCED}: a single condition that is a combination of
 * the number of embeddings and number of re-embedded neighbors (the numbers are
 * brought to a similar scale with the log function)
 * </UL>
 * 
 * <P>Like the {@link SequentialReembedding}, this algorithm actually scans
 * the chip several times and, since the OSPE never increases the amount of
 * conflicts (sum of border lengths or conflict indices), it is guaranteed to
 * never produce a worse embedding. However, this also means that the
 * improvements decrease with each iteration until no more improvements are
 * possible. In other words, the algorithm finds a local optimum solution.</P>
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
public class PriorityReembedding implements LayoutAlgorithm
{
	public static final int PRIORITY_NUM_OF_EMBEDDINGS = 0;
	
	public static final int PRIORITY_NUM_OF_NEIGHBORS = 1;
	
	public static final int PRIORITY_BALANCED = 2;
	
	private PriorityQueue<PendingSpot> queue;
	
	private Comparator<PendingSpot> comparator;
	
	private OptimumSingleProbeEmbedding embedder;
	
	private int mode;
	
	private double threshold ;
	
	public static final double DEFAULT_THRESHOLD = 0.001d;
	
	private boolean reset_first;
	
	private int spot_copy[][];
	
	private int num_rows;
	
	private int num_cols;
	
	private BitSet spot_ready;
	
	private BitSet spot_added;
	
	private int num_probes;
	
	private long num_embed[];
	
	private final int ADD_REGION_DIM = 1;

	/**
	 * Creates a new instance of the Priority Re-embedding algorithm with the
	 * default threshold and without the "reset first" feature.
	 * 
	 * <P>Mode of operation can be border length
	 * (@link OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN) or conflict index
	 * (@link OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN) minimization. The
	 * priority can be either {@link #PRIORITY_NUM_OF_EMBEDDINGS},
	 * {@link #PRIORITY_NUM_OF_NEIGHBORS} or {@link #PRIORITY_BALANCED}.</P>
	 * 
	 * @param mode conflict minimization mode
	 * @param priority pre-defined spot priority
	 */
	public PriorityReembedding (int mode, int priority)
	{
		this (mode, priority, false, DEFAULT_THRESHOLD);
	}

	/**
	 * Creates a new instance of the Priority Re-embedding algorithm with the
	 * default threshold.
	 * 
	 * <P>Mode of operation can be border length
	 * (@link OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN) or conflict index
	 * (@link OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN) minimization. The
	 * priority can be either {@link #PRIORITY_NUM_OF_EMBEDDINGS},
	 * {@link #PRIORITY_NUM_OF_NEIGHBORS} or {@link #PRIORITY_BALANCED}. The
	 * "reset first" feature allows the first pass to re-embed the probes
	 * considering only those that have already been re-embedded.</P>
	 * 
	 * @param mode conflict minimization mode
	 * @param priority pre-defined spot priority
	 * @param reset_first whether to consider only re-embedded probes on the
	 * first pass over the chip
	 */
	public PriorityReembedding (int mode, int priority,
			boolean reset_first)
	{
		this (mode, priority, reset_first, DEFAULT_THRESHOLD);
	}

	/**
	 * Creates a new instance of the Priority Re-embedding algorithm with the
	 * specified configuration.
	 * 
	 * <P>Mode of operation can be border length
	 * (@link OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN) or conflict index
	 * (@link OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN) minimization. The
	 * priority can be either {@link #PRIORITY_NUM_OF_EMBEDDINGS},
	 * {@link #PRIORITY_NUM_OF_NEIGHBORS} or {@link #PRIORITY_BALANCED}. The
	 * "reset first" feature allows the first pass to re-embed the probes
	 * considering only those that have already been re-embedded. The threshold
	 * sets a limit on the minimum reduction of conflicts that must be achieved
	 * in one pass in order to continue. In other words, when the improvement
	 * drops below the threshold, the algorithm stops.</P>
	 * 
	 * @param mode conflict minimization mode
	 * @param priority pre-defined spot priority
	 * @param reset_first whether to consider only re-embedded probes on the
	 * first pass over the chip
	 * @param threshold the minimum improvement that must be gained in on pass
	 * in order to continue the algorithm
	 */
	public PriorityReembedding (int mode, int priority,
			boolean reset_first, double threshold)
	{
		if (mode != OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN &&
			mode != OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN)
				throw new IllegalArgumentException
					("Unknown distance mode: " + mode);

		if (priority == PRIORITY_NUM_OF_EMBEDDINGS)
			comparator = new EmbeddingsPriority();
		else if (priority == PRIORITY_NUM_OF_NEIGHBORS)
			comparator = new NeighborsPriority();
		else if (priority == PRIORITY_BALANCED)
			comparator = new BalancedPriority();
		else
			throw new IllegalArgumentException
				("Unknown priority mode: " + priority);

		this.mode = mode;
		this.reset_first = reset_first;
		this.threshold = threshold;
	}

	/**
	 * Optimizes the layout of the chip by running the Priority Re-embedding
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
		this.num_probes = chip.getNumberOfProbes();
		
		this.spot_ready = new BitSet (num_rows * num_cols);
		this.spot_added = new BitSet (num_rows * num_cols);
		
		if (reset_first)
			this.spot_copy = new int [num_rows][num_cols];
		else
			this.spot_copy = null;
		
		this.num_embed = new long [chip.getNumberOfProbes()];
		
		// create the priority queue with an initial
		// capacity of 20% of the number of probes
		int capacity = (int) (.2d * this.num_probes);
		this.queue = new PriorityQueue<PendingSpot> (capacity, comparator);
		
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
		long pivot_threshold;
		double last_conf, curr_conf, impr = 1;
		boolean reset = this.reset_first;
		
		if (mode == OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN)
			last_conf = LayoutEvaluation.borderLength(chip);
		else // CONFLICT_INDEX_MIN
			last_conf = LayoutEvaluation.averageConflictIndex(chip);
		
		pivot_threshold = analyzeProbes();
		
		while (impr > threshold)
		{
			spot_added.clear();
			spot_ready.clear();
			
			fixPivots (chip, pivot_threshold, reset);
			processQueue (chip, reset);
			scanRemainingSpots (chip);
			processQueue (chip, reset);
			
			// reset the first iteration only
			reset = false;

			if (mode == OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN)
				curr_conf = LayoutEvaluation.borderLength(chip);
			else // CONFLICT_INDEX_MIN
				curr_conf = LayoutEvaluation.averageConflictIndex(chip);
			
			impr = (last_conf - curr_conf) / last_conf;
			
			last_conf = curr_conf;
		}
	}

	private void optimize (AffymetrixChip chip)
	{
		long pivot_threshold;
		double last_conf, curr_conf, impr = 1;
		boolean reset = this.reset_first;
		
		if (mode == OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN)
			last_conf = LayoutEvaluation.borderLength(chip);
		else // CONFLICT_INDEX_MIN
			last_conf = LayoutEvaluation.averageConflictIndex(chip);
		
		pivot_threshold = analyzeProbes(chip);
				
		while (impr > threshold)
		{
			spot_added.clear();
			spot_ready.clear();
			
			fixPivots (chip, pivot_threshold, reset);
			processQueue (chip, reset);
			scanRemainingSpots (chip);
			processQueue (chip, reset);
			
			// reset the first iteration only
			reset = false;

			if (mode == OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN)
				curr_conf = LayoutEvaluation.borderLength(chip);
			else // CONFLICT_INDEX_MIN
				curr_conf = LayoutEvaluation.averageConflictIndex(chip);
			
			impr = (last_conf - curr_conf) / last_conf;
			
			last_conf = curr_conf;
		}
	}

	private void fixPivots (SimpleChip chip, long limit, boolean reset)
	{
		int id, r, c;
		
		// first find pivots and mark their spots as ready
		for (r = 0; r < num_rows; r++)
			for (c = 0; c < num_cols; c++)
			{
				if ((id = chip.spot[r][c]) == Chip.EMPTY_SPOT)
				{
					setReadySpot(r, c);
					continue;
				}
				
				if (num_embed[id] <= limit)
				{
					if (!reset && num_embed[id] > 1)
						embedder.reembedSpot(r, c, id);
					
					setReadySpot(r, c);
				}
				else if (reset)
				{
					spot_copy[r][c] = id;
					chip.spot[r][c] = Chip.EMPTY_SPOT;
				}
			}
		
		// now add immediate neighbors of pivots to the queue
		for (r = 0; r < num_rows; r++)
			for (c = 0; c < num_cols; c++)
				if (isReadySpot(r, c) && chip.spot[r][c] != Chip.EMPTY_SPOT)
					addNeighbors (chip, r, c);
	}

	private void fixPivots (AffymetrixChip chip, long limit, boolean reset)
	{
		int id, r, c;
		
		// first find pivots and mark their spots as ready
		for (r = 0; r < num_rows; r++)
			for (c = 0; c < num_cols; c++)
			{
				if ((id = chip.spot[r][c]) == Chip.EMPTY_SPOT)
				{
					setReadySpot(r, c);
					continue;
				}
				
				if (!chip.isPMProbe(id)) continue;
				
				if (num_embed[id] <= limit)
				{
					if (!reset && num_embed[id] > 1)
						embedder.reembedSpot(r, c, id);
					
					setReadySpot(r, c);
					setReadySpot(r + 1, c);
				}
				else if (reset)
				{
					spot_copy[r][c] = id;
					spot_copy[r + 1][c] = chip.spot[r + 1][c];
					chip.spot[r][c] = Chip.EMPTY_SPOT;
					chip.spot[r + 1][c] = Chip.EMPTY_SPOT;
				}
			}
		
		// now add immediate neighbors of pivots to the queue
		for (r = 0; r < num_rows; r++)
			for (c = 0; c < num_cols; c++)
				if ((id = chip.spot[r][c]) != Chip.EMPTY_SPOT)
					if (isReadySpot(r,c) && chip.isPMProbe(id))
						addNeighbors (chip, r, c);
	}

	private long analyzeProbes ()
	{
		long min = Long.MAX_VALUE;
		
		for (int p = 0; p < num_probes; p++)
			if ((num_embed[p] = embedder.numberOfEmbeddings(p)) < min)
				min = num_embed[p];
		
		return min;
	}

	private long analyzeProbes (AffymetrixChip chip)
	{
		long min = Long.MAX_VALUE;
		
		for (int p = 0; p < num_probes; p++)
			if (chip.isPMProbe(p))
				if ((num_embed[p] = embedder.numberOfEmbeddings(p)) < min)
					min = num_embed[p];

		return min;
	}
	
	private void processQueue (Chip chip, boolean reset)
	{
		ArrayList<PendingSpot> updated;
		Iterator<PendingSpot> iter;
		PendingSpot s;
		int row, col;

		// list of spots that were updated and need
		// to be re-inserted into the queue
		updated = new ArrayList<PendingSpot> ();
		
		while((s = queue.poll()) != null)
		{
			row = s.row;
			col = s.col;
			
			if (chip instanceof SimpleChip)
				restoreSpot ((SimpleChip) chip, row, col, reset);
			else
				restoreSpot ((AffymetrixChip) chip, row, col, reset);
						
			// check spots that need to be updated due to
			// the restoration of the last spot
			for (iter = queue.iterator(); iter.hasNext(); )
			{
				s = iter.next();
				
				if (s.addNeighbor(row, col))
		        {
					// remove spots that were updated (they need to be
					// re-inserted in order to restore the queue's ordering)
		        	iter.remove();
		        	updated.add(s);
		        }
			}
			
			// re-insert updated spots
			queue.addAll(updated);
			updated.clear();

			if (chip instanceof SimpleChip)
				addNeighbors ((SimpleChip) chip, row, col);
			else
				addNeighbors ((AffymetrixChip) chip, row, col);
		}
	}

	private void addNeighbors (SimpleChip chip, int row, int col)
	{
		if (mode == OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN)
			addImmediateNeighbors(chip, row, col);
		else	// CONFLICT_INDEX_MIN
			addRegionNeighbors(chip, row, col);
	}

	private void addNeighbors (AffymetrixChip chip, int row, int col)
	{
		if (mode == OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN)
			addImmediateNeighbors(chip, row, col);
		else	// CONFLICT_INDEX_MIN
			addRegionNeighbors(chip, row, col);
	}

	private void addImmediateNeighbors (SimpleChip chip, int row, int col)
	{
		if (row > 0)
			addSpot (chip, row - 1, col);	// top
		if (row < chip.getNumberOfRows() - 1)
			addSpot (chip, row + 1, col);	// bottom
		if (col > 0)
			addSpot (chip, row, col - 1);	// left
		if (col < chip.getNumberOfColumns() - 1)
			addSpot (chip, row, col + 1);	// right
	}

	private void addImmediateNeighbors (AffymetrixChip chip, int row, int col)
	{
		// coordinates (row, col) are guaranteed
		// to be of a spot containing a PM probe
		
		// note that only spots with PM probes are added to the queue;
		// usually, a row contains exclusively PM or MM probes but sometimes
		// some regions might be out of sync and thus we cannot assume where
		// the PM neighboring probes are; if the spot contains a MM probe,
		// the addSpot method will skip it 
		
		if (row > +1)
			addSpot (chip, row - 2, col);		// top
		if (row > 0)
			addSpot (chip, row - 1, col);		// top (probably a MM probe)

		if (row < chip.getNumberOfRows() - 2)
			addSpot (chip, row + 2, col);		// bottom 
		if (row < chip.getNumberOfRows() - 1)
			addSpot (chip, row + 1, col);		// bottom (probably a MM probe)
		
		if (col > 0)
			addSpot (chip, row, col - 1);		// left
		if (col > 0)
			addSpot (chip, row + 1, col - 1);	// left (probably a MM probe)
		
		if (col < chip.getNumberOfColumns() - 1)
			addSpot (chip, row, col + 1);		// right
		if (col < chip.getNumberOfColumns() - 1)
			addSpot (chip, row + 1, col + 1);	// right (probably a MM probe)
}

	private void addRegionNeighbors (SimpleChip chip, int row, int col)
	{
		int min_row = Math.max(row - ADD_REGION_DIM, 0);
		int min_col = Math.max(col - ADD_REGION_DIM, 0);
		int max_row = Math.min(row + ADD_REGION_DIM, num_rows - 1);
		int max_col = Math.min(col + ADD_REGION_DIM, num_cols - 1);
		
		for (int r = min_row; r <= max_row; r++)
			for (int c = min_col; c <= max_col; c++)
				addSpot (chip, r, c);
	}

	private void addRegionNeighbors (AffymetrixChip chip, int row, int col)
	{
		// coordinates (row, col) are guaranteed
		// to be of a spot containing a PM probe
		
		// note that only spots with PM probes are added to the queue;
		// usually, a row contains exclusively PM or MM probes but sometimes
		// some regions might be out of sync and thus we cannot assume where
		// the PM neighboring probes are; if the spot contains a MM probe,
		// the addSpot method will skip it
		
		// TODO define exactly what the ADD_REGION_DIM really means
		// the current use only adds and subtracts from the row number
		// in order to try and reach a spot with a PM probe 

		int min_row = Math.max(row - 2 * ADD_REGION_DIM, 0);
		int min_col = Math.max(col - ADD_REGION_DIM, 0);
		int max_row = Math.min(row + 2 * ADD_REGION_DIM, num_rows - 1);
		int max_col = Math.min(col + ADD_REGION_DIM, num_cols - 1);
		
		for (int r = min_row; r <= max_row; r++)
			for (int c = min_col; c <= max_col; c++)
				addSpot (chip, r, c);
	}
	
	private void scanRemainingSpots (SimpleChip chip)
	{
		for (int r = 0; r < num_rows; r++)
			for (int c = 0; c < num_cols; c++)
				// only non-empty non-ready spots are added
				addSpot(chip, r, c);
	}

	private void scanRemainingSpots (AffymetrixChip chip)
	{
		for (int r = 0; r < num_rows; r++)
			for (int c = 0; c < num_cols; c++)
				// only non-empty non-ready spots are added
				addSpot(chip, r, c);
	}
	
	private void restoreSpot (SimpleChip chip, int row, int col, boolean reset)
	{
		int id;

		if (reset)
			chip.spot[row][col] = spot_copy[row][col];
		
		id = chip.spot[row][col];
		
		setReadySpot (row, col);
		
		embedder.reembedSpot(row, col, id);
	}

	private void restoreSpot (AffymetrixChip chip, int row, int col,
			boolean reset)
	{
		int id;
		
		if (reset)
		{
			chip.spot[row][col] = spot_copy[row][col];
			chip.spot[row + 1][col] = spot_copy[row + 1][col];
		}
		
		id = chip.spot[row][col];
		
		setReadySpot (row, col);
		setReadySpot (row + 1, col);
		
		embedder.reembedSpot(row, col, id);
	}
	
	private void addSpot (SimpleChip chip, int row, int col)
	{
		PendingSpot s;
		int id;
		
		if (isReadySpot(row, col) || isAddedSpot(row, col))
			return;
		
		if (spot_copy != null)
			id = spot_copy[row][col];
		else
			id = chip.spot[row][col];
		
		if (mode == OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN)
		{
			int neighbors = countImmediateNeighbors (chip, row, col); 
			s = new SimpleBL (row, col, num_embed[id], neighbors);
		}
		else // CONFLICT_INDEX_MIN
		{
			double neighbors = countRegionNeighbors (chip, row, col);
			s = new SimpleCI (row, col, num_embed[id], neighbors);
		}
		
		queue.offer(s);
		
		setAddedSpot(row, col);
	}

	private void addSpot (AffymetrixChip chip, int row, int col)
	{
		PendingSpot s;
		int id;
		
		if (isReadySpot(row, col) || isAddedSpot(row, col))
			return;
		
		if (spot_copy != null)
			id = spot_copy[row][col];
		else
			id = chip.spot[row][col];

		// only spots with PM probes are inserted into the queue
		if (!chip.isPMProbe(id))
			return;

		if (mode == OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN)
		{
			int neighbors = countImmediateNeighbors (chip, row, col);
			s = new AffyBL (row, col, num_embed[id], neighbors);
		}
		else // CONFLICT_INDEX_MIN
		{
			double neighbors = countRegionNeighbors (chip, row, col);
			s = new AffyCI (row, col, num_embed[id], neighbors);
		}
		
		queue.offer(s);
		
		setAddedSpot(row, col);
	}
	
	private int countImmediateNeighbors (SimpleChip chip, int row, int col)
	{
		int neighbors = 0;
		
		if (row > 0)
			if (chip.spot[row - 1][col] != Chip.EMPTY_SPOT &&
				isReadySpot(row - 1, col))
					neighbors++;
		
		if (row < chip.getNumberOfRows() - 1)
			if (chip.spot[row + 1][col] != Chip.EMPTY_SPOT &&
				isReadySpot(row + 1, col))
					neighbors++;
		
		if (col > 0)
			if (chip.spot[row][col - 1] != Chip.EMPTY_SPOT &&
				isReadySpot(row, col - 1))
					neighbors++;
		
		if (col < chip.getNumberOfColumns() - 1)
			if (chip.spot[row][col + 1] != Chip.EMPTY_SPOT &&
				isReadySpot(row, col + 1))
					neighbors++;
		
		return neighbors;
	}
	
	private int countImmediateNeighbors (AffymetrixChip chip, int row, int col)
	{
		// we assume that the spot has a PM probe
		int pm_row = row, mm_row = row + 1;
		int neighbors = 0;
		
		if (pm_row > 0)
			if (chip.spot[pm_row - 1][col] != Chip.EMPTY_SPOT &&
				isReadySpot(pm_row - 1, col))
					neighbors++;
		
		if (mm_row < chip.getNumberOfRows() - 1)
			if (chip.spot[mm_row + 1][col] != Chip.EMPTY_SPOT &&
				isReadySpot(mm_row + 1, col))
					neighbors++;
		
		if (col > 0)
			if (chip.spot[pm_row][col - 1] != Chip.EMPTY_SPOT &&
				isReadySpot(pm_row, col - 1))
					neighbors++;
			else if (chip.spot[mm_row][col - 1] != Chip.EMPTY_SPOT &&
					 isReadySpot(mm_row, col - 1))
					neighbors++;
		
		if (col < chip.getNumberOfColumns() - 1)
			if (chip.spot[pm_row][col + 1] != Chip.EMPTY_SPOT &&
				isReadySpot(pm_row, col + 1))
					neighbors++;
			else if (chip.spot[mm_row][col + 1] != Chip.EMPTY_SPOT &&
					 isReadySpot(mm_row, col + 1))
					neighbors++;
		
		return neighbors;
	}
	
	private double countRegionNeighbors (SimpleChip chip, int row, int col)
	{
		int dim;
		double n = 0;
		
		dim = ConflictIndex.dimConflictRegion();
		
		int min_row = Math.max (row - dim, 0);
		int max_row = Math.min (row + dim, chip.getNumberOfRows() - 1);
		int min_col = Math.max (col - dim, 0);
		int max_col = Math.min (col + dim, chip.getNumberOfColumns() - 1);
		
		for (int r = min_row; r <= max_row; r++)
		for (int c = min_col; c <= max_col; c++)
		{
			if (chip.spot[r][c] == Chip.EMPTY_SPOT)
				continue;
			
			if (!isReadySpot(r, c))
				continue;
			
			if (r == row && c == col)
				continue;

			n += ConflictIndex.distanceWeight(row, col, r, c);
			n += ConflictIndex.distanceWeight(r, c, row, col);
		}
		
		return n;
	}

	private double countRegionNeighbors (AffymetrixChip chip, int row, int col)
	{
		int dim, r_pm, r_mm;
		double n = 0;
		
		r_pm = row;
		r_mm = row + 1;
		dim = ConflictIndex.dimConflictRegion();
		
		int min_row = Math.max (r_pm - dim, 0);
		int max_row = Math.min (r_mm + dim, chip.getNumberOfRows() - 1);
		int min_col = Math.max (col - dim, 0);
		int max_col = Math.min (col + dim, chip.getNumberOfColumns() - 1);
		
		for (int r = min_row; r <= max_row; r++)
		for (int c = min_col; c <= max_col; c++)
		{
			if (chip.spot[r][c] == Chip.EMPTY_SPOT)
				continue;
			
			if (!isReadySpot(r, c))
				continue;
			
			if ((r == r_pm && c == col) || (r == r_mm && c == col))
				continue;
			
			if (r_pm >= r - dim && r_pm <= r + dim)
			{
				n += ConflictIndex.distanceWeight(r_pm, col, r, c);
				n += ConflictIndex.distanceWeight(r, c, r_pm, col);
			}
			
			if (r_mm >= r - dim && r_mm <= r + dim)
			{
				n += ConflictIndex.distanceWeight(r_mm, col, r, c);
				n += ConflictIndex.distanceWeight(r, c, r_mm, col);
			}
		}
		
		return n;
	}

	private void setReadySpot (int row, int col)
	{
		spot_ready.set(row * num_cols + col);
	}

	private boolean isReadySpot (int row, int col)
	{
		return spot_ready.get(row * num_cols + col);
	}

	private void setAddedSpot (int row, int col)
	{
		spot_added.set(row * num_cols + col);
	}

	private boolean isAddedSpot (int row, int col)
	{
		return spot_added.get(row * num_cols + col);
	}
	
	/**
	 * Returns the algorithm's name together with current options.
	 */
	@Override
	public String toString ()
	{
		String m, p, r;
		
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
		
		if (comparator instanceof EmbeddingsPriority)
			p = "-NumOfEmbed";
		else if (comparator instanceof NeighborsPriority)
			p = "-NumOfNeighbors";
		else if (comparator instanceof BalancedPriority)
			p = "-Balanced";
		else
			p = "-?";

		if (reset_first)
			r = "-Reset-";
		else
			r = "-NoReset-";
		
		return this.getClass().getSimpleName() + m + p + r + threshold;
	}

	/**
	 * Abstract class representing a pending spot, i.e. a spot whose probe still
	 * needs to be re-embedded. This class implements the Comparable interface
	 * so that its instances can be added to a priority queue. The order imposed
	 * allows the priority queue to return spots with probes that have the least
	 * degree of freedom in their embeddings. In order to achieve that, a spot
	 * <CODE>a</CODE> is considered less than spot <CODE>b</CODE> if and only
	 * if the probe of <CODE>a</CODE> has less valid of embeddings than the
	 * probe of <CODE>b</CODE>. In case of ties, <CODE>a</CODE> is considered
	 * less than <CODE>b</CODE> if <CODE>a</CODE> has less ready immediate
	 * neighbors (called borders). In case of ties...  
	 * 
	 * <P/>Note: this class has a natural ordering that is inconsistent with
	 * equals.</P> 
	 */
	private static abstract class PendingSpot
	{
		protected int row;
		protected int col;
		protected double neighbors;
		double embeds;
		
		private PendingSpot (int row, int col, double embeds, double neighbors)
		{
			this.row = row;
			this.col = col;
			this.embeds = embeds;
			this.neighbors = neighbors;
		}
		
		@Override
		public String toString ()
		{
			return "Spot " + row + "," + col + " (" + embeds + " embeddings, " +
						neighbors +	" finished neighbors)";
		}
		
		abstract boolean addNeighbor (int r, int c);
	}
	
	private static abstract class BorderLengthSpot extends PendingSpot
	{
		private BorderLengthSpot (int row, int col, long embeds,
				double neighbors)
		{
			super (row, col, Math.log10(embeds), neighbors);
			
			// the number of embeddings is scaled down to a number closer
			// to the number of embeddings so that they can be better
			// combined in case of PRIORITY_BALANCED
		}
	}

	private static abstract class ConflictIndexSpot extends PendingSpot
	{
		protected int ci_dim;
		
		private ConflictIndexSpot (int row, int col, long embeds,
				double neighbors)
		{
			// the number of embeddings is scaled down to a number closer
			// to the number of embeddings so that they can be better
			// combined in case of PRIORITY_BALANCED; log (base e) is used
			// since the number of neighbors here is greater than in the
			// case of BorderLengthSpot
			super (row, col, Math.log(embeds), neighbors);			
			
			ci_dim = ConflictIndex.dimConflictRegion();
		}		
	}
	
	private static class SimpleBL extends BorderLengthSpot
	{
		SimpleBL (int row, int col, long embeds, int neighbors)
		{
			super (row, col, embeds, neighbors);
		}
		
		@Override
		boolean addNeighbor (int r, int c)
		{
			int	add = 0;
			
			if (r == this.row)
			{
				if (c == this.col - 1)
					add = 1;
				else if (c == this.col + 1)
					add = 1;
			}
			else if (c == this.col)
			{
				if (r == this.row - 1)
					add = 1;
				else if (r == this.row + 1)
					add = 1;
			}
			
			if (add > 0)
			{
				this.neighbors += add;
				return true;
			}
			// else
				return false;
		}
	}

	private static class SimpleCI extends ConflictIndexSpot
	{
		SimpleCI (int row, int col, long embeds, double neighbors)
		{
			super (row, col, embeds, neighbors);
		}
		
		@Override
		boolean addNeighbor (int r, int c)
		{
			double	add = 0;
			
			if (this.row  < r - ci_dim || this.row > r + ci_dim)
				return false;
			
			if (this.col < c - ci_dim || this.col > c + ci_dim)
				return false;
			
			add += ConflictIndex.distanceWeight(this.row, this.col, r, c);
			add += ConflictIndex.distanceWeight(r, c, this.row, this.col);
			
			if (add > 0)
			{
				this.neighbors += add;
				return true;
			}
			// else
				return false;
		}
	}

	private static class AffyBL extends BorderLengthSpot
	{
		AffyBL (int row, int col, long embeds, int neighbors)
		{
			super (row, col, embeds, neighbors);
		}
		
		@Override
		boolean addNeighbor (int r, int c)
		{
			int	pm_row, mm_row, add = 0;
			
			// this class is intended to hold spots with PM probes only
			pm_row = this.row;
			mm_row = this.row + 1;
						
			if (r == pm_row || r == mm_row)
			{
				if (c == this.col - 1)
					add = 1;
				else if (c == this.col + 1)
					add = 1;
			}
			else if (c == this.col)
			{
				if (r == pm_row - 2 || r == pm_row + 2)
					add = 1;
			}
			
			if (add > 0)
			{
				this.neighbors += add;
				return true;
			}
			// else
				return false;
		}		
	}

	private static class AffyCI extends ConflictIndexSpot
	{
		AffyCI (int row, int col, long embeds, double neighbors)
		{
			super (row, col, embeds, neighbors);						
		}
		
		@Override
		boolean addNeighbor (int r, int c)
		{
			return addNeighborSingleProbe (r, c) ||
					addNeighborSingleProbe (r + 1, c);
		}
		
		private boolean addNeighborSingleProbe (int r, int c)
		{

			int		r_pm, r_mm;
			double	add = 0;
			
			r_pm = this.row;
			r_mm = r_pm + 1;
			
			if (this.col < c - ci_dim || this.col > c + ci_dim)
				return false;
			
			if (r_pm >= r - ci_dim && r_pm <= r + ci_dim)
			{
				add += ConflictIndex.distanceWeight(r_pm, this.col, r, c);
				add += ConflictIndex.distanceWeight(r, c, r_pm, this.col);
			}
			
			if (r_mm >= r - ci_dim && r_mm <= r + ci_dim)
			{
				add += ConflictIndex.distanceWeight(r_mm, this.col, r, c);
				add += ConflictIndex.distanceWeight(r, c, r_mm, this.col);
			}
			
			if (add > 0)
			{
				this.neighbors += add;
				return true;
			}
			// else
				return false;			
		}		
	}
	
	private static class EmbeddingsPriority implements Comparator<PendingSpot>
	{
		public int compare (PendingSpot s1, PendingSpot s2)
		{
			if (s1.equals(s2)) return 0;
			
			// a spot whose probe has less number of embeddings
			// has less flexibility and should come first
			if (s1.embeds < s2.embeds)
				return -1;
			
			if (s1.embeds > s2.embeds)
				return +1;
			
			// in case of ties, we look at the number of neighbors:
			// a spot with greater number of ready neighbors
			// has less flexibility and should come first
			if (s1.neighbors > s2.neighbors)
				return -1;
			
			if (s1.neighbors < s2.neighbors)
				return +1;

			// ideally, ties would not need to be broken since the priority
			// queue should be able to handle objects of equal order; however,
			// the PriorityQueue.remove(Object) is currently removing an object
			// based on its ordering and not based on Object.equals(Object),
			// thus another object may be deleted just because it has the same
			// ordering (compareTo returns zero) - this problem has been
			// reported as a Java bug (bug ID 6268068); moreover, by returning
			// a non-zero value, the comparator becomes consistent with equals
			return -1;
		}
	}

	private static class NeighborsPriority implements Comparator<PendingSpot>
	{
		public int compare (PendingSpot s1, PendingSpot s2)
		{
			if (s1.equals(s2)) return 0;

			// a spot with greater number of ready neighbors
			// has less flexibility and should come first
			if (s1.neighbors > s2.neighbors)
				return -1;
			
			if (s1.neighbors < s2.neighbors)
				return +1;
			
			// in case of ties, we look at the number of embeddings:
			// a spot whose probe has less number of embeddings
			// has less flexibility and should come first
			if (s1.embeds < s2.embeds)
				return -1;
			
			if (s1.embeds > s2.embeds)
				return +1;

			// CAUTION: see note on the return statement of the
			// EmbeddingsPriority.compare method
			return -1;
		}
	}

	private static class BalancedPriority implements Comparator<PendingSpot>
	{
		public int compare (PendingSpot s1, PendingSpot s2)
		{
			double cmp;
			
			if (s1.equals(s2)) return 0;
			
			cmp = (s1.embeds - s2.embeds);
			cmp += (s2.neighbors - s1.neighbors);
			
			if (cmp < 0)
				return -1;
			
			if (cmp > 0)
				return +1;
			
			// CAUTION: see note on the return statement of the
			// EmbeddingsPriority.compare method
			return -1;
		}
	}
}
