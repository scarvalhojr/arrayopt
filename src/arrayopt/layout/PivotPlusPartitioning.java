/*
 * PivotPlusPartitioning.java
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

import arrayopt.util.ArrayIndexedCollection;
import arrayopt.util.QuickSort;

/**
 * TODO document this
 * 
 * <P><B>Selection of pivots</B>. Normally, all probes that have at most 1
 * possible embedding into the deposition sequence (or 2 in the case of an
 * {@link AffymetrixChip}) are selected as pivots. These are called 'true'
 * pivots. However, some chips have no such probes, i.e. all probes have more
 * than 1 (or 2) possible embeddings. In this case, the probes with the minimum
 * number of embeddings are promoted to pivots (they are then called 'fake'
 * pivots).</P>
 * 
 * <P>This algorithm defines what is the minimum number of pivots that is needed
 * in order guarantee a good partitioning of the chip. If it decides that there
 * are not enough true pivots for a good partitioning, it promotes the probes
 * with the next minimum number of embeddings to pivots (as fake pivots). This
 * is repeated until the minimum number of pivots is reached.</P>
 * 
 * <P>The minimum number of pivots is usually 2 times 2 to the power of
 * {@link #max_depth}. This is the number of pivots needed to partition the
 * chip all the way down to the maximum partitioning depth. However, if the
 * {@link #max_depth} is set too high, this might force too many probes to
 * be promoted as fake probes, specially if the chip is not large enough.
 * Having too many fake probes can degrade the quality of the layout since
 * pivots are considered as having a single fixed embedding. To prevent the
 * promotion of too many fake probes, the minimum number of pivots cannot
 * exceed a percentage of the total number of probes defined by the
 * {@link #MIN_PERCENTAGE_PIVOTS} constant.</P>
 * 
 * @author Anna Domanski
 * @author Ronny Gaertner
 * @author Sergio A. de Carvalho Jr.
 */
public class PivotPlusPartitioning implements LayoutAlgorithm
{
	/**
	 * This variable stores the current minimization mode used by the algorithm.
	 * Possible values are {@link OptimumSingleProbeEmbedding#BORDER_LENGTH_MIN}
	 * and {@link OptimumSingleProbeEmbedding#CONFLICT_INDEX_MIN}. 
	 */
	private int mode;
	
	/**
	 * Filling algorithm used to place the probes in each final sub-region.
	 */
	private FillingAlgorithm filler;

	/**
	 * TODO document this
	 */
	private int max_depth;

	/**
	 * TODO document this
	 */
	public static final double MIN_PERCENTAGE_PIVOTS = 0.01;

	private OptimumSingleProbeEmbedding ospe;

	private Chip chip;

	private int pid[];

	private boolean fake_pivots;

	private long rank[];

	/**
	 * TODO use this to speed up the selection of pivots by computing the
	 * number of embedding of every probe at once, storing the results on this
	 * array sorting the probes by the number of embeddings
	 */
	private long num_embed[];

	private double dist[];

	private int offset;

	private RankSorting rank_sort;

	private DistanceSorting dist_sort;

	private int rows_per_probe;
	
	/**
	 * Creates an instance of the Pivot Partitioning algorithm with the
	 * specified filling algorithm, minimization mode and maximum partitioning
	 * depth.
	 * 
	 * @param filler filling algorithm used in the final sub-regions
	 * @param mode minimization mode
	 * @param max_depth maximum partitioning depth
	 */
	public PivotPlusPartitioning (FillingAlgorithm filler, int mode, int max_depth)
	{
		switch (mode)
		{
			case OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN:
			case OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN:
				this.mode = mode;
				break;
				
			default:
				throw new IllegalArgumentException
					("Illegal minimization mode: " + mode);
		}
		
		if (max_depth < 1)
			throw new IllegalArgumentException
				("Illegal maximum depth: " + max_depth);
		
		this.filler = filler;
		this.max_depth = max_depth;
	}

	/**
	 * Creates a new layout of a microarray chip using the Pivot Partitioning
	 * algorithm.
	 * 
	 * @param chip instance of a microarray chip
	 */
	public void changeLayout (Chip c)
	{
		int pivots, nonpivots;
		
		if (c instanceof SimpleChip)
			rows_per_probe = 1;
		else if (c instanceof AffymetrixChip)
			rows_per_probe = 2;
		else
			throw new IllegalArgumentException ("Unsupported chip type.");
		
		this.chip = c;
		this.ospe = OptimumSingleProbeEmbedding.createEmbedder(c, mode);
				
		// reset current layout (if any)
		chip.resetLayout();
		
		this.pid = chip.getMovableProbes ();
		if (pid.length < 1)
			throw new IllegalStateException
				("No probes available for placing.");

		// select probes which will serve as pivots
		pivots = selectPivots ();
		nonpivots = pid.length - pivots;
		
		// offset marks the index of the first non-pivot probe
		this.offset = pivots;
		
		// create probe ranking array
		rank = new long[nonpivots];
		rank_sort = new RankSorting (pid, rank, pivots);
		computeProbeRanks (pivots, pid.length - 1);
		
		// create probe distance array
		dist = new double[nonpivots];
		dist_sort = new DistanceSorting (pid, rank, dist, offset);
		
		horizontalDivide (1, chip.getChipRegion(), 0, pivots - 1, pivots,
			pid.length - 1);
	}
	
	private int selectPivots ()
	{
		int min_pivots, num_pivots, end, start_fake;
		
		// first select probes with the minimum number of embeddings
		num_pivots = findPivots (0, end = pid.length - 1);

		// check whether selected pivots are 'real' or 'fake' pivots:
		// 'real' pivots are those with only 1 possible embedding
		// (or 2 in case of Affymetrix chips since the middle bases of PM/MM
		// pairs usually have some degree of freedom); pivots with more than 2
		// possible embeddings are considered 'fake' pivots
		if (ospe.numberOfEmbeddings(pid[0]) > 2)
		{
			// yes: we call these probes as 'fake' pivots
			start_fake = 0;
			
			// TODO remove this
			// long noe = ospe.numberOfEmbeddings(pid[0]);
			// System.err.println(num_pivots + " fake pivots selected => noe: " + noe);
		}
		else
		{
			// no: the selected pivots are 'real' pivots, i.e.
			// they have only 1 embedding (or 2 in case of Affy chips)
			start_fake = num_pivots;
			
			// TODO remove this
			// long noe = ospe.numberOfEmbeddings(pid[0]);
			// System.err.println(num_pivots + " true pivots selected => noe: " + noe);
		}
		
		// set the minimum number of pivots as the maximum between:
		// a) 2 * 2 ^ max_depth
		// b) MIN_PERCENTAGE_PIVOTS * the total number of probes
		min_pivots = Math.max(
					2 * (int) Math.pow(2, max_depth),
					(int) (MIN_PERCENTAGE_PIVOTS * chip.getNumberOfProbes()));
		
		while (num_pivots < min_pivots)
		{
			// find more probes with the next minimum number of embeddings
			num_pivots = findPivots (num_pivots, end);
			
			// TODO remove this
			// long noe = ospe.numberOfEmbeddings(pid[num_pivots - 1]);
			// System.err.println("Pivot list extended to " + num_pivots + " => noe: " + noe);
		}
		
		if (start_fake < num_pivots)
		{
			// turn on the fake pivots flag
			this.fake_pivots = true;
			
			// cut the excess of fake pivots
			if (num_pivots > min_pivots)
				num_pivots = min_pivots;
			
			// TODO remove this
			// System.err.println("Pivot list trimmed to " + num_pivots);
			
			// TODO re-embed fake probes with the CenteredEmbedding algorithm
			// or implement a different one that move unmasked steps to the
			// extremities as far as possible
		}
		else
		{
			// turn off the fake pivots flag
			this.fake_pivots = false;
		}
		
		return num_pivots;
	}
	
	private int findPivots (int start, int end)
	{
		int i, border, tmp;
		long noe, min;
		
		min = ospe.numberOfEmbeddings(pid[start]);
		
		// select probes with the lowest number of embeddings as pivots
		for (border = start + 1, i = border; i <= end; i++)
		{
			if ((noe = ospe.numberOfEmbeddings(pid[i])) == min)
			{
				tmp = pid[border];
				pid[border] = pid[i];
				pid[i] = tmp;
				border++;
			}
			else if (noe < min)
			{
				min = noe;
				tmp = pid[start];
				pid[start] = pid[i];
				pid[i] = tmp;
				border = start + 1;
			}
		}
		
		return border;
	}
	
	private int horizontalDivide (int depth, RectangularRegion r,
			int f_pivot, int l_pivot, int f_probe, int l_probe)
	{
		RectangularRegion top, bottom;
		int part_pivot, part_probe, num_rows, num_cols, exceed, unplaced;
		int probes_1, total_1, rows_1, spots_1;
		int probes_2, total_2, rows_2, spots_2;
		
		/*
		System.err.println("\nhorizontalDivide (depth " + depth + "; max " + max_depth + ")");
		System.err.println("Region: " + (r.last_row - r.first_row + 1) +
				" x " + (r.last_col - r.first_col + 1) + " (" + r + ")");
		System.err.println("Pivots: " + (l_pivot - f_pivot + 1) + " pivots (" +
				f_pivot + " - " + l_pivot + ")");
		System.err.println("Probes: " + (l_probe - f_probe + 1) + " probes (" +
				f_probe + " - " + l_probe + ")");
		//*/
		
		// stop partitioning if reached maximum depth or
		// if the number of probes or pivots is insufficient
		if ((depth > max_depth) || (l_pivot <= f_pivot) || (l_probe <= f_probe))
			return fillRegion (r, f_pivot, l_pivot, f_probe, l_probe);
		
		num_rows = r.last_row - r.first_row + 1;
		num_cols = r.last_col - r.first_col + 1;
		
		if (num_rows <= rows_per_probe)
		{
			if (num_cols <= 1)
			{
				// region cannot be partitioned anymore
				return fillRegion (r, f_pivot, l_pivot, f_probe, l_probe);
			}

			// region can still be vertically partitined
			return verticalDivide (depth, r, f_pivot,l_pivot,f_probe,l_probe);
		}
		
		// select a pair of pivots, p1 and p2, with max hamming distance
		// between them, and assign the remaining pivots to the chosen pivot 
		// with minimum hamming distance 
		if (mode == OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN)
			part_pivot = choosePivotPair_bl (f_pivot, l_pivot);
		else // if (mode == OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN)
			part_pivot = choosePivotPair_ci (f_pivot, l_pivot);
		 
		part_probe = divideProbes (pid[f_pivot], pid[l_pivot], f_probe,l_probe);
		
		// check how many (non-pivot) probes were assigned to each region
		probes_1 = part_probe - f_probe;
		probes_2 = l_probe - part_probe + 1;
		
		// check the total number of probes of each region (including pivots)
		total_1 = (part_pivot - f_pivot) + probes_1;
		total_2 = (l_pivot - part_pivot + 1) + probes_2; 
		
		// divide the region proportionally
		rows_1 = rows_per_probe * (int) Math.round((double)probes_1 / num_cols);
		
		// make sure each region gets at least one row
		if (rows_1 == 0)
			rows_1 = rows_per_probe;
		else if (rows_1 == num_rows)
			rows_1 = num_rows - rows_per_probe;
		
		rows_2 = num_rows - rows_1;
		
		// TODO remove this
		// System.err.println("Partitioning (%): " + (100 * rows_1 / num_rows) + " x " + (100 * rows_2 / num_rows));
		
		// check how many spots each region has
		spots_1 = num_cols * rows_1 / rows_per_probe;
		spots_2 = num_cols * rows_2 / rows_per_probe;
		
		// check if each region can handle the assigned probes
		if ((exceed = total_1 - spots_1) > 0)
		{
			if (exceed <= probes_1)
			{
				part_probe -= exceed;
			}
			else
			{
				part_probe -= probes_1;
				exceed -= probes_1;
				part_pivot -= exceed;
			}
		}
		else if ((exceed = total_2 - spots_2) > 0)
		{
			if (exceed <= probes_2)
			{
				part_probe += exceed;
			}
			else
			{
				part_probe += probes_2;
				exceed -= probes_2;
				part_pivot += exceed;
			}
		}
		
		top = new RectangularRegion(r.first_row, r.first_row + rows_1 - 1,
									r.first_col, r.last_col);
		
		unplaced = verticalDivide (depth + 1, top, f_pivot, part_pivot - 1,
						f_probe, part_probe - 1);

		bottom = new RectangularRegion(r.first_row + rows_1, r.last_row,
										r.first_col, r.last_col);

		unplaced += verticalDivide (depth + 1, bottom, part_pivot, l_pivot,
						part_probe, l_probe);

		return unplaced;
	}

	private int verticalDivide (int depth, RectangularRegion r,
			int f_pivot, int l_pivot, int f_probe, int l_probe)
	{
		RectangularRegion left, right;
		int part_pivot, part_probe, num_rows, num_cols, exceed, unplaced;
		int probes_1, total_1, cols_1, spots_1;
		int probes_2, total_2, cols_2, spots_2;
		
		/*
		System.err.println("\nverticalDivide (depth " + depth + "; max " + max_depth + ")");
		System.err.println("Region: " + (r.last_row - r.first_row + 1) +
				" x " + (r.last_col - r.first_col + 1) + " (" + r + ")");
		System.err.println("Pivots: " + (l_pivot - f_pivot + 1) + " pivots (" +
				f_pivot + " - " + l_pivot + ")");
		System.err.println("Probes: " + (l_probe - f_probe + 1) + " probes (" +
				f_probe + " - " + l_probe + ")");
		//*/
		
		// stop partitioning if reached maximum depth or
		// if the number of probes or pivots is insufficient
		if ((depth > max_depth) || (l_pivot <= f_pivot) || (l_probe <= f_probe))
			return fillRegion (r, f_pivot, l_pivot, f_probe, l_probe);

		num_rows = r.last_row - r.first_row + 1;
		num_cols = r.last_col - r.first_col + 1;
		
		if (num_cols <= 1)
		{
			if (num_rows <= rows_per_probe)
			{
				// region cannot be partitioned anymore
				return fillRegion (r, f_pivot, l_pivot, f_probe, l_probe);
			}

			// region can still be horizontally partitined
			return horizontalDivide (depth, r, f_pivot,l_pivot,f_probe,l_probe);
		}

		// select a pair of pivots, p1 and p2, with max hamming distance
		// between them, and assign the remaining pivots to the chosen pivot 
		// with minimum hamming distance
		if (mode == OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN)
			part_pivot = choosePivotPair_bl (f_pivot, l_pivot);
		else // if (mode == OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN)
			part_pivot = choosePivotPair_ci (f_pivot, l_pivot);
		 
		part_probe = divideProbes (pid[f_pivot], pid[l_pivot], f_probe,l_probe);
		
		// check how many (non-pivot) probes were assigned to each region
		probes_1 = part_probe - f_probe;
		probes_2 = l_probe - part_probe + 1;
		
		// check the total number of probes of each region (including pivots)
		total_1 = (part_pivot - f_pivot) + probes_1;
		total_2 = (l_pivot - part_pivot + 1) + probes_2; 
		
		// divide the region proportionally
		cols_1 = (int) Math.round((double)probes_1 / (num_rows/rows_per_probe));
		
		// make sure each region gets at least one column
		if (cols_1 == 0)
			cols_1 = 1;
		else if (cols_1 == num_cols)
			cols_1 = num_cols - 1;
		
		cols_2 = num_cols - cols_1;
		
		// TODO remove this
		// System.err.println("Partitioning (%): " + (100 * cols_1 / num_cols) + " x " + (100 * cols_2 / num_cols));
		
		// check how many spots each region has
		spots_1 = cols_1 * num_rows / rows_per_probe;
		spots_2 = cols_2 * num_rows / rows_per_probe;
		
		// check if each region can handle the assigned probes
		if ((exceed = total_1 - spots_1) > 0)
		{
			if (exceed <= probes_1)
			{
				part_probe -= exceed;
			}
			else
			{
				part_probe -= probes_1;
				exceed -= probes_1;
				part_pivot -= exceed;
			}
		}
		else if ((exceed = total_2 - spots_2) > 0)
		{
			if (exceed <= probes_2)
			{
				part_probe += exceed;
			}
			else
			{
				part_probe += probes_2;
				exceed -= probes_2;
				part_pivot += exceed;
			}
		}
		
		left = new RectangularRegion(r.first_row, r.last_row,
									r.first_col, r.first_col + cols_1 - 1);
		
		unplaced = horizontalDivide (depth + 1, left, f_pivot, part_pivot - 1,
						f_probe, part_probe - 1);

		right = new RectangularRegion(r.first_row, r.last_row,
									r.first_col + cols_1, r.last_col);

		unplaced += horizontalDivide (depth + 1, right, part_pivot, l_pivot,
						part_probe, l_probe);

		return unplaced;
	}
	
	private int choosePivotPair_bl (int first, int last)
	{
		int i, j, p1, p2, d, maxdist, tmp;
		
		p1 = first;
		p2 = last;
		maxdist = -1;
		
		// find pair of pivots p1 and p2 with maximum Hamming distance
		for (i = first; i < last; i++)
			for (j = i + 1; j <= last; j++)
			{
				d = LayoutEvaluation.hammingDistance(chip, pid[i], pid[j]);
				if (d > maxdist)
				{
					maxdist = d;
					p1 = i;
					p2 = j;
				}
			}
		
		// move pivots to extremities
		tmp = pid[first];
		pid[first] = pid[p1];
		pid[p1] = tmp;
		
		tmp = pid[last];
		pid[last] = pid[p2];
		pid[p2] = tmp;
		
		// from now on, p1 and p2 store the ID of the chosen pivots
		p1 = pid[first];
		p2 = pid[last];
		
		// partition the remaining pivots into two groups
		// according to whether they are closer to p1 or p2
		int dist1, dist2, count1 = 0, count2 = 0; 
		
		for (i = first + 1, j = last - 1; i <= j;)
		{
			dist1 = LayoutEvaluation.hammingDistance(chip, pid[i], p1);
			dist2 = LayoutEvaluation.hammingDistance(chip, pid[i], p2);
			
			if ((dist1 < dist2) || (dist1 == dist2 && count1 < count2))
			{
				// assign pivot to p1
				i++;
				count1++;
			}
			else
			{
				// assign pivot to p2
				tmp = pid[i];
				pid[i] = pid[j];
				pid[j] = tmp;
				j--;
				count2++;
			}
		}
		
		// move p2 to the beginning of its own list
		pid[last] = pid[i];
		pid[i] = p2;

		return i;
	}

	private int choosePivotPair_ci (int first, int last)
	{
		int i, j, p1, p2, tmp;
		double d, maxdist;
		
		p1 = first;
		p2 = last;
		maxdist = -1;
		
		// find pair of pivots p1 and p2 with maximum Hamming distance
		for (i = first; i < last; i++)
			for (j = i + 1; j <= last; j++)
			{
				d = LayoutEvaluation.conflictDistance(chip, pid[i], pid[j]);
				if (d > maxdist)
				{
					maxdist = d;
					p1 = i;
					p2 = j;
				}
			}
		
		// move pivots to extremities
		tmp = pid[first];
		pid[first] = pid[p1];
		pid[p1] = tmp;
		
		tmp = pid[last];
		pid[last] = pid[p2];
		pid[p2] = tmp;
		
		// from now on, p1 and p2 store the ID of the chosen pivots
		p1 = pid[first];
		p2 = pid[last];
		
		// partition the remaining pivots into two groups
		// according to whether they are closer to p1 or p2
		int count1 = 0, count2 = 0;
		double dist1, dist2;
		
		for (i = first + 1, j = last - 1; i <= j;)
		{
			dist1 = LayoutEvaluation.conflictDistance(chip, pid[i], p1);
			dist2 = LayoutEvaluation.conflictDistance(chip, pid[i], p2);
			
			if ((dist1 < dist2) || (dist1 == dist2 && count1 < count2))
			{
				// assign pivot to p1
				i++;
				count1++;
			}
			else
			{
				// assign pivot to p2
				tmp = pid[i];
				pid[i] = pid[j];
				pid[j] = tmp;
				j--;
				count2++;
			}
		}
		
		// move p2 to the beginning of its own list
		pid[last] = pid[i];
		pid[i] = p2;

		return i;
	}

	private int divideProbes (int p1, int p2, int first, int last)
	{
		int i, total, count1 = 0, count2 = 0, count_any, delta;
		double d;
		
		total = last - first + 1;
		
		// sort probes lexicographically
		QuickSort.sort(rank_sort, first, total);
		
		// compute and save the minimum distance of every probe to pivot p1
		dist[first - offset] = ospe.minDistanceProbe(pid[first], p1);
		for (i = first + 1; i <= last; i++)
			dist[i - offset] = ospe.minDistanceProbe(pid[i]);
		
		// subtract the the min distance to p1 by the min distance to p2
		d = dist[first - offset] -= ospe.minDistanceProbe(pid[first], p2);
		if (d < 0) count1++; else if (d > 0) count2++;
		for (i = first + 1; i <= last; i++)
		{
			d = dist[i - offset] -= ospe.minDistanceProbe(pid[i]);
			if (d < 0) count1++; else if (d > 0) count2++;
		}
		
		// sort probes by the difference of the distances
		QuickSort.sort(dist_sort, first, total);
		
		// count how many probes have the same minimum distance to p1 and p2
		count_any = total - count1 - count2;

		// divide the set of probes as evenly as possible
		if ((delta = count1 - count2) < 0)
		{
			delta = -delta;
			if (count_any >= delta)
			{
				count_any -= delta;
				count1 += delta + count_any / 2;
			}
			else
			{
				count1 += count_any;
			}
		}
		else if (delta > 0)
		{
			if (count_any >= delta)
			{
				count_any -= delta;
				count1 += count_any / 2;
			}
		}
		
		// return the index of the first probe assigned to p2
		return first + count1;
	}
	
	private void computeProbeRanks (int start, int end)
	{
		int  i, word, step, bitmask = 0;
		long base_mask;
		
		for (i = start; i <= end; i++)
			rank[i - offset] = 0;
		
		for (word = -1, step = 0; step < chip.embed_len; step++)
		{
			if (step % Integer.SIZE == 0)
			{
				bitmask = 0x01 << (Integer.SIZE - 1);
				word++;
			}
			else
				bitmask >>>= 1;

			switch (chip.dep_seq[step])
			{
				case 'A':
					base_mask = 0x00;
					break;
					
				case 'C':
					base_mask = 0x01;
					break;

				case 'G':
					base_mask = 0x02;
					break;
					
				case 'T':
					base_mask = 0x03;
					break;
				
				default:
					throw new IllegalArgumentException
						("Illegal deposition sequence.");
			}
			
			for (i = start; i <= end; i++)
				if ((bitmask & chip.embed[pid[i]][word]) != 0)
				{
					rank[i - offset] <<= 2;
					rank[i - offset] |= base_mask;
				}
		}
	}

	private int fillRegion (RectangularRegion region, int f_pivot, int l_pivot,
			int f_probe, int l_probe)
	{
		int i, num_pivots, num_probes, all[];
		
		num_pivots = l_pivot - f_pivot + 1;
		num_probes = l_probe - f_probe + 1;
		
		// TODO restore this
		/*
		if (fake_pivots)
		{
			// reembed fake pivots optimally in regards to the main pivot
			if (f_pivot + 1 <= l_pivot)
				ospe.reembedProbe(pid[f_pivot + 1], pid[f_pivot]);
			for (i = f_pivot + 2; i <= l_pivot; i++)
				ospe.reembedProbe(pid[i]);
		}
		//*/

		// TODO restore this
		/*
		// sort probes lexicographically to speed up re-embeddings
		QuickSort.sort(rank_sort, f_probe, l_probe - f_probe + 1);

		// reembed non-pivots optimally in regards to all pivots
		if (f_probe <= l_probe)
			ospe.reembedProbe(pid[f_probe], pid, f_pivot, l_pivot);
		for (i = f_probe + 1; i <= l_probe; i++)
			ospe.reembedProbe(pid[i]);
		//*/
		
		// create an array with all probe IDs (pivots and non-pivots)
		all = new int [num_pivots + num_probes];
		System.arraycopy(pid, f_pivot, all, 0, num_pivots);
		System.arraycopy(pid, f_probe, all, num_pivots, num_probes);
		
		return filler.fillRegion(chip, region, all);
	}

	private class RankSorting implements ArrayIndexedCollection
	{
		private int probe_id[];
		
		private long probe_rank[];
		
		private int off;
		
		private long pivot;
		
		RankSorting (int probe_id[], long probe_rank[], int offset)
		{
			this.probe_id = probe_id;
			this.probe_rank = probe_rank;
			this.off = offset;
		}
		
		public int compare (int i, int j)
		{
			i -= off;
			j -= off;
			
			return probe_rank[i] < probe_rank[j] ? -1 :
					probe_rank[i] == probe_rank[j] ? 0 : +1;
		}
		
		public void swap (int i, int j)
		{
			int tmp1;
			tmp1 = probe_id[i];
			probe_id[i] = probe_id[j];
			probe_id[j] = tmp1;
			
			i -= off;
			j -= off;
			
			long tmp2;
			tmp2 = probe_rank[i];
			probe_rank[i] = probe_rank[j];
			probe_rank[j] = tmp2;
		}
		
		public void setPivot (int i)
		{
			this.pivot = probe_rank[i - off];
		}
		
		public int compareToPivot (int i)
		{
			i -= off;
			
			return probe_rank[i] < pivot ? -1 :
					probe_rank[i] == pivot ? 0 : +1;
		}
		
		public int medianOfThree (int i, int j, int k)
		{
			long rank_i = probe_rank[i - off];
			long rank_j = probe_rank[j - off];
			long rank_k = probe_rank[k - off];

			return rank_i < rank_j ?
					(rank_j < rank_k ? j : (rank_i < rank_k ? k : i)) :
					(rank_j > rank_k ? j : (rank_i > rank_k ? k : i));
		}
	}
	
	private class DistanceSorting implements ArrayIndexedCollection
	{
		private int probe_id[];
		
		private long probe_rank[];
		
		private double probe_dist[];
		
		private int off;
		
		private double pivot;
		
		DistanceSorting (int probe_id[], long probe_rank[], double probe_dist[],
				int offset)
		{
			this.probe_id = probe_id;
			this.probe_rank = probe_rank;
			this.probe_dist = probe_dist;
			this.off = offset;
		}
		
		public int compare (int i, int j)
		{
			i -= off;
			j -= off;
			
			return probe_dist[i] < probe_dist[j] ? -1 :
					probe_dist[i] == probe_dist[j] ? 0 : +1;
		}
		
		public void swap (int i, int j)
		{
			int tmp1;
			tmp1 = probe_id[i];
			probe_id[i] = probe_id[j];
			probe_id[j] = tmp1;
			
			i -= off;
			j -= off;
			
			long tmp2;
			tmp2 = probe_rank[i];
			probe_rank[i] = probe_rank[j];
			probe_rank[j] = tmp2;
			
			double tmp3;
			tmp3 = probe_dist[i];
			probe_dist[i] = probe_dist[j];
			probe_dist[j] = tmp3;
		}
		
		public void setPivot (int i)
		{
			this.pivot = probe_dist[i - off];
		}
		
		public int compareToPivot (int i)
		{
			i -= off;
			
			return probe_dist[i] < pivot ? -1 :
					probe_dist[i] == pivot ? 0 : +1;
		}
		
		public int medianOfThree (int i, int j, int k)
		{
			double dist_i = probe_dist[i - off];
			double dist_j = probe_dist[j - off];
			double dist_k = probe_dist[k - off];

			return dist_i <= dist_j ?
					(dist_j <= dist_k ? j : (dist_i <= dist_k ? k : i)) :
					(dist_j >= dist_k ? j : (dist_i >= dist_k ? k : i));
		}
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
			case OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN:
				m = "-BL-";
				break;
				
			case OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN:
				m = "-CI-";
				break;
			
			default:
				m = "-?";
				break;
		}
		
		return this.getClass().getSimpleName() + m + max_depth + "-" + filler;
	}
}