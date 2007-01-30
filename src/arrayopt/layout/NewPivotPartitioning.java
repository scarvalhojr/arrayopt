/*
 * NewPivotPartitioning.java
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
public class NewPivotPartitioning implements LayoutAlgorithm
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
	
	private static final int NUM_SEEDS = 100;

	private OptimumSingleProbeEmbedding ospe;

	private Chip chip;

	private int pid[];

	/**
	 * TODO use this to speed up the selection of pivots by computing the
	 * number of embedding of every probe at once, storing the results on this
	 * array sorting the probes by the number of embeddings
	 */
	// private long num_embed[];

	private double dist[];

	private int offset;

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
	public NewPivotPartitioning (FillingAlgorithm filler, int mode, int max_depth)
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
		
		if (max_depth < 0)
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
		
		// create probe distance array
		dist = new double[nonpivots];
		dist_sort = new DistanceSorting (pid, dist, offset);
		
		probe_len = chip.getProbeLength();
		pivot1 = new float[probe_len];
		pivot2 = new float[probe_len];
		
		horizontalDivide (1, chip.getChipRegion(), 0, pivots - 1, pivots,
			pid.length - 1);
	}
	
	private int selectPivots ()
	{
		int min_pivots, num_pivots, end, start_fake;
		
		// TODO remove this
		// System.err.println("Selecting pivots...");

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
			long noe = ospe.numberOfEmbeddings(pid[0]);
			System.err.println(num_pivots + " fake pivots selected => noe: " + noe);
		}
		else
		{
			// no: the selected pivots are 'real' pivots, i.e.
			// they have only 1 embedding (or 2 in case of Affy chips)
			start_fake = num_pivots;
			
			// TODO remove this
			long noe = ospe.numberOfEmbeddings(pid[0]);
			System.err.println(num_pivots + " true pivots selected => noe: " + noe);
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
			long noe = ospe.numberOfEmbeddings(pid[num_pivots - 1]);
			System.err.println("Pivot list extended to " + num_pivots + " => noe: " + noe);
		}
		
		if (start_fake < num_pivots)
		{
			// cut the excess of fake pivots
			if (num_pivots > min_pivots)
				num_pivots = min_pivots;
			
			// TODO remove this
			System.err.println("Pivot list trimmed to " + num_pivots);
			
			// TODO re-embed fake probes with the CenteredEmbedding algorithm
			// or implement a different one that move unmasked steps to the
			// extremities as far as possible
		}
		
		// TODO remove this
		/*
		System.err.println("List of pivots:");
		for (int p = 0; p < num_pivots; p++)
		{
			chip.printEmbedding(pid[p]);
			System.err.println("\tID " + pid[p]);
		}
		//*/
		
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
		
		// TODO remove this
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
			part_pivot = choosePivotPair_bl_new (f_pivot, l_pivot);
		else // if (mode == OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN)
			part_pivot = choosePivotPair_ci (f_pivot, l_pivot);
		
		// TODO restore this
		// part_probe = divideProbes (pid[f_pivot], pid[l_pivot], f_probe,l_probe);
		// makePseudoPivots(f_pivot, part_pivot, l_pivot);
		part_probe = divideProbes (f_probe, l_probe);
		
		// check how many (non-pivot) probes were assigned to each region
		probes_1 = part_probe - f_probe;
		probes_2 = l_probe - part_probe + 1;
		
		// check the total number of probes of each region (including pivots)
		total_1 = (part_pivot - f_pivot) + probes_1;
		total_2 = (l_pivot - part_pivot + 1) + probes_2; 
		
		// divide the region proportionally
		rows_1 = rows_per_probe * (int) Math.round((double)total_1 / num_cols);
		
		// make sure each region gets at least one row
		if (rows_1 == 0)
			rows_1 = rows_per_probe;
		else if (rows_1 == num_rows)
			rows_1 = num_rows - rows_per_probe;
		
		rows_2 = num_rows - rows_1;

		// check how many spots each region has
		spots_1 = num_cols * rows_1 / rows_per_probe;
		spots_2 = num_cols * rows_2 / rows_per_probe;

		// TODO remove this
		System.err.println("Partitioning (%): " + (100 * rows_1 / num_rows) + " x " + (100 * rows_2 / num_rows));
		/*
		System.err.println("Total1: " + total_1);
		System.err.println("Total2: " + total_2);
		System.err.println("Rows1: " + rows_1);
		System.err.println("Rows2: " + rows_2);
		System.err.println("Spots1: " + spots_1);
		System.err.println("Spots2: " + spots_2);
		System.err.println("Exceed: " + (total_1 - spots_1));
		//*/

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

		// TODO remove this
		/*
		int p, p1, p2, dd1, dd2, count1 = 0, count2 = 0;
		p1 = pid[f_pivot];
		p2 = pid[l_pivot];
		System.err.println("********************************************");
		System.err.println("Probes assigned to P1 (" + p1 + ")");
		chip.printEmbedding(p1);
		System.err.println("\n--------------------------------------------");
		for (p = f_probe; p < part_probe; p++)
		{
			dd1 = LayoutEvaluation.hammingDistance(chip, p1, pid[p]);
			dd2 = LayoutEvaluation.hammingDistance(chip, p2, pid[p]);
			chip.printEmbedding(pid[p]);
			if (dd1 < dd2)
				System.err.println("\t" + dd1 + " < " + dd2 + " => " + (dd1 - dd2));
			else if (dd1 == dd2)
				System.err.println("\t" + dd1 + " = " + dd2 + " => " + (dd1 - dd2));
			else
				System.err.println("\t" + dd1 + " > " + dd2 + " => " + (dd1 - dd2));
			count1++;
		}
		System.err.println(count1 + " probes");
		System.err.println("********************************************");
		System.err.println("Probes assigned to P2 (" + p2 + ")");
		chip.printEmbedding(p2);
		System.err.println("\n--------------------------------------------");
		for (; p <= l_probe; p++)
		{
			dd1 = LayoutEvaluation.hammingDistance(chip, p1, pid[p]);
			dd2 = LayoutEvaluation.hammingDistance(chip, p2, pid[p]);
			chip.printEmbedding(pid[p]);
			if (dd1 < dd2)
				System.err.println("\t" + dd1 + " < " + dd2 + " => " + (dd1 - dd2));
			else if (dd1 == dd2)
				System.err.println("\t" + dd1 + " = " + dd2 + " => " + (dd1 - dd2));
			else
				System.err.println("\t" + dd1 + " > " + dd2 + " => " + (dd1 - dd2));
			count2++;
		}
		System.err.println(count2 + " probes");
		System.err.println("********************************************");
		//*/
		
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
		
		// TODO remove this
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
			part_pivot = choosePivotPair_bl_new (f_pivot, l_pivot);
		else // if (mode == OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN)
			part_pivot = choosePivotPair_ci (f_pivot, l_pivot);
		
		// TODO remove this
		// System.err.println("Pivots chosen: " + pid[f_pivot] + " and " + pid[l_pivot]);
		
		// TODO restore this
		// part_probe = divideProbes (pid[f_pivot], pid[l_pivot], f_probe,l_probe);
		// makePseudoPivots(f_pivot, part_pivot, l_pivot);
		part_probe = divideProbes (f_probe, l_probe);

		// check how many (non-pivot) probes were assigned to each region
		probes_1 = part_probe - f_probe;
		probes_2 = l_probe - part_probe + 1;
		
		// check the total number of probes of each region (including pivots)
		total_1 = (part_pivot - f_pivot) + probes_1;
		total_2 = (l_pivot - part_pivot + 1) + probes_2; 
		
		// divide the region proportionally
		cols_1 = (int) Math.round((double)total_1 / (num_rows/rows_per_probe));
		
		// make sure each region gets at least one column
		if (cols_1 == 0)
			cols_1 = 1;
		else if (cols_1 == num_cols)
			cols_1 = num_cols - 1;
		
		cols_2 = num_cols - cols_1;
		
		// TODO remove this
		System.err.println("Partitioning (%): " + (100 * cols_1 / num_cols) + " x " + (100 * cols_2 / num_cols));
		
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
		int i, j, last_seed, p1, p2, d, maxdist, tmp;
		
		// TODO remove this
		// System.err.println("Choosing pivot pair");
		// long time = System.nanoTime();
		
		p1 = first;
		p2 = last;
		maxdist = -1;
		
		// find pair of pivots p1 and p2 with maximum Hamming distance
		last_seed = first + NUM_SEEDS - 1;
		if (last_seed > last) last_seed = last;
		for (i = first; i < last_seed; i++)
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
		
		// TODO remove this
		// System.err.println("P1 and P2 are pid[" + p1 + "]=" + pid[p1] + " and pid[" + p2 + "]=" + pid[p2] + " with distance " + maxdist);
		
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

		// TODO remove this
		/*
		// System.err.println("Chose pivot pair in " + ((System.nanoTime() - time)/Math.pow(10,9)) + " seconds with distance " + maxdist);
		System.err.println("P1: " + p1);
		chip.printEmbedding(p1);
		System.err.println("\nP2: " + p2);
		chip.printEmbedding(p2);
		System.err.println("\nDistance: " + LayoutEvaluation.hammingDistance(chip, p1, p2));
		//*/

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
		
		// TODO is this THE bug?
		// move p2 to the beginning of its own list
		// pid[last] = pid[i];
		// pid[i] = p2;
		
		// TODO remove this
		/*
		int p, dd1, dd2;
		System.err.println("********************************************");
		System.err.println("Pivots assigned to P1 (" + p1 + ")");
		chip.printEmbedding(p1);
		System.err.println("\n--------------------------------------------");
		for (p = first + 1; p < i; p++)
		{
			dd1 = LayoutEvaluation.hammingDistance(chip, p1, pid[p]);
			dd2 = LayoutEvaluation.hammingDistance(chip, p2, pid[p]);
			chip.printEmbedding(pid[p]);
			if (dd1 < dd2)
				System.err.println("\t" + dd1 + " < " + dd2);
			else if (dd1 == dd2)
				System.err.println("\t" + dd1 + " = " + dd2);
			else
				System.err.println("\t" + dd1 + " > " + dd2);
		}
		System.err.println(count1 + " pivots");
		System.err.println("********************************************");
		System.err.println("Pivots assigned to P2 (" + p2 + ")");
		chip.printEmbedding(p2);
		System.err.println("\n--------------------------------------------");
		for (; p < last; p++)
		{
			dd1 = LayoutEvaluation.hammingDistance(chip, p1, pid[p]);
			dd2 = LayoutEvaluation.hammingDistance(chip, p2, pid[p]);
			chip.printEmbedding(pid[p]);
			if (dd1 < dd2)
				System.err.println("\t" + dd1 + " < " + dd2);
			else if (dd1 == dd2)
				System.err.println("\t" + dd1 + " = " + dd2);
			else
				System.err.println("\t" + dd1 + " > " + dd2);
		}
		System.err.println(count2 + " pivots");
		System.err.println("********************************************");
		//*/

		return i;
	}

	private int choosePivotPair_bl_new (int first, int last)
	{
		int i, j, last_seed, p1, p2, d, maxdist, tmp;
		
		// TODO remove this
		// System.err.println("Choosing pivot pair");
		// long time = System.nanoTime();
		
		p1 = first;
		p2 = last;
		maxdist = -1;
		
		// find pair of pivots p1 and p2 with maximum Hamming distance
		last_seed = first + NUM_SEEDS - 1;
		if (last_seed > last) last_seed = last;
		for (i = first; i < last_seed; i++)
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
		
		// TODO remove this
		// System.err.println("P1 and P2 are pid[" + p1 + "]=" + pid[p1] + " and pid[" + p2 + "]=" + pid[p2] + " with distance " + maxdist);
		
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

		// TODO remove this
		//*
		// System.err.println("Chose pivot pair in " + ((System.nanoTime() - time)/Math.pow(10,9)) + " seconds with distance " + maxdist);
		System.err.println("P1: " + p1);
		chip.printEmbedding(p1);
		System.err.println("\nP2: " + p2);
		chip.printEmbedding(p2);
		System.err.println("\nDistance: " + LayoutEvaluation.hammingDistance(chip, p1, p2));
		//*/
		
		startPseudoPivots(p1, p2);

		// partition the remaining pivots into two groups
		// according to whether they are closer to p1 or p2
		int id;
		
		for (i = first + 1, j = last - 1; i <= j;)
		{
			id = pid[i];
			if (distanceToPivots (id) <= 0)
			{
				// assign pivot to p1
				i++;
				p1_count++;
				updatePseudoPivot(id, pivot1);
			}
			else
			{
				// assign pivot to p2
				tmp = pid[i];
				pid[i] = pid[j];
				pid[j] = tmp;
				j--;
				p2_count++;
				updatePseudoPivot(id, pivot2);
			}
		}
		
		// TODO remove this
		System.err.println((i - first) + " x " + (last - i + 1));
		
		// TODO is this THE bug?
		// move p2 to the beginning of its own list
		// pid[last] = pid[i];
		// pid[i] = p2;
		
		// TODO remove this
		/*
		int p;
		float dd;
		System.err.println("********************************************");
		System.err.println("Pivots assigned to P1 (" + p1 + ")");
		chip.printEmbedding(p1);
		System.err.println("\n--------------------------------------------");
		for (p = first + 1; p < i; p++)
		{
			dd = distanceToPivots(pid[p]);
			chip.printEmbedding(pid[p]);
			System.err.println("\t" + dd);
		}
		// System.err.println(count1 + " pivots");
		System.err.println("********************************************");
		System.err.println("Pivots assigned to P2 (" + p2 + ")");
		chip.printEmbedding(p2);
		System.err.println("\n--------------------------------------------");
		for (; p < last; p++)
		{
			dd = distanceToPivots(pid[p]);
			chip.printEmbedding(pid[p]);
			System.err.println("\t" + dd);
		}
		// System.err.println(count2 + " pivots");
		System.err.println("********************************************");
		//*/

		return i;
	}

	private int choosePivotPair_ci (int first, int last)
	{
		int i, j, last_seed, p1, p2, tmp;
		double d, maxdist;
		
		// TODO remove this
		// System.err.println("Choosing pivot pair");
		// long time = System.nanoTime();

		p1 = first;
		p2 = last;
		maxdist = -1;
		
		// find pair of pivots p1 and p2 with maximum Hamming distance
		last_seed = first + NUM_SEEDS - 1;
		if (last_seed > last) last_seed = last;
		for (i = first; i < last_seed; i++)
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
		
		// TODO is this THE bug?
		// move p2 to the beginning of its own list
		// pid[last] = pid[i];
		// pid[i] = p2;

		// TODO remove this
		// System.err.println("Chose pivot pair in " + ((System.nanoTime() - time)/Math.pow(10,9)) + " seconds with distance " + maxdist);

		return i;
	}

	private int divideProbes (int p1, int p2, int first, int last)
	{
		int i, total, count1 = 0, count2 = 0, count_any, delta;
		double d;
		
		total = last - first + 1;
		
		// TODO remove this
		// System.err.println("Dividing probes based on pivots " + p1 + " and " + p2);
		
		// compute and save the minimum distance of every probe to pivot p1
		for (i = first; i <= last; i++)
			dist[i - offset] = LayoutEvaluation.hammingDistance (chip, pid[i], p1);
		
		// subtract the the min distance to p1 by the min distance to p2
		for (i = first; i <= last; i++)
		{
			d = dist[i - offset] -= LayoutEvaluation.hammingDistance (chip, pid[i], p2);
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
		else // if (delta == 0)
		{
			count1 += count_any / 2;
		}
		
		// TODO remove this
		/*
		int p, dd1, dd2, mycount;
		System.err.println("********************************************");
		System.err.println("Probes assigned to pivot p1 (" + p1 + ")");
		chip.printEmbedding(p1);
		System.err.println("\n--------------------------------------------");
		for (mycount = 0, p = first; p < (first+count1); p++)
		{
			dd1 = LayoutEvaluation.hammingDistance(chip, p1, pid[p]);
			dd2 = LayoutEvaluation.hammingDistance(chip, p2, pid[p]);
			chip.printEmbedding(pid[p]);
			if (dd1 < dd2)
				System.err.println("\t" + dd1 + " < " + dd2 + " => " + (dd1 - dd2));
			else if (dd1 == dd2)
				System.err.println("\t" + dd1 + " = " + dd2 + " => " + (dd1 - dd2));
			else
				System.err.println("\t" + dd1 + " > " + dd2 + " => " + (dd1 - dd2));
			mycount++;
		}
		System.err.println(mycount + " probes");
		System.err.println("********************************************");
		System.err.println("Probes assigned to pivot p2 (" + p2 + ")");
		chip.printEmbedding(p2);
		System.err.println("\n--------------------------------------------");
		for (mycount = 0;p <= last; p++)
		{
			dd1 = LayoutEvaluation.hammingDistance(chip, p1, pid[p]);
			dd2 = LayoutEvaluation.hammingDistance(chip, p2, pid[p]);
			chip.printEmbedding(pid[p]);
			if (dd1 < dd2)
				System.err.println("\t" + dd1 + " < " + dd2 + " => " + (dd1 - dd2));
			else if (dd1 == dd2)
				System.err.println("\t" + dd1 + " = " + dd2 + " => " + (dd1 - dd2));
			else
				System.err.println("\t" + dd1 + " > " + dd2 + " => " + (dd1 - dd2));
			mycount++;
		}
		System.err.println(mycount + " probes");
		System.err.println("********************************************");
		//System.err.println("i = " + i + "; first=" + first + " + count1=" + count1 + " = " + (first + count1));
		//*/
		
		// return the index of the first probe assigned to p2
		return first + count1;
	}
	
	private float pivot1[];
	private int p1_count;
	private float pivot2[];
	private int p2_count;
	private int probe_len;
	
	private void makePseudoPivots (int f_pivot, int p_pivot, int l_pivot)
	{
		int i, pos, w, bitmask = 0;
		
		p1_count = p_pivot - f_pivot;
		p2_count = l_pivot - p_pivot + 1;
		
		for (w = -1, pos = 0; pos < probe_len; pos++)
		{
			if ((pos % Integer.SIZE) == 0)
			{
				w++;
				bitmask = 0x01 << (Integer.SIZE - 1);
			}
			else
				bitmask >>>= 1;
			
			pivot1[pos] = pivot1[pos] = 0;
			
			for (i = f_pivot; i < p_pivot; i++)
				if ((chip.embed[pid[i]][w] & bitmask) != 0)
					pivot1[pos]++;
			
			for (; i <= l_pivot; i++)
				if ((chip.embed[pid[i]][w] & bitmask) != 0)
					pivot2[pos]++;
		}
		
		// TODO remove this
		//*
		System.err.println("Pseudo pivot 1:");
		for (pos = 0; pos < probe_len; pos++)
			System.err.printf("%3.2f|", pivot1[pos]/p1_count);
		System.err.println("\nPseudo pivot 2:");
		for (pos = 0; pos < probe_len; pos++)
			System.err.printf("%3.2f|", pivot2[pos]/p2_count);
		System.err.println();
		//*/
	}
	
	private void startPseudoPivots (int p1, int p2)
	{
		int pos, w, bitmask = 0;
		
		p1_count = 1;
		p2_count = 1;
		
		for (w = -1, pos = 0; pos < probe_len; pos++)
		{
			if ((pos % Integer.SIZE) == 0)
			{
				w++;
				bitmask = 0x01 << (Integer.SIZE - 1);
			}
			else
				bitmask >>>= 1;
			
			pivot1[pos] = pivot2[pos] = 0;
			
			if ((chip.embed[p1][w] & bitmask) != 0)
				pivot1[pos]++;
			
			if ((chip.embed[p2][w] & bitmask) != 0)
				pivot2[pos]++;
		}
		
		// TODO remove this
		/*
		System.err.println("Pseudo pivot 1:");
		for (pos = 0; pos < probe_len; pos++)
			System.err.printf("%3.2f|", pivot1[pos]/p1_count);
		System.err.println("\nPseudo pivot 2:");
		for (pos = 0; pos < probe_len; pos++)
			System.err.printf("%3.2f|", pivot2[pos]/p2_count);
		System.err.println();
		//*/
	}
	
	private void updatePseudoPivot (int id, float pseudo[])
	{
		int pos, w, bitmask = 0;
		
		for (w = -1, pos = 0; pos < probe_len; pos++)
		{
			if ((pos % Integer.SIZE) == 0)
			{
				w++;
				bitmask = 0x01 << (Integer.SIZE - 1);
			}
			else
				bitmask >>>= 1;
			
			if ((chip.embed[id][w] & bitmask) != 0)
				pseudo[pos]++;
		}
		
		// TODO remove this
		/*
		System.err.println("********************************************");
		if (pseudo == pivot1)
		{
			System.err.println("Adding probe ID " + id + " to pseudo pivot p1");
			chip.printEmbedding(id);
			System.err.println("\nUpdated pseudo pivot p1:");
			for (pos = 0; pos < probe_len; pos++)
				System.err.printf("%3.2f|", pivot1[pos]/p1_count);
		}
		else
		{
			System.err.println("Adding probe ID " + id + " to pseudo pivot p2");
			chip.printEmbedding(id);
			System.err.println("\nUpdated pseudo pivot p2:");
			for (pos = 0; pos < probe_len; pos++)
				System.err.printf("%3.2f|", pivot2[pos]/p2_count);
		}
		System.err.println();
		//*/
	}
	
	private float distanceToPivots (int id)
	{
		int pos, w, bitmask = 0;
		float d = 0;
		
		for (w = -1, pos = 0; pos < probe_len; pos++)
		{
			if ((pos % Integer.SIZE) == 0)
			{
				w++;
				bitmask = 0x01 << (Integer.SIZE - 1);
			}
			else
				bitmask >>>= 1;
			
			if ((chip.embed[id][w] & bitmask) == 0)
				d += (pivot1[pos] / p1_count - pivot2[pos] / p2_count);
			else
				// d += (1 - pivot1[pos]/p1_count) - (1 - pivot2[pos]/p2_count);
				d += pivot2[pos]/p2_count - pivot1[pos]/p1_count;
		}

		return d;
	}
	
	private int divideProbes (int first, int last)
	{
		int id, total, count1 = 0, count2 = 0;
		
		total = last - first + 1;
		
		for (int p = first; p <= last; p++)
		{
			id = pid[p];
			if ((dist[p - offset] = distanceToPivots (id)) <= 0)
			{
				count1++;
				p1_count++;
				updatePseudoPivot (id, pivot1);
			}
			else
			{
				count2++;
				p2_count++;
				updatePseudoPivot (id, pivot2);
			}
		}
		
		// sort probes by the difference of the distances
		QuickSort.sort(dist_sort, first, total);
		
		// TODO remove this
		/*
		int p, mycount;
		float dd;
		System.err.println("********************************************");
		System.err.println("Probes assigned to pseudo pivot p1");
		for (mycount = 0, p = first; p < (first+count1); p++)
		{
			dd = distanceToPivots(pid[p]);
			chip.printEmbedding(pid[p]);
			System.err.printf("%6.4f\t%6.4f\n", dd, dist[p - offset]);
			mycount++;
		}
		System.err.println(mycount + " probes");
		System.err.println("********************************************");
		System.err.println("Probes assigned to pseudo pivot p2");
		for (mycount = 0;p <= last; p++)
		{
			dd = distanceToPivots(pid[p]);
			chip.printEmbedding(pid[p]);
			System.err.printf("%6.4f\t%6.4f\n", dd, dist[p - offset]);
			mycount++;
		}
		System.err.println(mycount + " probes");
		System.err.println("********************************************");
		//*/
		
		// return the index of the first probe assigned to p2
		return first + count1;
	}
	
	private int dividePivots (int first, int last)
	{
		int id, total, count1 = 0, count2 = 0;
		
		total = last - first + 1;
		
		for (int p = first; p <= last; p++)
		{
			id = pid[p];
			if ((dist[p - offset] = distanceToPivots (id)) <= 0)
			{
				count1++;
				p1_count++;
				updatePseudoPivot (id, pivot1);
			}
			else
			{
				count2++;
				p2_count++;
				updatePseudoPivot (id, pivot2);
			}
		}
		
		// sort probes by the difference of the distances
		QuickSort.sort(dist_sort, first, total);
		
		// TODO remove this
		/*
		int p, mycount;
		float dd;
		System.err.println("********************************************");
		System.err.println("Probes assigned to pseudo pivot p1");
		for (mycount = 0, p = first; p < (first+count1); p++)
		{
			dd = distanceToPivots(pid[p]);
			chip.printEmbedding(pid[p]);
			System.err.printf("%6.4f\t%6.4f\n", dd, dist[p - offset]);
			mycount++;
		}
		System.err.println(mycount + " probes");
		System.err.println("********************************************");
		System.err.println("Probes assigned to pseudo pivot p2");
		for (mycount = 0;p <= last; p++)
		{
			dd = distanceToPivots(pid[p]);
			chip.printEmbedding(pid[p]);
			System.err.printf("%6.4f\t%6.4f\n", dd, dist[p - offset]);
			mycount++;
		}
		System.err.println(mycount + " probes");
		System.err.println("********************************************");
		//*/
		
		// return the index of the first probe assigned to p2
		return first + count1;
	}
	
	private int fillRegion (RectangularRegion region, int f_pivot, int l_pivot,
			int f_probe, int l_probe)
	{
		int num_pivots, num_probes, all[];
		
		num_pivots = l_pivot - f_pivot + 1;
		num_probes = l_probe - f_probe + 1;		
		
		// create an array with all probe IDs (pivots and non-pivots)
		all = new int [num_pivots + num_probes];
		System.arraycopy(pid, f_pivot, all, 0, num_pivots);
		System.arraycopy(pid, f_probe, all, num_pivots, num_probes);
		
		return filler.fillRegion(chip, region, all);
	}
	
	private class DistanceSorting implements ArrayIndexedCollection
	{
		private int probe_id[];
		
		private double probe_dist[];
		
		private int off;
		
		private double pivot;
		
		DistanceSorting (int probe_id[], double probe_dist[], int offset)
		{
			this.probe_id = probe_id;
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