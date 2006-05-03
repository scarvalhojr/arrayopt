/*
 * PivotPartitioning.java
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
 * <P>The minimum number of pivots is usually 2 to the power of
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
 */
public class PivotPartitioning implements PlacementAlgorithm
{
	/**
	 * TODO document this
	 */
	private FillingAlgorithm filler;

	/**
	 * TODO document this
	 */
	public static final int DEFAULT_MAX_DEPTH = 6;

	/**
	 * TODO document this
	 */
	private int max_depth;

	/**
	 * TODO document this
	 */
	public static final double MIN_PERCENTAGE_PIVOTS = 0.02;

	/**
	 * TODO document this
	 */
	private OptimumSingleProbeEmbedding ospe;

	/**
	 * TODO document this
	 */
	private Chip chip;

	/**
	 * TODO document this
	 */
	private int pid[];

	/**
	 * TODO document this
	 */
	private int rows_per_probe;

	/**
	 * TODO document this
	 */
	public PivotPartitioning (FillingAlgorithm filler)
	{
		this(filler, DEFAULT_MAX_DEPTH);
	}

	/**
	 * TODO document this
	 */
	public PivotPartitioning (FillingAlgorithm filler, int max_depth)
	{
		this.filler = filler;
		this.max_depth = (max_depth < 1) ? 1 : max_depth;
	}

	/**
	 * TODO document this
	 */
	public int makeLayout (Chip c)
	{
		int pivots;
		
		if (c instanceof SimpleChip)
			rows_per_probe = 1;
		else if (c instanceof AffymetrixChip)
			rows_per_probe = 2;
		else
			throw new IllegalArgumentException ("Unsupported chip type.");
		
		this.chip = c;
		this.ospe = OptimumSingleProbeEmbedding.createEmbedder(c);
		
		// reset current layout (if any)
		chip.resetLayout();
		
		this.pid = chip.getMovableProbes ();
		if (pid.length < 1)
			throw new IllegalStateException
				("No probes available for placing.");

		// select probes which will serve as pivots
		pivots = selectPivots ();
		
		return horizontalDivide (1, chip.getChipRegion(), 0, pivots - 1, pivots,
				pid.length - 1);
	}
	
	private int selectPivots ()
	{
		int min_pivots, num_pivots, end, start_fake;

		// set the minimum number of pivots to as the minimum between:
		// a) 2 ^ max_depth
		// b) MIN_PERCENTAGE_PIVOTS * the total number of probes
		min_pivots = Math.min(
					(int) Math.pow(2, max_depth),
					(int) (MIN_PERCENTAGE_PIVOTS * chip.getNumberOfProbes()));
		
		// first find probes with the minimum number of embeddings
		num_pivots = findPivots (0, end = pid.length - 1);
		
		// check if the minimum number of embeddins is greater than 2
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
		
		while (num_pivots < min_pivots)
		{
			// find more probes with the next minimum number of embeddings
			num_pivots = findPivots (num_pivots, end);
			
			// TODO remove this
			long noe = ospe.numberOfEmbeddings(pid[num_pivots - 1]);
			System.err.println("pivot list extended to: " + num_pivots + " => noe: " + noe);
		}
		
		if (start_fake < num_pivots)
		{
			// TODO re-embed fake probes with the CenteredEmbedding algorithm
			// or implement a different one that move unmasked steps to the
			// extremities as far as possible
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
		System.err.println("\nhorizontalDivide (depth " + depth + ")");
		System.err.println("Region: " + r);
		System.err.println("Pivots: " + f_pivot + " - " + l_pivot);
		System.err.println("Probes: " + f_probe + " - " + l_probe);
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
		part_pivot = choosePivotPair (f_pivot, l_pivot);
		 
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
		// System.err.println("Partitioning: " + (100 * rows_1 / num_rows) + " x " + (100 * rows_2 / num_rows));
		
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
		System.err.println("\nverticalDivide (depth " + depth + ")");
		System.err.println("Region: " + r);
		System.err.println("Pivots: " + f_pivot + " - " + l_pivot);
		System.err.println("Probes: " + f_probe + " - " + l_probe);
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
		part_pivot = choosePivotPair (f_pivot, l_pivot);
		 
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
		// System.err.println("Partitioning: " + (100 * cols_1 / num_cols) + " x " + (100 * cols_2 / num_cols));
		
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
	
	private int choosePivotPair (int first, int last)
	{
		int i, j, p1, p2, dist, maxdist, tmp;
		
		p1 = first;
		p2 = last;
		maxdist = -1;
		
		// find pair of pivots p1 and p2 with maximum Hamming distance
		for (i = first; i < last; i++)
			for (j = i + 1; j <= last; j++)
			{
				dist = LayoutEvaluation.hammingDistance(chip, pid[i], pid[j]);
				if (dist > maxdist)
				{
					maxdist = dist;
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
		
		// TODO remove this
		System.err.println("Pivot IDs: " + p1 + ", " + p2);
		
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
		
		/*
		System.err.println(count1 + " pivots to p1");
		System.err.println(count2 + " pivots to p2");
		//*/
		
		// TODO remove this
		if (pid[last] != p2)
			throw new IllegalStateException ("Impossible");
			
		// move p2 to the beginning of its own list
		pid[last] = pid[i];
		pid[i] = p2;

		return i;
	}

	private int divideProbes (int p1, int p2, int first, int last)
	{
		int		i, count1 = 0, count2 = 0, tmp;
		double	dist1, dist2;
		
		for (i = first; i <= last;)
		{
			dist1 = ospe.minDistanceProbe(pid[i], p1);
			dist2 = ospe.minDistanceProbe(pid[i], p2);
			
			/*
			// TODO remove this
			if (d1 + 4 < d2)
				dwin1++;
			else if (d1 < d2)
				win1++;
			else if (d1 == d2)
				draw++;
			else if (d1 > d2 + 4)
				dwin2++;
			else if (d1 > d2)
				win2++;
			else
				throw new IllegalStateException ("What??");
			//*/

			if ((dist1 < dist2) || (dist1 == dist2 && count1 < count2))
			{
				// assign probe to p1
				i++;
				count1++;
			}
			else
			{
				// assign probe to p2
				tmp = pid[i];
				pid[i] = pid[last];
				pid[last] = tmp;
				last--;
				count2++;
			}
		}
		
		// TODO remove this
		/*
		System.err.println(dwin1 + "," + win1 + "," + draw + "," + win2 + "," + dwin2);
		System.err.println(count1 + " probes to p1, " + count2 + " probes to p2");
		//*/
		
		return i;
	}

	private int fillRegion (RectangularRegion region, int f_pivot, int l_pivot,
			int f_probe, int l_probe)
	{
		int i, pivot_id, unplaced;

		// ID of main pivot
		pivot_id = pid[f_pivot];
		
		// reembed pivots optimally in regards to the main pivot
		for (i = f_pivot + 1; i <= l_pivot; i++)
			ospe.reembedProbe(pid[i], pivot_id);

		// reembed probes optimally in regards to the main pivot
		for (i = f_probe; i <= l_probe; i++)
			ospe.reembedProbe(pid[i], pivot_id);
			// TODO remove this
			//               (pid[i], pid, f_pivot, l_pivot);
		
		// place pivots
		unplaced = filler.fillRegion(chip, region, pid, f_pivot, l_pivot);

		// place non-pivots
		unplaced += filler.fillRegion(chip, region, pid, f_probe, l_probe);

		return unplaced;
	}

	// TODO remove this
	/* private */ void assignProbes (int pivots)
	{
		int best, count[], assign[];
		double dist, mindist;
		
		System.err.println("Assigning probes to pivots...");
		
		count = new int[pivots];
		assign = new int[pid.length - pivots];
		
		for (int i = pivots; i < pid.length; i++)
		{
			if (i % 1000 == 0)
				System.err.println(i);
			
			best = 0;
			mindist = ospe.minDistanceProbe(pid[i], pid[best]);
			
			for (int j = 1; j < pivots; j++)
			{
				dist = ospe.minDistanceProbe(pid[i], pid[j]);
				
				if (dist > mindist)
					continue;
				
				if (dist == mindist)
					if (count[j] >= count[best])
						continue;
				
				mindist = dist;
				best = j;
			}
			
			// assign pid i to best pivot
			assign[i - pivots] = best;
			count[best]++;
		}
		
		System.exit(1);
	}
}
