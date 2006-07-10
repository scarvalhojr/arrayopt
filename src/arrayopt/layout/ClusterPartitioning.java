/*
 * ClusterPartitioning.java
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
 * @author Sergio A. de Carvalho Jr.
 */
public class ClusterPartitioning implements PlacementAlgorithm
{
	private int mode;
	
	public static final int MODE_BORDER_LENGTH = 0;
	
	public static final int MODE_CONFLICT_INDEX = 1;

	/**
	 * TODO document this
	 */
	private FillingAlgorithm filler;

	/**
	 * TODO document this
	 */
	public static final double MIN_PERCENTAGE_PIVOTS = 0.01;

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
	private boolean fake_pivots;

	/**
	 * TODO document this
	 * 
	 * This is used to sort probes lexicographically but also to store the index
	 * to the best pivot of each non-pivot
	 */
	private long rank[];

	/**
	 * TODO document this
	 */
	private double dist[];

	/**
	 * TODO document this
	 */
	private int offset;

	/**
	 * TODO document this
	 */
	private RankSorting rank_sort;

	/**
	 * TODO document this
	 */
	private int rows_per_probe;
	
	/**
	 * TODO document this
	 */
	public ClusterPartitioning (FillingAlgorithm filler)
	{
		this(filler, MODE_BORDER_LENGTH);
	}
	
	// TODO remove this
	private boolean sortprobes;

	/**
	 * TODO document this
	 */
	public ClusterPartitioning (FillingAlgorithm filler, int mode)
	{
		switch (mode)
		{
			case MODE_BORDER_LENGTH:
			case MODE_CONFLICT_INDEX:
				this.mode = mode;
				break;
				
			default:
				throw new IllegalArgumentException
					("Illegal value for argument 'mode': " + mode);
		}
		
		// TODO remove this!
		if (mode == MODE_BORDER_LENGTH)
			this.sortprobes = true;
		else
			this.sortprobes = false;
		this.mode = MODE_BORDER_LENGTH;
		
		this.filler = filler;
	}

	/**
	 * TODO document this
	 */
	public int makeLayout (Chip c)
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
		rank_sort = new RankSorting (pid, rank, offset);
		computeProbeRanks (offset, pid.length - 1);
		
		// TODO remove this if block
		if (sortprobes)
		{
			System.err.println("Sorting probes lexicographically!");
		
		// sort non-pivot probes lexicographically
		QuickSort.sort(rank_sort, offset, nonpivots);
		
		// TODO remove this if block
		} else {
			System.err.println("No Sorting!");
		}
		
		// create probe distance array
		dist = new double[nonpivots];
		for (int i = 0; i < dist.length; i++)
			dist[i] = Double.POSITIVE_INFINITY;
		
		makeClusters(0, pivots - 1, offset, pid.length - 1);
		
		return 0;
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
		
		// compute the minimum number of pivots
		min_pivots = (int) (MIN_PERCENTAGE_PIVOTS * chip.getNumberOfProbes());
		
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
			// turn on the fake pivots flag
			this.fake_pivots = true;
			
			// cut the excess of fake pivots
			if (num_pivots > min_pivots)
				num_pivots = min_pivots;
			
			// TODO remove this
			System.err.println("Pivot list trimmed to " + num_pivots);
			
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

	private void makeClusters (int f_pivot, int l_pivot, int f_probe, int l_probe)
	{
		double d;
		
		// System.err.println("Starting clustering...");
		
		for (int p = f_pivot; p <= l_pivot; p++)
		{
			d = ospe.minDistanceProbe(pid[f_probe], pid[p]);			

			/*
			System.err.println("\nPivot " + p + " (ID " + pid[p] + "): " + d);
			chip.printEmbedding(pid[p]);
			System.err.println();
			//*/
			
			if (d < dist[f_probe - offset])
			{
				/*
				System.err.println("New min dist for Probe ID " + pid[f_probe] + ":");
				ospe.reembedProbe(pid[f_probe]);
				chip.printEmbedding(pid[f_probe]);
				System.err.println();
				//*/
				
				dist[f_probe - offset] = d;
				rank[f_probe - offset] = p; 
			}
			
			for (int n = f_probe + 1; n <= l_probe; n++)
			{
				d = ospe.minDistanceProbe(pid[n]);
				if (d < dist[n - offset])
				{
					dist[n - offset] = d;
					rank[n - offset] = p; 
				}
			}
		}		
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
}
