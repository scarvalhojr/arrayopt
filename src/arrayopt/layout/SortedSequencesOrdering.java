/*
 * SortedSequencesOrdering.java
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
 *
 */
public class SortedSequencesOrdering implements ProbeOrderingAlgorithm	
{
	/**
	 *
	 */
	public void orderProbes (Chip chip, int[] id)
	{
		orderProbes (chip, id, 0, id.length - 1);
	}

	/**
	 *
	 */
	public void orderProbes (Chip chip, int[] id, int start, int end)
	{
		long probe_rank[] = chip.computeProbeRanks(id, 0, id.length - 1);
		QuickSort.sort(new SequenceSorting(id, probe_rank), start,
				end - start + 1);
	}
	
	private class SequenceSorting implements ArrayIndexedCollection
	{
		private int probe_id[];
		
		private long rank[];
		
		private long pivot;
		
		SequenceSorting (int probe_id[], long probe_rank[])
		{
			this.probe_id = probe_id;
			this.rank = probe_rank;
		}
		
		public int compare (int i, int j)
		{
			return rank[i] < rank[j] ? -1 :
					rank[i] == rank[j] ? 0 : +1;
		}
		
		public void swap (int i, int j)
		{
			int tmp1;
			tmp1 = probe_id[i];
			probe_id[i] = probe_id[j];
			probe_id[j] = tmp1;
			
			long tmp2;
			tmp2 = rank[i];
			rank[i] = rank[j];
			rank[j] = tmp2;
		}
		
		public void setPivot (int i)
		{
			this.pivot = rank[i];
		}
		
		public int compareToPivot (int i)
		{
			return rank[i] < pivot ? -1 :
					rank[i] == pivot ? 0 : +1;
		}
	}
}
