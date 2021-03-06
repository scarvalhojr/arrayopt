/*
 * SortedEmbeddingsOrdering.java
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
 * �rrayOpt is distributed in the hope that it will be useful, but WITHOUT ANY
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
public class SortedEmbeddingsOrdering implements ProbeOrderingAlgorithm	
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
		QuickSort.sort(new EmbeddingSort(chip, id), start, end - start + 1);
	}
	
	private class EmbeddingSort implements ArrayIndexedCollection
	{
		private Chip chip;
		
		private int words;
		
		private int probe_id[];
		
		private int pivot;
		
		EmbeddingSort (Chip chip, int probe_id[])
		{
			this.chip = chip;
			this.words = chip.embed[0].length;
			this.probe_id = probe_id;
		}
		
		public int compare (int i, int j)
		{
			int id_i = probe_id[i];
			int id_j = probe_id[j];
			
			for (int w = 0; w < words; w++)
			{
				// compare first bit (signal)
				// first bit is 1 => negative number
				// first bit is 0 => non-negative
				if (chip.embed[id_i][w] < 0)
				{
					if (chip.embed[id_j][w] >= 0)
						return 1;
				}
				else
				{
					if (chip.embed[id_j][w] < 0)
						return -1;
				}
				
				// compare remaining bits if both have same signal
				if (chip.embed[id_i][w] < chip.embed[id_j][w])
					return -1;
				else if (chip.embed[id_i][w] > chip.embed[id_j][w])
					return 1;
			}
			
			return 0;
		}
		
		public void swap (int i, int j)
		{
			int tmp;
			tmp = probe_id[i];
			probe_id[i] = probe_id[j];
			probe_id[j] = tmp;
		}
		
		public void setPivot (int i)
		{
			this.pivot = probe_id[i];
		}
		
		public int compareToPivot (int i)
		{
			int id_i = probe_id[i];

			for (int w = 0; w < words; w++)
			{
				// compare first bit (signal)
				// first bit is 1 => negative number
				// first bit is 0 => non-negative
				if (chip.embed[id_i][w] < 0)
				{
					if (chip.embed[pivot][w] >= 0)
						return 1;
				}
				else
				{
					if (chip.embed[pivot][w] < 0)
						return -1;
				}
				
				// compare remaining bits if both have same signal
				if (chip.embed[id_i][w] < chip.embed[pivot][w])
					return -1;
				else if (chip.embed[id_i][w] > chip.embed[pivot][w])
					return 1;
			}
			
			return 0;
		}
	}
}
