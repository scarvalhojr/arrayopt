/*
 * SNPPairRightmostEmbedding.java
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

import java.util.BitSet;

/**
 *
 * Find SNP pairs that are on adjacent spots and re-embedded them with a
 * pair-wise right-most synthesis strategy. If a probe is not adjacent to its
 * SNP mate, it is simply right-most embedded (not pair-wise).  
 */
public class SNPPairRightmostEmbedding implements LayoutAlgorithm
{
	private RightMostEmbedding right_emb = new RightMostEmbedding();
	
	/**
	 *
	 */
	public void changeLayout (Chip chip)
	{
		BitSet reembed;
		int num_rows, num_cols, embed_len, id1, id2;
		
		num_rows = chip.getNumberOfRows();
		num_cols = chip.getNumberOfColumns();
		embed_len = chip.getEmbeddingLength();
		
		reembed = new BitSet(num_rows * num_cols);
		 
		for (int r = 0; r < num_rows; r++)
			for (int c = 0; c < num_cols; c++)			
			{
				if ((id1 = chip.spot[r][c]) == Chip.EMPTY_SPOT)
					continue;
				
				if (reembed.get(r * num_cols + c))
					continue;
				
				if (c < num_cols - 1)
					if ((id2 = chip.spot[r][c + 1]) != Chip.EMPTY_SPOT)
						if (chip.isSNPProbePair(id1, id2))
						{
							reembedPair (chip, embed_len, id1, id2);
							reembed.set(r * num_cols + c);
							reembed.set(r * num_cols + (c + 1));
							continue;
						}
				
				if (r < num_rows - 1)
					if ((id2 = chip.spot[r + 1][c]) != Chip.EMPTY_SPOT)
						if (chip.isSNPProbePair(id1, id2))
						{
							reembedPair (chip, embed_len, id1, id2);
							reembed.set(r * num_cols + c);
							reembed.set((r + 1) * num_cols + c);
							continue;
						}
				
				// if probe is not adjacent to a SNP mate,
				// it is simply right-most embedded
				right_emb.reembedProbe(chip, id1, 0);
			}
	}
	
	private void reembedPair (Chip chip, int embed_len, int id1, int id2)
	{
		int pos1, pos2, w1, w2, bitmask1, bitmask2;
		int master_id, slave_id, master_pos, pos, w, bitmask;
		boolean mismatch = false, stop = false;
		
		// right-most re-embed both probes (without shifts)
		right_emb.reembedProbe(chip, id1, 0);
		right_emb.reembedProbe(chip, id2, 0);
		
		// proceed from right to left until mismatch is found
		pos1 = pos2 = embed_len;
		w1 = w2 = embed_len / Integer.SIZE;
		bitmask1 = bitmask2 = 0x01 << (Integer.SIZE -1 - (pos1 % Integer.SIZE));
		while (pos1 > 0 && !stop)
		{
			if ((pos1-- % Integer.SIZE) == 0)
			{
				w1--;
				bitmask1 = 0x01;
			}
			else
				bitmask1 <<= 1; 
			
			if ((chip.embed[id1][w1] & bitmask1) != 0)
			{
				while (pos2 > 0)
				{
					if ((pos2-- % Integer.SIZE) == 0)
					{
						w2--;
						bitmask2 = 0x01;
					}
					else
						bitmask2 <<= 1;
					
					if ((chip.embed[id2][w2] & bitmask2) != 0)
					{
						if (!mismatch)
						{
							// signal when mismatch position is found
							if (pos1 != pos2) mismatch = true;
						}
						else // (mismatch)
						{
							// if right-most productive step of the left part
							// is aligned, then all left part is aligned since
							// probes were right-most embedded beforehand:
							// nothing to do!
							if (pos1 == pos2)
								return;
							
							// otherwise, signal to quit "search loop"
							stop = true;
						}
						
						break;
					}
				}
			}
		}
		
		if (pos1 < pos2)
		{
			// left part of probe 1's embedding will be copied to probe 2 
			master_id = id1;
			master_pos = pos1;
			slave_id = id2;
			pos = pos2;
			w = w2;
			bitmask = bitmask2;
		}
		else // (pos2 < pos1)
		{
			// left part of probe 2's embedding will be copied to probe 1
			master_id = id2;
			master_pos = pos2;
			slave_id = id1;
			pos = pos1;
			w = w1;
			bitmask = bitmask1;
		}
		
		// first, erase rightmost of misaligned steps
		while (pos > master_pos)
		{
			chip.embed[slave_id][w] = chip.embed[slave_id][w] & ~bitmask;

			if ((pos-- % Integer.SIZE) == 0)
			{
				w--;
				bitmask = 0x01;
			}
			else
				bitmask <<= 1;
		}

		// copy all embedding steps from master to slave, until position zero
		while (pos >= 0)
		{
			if ((chip.embed[master_id][w] & bitmask) != 0)
				chip.embed[slave_id][w] = chip.embed[slave_id][w] | bitmask;
			else
				chip.embed[slave_id][w] = chip.embed[slave_id][w] & ~bitmask;

			if ((pos-- % Integer.SIZE) == 0)
			{
				w--;
				bitmask = 0x01;
			}
			else
				bitmask <<= 1;
		}
	}
	
	/**
	 * Returns the algorithm's name.
	 */
	@Override
	public String toString ()
	{
		return this.getClass().getSimpleName();
	}
}
