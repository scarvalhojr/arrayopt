/*
 * SNPRightmostEmbedding.java
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
public class SNPRightmostEmbedding implements LayoutAlgorithm
{
	private RightMostEmbedding right_emb = new RightMostEmbedding();
	
	private SNPLeftmostEmbedding snpleft_emb = new SNPLeftmostEmbedding();
	
	/**
	 *
	 */
	public void changeLayout (Chip chip)
	{
		BitSet reembed;
		int num_rows, num_cols, id1, id2, embed_len, cycle_len;
		
		num_rows = chip.getNumberOfRows();
		num_cols = chip.getNumberOfColumns();
		embed_len = chip.getEmbeddingLength();
		cycle_len = chip.depositionSequenceCycleLength();
		
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
							reembedPair (chip, embed_len, cycle_len, id1, id2);
							reembed.set(r * num_cols + c);
							reembed.set(r * num_cols + (c + 1));
							continue;
						}
				
				if (r < num_rows - 1)
					if ((id2 = chip.spot[r + 1][c]) != Chip.EMPTY_SPOT)
						if (chip.isSNPProbePair(id1, id2))
						{
							reembedPair (chip, embed_len, cycle_len, id1, id2);
							reembed.set(r * num_cols + c);
							reembed.set((r + 1) * num_cols + c);
							continue;
						}
				
				// if probe is not adjacent to a SNP mate,
				// it is simply right-most embedded
				right_emb.reembedProbe(chip, id1, 0);
			}
	}
	
	private void reembedPair (Chip chip, int embed_len, int cycle_len,
			int id1, int id2)
	{
		int shift, pos, w, bitmask = 0;
		
		// TODO remove this
		/*
		System.err.println("Probe " + id1);
		chip.printEmbedding(id1);
		System.err.println("\nProbe " + id2);
		chip.printEmbedding(id2);
		System.err.println();
		//*/
		
		// right-most re-embed both probes (without shifts)
		right_emb.reembedProbe(chip, id1, 0);
		right_emb.reembedProbe(chip, id2, 0);

		// TODO remove this
		/*
		System.err.println("\nRight-most embeddings:");
		System.err.println("Probe " + id1);
		chip.printEmbedding(id1);
		System.err.println("\nProbe " + id2);
		chip.printEmbedding(id2);
		System.err.println();
		//*/

		// check max number of free steps
		// (masked steps before first productive step)
		for (pos = -1, w = -1; pos < embed_len;)
		{
			if ((++pos % Integer.SIZE) == 0)
			{
				w++;
				bitmask = 0x01 << (Integer.SIZE - 1);
			}
			else
				bitmask >>>= 1; 
			
			if (((chip.embed[id1][w] & bitmask) != 0) ||
				((chip.embed[id2][w] & bitmask) != 0))
				break;
		}
		
		// compute shift
		shift = cycle_len * (int) Math.floor(pos / (double) cycle_len);
		
		// left-most re-embed probe pair with the calculated shift
		snpleft_emb.reembedPair(chip, embed_len, id1, id2, shift);
		
		// TODO remove this
		/*
		System.err.println("\nSNP right-most embeddings (shift = " + shift + "):");
		System.err.println("Probe " + id1);
		chip.printEmbedding(id1);
		System.err.println("\nProbe " + id2);
		chip.printEmbedding(id2);
		System.err.println();
		//*/
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
