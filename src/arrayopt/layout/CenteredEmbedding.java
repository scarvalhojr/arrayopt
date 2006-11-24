/*
 * CenteredEmbedding.java
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
 * ?rrayOpt is distributed in the hope that it will be useful, but WITHOUT ANY
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
 * This class implements the centered re-embedding algorithm. Probes are
 * re-embedded in such a way that the number of masked steps to the left of the
 * first productive step is approximately equal to the number of masked steps to
 * the right of the last productive step. This results in embeddings that
 * reduce the number of unmasked spots on the first and last masks. 
 *  
 * @author Anna Domanski & Ronny G?rtner
 * @author Sergio A. de Carvalho Jr.
 */
public class CenteredEmbedding implements SingleProbeEmbeddingAlgorithm,
	ProbeSetEmbeddingAlgorithm
{
	private LeftMostEmbedding embedder = new LeftMostEmbedding();
	
	/**
	 * Constant used internally to alternate between left- and right-most
	 * preference when both result in the embedding having the same difference
	 * of masked steps before the first productive step and after the last
	 * productive step. 
	 */
	private int alternate = - 1;

	/**
	 * embedd a complete set of probes of a given chip
	 */
	public void reembedProbeSet (Chip chip, int probe_id[])
	{
		reembedProbeSet (chip, probe_id, 0, probe_id.length - 1);
	}

	/**
	 * embedd a set of probes of a chip within a range: from first to last probe
	 */
	public void reembedProbeSet (Chip chip, int probe_id[], int first, int last)
	{
		for (int i = first; i <= last; i++)
			reembedProbe (chip, probe_id[i]);
	}
	
	/**
	 * embeds a single probe of a chip
	 */
	public void reembedProbe (Chip chip, int probe_id)
	{
		if (chip instanceof SimpleChip)
			reembedProbe ((SimpleChip) chip, probe_id);

		else if (chip instanceof AffymetrixChip)
			reembedProbe ((AffymetrixChip) chip, probe_id);

		else
			throw new IllegalArgumentException ("Unsupported chip type.");
	}

	/**
	 * Re-embeds a single probe in a centered fashion.
	 *  
	 * @param chip chip
	 * @param probe_id probe ID
	 */
	public void reembedProbe (SimpleChip chip, int probe_id)
	{
		int cycle_len, left, right, shift;
		int pos, idx, bitmask = 0;
		double s;

		if ((cycle_len = chip.depositionSequenceCycleLength()) <= 0)
			throw new IllegalArgumentException
				("Only chips with cyclical deposition sequences are supported.");
		
		// leftmost embedding with no shift
		embedder.reembedProbe(chip,probe_id,0);
		
		// count number of masked steps before first productive step
		for (left = 0, idx = -1, pos = 0; pos < chip.embed_len; pos++)
		{
			if ((pos % Integer.SIZE) == 0)
			{
				bitmask = 0x01 << (Integer.SIZE - 1);
				idx++;
			}
			else
				bitmask >>>= 1;
				
			if ((chip.embed[probe_id][idx] & bitmask) == 0)
				left++;
			else
				break;
		}
		
		// count number of masked steps after last productive step
		pos = chip.dep_seq.length - 1;
		idx = (int) Math.floor(pos / (double) Integer.SIZE);
		bitmask = 0x01 << (Integer.SIZE - (chip.dep_seq.length % Integer.SIZE));
		for (right = 0; pos >= 0; pos--)
		{
			if ((chip.embed[probe_id][idx] & bitmask) == 0)
				right++;
			else
				break;
			
			if ((pos % Integer.SIZE) == 0)
			{
				bitmask = 0x01;
				idx--;
			}
			else
				bitmask <<= 1;
		}

		// compute best shift
		if (right > left)
		{
			s = (right - left) / (float) (2 * cycle_len);
			
			if (((right - left) % (2 * cycle_len)) == cycle_len)
			{
				if ((alternate = - alternate) == 1)
					shift = (int) Math.floor(s);
				else
					shift = (int) Math.ceil(s);
			}
			else
				shift = (int) Math.round(s);

			embedder.reembedProbe(chip, probe_id, shift * cycle_len);
		}
	}

	public void reembedProbe (AffymetrixChip chip, int probe_id)
	{
		// TODO implement based on the SimpleChip code
		return;
	}
}
