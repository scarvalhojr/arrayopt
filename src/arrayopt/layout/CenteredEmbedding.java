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
 * This class contains centered algorithm for several microarray chips.
 * The goal of centered algorithm is to align a probe more in the center of the depostion sequence.
 * @author Anna Domanski & Ronny G?rtner
 */
public class CenteredEmbedding implements SingleProbeEmbeddingAlgorithm,
	ProbeSetEmbeddingAlgorithm
{
	private LeftMostEmbedding embedder = new LeftMostEmbedding();

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
	 * embeds a probe of a simple chip so that the bases of the probe will be 
	 * aligned more in the center of the deposition sequence
	 * @param chip simple chip
	 * @param probe_id probe ID
	 * @throws IllegalArgumentException if reembedding is not possible, e.g. because of former
	 * wrong embedding
	 */
	public void reembedProbe (SimpleChip chip, int probe_id)
	{
		// to do:
		int shift = 0, mask, int_index, pos, count_appearance = 0;
		
		pos = chip.dep_seq.length - 1;
		mask = 0x01 << ((Integer.SIZE-1) - (pos % Integer.SIZE));
		int_index = (int) Math.floor((double) pos/Integer.SIZE);

		// leftmost embedding with no shift
		embedder.reembedProbe(chip,probe_id,0);
		
		while((chip.embed[probe_id][int_index] & mask) == 0)
		{
			// throw exception if pos reaches area where bases must be for sure
			if (--pos < chip.probe_len-1)
			{
				throw new IllegalArgumentException ("Unable to reembed probe.");
			}
			if (pos%Integer.SIZE == 0)
			{
				int_index--;
				mask = 0x01;
			}
			else
				mask <<= 1;
		}
		
		// count appearance of last base in 'zero tail'
		for (int i = pos + 1; i < chip.dep_seq.length; i++)
		{
			if (chip.dep_seq[i] == chip.dep_seq[pos])
				count_appearance++;
		}
		//Attention! length of cycle is hardcoded --> here length of cycle = 4
		shift = ((int) Math.floor((double) count_appearance/2)) * 4;
	
		embedder.reembedProbe(chip, probe_id, shift);
		// do a left-most embedding and check how many masked steps
		// are left after the last synthesized probe, call it spaces

		// do a left-most embedding with a shift = spaces / 2
	}

	/**
	 * reembeds a probe of an affymetrix chip in such a way that the bases of the probe will be 
	 * aligned more to the center of the deposition sequence
	 * 'pair-awareness' is included
	 * @param chip affymetrix chip
	 * @param probe_id probe ID
	 * @throws IllegalArgumentException if reembedding is not possible, e.g. because of former
	 * wrong embedding
	 */
	public void reembedProbe (AffymetrixChip chip, int probe_id)
	{
		int shift = 0, mask, int_index, pos, count_appearance = 0;
		
		pos = chip.dep_seq.length - 1;
		mask = 0x01 << ((Integer.SIZE-1) - (pos % Integer.SIZE));
		int_index = (int) Math.floor((double) pos/Integer.SIZE);

		// leftmost embedding with no shift
		embedder.reembedProbe(chip,probe_id,0);
		
		while((chip.embed[probe_id][int_index] & mask) == 0)
		{
			// throw exception if pos reaches area where bases must be for sure
			if (--pos < chip.probe_len-1)
			{
				throw new IllegalArgumentException ("Unable to reembed probe.");
			}
			if (pos%Integer.SIZE == 0)
			{
				int_index--;
				mask = 0x01;
			}
			else
				mask <<= 1;
		}
		
		// count appearance of last base in 'zero tail'
		for (int i = pos + 1; i < chip.dep_seq.length; i++)
		{
			if (chip.dep_seq[i] == chip.dep_seq[pos])
				count_appearance++;
		}
		//Attention! length of cycle is hardcoded --> here length of cycle = 4
		shift = ((int) Math.floor((double) count_appearance/2)) * 4;
	
		embedder.reembedProbe(chip, probe_id, shift);
		// do a left-most embedding and check how many masked steps
		// are left after the last synthesized probe, call it spaces

		// do a left-most embedding with a shift = spaces / 2
	}
}
