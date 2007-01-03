/*
 * AlignedEmbedding.java
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
 * This class implements the aligned re-embedding algorithm. Probes are
 * re-embedded in such a way that the middle bases are aligned with the central
 * masks as far as possible. Such embeddings are good when combined with a
 * placement algorithm that optimizes the central masks such as the
 * {@link TwoDimensionalCentralPartitioning} because the middle bases will be
 * synthesized by masks that are highly optimized. 
 *  
 * @author Sergio A. de Carvalho Jr.
 */
public class AlignedEmbedding implements SingleProbeEmbeddingAlgorithm,
	LayoutAlgorithm
{
	private LeftMostEmbedding embedder = new LeftMostEmbedding();
	
	/**
	 * Creates a new layout of a microarray chip by re-embedding all probes
	 * using the aligned re-embedding algorithm.
	 * 
	 * @param chip instance of a microarray chip
	 */
	public void changeLayout (Chip chip)
	{
		reembedProbeSet (chip, chip.getAllProbes());
	}

	/**
	 * Re-embeds a set of probes of a given chip.
	 * 
	 * @param chip the chip containing the probes to be re-embedded
	 * @param id a list of probe IDs to be re-embedded
	 */
	public void reembedProbeSet (Chip chip, int id[])
	{
		reembedProbeSet (chip, id, 0, id.length - 1);
	}

	/**
	 * Re-embeds a selected sub-set of probes of a given chip. 
	 * 
	 * @param chip the chip containing the probes to be re-embedded
	 * @param id a list of probe IDs to be re-embedded
	 * @param first first element on the list to be re-embedded
	 * @param last last element on the list to be re-embedded
	 */
	public void reembedProbeSet (Chip chip, int id[], int first, int last)
	{
		SimpleChip schip;
		int embed_len, cycle, min_step, middle;
		
		if (chip instanceof SimpleChip)
			schip = (SimpleChip) chip;
		else
			throw new IllegalArgumentException ("Unsupported chip type.");
		
		embed_len = chip.getEmbeddingLength();
		cycle = chip.depositionSequenceCycleLength();
		
		if (cycle <= 0)
			throw new IllegalStateException
				("Deposition sequence must be cyclical.");
		
		// first step of the central masks
		min_step = (embed_len / 2) - (cycle / 2);
		
		// middle base
		middle = (int) Math.round(chip.getProbeLength() / (double) 2);
		
		for (int i = first; i <= last; i++)
			reembedProbe (schip, id[i], embed_len, cycle, middle, min_step);
	}
	
	/**
	 * Re-embeds a single probe of a given chip. 
	 * 
	 * @param chip the chip containing the probe to be re-embedded
	 * @param probe_id the ID of the probe to be re-embedded
	 */
	public void reembedProbe (Chip chip, int probe_id)
	{
		if (chip instanceof SimpleChip)
			reembedProbe ((SimpleChip) chip, probe_id);
		else
			throw new IllegalArgumentException ("Unsupported chip type.");
	}

	/**
	 * Re-embeds a single probe of a given chip. 
	 * 
	 * @param chip the chip containing the probe to be re-embedded
	 * @param probe_id the ID of the probe to be re-embedded
	 */
	public void reembedProbe (SimpleChip chip, int probe_id)
	{
		int embed_len, cycle, middle, min_step;
		
		embed_len = chip.getEmbeddingLength();
		cycle = chip.depositionSequenceCycleLength();
		
		if (cycle <= 0)
			throw new IllegalStateException
				("Deposition sequence must be cyclical.");
		
		// first step of the central masks
		min_step = (embed_len / 2) - (cycle / 2);
		
		// middle base
		middle = (int) Math.round(chip.getProbeLength() / (double) 2);
		
		reembedProbe (chip, probe_id, embed_len, cycle, middle, min_step);
	}
	
	private void reembedProbe (SimpleChip chip, int id, int embed_len,
			int cycle_len, int middle_base, int min_step)
	{
		int free, max_shift, shift;
		int base, pos, w, bitmask = 0;
		
		// leftmost embedding with no shift
		embedder.reembedProbe(chip, id, 0);
		
		// count number of masked steps after last productive step
		pos = embed_len - 1;
		w = (int) Math.floor(pos / (double) Integer.SIZE);
		bitmask = 0x01 << (Integer.SIZE - (embed_len % Integer.SIZE));
		for (free = 0; pos >= 0; pos--)
		{
			if ((chip.embed[id][w] & bitmask) == 0)
				free++;
			else
				break;
			
			if ((pos % Integer.SIZE) == 0)
			{
				bitmask = 0x01;
				w--;
			}
			else
				bitmask <<= 1;
		}
		
		// compute maximum shift
		max_shift = 4 * (int) Math.floor(free / (double) cycle_len);

		// find the middle base's current position
		pos = embed_len - 1;
		bitmask = 0x01 << (Integer.SIZE - (embed_len % Integer.SIZE));
		for (w = -1, pos = -1, base = 0; base < middle_base;)
		{
			if ((++pos % Integer.SIZE) == 0)
			{
				bitmask = 0x01 << (Integer.SIZE - 1);
				w++;
			}
			else
				bitmask >>>= 1;
			
			if ((chip.embed[id][w] & bitmask) != 0)
				base++;
		}
		
		if (pos >= min_step)
			// with a leftmost embedding, the middle base is already past
			// the central masks => it is already aligned or cannot be aligned
			return;
		
		// compute optimal shift
		shift = 4 * (int) Math.ceil((min_step - pos) / (double) cycle_len);

		// make sure optimal shift is possible 
		shift = shift > max_shift ? max_shift : shift;
		
		// re-embed with shift
		if (shift > 0) embedder.reembedProbe(chip, id, shift);
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
