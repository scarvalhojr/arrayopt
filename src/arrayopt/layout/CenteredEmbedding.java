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
 * This class contains centered algorithm for several microarray chips
 *
 * @author Anna Domanski & Ronny Gärtner
 */
public class CenteredEmbedding extends ProbeEmbeddingAlgorithm
{
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
	 */
	public void reembedProbe (SimpleChip chip, int probe_id)
	{
		// to do:
		LeftMostEmbedding embedder = new LeftMostEmbedding();
		boolean isshift = false;
		int shift = 0, spaces = 0, mask = 1; 
		int int_index = chip.embed[probe_id].length;
		int pos = chip.embed_len;	
	
		embedder.reembedProbe(chip,probe_id);
		while(!isshift)
		{
			if ((chip.embed[probe_id][int_index] & mask) !=0)
				isshift = true;
			else 
			{
				spaces++;
				pos--;
				mask <<= 1;
				if (pos % Integer.SIZE == 0)
				{
					int_index --;
					mask = 1;
				}
			}
		}
		//shift a whole cycle
		shift = (int) Math.floor((double)((spaces-2)/4));	
		embedder.reembedProbe(chip, probe_id, shift);
		// do a left-most embedding and check how many masked steps
		// are left after the last synthesized probe, call it spaces

		// do a left-most embedding with a shift = spaces / 2
	}

	/**
	 * reembeds a probe of a Affymetrix Chip in such a way that the bases of the probe will be aligned 
	 * in more to the center of the deposition sequence
	 * 'pair-awareness' is included
	 */
	public void reembedProbe (AffymetrixChip chip, int probe_id)
	{
		// to do:
		LeftMostEmbedding embedder = new LeftMostEmbedding();
		boolean isshift = false;
		int shift = 0, spaces = 0, mask = 1; 
		int int_index = chip.embed[probe_id].length;
		int pos = chip.embed_len;	
	
		embedder.reembedProbe(chip,probe_id);
		while(!isshift)
		{
			if ((chip.embed[probe_id][int_index] & mask) !=0)
				isshift = true;
			else 
			{
				spaces++;
				pos--;
				mask <<= 1;
				if (pos % Integer.SIZE == 0)
				{
					int_index --;
					mask = 1;
				}
			}
		}
		//shift a whole cycle
		shift = (int) Math.floor((double)((spaces-2)/4));	
		embedder.reembedProbe(chip, probe_id, shift);
		// do a left-most embedding and check how many masked steps
		// are left after the last synthesized probe, call it spaces

		// do a left-most embedding with a shift = spaces / 2
	}
}
