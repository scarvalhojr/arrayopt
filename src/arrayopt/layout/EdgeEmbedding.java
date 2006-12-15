/*
 * EdgeEmbedding.java
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
 */
public class EdgeEmbedding implements SingleProbeEmbeddingAlgorithm,
	LayoutAlgorithm
{
	private LeftMostEmbedding left_embed = new LeftMostEmbedding();
	
	private int embed[];
	
	public void changeLayout (Chip chip)
	{
		reembedProbeSet (chip, chip.getAllProbes());
	}

	public void reembedProbeSet (Chip chip, int probe_id[])
	{
		reembedProbeSet (chip, probe_id, 0, probe_id.length - 1);
	}

	public void reembedProbeSet (Chip chip, int probe_id[], int first, int last)
	{
		if (embed == null)
			embed = new int [chip.embed[0].length];
		else if (embed.length != chip.embed[0].length)
			embed = new int [chip.embed[0].length];
		
		if (chip instanceof SimpleChip)
		{
			for (int i = first; i <= last; i++)
				reembed ((SimpleChip) chip, probe_id[i]);
		}
		else if (chip instanceof AffymetrixChip)
		{
			for (int i = first; i <= last; i++)
				reembed ((AffymetrixChip) chip, probe_id[i]);
		}
		else
			throw new IllegalArgumentException ("Unsupported chip type.");
	}
	
	public void reembedProbe (Chip chip, int probe_id)
	{
		if (embed == null)
			embed = new int [chip.embed[0].length];
		else if (embed.length != chip.embed[0].length)
			embed = new int [chip.embed[0].length];

		if (chip instanceof SimpleChip)
			reembed ((SimpleChip) chip, probe_id);
		else if (chip instanceof AffymetrixChip)
			reembed ((AffymetrixChip) chip, probe_id);
		else
			throw new IllegalArgumentException ("Unsupported chip type.");
	}

	public void reembedProbe (SimpleChip chip, int probe_id)
	{
		if (embed == null)
			embed = new int [chip.embed[0].length];
		else if (embed.length != chip.embed[0].length)
			embed = new int [chip.embed[0].length];
		
		reembed (chip, probe_id);
	}

	public void reembedProbe (AffymetrixChip chip, int probe_id)
	{
		if (embed == null)
			embed = new int [chip.embed[0].length];
		else if (embed.length != chip.embed[0].length)
			embed = new int [chip.embed[0].length];
		
		reembed (chip, probe_id);
	}

	private void reembed (SimpleChip chip, int id)
	{
		int embed_len, step, rstep, mask, rmask, w, rw, len_right = 0;
		
		// first, a left-most embedding (with no shift)
		left_embed.reembedProbe(chip, id, 0);
		
		embed_len = chip.getEmbeddingLength();
		
		rstep = step = embed_len - 1;
		rmask = mask = 0x01 << Integer.SIZE - (embed_len % Integer.SIZE); 
		rw = w = chip.embed[id].length - 1;
		
		// find last unmasked step
		while ((chip.embed[id][w] & mask) == 0)
		{
			if ((step-- % Integer.SIZE) == 0)
			{
				w--;
				mask = 0x01;
			}
			else
				mask <<= 1;
		}
		
		// find end of right part
		while (chip.dep_seq[rstep] != chip.dep_seq[step])
		{
			if ((rstep-- % Integer.SIZE) == 0)
			{
				rw--;
				rmask = 0x01;
			}
			else
				rmask <<= 1;
		}
		
		// TODO this is not optimal!
		
		do
		{
			if ((chip.embed[id][w] & mask) != 0)
			{
				chip.embed[id][w] &= ~mask;
				chip.embed[id][rw] |= rmask;
				len_right = embed_len - rstep;
			}

			if ((step-- % Integer.SIZE) == 0)
			{
				w--;
				mask = 0x01;
			}
			else
				mask <<= 1;

			if ((rstep-- % Integer.SIZE) == 0)
			{
				rw--;
				rmask = 0x01;
			}
			else
				rmask <<= 1;
			
		} while (len_right < step);
	}

	private void reembed (AffymetrixChip chip, int id)
	{
		int embed_len, step, rstep, mask, rmask, w, rw, len_right = 0;
		
		// first, a left-most embedding (with no shift)
		left_embed.reembedProbe(chip, id, 0);
		
		embed_len = chip.getEmbeddingLength();
		
		rstep = step = embed_len - 1;
		rmask = mask = 0x01 << Integer.SIZE - (embed_len % Integer.SIZE); 
		rw = w = chip.embed[id].length - 1;
		
		// find last unmasked step
		while ((chip.embed[id][w] & mask) == 0)
		{
			if ((step-- % Integer.SIZE) == 0)
			{
				w--;
				mask = 0x01;
			}
			else
				mask <<= 1;
		}
		
		// find end of right part
		while (chip.dep_seq[rstep] != chip.dep_seq[step])
		{
			if ((rstep-- % Integer.SIZE) == 0)
			{
				rw--;
				rmask = 0x01;
			}
			else
				rmask <<= 1;
		}
		
		// TODO this is not optimal!
		
		do
		{
			if ((chip.embed[id][w] & mask) != 0)
			{
				chip.embed[id][w] &= ~mask;
				chip.embed[id][rw] |= rmask;
				len_right = embed_len - rstep;
			}

			if ((chip.embed[id + 1][w] & mask) != 0)
			{
				chip.embed[id + 1][w] &= ~mask;
				chip.embed[id + 1][rw] |= rmask;
				len_right = embed_len - rstep;
			}

			if ((step-- % Integer.SIZE) == 0)
			{
				w--;
				mask = 0x01;
			}
			else
				mask <<= 1;

			if ((rstep-- % Integer.SIZE) == 0)
			{
				rw--;
				rmask = 0x01;
			}
			else
				rmask <<= 1;
			
		} while (len_right < step);
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
