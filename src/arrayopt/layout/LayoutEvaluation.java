/*
 * LayoutEvaluation.java
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
 *
 */
public class LayoutEvaluation
{
	// package access
	// distance-based weighting function
	static final double WEIGHT_DIST[][] = {
				{	0,		0,		0.1,	0.1111,	0.1,	0,		0		},
				{	0,		0.125,	0.2,	0.25,	0.2,	0.125,	0		},
				{	0.1,	0.2,	0.5,	1,		0.5,	0.2,	0.1		},
				{	0.1111,	0.25,	1,		0,		1,		0.25,	0.1111	},
				{	0.1,	0.2,	0.5,	1,		0.5,	0.2,	0.1		},
				{	0,		0.125,	0.2,	0.25,	0.2,	0.125,	0		},
				{	0,		0,		0.1,	0.1111,	0.1,	0,		0		}};

	// package access
	// compute the hamming distance between two probes' embeddings
	static int hammingDistance (Chip chip, int id_1, int id_2)
	{
		if (chip instanceof SimpleChip)
			return hammingDistance ((SimpleChip) chip, id_1, id_2);

		else if (chip instanceof AffymetrixChip)
			return hammingDistance ((AffymetrixChip) chip, id_1, id_2);

		else
			throw new IllegalArgumentException ("Unsupported chip type.");
	}

	// package access
	// compute the distance between two probes' embeddings with
	// distance-dependent weights (as defined for the conflict index)
	static double weightedDistance (Chip chip, int id_1, int id_2)
	{
		if (chip instanceof SimpleChip)
			return weightedDistance ((SimpleChip) chip, id_1, id_2);

		else if (chip instanceof AffymetrixChip)
			return weightedDistance ((AffymetrixChip) chip, id_1, id_2);

		else
			throw new IllegalArgumentException ("Unsupported chip type.");
	}

	/**
	 * TO DO: use Integer.bitCount(int i) instead of this ugly bit shifting stuff
	 */
	protected static int hammingDistance (SimpleChip chip, int id_1, int id_2)
	{
		char ch;
		int  bits = 0, hd = 0, mask = 0, w, pos;

		for (w = -1, pos = 0; pos < chip.getEmbeddingLength(); pos++)
		{
			if (pos % Integer.SIZE == 0)
			{
				// next 4-byte word
				w++;

				// xor
				bits = chip.embed[id_1][w] ^ chip.embed[id_2][w];

				// turn on very last bit of mask only
				mask = 0x01 << (Integer.SIZE - 1);
			}

			if ((mask & bits) != 0) hd++;

			// shift interrogating bit to the right
			// ('>>>' means a unsigned shift)
			mask >>>= 1;
		}

		return hd;
	}

	/**
	 * TO DO: use Integer.bitCount(int i) instead of this ugly bit shifting stuff
	 */
	protected static int hammingDistance (AffymetrixChip chip, int id_1, int id_2)
	{
		char ch;
		int  bits = 0, hd = 0, mask = 0, w, pos;

		for (w = -1, pos = 0; pos < chip.getEmbeddingLength(); pos++)
		{
			if (pos % Integer.SIZE == 0)
			{
				// next 4-byte word
				w++;

				// xor
				bits = chip.embed[id_1][w] | chip.embed[id_1 + 1][w];
				bits = bits ^ (chip.embed[id_2][w] | chip.embed[id_2 + 1][w]);

				// turn on very last bit of mask only
				mask = 0x01 << (Integer.SIZE - 1);
			}

			if ((mask & bits) != 0) hd++;

			// shift interrogating bit to the right
			// ('>>>' means a unsigned shift)
			mask >>>= 1;
		}

		return hd;
	}

	protected static double weightedDistance (SimpleChip chip, int id_1, int id_2)
	{
		// to do

		return 0;
	}

	protected static double weightedDistance (AffymetrixChip chip, int id_1, int id_2)
	{
		// to do

		return 0;
	}

	protected static double spotConflict (Chip chip, int row, int col)
	{
		if (chip.spot[row][col] == Chip.EMPTY_SPOT)
			return 0;
		else
			return spotConflict (chip, row, col, chip.spot[row][col]);
	}

	protected static double spotConflict (Chip chip, int row, int col, int probe_id)
	{
		double	conflict, pos_mult;
		int		base, r, c, w, mask, embed_len, probe_len, num_rows, num_cols;

		num_rows = chip.getNumberOfRows();
		num_cols = chip.getNumberOfColumns();
		embed_len = chip.getEmbeddingLength();
		probe_len = chip.getProbeLength();

		conflict = 0;

		base = 0;
		w = -1;
		mask = 0;

		for (int step = 0; step < embed_len; step++)
		{
			// prepare mask to interrogate the bit
			// corresponding to the current masking step
			if (step % Integer.SIZE == 0)
			{
				// jump to the next 4-byte word
				w++;

				// prepare mask to interrogate first bit
				mask = 0x01 << (Integer.SIZE - 1);
			}
			else
			{
				// shift mask to the right
				mask >>>= 1;
			}

			// check the state of the proposed
			// embedding at the current masking step
			if ((chip.embed[probe_id][w] & mask) != 0)
			{
				// spot is in an unmasked step
				// (a nucleotide is being synthesized)

				// increment the number of synthesized bases
				base++;

				// no conflict in this step (light directed to
				// neighboring spots cannot cause any damage here)
				continue;
			}

			// compute position multiplier: conflicts
			// are more harmful in the middle of a probe
			// (a conflict would harm the next nucleotide
			// to be synthesized, if there is any)
			// pos_mult = sqrt ( min (base + 1, probe_len - base + 1) )
			pos_mult = (base + 1 <= probe_len - base + 1) ?
						base + 1 :  probe_len - base + 1;

			pos_mult = Math.sqrt(pos_mult);

			r = (row >= 3) ? row - 3 : 0;
			for (; r <= row + 3 && r < num_rows; r++)
			{
				c = (col >= 3) ? col - 3 : 0;
				for (; c <= col + 3 && c < num_cols; c++)
				{
					// skip if neighbor is empty
					if (chip.spot[r][c] == chip.EMPTY_SPOT) continue;

					// conflict only when neighbor is unmasked
					if ((chip.embed[chip.spot[r][c]][w] & mask) == 0) continue;

					conflict += pos_mult * WEIGHT_DIST[r - row + 3][c - col + 3];
				}
			}
		}

		return conflict;
	}
}
