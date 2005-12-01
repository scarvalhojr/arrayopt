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
	private static final int POS_INFINITY = 1000; // Integer.MAX_VALUE;   <= FIX THIS!!!!

	// distance-based conflict weight function
	public static final double WEIGHT_DIST[][] = {	{	0,		0,		0.1,	0.1111,	0.1,	0,		0		},
													{	0,		0.125,	0.2,	0.25,	0.2,	0.125,	0		},
													{	0.1,	0.2,	0.5,	1,		0.5,	0.2,	0.1		},
													{	0.1111,	0.25,	1,		0,		1,		0.25,	0.1111	},
													{	0.1,	0.2,	0.5,	1,		0.5,	0.2,	0.1		},
													{	0,		0.125,	0.2,	0.25,	0.2,	0.125,	0		},
													{	0,		0,		0.1,	0.1111,	0.1,	0,		0		}};

	private static int matrix[][];

	private static int emb_1[];

	private static int emb_2[];

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
	// compute the minimum hamming distance between all possible embeddings
	// of probe 1 and a fixed embedding of probe 2
	static int minHammingDistance (Chip chip, int id_1, int id_2)
	{
		if (chip instanceof SimpleChip)
			return minHammingDistance ((SimpleChip) chip, id_1, id_2);

		else if (chip instanceof AffymetrixChip)
			return minHammingDistance ((AffymetrixChip) chip, id_1, id_2);

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

	// package access
	// reembed probe 1 so that it has the minimum hamming distance to
	// a fixed embedding of probe 2
	static void optimalReembedding (Chip chip, int id_1, int id_2)
	{
		if (chip instanceof SimpleChip)
			optimalReembedding ((SimpleChip) chip, id_1, id_2);

		else if (chip instanceof AffymetrixChip)
			optimalReembedding ((AffymetrixChip) chip, id_1, id_2);

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

	// compute the minimum hamming distance between all possible embeddings
	// of probe 1 and a fixed embedding of probe 2
	protected static int minHammingDistance (SimpleChip chip, int id_1, int id_2)
	{
		int 	n, m;

		// n is the number of rows (length
		// of probe 1, i.e. number of bases),
		// m is the number of columns
		n = 1 + chip.getProbeLength();
		m = 1 + chip.getEmbeddingLength();

		// compute matrix values
		computeMinDistanceMatrix (n, m, chip.embed[id_1], chip.embed[id_2], chip.dep_seq);

		// return the minimum cost
		return matrix[n - 1][m - 1];
	}

	protected static int minHammingDistance (AffymetrixChip chip, int id_1, int id_2)
	{
		int 	n, m, words;

		// n is the number of rows (length
		// of combined probe pair of id 1),
		// m is the number of columns
		n = 1 + chip.getProbeLength() + 1;
		m = 1 + chip.getEmbeddingLength();

		// assume both have the same number of words
		words = chip.embed[id_1].length;

		// create auxiliary arrays if necessary
		if (emb_1 == null || emb_2 == null)
		{
			emb_1 = new int[words];
			emb_2 = new int[words];
		}
		else if (emb_1.length < words || emb_2.length < words)
		{
			emb_1 = new int[words];
			emb_2 = new int[words];
		}

		// combine embeddings
		for (int i = 0; i < words; i++)
		{
			emb_1[i] = chip.embed[id_1][i] | chip.embed[id_1 + 1][i];
			emb_2[i] = chip.embed[id_2][i] | chip.embed[id_2 + 1][i];
		}

		// BEWARE: this is not optimal!! in this way, we are only considering
		// one possible combined sequence. to do it right, we first need to
		// to change how the computeMinDistanceMatrix method works. it must
		// not take an embedding, but a probe as a sequence of bases. then,
		// we would need to compute two possibilites, on for each possible
		// combined sequence, only to choose the best one

		// compute matrix values
		computeMinDistanceMatrix (n, m, emb_1, emb_2, chip.dep_seq);

		// return the minimum cost
		return matrix[n - 1][m - 1];
	}

	protected static void optimalReembedding (SimpleChip chip, int id_1, int id_2)
	{
		int 	n, m, i, j, cost_unpr;
		int		pos_2, word_2, mask_2 = 0;

		// n is the number of rows (length
		// of probe 1, i.e. number of bases),
		// m is the number of columns
		n = 1 + chip.getProbeLength();
		m = 1 + chip.getEmbeddingLength();

		// compute matrix values
		computeMinDistanceMatrix (n, m, chip.embed[id_1], chip.embed[id_2], chip.dep_seq);

		// since all we know about the probe is
		// a given embedding, we first need to
		// find out which position of the given
		// embedding contains the last base
		pos_2  = chip.getEmbeddingLength() - 1;
		word_2 = pos_2 / Integer.SIZE;
		mask_2 = 0x01 << ((word_2 + 1) * Integer.SIZE - pos_2 - 1);

		for (i = n - 1, j = m - 1; j > 0; j--)
		{
			if ((chip.embed[id_2][word_2] & mask_2) == 0)
				cost_unpr = 0;
			else
				cost_unpr = 1;

			if (matrix[i][j] != matrix[i][j - 1] + cost_unpr)
			{
				// set bit to 1
				chip.embed[id_1][word_2] |= mask_2;

				i--;
			}
			else
			{
				// set bit to 0
				chip.embed[id_1][word_2] &= (~ mask_2);
			}

			// prepare mask to interrogate id_2's previous position
			if ((--pos_2 + 1) % Integer.SIZE == 0)
			{
				// previous 4-byte word
				word_2--;

				// turn on very first bit of mask
				mask_2 = 0x01;
			}
			else
			{
				// shift interrogating bit to the left
				mask_2 <<= 1;
			}
		}
	}

	protected static void optimalReembedding (AffymetrixChip chip, int id_1, int id_2)
	{
		int 	n, m, i, j, cost_unpr, middle_pm, middle_mm;
		int		pos_2, words, word_2, mask_2 = 0;

		// n is the number of rows (length
		// of combined probe pair of id 1),
		// m is the number of columns
		n = 1 + chip.getProbeLength() + 1;
		m = 1 + chip.getEmbeddingLength();

		// THIS IS WRONG!! The order in which the middle bases
		// appear in the combined probe might be different!!
		middle_pm = n / 2;
		middle_mm = middle_pm + 1;

		// assume both have the same number of words
		words = chip.embed[id_1].length;

		// create auxiliary arrays if necessary
		if (emb_1 == null || emb_2 == null)
		{
			emb_1 = new int[words];
			emb_2 = new int[words];
		}
		else if (emb_1.length < words || emb_2.length < words)
		{
			emb_1 = new int[words];
			emb_2 = new int[words];
		}

		// combine embeddings
		for (i = 0; i < words; i++)
		{
			emb_1[i] = chip.embed[id_1][i] | chip.embed[id_1 + 1][i];
			emb_2[i] = chip.embed[id_2][i] | chip.embed[id_2 + 1][i];
		}

		// BEWARE: this is not optimal!! in this way, we are only considering
		// one possible combined sequence. to do it right, we first need to
		// to change how the computeMinDistanceMatrix method works. it must
		// not take an embedding, but a probe as a sequence of bases. then,
		// we would need to compute two possibilites, on for each possible
		// combined sequence, only to choose the best one

		// compute matrix values
		computeMinDistanceMatrix (n, m, emb_1, emb_2, chip.dep_seq);

		// since all we know about the probe is
		// a given embedding, we first need to
		// find out which position of the given
		// embedding contains the last base
		pos_2  = chip.getEmbeddingLength() - 1;
		word_2 = pos_2 / Integer.SIZE;
		mask_2 = 0x01 << ((word_2 + 1) * Integer.SIZE - pos_2 - 1);

		for (i = n - 1, j = m - 1; j > 0; j--)
		{
			if ((emb_2[word_2] & mask_2) == 0)
				cost_unpr = 0;
			else
				cost_unpr = 1;

			if (matrix[i][j] == matrix[i][j - 1] + cost_unpr)
			{
				// set bits of both PM and MM probes to 0
				chip.embed[id_1][word_2] &= (~ mask_2);
				chip.embed[id_1 + 1][word_2] &= (~ mask_2);
			}
			else
			{
				if (i == middle_pm)
				{
					// set bit of PM probe to 1
					chip.embed[id_1][word_2] |= mask_2;

					// and set bit of MM probe to 0
					chip.embed[id_1 + 1][word_2] &= (~ mask_2);
				}
				else if (i == middle_mm)
				{
					// set bit of PM probe to 0
					chip.embed[id_1][word_2] &= (~ mask_2);

					// and set bit of MM probe to 1
					chip.embed[id_1 + 1][word_2] |= mask_2;
				}
				else
				{
					// set bits of both PM and MM probes to 1
					chip.embed[id_1][word_2] |= mask_2;
					chip.embed[id_1 + 1][word_2] |= mask_2;
				}

				i--;
			}

			// prepare mask to interrogate id_2's previous position
			if ((--pos_2 + 1) % Integer.SIZE == 0)
			{
				// previous 4-byte word
				word_2--;

				// turn on very first bit of mask
				mask_2 = 0x01;
			}
			else
			{
				// shift interrogating bit to the left
				mask_2 <<= 1;
			}
		}
	}

	protected static void computeMinDistanceMatrix (int n, int m, int e1[], int e2[], char dep_seq[])
	{
		int 	i, j, cost_prod, cost_unpr, min_j, max_j;
		int		pos_1, pos_2, word_1, word_2, mask_1 = 0, mask_2 = 0;
		boolean	started;

		// check if matrix has already been instantiated
		if (matrix == null)
		{
			// no, so allocate one
			matrix = new int[n][m];
		}
		else
		{
			// yes: check dimensions
			if (matrix.length < n || matrix[0].length < m)
			{
				// insufficient dimensions: create a new one
				matrix = new int[n][m];
			}
		}

		matrix[0][0] = 0;

		min_j = 1;
		max_j = m - n;

		// initialize first row
		for (j = min_j, pos_2 = -1, word_2 = -1; j <= max_j; j++)
		{
			// prepare mask to interrogate e2's next position
			if (++pos_2 % Integer.SIZE == 0)
			{
				// next 4-byte word
				word_2++;

				// turn on very last bit of mask only
				mask_2 = 0x01 << (Integer.SIZE - 1);
			}
			else
			{
				// shift interrogating bit to the right
				// ('>>>' means a unsigned shift)
				mask_2 >>>= 1;
			}

			// the cost depends on whether e2's embedding
			// is productive or not at this step
			if ((e2[word_2] & mask_2) != 0)
			{
				// e2 is in a productive step:
				// there will be a conflict if e1
				// is put in an unproductive step
				matrix[0][j] = matrix[0][j - 1] + 1;
			}
			else
			{
				// e2 is in an unproductive step:
				// there will be no conflict if e1
				// is put in a productive step
				matrix[0][j] = matrix[0][j - 1];
			}
		}

		// compute matrix (row-wise)
		for (i = 1, pos_1 = -1, word_1 = -1; i < n; i++)
		{
			// since all we know about the probe is
			// a given embedding, we first need to
			// find out which position of the embedding
			// contains the i-th base
			do
			{
				// prepare mask to interrogate e1's next position
				if (++pos_1 % Integer.SIZE == 0)
				{
					// next 4-byte word
					word_1++;

					// turn on very last bit of mask only
					mask_1 = 0x01 << (Integer.SIZE - 1);
				}
				else
				{
					// shift interrogating bit to the right
					// ('>>>' means a unsigned shift)
					mask_1 >>>= 1;
				}

				// continue until the next productive step is found
			}
			while ((e1[word_1] & mask_1) == 0);

			if (max_j < m - 1) max_j ++;

			pos_2 = min_j - 1;
			word_2 = pos_2 / Integer.SIZE;
			mask_2 = 0x01 << ((word_2 + 1) * Integer.SIZE - pos_2 - 1);

			matrix[i][min_j - 1] = POS_INFINITY;

			started = false;

			// compute the values on the i-th row of the matrix
			for (j = min_j; j <= max_j; j++)
			{
				// two options need to be evaluated:
				// => e1's is productive (i-th base is synthesized at j-th step)
				// => e1's is unproductive (i-th base is *not* synthesized now)
				//
				// the cost of each option depends on:
				// a) whether the deposition sequence allows e1's
				//    current base to be synthesized at this step
				// b) whether e2's embedding is productive or not
				//    at this step

				// condition a)
				if (dep_seq[pos_1] != dep_seq[pos_2])
				{
					// i-th base of e1 cannot be
					// synthesized at j-th step
					cost_prod = POS_INFINITY;
				}
				else
				{
					cost_prod = matrix[i - 1][j - 1];
				}

				// condition b)

				cost_unpr = matrix[i][j - 1];

				// check the state of e2's embedding
				// at the current position
				if ((e2[word_2] & mask_2) == 0)
				{
					// e2 is in an unproductive step:
					// there will be a conflict if e1
					// is put in a productive step
					if (cost_prod != POS_INFINITY)
						cost_prod++;
				}
				else
				{
					// e2 is in a productive step:
					// there will be a conflict if e1
					// is put in an unproductive step
					if (cost_unpr != POS_INFINITY)
						cost_unpr++;
				}

				// choose the option with minimum cost
				matrix[i][j] = (cost_prod < cost_unpr ? cost_prod : cost_unpr);

				if (!started)
					if (matrix[i][j] != POS_INFINITY)
					{
						started = true;
						min_j = j + 1;
					}

				// prepare mask to interrogate e2's next position
				if (++pos_2 % Integer.SIZE == 0)
				{
					// next 4-byte word
					word_2++;

					// turn on very last bit of mask only
					mask_2 = 0x01 << (Integer.SIZE - 1);
				}
				else
				{
					// shift interrogating bit to the right
					// ('>>>' means a unsigned shift)
					mask_2 >>>= 1;
				}
			}
		}
	}

	protected static void computeMinDistanceMatrix_old (int n, int m, int e1[], int e2[], char dep_seq[])
	{
		int 	i, j, cost_prod, cost_unpr;
		int		pos_1, pos_2, word_1, word_2, mask_1 = 0, mask_2 = 0;
		char	base;

		// check if matrix has already been instantiated
		if (matrix == null)
		{
			// no, so allocate one
			matrix = new int[n][m];
		}
		else
		{
			// yes: check dimensions
			if (matrix.length < n || matrix[0].length < m)
			{
				// insufficient dimensions: create a new one
				matrix = new int[n][m];
			}
		}

		// initialize first column
		matrix[0][0] = 0;
		for (i = 1; i < n; i++)
			matrix[i][0] = POS_INFINITY;

		// initialize first row
		for (j = 1, pos_2 = -1, word_2 = -1; j < m; j++)
		{
			// prepare mask to interrogate e2's next position
			if (++pos_2 % Integer.SIZE == 0)
			{
				// next 4-byte word
				word_2++;

				// turn on very last bit of mask only
				mask_2 = 0x01 << (Integer.SIZE - 1);
			}
			else
			{
				// shift interrogating bit to the right
				// ('>>>' means a unsigned shift)
				mask_2 >>>= 1;
			}

			matrix[0][j] = matrix[0][j - 1];

			// the cost depends on whether e2's embedding
			// is productive or not at this step
			if ((e2[word_2] & mask_2) != 0)
			{
				// e2 is in a productive step:
				// there will be a conflict if e1
				// is put in an unproductive step
				matrix[0][j] ++;
			}
		}

		// compute matrix (row-wise)
		for (i = 1, pos_1 = -1, word_1 = -1; i < n; i++)
		{
			// since all we know about the probe is
			// a given embedding, we first need to
			// find out which position of the embedding
			// contains the i-th base
			do
			{
				// prepare mask to interrogate e1's next position
				if (++pos_1 % Integer.SIZE == 0)
				{
					// next 4-byte word
					word_1++;

					// turn on very last bit of mask only
					mask_1 = 0x01 << (Integer.SIZE - 1);
				}
				else
				{
					// shift interrogating bit to the right
					// ('>>>' means a unsigned shift)
					mask_1 >>>= 1;
				}

				// continue until the next productive step is found
			}
			while ((e1[word_1] & mask_1) == 0);

			// find out which base is at the current position
			base = dep_seq[pos_1];

			// compute the values on the i-th row of the matrix
			for (pos_2 = -1, word_2 = -1, j = 1; j < m; j++)
			{
				// prepare mask to interrogate e2's next position
				if (++pos_2 % Integer.SIZE == 0)
				{
					// next 4-byte word
					word_2++;

					// turn on very last bit of mask only
					mask_2 = 0x01 << (Integer.SIZE - 1);
				}
				else
				{
					// shift interrogating bit to the right
					// ('>>>' means a unsigned shift)
					mask_2 >>>= 1;
				}

				// two options need to be evaluated:
				// => e1's is productive (i-th base is synthesized at j-th step)
				// => e1's is unproductive (i-th base is *not* synthesized now)
				cost_prod = matrix[i - 1][j - 1];
				cost_unpr = matrix[i][j - 1];

				// the cost of each option depends on:
				// a) whether the deposition sequence allows e1's
				//    current base to be synthesized at this step
				// b) whether e2's embedding is productive or not
				//    at this step

				// condition a)
				if (dep_seq[pos_2] != base)
				{
					// i-th base of e1 cannot be
					// synthesized at j-th step
					cost_prod = POS_INFINITY;
				}

				// condition b)

				// check the state of e2's embedding
				// at the current position
				if ((e2[word_2] & mask_2) == 0)
				{
					// e2 is in an unproductive step:
					// there will be a conflict if e1
					// is put in a productive step
					if (cost_prod != POS_INFINITY)
						cost_prod++;
				}
				else
				{
					// e2 is in a productive step:
					// there will be a conflict if e1
					// is put in an unproductive step
					if (cost_unpr != POS_INFINITY)
						cost_unpr++;
				}

				// choose the option with minimum cost
				matrix[i][j] = (cost_prod < cost_unpr ? cost_prod : cost_unpr);
			}
		}
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
