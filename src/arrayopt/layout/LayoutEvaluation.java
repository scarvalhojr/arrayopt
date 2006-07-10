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
	/**
	 * Computes the Hamming distance between the embeddings of two probes. The
	 * Hamming distance gives the number of masking steps that the two
	 * embeddings differ.
	 * 
	 * @param chip the chip containing the probes
	 * @param id_1 ID of first probe
	 * @param id_2 ID of second probe 
	 * @return the Hamming distance between their embeddings 
	 */
	public static int hammingDistance (Chip chip, int id_1, int id_2)
	{
		if (chip instanceof SimpleChip)
			return hammingDistance ((SimpleChip) chip, id_1, id_2);

		if (chip instanceof AffymetrixChip)
			return hammingDistance ((AffymetrixChip) chip, id_1, id_2);

		throw new IllegalArgumentException ("Unsupported chip type.");
	}
	
	public static int hammingDistance (SimpleChip chip, int id_1, int id_2)
	{
		int w, lastw, hd = 0, shift, bits;
		
		// index of the last word
		lastw = chip.embed[id_1].length - 1;

		// count the differences in the first words
		for (w = 0; w < lastw; w++)
			hd += Integer.bitCount(chip.embed[id_1][w] ^ chip.embed[id_2][w]);
		
		// bitwise xor of the last word
		bits = chip.embed[id_1][lastw] ^ chip.embed[id_2][lastw];
		
		// clear any unused bits
		shift = (1 + lastw) * Integer.SIZE - chip.embed_len;
		bits >>>= shift;
		bits <<= shift;

		// count the differences in the last word
		hd += Integer.bitCount(bits);

		return hd;
	}

	public static int hammingDistanceSpots (Chip chip, int id_1, int id_2)
	{
		int w, lastw, hd = 0, shift, bits;
		
		// index of the last word
		lastw = chip.embed[id_1].length - 1;

		// count the differences in the first words
		for (w = 0; w < lastw; w++)
			hd += Integer.bitCount(chip.embed[id_1][w] ^ chip.embed[id_2][w]);
		
		// bitwise xor of the last word
		bits = chip.embed[id_1][lastw] ^ chip.embed[id_2][lastw];
		
		// clear any unused bits
		shift = (1 + lastw) * Integer.SIZE - chip.embed_len;
		bits >>>= shift;
		bits <<= shift;

		// count the differences in the last word
		hd += Integer.bitCount(bits);

		return hd;
	}

	public static int hammingDistance (AffymetrixChip chip, int id_1, int id_2)
	{
		int w, lastw, hd = 0, shift, bits;
		
		// index of the last word
		lastw = chip.embed[id_1].length - 1;

		// count the differences in the first words
		for (w = 0; w < lastw; w++)
		{
			bits = chip.embed[id_1][w] | chip.embed[id_1 + 1][w];
			bits = bits ^ (chip.embed[id_2][w] | chip.embed[id_2 + 1][w]);
			hd += Integer.bitCount(bits);
		}
		
		// bitwise xor of the last word
		bits = chip.embed[id_1][lastw] | chip.embed[id_1 + 1][lastw];
		bits = bits ^ (chip.embed[id_2][lastw] | chip.embed[id_2 + 1][lastw]);
		
		// clear any unused bits
		shift = (1 + lastw) * Integer.SIZE - chip.embed_len;
		bits >>>= shift;
		bits <<= shift;

		// count the differences in the last word
		hd += Integer.bitCount(bits);
		
		return hd;
	}
	
	/**
	 * Computes the distance between the embeddings of two probes using the
	 * distance-dependent weights given by the current definition of conflict
	 * index (see {@link ConflictIndex}).
	 * 
	 * @param chip chip containing the probes
	 * @param id_1 ID of first probe
	 * @param id_2 ID of second probe
	 * @return the weighted distance between their embeddings
	 */
	public static double weightedDistance (Chip chip, int id_1, int id_2)
	{
		if (chip instanceof SimpleChip)
			return weightedDistance ((SimpleChip) chip, id_1, id_2);

		if (chip instanceof AffymetrixChip)
			return weightedDistance ((AffymetrixChip) chip, id_1, id_2);

		throw new IllegalArgumentException ("Unsupported chip type.");
	}
	
	public static double weightedDistance (SimpleChip chip, int id_1,
			int id_2)
	{
		int base = 0, bitmask = 0, w, pos, probe_len;
		double dist = 0;
		
		probe_len = chip.getProbeLength();

		for (w = -1, pos = 0; pos < chip.getEmbeddingLength(); pos++)
		{
			if (pos % Integer.SIZE == 0)
			{
				// next 4-byte word
				w++;

				// turn on very last bit of mask only
				bitmask = 0x01 << (Integer.SIZE - 1);
			}
			
			if ((bitmask & chip.embed[id_1][w]) != 0)
			{
				// probe id_1 is unmasked: no conflict comming from id_2
				base++;
			}
			else if ((bitmask & chip.embed[id_2][w]) != 0)
			{
				// probe id_2 is unmasked: there is a conflict!
				dist += ConflictIndex.positionWeight(base, probe_len);
			}

			// shift interrogating bit to the right
			// ('>>>' means a unsigned shift)
			bitmask >>>= 1;
		}

		return dist;
	}

	public static double weightedDistance (AffymetrixChip chip, int id_1,
			int id_2)
	{
		// TODO implement this

		return hammingDistance(chip, id_1, id_2);
	}
	
	public static double analyzeBorderLength (Chip chip)
	{
		RectangularRegion region;
		int embed_len, num_probes, r, c, id1, id2, step, word, bitmask = 0;
		int border[], total_bl = 0;
		double norm_bl[], total_norm_bl = 0;
		
		region = chip.getChipRegion();
		embed_len = chip.getEmbeddingLength();
		num_probes = chip.getNumberOfProbes();
		
		border = new int [embed_len];
		norm_bl = new double [embed_len];
		
		for (word = -1, step = 0; step < embed_len; step++)
		{
			if (step % Integer.SIZE == 0)
			{
				word++;
				bitmask = 0x01 << (Integer.SIZE - 1);
			}
			else
				bitmask >>>= 1;
			
			// compute current mask's border length
			border[step] = 0;
			
			// first, conflicts in the vertical borders between spots
			for (r = region.first_row; r <= region.last_row; r ++)
				for (c = region.first_col; c < region.last_col; c++)
				{
					if ((id1 = chip.spot[r][c]) == Chip.EMPTY_SPOT)
						continue;
					
					if ((id2 = chip.spot[r][c + 1]) == Chip.EMPTY_SPOT)
						continue;
					
					if (((chip.embed[id1][word] ^ chip.embed[id2][word])
							& bitmask) != 0)
						border[step]++;
				}

			// now, conflicts in the horizontal borders between spots
			for (c = region.first_col; c <= region.last_col; c ++)
				for (r = region.first_row; r < region.last_row; r++)
				{
					if ((id1 = chip.spot[r][c]) == Chip.EMPTY_SPOT)
						continue;
					
					if ((id2 = chip.spot[r + 1][c]) == Chip.EMPTY_SPOT)
						continue;
					
					if (((chip.embed[id1][word] ^ chip.embed[id2][word])
							& bitmask) != 0)
						border[step]++;
				}
			
			norm_bl[step] = border[step] / (double) num_probes;
			total_bl += border[step];
			total_norm_bl += norm_bl[step]; 
			
			System.err.println("Step " + step + " => border length: " +
					border[step] + " (normalized: " + norm_bl[step] +
					")");
		}
		
		System.err.println("Total border length: " + total_bl +
				" (normalized: " + total_norm_bl + ")");
		
		return total_norm_bl;
	}

	public static long borderLength (Chip chip)
	{
		return borderLength (chip, chip.getChipRegion());
	}

	public static double normalizedBorderLength (Chip chip)
	{
		return borderLength (chip, chip.getChipRegion()) /
								chip.getNumberOfProbes();
	}

	public static long borderLength (Chip chip, RectangularRegion region)
	{
		long border = 0;
		int r, c, id1, id2;
		
		for (r = region.first_row; r <= region.last_row; r ++)
			for (c = region.first_col; c < region.last_col; c++)
			{
				if ((id1 = chip.spot[r][c]) == Chip.EMPTY_SPOT)
					continue;
				if ((id2 = chip.spot[r][c + 1]) == Chip.EMPTY_SPOT)
					continue;
				border += hammingDistanceSpots(chip, id1, id2);
			}
			
		for (c = region.first_col; c <= region.last_col; c++)
			for (r = region.first_row; r < region.last_row; r ++)
			{
				if ((id1 = chip.spot[r][c]) == Chip.EMPTY_SPOT)
					continue;
				if ((id2 = chip.spot[r + 1][c]) == Chip.EMPTY_SPOT)
					continue;
				border += hammingDistanceSpots(chip, id1, id2);
			}

		return border;
	}

	public static long borderLength (Chip chip, int row, int col)
	{
		int id;
		
		if ((id = chip.spot[row][col]) == Chip.EMPTY_SPOT)
			return 0;
		
		if (chip instanceof SimpleChip)
			return borderLength ((SimpleChip) chip, row, col, id);

		if (chip instanceof AffymetrixChip)
			return borderLength ((AffymetrixChip) chip, row, col, id);

		throw new IllegalArgumentException ("Unsupported chip type.");
	}

	public static long borderLength (SimpleChip chip, int row, int col,
			int id)
	{
		RectangularRegion region;
		int id2, border = 0;

		region = chip.getChipRegion();
		
		if (row - 1 >= region.first_row)
			if ((id2 = chip.spot[row - 1][col]) != Chip.EMPTY_SPOT)
				border += hammingDistance(chip, id, id2);

		if (row + 1 <= region.last_row)
			if ((id2 = chip.spot[row + 1][col]) != Chip.EMPTY_SPOT)
				border += hammingDistance(chip, id, id2);

		if (col - 1 >= region.first_col)
			if ((id2 = chip.spot[row][col - 1]) != Chip.EMPTY_SPOT)
				border += hammingDistance(chip, id, id2);

		if (col + 1 <= region.last_col)
			if ((id2 = chip.spot[row][col + 1]) != Chip.EMPTY_SPOT)
				border += hammingDistance(chip, id, id2);

		return border;
	}

	public static long borderLength (AffymetrixChip chip, int row, int col,
			int id)
	{
		// TODO implement this
		
		return -1;
	}

	public static double analyzeConflictIndex (Chip chip)
	{
		RectangularRegion region;
		double	conf, avg_conf, min_conf, max_conf;
		int		id, num_probes, min_row, max_row, min_col, max_col;
		
		min_conf = Double.POSITIVE_INFINITY;
		avg_conf = max_conf = 0;
		min_row = max_row = min_col = max_col = -1;
		
		region = chip.getChipRegion();
		num_probes = chip.getNumberOfProbes();
		
		for (int r = region.first_row; r <= region.last_row; r++)
			for (int c = region.first_col; c <= region.last_col; c++)
				if ((id = chip.spot[r][c]) != Chip.EMPTY_SPOT)
				{
					conf = conflictIndex(chip, r, c, id);
					
					if (conf < min_conf)
					{
						min_conf = conf;
						min_row = r;
						min_col = c;
					}
					else if (conf > max_conf)
					{
						max_conf = conf;
						max_row = r;
						max_col = c;
					}
					
					avg_conf += conf / num_probes; 
				}
		
		System.err.println("Min conflict: " + min_conf + "(" + min_row + "," +
				min_col + ")");
		System.err.println("Avg conflict: " + avg_conf);
		System.err.println("Max conflict: " + max_conf + "(" + max_row + "," +
				max_col + ")");
		
		return avg_conf;
	}

	public static double averageConflictIndex (Chip chip)
	{
		RectangularRegion region;
		double	conf = 0;
		int		id, num_probes;
		
		region = chip.getChipRegion();
		num_probes = chip.getNumberOfProbes();
		
		for (int r = region.first_row; r <= region.last_row; r++)
			for (int c = region.first_col; c <= region.last_col; c++)
				if ((id = chip.spot[r][c]) != Chip.EMPTY_SPOT)
					conf += conflictIndex(chip, r, c, id) / num_probes;
		
		return conf;
	}

	public static double totalConflictIndex (Chip chip, RectangularRegion
			region)
	{
		double	conf = 0;
		int		id;
		
		for (int r = region.first_row; r <= region.last_row; r++)
			for (int c = region.first_col; c <= region.last_col; c++)
				if ((id = chip.spot[r][c]) != Chip.EMPTY_SPOT)
					conf += conflictIndex(chip, r, c, id);
		
		return conf;
	}

	/**
	 * Computes the conflict index of a spot. If the spot is empty, its confict
	 * index is zero. Otherwise, this method calls the
	 * {@link #conflictIndex(Chip, int, int, int)} method.
	 * 
	 * @param chip chip containing the spot
	 * @param row row coordinate of the spot
	 * @param col column coordinate of the spot
	 * @return conflict index of the spot
	 */
	public static double conflictIndex (Chip chip, int row, int col)
	{
		int id;
		
		if ((id = chip.spot[row][col]) == Chip.EMPTY_SPOT)
			return 0;
		
		if (chip instanceof SimpleChip)
			return conflictIndex ((SimpleChip) chip, row, col, id);

		if (chip instanceof AffymetrixChip)
			return conflictIndex ((AffymetrixChip) chip, row, col, id);

		throw new IllegalArgumentException ("Unsupported chip type.");
	}

	/**
	 * 
	 */
	public static double conflictIndex (Chip chip, int row, int col, int id)
	{
		if (chip instanceof SimpleChip)
			return conflictIndex ((SimpleChip) chip, row, col, id);

		if (chip instanceof AffymetrixChip)
			return conflictIndex ((AffymetrixChip) chip, row, col, id);

		throw new IllegalArgumentException ("Unsupported chip type.");
	}

	/**
	 * Computes the conflict index of probe when it is placed on a given spot.
	 * This method analyzes the probes in the region around the spot and uses
	 * the current definition of conflict index (see {@link ConflictIndex}).
	 * 
	 * @param chip chip containing the spot
	 * @param row row coordinate of the spot
	 * @param col column coordinate of the spot
	 * @param pid probe ID 
	 * @return conflict index of the spot
	 */
	public static double conflictIndex (SimpleChip chip, int row, int col,
			int pid)
	{
		RectangularRegion	region;
		double	conf, posw;
		int		ci_dim, embed_len, probe_len, r, c, step;
		int		word, base, bitmask = 0;

		ci_dim = ConflictIndex.dimConflictRegion();
		region = chip.getChipRegion();
		embed_len = chip.getEmbeddingLength();
		probe_len = chip.getProbeLength();

		// define region around the spot that needs to be examined
		int min_row = Math.max(row - ci_dim, region.first_row);
		int max_row = Math.min(row + ci_dim, region.last_row);
		int min_col = Math.max(col - ci_dim, region.first_col);
		int max_col = Math.min(col + ci_dim, region.last_col);
		
		conf = 0;

		for (base = 0, step = 0, word = - 1; step < embed_len; step++)
		{
			if (step % Integer.SIZE == 0)
			{
				word++;
				bitmask = 0x01 << (Integer.SIZE - 1);
			}
			else
				bitmask >>>= 1;

			// check the state of the embedding
			// at the current masking step
			if ((chip.embed[pid][word] & bitmask) != 0)
			{
				// spot is in an unmasked step
				// (a nucleotide is being synthesized)

				// increment the number of synthesized bases
				base++;

				// no conflict in this step (light directed to
				// neighboring spots cannot cause any damage here)
				continue;
			}

			// masked step: compute position multiplier (a conflict
			// would harm the next nucleotide to be synthesized)
			posw = ConflictIndex.positionWeight (base, probe_len);

			for (r = min_row; r <= max_row; r++)
			{
				for (c = min_col; c <= max_col; c++)
				{
					// skip if neighbor is empty
					if (chip.spot[r][c] == Chip.EMPTY_SPOT)
						continue;

					// conflict only when neighbor is unmasked
					if ((chip.embed[chip.spot[r][c]][word] & bitmask) == 0)
						continue;

					conf += posw * ConflictIndex.distanceWeight(r,c, row, col);
				}
			}
		}

		return conf;
	}

	/**
	 * Computes the conflict index of probe when it is placed on a given spot.
	 * This method analyzes the probes in the region around the spot and uses
	 * the current definition of conflict index (see {@link ConflictIndex}).
	 * 
	 * @param chip chip containing the spot
	 * @param row row coordinate of the spot
	 * @param col column coordinate of the spot
	 * @param pid probe ID 
	 * @return conflict index of the spot
	 */
	public static double conflictIndex (AffymetrixChip chip, int row, int col,
			int pid)
	{
		// TODO implement this
		
		return -1;
	}
}
