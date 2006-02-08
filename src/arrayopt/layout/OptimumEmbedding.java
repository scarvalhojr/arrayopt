/*
 * OptimumEmbedding.java
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
 * This class implements several methods for computing an optimum embedding of a
 * probe. It can work with two types of distance measures:
 * {@link #MODE_BORDER_LENGTH} and {@link #MODE_CONFLICT_INDEX}.
 * 
 * <P>Two types of methods are provided by this class: 1) compute the minimum
 * distance that any embedding of a probe (called source) can have to the
 * current embedding of another probe (called destination); 2) re-embed the
 * source probe so that its distance to the destination probe is minimal. The
 * first type contains the <CODE>minDistanceProbe</CODE> and
 * <CODE>minDistanceSpot</CODE> methods while the second class contains the
 * <CODE>reembedProbe</CODE> and <CODE>reembedSpot</CODE> members.></P>
 * 
 * <P>Note that the current embedding of the destination probe is considered
 * when computing the distance (and this embedding is considered as fixed). This
 * fact allows the use of a simple and efficient dynamic programming approach
 * with quadratic space and time complexities. The actual implementation aims at
 * maximizing the reuse of common routines by the sub-classes rather than
 * achiving maximum speed.</P>
 *  
 * <P>In fact, each of two types of methods have several method calls that
 * allows the computed distance to take into account not only a single
 * destination probe but also a set of probes or a neighboring region. Moreover,
 * the source probe can be specified by a probe ID or its current location on
 * the chip.</P>
 * 
 * <P>This class is abstract and the actual implementation details are provided
 * by (inner) sub-classes for each mode (as well as for different types of
 * chips). Note, however, that neither this class or any of its subclasses can
 * be instantiated directly. Instead, the {@link #createEmbedder} factory method
 * should be called with a chip instance and the desired mode of operation
 * passed as arguments. The returned object will provide the functionality
 * specified here as needed.</P>

 * <P>This class is intended for use by other classes inside the layout package
 * and, therefore, has no public methods.</P>
 *
 * <P><B>Note that this implementation is not synchronized.</B> If multiple
 * threads access an OptimumEmbedding instance concurrently, it must be
 * synchronized externally.</P>
 * 
 * @see OptimumEmbedding.Simple.BorderLength
 * @see OptimumEmbedding.Simple.ConflictIndex
 * @see OptimumEmbedding.Affy.BorderLength
 * @see OptimumEmbedding.Affy.ConflictIndex
 */
public abstract class OptimumEmbedding
{
	/**
	 * Reference to the chip instance for which this object is configured to
	 * work with.
	 */
	protected Chip chip;
	
	/**
	 * Reference to the total region of the {@link #chip}. 
	 */
	protected RectangularRegion chip_region;

	/**
	 * An internal copy of the {@link #chip}'s deposition sequence.
	 */
	protected char dep_seq[];
	
	/**
	 * An array holding the costs of masking a spot (the total conflict that
	 * such an operation would create in that spot).
	 */
	protected double mask_cost[];
	
	/**
	 * An array holding the costs of leaving a spot unmasked (the total conflict
	 * that such an operation would create in the neighboring spots).
	 */
	protected double unmask_cost[];
	
	/**
	 * An array with the position multipliers (mainly used when configured with
	 * {@link #MODE_CONFLICT_INDEX}; with {@link #MODE_BORDER_LENGTH}, all
	 * positions of the array contain a 1.
	 */
	protected double pos_mult[];

	/**
	 * The length of the embeddings found on the {@link #chip} (i.e the length
	 * of the deposition sequente). 
	 */
	protected int embed_len;
	
	/**
	 * The length of probe sequence found on the {@link #chip}. 
	 */
	protected int probe_len;

	/**
	 * Constant that indicates that the border length should be considered when
	 * computing the distance between embeddings/spots.
	 */
	static final int MODE_BORDER_LENGTH = 0;
	
	/**
	 * Constant that indicates that the conflict index should be used when
	 * computing the distance between embeddings/spots.
	 */
	static final int MODE_CONFLICT_INDEX = 1;

	/**
	 * Create an instance of the appropriate sub-class according to the type
	 * of chip and the desired distance mode ({@link #MODE_BORDER_LENGTH} or
	 * {@link #MODE_BORDER_LENGTH}).
	 * 
	 * @param c a chip instance  
	 * @param mode type of the desired distance measure
	 * @return an instance of an OptimumEmbedding for the type of chip and the
	 * selected mode 
	 */
	static OptimumEmbedding createEmbedder (Chip c, int mode)
	{
		if (c instanceof SimpleChip)
		{
			if (mode == MODE_BORDER_LENGTH)
				return new Simple.BorderLength ((SimpleChip) c);
			
			if (mode == MODE_CONFLICT_INDEX)
				return new Simple.ConflictIndex ((SimpleChip) c);
			
			throw new IllegalArgumentException
				("Illegal value for argument 'mode'.");
		}
		
		if (c instanceof AffymetrixChip)
		{
			if (mode == MODE_BORDER_LENGTH)
				return new Affy.BorderLength ((AffymetrixChip) c);
			
			if (mode == MODE_CONFLICT_INDEX)
				return new Affy.ConflictIndex ((AffymetrixChip) c);
			
			throw new IllegalArgumentException
				("Illegal value for argument 'mode'.");
		}
		
		// else
		throw new IllegalArgumentException ("Unsupported chip type.");
	}
	
	/**
	 * This constructor can only be accessed internally. For instantiating
	 * objects of this class, please use {@link createEmbedder}.  
	 * @param chip a chip instance
	 * @param probe_len the length of the probes 
	 */
	protected OptimumEmbedding (Chip chip, int probe_len)
	{
		this.chip = chip;
		this.probe_len = probe_len;
		this.embed_len = chip.getEmbeddingLength();
		
		// create a local copy of the deposition sequence
		// shifted once to the right (aligned with the matrix)
		this.dep_seq = new char[embed_len + 1];
		this.dep_seq[0] = ' ';
		System.arraycopy(chip.dep_seq, 0, this.dep_seq, 1, embed_len);
		
		// save a local copy of the chip dimensions
		this.chip_region = chip.getChipRegion();
				
		// create auxiliary arrays and matrices
		this.mask_cost = new double [embed_len];
		this.unmask_cost = new double [embed_len];
		this.pos_mult = new double [probe_len + 1];
	}
	
	/**
	 * Computes the minimum distance between any valid embedding of a probe
	 * (<CODE>id_1</CODE>) and the current embedding of <CODE>id_2</CODE>.
	 * @param id_1 the ID of the first probe
	 * @param id_2 the ID of the probe with a fixed embedding
	 * @return the minimum distance between id_1 and id_2
	 */
	double minDistanceProbe (int id_1, int id_2)
	{
		resetCosts();
		addProbeCost (id_2);
		return computeMinDistance(id_1);
	}

	/**
	 * Computes the minimum distance between any valid embedding of a probe
	 * (<CODE>id_1</CODE>) and the current embeddings of a given set of probes.
	 * @param id_1 the ID of the first probe
	 * @param id the IDs of the set of probe with a fixed embedding
	 * @return the minimum distance between id_1 and the set of probes
	 */
	double minDistanceProbe (int id_1, int id[])
	{
		return minDistanceProbe (id_1, id, 0, id.length - 1);
	}

	/**
	 * Computes the minimum distance between any valid embedding of a probe
	 * (<CODE>id_1</CODE>) and the current embeddings of a given set of probes.
	 * @param id_1 the ID of the first probe
	 * @param id the IDs of the set of probe with a fixed embedding
	 * @param start the index of the first element on the list
	 * @param end the index of the last element on the list 
	 * @return the minimum distance between id_1 and the set of probes
	 */
	double minDistanceProbe (int id_1, int id[], int start, int end)
	{
		resetCosts();
		addProbeCost (id, start, end);
		return computeMinDistance(id_1);
	}

	/**
	 * Computes the minimum distance between any valid embedding of the probe
	 * found on a spot and the current embeddings of the neighboring probes.
	 * @param row spot's row coordinate
	 * @param col spot's column coordinate
	 * @return the minimum conflict that the probe on the given spot can cause
	 * to the neighboring probes
	 */
	double minDistanceSpot (int row, int col)
	{
		int id;
		
		if ((id = chip.spot[row][col]) == Chip.EMPTY_SPOT) return 0;
		
		resetCosts();
		addSpotCost (row, col);
		return computeMinDistance(id);
	}

	/**
	 * Reembeds a probe (<CODE>id_1</CODE>) so that its distance to the current
	 * embeddings of probe <CODE>id_2</CODE> is minimal.
	 * @param id_1 the ID of the first probe
	 * @param id_2 the ID of the probe with a fixed embedding
	 * @return the minimum distance between id_1 and id_2
	 */
	double reembedProbe (int id_1, int id_2)
	{
		resetCosts();
		addProbeCost (id_2);
		return reembedOptimally (id_1);
	}

	/**
	 * Reembeds a probe (<CODE>id_1</CODE>) so that its distance to the current
	 * embeddings of the given set of probes is minimal.
	 * @param id_1 the ID of the first probe
	 * @param id the IDs of the set of probe with a fixed embedding
	 * @return the minimum distance between id_1 and the set of probes
	 */
	double reembedProbe (int id_1, int id[])
	{
		return reembedProbe (id_1, id, 0, id.length - 1); 
	}

	/**
	 * Reembeds a probe (<CODE>id_1</CODE>) so that its distance to the current
	 * embeddings of the given set of probes is minimal.
	 * @param id_1 the ID of the first probe
	 * @param id the IDs of the set of probe with a fixed embedding
	 * @param start the index of the first element on the list
	 * @param end the index of the last element on the list 
	 * @return the minimum distance between id_1 and the set of probes
	 */
	double reembedProbe (int id_1, int id[], int start, int end)
	{
		resetCosts();
		addProbeCost (id, start, end);
		return reembedOptimally (id_1);
	}

	/**
	 * Reembeds a probe of the given spot so that its distance to its neighbors
	 * is minimal.
	 * @param row spot's row coordinate
	 * @param col spot's column coordinate
	 * @return the minimum conflict that the probe on the given spot can cause
	 * to the neighboring probes or zero if the spot is empty
	 */
	double reembedSpot (int row, int col)
	{
		int id;
		
		if ((id = chip.spot[row][col]) == Chip.EMPTY_SPOT) return 0;

		resetCosts();
		addSpotCost (row, col);
		return reembedOptimally (id);
	}
	
	/**
	 * Reset cost arrays before a new distance is computed.
	 */
	protected void resetCosts ()
	{
		for (int i = 0; i < embed_len; i++)
			mask_cost[i] = unmask_cost[i] = 0;
	}
	
	protected abstract void addProbeCost (int id);
	
	protected abstract void addProbeCost (int id[], int start, int end);
	
	protected abstract void addSpotCost (int row, int col);
	
	protected abstract double computeMinDistance (int id);
	
	protected abstract double reembedOptimally (int id);
	
	protected double computeMatrix (double matrix[][], char probe[])
	{
		double	cost, un_cost;
		int 	r, c, start;
		
		matrix[0][0] = 0;
		for (c = (start = 1); c <= embed_len; c++)
			matrix[0][c] = matrix[0][c - 1] + pos_mult[0] * mask_cost[c -1];
		
		for (r = 1; r <= probe_len; r++)
		{
			matrix[r][start - 1] = Double.POSITIVE_INFINITY;
			
			for (c = start; c <= embed_len; c++)
			{
				cost = matrix[r][c - 1] + pos_mult[r] * mask_cost[c -1];
				
				if (probe[r] == dep_seq[c])
				{
					un_cost = matrix[r - 1][c - 1] + unmask_cost[c -1];
					
					if (un_cost < cost)
						cost = un_cost;
				}
				
				matrix[r][c] = cost;
				
				if (Double.isInfinite(cost)) start++;
			}
		}
		
		return matrix[probe_len][embed_len];
	}

	protected static abstract class Simple extends OptimumEmbedding
	{
		protected double matrix[][];
		
		protected char probe[];

		protected Simple (SimpleChip chip)
		{
			super (chip, chip.getProbeLength());

			this.probe = new char [probe_len + 1];
			this.matrix = new double [probe_len + 1][embed_len + 1];
		}

		@Override
		protected double computeMinDistance (int id)
		{
			decodeEmbedding (id);
						
			return computeMatrix (matrix, probe);
		}

		@Override
		protected double reembedOptimally (int id)
		{
			decodeEmbedding (id);
			
			double d = computeMatrix (matrix, probe);
			
			encodeEmbedding (id);
			
			return d;
		}

		protected void decodeEmbedding (int id)
		{
			int i, pos, word, mask;
			
			for (i = 1, mask = 0, word = - 1, pos = 0; pos < embed_len; pos++)
			{
				if (pos % Integer.SIZE == 0)
				{
					word++;
					mask = 0x01 << (Integer.SIZE - 1);
				}
				else
					mask >>>= 1;
				
				if ((chip.embed[id][word] & mask) != 0)
					probe[i++] = chip.dep_seq[pos];
			}
		}

		protected void encodeEmbedding (int id)
		{
			int r, c, pos, word, mask;
			
			pos = embed_len;
			word = pos / Integer.SIZE;
			mask = 0x01 << (Integer.SIZE - 1 - (pos % Integer.SIZE));
			
			for (r = probe_len, c = embed_len; pos > 0; c--)
			{
				if ((pos-- % Integer.SIZE) == 0)
				{
					word--;
					mask = 0x01;
				}
				else
					mask <<= 1;
				
				if (r == 0)
				{
					chip.embed[id][word] &= ~mask;
				}
				else if (probe[r] == dep_seq[c] &&
						matrix[r][c] == matrix[r-1][c-1] + unmask_cost[c -1])
				{
					chip.embed[id][word] |= mask;
					r--;
				}
				else
				{
					chip.embed[id][word] &= ~mask;
				}
			}
		}

		protected static class BorderLength extends Simple
		{
			public BorderLength (SimpleChip chip)
			{
				super (chip);
				
				// the position multipliers (only) depend on the length
				// of the probes (e.g., it needs to be computed only once)
				for (int b = 0; b <= probe_len; b++)
					pos_mult[b] = 1;
			}

			@Override
			protected void addProbeCost (int id)
			{
				int word, pos, mask = 0;
				
				for (word = -1, pos = 0; pos < embed_len; pos++)
				{
					if (pos % Integer.SIZE == 0)
					{
						word++;
						mask = 0x01 << (Integer.SIZE-1);
					}
					else
						mask >>>= 1;
					
					if ((chip.embed[id][word] & mask) != 0)
						mask_cost[pos] += 1;
					else
						unmask_cost[pos] += 1;
				}
			}

			@Override
			protected void addProbeCost (int id[], int start, int end)
			{
				int i, word, pos, mask = 0;
				
				for (word = -1, pos = 0; pos < embed_len; pos++)
				{
					if (pos % Integer.SIZE == 0)
					{
						word++;
						mask = 0x01 << (Integer.SIZE-1);
					}
					else
						mask >>>= 1;
					
					for (i = start; i <= end; i++)
					{
						if ((chip.embed[id[i]][word] & mask) != 0)
							mask_cost[pos] += 1;
						else
							unmask_cost[pos] += 1;
					}
				}
			}

			@Override
			protected void addSpotCost (int row, int col)
			{
				int r_2, c_2, id_2;
				
				// top
				r_2 = row - 1;
				if (r_2 >= chip_region.first_row)
					if ((id_2 = chip.spot[r_2][col]) != Chip.EMPTY_SPOT)
						addProbeCost(id_2);
				
				// bottom
				r_2 = row + 1;
				if (r_2 <= chip_region.last_row)
					if ((id_2 = chip.spot[r_2][col]) != Chip.EMPTY_SPOT)
						addProbeCost(id_2);

				// left
				c_2 = col - 1;
				if (c_2 >= chip_region.first_col)
					if ((id_2 = chip.spot[row][c_2]) != Chip.EMPTY_SPOT)
						addProbeCost(id_2);

				// right
				c_2 = col + 1;
				if (c_2 <= chip_region.last_col)
					if ((id_2 = chip.spot[row][c_2]) != Chip.EMPTY_SPOT)
						addProbeCost(id_2);
			}
		}

		protected static class ConflictIndex extends Simple
		{
			protected ConflictIndex (SimpleChip chip)
			{
				super (chip);

				// the position multipliers (only) depend on the length
				// of the probes (e.g., it needs to be computed only once) 
				for (int b = 0; b <= probe_len; b++)
					pos_mult[b] = LayoutEvaluation.positionMultiplier (
															probe_len, b);
			}
			
			@Override
			protected void addProbeCost (int id)
			{
				addProbeCost (id, 1);
			}

			@Override
			protected void addProbeCost (int id[], int start, int end)
			{
				for (int p = start; p <= end; p++)
					addProbeCost (id[p], 1);
			}

			@Override
			protected void addSpotCost (int r_1, int c_1)
			{
				int		dim, r_2, c_2, r_min, c_min, r_max, c_max, id;
				double	w;
				
				dim = LayoutEvaluation.DIM_CONFLICT_REGION;
				
				r_min = Math.max(r_1 - dim, chip_region.first_row);
				c_min = Math.max(c_1 - dim, chip_region.first_col);
				r_max = Math.min(r_1 + dim, chip_region.last_row);
				c_max = Math.min(c_1 + dim, chip_region.last_col);
				
				for (r_2 = r_min; r_2 <= r_max; r_2++)
					for (c_2 = c_min; c_2 <= c_max; c_2++)
					{
						if (r_1 == r_2 && c_1 == c_2)
							continue;
						
						if ((id = chip.spot[r_2][c_2]) == Chip.EMPTY_SPOT)
							continue;
						
						w = LayoutEvaluation.WEIGHT_DIST[dim + r_1 - r_2]
						                                 [dim + c_1 - c_2];
						if (w > 0) addProbeCost (id, w);
					}
			}
			
			private void addProbeCost (int id, double weight_dist)
			{
				int		b, pos, word, mask = 0;
				
				for (b = 0, word = -1, pos = 0; pos < embed_len; pos++)
				{
					if (pos % Integer.SIZE == 0)
					{
						word++;
						mask = 0x01 << (Integer.SIZE-1);
					}
					else
						mask >>>= 1;
					
					if ((chip.embed[id][word] & mask) != 0)
					{
						mask_cost[pos] += weight_dist;
						b++;
					}
					else
						unmask_cost[pos] += weight_dist *
							LayoutEvaluation.positionMultiplier(probe_len, b);
				}
			}
		}
	}

	protected static abstract class Affy extends OptimumEmbedding
	{
		protected double matrix_1[][];
		
		protected double matrix_2[][];

		protected char probe_1[];

		protected char probe_2[];
		
		protected Affy (AffymetrixChip chip)
		{
			// with AffymetrixChips we work with "combined" PM-MM probes 
			// whose length is the length of the individual probes plus 1
			// because they differ (only) in the middle base
			super (chip, chip.getProbeLength() + 1);
			
			// since there are two possible combined probes, we need
			// to create two probe arrays and two matrices
			this.probe_1 = new char [probe_len + 1];
			this.probe_2 = new char [probe_len + 1];
			this.matrix_1 = new double [probe_len + 1][embed_len + 1];
			this.matrix_2 = new double [probe_len + 1][embed_len + 1];
		}

		@Override
		protected double computeMinDistance (int id)
		{
			double	d1, d2;
			
			decodeEmbedding (id);
			
			d1 = computeMatrix (matrix_1, probe_1);
			d2 = computeMatrix (matrix_2, probe_2);
			
			return (d1 <= d2 ? d1 : d2);
		}

		@Override
		protected double reembedOptimally (int id)
		{
			double	d1, d2;
			
			decodeEmbedding (id);

			d1 = computeMatrix (matrix_1, probe_1);
			d2 = computeMatrix (matrix_2, probe_2);
			
			if (d1 <= d2)
			{
				encodeEmbedding (id, id + 1, matrix_1, probe_1);
				return d1;
			}
			
			// else
			encodeEmbedding (id + 1, id, matrix_2, probe_2);
			return d2;
		}

		protected void decodeEmbedding (int id)
		{	
			int		i, pos, word, mask, middle_base;
			char	base, comp;

			// position where the middle base appears
			middle_base = AffymetrixChip.AFFY_MIDDLE_BASE;

			// decode embedding up to (but not including middle base) 
			for (i = 1, mask = 0, word = - 1, pos = 0; i < middle_base; pos++)
			{
				if (pos % Integer.SIZE == 0)
				{
					word++;
					mask = 0x01 << (Integer.SIZE - 1);
				}
				else
					mask >>>= 1;

				if ((chip.embed[id][word] & mask) != 0)
				{
					probe_1[i]   = chip.dep_seq[pos];
					probe_2[i++] = chip.dep_seq[pos];
				}
			}

			// decode middle bases 
			for (; i == middle_base; pos++)
			{
				if (pos % Integer.SIZE == 0)
				{
					word++;
					mask = 0x01 << (Integer.SIZE - 1);
				}
				else
					mask >>>= 1;

				if ((chip.embed[id][word] & mask) != 0)
				{
					base = chip.dep_seq[pos];
					comp = AffymetrixChip.getBaseComplement(base);
					
					probe_1[i] = base;
					probe_2[i++] = comp;
					probe_1[i] = comp;
					probe_2[i++] = base;
				}
			}

			// decode the remaining part 
			for (; pos < embed_len; pos++)
			{
				if (pos % Integer.SIZE == 0)
				{
					word++;
					mask = 0x01 << (Integer.SIZE - 1);
				}
				else
					mask >>>= 1;

				if ((chip.embed[id][word] & mask) != 0)
				{
					probe_1[i]   = chip.dep_seq[pos];
					probe_2[i++] = chip.dep_seq[pos];
				}
			}
		}

		protected void encodeEmbedding (int id_1, int id_2, double matrix[][],
				char probe[])
		{
			int r, c, pos, word, mask;
			
			pos = embed_len;
			word = pos / Integer.SIZE;
			mask = 0x01 << (Integer.SIZE - 1 - (pos % Integer.SIZE));
			
			for (r = probe_len, c = embed_len; pos > 0; c--)
			{
				if ((pos-- % Integer.SIZE) == 0)
				{
					word--;
					mask = 0x01;
				}
				else
					mask <<= 1;
				
				if (r == 0)
				{
					chip.embed[id_1][word] &= ~mask;
					chip.embed[id_2][word] &= ~mask;
				}
				else if (probe[r] == dep_seq[c] &&
						matrix[r][c] == matrix[r-1][c-1] + unmask_cost[c -1])
				{
					chip.embed[id_2][word] |= mask;
					chip.embed[id_1][word] |= mask;
					
					if (r == AffymetrixChip.AFFY_MIDDLE_BASE)
						chip.embed[id_1][word] &= ~mask;
					
					else if (r == AffymetrixChip.AFFY_MIDDLE_BASE + 1)
						chip.embed[id_2][word] &= ~mask;
					
					r--;
				}
				else
				{
					chip.embed[id_1][word] &= ~mask;
					chip.embed[id_2][word] &= ~mask;
				}
			}
		}

		protected static class BorderLength extends Affy
		{
			protected BorderLength (AffymetrixChip chip)
			{
				super (chip);
				
				// the position multipliers (only) depend on the length
				// of the probes (e.g., it needs to be computed only once)
				for (int b = 0; b <= probe_len; b++)
					pos_mult[b] = 1;
			}
			
			@Override
			protected void addProbeCost (int id_1)
			{
				int id_2, word, pos, mask = 0;
				
				if (((AffymetrixChip) chip).isPMProbe(id_1))
					id_2 = id_1 + 1;
				else
					id_2 = id_1 - 1;
				
				for (word = -1, pos = 0; pos < embed_len; pos++)
				{
					if (pos % Integer.SIZE == 0)
					{
						word++;
						mask = 0x01 << (Integer.SIZE-1);
					}
					else
						mask >>>= 1;
					
					if ((chip.embed[id_1][word] & mask) != 0 ||
						(chip.embed[id_2][word] & mask) != 0)
						mask_cost[pos] += 1;
					else
						unmask_cost[pos] += 1;
				}
			}

			@Override
			protected void addProbeCost (int id[], int start, int end)
			{
				int i, id_1, id_2, word, pos, mask = 0;
				
				for (word = -1, pos = 0; pos < embed_len; pos++)
				{
					if (pos % Integer.SIZE == 0)
					{
						word++;
						mask = 0x01 << (Integer.SIZE-1);
					}
					else
						mask >>>= 1;
					
					for (i = start; i <= end; i++)
					{
						id_1 = id[i];
						
						if (((AffymetrixChip) chip).isPMProbe(id_1))
							id_2 = id_1 + 1;
						else
							id_2 = id_1 - 1;

						if ((chip.embed[id_1][word] & mask) != 0 ||
							(chip.embed[id_2][word] & mask) != 0)
							mask_cost[pos] += 1;
						else
							unmask_cost[pos] += 1;
					}
				}
			}

			@Override
			protected void addSpotCost (int row, int col)
			{
				int r_2, c_2, id_2;
				
				// top
				r_2 = row - 2;
				if (r_2 >= chip_region.first_row)
					if ((id_2 = chip.spot[r_2][col]) != Chip.EMPTY_SPOT)
						addProbeCost(id_2);
				
				// bottom
				r_2 = row + 2;
				if (r_2 <= chip_region.last_row)
					if ((id_2 = chip.spot[r_2][col]) != Chip.EMPTY_SPOT)
						addProbeCost(id_2);

				// left
				c_2 = col - 1;
				if (c_2 >= chip_region.first_col)
					if ((id_2 = chip.spot[row][c_2]) != Chip.EMPTY_SPOT)
						addProbeCost(id_2);

				// right
				c_2 = col + 1;
				if (c_2 <= chip_region.last_col)
					if ((id_2 = chip.spot[row][c_2]) != Chip.EMPTY_SPOT)
						addProbeCost(id_2);
			}
		}

		protected static class ConflictIndex extends Affy
		{
			protected ConflictIndex (AffymetrixChip chip)
			{
				super (chip);
				
				int b;
				
				// the position multipliers (only) depend on the length
				// of the probes (e.g., it needs to be computed only once)
				
				// first, bases up to (but not including) the middle base
				for (b = 0; b < AffymetrixChip.AFFY_MIDDLE_BASE; b++)
					pos_mult[b] = LayoutEvaluation.positionMultiplier (
															probe_len, b);
				
				// then, the middle bases (the combined probe has two) 
				pos_mult[b] = LayoutEvaluation.positionMultiplier (
															probe_len, b);
				pos_mult[b + 1] = pos_mult[b]; 
				
				// finally, the remaining bases
				for (b += 2; b <= probe_len; b++)
					pos_mult[b] = LayoutEvaluation.positionMultiplier (
															probe_len, b - 1);
			}
			
			@Override
			protected void addProbeCost (int id)
			{
				addProbeCost (id, 1);
			}

			@Override
			protected void addProbeCost (int id[], int start, int end)
			{
				for (int p = start; p <= end; p++)
					addProbeCost (id[p], 1);
			}

			@Override
			protected void addSpotCost (int r_1, int c_1)
			{
				int		dim, r_2, c_2, r_min, c_min, r_max, c_max, id;
				double	w;
				
				dim = LayoutEvaluation.DIM_CONFLICT_REGION;
				
				r_min = Math.max(r_1 - 2 * dim, chip_region.first_row);
				r_max = Math.min(r_1 + 2 * dim, chip_region.last_row);
				c_min = Math.max(c_1 - dim, chip_region.first_col);
				c_max = Math.min(c_1 + dim, chip_region.last_col);
				
				for (r_2 = r_min; r_2 <= r_max; r_2 += 2)
					for (c_2 = c_min; c_2 <= c_max; c_2++)
					{
						if (r_1 == r_2 && c_1 == c_2)
							continue;
						
						if ((id = chip.spot[r_2][c_2]) == Chip.EMPTY_SPOT)
							continue;
						
						w = LayoutEvaluation.WEIGHT_DIST[2 * dim + r_1 - r_2]
						                                 [dim + c_1 - c_2];
						if (w > 0) addProbeCost (id, w);
					}
			}
			
			private void addProbeCost (int id_1, double weight_dist)
			{
				int		id_2, b, pos, word, mask = 0;
				boolean	pm, mm;
				
				if (((AffymetrixChip) chip).isPMProbe(id_1))
					id_2 = id_1 + 1;
				else
					id_2 = id_1 - 1;
				
				for (b = 0, word = -1, pos = 0; pos < embed_len; pos++)
				{
					if (pos % Integer.SIZE == 0)
					{
						word++;
						mask = 0x01 << (Integer.SIZE-1);
					}
					else
						mask >>>= 1;
					
					pm = (chip.embed[id_1][word] & mask) != 0;
					mm = (chip.embed[id_2][word] & mask) != 0;
					
					if (pm || mm)
					{
						mask_cost[pos] += weight_dist;
						
						if (pm) b++;
					}
					else
					{
						unmask_cost[pos] += weight_dist *
						   LayoutEvaluation.positionMultiplier(probe_len -1, b);
					}
				}
			}
		}
	}
}

