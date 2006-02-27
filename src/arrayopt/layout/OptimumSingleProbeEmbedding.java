/*
 * OptimumSingleProbeEmbedding.java
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
 * �rrayOpt is distributed in the hope that it will be useful, but WITHOUT ANY
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
 * probe. The implementation is based on algorithms described in the following
 * papers:
 * 
 * <P><UL>
 * <LI>"Border Length Minimization in DNA Array Design", A. Kahng, I. Mandoiu,
 * P. Pevzner, S. Reda, and A. Zelikovsky, Proc. 2nd Int. Workshop on Algorithms
 * in Bioinformatics (WABI 2002), pp. 435�448.
 * <LI>"Engineering a scalable placement heuristic for DNA probe arrays",
 * A.B. Kahng, I. Mandoiu, P. Pevzner, S. Reda, and A. Zelikovsky, Proc. 7th
 * Annual Int. Conference on Research in Computational Molecular Biology
 * (RECOMB 2003), pp. 148-83.
 * </UL></P>
 * 
 * <P>This implementaion can work with two types of conflict minimization:
 * border length ({@link #BORDER_LENGTH_MIN}) or conflict index
 * ({@link #CONFLICT_INDEX_MIN}). These are refered to as distance or quality
 * measures.</P>
 * 
 * <P>Two main functionalities are provided. Firstly, it can compute the minimum
 * distance that any embedding of a probe can have to the current embedding of
 * another probe. (The distance depends on the chosen quality measure and whether
 * the location of the probes on the chip is taken into account.) This is
 * accomplished by the {@link minDistanceProbe} and {@link minDistanceSpot}
 * methods (the first only considers the embeddings, while the latter also
 * considers the location of the probes on the chip).</P>
 * 
 * <P>The other main feature of this class is provided by the
 * {@link reembedProbe} and {@link reembedSpot} methods, which re-embeds a
 * probe with the minimum distance.</P>
 *
 * <P>In fact, these methods have several signatures that allow the computed
 * distance to take into account not only a single probe but also a set of
 * probes or all probes around a given spot (which probes are considered depends
 * on the quality measure).</P>
 * 
 * <P>Note that the current embeddings of the probes towards which the minimum
 * distance is sought are considered fixed. This allows the use of a simple and
 * efficient dynamic programming approach with quadratic space and time
 * complexities. The implementation aimed at maximizing the reuse of common
 * routines by the sub-classes rather than achiving maximum speed.</P>
 *  
 * <P>This class is abstract and the implementation details are provided by
 * (inner) sub-classes for each quality measure (as well as for different types
 * of chips). These sub-classes, however, cannot be instantiated directly.
 * Instead, the {@link #createEmbedder} factory method should be called with a
 * chip instance and the desired type of minimization (either
 * {@link #BORDER_LENGTH_MIN} or {@link #CONFLICT_INDEX_MIN}) passed as
 * arguments. The returned object will provide the functionality specified
 * here.</P>
 *
 * <P>This class is intended for use by other classes inside the layout package
 * and, therefore, has no public methods.</P>
 * 
 * <P><B>Implementation note:</B> after the dynamic programming matrix is
 * computed, an optimal embedding of a probe is retrieved by tracing back a path
 * in the matrix from the bottom right cell to the top left one. In some cases
 * several embeddings might be optimal which means that several paths can be
 * traced in the matrix. The current implementation retrieves the left-most
 * optimum embedding. It may be interesting to allow the user to configure the
 * algorithm to return a different type of embedding, for instance, a random or
 * the right-most optimum embedding. The random case might be useful to generate
 * a better distribution of the unmasked steps. The left-most and right-most
 * cases might give better results depending on the layout of the chip.</B></P>   
 *
 * @author Sergio A. de Carvalho Jr.
 */
public abstract class OptimumSingleProbeEmbedding
{
	/**
	 * Reference to the region of the chip for which this object is configured
	 * to work with.
	 */
	protected RectangularRegion chip_region;

	/**
	 * The length of the embeddings found on the chip (i.e the length of the
	 * deposition sequente). 
	 */
	protected int embed_len;
	
	/**
	 * The length of probe sequences found on the chip. 
	 */
	protected int probe_len;
	
	/**
	 * Constant that indicates that the border length should be considered when
	 * computing the distance between embeddings/spots. It is used to create a
	 * new instance of this class with the {@link #createEmbedder(Chip, int)}
	 * method.
	 * 
	 * <P>Border length minimization means that only immediate neighbors of a
	 * spot are considered when {@link minDistanceSpot} or {@link reembedSpot}
	 * are called. Also, any conflicts between the probes have the same
	 * weight.</P>
	 */
	public static final int BORDER_LENGTH_MIN = 0;
	
	/**
	 * Constant that indicates that the conflict index should be used when
	 * computing the distance between embeddings/spots. It is used to create a
	 * new instance of this class with the {@link #createEmbedder(Chip, int)}
	 * method.
	 * 
	 * <P>Conflict index minimization means that not only immediate neighbors
	 * but also neighbors "close to" the spot are considered when
	 * {@link minDistanceSpot} or {@link reembedSpot} are called. Also, the
	 * conflicts between the probes are weighted according to the position of
	 * the base in the probe sequence which could be damaged by the conflict.
	 * The exact values such weights are defined by the
	 * {@linkplain LayoutEvaluation} class.
	 */
	public static final int CONFLICT_INDEX_MIN = 1;
	
	/**
	 * This constructor can only be accessed internally. For instantiating
	 * objects of this class and sub-classes, use {@link createEmbedder}.  
	 * @param chip a chip instance
	 * @param probe_len the length of the probes 
	 */
	protected OptimumSingleProbeEmbedding (Chip chip)
	{
		this.probe_len = chip.getProbeLength();
		this.embed_len = chip.getEmbeddingLength();
		this.chip_region = chip.getChipRegion();
	}

	/**
	 * Create an instance of the appropriate sub-class according to the type
	 * of chip and the default minimization function
	 * ({@link #BORDER_LENGTH_MIN}).
	 * 
	 * @param c a chip instance  
	 * @param mode type of the desired distance measure
	 * @return an instance of an OptimumSingleProbeEmbedding for the type of
	 * chip and the selected mode 
	 */
	static OptimumSingleProbeEmbedding createEmbedder (Chip c)
	{
		return createEmbedder (c, BORDER_LENGTH_MIN);
	}

	/**
	 * Create an instance of the Optimum Single-Probe Embedding (OSPE) algorithm
	 * with the desired type of minimization function
	 * ({@link #BORDER_LENGTH_MIN} or {@link #CONFLICT_INDEX_MIN}).
	 * 
	 * @param c a chip instance  
	 * @param mode type of the desired distance measure
	 * @return an instance of an OptimumSingleProbeEmbedding for the type of
	 * chip and the selected mode 
	 */
	static OptimumSingleProbeEmbedding createEmbedder (Chip c, int mode)
	{
		if (mode != BORDER_LENGTH_MIN && mode != CONFLICT_INDEX_MIN)
			throw new IllegalArgumentException
				("Unknown distance mode: " + mode);
			
		if (c instanceof SimpleChip)
		{
			if (mode == BORDER_LENGTH_MIN)
				return new Simple.BorderLength ((SimpleChip) c);
			
			// else: CONFLICT_INDEX_MIN
				return new Simple.ConflictIndex ((SimpleChip) c);
		}
		
		if (c instanceof AffymetrixChip)
		{
			if (mode == BORDER_LENGTH_MIN)
				return new Affymetrix.BorderLength ((AffymetrixChip) c);
			
			// else: CONFLICT_INDEX_MIN
				return new Affymetrix.ConflictIndex ((AffymetrixChip) c);
		}
		
		// else
		throw new IllegalArgumentException ("Unsupported chip type.");
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
	 * @param id_2 the IDs of the set of probe with a fixed embedding
	 * @return the minimum distance between id_1 and the set of probes
	 */
	double minDistanceProbe (int id_1, int id_2[])
	{
		return minDistanceProbe (id_1, id_2, 0, id_2.length - 1);
	}

	/**
	 * Computes the minimum distance between any valid embedding of a probe
	 * (<CODE>id_1</CODE>) and the current embeddings of a given set of probes.
	 * @param id_1 the ID of the first probe
	 * @param id_2 the IDs of the set of probe with a fixed embedding
	 * @param start the index of the first element on the list
	 * @param end the index of the last element on the list 
	 * @return the minimum distance between id_1 and the set of probes
	 */
	double minDistanceProbe (int id_1, int id_2[], int start, int end)
	{
		resetCosts();
		addProbeCost (id_2, start, end);
		return computeMinDistance(id_1);
	}

	/**
	 * Computes the minimum distance between any valid embedding of the current
	 * probe of a spot and the embeddings of the neighboring probes.
	 * @param row spot's row coordinate
	 * @param col spot's column coordinate
	 * @return the minimum conflict that any embedding of the probe on the spot
	 * can cause
	 */
	double minDistanceSpot (int row, int col)
	{
		int id;
		
		if ((id = getProbeID(row, col)) == Chip.EMPTY_SPOT)
			return 0;
		
		resetCosts();
		addSpotCost (row, col);
		return computeMinDistance(id);
	}

	/**
	 * Computes the minimum distance between any valid embedding of a probe,
	 * when placed on a particular spot, and the current embeddings of the
	 * neighboring probes.
	 * @param row spot's row coordinate
	 * @param col spot's column coordinate
	 * @param id probe ID
	 * @return the minimum conflict that the probe spot can cause when placed
	 * on the spot
	 */
	double minDistanceSpot (int row, int col, int id)
	{
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
	 * @param id_2 the IDs of the set of probe with a fixed embedding
	 * @return the minimum distance between id_1 and the set of probes
	 */
	double reembedProbe (int id_1, int id_2[])
	{
		return reembedProbe (id_1, id_2, 0, id_2.length - 1); 
	}

	/**
	 * Reembeds a probe (<CODE>id_1</CODE>) so that its distance to the current
	 * embeddings of the given set of probes is minimal.
	 * @param id_1 the ID of the first probe
	 * @param id_2 the IDs of the set of probe with a fixed embedding
	 * @param start the index of the first element on the list
	 * @param end the index of the last element on the list 
	 * @return the minimum distance between id_1 and the set of probes
	 */
	double reembedProbe (int id_1, int id_2[], int start, int end)
	{
		resetCosts();
		addProbeCost (id_2, start, end);
		return reembedOptimally (id_1);
	}

	/**
	 * Reembeds a spot's current probe so that the conflicts that it generates
	 * in the region are minimized.
	 * is minimal.
	 * @param row spot's row coordinate
	 * @param col spot's column coordinate
	 * @return the minimum conflict that the probe on the given spot can cause
	 * to the neighboring probes or zero if the spot is empty
	 */
	double reembedSpot (int row, int col)
	{
		int id;
		
		if ((id = getProbeID(row, col)) == Chip.EMPTY_SPOT)
			return 0;

		resetCosts();
		addSpotCost (row, col);
		return reembedOptimally (id);
	}

	/**
	 * Reembeds a probe so that the conflicts that it generates when placed on
	 * a given spot are minimized.
	 * is minimal.
	 * @param id probe ID
	 * @param row spot's row coordinate
	 * @param col spot's column coordinate
	 * @return the minimum conflict that the probe can cause to the neighboring
	 * probes when placed on the spot
	 */
	double reembedSpot (int row, int col, int id)
	{
		resetCosts();
		addSpotCost (row, col);
		return reembedOptimally (id);
	}

	/**
	 * TODO: move this function to another class.
	 * Computes the number of ways in which a probe can be embedded into the
	 * deposition sequence. 
	 * @param id probe ID
	 * @return the number of different embeddings the probe can have
	 */
	abstract int numberOfEmbeddings (int id);
	
	/**
	 * Reset cost arrays before a new distance is computed.
	 */
	protected abstract void resetCosts ();
	
	protected abstract int getProbeID (int row, int col);
	
	protected abstract void addProbeCost (int id);
	
	protected abstract void addProbeCost (int id[], int start, int end);
	
	protected abstract void addSpotCost (int row, int col);
	
	protected abstract double computeMinDistance (int id);
	
	protected abstract double reembedOptimally (int id);

	protected int numberOfEmbeddings (char probe[], char dep_seq[], int m[])
	{
		int r, c, top, last_row, tmp;
		
		for (r = 0; r < probe.length; r++)
			m[r] = 0;
		
		last_row = 0;
		
		for (c = 0; c < dep_seq.length; c++)
		{
			top = 1;
			
			for (r = 0; r <= last_row; r++)
			{
				if (probe[r] == dep_seq[c])
				{
					tmp = m[r];
					m[r] += top;
					top = tmp;
				}
				else
					top = m[r];
			}
			
			if (m[last_row] > 0)
				if (last_row < probe.length - 1)
					last_row++;
		}
		
		return m[probe.length - 1];
	}

	protected static abstract class Simple extends OptimumSingleProbeEmbedding
	{
		protected SimpleChip chip;

		protected double matrix[][];
		
		protected char probe[];

		protected double mask_cost[];
		
		protected double unmask_cost[];
		
		protected int num_embed[];
		
		protected Simple (SimpleChip chip)
		{
			super (chip);

			this.chip = chip;
			this.matrix = new double [probe_len + 1][embed_len + 1];
			
			this.probe = new char [probe_len];
			this.num_embed = new int [probe_len];
			
			this.mask_cost = new double [embed_len];
			this.unmask_cost = new double [embed_len];
		}
		
		@Override
		protected void resetCosts ()
		{
			for (int i = 0; i < embed_len; i++)
				mask_cost[i] = unmask_cost[i] = 0;
		}
		
		@Override
		protected int getProbeID (int row, int col)
		{
			return chip.spot[row][col];
		}

		@Override
		protected double computeMinDistance (int id)
		{
			decodeEmbedding (id);
			return computeMatrix ();
		}

		@Override
		protected double reembedOptimally (int id)
		{
			double d;
			
			decodeEmbedding (id);
			d = computeMatrix ();
			encodeEmbedding (id);
			return d;
		}

		@Override
		int numberOfEmbeddings (int id)
		{
			decodeEmbedding (id);
			
			return numberOfEmbeddings (probe, chip.dep_seq, num_embed);
		}

		protected void decodeEmbedding (int id)
		{
			int i, pos, word, bitmask = 0;
			
			for (i = 0, word = - 1, pos = 0; pos < embed_len; pos++)
			{
				if (pos % Integer.SIZE == 0)
				{
					word++;
					bitmask = 0x01 << (Integer.SIZE - 1);
				}
				else
					bitmask >>>= 1;
				
				if ((chip.embed[id][word] & bitmask) != 0)
					probe[i++] = chip.dep_seq[pos];
			}
		}
		
		protected abstract double computeMatrix ();
		
		protected abstract void encodeEmbedding (int id);

		private static class BorderLength extends Simple
		{
			public BorderLength (SimpleChip chip)
			{
				super (chip);
			}

			@Override
			protected void addProbeCost (int id)
			{
				int word, pos, bitmask = 0;
				
				for (word = -1, pos = 0; pos < embed_len; pos++)
				{
					if (pos % Integer.SIZE == 0)
					{
						word++;
						bitmask = 0x01 << (Integer.SIZE-1);
					}
					else
						bitmask >>>= 1;
					
					if ((chip.embed[id][word] & bitmask) != 0)
						mask_cost[pos] += 1;
					else
						unmask_cost[pos] += 1;
				}
			}

			@Override
			protected void addProbeCost (int id[], int start, int end)
			{
				int i, word, pos, bitmask = 0;
				
				for (word = -1, pos = 0; pos < embed_len; pos++)
				{
					if (pos % Integer.SIZE == 0)
					{
						word++;
						bitmask = 0x01 << (Integer.SIZE-1);
					}
					else
						bitmask >>>= 1;
					
					for (i = start; i <= end; i++)
					{
						if ((chip.embed[id[i]][word] & bitmask) != 0)
							mask_cost[pos] += 1;
						else
							unmask_cost[pos] += 1;
					}
				}
			}

			@Override
			protected void addSpotCost (int row, int col)
			{
				int r, c, id;
				
				// top
				r = row - 1;
				if (r >= chip_region.first_row)
					if ((id = chip.spot[r][col]) != Chip.EMPTY_SPOT)
						addProbeCost(id);
				
				// bottom
				r = row + 1;
				if (r <= chip_region.last_row)
					if ((id = chip.spot[r][col]) != Chip.EMPTY_SPOT)
						addProbeCost(id);

				// left
				c = col - 1;
				if (c >= chip_region.first_col)
					if ((id = chip.spot[row][c]) != Chip.EMPTY_SPOT)
						addProbeCost(id);

				// right
				c = col + 1;
				if (c <= chip_region.last_col)
					if ((id = chip.spot[row][c]) != Chip.EMPTY_SPOT)
						addProbeCost(id);
			}
			
			@Override
			protected double computeMatrix ()
			{
				double	mask, unmask;
				int 	r, c, start;
				
				matrix[0][0] = 0;
				for (c = (start = 1); c <= embed_len; c++)
					matrix[0][c] = matrix[0][c - 1] + mask_cost[c -1];
				
				for (r = 1; r <= probe_len; r++, start++)
				{
					matrix[r][start - 1] = Double.POSITIVE_INFINITY;
					
					for (c = start; c <= embed_len; c++)
					{
						mask = matrix[r][c - 1] + mask_cost[c -1];
						
						if (probe[r - 1] == chip.dep_seq[c - 1])
							unmask = matrix[r - 1][c - 1] + unmask_cost[c -1];
						else
							unmask = Double.POSITIVE_INFINITY;
						
						matrix[r][c] = Math.min(mask, unmask);
						
						if (Double.isInfinite(matrix[r][c])) start++;
					}
				}
				
				return matrix[probe_len][embed_len];
			}

			@Override
			protected void encodeEmbedding (int id)
			{
				int r, c, pos, word, bitmask;
				
				pos = embed_len;
				word = pos / Integer.SIZE;
				bitmask = 0x01 << (Integer.SIZE - 1 - (pos % Integer.SIZE));
				
				for (r = probe_len, c = embed_len; pos > 0; c--)
				{
					if ((pos-- % Integer.SIZE) == 0)
					{
						word--;
						bitmask = 0x01;
					}
					else
						bitmask <<= 1;
					
					chip.embed[id][word] &= ~bitmask;
					
					if (r == 0) continue;

					if (matrix[r][c] == matrix[r][c - 1] + mask_cost[c -1])
						continue;
					
					chip.embed[id][word] |= bitmask;
					r--;
				}
			}
		}
		
		private static class ConflictIndex extends Simple
		{
			private double pos_mult[];
			
			protected ConflictIndex (SimpleChip chip)
			{
				super (chip);
				
				this.pos_mult = new double [probe_len + 1];
				
				for (int b = 0; b <= probe_len; b++)
					pos_mult[b] = LayoutEvaluation.positionMultiplier (
															probe_len, b);
			}
			
			@Override
			protected void addProbeCost (int id)
			{
				addSingleProbeCost (id, 1, 1);
			}
		
			@Override
			protected void addProbeCost (int id[], int start, int end)
			{
				for (int p = start; p <= end; p++)
					addSingleProbeCost (id[p], 1, 1);
			}
		
			@Override
			protected void addSpotCost (int row, int col)
			{
				int		r, c, id, dim, r_min, c_min, r_max, c_max;
				double	mask_w, unmask_w;
				
				dim = LayoutEvaluation.DIM_CONFLICT_REGION;
				
				r_min = Math.max(row - dim, chip_region.first_row);
				c_min = Math.max(col - dim, chip_region.first_col);
				r_max = Math.min(row + dim, chip_region.last_row);
				c_max = Math.min(col + dim, chip_region.last_col);
				
				for (r = r_min; r <= r_max; r++)
					for (c = c_min; c <= c_max; c++)
					{
						if (row == r && col == c)
							continue;
						
						if ((id = chip.spot[r][c]) == Chip.EMPTY_SPOT)
							continue;
						
						mask_w = LayoutEvaluation.WEIGHT_DIST
									[dim + r - row][dim + c - col];
						
						unmask_w = LayoutEvaluation.WEIGHT_DIST
									[dim + row - r][dim + col - c];
						
						if (mask_w > 0 || unmask_w > 0)
							addSingleProbeCost (id, mask_w, unmask_w);
					}
			}
			
			private void addSingleProbeCost (int id, double mask_weight,
					double unmask_weight)
			{
				int		b, pos, word, bitmask = 0;
				
				for (b = 0, word = -1, pos = 0; pos < embed_len; pos++)
				{
					if (pos % Integer.SIZE == 0)
					{
						word++;
						bitmask = 0x01 << (Integer.SIZE-1);
					}
					else
						bitmask >>>= 1;
					
					if ((chip.embed[id][word] & bitmask) != 0)
					{
						mask_cost[pos] += mask_weight;
						b++;
					}
					else
						unmask_cost[pos] += unmask_weight *
							LayoutEvaluation.positionMultiplier(probe_len, b);
				}
			}
			
			@Override
			protected double computeMatrix ()
			{
				double	mask, unmask;
				int 	r, c, start;
				
				matrix[0][0] = 0;
				for (c = (start = 1); c <= embed_len; c++)
					matrix[0][c] = matrix[0][c - 1] +
									pos_mult[0] * mask_cost[c -1];
				
				for (r = 1; r <= probe_len; r++, start++)
				{
					matrix[r][start - 1] = Double.POSITIVE_INFINITY;
					
					for (c = start; c <= embed_len; c++)
					{
						mask = matrix[r][c - 1] + pos_mult[r] * mask_cost[c -1];
						
						if (probe[r - 1] == chip.dep_seq[c - 1])
							unmask = matrix[r - 1][c - 1] + unmask_cost[c -1];
						else
							unmask = Double.POSITIVE_INFINITY;
						
						matrix[r][c] = Math.min(mask, unmask);
						
						if (Double.isInfinite(matrix[r][c])) start++;
					}
				}
				
				return matrix[probe_len][embed_len];
			}
			
			@Override
			protected void encodeEmbedding (int id)
			{
				int r, c, pos, word, bitmask;
				
				pos = embed_len;
				word = pos / Integer.SIZE;
				bitmask = 0x01 << (Integer.SIZE - 1 - (pos % Integer.SIZE));
				
				for (r = probe_len, c = embed_len; pos > 0; c--)
				{
					if ((pos-- % Integer.SIZE) == 0)
					{
						word--;
						bitmask = 0x01;
					}
					else
						bitmask <<= 1;
					
					chip.embed[id][word] &= ~bitmask;
					
					if (r == 0) continue;

					if (matrix[r][c] ==
						matrix[r][c - 1] + pos_mult[r] * mask_cost[c -1])
						continue;
					
					chip.embed[id][word] |= bitmask;
					r--;
				}
			}
		}
	}

	protected static abstract class Affymetrix extends OptimumSingleProbeEmbedding
	{
		protected AffymetrixChip chip;

		protected double matrix_1[][];
		
		protected double matrix_2[][];

		protected char probe_1[];

		protected char probe_2[];
		
		protected double mask_cost_pm[];
		
		protected double mask_cost_mm[];
		
		protected double unmask_cost_pm[];
		
		protected double unmask_cost_mm[];
		
		private int num_embed[];
		
		protected Affymetrix (AffymetrixChip chip)
		{
			super (chip);
			
			this.chip = chip;

			// with AffymetrixChips we work with "combined" PM-MM probes 
			// whose length is the length of the individual probes plus 1
			// because they differ (only) in the middle bases; also,
			// since there are two possible combined probes, we need
			// to create two probe arrays and two matrices
			this.probe_1 = new char [probe_len + 1];
			this.probe_2 = new char [probe_len + 1];			
			this.matrix_1 = new double [probe_len + 2][embed_len + 1];
			this.matrix_2 = new double [probe_len + 2][embed_len + 1];

			// the conflict costs are divided into those
			// generated/suffered by the PM and MM probe 
			this.mask_cost_pm = new double [embed_len];
			this.mask_cost_mm = new double [embed_len];
			this.unmask_cost_pm = new double [embed_len];
			this.unmask_cost_mm = new double [embed_len];

			this.num_embed = new int [probe_len + 1];
		}
		
		@Override
		protected void resetCosts ()
		{
			for (int i = 0; i < embed_len; i++)
			{
				mask_cost_pm[i] = mask_cost_mm[i] = 0;
				unmask_cost_pm[i] = unmask_cost_mm[i] = 0;
			}
		}
		
		@Override
		protected int getProbeID (int row, int col)
		{
			return chip.spot[row][col];
		}
		
		@Override
		protected double computeMinDistance (int id)
		{
			double	d1, d2;
			int		mid = AffymetrixChip.AFFY_MIDDLE_BASE;
			
			decodeEmbedding (id);
			d1 = computeMatrix (matrix_1, probe_1, mid, mid + 1);
			d2 = computeMatrix (matrix_2, probe_2, mid + 1, mid);
			
			return (d1 <= d2 ? d1 : d2);
		}

		@Override
		protected double reembedOptimally (int id)
		{
			double	d1, d2;
			int		mid = AffymetrixChip.AFFY_MIDDLE_BASE;
			
			if (!chip.isPMProbe(id)) id = id -1;

			decodeEmbedding (id);
			d1 = computeMatrix (matrix_1, probe_1, mid, mid + 1);
			d2 = computeMatrix (matrix_2, probe_2, mid + 1, mid);
			
			if (d1 <= d2)
			{
				encodeEmbedding (id, id + 1, matrix_1);
				return d1;
			}
			// else
			encodeEmbedding (id + 1, id, matrix_2);
			return d2;
		}
		
		@Override
		int numberOfEmbeddings (int id)
		{
			decodeEmbedding (id);
			
			return numberOfEmbeddings (probe_1, chip.dep_seq, num_embed) +
				numberOfEmbeddings (probe_2, chip.dep_seq, num_embed);
		}
		
		protected void decodeEmbedding (int id)
		{	
			int		i, pos, word, bitmask = 0;
			char	base, comp;

			for (i = 0, word = - 1, pos = 0; pos < embed_len; pos++)
			{
				if (pos % Integer.SIZE == 0)
				{
					word++;
					bitmask = 0x01 << (Integer.SIZE - 1);
				}
				else
					bitmask >>>= 1;

				if ((chip.embed[id][word] & bitmask) == 0)
					continue;
				
				if (i != AffymetrixChip.AFFY_MIDDLE_BASE - 1)
				{
					probe_1[i]   = chip.dep_seq[pos];
					probe_2[i++] = chip.dep_seq[pos];
				}
				else
				{
					base = chip.dep_seq[pos];
					comp = AffymetrixChip.getBaseComplement(base);
					
					probe_1[i] = base;
					probe_2[i++] = comp;
					probe_1[i] = comp;
					probe_2[i++] = base;
				}
			}
		}
		
		protected abstract double computeMatrix (double matrix[][],
				char probe[], int mid_pm, int mid_mm);

		protected abstract void encodeEmbedding (int id_1, int id_2,
				double matrix[][]);
		
		private static class BorderLength extends Affymetrix
		{
			protected BorderLength (AffymetrixChip chip)
			{
				super (chip);
			}
			
			@Override
			protected void addProbeCost (int id_1)
			{
				int id_2, word, pos, bitmask = 0;
				
				if (chip.isPMProbe(id_1))
					id_2 = id_1 + 1;
				else
					id_2 = id_1 - 1;
				
				for (word = -1, pos = 0; pos < embed_len; pos++)
				{
					if (pos % Integer.SIZE == 0)
					{
						word++;
						bitmask = 0x01 << (Integer.SIZE-1);
					}
					else
						bitmask >>>= 1;
					
					if ((chip.embed[id_1][word] & bitmask) != 0 ||
						(chip.embed[id_2][word] & bitmask) != 0)
					{
						mask_cost_pm[pos] += 1;
						mask_cost_mm[pos] += 1;
					}
					else
					{
						unmask_cost_pm[pos] += 1;
						unmask_cost_mm[pos] += 1;
					}
				}
			}

			@Override
			protected void addProbeCost (int id[], int start, int end)
			{
				int i, id_1, id_2, word, pos, bitmask = 0;
				
				for (word = -1, pos = 0; pos < embed_len; pos++)
				{
					if (pos % Integer.SIZE == 0)
					{
						word++;
						bitmask = 0x01 << (Integer.SIZE-1);
					}
					else
						bitmask >>>= 1;
					
					for (i = start; i <= end; i++)
					{
						id_1 = id[i];
						
						if (chip.isPMProbe(id_1))
							id_2 = id_1 + 1;
						else
							id_2 = id_1 - 1;

						if ((chip.embed[id_1][word] & bitmask) != 0 ||
							(chip.embed[id_2][word] & bitmask) != 0)
						{
							mask_cost_pm[pos] += 1;
							mask_cost_mm[pos] += 1;
						}
						else
						{
							unmask_cost_pm[pos] += 1;
							unmask_cost_mm[pos] += 1;
						}
					}
				}
			}

			@Override
			protected void addSpotCost (int row, int col)
			{
				int r, c, id, pm_row, mm_row;
				
				if (chip.isPMProbe(chip.spot[row][col]))
				{
					pm_row = row;
					mm_row = row + 1;
				}	
				else
				{
					mm_row = row;
					pm_row = row - 1;
				}
				
				// top
				r = pm_row - 1;
				if (r >= chip_region.first_row)
					if ((id = chip.spot[r][col]) != Chip.EMPTY_SPOT)
						addSingleProbeCost(id, 1, 0);
				
				// bottom
				r = mm_row + 1;
				if (r <= chip_region.last_row)
					if ((id = chip.spot[r][col]) != Chip.EMPTY_SPOT)
						addSingleProbeCost(id, 0, 1);

				// left
				c = col - 1;
				if (c >= chip_region.first_col)
				{
					// PM row
					if ((id = chip.spot[pm_row][c]) != Chip.EMPTY_SPOT)
						addSingleProbeCost(id, 1, 0);
					// MM row
					if ((id = chip.spot[mm_row][c]) != Chip.EMPTY_SPOT)
						addSingleProbeCost(id, 0, 1);
				}

				// right
				c = col + 1;
				if (c <= chip_region.last_col)
				{
					// PM row
					if ((id = chip.spot[pm_row][c]) != Chip.EMPTY_SPOT)
						addSingleProbeCost(id, 1, 0);
					if ((id = chip.spot[mm_row][c]) != Chip.EMPTY_SPOT)
						addSingleProbeCost(id, 0, 1);
				}
			}

			private void addSingleProbeCost (int id, int pm_mult, int mm_mult)
			{
				int word, pos, bitmask = 0;
				
				for (word = -1, pos = 0; pos < embed_len; pos++)
				{
					if (pos % Integer.SIZE == 0)
					{
						word++;
						bitmask = 0x01 << (Integer.SIZE-1);
					}
					else
						bitmask >>>= 1;
					
					if ((chip.embed[id][word] & bitmask) != 0)
					{
						mask_cost_pm[pos] += 1 * pm_mult;
						mask_cost_mm[pos] += 1 * mm_mult;
					}
					else
					{
						unmask_cost_pm[pos] += 1 * pm_mult;
						unmask_cost_mm[pos] += 1 * mm_mult;
					}
				}
			}

			@Override
			protected double computeMatrix (double matrix[][], char probe[],
					int mid_pm, int mid_mm)
			{
				double	mask, unmask;
				int 	r, c, start;
				
				matrix[0][0] = 0;
				for (c = (start = 1); c <= embed_len; c++)
					matrix[0][c] = matrix[0][c - 1] + mask_cost_pm[c -1]
									+ mask_cost_mm[c -1];
				
				for (r = 1; r <= probe_len + 1; r++, start++)
				{
					matrix[r][start - 1] = Double.POSITIVE_INFINITY;
					
					for (c = start; c <= embed_len; c++)
					{
						mask = matrix[r][c - 1] + mask_cost_pm[c -1]
								+ mask_cost_mm[c -1];
						
						if (probe[r - 1] != chip.dep_seq[c -1])
						{
							unmask = Double.POSITIVE_INFINITY;
						}
						else
						{
							unmask = matrix[r - 1][c - 1];
							
							if (r == mid_pm)
								unmask += mask_cost_mm[c -1];
							else
								unmask += unmask_cost_mm[c - 1];

							if (r == mid_mm)
								unmask += mask_cost_pm[c -1];
							else
								unmask += unmask_cost_pm[c - 1];
						}
						
						matrix[r][c] = Math.min(mask, unmask);
						
						if (Double.isInfinite(matrix[r][c])) start++;
					}
				}
				
				return matrix[probe_len + 1][embed_len];
			}
			
			@Override
			protected void encodeEmbedding (int id_1, int id_2, double matrix[][])
			{
				int		r, c, pos, word, bitmask;
				double	mask_cost;
				
				pos = embed_len;
				word = pos / Integer.SIZE;
				bitmask = 0x01 << (Integer.SIZE - 1 - (pos % Integer.SIZE));
				
				for (r = probe_len + 1, c = embed_len; pos > 0; c--)
				{
					if ((pos-- % Integer.SIZE) == 0)
					{
						word--;
						bitmask = 0x01;
					}
					else
						bitmask <<= 1;
					
					chip.embed[id_1][word] &= ~bitmask;
					chip.embed[id_2][word] &= ~bitmask;
					
					if (r == 0) continue;
					
					mask_cost = mask_cost_mm[c - 1] + mask_cost_pm[c - 1];
					
					if (matrix[r][c] == matrix[r][c - 1] + mask_cost)
						continue;
					
					if (r != AffymetrixChip.AFFY_MIDDLE_BASE)
						chip.embed[id_2][word] |= bitmask;
					
					if (r != AffymetrixChip.AFFY_MIDDLE_BASE + 1)
						chip.embed[id_1][word] |= bitmask;

					r--;
				}
			}
		}
		
		private static class ConflictIndex extends Affymetrix
		{
			private double pos_mult[];
			
			protected ConflictIndex (AffymetrixChip chip)
			{
				super (chip);
				
				int b;
				
				this.pos_mult = new double [probe_len + 2];
				
				// first, bases up to (but not including) the middle base
				for (b = 0; b < AffymetrixChip.AFFY_MIDDLE_BASE; b++)
					pos_mult[b] = LayoutEvaluation.positionMultiplier (
															probe_len, b);
				
				// then, the middle bases (the combined probe has two) 
				pos_mult[b] = LayoutEvaluation.positionMultiplier (
															probe_len, b);
				pos_mult[b + 1] = pos_mult[b]; 
				
				// finally, the remaining bases
				for (b += 2; b <= probe_len + 1; b++)
					pos_mult[b] = LayoutEvaluation.positionMultiplier (
														probe_len, b - 1);
			}			
		
			@Override
			protected void addProbeCost (int id_1)
			{
				int		id_2, base, pos, word, bitmask = 0;
				boolean	p_1, p_2, middle_base = false;
				double	pos_mul;
				
				if (chip.isPMProbe(id_1))
					id_2 = id_1 + 1;
				else
					id_2 = id_1 - 1;
				
				for (base = 0, word = -1, pos = 0; pos < embed_len; pos++)
				{
					if (pos % Integer.SIZE == 0)
					{
						word++;
						bitmask = 0x01 << (Integer.SIZE-1);
					}
					else
						bitmask >>>= 1;
					
					p_1 = (chip.embed[id_1][word] & bitmask) != 0;
					p_2 = (chip.embed[id_2][word] & bitmask) != 0;
					
					if (p_1 || p_2)
					{
						mask_cost_pm[pos] += 1;
						mask_cost_pm[pos] += 1;
						
						if (p_1 && p_2)
							base++;
						else if (!middle_base)
						{
							base++;
							middle_base = true;
						}
							
					}
					else
					{
						pos_mul = LayoutEvaluation.positionMultiplier(probe_len,
									base);
						unmask_cost_pm[pos] += pos_mul;
						unmask_cost_mm[pos] += pos_mul;
					}
				}
			}

			@Override
			protected void addProbeCost (int id[], int start, int end)
			{
				for (int p = start; p <= end; p++)
					addProbeCost (id[p]);
			}

			@Override
			protected void addSpotCost (int row, int col)
			{
				int		id, pm_row, mm_row, dim, r_min, r_max, c_min, c_max;
				double	mask_pm, mask_mm, unmask_pm, unmask_mm;
				boolean	zero;
				
				if (chip.isPMProbe(chip.spot[row][col]))
				{
					pm_row = row;
					mm_row = row + 1;
				}	
				else
				{
					mm_row = row;
					pm_row = row - 1;
				}

				dim = LayoutEvaluation.DIM_CONFLICT_REGION;
				r_min = Math.max(pm_row - dim, chip_region.first_row);
				r_max = Math.min(mm_row + dim, chip_region.last_row);
				c_min = Math.max(col - dim, chip_region.first_col);
				c_max = Math.min(col + dim, chip_region.last_col);
				
				for (int r = r_min; r <= r_max; r++)
					for (int c = c_min; c <= c_max; c++)
					{
						if ((id = chip.spot[r][c]) == Chip.EMPTY_SPOT)
							continue;
						
						if ((r == pm_row && c == col) ||
							(r == mm_row && c == col))
							continue;
						
						zero = true; 
						
						if (r <= pm_row + dim)
						{
							mask_pm = LayoutEvaluation.WEIGHT_DIST
										[dim + r - pm_row][dim + c - col];
						
							unmask_pm = LayoutEvaluation.WEIGHT_DIST
										[dim + pm_row - r][dim + col - c];
							
							if (mask_pm > 0 || unmask_pm > 0)
								zero = false;
						}
						else
							mask_pm = unmask_pm = 0;

						if (r >= mm_row - dim)
						{
							mask_mm = LayoutEvaluation.WEIGHT_DIST
										[dim + r - mm_row][dim + c - col];

							unmask_mm = LayoutEvaluation.WEIGHT_DIST
										[dim + mm_row - r][dim + col - c];
							
							if (mask_mm > 0 || unmask_mm > 0)
								zero = false;
						}
						else
							mask_mm = unmask_mm = 0;

						if (!zero)
							addSingleProbeCost (id, mask_pm, unmask_pm,
														mask_mm, unmask_mm);
					}
			}
			
			private void addSingleProbeCost (int id, double mask_pm,
					double unmask_pm, double mask_mm, double unmask_mm)
			{
				int		b, pos, word, bitmask = 0;
				double	m;
				
				for (b = 0, word = -1, pos = 0; pos < embed_len; pos++)
				{
					if (pos % Integer.SIZE == 0)
					{
						word++;
						bitmask = 0x01 << (Integer.SIZE-1);
					}
					else
						bitmask >>>= 1;
					
					if ((chip.embed[id][word] & bitmask) != 0)
					{
						mask_cost_pm[pos] += mask_pm;
						mask_cost_mm[pos] += mask_mm;
						b++;
					}
					else
					{
						m = LayoutEvaluation.positionMultiplier (probe_len, b);
						unmask_cost_pm[pos] += unmask_pm * m;
						unmask_cost_mm[pos] += unmask_mm * m;							
					}
				}
			}
			
			@Override
			protected double computeMatrix (double matrix[][], char probe[],
					int mid_pm, int mid_mm)
			{
				double	mask, unmask;
				int 	r, c, start;
				
				matrix[0][0] = 0;
				for (c = (start = 1); c <= embed_len; c++)
					matrix[0][c] = matrix[0][c - 1] + pos_mult[0] *
									(mask_cost_pm[c -1] + mask_cost_mm[c -1]);
				
				for (r = 1; r <= probe_len + 1; r++, start++)
				{
					matrix[r][start - 1] = Double.POSITIVE_INFINITY;
					
					for (c = start; c <= embed_len; c++)
					{
						mask = mask_cost_mm[c - 1] + mask_cost_pm[c - 1];
						mask = matrix[r][c - 1] + pos_mult[r] * mask; 
						
						if (probe[r - 1] != chip.dep_seq[c -1])
						{
							unmask = Double.POSITIVE_INFINITY;
						}
						else
						{
							unmask = matrix[r - 1][c - 1];
							
							if (r == mid_pm)
								unmask += pos_mult[mid_mm] * mask_cost_mm[c -1];
							else
								unmask += unmask_cost_mm[c - 1];

							if (r == mid_mm)
								unmask += pos_mult[mid_pm] * mask_cost_pm[c -1];
							else
								unmask += unmask_cost_pm[c - 1];
						}
						
						matrix[r][c] = Math.min(mask, unmask);
						
						if (Double.isInfinite(matrix[r][c])) start++;
					}
				}
				
				return matrix[probe_len + 1][embed_len];
			}
			
			@Override
			protected void encodeEmbedding (int id_1, int id_2,
					double matrix[][])
			{
				int		r, c, pos, word, bitmask;
				double	mask_cost;
				
				pos = embed_len;
				word = pos / Integer.SIZE;
				bitmask = 0x01 << (Integer.SIZE - 1 - (pos % Integer.SIZE));
				
				for (r = probe_len + 1, c = embed_len; pos > 0; c--)
				{
					if ((pos-- % Integer.SIZE) == 0)
					{
						word--;
						bitmask = 0x01;
					}
					else
						bitmask <<= 1;

					chip.embed[id_1][word] &= ~bitmask;
					chip.embed[id_2][word] &= ~bitmask;
					
					if (r == 0) continue;
					
					mask_cost = pos_mult[r] * (mask_cost_mm[c - 1] +
													mask_cost_pm[c - 1]);
					
					if (matrix[r][c] == matrix[r][c - 1] + mask_cost)
						continue;
					
					if (r != AffymetrixChip.AFFY_MIDDLE_BASE)
						chip.embed[id_2][word] |= bitmask;
					
					if (r != AffymetrixChip.AFFY_MIDDLE_BASE + 1)
						chip.embed[id_1][word] |= bitmask;
					
					r--;
				}
			}
		}
	}
}
