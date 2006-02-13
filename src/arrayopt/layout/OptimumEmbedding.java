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
 * probe. It can work with two types of quality measures: border length
 * ({@link #MODE_BORDER_LENGTH}) and conflict index
 * ({@link #MODE_CONFLICT_INDEX}).
 * 
 * <P>This class provides two main functionalities. Firstly, it can compute
 * the minimum distance that any embedding of a probe can have to the current
 * embedding of another probe. (The distance depends on the chosen quality
 * measure and whether the location of the probes on the chip is taken into
 * account.) This is accomplished by the {@link minDistanceProbe} and
 * {@link minDistanceSpot} methods (the first only considers the embeddings,
 * while the latter also includes the location of the probes on the chip).
 * 
 * <P>The other main feature of this class is provided by the
 * {@link reembedProbe} and {@link reembedSpot} methods, which re-embeds a
 * probe with the minimum distance.</P>
 *
 * <P>In fact, some of these methods have several signatures that allow the
 * computed distance to take into account not only a single probe but a set of
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
 * chip instance and the desired mode of operation (either
 * {@link #MODE_BORDER_LENGTH} or {@link #MODE_CONFLICT_INDEX}) passed as
 * arguments. The returned object will provide the functionality specified
 * here.</P>

 * <P>This class is intended for use by other classes inside the layout package
 * and, therefore, has no public methods.</P>
 *
 * @see OptimumEmbedding.Simple.BorderLength
 * @see OptimumEmbedding.Simple.ConflictIndex
 * @see OptimumEmbedding.Affymetrix.BorderLength
 * @see OptimumEmbedding.Affymetrix.ConflictIndex
 */
public abstract class OptimumEmbedding
{
	// TODO Solve the problem with rounding of floating point numbers (when
	// comparing doubles in the encodeEmbedding methods)
	
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
	 */
	static final int MODE_BORDER_LENGTH = 0;
	
	/**
	 * Constant that indicates that the conflict index should be used when
	 * computing the distance between embeddings/spots. It is used to create a
	 * new instance of this class with the {@link #createEmbedder(Chip, int)}
	 * method. 
	 */
	static final int MODE_CONFLICT_INDEX = 1;
	
	/**
	 * This constructor can only be accessed internally. For instantiating
	 * objects of this class, please use {@link createEmbedder}.  
	 * @param chip a chip instance
	 * @param probe_len the length of the probes 
	 */
	protected OptimumEmbedding (Chip chip)
	{
		this.probe_len = chip.getProbeLength();
		this.embed_len = chip.getEmbeddingLength();
		this.chip_region = chip.getChipRegion();
	}
	
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
				return new Affymetrix.BorderLength ((AffymetrixChip) c);
			
			if (mode == MODE_CONFLICT_INDEX)
				return new Affymetrix.ConflictIndex ((AffymetrixChip) c);
			
			throw new IllegalArgumentException
				("Illegal value for argument 'mode'.");
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
	 * @param id probe ID
	 * @param row spot's row coordinate
	 * @param col spot's column coordinate
	 * @return the minimum conflict that the probe spot can cause when placed
	 * on the spot
	 */
	double minDistanceSpot (int id, int row, int col)
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
	double reembedSpot (int id, int row, int col)
	{
		resetCosts();
		addSpotCost (row, col);
		return reembedOptimally (id);
	}

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

	protected static abstract class Simple extends OptimumEmbedding
	{
		protected SimpleChip chip;

		protected double matrix[][];
		
		protected char probe[];

		protected double mask_cost[];
		
		protected double unmask_cost[];
		
		protected Simple (SimpleChip chip)
		{
			super (chip);

			this.chip = chip;
			this.matrix = new double [probe_len + 1][embed_len + 1];
			
			this.probe = new char [probe_len];
			
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
		
		protected abstract void encodeEmbedding (int id);
		
		protected abstract double computeMatrix ();

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
				addProbeCost (id, 1, 1);
			}
		
			@Override
			protected void addProbeCost (int id[], int start, int end)
			{
				for (int p = start; p <= end; p++)
					addProbeCost (id[p], 1, 1);
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
							addProbeCost (id, mask_w, unmask_w);
					}
			}
			
			private void addProbeCost (int id, double mask_weight,
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

	protected static abstract class Affymetrix extends OptimumEmbedding
	{
		protected AffymetrixChip chip;

		protected double matrix_1[][];
		
		protected double matrix_2[][];

		protected char probe_1[];

		protected char probe_2[];
		
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
		}
		
		@Override
		protected int getProbeID (int row, int col)
		{
			return chip.spot[row][col];
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

		private static class BorderLength extends Affymetrix
		{
			private double mask_cost[];
			
			private double unmask_cost[];
			
			protected BorderLength (AffymetrixChip chip)
			{
				super (chip);
				
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
				
				if (!chip.isPMProbe(id)) id = id -1;
				
				decodeEmbedding (id);
				d1 = computeMatrix (matrix_1, probe_1);
				d2 = computeMatrix (matrix_2, probe_2);
				
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
						mask_cost[pos] += 1;
					}
					else
						unmask_cost[pos] += 1;
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
							mask_cost[pos] += 1;
						}
						else
							unmask_cost[pos] += 1;
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
						addSingleProbeCost(id);
				
				// bottom
				r = mm_row + 1;
				if (r <= chip_region.last_row)
					if ((id = chip.spot[r][col]) != Chip.EMPTY_SPOT)
						addSingleProbeCost(id);

				// left
				c = col - 1;
				if (c >= chip_region.first_col)
				{
					// PM row
					if ((id = chip.spot[pm_row][c]) != Chip.EMPTY_SPOT)
						addSingleProbeCost(id);
					// MM row
					if ((id = chip.spot[mm_row][c]) != Chip.EMPTY_SPOT)
						addSingleProbeCost(id);
				}

				// right
				c = col + 1;
				if (c <= chip_region.last_col)
				{
					// PM row
					if ((id = chip.spot[pm_row][c]) != Chip.EMPTY_SPOT)
						addSingleProbeCost(id);
					if ((id = chip.spot[mm_row][c]) != Chip.EMPTY_SPOT)
						addSingleProbeCost(id);
				}
			}

			private void addSingleProbeCost (int id)
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

			private double computeMatrix (double matrix[][], char probe[])
			{
				double	mask, unmask;
				int 	r, c, start;
				
				matrix[0][0] = 0;
				for (c = (start = 1); c <= embed_len; c++)
					matrix[0][c] = matrix[0][c - 1] + mask_cost[c -1];
				
				for (r = 1; r <= probe_len + 1; r++, start++)
				{
					matrix[r][start - 1] = Double.POSITIVE_INFINITY;
					
					for (c = start; c <= embed_len; c++)
					{
						mask = matrix[r][c - 1] + mask_cost[c -1];
						
						if (probe[r - 1] == chip.dep_seq[c -1])
							unmask = matrix[r - 1][c - 1] + unmask_cost[c -1];
						else
							unmask = Double.POSITIVE_INFINITY;
						
						matrix[r][c] = Math.min(mask, unmask);
						
						if (Double.isInfinite(matrix[r][c])) start++;
					}
				}
				
				return matrix[probe_len + 1][embed_len];
			}
			
			private void encodeEmbedding (int id_1, int id_2, double matrix[][])
			{
				int r, c, pos, word, bitmask;
				
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
					
					if (r == 0)	continue;

					if (matrix[r][c] == matrix[r][c -1] + mask_cost[c -1])
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
			private double mask_cost_pm[];
			
			private double mask_cost_mm[];
			
			private double unmask_cost_pm[];
			
			private double unmask_cost_mm[];
			
			private double pos_mult[];
			
			protected ConflictIndex (AffymetrixChip chip)
			{
				super (chip);
				
				int b;
				
				this.mask_cost_pm = new double [embed_len];
				this.mask_cost_mm = new double [embed_len];
				this.unmask_cost_pm = new double [embed_len];
				this.unmask_cost_mm = new double [embed_len];
				
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
			protected void resetCosts ()
			{
				for (int i = 0; i < embed_len; i++)
				{
					mask_cost_pm[i] = mask_cost_mm[i] = 0;
					unmask_cost_pm[i] = unmask_cost_mm[i] = 0;
				}
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
			
			private double computeMatrix (double matrix[][], char probe[],
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
			
			private void encodeEmbedding (int id_1, int id_2, double matrix[][])
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
					
					if (Math.abs(matrix[r][c] - matrix[r][c - 1] - mask_cost) < 0.0001d)
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

