/*
 * Chip.java
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

import java.util.*;
import java.io.*;

/**
 * This class contains the layout specification of a (high-density) microarray
 * chip. It consists mainly of a rectangular grid of sites (or <EM>spots</EM>)
 * and a set of probes that are mapped to sites on a one-to-one basis.
 *
 * <P>Probes can appear in single instances, pairs or, in the general case, in
 * tuples. The exact scheme, however, is not determined by this class. Instead,
 * this must be specified by sub-classes (see known sub-classes below for
 * available options).</P>
 *
 * <P>Probes are typically synthesized in parallel, on the chip, in a series of
 * repetitive steps. Each step appends the same nucleotide to probes positioned
 * in selected regions of the chip. Selection occurs by exposure to light with
 * the help of a photolithographic mask that allows or obstructs the passage of
 * light accordingly. Thus, each masking step induces the addition of a
 * particular nucleotide to a group of probes. The particular sequence of steps
 * used to produce the chip is called the <EM>deposition sequence</EM>.</P>
 *
 * <P>Obviously, the deposition sequence must be a supersequence of all probe
 * sequences on the chip. In general, a probe can be embedded within the
 * deposition sequence in several ways. An embedding can be seen as a binary
 * string with zeros corresponding to masked and ones corresponding to unmasked
 * steps. In order to save space, probes sequences are stored as a binary
 * sequence according to one particular embedding within the deposition
 * sequence. This embedding can be modified at any time, provided that the
 * resulting probe sequence remains the same.</P>
 *
 * <P>The main purpose of this class is to hold the full specification of a
 * microarray chip that can be used by a {@linkplain PlacementAlgorithm} to
 * produce a new layout containing the same set of probes.</P>
 *
 * @author Sergio A. de Carvalho Jr.
 * @see PlacementAlgorithm
 */
public abstract class Chip implements Cloneable
{
	/**
	 * Total number of rows of sites on the chip.
	 * The rows will be numbered from 0 to num_rows - 1.
	 */
	protected int num_rows;

	/**
	 * Total number of columns of sites on the chip.
	 * The columns will be numbered from 0 to num_cols - 1.
	 */
	protected int num_cols;

	/**
	 * The total number of probe instances (single probes, pairs or tuples) the
	 * chip contains.
	 */
	protected int num_probes;

	/**
	 * The embedding length of any probe on the chip. This is precisely the
	 * length of the deposition sequence used to synthesize the probes.
	 */
	protected int embed_len;

	/**
	 * The length of the probes (number of bases) on this chip.
	 */
	protected int probe_len;

	/**
	 * The deposition sequence, i.e. the actual sequence of bases used to
	 * synthesize the probes on the chip. <B>Implementation note:</B> for
	 * ease of access, this is a public variable.
	 */
	public char dep_seq[];

	/**
	 * The length of the basic permutation that forms the deposition sequence.
	 */
	public int dep_seq_cycle;

	/**
	 * The matrix of sites (or spots) on the chip. The matrix's first dimension
	 * corresponds to the spot's row number while the second dimension is the
	 * column number (<CODE>spot[row][col]</CODE>). Each spot can either be
	 * uninitialized ({@link #UNINITIALIZED_SPOT}), empty ({@link #EMPTY_SPOT})
	 * or associated with a unique probe ID (see {@link #embed}).
	 * <B>Implementation note:</B> for ease of access, this is a public
	 * variable.
	 */
	public int spot[][];

	/**
	 * Constant that flags that a certain spot is empty (does not contain a
	 * probe).
	 */
	public static final int EMPTY_SPOT = -1;

	/**
	 * Constant that flags that a certain spot has not been initialized (has
	 * been neither associated with a probe, nor marked as empty).
	 */
	protected static final int UNINITIALIZED_SPOT = -2;

	/**
	 * This is an array containing the embeddings of all probes on the chip
	 * (one embedding for every single probe of a probe pair or tuple). Each
	 * probe has a unique probe ID which is used to index this array. To reduce
	 * memory requirements, all embeddings are encoded as binary strings (whose
	 * length equals {@link #embed_len}) that are stored in an array of
	 * integers. This binary string tells how the probe is embedded into the
	 * deposition sequence. A <CODE>'1'</CODE> means that the base at the
	 * corresponding position of the deposition sequence is appended to the
	 * probe during the synthesis process (the spot is unmasked), while a
	 * <CODE>'0'</CODE> means that the probe does not receive that base (the
	 * spot is masked). For example, if the deposition sequence is
	 * <CODE>'TGCATGCATGCA...'</CODE>, the embedding
	 * <CODE>'0010111001001'</CODE> corresponds to the probe
	 * <CODE>'CTGTA...'</CODE>. <B>Implementation notes:</B> 1) for ease of
	 * access, this is a public variable. <B>To do:</B> use an array of bites
	 * instead of an array of integers to reduce wasted memory space.
	 */
	public int embed[][];

	/**
	 * Flags which spots are fixed. Fixed spots cannot be changed, i.e. if they
	 * are empty, they must remain empty; if they have a probe, this probe
	 * cannot be moved anywhere else. <B>Implementation note:</B> this is
	 * implemented as a BitSet to reduce memory usage.
	 */
	protected BitSet fixed_spots;

	/**
	 * Indicates that the chip layout specification (list of probes, fixed
	 * spots, etc.) has already been loaded.
	 */
	protected boolean input_done;

	/**
	 * A list of fixed probe IDs. Fixed probes are those that belong to fixed
	 * spots and thus cannot be relocated by a {@linkplain PlacementAlgorithm}.
	 * This list is assembled by the {@link #readLayout} abstract method. Its
	 * contents, therefore, depends on the specific implementations provided by
	 * the subclasses, but it is generally encouraged that this list should
	 * contain only the main probe ID in case of a multi-probe type of chip
	 * (such as {@linkplain AffymetrixChip}). This list is mainly used by the
	 * {link getMovableProbes} method to generate a list of movable probe IDs
	 * (any probe ID not found on this list is considered <EM>movable</EM>.
	 */
	protected int fixed_probe[];

	/**
	 * Creates a new instance of a chip.
	 *
	 * @param num_rows number of rows of sites
	 * @param num_cols number of columns of sites
	 * @param num_probes total number of probes (not probe pairs or tuples)
	 * @param probe_len length of the probes
	 * @param dep_seq deposition sequence
	 */
	public Chip (int num_rows, int num_cols, int num_probes, int probe_len,
		String dep_seq)
	{
		// validate parameters
		if (num_rows < 1 || num_cols < 1)
			throw new IllegalArgumentException
				("Invalid chip dimensions.");

		if (num_probes < 1 || num_probes > num_rows * num_cols)
			throw new IllegalArgumentException
				("Invalid number of probes.");

		if (probe_len < 1)
			throw new IllegalArgumentException
				("Invalid probe length.");
		
		checkDepositionSequence (dep_seq);
		
		// store chip parameters
		this.num_rows = num_rows;
		this.num_cols = num_cols;
		this.num_probes = num_probes;
		this.embed_len = dep_seq.length();
		this.probe_len = probe_len;

		// how many integers are needed to encode each embedding?
		int b = (int) Math.ceil((double) embed_len / Integer.SIZE);

		// allocate space for probe embeddings
		this.embed = new int [num_probes][b];

		// allocate space for spots
		this.spot = new int [num_rows][num_cols];

		// create BitSet of flags for fixed spots
		// (initially, all bits in the set are set to false)
		this.fixed_spots = new BitSet (num_rows * num_cols);

		// flag that the layout spec must still be loaded
		input_done = false;
	}
	
	/**
	 * Checks whether the deposition sequence is valid. This method is called
	 * by the constructor and, in the case of an invalid deposition sequence,
	 * throws a IllegalArgumentException which aborts the instantiation of
	 * the Chip object. This method also checks if the deposition sequence is
	 * cyclical (and stores the cycle length in {@link #dep_seq_cycle}.
	 * 
	 * @param seq the deposition sequence to be checked
	 */
	protected void checkDepositionSequence (String seq)
	{
		int		i, c;
		char	base;
		boolean	cyclical;
		
		embed_len = seq.length();
		
		// check minimum sequence length
		if (embed_len <= 0)
			throw new IllegalArgumentException
				("Invalid deposition sequence: invalid length.");
		
		dep_seq = new char [embed_len];

		// check if fully determined and composed of valid bases
		for (i = 0; i < embed_len; i++)
		{
			base = seq.charAt(i);
			
			if (base == 'A' || base == 'C' || base == 'G' || base == 'T')
				dep_seq[i] = base;
			else
				throw new IllegalArgumentException
					("Invalid deposition sequence: invalid base at position " +
						(i + 1) + ".");
		}
		
		// check if dep. seq. is cyclical
		// max cycle length is half the length of dep. seq.
		c = 1;
		while (c <= embed_len / 2)
		{
			cyclical = true;
			
			for (i = 0; cyclical && i + c < embed_len; i++)
				if (dep_seq[i] != dep_seq[i + c])
					cyclical = false;

			if (cyclical)
			{
				// save cylce length and return
				dep_seq_cycle = c;
				return;
			}
			
			c++;
		}
		
		// zero means dep. seq. is not cyclical
		dep_seq_cycle = 0;
	}

	/**
	 * Returns the lenght of basic permutation that forms the deposition
	 * sequence or zero if the deposition sequence is not cyclical.
	 *
	 * @return deposition sequence's cycle lenght or zero if non-cyclical
	 */
	protected int depositionSequenceCycleLength ()
	{
		return dep_seq_cycle;
	}

	/**
	 * Returns the number of rows of sites this chip has.
	 *
	 * @return number of rows of sites on the chip
	 */
	public int getNumberOfRows ()
	{
		return this.num_rows;
	}

	/**
	 * Returns the number of columns of sites this chip has.
	 *
	 * @return number of columns of sites on the chip
	 */
	public int getNumberOfColumns ()
	{
		return this.num_cols;
	}

	/**
	 * Returns the number of probes this chip has.
	 *
	 * @return number of probes on the chip
	 */
	public int getNumberOfProbes ()
	{
		return this.num_probes;
	}

	/**
	 * Returns the number of internal borders this chip has. An internal border
	 * separates two spots of the chip. Since the chip is a rectangular grid,
	 * the number of internal borders equals
	 * <CODE>r * (c - 1) + c * (r - 1)</CODE>, where <CODE>r</CODE> and
	 * <CODE>c</CODE> are the number of rows and columns in the chip,
	 * respectively. 
	 *
	 * @return number of probes on the chip
	 */
	public int getNumberOfBorders ()
	{
		return num_rows * (num_cols - 1) + num_cols * (num_rows - 1);
	}

	/**
	 * Returns the length of the probe embeddings (which equals the length of
	 * the deposition sequence).
	 *
	 * @return length of the probes' embeddings
	 */
	public int getEmbeddingLength ()
	{
		return this.embed_len;
	}

	/**
	 * Returns the length of the probe sequences (number of nucleotides).
	 *
	 * @return length of the probe sequences
	 */
	public int getProbeLength ()
	{
		return this.probe_len;
	}

	/**
	 * Returns a rectangular region representing the grid of spots of this
	 * chip.
	 *
	 * @return rectangular region representing the grid of spots
	 */
	public RectangularRegion getChipRegion ()
	{
		return new RectangularRegion (0, getNumberOfRows() - 1,
										0, getNumberOfColumns() - 1);
	}

	/**
	 * Returns a list of non-fixed probe IDs. These are probes not located on
	 * fixed spots and which, therefore, can be relocated by a
	 * {@linkplain PlacementAlgorithm}. This is an abstract method. Thus, the
	 * exact result depends on the specific implementations provided by the
	 * subclasses, but it is generally encouraged that this list should contain
	 * only the main probe ID in case of a multi-probe type of chip (such as
	 * {@linkplain AffymetrixChip}). This list is assembled from the list of
	 * fixed probes ({@link #fixed_probe}).
	 *
	 * @return a list of movable (not fixed) probe IDs
	 */
	public abstract int[] getMovableProbes ();

	/**
	 * Returns a list of ALL probe IDs (including those of fixed probes and all
	 * probes instances in case of a multi-probe type of chip).
	 *
	 * @return a list of all probe IDs
	 */
	public int[] getAllProbes ()
	{
		int id[] = new int [this.num_probes];
		
		for (int i = 0; i < this.num_probes; i++)
			id[i] = i;
		
		return id;
	}

	/**
	 * Resets the layout of this microarray chip by marking all non-fixed spots
	 * as empty. This allows a new layout to be produced by a
	 * {@linkplain PlacementAlgorithm}.
	 */
	public void resetLayout ()
	{
		for (int r = 0; r < num_rows; r++)
			for (int c = 0; c < num_cols; c++)
				if (!isFixedSpot(r, c))
					spot[r][c] = EMPTY_SPOT;
	}

	/**
	 * Returns true if the spot at position (row, col) is fixed (cannot be
	 * changed), or false otherwise.
	 *
	 * @param row spot's row number
	 * @param col spot's column number
	 * @return boolean indicating whether the given spot is fixed or not
	 */
	public boolean isFixedSpot (int row, int col)
	{
		return this.fixed_spots.get (row * this.num_cols + col);
	}

	/**
	 * Marks the spot at position (row, col) as fixed or non-fixed. This method
	 * can only be called internally.
	 *
	 * @param row spot's row number
	 * @param col spot's column number
	 * @param fixed true is the spot must be fixed, false otherwise
	 */
	protected void setFixedSpot (int row, int col, boolean fixed)
	{
		this.fixed_spots.set (row * this.num_cols + col, fixed);
	}

	/**
	 * Returns a BitSet containing a set of flags indicating which spots are
	 * fixed and which spots are not. Fixed spots cannot be changed. This might
	 * be useful in several cases, e.g. for a {@linkplain PlacementAlgorithm}
	 * to know how many spots are fixed. Spot at position (row, col) has index
	 * (row * num_cols + col) of the BitSet, where num_cols is the number of
	 * columns this chip has. This method returns a clone of the current
	 * BitSet, so changes will not affect the chip's specification.
	 *
	 * @return BitSet a set of bits as a BitSet instance indicating which spots
	 * are fixed
	 */
	public BitSet getFixedSpots ()
	{
		return (BitSet) fixed_spots.clone();
	}

	/**
	 * Read a chip layout specification from a character input stream. This
	 * method calls the other {@link #readLayout(Reader,boolean)} method (that
	 * must be provided by sub-classes) with the default behaviour
	 * (ignore_fixed = false).
	 *
	 * @param input an character input stream
	 * @throws IOException if an I/O error occurrs or input does not comply
	 * with the format rules
	 */
	public void readLayout (Reader input) throws IOException
	{
		readLayout (input, false);
	}

	/**
	 * Read a chip layout specification from a character input stream. This
	 * method must be provided by sub-classes who should specify the input's
	 * format. There is an option on whether fixed spots must have its
	 * contents preserved (ignore_fixed = false) or not (ignore_fixed = true).
	 * In the first case, the chip specification will contain a list of fixed
	 * spots and fixed probes (those found on fixed spots). Otherwise, all
	 * fixed spots will be treated as non-fixed (which may result, for
	 * instance, in their probes being moved to a new location).
	 *
	 * @param input an character input stream
	 * @param ignore_fixed true if fixed status should be ignored, false
	 * otherwise
	 * @throws IOException if an I/O error occurrs or input does not comply
	 * with the format rules
	 */
	public abstract void readLayout (Reader input, boolean ignore_fixed)
		throws IOException;

	/**
	 * Create a random set of probes and a random layout for this chip. This
	 * method is an alternative to reading a layout from an input stream
	 * ({@link #readLayout(Reader)}), and is specially useful for evaluating
	 * algorithms. This method must be provided by sub-classes according to
	 * their specific implementation details.
	 */
	public abstract void createRandomLayout ();
	
	/**
	 * This method stores a probe sequence as a binary string encoded in
	 * integers representing its embedding into the deposition sequnece. The
	 * input consists of the probe's ID and the embedding specified as a string
	 * of letters (representing the nucleotide bases synthesized at a
	 * particular step) and white spaces (representing masked steps).
	 *
	 * @param embedding embedding specified as a string of letters and white
	 * spaces
	 * @param probe_id probe ID
	 */
	protected void encodeEmbedding (String probe, String embedding, int probe_id)
	{
		char ch;
		int  mask = 0, w, pos, len = 0;

		// validate embedding length
		if (embedding.length() != embed_len)
			throw new IllegalArgumentException ("invalid length");

		// turn all bits off
		for (w = 0; w < embed[probe_id].length; w++)
			embed[probe_id][w] = 0;

		for (w = -1, pos = 0; pos < embed_len; pos++)
		{
			if (pos % Integer.SIZE == 0)
			{
				// next 4-byte word
				w++;

				// turn on very first bit of mask only
				mask = 0x01 << (Integer.SIZE - 1);
			}

			// if step is not masked
			if ((ch = embedding.charAt(pos)) != ' ')
			{
				// check that the embedding "agree" with the
				// probe's base at this postition
				if (probe.charAt(len) != ch)
					throw new IllegalArgumentException ("probe sequence and " +
						"embedding differ at step " + pos);
				
				// check that the embedding "agree" with the
				// deposition sequence at this postition
				if (dep_seq[pos] != ch)
					throw new IllegalArgumentException ("base at step " + pos +
						" is not synchronized with the deposition sequence");
				
				// turn on bit to indicate productive step
				embed[probe_id][w] |= mask;
				
				len++;
			}

			// shift bit to the right
			// ('>>>' means unsigned shift)
			mask >>>= 1;
		}
		
		// check probe length
		if (len != probe_len)
			throw new IllegalArgumentException ("unexpected probe length: " +
				len);
	}
	
	/**
	 * Computes a rank for each probe based on the probe's sequence, which can
	 * be used to sort the probes lexicographically.
	 * 
	 * @param pid array of probe IDs
	 * @param start starting position in the array 
	 * @param end last position in the array
	 * @return an array of numbers corresponding to each input probe
	 */
	public long[] computeProbeRanks (int pid[], int start, int end)
	{
		int  i, word, step, bitmask = 0;
		long base_mask, rank[];
		
		rank = new long[end - start + 1];
		
		for (i = start; i <= end; i++)
			rank[i - start] = 0;
		
		for (word = -1, step = 0; step < embed_len; step++)
		{
			if (step % Integer.SIZE == 0)
			{
				bitmask = 0x01 << (Integer.SIZE - 1);
				word++;
			}
			else
				bitmask >>>= 1;

			switch (dep_seq[step])
			{
				case 'A':
					base_mask = 0x00;
					break;
					
				case 'C':
					base_mask = 0x01;
					break;

				case 'G':
					base_mask = 0x02;
					break;
					
				case 'T':
					base_mask = 0x03;
					break;
				
				default:
					throw new IllegalArgumentException
						("Illegal deposition sequence.");
			}
			
			for (i = start; i <= end; i++)
				if ((bitmask & embed[pid[i]][word]) != 0)
				{
					rank[i - start] <<= 2;
					rank[i - start] |= base_mask;
				}
		}
		
		return rank;
	}

	/**
	 * Print the specification of the chip's current layout. This method should
	 * be provided by the sub-classes depending on the specific probe scheme
	 * used (single probes, probe pairs, etc.).
	 *
	 * @param out a PrintWriter stream
	 * @throws IOException if an error occurs while writing on the stream
	 */
	public abstract void writeLayout (PrintWriter out) throws IOException;

	/**
	 * Prints the stored embedding of a probe on the standard error output
	 * stream. This method serves for debugging purposes only.
	 *
	 * @param probe_id probe ID
	 */
	public void printEmbedding (int probe_id)
	{
		int mask = 0, w, pos;
		
		for (w = -1, pos = 0; pos < embed_len; pos++)
		{
			if (pos % Integer.SIZE == 0)
			{
				// use next 4-byte word
				w++;

				// turn on very first bit of mask only
				mask = 0x01 << (Integer.SIZE - 1);
			}

			if ((mask & embed[probe_id][w]) == 0)
				System.err.print (" ");
			else
				System.err.print (dep_seq[pos]);

			mask >>>= 1;
		}
	}

	/**
	 * Prints the stored embedding of a probe as a binary string on the standard
	 * error output stream. This method serves for debugging purposes only.
	 *
	 * @param probe_id probe ID
	 */
	public void printBinaryEmbedding (int probe_id)
	{
		int mask = 0, w, pos;
		
		for (w = -1, pos = 0; pos < embed_len; pos++)
		{
			if (pos % Integer.SIZE == 0)
			{
				// use next 4-byte word
				w++;

				// turn on very first bit of mask only
				mask = 0x01 << (Integer.SIZE - 1);
			}

			if ((mask & embed[probe_id][w]) == 0)
				System.err.print ('0');
			else
				System.err.print ('1');

			mask >>>= 1;
		}
	}

	/**
	 * Prints the stored probe sequence on the standard error output stream.
	 * This method serves for debugging purposes only.
	 *
	 * @param probe_id probe ID
	 */
	public void printProbe (int probe_id)
	{
		int mask = 0, w, pos;
		
		for (w = -1, pos = 0; pos < embed_len; pos++)
		{
			if (pos % Integer.SIZE == 0)
			{
				// use next 4-byte word
				w++;

				// turn on very first bit of mask only
				mask = 0x01 << (Integer.SIZE - 1);
			}
			else
				mask >>>= 1;

			if ((mask & embed[probe_id][w]) != 0)
				System.err.print (dep_seq[pos]);
		}
	}

	/**
	 * Creates and returns a copy of this Chip object. The new object will
	 * contain the same chip specification as the cloned object, i.e. same
	 * chip dimension, same probes on the same locations and with the same
	 * embeddings, same deposition sequence and so on. Specific
	 * sub-classes may need to extend this method to account for added
	 * member variables.
	 *
	 * @return a clone of this instance
	 */
	@Override
	protected Chip clone ()
	{
		Chip c;
		int i;
		
		try
		{
			// clone the chip
			c = (Chip) super.clone();
		}
		catch (CloneNotSupportedException e)
		{
			// this shouldn't happen anyway...
			throw new AssertionError();
		}
		
		// now we need to clone members that are
		// not of primitive types (including arrays)

		c.dep_seq = this.dep_seq.clone();
		c.fixed_probe = this.fixed_probe.clone();
		c.fixed_spots = (BitSet) this.fixed_spots.clone();

		// special care is needed with multidimensional arrays since
		// their clones are shallow (subarrays are still shared)
		// see the Java Language Specification's chapter on arrays:
		// java.sun.com/docs/books/jls/second_edition/html/arrays.doc.html

		// spots
		c.spot = this.spot.clone();
		for (i = 0; i < this.spot.length; i++)
			c.spot[i] = this.spot[i].clone();
		
		// embeddings
		c.embed = this.embed.clone();
		for (i = 0; i < this.embed.length; i++)
			c.embed[i] = this.embed[i].clone();
		
		return c;
	}

	/**
	 * Indicates whether some other Chip instance is "equal to" this one.
	 * Equals here means that both chips have the same specification, i.e. both
	 * have same dimension, same probes on the same locations and with the same
	 * embeddings, same deposition sequence and so on. Specific sub-classes may
	 * need to extend this method to account for added member variables. The
	 * {link #compatible} method provides another way of comparing two chip
	 * instances.
	 *
	 * @param obj the reference object with which to compare
	 * @return true if this Chip is "equal to" the argument; false otherwise
	 */
	@Override
	public boolean equals (Object obj)
	{
		Chip other;
		int i, j;
		
		// equal if references point
		// to the same instance
		if (this == obj) return true;
		
		// cannot compare objects that are not Chip instances
		if (!(obj instanceof Chip)) return false;
		
		other = (Chip) obj;
		
		// check basic properties
		if (this.num_rows != other.num_rows) return false;
		if (this.num_cols != other.num_cols) return false;
		if (this.num_probes != other.num_probes) return false;
		if (this.probe_len != other.probe_len) return false;
		
		// check if have the same deposition sequence
		if (this.embed_len != other.embed_len) return false;
		
		for (i = 0; i < dep_seq.length; i++)
			if (this.dep_seq[i] != other.dep_seq[i])
				return false;
		
		// check if specification (probes, spots, etc.)
		// of both chips have been loaded (or not)
		if (this.input_done != other.input_done) return false;
		
		// if specifications have not been loaded yet,
		// then the chips are considered equal
		if (!input_done) return true;
		
		// check if spots have the same contents
		if (this.spot.length != other.spot.length) return false;
		
		for (i = 0; i < this.spot.length; i++)
		{
			if (this.spot[i].length != other.spot[i].length) return false;
			
			for (j = 0; j < this.spot[i].length; j++)
				if (this.spot[i][j] != other.spot[i][j])
					return false;
		}
							
		// check if embeddings have the same contents
		if (this.embed.length != other.embed.length) return false;
		
		for (i = 0; i < this.embed.length; i++)
		{
			if (this.embed[i].length != other.embed[i].length) return false;
			
			for (j = 0; j < this.embed[i].length; j++)
				if (this.embed[i][j] != other.embed[i][j])
					return false;
		}

		// check the list of fixed probes
		if (this.fixed_probe.length != other.fixed_probe.length) return false;
		
		for (i = 0; i < fixed_probe.length; i++)
			if (this.fixed_probe[i] != other.fixed_probe[i])
				return false;
		
		// check the list of fixed spots
		if (this.fixed_spots == null)
		{
			if(!(other.fixed_spots == null)) return false;
		}
		else
		{
			if (!(this.fixed_spots.equals(other.fixed_spots))) return false;
		}
		
		// if passed all tests, they are considered equal
		return true;
	}	

	/**
	 * Indicates whether some other Chip instance is "compatible with" this
	 * one. Compatible here means that both chips have almost the same
	 * specification, i.e. both have same dimension, same probes, same
	 * deposition sequence and so on. However, two differences are allowed:
	 * 1) the probes might be located on different spots; and 2) the same
	 * probe on one chip might have a different but "compatible" embeddings
	 * on the other chip). Compatible embeddings means that they produce the
	 * same probe sequence under their own deposition sequences (which must
	 * not necessarily be the same for both chips); see
	 * {link #compatibleEmbedding}.
	 *
	 * <P>In normal circunstances, compatible chips can only be produced by
	 * cloning a chip (and optionally changing its layout). This is mainly
	 * because the compatibility test assumes that a given probe has the same
	 * ID on both chips. Thus, it is possible that two chips containing
	 * identical set of probes may be regarded as not compatible just because
	 * the lists of probes were input in a different order.</P>
	 *
	 * <P>The {link #validateLayout} method is called on both chips before
	 * checking their compatibility to make sure that both are in a correct
	 * state (and thus avoiding errors such as NullPointerException).</P>
	 *
	 * <P>Sub-classes may need to extend this method to account for their own
	 * member variables and specific requirements. Note that this method is
	 * similar but conceptually different from the {link #equals} method.</P>
	 *
	 * @param other the reference object with which to compare
	 * @return true if this Chip is "compatible with" the argument; false
	 * otherwise
	 */
	public boolean compatible (Chip other)
	{
		int	i, j;
		
		// first check that both chips have valid layouts;
		if (!this.validateLayout()) return false;
		if (!other.validateLayout()) return false;
		
		// if references point to the same instance
		// there is no need to check for compatibility
		if (this == other) return true;
		
		// check basic properties
		if (this.num_rows != other.num_rows) return false;
		if (this.num_cols != other.num_cols) return false;
		if (this.num_probes != other.num_probes) return false;
		if (this.probe_len != other.probe_len) return false;
		
		// check if specification (probes, spots, etc.)
		// of both chips have been loaded (or not)
		if (this.input_done != other.input_done) return false;
		
		// if specifications have not been loaded yet,
		// then the chips are considered compatible
		if (!input_done) return true;
		
		// the list of fixed probes must be equal
		if (this.fixed_probe.length != other.fixed_probe.length) return false;
		for (i = 0; i < fixed_probe.length; i++)
			if (this.fixed_probe[i] != other.fixed_probe[i])
				return false;
		
		// the list of fixed spots must also be equal
		if (!(this.fixed_spots.equals(other.fixed_spots))) return false;
		
		// check if embeddings are compatible
		for (i = 0; i < num_probes; i++)
			if (!compatibleEmbedding(i, other))
				return false;

		// check that the chips have the same set of fixed spots
		for (i = 0; i < num_rows; i++)
			for (j = 0; j < num_cols; j++)
			{
				if (this.isFixedSpot (i, j))
				{
					if (!other.isFixedSpot(i, j)) return false;
					
					// also check that fixed spots have the same contents
					if (this.spot[i][j] != other.spot[i][j]) return false;
				}
				else
				{
					if (other.isFixedSpot(i, j)) return false;
				}
			}
		
		// if passed all tests, they are considered compatible
		return true;
	}

	/**
	 * Indicates whether a given probe has a "compatible" embedding on other
	 * Chip instance. Compatible here means that the embeddings produce the
	 * same probe sequence under their own deposition sequence (which must not
	 * necessarily be the same for both chips). This method is mainly used by
	 * the {link #compatible} method.
	 *
	 * @param id the ID of the probe whose embeddings are to be verified
	 * @param other the other chip where the embedding is to be checked
	 * @return true if the embeddings are "compatible"; false otherwise
	 */
	protected boolean compatibleEmbedding (int id, Chip other)
	{
		int 	this_pos, other_pos, this_w, other_w, this_mask, other_mask;
		char	ch;
		
		this_pos = other_pos = -1;
		this_w = other_w = -1;
		this_mask = other_mask = 0;
		
		while (++this_pos < this.embed_len)
		{
			if (this_pos % Integer.SIZE == 0)
			{
				// use next 4-byte word
				this_w++;
				this_mask = 0x01 << (Integer.SIZE - 1);
			}
			else
			{
				this_mask >>>= 1;
			}

			// when the embedding has a set bit, that
			// means we found the probe's next base
			if ((this_mask & this.embed[id][this_w]) != 0)
			{
				ch = this.dep_seq[this_pos];
				
				// now we need to find the next
				// base on the 'other' chip
				while (++other_pos < other.embed_len)
				{
					if (other_pos % Integer.SIZE == 0)
					{
						// use next 4-byte word
						other_w++;
						other_mask = 0x01 << (Integer.SIZE - 1);
					}
					else
					{
						other_mask >>>= 1;
					}

					if ((other_mask & other.embed[id][other_w]) != 0)
					{
						// check if bases are equal
						if (other.dep_seq[other_pos] != ch)
							// if not, the embeddings produce
							// different probe sequences
							return false;
							
						break;
					}
				}
			}
		}
		
		// if reached here,
		// embeddings are "compatible"
		return true;
	}

	/**
	 * Checks that the layout specification is valid. This method is useful for
	 * checking that a {@linkplain PlacementAlgorithm} has produced a new but
	 * still valid layout (according to the original specification). For
	 * instance, it checks whether the fixed spots were respected (their
	 * contents remain unchanged).
	 * 
	 * <P>Sub-classes may need to extend this method to account for their own
	 * member variables and specific requirements.</P>
	 *
	 * <P>Note that this method is called by the {link #compatible} method to
	 * ensure that both chips being compared have a valid layout. Nonetheless,
	 * this method can be called independently to check that the chip is in a
	 * valid state after, for instance, a new layout is created.<P>
	 * 
	 * @return true if the current layout specification is valid; false
	 * otherwise
	 */
	public abstract boolean validateLayout ();
}
