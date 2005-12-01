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
public abstract class Chip
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
	 * The matrix of sites (or spots) on the chip. Each spot can either be
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
	 * @param embed_len length of the embeddings (deposition sequence)
	 */
	public Chip (int num_rows, int num_cols, int num_probes, int probe_len,
		int embed_len)
	{
		// store chip parameters
		this.num_rows = num_rows;
		this.num_cols = num_cols;
		this.num_probes = num_probes;
		this.embed_len = embed_len;
		this.probe_len = probe_len;

		// how many integers are needed to encode each embedding?
		int b = (int) Math.ceil((double) embed_len / Integer.SIZE);

		// allocate space for probe embeddings
		this.embed = new int [num_probes][b];

		// allocate space for the deposition sequence
		this.dep_seq = new char [embed_len];

		// initialize the deposition sequence as a string of white spaces
		for (b = 0; b < embed_len; b++)
			dep_seq[b] = ' ';

		// allocate space for spots
		this.spot = new int [num_rows][num_cols];

		// create BitSet of flags for fixed spots
		// (initially, all bits in the set are set to false)
		this.fixed_spots = new BitSet (num_rows * num_cols);

		// flag that the layout spec must still be loaded
		input_done = false;
	}

	/**
	 * Creates a new instance of a chip with a pre-defined deposition sequence.
	 *
	 * @param num_rows number of rows of sites
	 * @param num_cols number of columns of sites
	 * @param num_probes total number of probes (not probe pairs or tuples)
	 * @param probe_len length of the probes
	 * @param dep_seq deposition sequence as a String
	 */
	public Chip (int num_rows, int num_cols, int num_probes, int probe_len,
		String dep_seq)
	{
		this (num_rows, num_cols, num_probes, probe_len, dep_seq.length());

		// set deposition sequence
		for (int i = 0; i < dep_seq.length(); i++)
			this.dep_seq[i] = dep_seq.charAt(i);
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
		return new RectangularRegion (0, getNumberOfRows() - 1, 0, getNumberOfColumns() - 1);
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
	 * Marks the spot at position (row, col) as fixed or non-fixed.
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
	 * method must be provided by sub-classes who should specify the input's
	 * format.
	 *
	 * @param input an character input stream
	 * @throws IOException if an I/O error occurrs or input does not comply
	 * with the format rules
	 */
	public abstract void readLayout (Reader input) throws IOException;

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
	protected void encodeEmbedding (String embedding, int probe_id)
	{
		char ch;
		int  mask = 0, w, pos;

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
				if (dep_seq[pos] == ' ')
				{
					// save the base at the corresponding position in
					// the deposition sequence if it is still unknown
					dep_seq[pos] = ch;
				}
				else if (dep_seq[pos] != ch)
				{
					// the embedding does not "agree" with the
					// deposition sequence at this postition
					throw new IllegalArgumentException ("unexpected base at step " + pos);
				}

				// turn on bit to indicate productive step
				embed[probe_id][w] |= mask;
			}

			// shift bit to the right
			// ('>>>' means unsigned shift)
			mask >>>= 1;
		}
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
	 * Prints the stored embedding of a probe on the standard output. This
	 * method serves for debugging purposes only.
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
				System.out.print (" ");
			else
				System.out.print (dep_seq[pos]);

			mask >>>= 1;
		}
	}
}
