/*
 * AffymetrixChip.java
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

import java.util.regex.*;
import java.util.*;
import java.io.*;

/**
 * This class represents Affymetrix chips. The main particularity of this type
 * of chip is that probes appear in pairs: one probe is called the
 * perfect-match probe (PM) while the other is called the mismatch (MM) probe.
 * The MM probe is a copy of the PM probe where the middle base is exchanged
 * according to the DNA complementarity rules (<CODE>A <-> T</CODE>,
 * <CODE>C <-> G</CODE>).
 *
 * <P>Moreover, probe length if fixed to 25 ({@link #AFFY_PROBE_LENGTH}) and
 * there is a default deposition sequence, which is a truncated repetition of
 * <CODE>TGCA</CODE> of length 74 (all known Affymetrix chips use this
 * sequence).</P>
 *
 * @author Sergio A. de Carvalho Jr.
 */
public class AffymetrixChip extends Chip
{
	/**
	 * This constant holds the length of Affymetrix probes (which cannot be
	 * changed).
	 */
	public static final int AFFY_PROBE_LENGTH = 25;

	/**
	 * This constant holds the middle base position of a probe, where the PM
	 * probes differ from the MM probes.
	 */
	public static final int AFFY_MIDDLE_BASE = 12;

	/**
	 * This constant holds the default deposition sequence of Affymetrix chips,
	 * which is a truncated repetition of <CODE>TGCA</CODE> of length 74.
	 */
	public static final String AFFY_DEFAULT_DEPSEQ = "TGCATGCATGCATGCATGCA" +
					"TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG";

	/**
	 * This constant holds the character used to indicate if a probe is a PM
	 * (perfect match) probe (in the input stream).
	 */
	public static final char AFFY_PM_PROBE = 'P';

	/**
	 * This constant holds the character used to indicate if a probe is a MM
	 * (mismatch) probe (in the input stream).
	 */
	public static final char AFFY_MM_PROBE = 'M';

	/**
	 * A BitSet flagging which probes are of tyoe PM (perfect match). If a
	 * probe is not a PM probe, it must be a MM (mismatch) probe.
	 * <B>Implementation note:</B> this is implemented as a BitSet to reduce
	 * memory usage.
	 */
	protected BitSet pm_probe;

	/**
	 * Creates a new AffymetrixChip instance with the default deposition
	 * sequence (see {@link #AFFY_DEFAULT_DEPSEQ}).
	 *
	 * @param num_rows number of rows of sites
	 * @param num_cols number of columns of sites
	 * @param num_probes total number of probes (not probe pairs)
	 */
	public AffymetrixChip (int num_rows, int num_cols, int num_probes)
	{
		super (num_rows, num_cols, num_probes, AFFY_PROBE_LENGTH,
			AFFY_DEFAULT_DEPSEQ);

		// create BitSet of flags for PM probes
		this.pm_probe = new BitSet (num_probes);
	}

	/**
	 * Creates a new AffymetrixChip instance with the specified deposition
	 * sequence.
	 *
	 * @param num_rows number of rows of sites
	 * @param num_cols number of columns of sites
	 * @param num_probes total number of probes (not probe pairs)
	 * @param dep_seq deposition sequence as a String
	 */
	public AffymetrixChip (int num_rows, int num_cols, int num_probes,
		String dep_seq)
	{
		super (num_rows, num_cols, num_probes, AFFY_PROBE_LENGTH,
			dep_seq);

		// create BitSet of flags for PM probes
		this.pm_probe = new BitSet (num_probes);
	}

	/**
	 * Returns true if a the probe with the given ID is a PM (perfect match)
	 * probe, or false if the probe is a MM (mismatch) probe.
	 *
	 * @param probe_id probe ID
	 * @return boolean indicating whether the probe is a PM probe
	 */
	public boolean isPMProbe (int probe_id)
	{
		return this.pm_probe.get (probe_id);
	}

	/**
	 * Marks the probe with a given ID as a PM probe (true) or as a MM probe
	 * (false).
	 *
	 * @param probe_id probe ID
	 * @param flag true if the probe is a PM probe, false otherwise
	 */
	protected void setPMProbe (int probe_id, boolean flag)
	{
		this.pm_probe.set (probe_id, flag);
	}

	/**
	 * Returns the base complement of a given base. This can be used to
	 * identify the middle base of a MM probe based on the middle base of the
	 * corresponding PM probe (or vice-versa). The complementarity rules are
	 * as follows: <CODE>A <-> T</CODE> and <CODE>C <-> G</CODE>.
	 *
	 * @param base a base whose complement is wanted
	 * @return the base complement
	 */
	public char getBaseComplement (char base)
	{
		switch (base)
		{
			case 'A':
				return 'T';

			case 'C':
				return 'G';

			case 'G':
				return 'C';

			case 'T':
				return 'A';

			default:
				throw new IllegalArgumentException ("Invalid base '" + base +
													"'.");
		}
	}

	/**
	 * Read a chip layout specification from a character input stream. The
	 * input must consist of lines of text with the following TAB-delimited
	 * fields describing a single spot of the chip:<BR>
	 * <BR>
	 * 1. Column number (X coordinate)<BR>
	 * 2. Row number (Y coordinate)<BR>
	 * 3. Group<BR>
	 * 4. Fixed: <CODE>Y</CODE> or <CODE>N</CODE><BR>
	 * 5. PM/MM: <CODE>P</CODE> or <CODE>M</CODE><BR>
	 * 6. Probe sequence<BR>
	 * 7. Embedding<BR>
	 *
	 * <P>The first two fields are the row and column coordinates of the spot
	 * on the chip. Every spot can be either empty or associated with a single
	 * probe. If it is empty, only fields 1, 2 and 4 are relevant, while fields
	 * 5, 6 and 7 must contain a dash (<CODE>-</CODE>) and field 3 is ignored.
	 * The 4th field must contain a <CODE>Y</CODE> if the spot is fixed (its
	 * content cannot be changed; if its empty it must remain empty, if it
	 * contains a probe, this probe cannot be moved elsewhere) or
	 * <CODE>N</CODE> otherwise. If the spot contains a probe, a text tag can
	 * appear on the third field to "group" a set of probes under a single name
	 * (this grouping may or may not be explored by a placement algorithm). The
	 * 5th field must either contain a <CODE>P</CODE> for PM (perfect match)
	 * probes or a <CODE>M</CODE> for MM (mismatch) probes.<P>
	 *
	 * <P>The exact probe sequence appears at field 6th as a sequence of the
	 * characters <CODE>A</CODE>, <CODE>C</CODE>, <CODE>G</CODE> and
	 * <CODE>T</CODE>. The last field must contain a valid embedding of the
	 * probe into the deposition sequence formed by inserting a number of
	 * spaces into the probe sequence so that the resulting string is
	 * "synchronized" with the deposition sequence (see
	 * {@linkplain #encodeEmbedding}).</P>
	 *
	 * <P>Every probe in the input receives a unique probe ID. Any coordinate
	 * pair can appear only once in the input (a spot listed twice will throw
	 * an IOException) but spots can be listed in any order. However, probes of
	 * a pair must be listed one after the other (PM, then MM probe) so that
	 * they are assigned consecutive probe IDs. This is important as the list
	 * of "moveable" probes (see {link getMovableProbes}) will contain only the
	 * ID of the PM probes. Moreover, PM and MM probes must be located on
	 * consecutive rows of the same column of the grid.<P>
	 *
	 * <P>A pair of probes must have the same fixed status (both fixed or both
	 * not fixed). Their embeddings must also contains sequences that differ
	 * only in the middle base ({link AFFY_MIDDLE_BASE})</P>.
	 *
	 * @param input a character input stream
	 * @throws IOException if an I/O error occurrs or input is not compliantan
	 */
	public void readLayout (Reader input) throws IOException
	{
		ArrayList<Integer>	fixed_list;
		BufferedReader		in;
		Pattern				parser;
		String				line, field[];
		int					ln = 0, r, c, probe_id = -1, i = 0;
		boolean				empty, fixed;
		char				type, last_type;

		// check if chip spec has already been input
		if (input_done)
			throw new IllegalStateException
				("Layout specification has already been loaded.");

		// mark all spots as unitialized
		for (r = 0; r < num_rows; r++)
			for (c = 0; c < num_cols; c++)
				this.spot[r][c] = UNINITIALIZED_SPOT;

		// create a list of IDs of fixed probes with an initial
		// capacity of about 5% of the number of probes
		fixed_list = new ArrayList<Integer> ((int) (.05 * num_probes));

		// create a buffered reader to read lines
		in = new BufferedReader (input);

		// lines will be parsed into fields...
		parser = Pattern.compile ("\t");

		// ...and stored in this array of strings
		field = new String [7];

		// probe pairs must be listed together,
		// with the PM probe followed by a MM probe
		last_type = AFFY_MM_PROBE;

		while ((line = in.readLine()) != null)
		{
			// line number
			ln++;

			// parse fields
			field = parser.split(line, field.length);

			try
			{
				// spot coordinates
				// field 0: X coordinate -> column
				// field 1: Y coordinate -> row
				c = Integer.parseInt (field[0]);
				r = Integer.parseInt (field[1]);

				// fixed spot?
				if (field[3].equals("Y"))
					fixed = true;
				else if (field[3].equals("N"))
					fixed = false;
				else
					throw new IOException ("Invalid fixed flag at line " +
											ln + ".");
				// empty spot?
				empty = (field[6].equals("-")) ? true : false;

				// probe type (only first character is significant)
				type = field[4].charAt(0);
			}
			catch (ArrayIndexOutOfBoundsException e)
			{
				// invalid file format
				throw new IOException ("Unable to parse input file at line " +
										ln + ".");
			}
			catch (NumberFormatException e)
			{
				// invalid file format
				throw new IOException ("Invalid spot coordinates at line " +
										ln + ".");
			}

			// validate row and column numbers
			if (r < 0 || r >= num_rows || c < 0 || c >= num_cols)
				throw new IOException ("Invalid spot coordinates at line " +
										ln + ".");

			// check for spot conflict
			if (spot[r][c] != UNINITIALIZED_SPOT)
				throw new IOException ("Spot conflict at row " + r +
										", column " + c + ".");

			// mark spot as fixed or non-fixed
			setFixedSpot(r, c, fixed);

			if (empty)
			{
				// mark spot as empty
				spot[r][c] = EMPTY_SPOT;

				// read next line
				continue;
			}

			// new probe found
			probe_id++;

			// check if number of probes has been exceeded
			if (probe_id >= num_probes)
				throw new IOException
					("Found more probes in the input than expected.");

			// place probe on the spot (mark spot as used)
			this.spot[r][c] = probe_id;

			// check if input alternates between PM and MM probes
			if (type == last_type)
				throw new IOException ("Unexpected probe type '" + type +
										"' at line " + ln + ".");
			try
			{
				// encode probe embedding
				encodeEmbedding (field[5], field[6], probe_id);
			}
			catch (IllegalArgumentException e)
			{
				throw new IOException ("Invalid embedding at line " + ln +
										" (" + e.getMessage() + ").");
			}

			// save probe type
			if (type == AFFY_PM_PROBE)
			{
				setPMProbe (probe_id, true);

				if (fixed)
					// add probe ID to the list of fixed probes
					fixed_list.add(probe_id);
			}
			else if (type == AFFY_MM_PROBE)
			{
				// check if PM and MM embeddings are compatible
				if (!validateEmbeddings(probe_id - 1, probe_id))
					throw new IOException ("Embedding of MM probe at line " +
						ln + " is not compatible with PM probe.");

				// check that PM probe is a located on row (r-1), column c
				if (spot[r-1][c] != probe_id - 1)
					throw new IOException ("MM probe on line " + ln +
						" does not correspond to PM probe on previous line.");

				// check that PM and MM probes are both fixed or both non-fixed
				if (fixed ^ isFixedSpot (r - 1, c))
					throw new IOException ("Probe pair at lines " + (ln - 1) +
						" and " + ln + " have different fixed flags.");

				setPMProbe (probe_id, false);
			}
			else
				throw new IOException ("Invalid probe type '" + type +
									"' at row " + r + ", column " + c + ".");

			last_type = type;
		}

		// check number of probes
		if (probe_id + 1 != num_probes)
			throw new IOException ("Only " + (probe_id + 1) + " of the " +
									num_probes + " probes were found.");

		// save list of fixed probes as a normal int array
		this.fixed_probe = new int [fixed_list.size()];
		for (Integer id : fixed_list)
			this.fixed_probe[i++] = id;

		// set uninitialized spots as empty
		for (r = 0; r < num_rows; r++)
			for (c = 0; c < num_cols; c++)
				if (spot[r][c] == UNINITIALIZED_SPOT)
					spot[r][c] = EMPTY_SPOT;

		// reading successful
		input_done = true;
	}

	/**
	 * Checks that the embeddings of a pair of probes are valid, i.e. the
	 * resulting probe sequences only differ in the middle base
	 * ({@link #AFFY_MIDDLE_BASE}.
	 *
	 * @param p1 ID of one probe instance of a probe pair
	 * @param p2 ID of the other probe instance of the pair
	 * @return true if embeddings are valid, false otherwise
	 */
	protected boolean validateEmbeddings (int p1, int p2)
	{
		int diff = 0;

		// count the number of differing bits
		for (int i = 0; i < embed[p1].length; i++)
			diff += Integer.bitCount(embed[p1][i] ^ embed[p2][i]);

		// there must be two differences
		// (due to the middle bases)
		if (diff != 2)
			return false;
		else
			return true;

		// TO DO:
		// In reality, this test is not enough as we do not check where
		// the difference is occurring. A correct implementation would also
		// check if the number of differences before and after the middle
		// base is zero.
	}

	/**
	 * Returns a list of non-fixed probe IDs. These are probes not located
	 * on fixed spots and which, therefore, can be relocated by a
	 * {@linkplain PlacementAlgorithm}. The resulting list contain only the
	 * PM probe IDs since the MM IDs can be easily deduced.
	 *
	 * @return a list of movable (not fixed) PM probe IDs
	 */
	public int[] getMovableProbes ()
	{
		int i, f, m, movable[];

		// check if chip spec has already been input
		if (!input_done)
			throw new IllegalStateException
				("Layout specification has not been loaded yet.");

		// only ID of non-fixed PM probes are stored
		movable = new int [(num_probes / 2) - fixed_probe.length];

		// copy (PM) probe IDs which are not
		// in the list of fixed probes
		for (i = f = m = 0; f < fixed_probe.length; i += 2)
		{
			if (i < fixed_probe[f])
				// probe is not fixed: add ID to the list
				movable[m++] = i;

			else
				// means that (i == fixed_probe[f])
				// probe is fixed: ignore ID and move f pointer
				f++;
		}

		// end of fixed list: from now on, all
		// probes are guaranteed to be movable
		for (; i < num_probes; i += 2)
			movable[m++] = i;

		return movable;
	}

	/**
	 * Print the specification of the chip's current layout. The layout is
	 * written in accordance with the format specified by the
	 * {@link #readLayout} method.
	 *
	 * @param out a PrintWriter stream
	 * @throws IOException if an error occurs while writing on the stream
	 */
	public void writeLayout (PrintWriter out) throws IOException
	{
		int		mask = 0, w, pos;
		char	fix, typ;

		for (int c = 0; c < num_cols; c++)
		for (int r = 0; r < num_rows; r++)
		{
			fix = isFixedSpot(r, c) ? 'Y' : 'N';

			if (spot[r][c] == EMPTY_SPOT)
			{
				// print empty spot
				out.println (c + "\t" + r + "\tEMPTY\t" + fix + "\t-\t-\t-");
			}
			else
			{
				typ = isPMProbe(spot[r][c]) ? AFFY_PM_PROBE : AFFY_MM_PROBE;

				// print spot coordinates
				out.print (c + "\t" + r + "\t-\t" + fix + "\t" + typ + "\t");

				// print probe
				for (w = -1, pos = 0; pos < embed_len; pos++)
				{
					if (pos % Integer.SIZE == 0)
					{
						// use next 4-byte word
						w++;

						// turn on very first bit of mask only
						mask = 0x01 << (Integer.SIZE - 1);
					}

					if ((mask & embed[spot[r][c]][w]) != 0)
						out.print (dep_seq[pos]);

					mask >>>= 1;
				}

				out.print ('\t');

				// print embedding
				for (w = -1, pos = 0; pos < embed_len; pos++)
				{
					if (pos % Integer.SIZE == 0)
					{
						// use next 4-byte word
						w++;

						// turn on very first bit of mask only
						mask = 0x01 << (Integer.SIZE - 1);
					}

					if ((mask & embed[spot[r][c]][w]) != 0)
						out.print (dep_seq[pos]);
					else
						out.print (' ');

					mask >>>= 1;
				}

				// end of line
				out.println ();
			}
		}

		// end of output
		out.flush();
	}

	/**
	 * Creates and returns a copy of this AffymetrixChip object. The new object
	 * will contain the same chip specification as the cloned object, i.e.
	 * same chip dimension, same PM/MM pairs of probes on the same locations
	 * and with the same embeddings, same deposition sequence and so on.
	 *
	 * @return a clone of this instance
	 */
	public AffymetrixChip clone ()
	{
		AffymetrixChip c;
		
		// use superclass' clone method
		c = (AffymetrixChip) super.clone();
		
		// now clone the BitSet containing the PM flags
		c.pm_probe = (BitSet) this.pm_probe.clone();
		
		return c;
	}

	/**
	 * Indicates whether some other AffymetrixChip instance is "equal to" this
	 * one. Equals here means that both chips have the same specification, i.e.
	 * both have same dimension, same probes on the same locations and with the
	 * same embeddings, same deposition sequence and so on. The
	 * {link #compatible} method provides another way of comparing two chip
	 * instances.
	 *
	 * @param obj the reference object with which to compare
	 * @return true if this Chip is "equal to" the argument; false otherwise
	 */
	public boolean equals (Object obj)
	{
		AffymetrixChip other;

		// check if superclass members are equal
		if (!super.equals(obj)) return false;
		
		// cannot compare objects that are not AffymetrixChip instances
		if (!(obj instanceof AffymetrixChip)) return false;
		
		other = (AffymetrixChip) obj;
		
		// if specifications have not been loaded yet,
		// then the chips are considered equal
		if (!input_done) return true;
		
		// check the list of PM probe IDs
		if (this.pm_probe == null)
		{
			if(!(other.pm_probe == null)) return false;
		}
		else
		{
			if (!(this.pm_probe.equals(other.pm_probe))) return false;
		}
		
		// if passed all tests, they are considered equal
		return true;
	}

	/**
	 * Indicates whether some other AffymetrixChip instance is
	 * "compatible with" this one. Compatible here means that both chips have
	 * almost the same specification, i.e. both have same dimension, same
	 * probes, same deposition sequence and so on. However, two differences
	 * are allowed: 1) the probes might be located on different spots; and 2)
	 * the same probe on one chip might have a different but "compatible"
	 * embeddings on the other chip). Compatible embeddings means that they
	 * produce the same probe sequence under the given deposition sequence
	 * (which must 	 * be the same for both chips); see
	 * {link #compatibleEmbedding}. Note that this method is similar but
	 * conceptually different from the {link #equals} method.
	 *
	 * @param other the reference object with which to compare
	 * @return true if this Chip is "compatible with" the argument; false
	 * otherwise
	 */
	public boolean compatible (AffymetrixChip other)
	{
		int i, j;

		// check the superclass members for compatibility
		if (!super.compatible(other)) return false;
		
		// if specifications have not been loaded yet,
		// then the chips are considered compatible
		if (!input_done) return true;
		
		// check the list of PM probes
		if (this.pm_probe == null)
		{
			if(!(other.fixed_spots == null)) return false;
		}
		else
		{
			if (!(this.fixed_spots.equals(other.fixed_spots))) return false;
		}
		
		// if passed all test, they are considered equal
		return true;
	}

	/**
	 * Checks that the layout specification is valid. Several test are
	 * performed. For instace, it is checked whether all probe IDs are assigned
	 * to a (single) spot, and whether all probes on fixed spots have their IDs
	 * on the list of fixed probes.
	 *
	 * A chip whose specification has not been loaded yet (by the
	 * {link #readSpecification} method) is NOT considered to be in a valid
	 * state.
	 *
	 * @return true if the current layout specification is valid; false
	 * otherwise
	 */
	public boolean validateLayout ()
	{
		HashSet<Integer>	fixed_id;
		BitSet				visited_id;
		int					i, j, probe_id;
		
		// check basic properties
		if (spot == null || fixed_spots == null || fixed_probe == null)
			return false;
		
		// if specification has not been loaded yet,
		// the chip is considered in an INVALID state
		if (!input_done) return false;
		
		// create a BitSet object to check whether all probe IDs are placed (and that
		// they are placed only once); a zero means the corresponding ID has not been
		// seen yet while a 1 flags a "visited" ID
		visited_id = new BitSet (num_probes);
		
		// create a HashSet to check whether all probes on fixed spots are
		// marked as fixed (and that only these probes are marked as fixed);
		// its size must be twice the number of elements on the fixed list
		// since this list only stores the main ID of a probe pair
		fixed_id = new HashSet<Integer> (2 * fixed_probe.length);
		
		// the set is initially populated with all IDs found on the list of
		// fixed probes; then, while the spot matrix is scanned, every probe
		// found on a fixed spot will have its ID removed from the set; in
		// the end, the set must be empty (otherwise, an ID in the fixed list
		// is not placed on a fixed spot)
		for (i = 0; i < fixed_probe.length; i++)
		{
			// add the PM and MM probe's IDs
			fixed_id.add (fixed_probe[i]);
			fixed_id.add (fixed_probe[i] + 1);
		}

		// check spot matrix			
		if (spot.length != num_rows) return false;
		
		for (i = 0; i < num_rows; i++)
		{
			if (spot[i].length != num_cols)
			{
				System.err.println ("\tbad num_cols!");
				return false;
			}
			
			for (j = 0; j < num_cols; j++)
			{
				probe_id = spot[i][j];
				
				// check if ID is in the valid range
				if (probe_id >= 0 && probe_id < num_probes)
				{
					// check that ID has not been visited already (i.e. that
					// the same probe does not appear in more than one spot)
					if (visited_id.get(probe_id))
						return false;
					
					// mark probe ID as "visited"
					visited_id.set(probe_id);
					
					// if spot is fixed...
					if (isFixedSpot(i, j))
					{
						// ...remove probe ID from the hash set
						if (!fixed_id.remove(probe_id))
						{
							// if it fails, then there is something wrong as
							// the probe placed on the fixed spot does not have
							// its ID in the fixed list
							return false;
						}
					}
				}
				else
				{
					// the only other possible value is the EMPTY flag
					if (probe_id != EMPTY_SPOT)
						return false;
				}
			}
		}
		
		// if the cardinality is not equal to the number of probes,
		// one or more probes are not placed
		if (visited_id.cardinality() != num_probes) return false;
		
		// if the fixed list is not empty, then one or more fixed
		// probes were not found on fixed spots
		if (!fixed_id.isEmpty()) return false;

		// valid layout if passed all tests
		return true;
	}
}
