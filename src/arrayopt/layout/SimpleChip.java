/*
 * SimpleChip.java
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
 * This class represents the simplest type of chips. Probes appear in single
 * instances, in contrast to Affymetrix's pairs (see
 * {@linkplain AffymetrixChip}).
 *
 * @author Sergio A. de Carvalho Jr.
 */
public class SimpleChip extends Chip
{
	/**
	 * Creates a new SimpleChip instance.
	 *
	 * @param num_rows number of rows of sites
	 * @param num_cols number of columns of sites
	 * @param num_probes total number of probes (not probe pairs or tuples)
	 * @param probe_len length of the probes
	 * @param dep_seq deposition sequence as a String
	 */
	public SimpleChip (int num_rows, int num_cols, int num_probes,
		int probe_len, String dep_seq)
	{
		super (num_rows, num_cols, num_probes, probe_len, dep_seq);
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
	 * 5. PM/MM: must contain a dash (<CODE>-</CODE>)<BR>
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
	 * 5th field is not used but it must contain a dash.<P>
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
	 * an IOException) but spots can be listed in any order.<P>
	 *
	 * <P>There is an option on whether fixed spots must have its contents
	 * preserved (ignore_fixed = false) or not (ignore_fixed = true). In the
	 * first case, the chip specification will contain a list of fixed spots
	 * and fixed probes (those found on fixed spots). Otherwise, all fixed spots
	 * will be treated as non-fixed (which may result, for instance, in their
	 * their probes being moved to a new location).</P>
	 *
	 * @param input an character input stream
	 * @param ignore_fixed true if fixed status should be ignored, false
	 * otherwise
	 * @throws IOException if an I/O error occurrs or input is not compliantan
	 */
	@Override
	public void readLayout (Reader input, boolean ignore_fixed)
		throws IOException
	{
		ArrayList<Integer>	fixed_list;
		BufferedReader		in;
		Pattern				parser;
		String				line, field[];
		int					ln = 0, r, c, probe_id = -1, i = 0;
		boolean				empty, fixed;

		// check if chip spec has already been input
		if (input_done)
			throw new IllegalStateException
				("Layout specification has already been loaded.");

		// mark all spots as unitialized and not fixed
		for (r = 0; r < num_rows; r++)
			for (c = 0; c < num_cols; c++)
			{
				this.spot[r][c] = UNINITIALIZED_SPOT;
			}

		// create a list of IDs of fixed probes with an initial
		// capacity of about 5% of the number of probes
		fixed_list = new ArrayList<Integer> ((int) (.05 * num_probes));

		// create a buffered reader to read lines
		in = new BufferedReader (input);

		// lines will be parsed into fields...
		parser = Pattern.compile ("\t");

		// ...and stored in this array of strings
		field = new String [7];

		while ((line = in.readLine()) != null)
		{
			// line number
			ln++;
			
			// skip blank lines
			if (line.length() <= 0) continue;

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
				if (ignore_fixed)
				{
					fixed = false;
				}
				else
				{
					if (field[3].equals("Y"))
						fixed = true;
					else if (field[3].equals("N"))
						fixed = false;
					else
						throw new IOException ("Invalid fixed flag at line " +
												ln + ".");
				}
				
				if (!field[4].equals("-"))
					throw new IOException ("Non-empty probe type flag at line "
											+ ln + " (must contain a dash).");
				
				// empty spot?
				empty = (field[6].equals("-")) ? true : false;
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

			if (fixed)
				// add probe ID to the list of fixed probes
				fixed_list.add(probe_id);

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
	 * Returns a list of non-fixed probe IDs. These are the probes not located
	 * on fixed spots and which, therefore, can be relocated by a
	 * {@linkplain PlacementAlgorithm}.
	 *
	 * @return a list of movable (not fixed) probe IDs
	 */
	@Override
	public int[] getMovableProbes ()
	{
		int i, f, m, movable[];

		// check if chip spec has already been input
		if (!input_done)
			throw new IllegalStateException
				("Layout specification has not been loaded yet.");

		movable = new int [num_probes - fixed_probe.length];

		// copy probe IDs which are not
		// in the list of fixed probes
		for (i = f = m = 0; f < fixed_probe.length; i++)
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
		for (; i < num_probes; i++)
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
	@Override
	public void writeLayout (PrintWriter out)
	{
		int		mask = 0, w, pos;
		char	fix;

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
				// print spot coordinates
				out.print (c + "\t" + r + "\t-\t" + fix + "\t-\t");

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
	 * Creates and returns a copy of this SimpleChip object. The new object
	 * will contain the same chip specification as the cloned object, i.e.
	 * same chip dimension, same probes on the same locations and with the
	 * same embeddings, same deposition sequence and so on.
	 *
	 * @return a clone of this instance
	 */
	@Override
	public SimpleChip clone ()
	{
		// use superclass' clone method as
		// this class has no extra member variables
		return (SimpleChip) super.clone();
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
	@Override
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
		// marked as fixed (and that only these probes are marked as fixed)
		fixed_id = new HashSet<Integer> (fixed_probe.length);
		
		// the set is initially populated with all IDs found on the list of
		// fixed probes; then, while the spot matrix is scanned, every probe
		// found on a fixed spot will have its ID removed from the set; in
		// the end, the set must be empty (otherwise, an ID in the fixed list
		// is not placed on a fixed spot)
		for (i = 0; i < fixed_probe.length; i++)
			fixed_id.add (fixed_probe[i]);
			
		// check spot matrix
		if (spot.length != num_rows) return false;
		
		for (i = 0; i < num_rows; i++)
		{
			if (spot[i].length != num_cols)
				return false;
			
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
							// if it fails, then there is something wrong as
							// the probe placed on the fixed spot does not have
							// its ID in the fixed list
							return false;
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
