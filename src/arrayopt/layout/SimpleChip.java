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

import java.io.*;
import java.util.regex.*;

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
	 * Creates a new chip instance.
	 *
	 * @param num_rows number of rows of sites
	 * @param num_cols number of columns of sites
	 * @param num_probes total number of probes (not probe pairs)
	 * @param probe_len length of probe sequences
	 * @param embed_len length of embeddings (deposition sequence)
	 */
	public SimpleChip (int num_rows, int num_cols, int num_probes, int probe_len, int embed_len)
	{
		super (num_rows, num_cols, num_probes, probe_len, embed_len);

		// allocate space for a list of "moveable"-probe IDs
		this.probe_list = new int [num_probes];
	}

	/**
	 * Read a chip specification from a character input stream. The input must
	 * consist of lines of text with the following TAB-delimited fields
	 * describing a single spot of the chip:<BR>
	 * <BR>
	 * 1. Row number<BR>
	 * 2. Column number<BR>
	 * 3. Group<BR>
	 * 4. Fixed: <CODE>Y</CODE> or <CODE>N</CODE><BR>
	 * 5. PM/MM: not used
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
	 * 5th field is not used and its content is just ignored.<P>
	 *
	 * <P>The exact probe sequence appears at field 6th as a sequence of the
	 * characters <CODE>A</CODE>, <CODE>C</CODE>, <CODE>G</CODE> and
	 * <CODE>T</CODE>. The last field must contain a valid embedding of the
	 * probe into the deposition sequence formed by inserting a number of
	 * spaces into the probe sequence so that the resulting string is
	 * "synchronozed" with the deposition sequence (see
	 * {@linkplain #encodeEmbedding}).</P>
	 *
	 * <P>Spots can be listed in any order. Any coordinate pair can appear only
	 * once in the input (a spot listed twice is considered ambiguous and will
	 * be rejected). Every probe in the input receives a unique probe ID.<P>
	 *
	 * @param input an character input stream
	 * @throws IOException if an I/O error occurrs or input is not compliantan
	 */
	public void readSpecification (Reader input) throws IOException
	{
		BufferedReader	in;
		Pattern			parser;
		String			line, field[];
		int				num_fields = 7, ln = 0, r, c, probe_id = -1, state = 0;
		boolean			empty, fixed;

		// check if chip spec has already been input
		if (input_done)
			throw new IllegalStateException ("Chip specification has already been input.");

		// mark all spots as uninitialized and not fixed
		for (r = 0; r < num_rows; r++)
			for (c = 0; c < num_cols; c++)
			{
				this.spot[r][c] = UNINITIALIZED_SPOT;
				this.fixed[r][c] = false;
			}

		// crete a buffered reader to read lines
		in = new BufferedReader (input);

		// lines will be parsed into fields...
		parser = Pattern.compile ("\t");

		// ...and stored in this array of strings
		field = new String [num_fields];

		while ((line = in.readLine()) != null)
		{
			// line number
			ln++;

			// parse fields
			field = parser.split(line, num_fields);

			try
			{
				// spot coordinates
				// field 0: X coordinate -> column
				// field 1: Y coordinate -> row
				c = Integer.parseInt (field[0]);
				r = Integer.parseInt (field[1]);

				// is it fixed?
				fixed = (field[3].equals("Y")) ? true : false;

				// is spot empty?
				empty = (field[6].equals("-")) ? true : false;
			}
			catch (ArrayIndexOutOfBoundsException e)
			{
				// invalid file format
				throw new IOException ("Unable to parse input file at line " + ln + ".");
			}
			catch (NumberFormatException e)
			{
				// invalid file format
				throw new IOException ("Invalid spot coordinates at line " + ln + ".");
			}

			// validate row and column numbers
			if (r < 0 || r >= num_rows || c < 0 || c >= num_cols)
				throw new IOException ("Invalid spot coordinates at line " + ln + ".");

			// check for spot conflict
			if (spot[r][c] != UNINITIALIZED_SPOT)
				throw new IOException ("Spot conflict at row " + r + ", column " + c + ".");

			// check if number of probes has been exceeded
			if (probe_id >= num_probes)
				throw new IOException ("Found more probes than expected.");

			if (fixed)
			{
				// mark spot as fixed
				this.fixed[r][c] = true;
			}

			if (empty)
			{
				// mark spot as empty
				spot[r][c] = EMPTY_SPOT;

				// read next line
				continue;
			}

			// new probe found: mark spot as used
			probe_id++;
			this.spot[r][c] = probe_id;

			try
			{
				// encode probe embedding
				encodeEmbedding (field[6], probe_id);
			}
			catch (IllegalArgumentException e)
			{
				throw new IOException ("Invalid embedding at line " + ln + " (" + e.getMessage() + ").");
			}
		}

		// check number of probes
		if (probe_id + 1 != num_probes)
			throw new IOException ("Only " + (probe_id + 1) + " of the " + num_probes + " probes were found.");

		// reading successful
		input_done = true;
	}

	/**
	 * Build a list of probe IDs whose location are not fixed. This list is
	 * stored in {@link #probe_list} and can be later used by a
	 * {@linkplain PlacementAlgorithm}.
	 */
	protected void buildProbeList ()
	{
		int i = 0;

		for (int r = 0; r < num_rows; r++)
			for (int c = 0; c < num_cols; c++)
				if (spot[r][c] != EMPTY_SPOT && spot[r][c] != UNINITIALIZED_SPOT)
				{
					if (!fixed[r][c])
					{
						// all IDs are added to the list
						probe_list[i++] = spot[r][c];
					}
				}

		// store the number of elements
		this.probe_list_len = i;
	}

	// x-pos (column), y-pos (row), group, fixed?, PM/MM, embedding
	public void printLayout (PrintWriter out) throws IOException
	{
		int mask = 0, w, pos;

		for (int r = 0; r < num_rows; r++)
			for (int c = 0; c < num_cols; c++)
				if (spot[r][c] == EMPTY_SPOT)
				{
					// print empty spot
					out.println(c + "\t" + r + "\tEMPTY\t" + (fixed[r][c] ? 'Y' : 'N') + "\t-\t-\t-");
				}
				else
				{
					// print spot coordinates
					out.print(c + "\t" + r + "\t-\t" + (fixed[r][c] ? 'Y' : 'N') + "\t-\t");

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

					out.println ();
				}

		out.flush();
	}
}
