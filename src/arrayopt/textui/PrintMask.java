/*
 * PrintMask.java
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

package arrayopt.textui;

import java.io.*;

import arrayopt.layout.*;

/**
 * A command-line utility for generating a representation of photolithographic
 * masks of a chip layout in bitmap images (BMP).
 *   
 * @author Sergio A. de Carvalho Jr.
 */
public class PrintMask
{
	private static final String AFFY_DEP_SEQ = "TGCATGCATGCATGCATGCATGCA" +
	"TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG";

	private static final String SYNC_DEP_SEQ = "ACGTACGTACGTACGTACGTACGTACGT" +
	"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT" +
	"ACGTACGTACGTACGTACGT";
	
	private static int SIMPLE_CHIP = 0;
	
	private static int AFFY_CHIP = 1;
	
	public static void main (String args[])
	{
		OutputStream out;
		Chip	chip;
		String	filename, dep_seq, outfile;
		int		type, rows, cols, probes, probe_len, start, end;

		try
		{
			// get command-line arguments
			
			// chip type: 'simple' or 'affy' 
			if (args[0].equalsIgnoreCase("simple"))
			{
				type = SIMPLE_CHIP;
			}
			else if (args[0].equalsIgnoreCase("affy"))
			{
				type = AFFY_CHIP;
			}
			else
				throw new IllegalArgumentException ("unknown chip type '" +
					args[0] + "'");
						
			filename     = args[1];
			rows         = Integer.parseInt(args[2]);
			cols         = Integer.parseInt(args[3]);
			probes       = Integer.parseInt(args[4]);
			probe_len    = Integer.parseInt(args[5]);
			dep_seq	     = args[6];
			
			// check use of standard deposition sequences
			if (dep_seq.equalsIgnoreCase("AFFY"))
				dep_seq = AFFY_DEP_SEQ;
			else if (dep_seq.equalsIgnoreCase("SYNC"))
				dep_seq = SYNC_DEP_SEQ;
			
			if (args[7].equalsIgnoreCase("ALL"))
			{
				start = 0;
				end = dep_seq.length() - 1;
			}
			else if (args[7].matches("\\A\\d+\\z"))
			{
				start = end = Integer.parseInt(args[7]);
			}
			else if (args[7].matches("\\A\\d+-\\d+\\z"))
			{
				int idx = args[7].indexOf("-");
				start = Integer.parseInt(args[7].substring(0, idx));
				end = Integer.parseInt(args[7].substring(idx + 1));
			}
			else
				throw new IllegalArgumentException ("invalid mask range '" +
						args[7] + "'.");
			
			if (start > end || start < 0 || end >= dep_seq.length())
				throw new IllegalArgumentException
					("mask range out of bounds.");
		}
		catch (NumberFormatException e)
		{
			usage();
			System.err.println("ERROR: Invalid numeric argument " +
					e.getMessage());
			System.exit(1);
			return;
		}
		catch (ArrayIndexOutOfBoundsException e)
		{
			usage();
			System.err.println("ERROR: Missing mandatory argument(s)");
			System.exit(1);
			return;
		}
		catch (IllegalArgumentException e)
		{
			usage();
			System.err.println("Illegal argument: " + e.getMessage());
			System.exit(1);
			return;
		}
		
		// create chip
		if (type == SIMPLE_CHIP)
		{
			chip = new SimpleChip (rows, cols, probes, probe_len, dep_seq);
		}
		else // if (type == AFFY_CHIP)
		{
			if (probe_len != AffymetrixChip.AFFY_PROBE_LENGTH)
				throw new IllegalArgumentException
					("Invalid probe length (Affymetrix probes must be " +
					 AffymetrixChip.AFFY_PROBE_LENGTH + " base-long).");

			chip = new AffymetrixChip (rows, cols, probes, dep_seq);
		}
		
		System.err.println("Reading input file '" + filename + "'...");
		
		try
		{
			FileReader file = new FileReader(filename);
			chip.readLayout (file);
			file.close();
		}
		catch (Exception e)
		{
			e.printStackTrace ();
			System.exit(1);
			return;
		}
		
		for (int m = start; m <= end; m++)
		{
			if (m < 10)
				outfile = filename + "_mask0" + m + ".bmp";
			else
				outfile = filename + "_mask" + m + ".bmp";
			
			System.err.println("Writing mask " + m + " to " + outfile);
			
			try
			{
				out = new BufferedOutputStream (new FileOutputStream(outfile));
				
				chip.writeMaskBMP(m, out);
			}
			catch (IOException e)
			{
				e.printStackTrace();
				System.exit(1);
				return;
			}
		}
		
		System.exit(0);
	}
	
	private static void usage ()
	{
		System.err.println (
	"--------------------------\n" +
	"ArrayOpt Microarray Design\n" +
	"--------------------------\n\n" +
	"Usage: PrintMask (affy | simple) <input> <rows> <columns>\n" +
	"          <probes> <length> <dep-seq> [<masks>]\n");
		System.err.println (
	"where: 'affy'      indicates an Affymetrix chip type\n" +
	"       'simple'    indicates a simple chip type\n" +
	"       <input>     is a file name\n" +
	"       <rows>      is the number of rows in the chip\n" +
	"       <columns>   is the number of columns in the chip\n" +
	"       <probes>    is the number of probe in the chip\n" +
	"       <length>    is the length of probes sequences\n" +
	"       <dep-seq>   is the deposition sequence or\n" +
	"                      AFFY for Affymetrix's sequence or\n" +
	"                      SYNC for a 100-step ACGT repetition\n" +
	"       <masks>     is a mask range (e.g. 2-5), a single\n" +
	"                      mask number or ALL for all masks\n");
	}
}
