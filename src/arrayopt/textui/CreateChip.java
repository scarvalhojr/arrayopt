/*
 * CreateChip.java
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
 * ArrayOpt is distributed in the hope that it will be useful, but WITHOUT ANY
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

import arrayopt.layout.*;
import arrayopt.util.*;
import java.io.*;

public class CreateChip
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
		ChipGenerator	generator;
		Chip			chip;
		String			dep_seq;
		int				type, generation, rows, cols, probes, probe_len;

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
			
			// chip type: 'simple' or 'affy' 
			if (args[1].equalsIgnoreCase("ALL_ASYNC_EMBED"))
			{
				generation = ChipGenerator.ALL_ASYNC_EMBEDDINGS;
			}
			else if (args[1].equalsIgnoreCase("ALL_SYNC_EMBED"))
			{
				generation = ChipGenerator.ALL_SYNC_EMBEDDINGS;
			}
			else
				throw new IllegalArgumentException ("unknown generation '" +
					args[1] + "'");
			
			rows       = Integer.parseInt(args[2]);
			cols       = Integer.parseInt(args[3]);
			probes     = Integer.parseInt(args[4]);
			probe_len  = Integer.parseInt(args[5]);
			dep_seq	   = args[6];
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
		
		// check use of standard deposition sequences
		if (dep_seq.equalsIgnoreCase("AFFY"))
			dep_seq = AFFY_DEP_SEQ;
		else if (dep_seq.equalsIgnoreCase("SYNC"))
			dep_seq = SYNC_DEP_SEQ;		

		try
		{
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
		}
		catch (IllegalArgumentException e)
		{
			System.err.println("Unable to instantiate chip: " + e.getMessage());
			System.exit(1);
			return;
		}
		
		System.err.println("Generating chip...");
		generator = new ChipGenerator (chip);
		
		try
		{
			generator.generateLayout(generation);
		}
		catch (IllegalArgumentException e)
		{
			System.err.println("Unable to generate chip: " +
					e.getMessage());
			System.exit(1);
			return;
		}
		
		// print chip layout
		chip.writeLayout(new PrintWriter(System.out));
		
		System.exit(0);
	}
	
	private static void usage ()
	{
		System.err.println (
	"--------------------------\n" +
	"ArrayOpt Microarray Design\n" +
	"--------------------------\n\n" +
	"Usage: CreateChip (affy | simple) <type> <rows> <columns> <probes> <length> <dep-seq>\n" +
	"where: 'affy'      indicates an Affymetrix chip type\n" +
	"       'simple'    indicates a simple chip type\n" +
	"       <type>      is ALLEMBED for all embeddings of a given length or\n" +
	"                      or ALLPROBES for all probes of a given length\n" +
	"                      (with synchronous embeddings)\n" +
	"       <rows>      is the number of rows in the chip\n" +
	"       <columns>   is the number of columns in the chip\n" +
	"       <probes>    is the number of probe in the chip\n" +
	"       <length>    is the length of probes sequences\n" +
	"       <dep-seq>   is the deposition sequence or\n" +
	"                      AFFY for Affymetrix's sequence or\n" +
	"                      SYNC for a 100-step ACGT repetition\n");	
	}
}
