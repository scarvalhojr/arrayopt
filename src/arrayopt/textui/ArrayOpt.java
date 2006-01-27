/*
 * ArrayOpt.java
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
import arrayopt.qap.*;
import java.io.*;

/**
 *
 */
public class ArrayOpt
{
	private static final String DEFAULT_DEP_SEQ = "TGCATGCATGCATGCATGCATGCA" +
						"TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG";
	
	public static void main (String args[])
	{
		Chip					chip, copy;
		PlacementAlgorithm		placer = null;
		PostPlacementAlgorithm	optimizer = null;
		int						rows, cols, probes, probe_len, unplaced;
		String					chip_type, filename, algorithm, dep_seq;

		try
		{
			// get mandatory command-line arguments
			chip_type = args[0];
			filename  = args[1];
			rows      = Integer.parseInt(args[2]);
			cols      = Integer.parseInt(args[3]);
			probes    = Integer.parseInt(args[4]);
			probe_len = Integer.parseInt(args[5]);
			algorithm = args[6];
		}
		catch (Exception e)
		{
			usage();
			System.exit(1);
			return;
		}

		try
		{
			// get optional command-line arguments
			dep_seq = args[7];
		}
		catch (ArrayIndexOutOfBoundsException e)
		{
			dep_seq = DEFAULT_DEP_SEQ;
		}

		try
		{
			// algorithm selection

			if (algorithm.equalsIgnoreCase("RANDOM"))
				placer = new RandomFiller();

			else if (algorithm.equalsIgnoreCase("GREEDY10"))
				placer = new GreedyFiller(10);
			else if (algorithm.equalsIgnoreCase("GREEDY100"))
				placer = new GreedyFiller(100);

			else if (algorithm.equalsIgnoreCase("GREEDY10R"))
				placer = new GreedyFiller(10, true);
			else if (algorithm.equalsIgnoreCase("GREEDY100R"))
				placer = new GreedyFiller(100, true);

			else if (algorithm.equalsIgnoreCase("2D1-SEQ"))
				placer = new TwoDimensionalPartitioning(new SequentialFiller(), 1);

			else if (algorithm.equalsIgnoreCase("2D4-SEQ"))
				placer = new TwoDimensionalPartitioning(new SequentialFiller(), 4);

			else if (algorithm.equalsIgnoreCase("2D1-GREEDY"))
				placer = new TwoDimensionalPartitioning(new GreedyFiller(), 1);

			else if (algorithm.equalsIgnoreCase("2D4-GREEDY"))
				placer = new TwoDimensionalPartitioning(new GreedyFiller(), 4);

			/*
			else if (algorithm.equalsIgnoreCase("PIVOT-GREEDY"))
				placer = new PivotPartitioning(new GreedyFiller());

			else if (algorithm.equalsIgnoreCase("BESTPIVOT-GREEDY"))
				placer = new BestPivotPartitioning(new GreedyFiller());

			else if (algorithm.equalsIgnoreCase("BESTPIVOT-SEQ"))
				placer = new BestPivotPartitioning(new SequentialFiller());

			else if (algorithm.equalsIgnoreCase("PREPIVOT-SEQ"))
				placer = new PrePivotPartitioning(new SequentialFiller());

			else if (algorithm.equalsIgnoreCase("IMAGINARY1-SEQ"))
				placer = new ImaginaryChipPartitioning(new SequentialFiller(), 1);

			else if (algorithm.equalsIgnoreCase("LOCALIZED7-SEQ"))
				placer = new LocalizedPartitioning(new SequentialFiller(), 7);
			//*/

			else if (algorithm.equalsIgnoreCase("2D4-QAPBB"))
				placer = new TwoDimensionalPartitioning(new QAPOptimization(new BranchAndBound()), 4);

			else if (algorithm.equalsIgnoreCase("QAPGRASPD"))
				placer = new QAPOptimization(new GraspDense());

			else if (algorithm.equalsIgnoreCase("LEFT"))
				placer = new ProbeSetEmbeddingWrapper(new LeftMostEmbedding());
				
			else if (algorithm.equalsIgnoreCase("RIGHT"))
				placer = new ProbeSetEmbeddingWrapper(new RightMostEmbedding());
				
			else if (algorithm.equalsIgnoreCase("CENTER"))
				placer = new ProbeSetEmbeddingWrapper(new CenteredEmbedding());
				
			else if (algorithm.equalsIgnoreCase("PIVOT"))
				placer = new ProbeSetEmbeddingWrapper(new PivotEmbedding());

			else if (algorithm.equalsIgnoreCase("SWIN-GRASPD"))
				optimizer = new SlidingWindowOptimization (
								new QAPOptimization (
									new GraspDense()
							), 6, 2);

			else
				throw new IllegalArgumentException
					("Unknown placement algorithm.");

			// chip type

			if (chip_type.equalsIgnoreCase("simple"))
			{
				chip = new SimpleChip (rows, cols, probes, probe_len, dep_seq);
			}
			else if (chip_type.equalsIgnoreCase("affy"))
			{
				if (probe_len != AffymetrixChip.AFFY_PROBE_LENGTH)
					throw new IllegalArgumentException
						("Invalid probe length (Affymetrix probes must be " +
						 AffymetrixChip.AFFY_PROBE_LENGTH + " base-long).");

				chip = new AffymetrixChip (rows, cols, probes, dep_seq);
			}

			else
				throw new IllegalArgumentException
					("Unknown chip type.");
		}
		catch (Exception e)
		{
			e.printStackTrace ();
			System.exit(1);
			return;
		}
		
		try
		{
			// create a file reader for the input
			FileReader file = new FileReader(filename);

			// read input
			chip.readLayout (file);

			// close reader
			file.close();

		}
		catch (Exception e)
		{
			e.printStackTrace ();
			System.exit(1);
			return;
		}

		// make a copy of the chip before modifying its layout
		if (chip instanceof SimpleChip)
		{
			copy = ((SimpleChip) chip).clone();
		}
		else
		{
			copy = ((AffymetrixChip) chip).clone();
		}

		if (placer != null)
		{
			// re-place probes on the chip
			unplaced = placer.makeLayout(chip);

			if (unplaced > 0)
				System.err.println(unplaced + " unplaced probe(s)");
		}
		
		if (optimizer != null)
		{
			optimizer.optimizeLayout (chip);
		}

		if (!chip.compatible(copy))
			System.err.println("WARNING: new layout is not compatible with the original specification.");			

		try
		{
			// print chip layout
			chip.writeLayout(new PrintWriter(System.out));
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}

		System.exit(0);
	}

	private static void usage ()
	{
		System.out.println (
			"Usage: ArrayOpt [affy | simple] <input file> <# of rows>" +
				" <# of columns> <# of probes> <probe length> " +
				" <placement algorithm> [<deposition sequence>]");
	}
}
