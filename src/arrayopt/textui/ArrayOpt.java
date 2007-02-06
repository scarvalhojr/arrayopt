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
	private static final String AFFY_DEP_SEQ = "TGCATGCATGCATGCATGCATGCA" +
						"TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG";

	private static final String SYNC_DEP_SEQ = "ACGTACGTACGTACGTACGTACGTACGT" +
						"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT" +
						"ACGTACGTACGTACGTACGT";

	private static int SIMPLE_CHIP = 0;
	
	private static int AFFY_CHIP = 1;
	
	public static void main (String args[])
	{
		LayoutAlgorithm alg[];
		Chip	chip, copy = null;
		String	filename, dep_seq;
		int		i, type, rows, cols, probes, probe_len, num_alg, a;
		boolean	ignore_fixed, check, calc_bl, calc_blm, calc_ci, print_chip;
		long	bl, start, end, total = 0;
		double	norm_bl;

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
			
			// ignore fixed spots?
			if (args[1].equalsIgnoreCase("fix"))
				ignore_fixed = false;
			else if (args[1].equalsIgnoreCase("nofix"))
				ignore_fixed = true;
			else
				throw new IllegalArgumentException ("'" + args[1] +
					"' (expected 'fix' or 'nofix')");
			
			filename  = args[2];
			rows      = Integer.parseInt(args[3]);
			cols      = Integer.parseInt(args[4]);
			probes    = Integer.parseInt(args[5]);
			probe_len = Integer.parseInt(args[6]);
			dep_seq	  = args[7];

			// perform validation checks?
			if (args[8].equalsIgnoreCase("check"))
				check = true;
			else if (args[8].equalsIgnoreCase("no-check"))
				check = false;
			else
				throw new IllegalArgumentException ("'" + args[8] +
					"' (expected 'check' or 'no-check')");

			// compute and print total border
			// length or average conflict indices? 
			if (args[9].equalsIgnoreCase("calc-bl"))
			{
				calc_bl  = true;
				calc_blm = false;
				calc_ci  = false;
			}
			else if (args[9].equalsIgnoreCase("calc-blm"))
			{
				calc_bl  = false;
				calc_blm = true;
				calc_ci  = false;
			}
			else if (args[9].equalsIgnoreCase("calc-ci"))
			{
				calc_bl  = false;
				calc_blm = false;
				calc_ci  = true;
			}
			else if (args[9].equalsIgnoreCase("no-calc"))
			{
				calc_bl  = false;
				calc_blm = false;
				calc_ci  = false;
			}
			else
				throw new IllegalArgumentException ("'" + args[9] +
				"' (expected 'calc-bl', 'calc-blm', 'calc-ci' or 'no-calc')");
			
			// print produced layout?
			if (args[10].equalsIgnoreCase("print"))
				print_chip = true;
			else if (args[10].equalsIgnoreCase("no-print"))
				print_chip = false;
			else
				throw new IllegalArgumentException ("'" + args[10] +
					"' (expected 'print' or 'no-print')");
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
		
		// create array of algorithms
		num_alg = args.length - 11;
		if (num_alg > 0)
		{
			alg = new LayoutAlgorithm[num_alg];
			
			try
			{
				// create an instance of each algorithm
				for (i = 11, a = 0; i < args.length; a++, i++)
					alg[a] = parseAlgorithmName(args[i]);
			}
			catch (Exception e)
			{
				usage();
				System.err.println("Illegal argument: " + e.getMessage());
				System.exit(1);
				return;
			}
		}
		else
		{
			alg = null;
			check = false;
		}
		
		if (filename.equalsIgnoreCase("RANDOM"))
		{
			System.err.println("Generating random chip...");
			
			chip.createRandomLayout();
		}
		else
		{
			System.err.println("Reading input file '" + filename + "'...");
			
			try
			{
				FileReader file = new FileReader(filename);
				chip.readLayout (file, ignore_fixed);
				file.close();
	
			}
			catch (FileNotFoundException e)
			{
				System.err.println("Unable to open input file for reading: " +
						e.getMessage());
				System.exit(1);
				return;
			}
			catch (IOException e)
			{
				System.err.println("Error while reading input file: " +
						e.getMessage());
				System.exit(1);
				return;
			}
		}

		if (check)
		{
			System.err.println("Cloning chip for later checks...");
			
			// make a copy of the chip before modifying its layout
			if (chip instanceof SimpleChip)
			{
				copy = ((SimpleChip) chip).clone();
			}
			else
			{
				copy = ((AffymetrixChip) chip).clone();
			}
		}
		
		for (a = 0; a < num_alg; a++)
		{
			if (calc_bl)
			{
				bl = LayoutEvaluation.borderLength(chip);
				norm_bl = bl / (double) chip.getNumberOfBorders();
				System.err.println("Total border length: " + bl +
						" -> normalized: " + norm_bl);
			}
			else if (calc_ci)
			{
				System.err.println("Average conflict index: " +
						LayoutEvaluation.averageConflictIndex(chip));
			}
			
			System.err.println("Running " + alg[a] + "...");
			
			start = System.nanoTime();

			// re-place probes on the chip
			alg[a].changeLayout(chip);
			
			end = System.nanoTime();
			
			total += end - start;
			
			System.err.println("Elapsed time: " +
					(end - start)/Math.pow(10,9) + " sec");
		}

		if (calc_bl)
		{
			bl = LayoutEvaluation.borderLength(chip);
			norm_bl = bl / (double) chip.getNumberOfBorders();
			System.err.println("Total border length: " + bl +
					" -> normalized: " + norm_bl);
		}
		else if (calc_ci)
		{
			System.err.println("Average conflict index: " +
					LayoutEvaluation.averageConflictIndex(chip));
		}

		if (total > 0)
			System.err.println("Total time: " + total/Math.pow(10,9) + " sec");

		if (copy != null)
		{
			if (chip.equals(copy))
				System.err.println("WARNING: new layout is equal to the "+
						"original specification.");
			else
				System.err.println("New layout is not equal to the " +
						"original specification.");

			if (!chip.compatible(copy))
				System.err.println("WARNING: new layout is NOT compatible " +
						"with the original specification.");
			else
				System.err.println("New layout is compatible " +
						"with the original specification.");
		}
		
		if (print_chip)
		{
			// print chip layout
			chip.writeLayout(new PrintWriter(System.out));
		}
		
		if (calc_blm)
		{
			// print border length per masking step
			LayoutEvaluation.analyzeBorderLength(chip,
					new PrintWriter(System.out));
		}
		
		System.exit(0);
	}
	
	private static LayoutAlgorithm parseAlgorithmName (String name)
	{
		LayoutAlgorithm alg;
		String args[];
		
		// parse algorithm options
		args = name.split("-");

		// *****************
		// Sequential placer
		// *****************
		if (args[0].equalsIgnoreCase("SEQPLACER"))
		{
			int order;
			
			if (args.length != 2)
				throw new IllegalArgumentException
					("Missing arguments for Sequential placement algorithm.");
			
			if (args[1].equalsIgnoreCase("KEEP"))
				order = SequentialPlacer.KEEP_ORDER;
			else if (args[1].equalsIgnoreCase("RANDOM"))
				order = SequentialPlacer.RANDOM_ORDER;
			else if (args[1].equalsIgnoreCase("SORTSEQ"))
				order = SequentialPlacer.SORT_SEQUENCES;
			else if (args[1].equalsIgnoreCase("SORTEMBED"))
				order = SequentialPlacer.SORT_EMBEDDINGS;
			else if (args[1].equalsIgnoreCase("TSP"))
				order = SequentialPlacer.TSP_ORDER;
			else
				throw new IllegalArgumentException ("Unknown ordering '" +
						args[1] + "' for Sequential placement algorithm.");
			
			alg = new SequentialPlacer(order);
		}		

		// ******************
		// k-threading placer
		// ******************
		else if (args[0].equalsIgnoreCase("KTHREAD"))
		{
			int kvalue, order;
			
			if (args.length != 3)
				throw new IllegalArgumentException
					("Missing arguments for k-threading placement algorithm.");
			
			if (args[1].matches("\\A\\d+\\z"))
				kvalue = Integer.parseInt(args[1]);
			else
				throw new IllegalArgumentException ("Invalid k-value '" +
						args[1] + "' for k-threading placement algorithm.");
			
			if (args[2].equalsIgnoreCase("KEEP"))
				order = KThreadingPlacer.KEEP_ORDER;
			else if (args[2].equalsIgnoreCase("RANDOM"))
				order = KThreadingPlacer.RANDOM_ORDER;
			else if (args[2].equalsIgnoreCase("SORTSEQ"))
				order = KThreadingPlacer.SORT_SEQUENCES;
			else if (args[2].equalsIgnoreCase("SORTEMBED"))
				order = KThreadingPlacer.SORT_EMBEDDINGS;
			else if (args[2].equalsIgnoreCase("TSP"))
				order = KThreadingPlacer.TSP_ORDER;
			else
				throw new IllegalArgumentException ("Unknown ordering '" +
						args[2] + "' for k-threading placement algorithm.");
			
			alg = new KThreadingPlacer(kvalue, order);
		}

		// ****************
		// Greedy placement
		// ****************
		else if (args[0].equalsIgnoreCase("GREEDYPLACER"))
		{
			int mode, window, order, kvalue;
			
			if (args.length != 5)
				throw new IllegalArgumentException
					("Missing arguments for Greedy placement algorithm.");
			
			if (args[1].equalsIgnoreCase("BL"))
				mode = GreedyPlacer.BORDER_LENGTH_MIN;
			else if (args[1].equalsIgnoreCase("CI"))
				mode = GreedyPlacer.CONFLICT_INDEX_MIN;
			else
				throw new IllegalArgumentException ("Unknown '" + args[1] +
						"' mode for Greedy placement algorithm.");
			
			if (args[2].matches("\\A\\d+\\z"))
				window = Integer.parseInt(args[2]);
			else if (args[2].matches("\\A\\d+[Kk]\\z"))
				window = 1000 * Integer.parseInt(
									args[2].substring(0,args[2].length()-1));
			else
				throw new IllegalArgumentException
					("Invalid window-size value '" + args[2] +
						"' for Greedy placement algorithm.");
			
			if (args[3].matches("\\A\\d+\\z"))
				kvalue = Integer.parseInt(args[3]);
			else
				throw new IllegalArgumentException ("Invalid k-value '" +
						args[3] + "' for Greedy placement algorithm.");
			
			if (args[4].equalsIgnoreCase("KEEP"))
				order = GreedyPlacer.KEEP_ORDER;
			else if (args[4].equalsIgnoreCase("RANDOM"))
				order = GreedyPlacer.RANDOM_ORDER;
			else if (args[4].equalsIgnoreCase("SORTSEQ"))
				order = GreedyPlacer.SORT_SEQUENCES;
			else if (args[4].equalsIgnoreCase("SORTEMBED"))
				order = GreedyPlacer.SORT_EMBEDDINGS;
			else if (args[4].equalsIgnoreCase("TSP"))
				order = GreedyPlacer.TSP_ORDER;
			else
				throw new IllegalArgumentException ("Unknown ordering '" +
						args[4] + "' for Greedy placement algorithm.");
			
			alg = new GreedyPlacer(mode, window, kvalue, order);
		}		
		
		// *********************
		// Greedy Plus placement
		// *********************
		else if (args[0].equalsIgnoreCase("GREEDYPLUS"))
		{
			int mode, window, kvalue;
			
			if (args.length != 4)
				throw new IllegalArgumentException
					("Missing arguments for Greedy Plus placement algorithm.");
			
			if (args[1].equalsIgnoreCase("BL"))
				mode = OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN;
			else if (args[1].equalsIgnoreCase("CI"))
				mode = OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN;
			else
				throw new IllegalArgumentException ("Unknown '" + args[1] +
						"' mode for Greedy Plus placement algorithm.");
			
			if (args[2].matches("\\A\\d+\\z"))
				window = Integer.parseInt(args[2]);
			else if (args[2].matches("[0-9]+[Kk]"))
				window = 1000 * Integer.parseInt(
									args[2].substring(0,args[2].length()-1));
			else
				throw new IllegalArgumentException
					("Invalid window-size value '" + args[2] +
						"' for Greedy Plus placement algorithm.");
			
			if (args[3].matches("[0-9]+"))
				kvalue = Integer.parseInt(args[3]);
			else
				throw new IllegalArgumentException ("Invalid k-value '" +
						args[3] + "' for Greedy PLus placement algorithm.");
			
			alg = new GreedyPlusPlacer(mode, window, kvalue);
		}		
		
		// *************
		// Row-epitaxial
		// *************
		else if (args[0].equalsIgnoreCase("REPTX"))
		{
			int mode, look_ahead;
			
			if (args.length != 3)
				throw new IllegalArgumentException
				("Missing arguments for Row-epitaxial placement algorithm.");
			
			if (args[1].equalsIgnoreCase("BL"))
				mode = RowEpitaxial.BORDER_LENGTH_MIN;
			else if (args[1].equalsIgnoreCase("CI"))
				mode = RowEpitaxial.CONFLICT_INDEX_MIN;
			else
				throw new IllegalArgumentException ("Unknown '" + args[1] +
						"' mode for Row-epitaxial placement algorithm.");
			
			if (args[2].matches("[0-9]+"))
				look_ahead = Integer.parseInt(args[2]);
			else if (args[2].matches("[0-9]+[Kk]"))
				look_ahead = 1000 * Integer.parseInt(
									args[2].substring(0,args[2].length()-1));
			else
				throw new IllegalArgumentException
					("Invalid look-ahead value '" + args[2] +
						"' for Row-epitaxial placement algorithm.");
			
			alg = new RowEpitaxial(mode, look_ahead);
		}		

		// ***************
		// 1D-Partitioning
		// ***************
		else if (args[0].equalsIgnoreCase("1DPART"))
		{
			FillingAlgorithm filler;
			String filler_name;
			int stop_dim, idx;
			
			if (args.length < 3)
				throw new IllegalArgumentException
					("Missing arguments for 1D Partitioning.");
			
			if (args[1].matches("[0-9]+"))
				stop_dim = Integer.parseInt(args[1]);
			else
				throw new IllegalArgumentException ("Invalid stop dimension " +
						args[1] + " for 1D Partitioning.");
			
			// get filling algorithm's name
			idx = args[0].length() + args[1].length() + 2;
			filler_name = name.substring(idx);
			
			try
			{
				filler = (FillingAlgorithm) parseAlgorithmName(filler_name);
			}
			catch (ClassCastException e)
			{
				throw new IllegalArgumentException ("Cannot use '" + filler_name
						+ "' as a filling algorithm for 1D Partitioning.");
			}
			catch (Exception e)
			{
				throw new IllegalArgumentException ("Unable to instantiate " +
						"filling algorithm for 1D Partitioning: " +
						e.getMessage());				
			}
			
			alg = new OneDimensionalPartitioning (filler, stop_dim);
		}
		
		// ***************
		// 2D-Partitioning
		// ***************
		else if (args[0].equalsIgnoreCase("2DPART"))
		{
			FillingAlgorithm filler;
			String filler_name;
			int stop_dim, idx;
			
			if (args.length < 3)
				throw new IllegalArgumentException
					("Missing arguments for 2D Partitioning.");
			
			if (args[1].matches("[0-9]+"))
				stop_dim = Integer.parseInt(args[1]);
			else
				throw new IllegalArgumentException ("Invalid stop dimension " +
						args[1] + " for 2D Partitioning.");
			
			// get filling algorithm's name
			idx = args[0].length() + args[1].length() + 2;
			filler_name = name.substring(idx);
			
			try
			{
				filler = (FillingAlgorithm) parseAlgorithmName(filler_name);
			}
			catch (ClassCastException e)
			{
				throw new IllegalArgumentException ("Cannot use '" + filler_name
						+ "' as a filling algorithm for 2D Partitioning.");
			}
			catch (Exception e)
			{
				throw new IllegalArgumentException ("Unable to instantiate " +
						"filling algorithm for 2D Partitioning: " +
						e.getMessage());				
			}
			
			alg = new TwoDimensionalPartitioning (filler, stop_dim);
		}
		
		// **********************************************
		// 2D-Partitioning with central mask optimization
		// **********************************************
		else if (args[0].equalsIgnoreCase("2DCENTRALPART"))
		{
			FillingAlgorithm filler;
			String filler_name;
			int stop_dim, idx;
			
			if (args.length < 3)
				throw new IllegalArgumentException
					("Missing arguments for 2D Partitioning.");
			
			if (args[1].matches("[0-9]+"))
				stop_dim = Integer.parseInt(args[1]);
			else
				throw new IllegalArgumentException ("Invalid stop dimension " +
						args[1] + " for 2D Partitioning.");
			
			// get filling algorithm's name
			idx = args[0].length() + args[1].length() + 2;
			filler_name = name.substring(idx);
			
			try
			{
				filler = (FillingAlgorithm) parseAlgorithmName(filler_name);
			}
			catch (ClassCastException e)
			{
				throw new IllegalArgumentException ("Cannot use '" + filler_name
						+ "' as a filling algorithm for 2D Partitioning.");
			}
			catch (Exception e)
			{
				throw new IllegalArgumentException ("Unable to instantiate " +
						"filling algorithm for 2D Partitioning: " +
						e.getMessage());				
			}
			
			alg = new TwoDimensionalCentralPartitioning (filler, stop_dim);
		}
		
		// ******************
		// Pivot Partitioning
		// ******************
		else if (args[0].equalsIgnoreCase("PIVOTPART"))
		{
			FillingAlgorithm filler;
			String filler_name;
			
			int mode, max_depth, idx;
			
			if (args.length < 4)
				throw new IllegalArgumentException
					("Missing arguments for Pivot Partitioning.");
			
			if (args[1].equalsIgnoreCase("BL"))
				mode = OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN;
			else if (args[1].equalsIgnoreCase("CI"))
				mode = OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN;
			else
				throw new IllegalArgumentException ("Unknown '" + args[1] +
						"' mode for Pivot Partitioning algorithm.");
			
			if (args[2].matches("[0-9]+"))
				max_depth = Integer.parseInt(args[2]);
			else
				throw new IllegalArgumentException ("Invalid partitioning " +
						"depth " + args[2] + " for Pivot Partitioning.");
			
			// get filling algorithm's name
			idx = args[0].length() + args[1].length() + args[2].length() + 3;
			filler_name = name.substring(idx);
			
			try
			{
				filler = (FillingAlgorithm) parseAlgorithmName(filler_name);
			}
			catch (ClassCastException e)
			{
				throw new IllegalArgumentException ("Cannot use '" + filler_name
						+ "' as a filling algorithm for Pivot Partitioning.");
			}
			catch (Exception e)
			{
				throw new IllegalArgumentException ("Unable to instantiate " +
						"filling algorithm for Pivot Partitioning: " +
						e.getMessage());				
			}
			
			alg = new PivotPartitioning (filler, mode, max_depth);
		}
		
		// ******************************
		// Simple re-embedding algorithms
		// ******************************
		else if (name.equalsIgnoreCase("LEFTMOST"))
			alg = new LeftMostEmbedding();
		
		else if (name.equalsIgnoreCase("RIGHTMOST"))
			alg = new RightMostEmbedding();
		
		else if (name.equalsIgnoreCase("CENTERED"))
			alg = new CenteredEmbedding();
		
		else if (name.equalsIgnoreCase("ALIGNED"))
			alg = new AlignedEmbedding();
		
		// ***********************
		// Sequential re-embedding
		// ***********************
		else if (args[0].equalsIgnoreCase("SEQREEMBED"))
		{
			int mode;
			boolean reset;
			double threshold;
			
			if (args.length != 4)
				throw new IllegalArgumentException
					("Missing arguments for Sequential re-embedding.");
			
			if (args[1].equalsIgnoreCase("BL"))
				mode = OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN;
			else if (args[1].equalsIgnoreCase("CI"))
				mode = OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN;
			else
				throw new IllegalArgumentException ("Unknown '" + args[1] +
						"' mode for Sequential re-embedding algorithm.");
			
			if (args[2].equalsIgnoreCase("RESET"))
				reset = true;
			else if (args[2].equalsIgnoreCase("NORESET"))
				reset = false;
			else
				throw new IllegalArgumentException ("Invalid argument '" +
						args[2] + "' for Sequential re-embedding algorithm.");

			if (args[3].matches("([0-9]+(\\.[0-9]*)?)|(\\.[0-9]+)"))
				threshold = Double.parseDouble(args[3]);
			else
				throw new IllegalArgumentException ("Invalid threshold/passes '"
						+ args[3] + "' for Sequential re-embedding algorithm.");
			
			alg = new SequentialReembedding(mode, reset, threshold);
		}
		
		// *********************
		// Priority re-embedding
		// *********************
		else if (args[0].equalsIgnoreCase("PRIORITY"))
		{
			int mode, priority;
			boolean reset;
			double threshold;
			
			if (args.length != 5)
				throw new IllegalArgumentException
					("Missing arguments for Priority re-embedding.");
			
			if (args[1].equalsIgnoreCase("BL"))
				mode = OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN;
			else if (args[1].equalsIgnoreCase("CI"))
				mode = OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN;
			else
				throw new IllegalArgumentException ("Unknown '" + args[1] +
						"' mode for Priority re-embedding algorithm.");
			
			if (args[2].equalsIgnoreCase("NUMEMBED"))
				priority = PriorityReembedding.PRIORITY_NUM_OF_EMBEDDINGS;
			else if (args[2].equalsIgnoreCase("NUMNEIGHBORS"))
				priority =  PriorityReembedding.PRIORITY_NUM_OF_NEIGHBORS;
			else if (args[2].equalsIgnoreCase("BALANCED"))
				priority = PriorityReembedding.PRIORITY_BALANCED; 
			else
				throw new IllegalArgumentException ("Unknown priority '" +
						args[2] + "' for Priority re-embedding algorithm.");

			if (args[3].equalsIgnoreCase("RESET"))
				reset = true;
			else if (args[3].equalsIgnoreCase("NORESET"))
				reset = false;
			else
				throw new IllegalArgumentException ("Invalid argument '" +
						args[3] + "' for Priority re-embedding algorithm.");

			if (args[4].matches("([0-9]+(\\.[0-9]*)?)|(\\.[0-9]+)"))
				threshold = Double.parseDouble(args[4]);
			else
				throw new IllegalArgumentException ("Invalid threshold/passes '"
						+ args[4] + "' for Priority re-embedding algorithm.");
			
			alg = new PriorityReembedding(mode, priority, reset, threshold);
		}		

		else
			throw new IllegalArgumentException
				("Unknown algorithm '" + name + "'");
		
		return alg;
	}

	private static void usage ()
	{
		System.err.println (
	"--------------------------\n" +
	"ArrayOpt Microarray Design\n" +
	"--------------------------\n\n" +
	"Usage: ArrayOpt (affy | simple) (fix | nofix) <input> <rows> <columns> " +
	                                          "<probes> <length> <dep-seq>\n" +
	"          (check | no-check) (calc-bl | calc-blm | calc-ci | no-calc) "  +
	                                          "(print | no-print) <alg>*\n\n" +
	"where: 'affy'      indicates an Affymetrix chip type\n" +
	"       'simple'    indicates a simple chip type\n" +
	"       'fix'       considers fixed spots\n" +
	"       'no-fix'    ignores fixed spots in the input\n" +
	"       <input>     is a file name or\n" +
	"                      RANDOM for a randomly generated chip\n" +
	"       <rows>      is the number of rows in the chip\n" +
	"       <columns>   is the number of columns in the chip\n" +
	"       <probes>    is the number of probe in the chip\n" +
	"       <length>    is the length of probes sequences\n" +
	"       <dep-seq>   is the deposition sequence or\n" +
	"                      AFFY for Affymetrix's sequence or\n" +
	"                      SYNC for a 100-step ACGT repetition\n" +
	"       'check'     performs validation checks\n" +
	"       'no-check'  does not perform validation checks\n" +
	"       'calc-bl'   prints total and normalized border length\n" +
	"       'calc-blm'  prints border length per masking step\n" +
	"       'calc-ci'   prints average conflict index\n" +
	"       'no-calc'   does not print any quality measure\n" +
	"       'print'     prints the resulting layout on standard output\n" +
	"       'no-print'  does not print the resulting layout\n" +
	"       <alg>*      zero or more layout algorithms\n");	
	}
}
