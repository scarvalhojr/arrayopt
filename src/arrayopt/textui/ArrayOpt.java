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
		boolean	ignore_fixed, check, calc_bl, calc_ci, print_chip;
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
			
			filename     = args[2];
			rows         = Integer.parseInt(args[3]);
			cols         = Integer.parseInt(args[4]);
			probes       = Integer.parseInt(args[5]);
			probe_len    = Integer.parseInt(args[6]);
			dep_seq	     = args[7];

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
				calc_bl = true;
				calc_ci = false;
			}
			else if (args[9].equalsIgnoreCase("calc-ci"))
			{
				calc_bl = false;
				calc_ci = true;
			}
			else if (args[9].equalsIgnoreCase("no-calc"))
			{
				calc_bl = false;
				calc_ci = false;
			}
			else
				throw new IllegalArgumentException ("'" + args[9] +
					"' (expected 'calc-bl' or 'calc-ci' or 'no-calc')");
			
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
			catch (Exception e)
			{
				e.printStackTrace ();
				System.exit(1);
				return;
			}
		}

		if (check)
		{
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
			try
			{
				// print chip layout
				chip.writeLayout(new PrintWriter(System.out));
			}
			catch (IOException e)
			{
				e.printStackTrace();
			}
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
		if (args[0].equalsIgnoreCase("SEQPLACER") && args.length == 2)
		{
			int order;
			
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
		else if (args[0].equalsIgnoreCase("KTHREAD") && args.length == 3)
		{
			int kvalue, order;
			
			if (args[1].matches("[0-9]+"))
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
		else if (args[0].equalsIgnoreCase("GREEDYPLACER") && args.length == 5)
		{
			int mode, window, order, kvalue;
			
			if (args[1].equalsIgnoreCase("BL"))
				mode = GreedyPlacer.BORDER_LENGTH_MIN;
			else if (args[1].equalsIgnoreCase("CI"))
				mode = GreedyPlacer.CONFLICT_INDEX_MIN;
			else
				throw new IllegalArgumentException ("Unknown '" + args[1] +
						"' mode for Greedy placement algorithm.");
			
			if (args[2].matches("[0-9]+"))
				window = Integer.parseInt(args[2]);
			else if (args[2].matches("[0-9]+[Kk]"))
				window = 1000 * Integer.parseInt(
									args[2].substring(0,args[2].length()-1));
			else
				throw new IllegalArgumentException
					("Invalid window-size value '" + args[2] +
						"' for Greedy placement algorithm.");
			
			if (args[3].matches("[0-9]+"))
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
		
		// *************
		// Row-epitaxial
		// *************
		else if (args[0].equalsIgnoreCase("REPTX") && args.length == 3)
		{
			int mode, look_ahead;
			
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

		// ********************************************************
		// Leftmost, Rightmost and Centered re-embedding algorithms
		// ********************************************************
		else if (name.equalsIgnoreCase("LEFTMOST"))
			alg = new LeftMostEmbedding();
		
		else if (name.equalsIgnoreCase("RIGHTMOST"))
			alg = new RightMostEmbedding();
		
		else if (name.equalsIgnoreCase("CENTERED"))
			alg = new CenteredEmbedding();
		
		// ***********************
		// Sequential re-embedding
		// ***********************
		else if (name.equalsIgnoreCase("SEQREEMBED-BL-RESET"))
			alg = new SequentialReembedding(
				OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN, true);
		
		else if (name.equalsIgnoreCase("SEQREEMBED-BL-NORESET"))
			alg = new SequentialReembedding(
				OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN, false);
		
		else if (name.equalsIgnoreCase("SEQREEMBED-CI-RESET"))
			alg = new SequentialReembedding(
				OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN, true);
		
		else if (name.equalsIgnoreCase("SEQREEMBED-CI-NORESET"))
			alg = new SequentialReembedding(
				OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN, false);

		// *********************
		// Priority re-embedding
		// *********************
		else if (name.equalsIgnoreCase("PRIORITYREEMBED-BL-NUMEMBED-RESET"))
			alg = new PriorityReembedding(
					OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN,
					PriorityReembedding.PRIORITY_NUM_OF_EMBEDDINGS, true);

		else if (name.equalsIgnoreCase("PRIORITYREEMBED-BL-NUMEMBED-NORESET"))
			alg = new PriorityReembedding(
					OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN,
					PriorityReembedding.PRIORITY_NUM_OF_EMBEDDINGS, false);

		else if (name.equalsIgnoreCase("PRIORITYREEMBED-BL-NUMNEIGHBORS-RESET"))
			alg = new PriorityReembedding(
					OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN,
					PriorityReembedding.PRIORITY_NUM_OF_NEIGHBORS, true);

		else if (name.equalsIgnoreCase("PRIORITYREEMBED-BL-NUMNEIGHBORS-NORESET"))
			alg = new PriorityReembedding(
					OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN,
					PriorityReembedding.PRIORITY_NUM_OF_NEIGHBORS, false);

		else if (name.equalsIgnoreCase("PRIORITYREEMBED-BL-BALANCED-RESET"))
			alg = new PriorityReembedding(
					OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN,
					PriorityReembedding.PRIORITY_BALANCED, true);

		else if (name.equalsIgnoreCase("PRIORITYREEMBED-BL-BALANCED-NORESET"))
			alg = new PriorityReembedding(
					OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN,
					PriorityReembedding.PRIORITY_BALANCED, false);

		else if (name.equalsIgnoreCase("PRIORITYREEMBED-CI-NUMEMBED-RESET"))
			alg = new PriorityReembedding(
					OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN,
					PriorityReembedding.PRIORITY_NUM_OF_EMBEDDINGS, true);

		else if (name.equalsIgnoreCase("PRIORITYREEMBED-CI-NUMEMBED-NORESET"))
			alg = new PriorityReembedding(
					OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN,
					PriorityReembedding.PRIORITY_NUM_OF_EMBEDDINGS, false);

		else if (name.equalsIgnoreCase("PRIORITYREEMBED-CI-NUMNEIGHBORS-RESET"))
			alg = new PriorityReembedding(
					OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN,
					PriorityReembedding.PRIORITY_NUM_OF_NEIGHBORS, true);

		else if (name.equalsIgnoreCase("PRIORITYREEMBED-CI-NUMNEIGHBORS-NORESET"))
			alg = new PriorityReembedding(
					OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN,
					PriorityReembedding.PRIORITY_NUM_OF_NEIGHBORS, false);

		else if (name.equalsIgnoreCase("PRIORITYREEMBED-CI-BALANCED-RESET"))
			alg = new PriorityReembedding(
					OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN,
					PriorityReembedding.PRIORITY_BALANCED, true);

		else if (name.equalsIgnoreCase("PRIORITYREEMBED-CI-BALANCED-NORESET"))
			alg = new PriorityReembedding(
					OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN,
					PriorityReembedding.PRIORITY_BALANCED, false);


		else
			throw new IllegalArgumentException
				("Unknown layout algorithm '" + name + "'");
		
		return alg;
	}

	private static void usage ()
	{
		System.err.println (
	"--------------------------\n" +
	"ArrayOpt Microarray Design\n" +
	"--------------------------\n\n" +
	"Usage: ArrayOpt (affy | simple) (fix | nofix) <input> <rows> <columns>\n" +
	"          <probes> <length> <dep-seq> (check | no-check)\n" +
	"          (calc-bl | calc-ci | no-calc) (print | no-print) <alg>*\n");
		System.err.println (
	"where: 'affy'      indicates an Affymetrix chip type\n" +
	"       'simple'    indicates a simple chip type\n" +
	"       'fix'       considers fixed spots\n" +
	"       'no-fix'    ignores fixed spots in the input\n" +
	"       <input>     is a file name or RANDOM for a randomly generated\n" +
	"                      chip\n" +
	"       <rows>      is the number of rows in the chip\n" +
	"       <columns>   is the number of columns in the chip\n" +
	"       <probes>    is the number of probe in the chip\n" +
	"       <length>    is the length of probes sequences\n" +
	"       <dep-seq>   is the deposition sequence or\n" +
	"                      AFFY for Affymetrix's sequence or\n" +
	"                      SYNC for a 100-step ACGT repetition\n" +
	"       'check'     performs validation checks\n" +
	"       'no-check'  does not perform validation checks\n" +
	"       'calc-bl'   computes and prints total border length\n" +
	"       'calc-ci'   computes and prints average conflict indices\n" +
	"       'no-calc'   does not compute/print any quality measure\n" +
	"       'print'     prints the resulting layout on standard output\n" +
	"       'no-print'  does not print the resulting layout\n" +
	"       <alg>*      zero or more layout algorithms\n");
	
	}
	
	/*

		// ******************************************************
		// placement algorithm selection
		// ******************************************************

		// Greedy Embeddings Placer

		// TODO delete this
		else if (placement.equalsIgnoreCase("GREEDYEMB-BL-5-S"))
			placer = new GreedyEmbeddingsPlacer(GreedyEmbeddingsPlacer.BORDER_LENGTH_MIN, 5, GreedyEmbeddingsPlacer.SORT_PROBES);

		else if (placement.equalsIgnoreCase("GREEDYEMB-BL-50-S"))
			placer = new GreedyEmbeddingsPlacer(GreedyEmbeddingsPlacer.BORDER_LENGTH_MIN, 50, GreedyEmbeddingsPlacer.SORT_PROBES);

		else if (placement.equalsIgnoreCase("GREEDYEMB-BL-100-S"))
			placer = new GreedyEmbeddingsPlacer(GreedyEmbeddingsPlacer.BORDER_LENGTH_MIN, 100, GreedyEmbeddingsPlacer.SORT_PROBES);

		else if (placement.equalsIgnoreCase("GREEDYEMB-BL-200-S"))
			placer = new GreedyEmbeddingsPlacer(GreedyEmbeddingsPlacer.BORDER_LENGTH_MIN, 200, GreedyEmbeddingsPlacer.SORT_PROBES);

		else if (placement.equalsIgnoreCase("GREEDYEMB-BL-250-S"))
			placer = new GreedyEmbeddingsPlacer(GreedyEmbeddingsPlacer.BORDER_LENGTH_MIN, 250, GreedyEmbeddingsPlacer.SORT_PROBES);

		else if (placement.equalsIgnoreCase("GREEDYEMB-BL-350-S"))
			placer = new GreedyEmbeddingsPlacer(GreedyEmbeddingsPlacer.BORDER_LENGTH_MIN, 350, GreedyEmbeddingsPlacer.SORT_PROBES);

		else if (placement.equalsIgnoreCase("GREEDYEMB-BL-400-S"))
			placer = new GreedyEmbeddingsPlacer(GreedyEmbeddingsPlacer.BORDER_LENGTH_MIN, 400, GreedyEmbeddingsPlacer.SORT_PROBES);

		else if (placement.equalsIgnoreCase("GREEDYEMB-BL-500-S"))
			placer = new GreedyEmbeddingsPlacer(GreedyEmbeddingsPlacer.BORDER_LENGTH_MIN, 500, GreedyEmbeddingsPlacer.SORT_PROBES);

		else if (placement.equalsIgnoreCase("GREEDYEMB-BL-700-S"))
			placer = new GreedyEmbeddingsPlacer(GreedyEmbeddingsPlacer.BORDER_LENGTH_MIN, 700, GreedyEmbeddingsPlacer.SORT_PROBES);

		else if (placement.equalsIgnoreCase("GREEDYEMB-BL-800-S"))
			placer = new GreedyEmbeddingsPlacer(GreedyEmbeddingsPlacer.BORDER_LENGTH_MIN, 800, GreedyEmbeddingsPlacer.SORT_PROBES);

		else if (placement.equalsIgnoreCase("GREEDYEMB-BL-1K-S"))
			placer = new GreedyEmbeddingsPlacer(GreedyEmbeddingsPlacer.BORDER_LENGTH_MIN, 1000, GreedyEmbeddingsPlacer.SORT_PROBES);
		
		else if (placement.equalsIgnoreCase("GREEDYEMB-BL-10K"))
			placer = new GreedyEmbeddingsPlacer(GreedyEmbeddingsPlacer.BORDER_LENGTH_MIN, 10000);

		else if (placement.equalsIgnoreCase("GREEDYEMB-BL-5K"))
			placer = new GreedyEmbeddingsPlacer(GreedyEmbeddingsPlacer.BORDER_LENGTH_MIN, 5000);

		else if (placement.equalsIgnoreCase("GREEDYEMB-BL-2K"))
			placer = new GreedyEmbeddingsPlacer(GreedyEmbeddingsPlacer.BORDER_LENGTH_MIN, 2000);

		else if (placement.equalsIgnoreCase("GREEDYEMB-BL-1K"))
			placer = new GreedyEmbeddingsPlacer(GreedyEmbeddingsPlacer.BORDER_LENGTH_MIN, 1000);

		else if (placement.equalsIgnoreCase("GREEDYEMB-CI-50-S"))
			placer = new GreedyEmbeddingsPlacer(GreedyEmbeddingsPlacer.CONFLICT_INDEX_MIN, 50, GreedyEmbeddingsPlacer.SORT_PROBES);

		else if (placement.equalsIgnoreCase("GREEDYEMB-CI-100-S"))
			placer = new GreedyEmbeddingsPlacer(GreedyEmbeddingsPlacer.CONFLICT_INDEX_MIN, 100, GreedyEmbeddingsPlacer.SORT_PROBES);

		else if (placement.equalsIgnoreCase("GREEDYEMB-CI-200-S"))
			placer = new GreedyEmbeddingsPlacer(GreedyEmbeddingsPlacer.CONFLICT_INDEX_MIN, 200, GreedyEmbeddingsPlacer.SORT_PROBES);

		else if (placement.equalsIgnoreCase("GREEDYEMB-CI-300-S"))
			placer = new GreedyEmbeddingsPlacer(GreedyEmbeddingsPlacer.CONFLICT_INDEX_MIN, 300, GreedyEmbeddingsPlacer.SORT_PROBES);

		else if (placement.equalsIgnoreCase("GREEDYEMB-CI-100"))
			placer = new GreedyEmbeddingsPlacer(GreedyEmbeddingsPlacer.CONFLICT_INDEX_MIN, 100);


		else if (placement.equalsIgnoreCase("2D1-LEFT-GREEDY"))
			placer = new TwoDimensionalPartitioning(new GreedyPlacer(), 1);

		else if (placement.equalsIgnoreCase("2D4-LEFT-GREEDY"))
			placer = new TwoDimensionalPartitioning(new GreedyPlacer(), 4);

		else if (placement.equalsIgnoreCase("2D8-LEFT-GREEDY"))
			placer = new TwoDimensionalPartitioning(new GreedyPlacer(), 8);

		else if (placement.equalsIgnoreCase("2D16-LEFT-GREEDY"))
			placer = new TwoDimensionalPartitioning(new GreedyPlacer(), 16);

		else if (placement.equalsIgnoreCase("2D32-LEFT-GREEDY"))
			placer = new TwoDimensionalPartitioning(new GreedyPlacer(), 32);

		else if (placement.equalsIgnoreCase("2D64-LEFT-GREEDY"))
			placer = new TwoDimensionalPartitioning(new GreedyPlacer(), 64);

		else if (placement.equalsIgnoreCase("2D1-CENTER-GREEDY"))
			placer = new TwoDimensionalCenteredPartitioning(new GreedyPlacer(), 1);

		else if (placement.equalsIgnoreCase("2D4-CENTER-GREEDY"))
			placer = new TwoDimensionalCenteredPartitioning(new GreedyPlacer(), 4);

		else if (placement.equalsIgnoreCase("2D8-CENTER-GREEDY"))
			placer = new TwoDimensionalCenteredPartitioning(new GreedyPlacer(), 8);

		else if (placement.equalsIgnoreCase("2D16-CENTER-GREEDY"))
			placer = new TwoDimensionalCenteredPartitioning(new GreedyPlacer(), 16);

		else if (placement.equalsIgnoreCase("2D32-CENTER-GREEDY"))
			placer = new TwoDimensionalCenteredPartitioning(new GreedyPlacer(), 32);

		else if (placement.equalsIgnoreCase("2D64-CENTER-GREEDY"))
			placer = new TwoDimensionalCenteredPartitioning(new GreedyPlacer(), 64);

		else if (placement.equalsIgnoreCase("2D8-GRASPD"))
			placer = new TwoDimensionalPartitioning(
						new QuadraticAssignmentPlacer(
							new GraspDense(), QuadraticAssignmentPlacer.MODE_CONFLICT_INDEX),
						8);


		else if (placement.equalsIgnoreCase("PIVOTPART"))
			placer = new PivotPartitioning(new SequentialPlacer());

		else if (placement.equalsIgnoreCase("PIVOTPART1+REPTX20K-BL"))
			placer = new PivotPartitioning(new RowEpitaxial(20000),
							PivotPartitioning.MODE_BORDER_LENGTH,
							1);

		else if (placement.equalsIgnoreCase("PIVOTPART2+REPTX20K-BL"))
			placer = new PivotPartitioning(new RowEpitaxial(20000),
							PivotPartitioning.MODE_BORDER_LENGTH, 2);

		else if (placement.equalsIgnoreCase("PIVOTPART3+REPTX20K-BL"))
			placer = new PivotPartitioning(new RowEpitaxial(20000),
							PivotPartitioning.MODE_BORDER_LENGTH, 3);

		else if (placement.equalsIgnoreCase("PIVOTPART4+REPTX20K-BL"))
			placer = new PivotPartitioning(new RowEpitaxial(20000),
							PivotPartitioning.MODE_BORDER_LENGTH, 4);

		else if (placement.equalsIgnoreCase("PIVOTPART5+REPTX20K-BL"))
			placer = new PivotPartitioning(new RowEpitaxial(20000),
							PivotPartitioning.MODE_BORDER_LENGTH, 5);

		else if (placement.equalsIgnoreCase("PIVOTPART6+REPTX20K-BL"))
			placer = new PivotPartitioning(new RowEpitaxial(20000),
							PivotPartitioning.MODE_BORDER_LENGTH, 6);
		
		// *************************

		else if (placement.equalsIgnoreCase("PIVOTPART2+GREEDY20K-BL"))
			placer = new PivotPartitioning(new GreedyPlacer(
							GreedyPlacer.MODE_BORDER_LENGTH, 20000),
							PivotPartitioning.MODE_BORDER_LENGTH, 2);

		else if (placement.equalsIgnoreCase("PIVOTPART4+GREEDY20K-BL"))
			placer = new PivotPartitioning(new GreedyPlacer(
							GreedyPlacer.MODE_BORDER_LENGTH, 20000),
							PivotPartitioning.MODE_BORDER_LENGTH, 4);

		else if (placement.equalsIgnoreCase("PIVOTPART6+GREEDY20K-BL"))
			placer = new PivotPartitioning(new GreedyPlacer(
							GreedyPlacer.MODE_BORDER_LENGTH, 20000),
							PivotPartitioning.MODE_BORDER_LENGTH, 6);

		// **************************
		
		else if (placement.equalsIgnoreCase("PIVOTPART2+GREEDY2K-BL"))
			placer = new PivotPartitioning(new GreedyPlacer(
							GreedyPlacer.MODE_BORDER_LENGTH, 2000),
							PivotPartitioning.MODE_BORDER_LENGTH, 2);

		else if (placement.equalsIgnoreCase("PIVOTPART4+GREEDY2K-BL"))
			placer = new PivotPartitioning(new GreedyPlacer(
							GreedyPlacer.MODE_BORDER_LENGTH, 2000),
							PivotPartitioning.MODE_BORDER_LENGTH, 4);

		else if (placement.equalsIgnoreCase("PIVOTPART6+GREEDY2K-BL"))
			placer = new PivotPartitioning(new GreedyPlacer(
							GreedyPlacer.MODE_BORDER_LENGTH, 2000),
							PivotPartitioning.MODE_BORDER_LENGTH, 6);

		else if (placement.equalsIgnoreCase("PIVOTPART8+GREEDY2K-BL"))
			placer = new PivotPartitioning(new GreedyPlacer(
							GreedyPlacer.MODE_BORDER_LENGTH, 2000),
							PivotPartitioning.MODE_BORDER_LENGTH, 8);

		else if (placement.equalsIgnoreCase("PIVOTPART2+GREEDY2K-CI"))
			placer = new PivotPartitioning(new GreedyPlacer(
							GreedyPlacer.MODE_CONFLICT_INDEX, 2000),
							PivotPartitioning.MODE_CONFLICT_INDEX, 2);

		else if (placement.equalsIgnoreCase("PIVOTPART4+GREEDY2K-CI"))
			placer = new PivotPartitioning(new GreedyPlacer(
							GreedyPlacer.MODE_CONFLICT_INDEX, 2000),
							PivotPartitioning.MODE_CONFLICT_INDEX, 4);

		else if (placement.equalsIgnoreCase("PIVOTPART6+GREEDY2K-CI"))
			placer = new PivotPartitioning(new GreedyPlacer(
							GreedyPlacer.MODE_CONFLICT_INDEX, 2000),
							PivotPartitioning.MODE_CONFLICT_INDEX, 6);

		else if (placement.equalsIgnoreCase("PIVOTPART8+GREEDY2K-CI"))
			placer = new PivotPartitioning(new GreedyPlacer(
							GreedyPlacer.MODE_CONFLICT_INDEX, 2000),
							PivotPartitioning.MODE_CONFLICT_INDEX, 8);
		
		// **************************

		else if (placement.equalsIgnoreCase("PIVOTPART2+GREEDY1K-BL"))
			placer = new PivotPartitioning(new GreedyPlacer(
							GreedyPlacer.MODE_BORDER_LENGTH, 1000),
							PivotPartitioning.MODE_BORDER_LENGTH, 2);

		else if (placement.equalsIgnoreCase("PIVOTPART4+GREEDY1K-BL"))
			placer = new PivotPartitioning(new GreedyPlacer(
							GreedyPlacer.MODE_BORDER_LENGTH, 1000),
							PivotPartitioning.MODE_BORDER_LENGTH, 4);

		else if (placement.equalsIgnoreCase("PIVOTPART6+GREEDY1K-BL"))
			placer = new PivotPartitioning(new GreedyPlacer(
							GreedyPlacer.MODE_BORDER_LENGTH, 1000),
							PivotPartitioning.MODE_BORDER_LENGTH, 6);

		else if (placement.equalsIgnoreCase("PIVOTPART8+GREEDY1K-BL"))
			placer = new PivotPartitioning(new GreedyPlacer(
							GreedyPlacer.MODE_BORDER_LENGTH, 1000),
							PivotPartitioning.MODE_BORDER_LENGTH, 8);

		else if (placement.equalsIgnoreCase("PIVOTPART2+GREEDY1K-CI"))
			placer = new PivotPartitioning(new GreedyPlacer(
							GreedyPlacer.MODE_CONFLICT_INDEX, 1000),
							PivotPartitioning.MODE_CONFLICT_INDEX, 2);

		else if (placement.equalsIgnoreCase("PIVOTPART4+GREEDY1K-CI"))
			placer = new PivotPartitioning(new GreedyPlacer(
							GreedyPlacer.MODE_CONFLICT_INDEX, 1000),
							PivotPartitioning.MODE_CONFLICT_INDEX, 4);

		else if (placement.equalsIgnoreCase("PIVOTPART6+GREEDY1K-CI"))
			placer = new PivotPartitioning(new GreedyPlacer(
							GreedyPlacer.MODE_CONFLICT_INDEX, 1000),
							PivotPartitioning.MODE_CONFLICT_INDEX, 6);

		else if (placement.equalsIgnoreCase("PIVOTPART8+GREEDY1K-CI"))
			placer = new PivotPartitioning(new GreedyPlacer(
							GreedyPlacer.MODE_CONFLICT_INDEX, 1000),
							PivotPartitioning.MODE_CONFLICT_INDEX, 8);
		
		// **************************

		else if (placement.equalsIgnoreCase("PIVOTPART2+GREEDY100-BL"))
			placer = new PivotPartitioning(new GreedyPlacer(
							GreedyPlacer.MODE_BORDER_LENGTH, 100),
							PivotPartitioning.MODE_BORDER_LENGTH, 2);

		else if (placement.equalsIgnoreCase("PIVOTPART4+GREEDY100-BL"))
			placer = new PivotPartitioning(new GreedyPlacer(
							GreedyPlacer.MODE_BORDER_LENGTH, 100),
							PivotPartitioning.MODE_BORDER_LENGTH, 4);

		else if (placement.equalsIgnoreCase("PIVOTPART6+GREEDY100-BL"))
			placer = new PivotPartitioning(new GreedyPlacer(
							GreedyPlacer.MODE_BORDER_LENGTH, 100),
							PivotPartitioning.MODE_BORDER_LENGTH, 6);

		else if (placement.equalsIgnoreCase("PIVOTPART8+GREEDY100-BL"))
			placer = new PivotPartitioning(new GreedyPlacer(
							GreedyPlacer.MODE_BORDER_LENGTH, 100),
							PivotPartitioning.MODE_BORDER_LENGTH, 8);

		else if (placement.equalsIgnoreCase("PIVOTPART2+GREEDY100-CI"))
			placer = new PivotPartitioning(new GreedyPlacer(
							GreedyPlacer.MODE_CONFLICT_INDEX, 100),
							PivotPartitioning.MODE_CONFLICT_INDEX, 2);

		else if (placement.equalsIgnoreCase("PIVOTPART4+GREEDY100-CI"))
			placer = new PivotPartitioning(new GreedyPlacer(
							GreedyPlacer.MODE_CONFLICT_INDEX, 100),
							PivotPartitioning.MODE_CONFLICT_INDEX, 4);

		else if (placement.equalsIgnoreCase("PIVOTPART6+GREEDY100-CI"))
			placer = new PivotPartitioning(new GreedyPlacer(
							GreedyPlacer.MODE_CONFLICT_INDEX, 100),
							PivotPartitioning.MODE_CONFLICT_INDEX, 6);

		else if (placement.equalsIgnoreCase("PIVOTPART8+GREEDY100-CI"))
			placer = new PivotPartitioning(new GreedyPlacer(
							GreedyPlacer.MODE_CONFLICT_INDEX, 100),
							PivotPartitioning.MODE_CONFLICT_INDEX, 8);
		
		// **************************

		else if (placement.equalsIgnoreCase("CLUSTER-SEQ-BL"))
			placer = new ClusterPartitioning(new SequentialPlacer(),
							PivotPartitioning.MODE_BORDER_LENGTH);

		else if (placement.equalsIgnoreCase("CLUSTER-SEQ-CI"))
			placer = new ClusterPartitioning(new SequentialPlacer(),
							PivotPartitioning.MODE_CONFLICT_INDEX);

		// **************************

		else if (placement.equalsIgnoreCase("PIVOTPART8+GREEDY50-CI"))
			placer = new PivotPartitioning(new GreedyPlacer(
							GreedyPlacer.MODE_CONFLICT_INDEX, 50),
							PivotPartitioning.MODE_CONFLICT_INDEX, 8);

		
		// ****************** GRASP Sparse ******************
		
		else if (placement.equalsIgnoreCase("GRASPS-BL"))
			placer = new QuadraticAssignmentPlacer (
							new GraspSparse(), QuadraticAssignmentPlacer.MODE_BORDER_LENGTH
						);

		else if (placement.equalsIgnoreCase("GRASPS-INV-BL"))
		{
			QAPSolverAlgorithm qapsolver = new GraspSparse();
			qapsolver.setSwapMatrices(true);
			placer = new QuadraticAssignmentPlacer (
							qapsolver, QuadraticAssignmentPlacer.MODE_BORDER_LENGTH
						);
		}

		else if (placement.equalsIgnoreCase("GRASPS-CI"))
			placer = new QuadraticAssignmentPlacer (
							new GraspSparse(), QuadraticAssignmentPlacer.MODE_CONFLICT_INDEX
						);

		else if (placement.equalsIgnoreCase("GRASPS-INV-CI"))
		{
			QAPSolverAlgorithm qapsolver = new GraspSparse();
			qapsolver.setSwapMatrices(true);
			placer = new QuadraticAssignmentPlacer (
							qapsolver, QuadraticAssignmentPlacer.MODE_CONFLICT_INDEX
						);
		}

		// ****************** GRASP Dense ******************
		
		else if (placement.equalsIgnoreCase("GRASPD-BL"))
			placer = new QuadraticAssignmentPlacer (
							new GraspDense(), QuadraticAssignmentPlacer.MODE_BORDER_LENGTH
						);

		else if (placement.equalsIgnoreCase("GRASPD-INV-BL"))
		{
			QAPSolverAlgorithm qapsolver = new GraspDense();
			qapsolver.setSwapMatrices(true);
			placer = new QuadraticAssignmentPlacer (
							qapsolver, QuadraticAssignmentPlacer.MODE_BORDER_LENGTH
						);
		}

		else if (placement.equalsIgnoreCase("GRASPD-CI"))
			placer = new QuadraticAssignmentPlacer (
							new GraspDense(), QuadraticAssignmentPlacer.MODE_CONFLICT_INDEX
						);

		else if (placement.equalsIgnoreCase("GRASPD-INV-CI"))
		{
			QAPSolverAlgorithm qapsolver = new GraspDense();
			qapsolver.setSwapMatrices(true);
			placer = new QuadraticAssignmentPlacer (
							qapsolver, QuadraticAssignmentPlacer.MODE_CONFLICT_INDEX
						);
		}
		
		else
			throw new IllegalArgumentException
				("Unknown placement algorithm: '" + placement + "'");

		// ******************************************************
		// post-placement algorithm selection
		// ******************************************************
		
		if (optimization.equalsIgnoreCase("NONE"))
			optimizer = null;

		else if (optimization.equalsIgnoreCase("PRIORITY-BL-E-INC"))
			optimizer = new PriorityReembedding(
					OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN,
					PriorityReembedding.PRIORITY_NUM_OF_EMBEDDINGS,
					true);

		else if (optimization.equalsIgnoreCase("PRIORITY-BL-E-TOT"))
			optimizer = new PriorityReembedding(
					OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN,
					PriorityReembedding.PRIORITY_NUM_OF_EMBEDDINGS,
					false);

		else if (optimization.equalsIgnoreCase("PRIORITY-BL-N-INC"))
			optimizer = new PriorityReembedding(
					OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN,
					PriorityReembedding.PRIORITY_NUM_OF_NEIGHBORS,
					true);

		else if (optimization.equalsIgnoreCase("PRIORITY-BL-N-TOT"))
			optimizer = new PriorityReembedding(
					OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN,
					PriorityReembedding.PRIORITY_NUM_OF_NEIGHBORS,
					false);

		else if (optimization.equalsIgnoreCase("PRIORITY-BL-B-INC"))
			optimizer = new PriorityReembedding(
					OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN,
					PriorityReembedding.PRIORITY_BALANCED,
					true);

		else if (optimization.equalsIgnoreCase("PRIORITY-BL-B-TOT"))
			optimizer = new PriorityReembedding(
					OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN,
					PriorityReembedding.PRIORITY_BALANCED,
					false);
		
		else if (optimization.equalsIgnoreCase("PRIORITY-CI-E-INC"))
			optimizer = new PriorityReembedding(
					OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN,
					PriorityReembedding.PRIORITY_NUM_OF_EMBEDDINGS,
					true);

		else if (optimization.equalsIgnoreCase("PRIORITY-CI-E-TOT"))
			optimizer = new PriorityReembedding(
					OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN,
					PriorityReembedding.PRIORITY_NUM_OF_EMBEDDINGS,
					false);

		else if (optimization.equalsIgnoreCase("PRIORITY-CI-N-INC"))
			optimizer = new PriorityReembedding(
					OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN,
					PriorityReembedding.PRIORITY_NUM_OF_NEIGHBORS,
					true);

		else if (optimization.equalsIgnoreCase("PRIORITY-CI-N-TOT"))
			optimizer = new PriorityReembedding(
					OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN,
					PriorityReembedding.PRIORITY_NUM_OF_NEIGHBORS,
					false);

		else if (optimization.equalsIgnoreCase("PRIORITY-CI-B-INC"))
			optimizer = new PriorityReembedding(
					OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN,
					PriorityReembedding.PRIORITY_BALANCED,
					true);

		else if (optimization.equalsIgnoreCase("PRIORITY-CI-B-TOT"))
			optimizer = new PriorityReembedding(
					OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN,
					PriorityReembedding.PRIORITY_BALANCED,
					false);

		else if (optimization.equalsIgnoreCase("SEQREEMBED-BL-RESET"))
			optimizer = new SequentialReembedding(
							OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN, true);
		
		else if (optimization.equalsIgnoreCase("SEQREEMBED-CI-RESET"))
			optimizer = new SequentialReembedding(
							OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN, true);

		else if (optimization.equalsIgnoreCase("SEQREEMBED-BL-NORESET"))
			optimizer = new SequentialReembedding(
							OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN, false);
		
		else if (optimization.equalsIgnoreCase("SEQREEMBED-CI-NORESET"))
			optimizer = new SequentialReembedding(
							OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN, false, 0.005);

		// ****************** GRASP Path-relinking ******************

		else if (optimization.equalsIgnoreCase("GRASPPR32-BL"))
		{
			optimizer = new QuadraticAssignmentPlacer (
							new GraspPathRelinking(1,32), QuadraticAssignmentPlacer.MODE_BORDER_LENGTH
						);
		}

		else if (optimization.equalsIgnoreCase("GRASPPR32-INV-BL"))
		{
			QAPSolverAlgorithm qapsolver = new GraspPathRelinking(1,32);
			qapsolver.setSwapMatrices(true);
			optimizer = new QuadraticAssignmentPlacer (
							qapsolver, QuadraticAssignmentPlacer.MODE_BORDER_LENGTH
						);
		}

		else if (optimization.equalsIgnoreCase("GRASPPR32-CI"))
		{
			optimizer = new QuadraticAssignmentPlacer (
							new GraspPathRelinking(1,32), QuadraticAssignmentPlacer.MODE_CONFLICT_INDEX
						);
		}

		else if (optimization.equalsIgnoreCase("GRASPPR32-INV-CI"))
		{
			QAPSolverAlgorithm qapsolver = new GraspPathRelinking(1,32);
			qapsolver.setSwapMatrices(true);
			optimizer = new QuadraticAssignmentPlacer (
							qapsolver, QuadraticAssignmentPlacer.MODE_CONFLICT_INDEX
						);
		}

		// ****************** GRASP 128 Path-relinking ******************

		else if (optimization.equalsIgnoreCase("GRASPPR128-BL"))
		{
			optimizer = new QuadraticAssignmentPlacer (
							new GraspPathRelinking(1, 128), QuadraticAssignmentPlacer.MODE_BORDER_LENGTH
						);
		}

		else if (optimization.equalsIgnoreCase("GRASPPR128-INV-BL"))
		{
			QAPSolverAlgorithm qapsolver = new GraspPathRelinking(1, 128);
			qapsolver.setSwapMatrices(true);
			optimizer = new QuadraticAssignmentPlacer (
							qapsolver, QuadraticAssignmentPlacer.MODE_BORDER_LENGTH
						);
		}

		else if (optimization.equalsIgnoreCase("GRASPPR128-CI"))
		{
			optimizer = new QuadraticAssignmentPlacer (
							new GraspPathRelinking(1, 128), QuadraticAssignmentPlacer.MODE_CONFLICT_INDEX
						);
		}

		else if (optimization.equalsIgnoreCase("GRASPPR128-INV-CI"))
		{
			QAPSolverAlgorithm qapsolver = new GraspPathRelinking(1, 128);
			qapsolver.setSwapMatrices(true);
			optimizer = new QuadraticAssignmentPlacer (
							qapsolver, QuadraticAssignmentPlacer.MODE_CONFLICT_INDEX
						);
		}

		else
			throw new IllegalArgumentException
				("Unknown post-placement algorithm '" + optimization + "'");

	 */
}
