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

	public static void main (String args[])
	{
		Chip					chip, copy;
		PlacementAlgorithm		placer = null;
		PostPlacementAlgorithm	optimizer = null;
		String					chip_type, filename, algorithm, dep_seq;
		int						rows, cols, probes, probe_len, unplaced;
		boolean					ignore_fixed;
		
		try
		{
			// get mandatory command-line arguments
			chip_type = args[0];
			filename  = args[2];
			rows      = Integer.parseInt(args[3]);
			cols      = Integer.parseInt(args[4]);
			probes    = Integer.parseInt(args[5]);
			probe_len = Integer.parseInt(args[6]);
			dep_seq	  = args[7];
			algorithm = args[8];
		}
		catch (Exception e)
		{
			usage();
			System.exit(1);
			return;
		}

		// ignore fixed spots?
		if (args[1].equalsIgnoreCase("fix"))
		{
			ignore_fixed = false;
		}
		else if (args[1].equalsIgnoreCase("nofix"))
		{
			ignore_fixed = true;
		}
		else
		{
			usage ();
			System.exit(1);
			return;
		}

		// standard deposition sequences
		if (dep_seq.equalsIgnoreCase("AFFY"))
			dep_seq = AFFY_DEP_SEQ;
		else if (dep_seq.equalsIgnoreCase("SYNC"))
			dep_seq = SYNC_DEP_SEQ;		
		
		try
		{
			// algorithm selection

			if (algorithm.equalsIgnoreCase("RANDOM"))
				placer = new RandomFiller();

			else if (algorithm.equalsIgnoreCase("SEQUENTIAL"))
				placer = new SequentialFiller();

			else if (algorithm.equalsIgnoreCase("GREEDY100-BL"))
				placer = new GreedyFiller(GreedyFiller.MODE_BORDER_LENGTH,
											100, false);
			
			else if (algorithm.equalsIgnoreCase("GREEDY20K-BL"))
				placer = new GreedyFiller(GreedyFiller.MODE_BORDER_LENGTH,
											20000, false);

			else if (algorithm.equalsIgnoreCase("REPTX20K-BL"))
				placer = new RowEpitaxial(20000, false);

			else if (algorithm.equalsIgnoreCase("2D1-SEQ"))
				placer = new TwoDimensionalPartitioning(new SequentialFiller(), 1);

			else if (algorithm.equalsIgnoreCase("2D4-SEQ"))
				placer = new TwoDimensionalPartitioning(new SequentialFiller(), 4);

			else if (algorithm.equalsIgnoreCase("2D1-GREEDY"))
				placer = new TwoDimensionalPartitioning(new GreedyFiller(), 1);

			else if (algorithm.equalsIgnoreCase("2D4-GREEDY"))
				placer = new TwoDimensionalPartitioning(new GreedyFiller(), 4);

			else if (algorithm.equalsIgnoreCase("2D8-GREEDY"))
				placer = new TwoDimensionalPartitioning(new GreedyFiller(), 8);

			else if (algorithm.equalsIgnoreCase("2D8-GRASPD"))
				placer = new TwoDimensionalPartitioning(
							new QAPOptimization(
								new GraspDense(), QAPOptimization.MODE_CONFLICT_INDEX),
							8);

			else if (algorithm.equalsIgnoreCase("LEFT"))
				placer = new ProbeSetEmbeddingWrapper(new LeftMostEmbedding());
				
			else if (algorithm.equalsIgnoreCase("RIGHT"))
				placer = new ProbeSetEmbeddingWrapper(new RightMostEmbedding());
				
			else if (algorithm.equalsIgnoreCase("CENTER"))
				placer = new ProbeSetEmbeddingWrapper(new CenteredEmbedding());
				
			else if (algorithm.equalsIgnoreCase("PIVOT"))
				placer = new ProbeSetEmbeddingWrapper(new PivotEmbedding());

			else if (algorithm.equalsIgnoreCase("PIVOTPART"))
				placer = new PivotPartitioning(new SequentialFiller());

			else if (algorithm.equalsIgnoreCase("PIVOTPART1+REPTX20K-BL"))
				placer = new PivotPartitioning(new RowEpitaxial(20000), 1);

			else if (algorithm.equalsIgnoreCase("PIVOTPART2+REPTX20K-BL"))
				placer = new PivotPartitioning(new RowEpitaxial(20000), 2);

			else if (algorithm.equalsIgnoreCase("PIVOTPART3+REPTX20K-BL"))
				placer = new PivotPartitioning(new RowEpitaxial(20000), 3);

			else if (algorithm.equalsIgnoreCase("PIVOTPART4+REPTX20K-BL"))
				placer = new PivotPartitioning(new RowEpitaxial(20000), 4);

			else if (algorithm.equalsIgnoreCase("PRIORITY-BL-E-INC"))
				optimizer = new PriorityReembedding(
						OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN,
						PriorityReembedding.PRIORITY_NUM_OF_EMBEDDINGS,
						true);

			else if (algorithm.equalsIgnoreCase("PRIORITY-BL-E-TOT"))
				optimizer = new PriorityReembedding(
						OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN,
						PriorityReembedding.PRIORITY_NUM_OF_EMBEDDINGS,
						false);

			else if (algorithm.equalsIgnoreCase("PRIORITY-BL-N-INC"))
				optimizer = new PriorityReembedding(
						OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN,
						PriorityReembedding.PRIORITY_NUM_OF_NEIGHBORS,
						true);

			else if (algorithm.equalsIgnoreCase("PRIORITY-BL-N-TOT"))
				optimizer = new PriorityReembedding(
						OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN,
						PriorityReembedding.PRIORITY_NUM_OF_NEIGHBORS,
						false);

			else if (algorithm.equalsIgnoreCase("PRIORITY-BL-B-INC"))
				optimizer = new PriorityReembedding(
						OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN,
						PriorityReembedding.PRIORITY_BALANCED,
						true);

			else if (algorithm.equalsIgnoreCase("PRIORITY-BL-B-TOT"))
				optimizer = new PriorityReembedding(
						OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN,
						PriorityReembedding.PRIORITY_BALANCED,
						false);
			
			else if (algorithm.equalsIgnoreCase("PRIORITY-CI-E-INC"))
				optimizer = new PriorityReembedding(
						OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN,
						PriorityReembedding.PRIORITY_NUM_OF_EMBEDDINGS,
						true);

			else if (algorithm.equalsIgnoreCase("PRIORITY-CI-E-TOT"))
				optimizer = new PriorityReembedding(
						OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN,
						PriorityReembedding.PRIORITY_NUM_OF_EMBEDDINGS,
						false);

			else if (algorithm.equalsIgnoreCase("PRIORITY-CI-N-INC"))
				optimizer = new PriorityReembedding(
						OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN,
						PriorityReembedding.PRIORITY_NUM_OF_NEIGHBORS,
						true);

			else if (algorithm.equalsIgnoreCase("PRIORITY-CI-N-TOT"))
				optimizer = new PriorityReembedding(
						OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN,
						PriorityReembedding.PRIORITY_NUM_OF_NEIGHBORS,
						false);

			else if (algorithm.equalsIgnoreCase("PRIORITY-CI-B-INC"))
				optimizer = new PriorityReembedding(
						OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN,
						PriorityReembedding.PRIORITY_BALANCED,
						true);

			else if (algorithm.equalsIgnoreCase("PRIORITY-CI-B-TOT"))
				optimizer = new PriorityReembedding(
						OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN,
						PriorityReembedding.PRIORITY_BALANCED,
						false);

			else if (algorithm.equalsIgnoreCase("SEQREEMBED-BL-INC"))
				optimizer = new SequentialReembedding(
								OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN, true);
			
			else if (algorithm.equalsIgnoreCase("SEQREEMBED-CI-INC"))
				optimizer = new SequentialReembedding(
								OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN, true);

			else if (algorithm.equalsIgnoreCase("SEQREEMBED-BL-TOT"))
				optimizer = new SequentialReembedding(
								OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN, false);
			
			else if (algorithm.equalsIgnoreCase("SEQREEMBED-CI-TOT"))
				optimizer = new SequentialReembedding(
								OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN, false);

			// ****************** GRASP Sparse ******************
			
			else if (algorithm.equalsIgnoreCase("GRASPS-BL"))
				placer = new QAPOptimization (
								new GraspSparse(), QAPOptimization.MODE_BORDER_LENGTH
							);

			else if (algorithm.equalsIgnoreCase("GRASPS-INV-BL"))
			{
				QAPSolverAlgorithm qapsolver = new GraspSparse();
				qapsolver.setSwapMatrices(true);
				placer = new QAPOptimization (
								qapsolver, QAPOptimization.MODE_BORDER_LENGTH
							);
			}

			else if (algorithm.equalsIgnoreCase("GRASPS-CI"))
				placer = new QAPOptimization (
								new GraspSparse(), QAPOptimization.MODE_CONFLICT_INDEX
							);

			else if (algorithm.equalsIgnoreCase("GRASPS-INV-CI"))
			{
				QAPSolverAlgorithm qapsolver = new GraspSparse();
				qapsolver.setSwapMatrices(true);
				placer = new QAPOptimization (
								qapsolver, QAPOptimization.MODE_CONFLICT_INDEX
							);
			}

			// ****************** GRASP Dense ******************
			
			else if (algorithm.equalsIgnoreCase("GRASPD-BL"))
				placer = new QAPOptimization (
								new GraspDense(), QAPOptimization.MODE_BORDER_LENGTH
							);

			else if (algorithm.equalsIgnoreCase("GRASPD-INV-BL"))
			{
				QAPSolverAlgorithm qapsolver = new GraspDense();
				qapsolver.setSwapMatrices(true);
				placer = new QAPOptimization (
								qapsolver, QAPOptimization.MODE_BORDER_LENGTH
							);
			}

			else if (algorithm.equalsIgnoreCase("GRASPD-CI"))
				placer = new QAPOptimization (
								new GraspDense(), QAPOptimization.MODE_CONFLICT_INDEX
							);

			else if (algorithm.equalsIgnoreCase("GRASPD-INV-CI"))
			{
				QAPSolverAlgorithm qapsolver = new GraspDense();
				qapsolver.setSwapMatrices(true);
				placer = new QAPOptimization (
								qapsolver, QAPOptimization.MODE_CONFLICT_INDEX
							);
			}
			
			// ****************** GRASP Path-relinking ******************

			else if (algorithm.equalsIgnoreCase("GRASPPR-BL"))
			{
				optimizer = new QAPOptimization (
								new GraspPathRelinking(), QAPOptimization.MODE_BORDER_LENGTH
							);
			}

			else if (algorithm.equalsIgnoreCase("GRASPPR-INV-BL"))
			{
				QAPSolverAlgorithm qapsolver = new GraspPathRelinking();
				qapsolver.setSwapMatrices(true);
				optimizer = new QAPOptimization (
								qapsolver, QAPOptimization.MODE_BORDER_LENGTH
							);
			}

			else if (algorithm.equalsIgnoreCase("GRASPPR-CI"))
			{
				optimizer = new QAPOptimization (
								new GraspPathRelinking(), QAPOptimization.MODE_CONFLICT_INDEX
							);
			}

			else if (algorithm.equalsIgnoreCase("GRASPPR-INV-CI"))
			{
				QAPSolverAlgorithm qapsolver = new GraspPathRelinking();
				qapsolver.setSwapMatrices(true);
				optimizer = new QAPOptimization (
								qapsolver, QAPOptimization.MODE_CONFLICT_INDEX
							);
			}

			// ****************** GRASP 128 Path-relinking ******************

			else if (algorithm.equalsIgnoreCase("GRASPPR128-BL"))
			{
				optimizer = new QAPOptimization (
								new GraspPathRelinking(1, 128), QAPOptimization.MODE_BORDER_LENGTH
							);
			}

			else if (algorithm.equalsIgnoreCase("GRASPPR128-INV-BL"))
			{
				QAPSolverAlgorithm qapsolver = new GraspPathRelinking(1, 128);
				qapsolver.setSwapMatrices(true);
				optimizer = new QAPOptimization (
								qapsolver, QAPOptimization.MODE_BORDER_LENGTH
							);
			}

			else if (algorithm.equalsIgnoreCase("GRASPPR128-CI"))
			{
				optimizer = new QAPOptimization (
								new GraspPathRelinking(1, 128), QAPOptimization.MODE_CONFLICT_INDEX
							);
			}

			else if (algorithm.equalsIgnoreCase("GRASPPR128-INV-CI"))
			{
				QAPSolverAlgorithm qapsolver = new GraspPathRelinking(1, 128);
				qapsolver.setSwapMatrices(true);
				optimizer = new QAPOptimization (
								qapsolver, QAPOptimization.MODE_CONFLICT_INDEX
							);
			}

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
		
		if (filename.equalsIgnoreCase("RANDOM"))
		{
			chip.createRandomLayout();
		}
		else
		{
			try
			{
				// create a file reader for the input
				FileReader file = new FileReader(filename);
	
				// read input
				chip.readLayout (file, ignore_fixed);
	
				// close reader
				file.close();
	
			}
			catch (Exception e)
			{
				e.printStackTrace ();
				System.exit(1);
				return;
			}
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
			long start = System.nanoTime();

			// re-place probes on the chip
			unplaced = placer.makeLayout(chip);
			
			long end = System.nanoTime();
			System.err.println("Elapsed time: " + (end - start)/Math.pow(10,9) + " sec");
			
			if (unplaced > 0)
				System.err.println("WARNING: " + unplaced + " unplaced probe(s)");
		}

		if (optimizer != null)
		{
			long start = System.nanoTime();
			
			optimizer.optimizeLayout (chip);

			long end = System.nanoTime();
			System.err.println("Elapsed time: " + (end - start)/Math.pow(10,9) + " sec");
		}

		// TODO remove this
		if (algorithm.endsWith("BL"))
			System.err.println("Border length after: " + LayoutEvaluation.borderLength(chip));
		else if (algorithm.endsWith("CI"))
			System.err.println("Average conflict index after: " + LayoutEvaluation.averageConflictIndex(chip));

		if (chip.equals(copy))
			System.err.println("WARNING: new layout is equal to the original specification.");
		else
			System.err.println("New layout is not equal to the original specification.");

		if (!chip.compatible(copy))
			System.err.println("WARNING: new layout is not compatible with the original specification.");
		else
			System.err.println("New layout is compatible with the original specification.");

		//*
		try
		{
			// print chip layout
			chip.writeLayout(new PrintWriter(System.out));
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
		//*/

		System.exit(0);
	}

	private static void usage ()
	{
		System.err.println (
			"Usage: ArrayOpt [affy | simple] [fix | nofix] <input> <rows> " +
				"<columns> <probes> <probe length> <dep seq> <algo>");
		System.err.println (
			"where: <input> is a file name or 'RANDOM' for a randomly " +
			"generated chip\n" +
			"       <dep seq> is a deposition sequence or AFFY for " +
			"Affymetrix's sequence or SYNC for a 100-step 'ACGT' repetition");
	}
}
