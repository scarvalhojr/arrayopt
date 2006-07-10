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
		PostPlacementAlgorithm	optimizer = null;
		PlacementAlgorithm		placer = null;
		Chip	chip, copy = null;
		String	type, filename, dep_seq, placement, optimization;
		int		rows, cols, probes, probe_len, unplaced;
		boolean	ignore_fixed, check, calc_bl, calc_ci, print_chip;
		long	start, end;

		try
		{
			// get command-line arguments
			type = args[0];
			
			// ignore fixed spots?
			if (args[1].equalsIgnoreCase("fix"))
				ignore_fixed = false;
			else if (args[1].equalsIgnoreCase("nofix"))
				ignore_fixed = true;
			else
				throw new IllegalArgumentException
					("Illegal argument: " + args[1]); 

			filename     = args[2];
			rows         = Integer.parseInt(args[3]);
			cols         = Integer.parseInt(args[4]);
			probes       = Integer.parseInt(args[5]);
			probe_len    = Integer.parseInt(args[6]);
			dep_seq	     = args[7];
			placement    = args[8];
			optimization = args[9];

			// perform validation checks?
			if (args[10].equalsIgnoreCase("check"))
				check = true;
			else if (args[10].equalsIgnoreCase("no-check"))
				check = false;
			else
				throw new IllegalArgumentException
					("Illegal argument '" + args[10] + "'"); 

			// compute and print total border
			// length or average conflict indices? 
			if (args[11].equalsIgnoreCase("calc-bl"))
			{
				calc_bl = true;
				calc_ci = false;
			}
			else if (args[11].equalsIgnoreCase("calc-ci"))
			{
				calc_bl = false;
				calc_ci = true;
			}
			else if (args[11].equalsIgnoreCase("no-calc"))
			{
				calc_bl = false;
				calc_ci = false;
			}
			else
				throw new IllegalArgumentException
					("Illegal argument '" + args[11] + "'"); 
			
			// print produced layout?
			if (args[12].equalsIgnoreCase("print"))
				print_chip = true;
			else if (args[12].equalsIgnoreCase("no-print"))
				print_chip = false;
			else
				throw new IllegalArgumentException
					("Illegal argument '" + args[12] + "'"); 
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
			System.err.println("ERROR: All arguments are mandatory");
			System.exit(1);
			return;
		}
		catch (IllegalArgumentException e)
		{
			usage();
			System.err.println("ERROR: " + e.getMessage());
			System.exit(1);
			return;
		}
		
		// standard deposition sequences
		if (dep_seq.equalsIgnoreCase("AFFY"))
			dep_seq = AFFY_DEP_SEQ;
		else if (dep_seq.equalsIgnoreCase("SYNC"))
			dep_seq = SYNC_DEP_SEQ;		

		// create chip
		if (type.equalsIgnoreCase("simple"))
		{
			chip = new SimpleChip (rows, cols, probes, probe_len, dep_seq);
		}
		else if (type.equalsIgnoreCase("affy"))
		{
			if (probe_len != AffymetrixChip.AFFY_PROBE_LENGTH)
				throw new IllegalArgumentException
					("Invalid probe length (Affymetrix probes must be " +
					 AffymetrixChip.AFFY_PROBE_LENGTH + " base-long).");

			chip = new AffymetrixChip (rows, cols, probes, dep_seq);
		}
		else
			throw new IllegalArgumentException
				("Unknown chip type '" + type + "'");

		// ******************************************************
		// placement algorithm selection
		// ******************************************************

		if (placement.equalsIgnoreCase("NONE"))
			placer = null;

		else if (placement.equalsIgnoreCase("RANDOM"))
			placer = new RandomFiller();

		else if (placement.equalsIgnoreCase("SEQUENTIAL"))
			placer = new SequentialFiller();

		else if (placement.equalsIgnoreCase("GREEDY-BL"))
			placer = new GreedyFiller(GreedyFiller.MODE_BORDER_LENGTH);

		else if (placement.equalsIgnoreCase("GREEDY-CI"))
			placer = new GreedyFiller(GreedyFiller.MODE_CONFLICT_INDEX);

		else if (placement.equalsIgnoreCase("GREEDY-BL-20K"))
			placer = new GreedyFiller(GreedyFiller.MODE_BORDER_LENGTH, 20000);

		else if (placement.equalsIgnoreCase("GREEDY-BL-10K"))
			placer = new GreedyFiller(GreedyFiller.MODE_BORDER_LENGTH, 10000);

		else if (placement.equalsIgnoreCase("GREEDY-BL-5K"))
			placer = new GreedyFiller(GreedyFiller.MODE_BORDER_LENGTH, 5000);

		else if (placement.equalsIgnoreCase("GREEDY-BL-2K"))
			placer = new GreedyFiller(GreedyFiller.MODE_BORDER_LENGTH, 2000);

		else if (placement.equalsIgnoreCase("GREEDY-BL-2K-S"))
			placer = new GreedyFiller(GreedyFiller.MODE_BORDER_LENGTH, 2000,
										GreedyFiller.SORT_EMBEDDINGS);

		else if (placement.equalsIgnoreCase("GREEDY-BL-2K-R"))
			placer = new GreedyFiller(GreedyFiller.MODE_BORDER_LENGTH, 2000,
										GreedyFiller.RANDOMIZE_INPUT);

		else if (placement.equalsIgnoreCase("GREEDY-BL-1K"))
			placer = new GreedyFiller(GreedyFiller.MODE_BORDER_LENGTH, 1000);

		else if (placement.equalsIgnoreCase("GREEDY-BL-100"))
			placer = new GreedyFiller(GreedyFiller.MODE_BORDER_LENGTH, 100);		

		else if (placement.equalsIgnoreCase("GREEDYBL-50"))
			placer = new GreedyFiller(GreedyFiller.MODE_BORDER_LENGTH, 50);		

		else if (placement.equalsIgnoreCase("GREEDY-CI-10K"))
			placer = new GreedyFiller(GreedyFiller.MODE_CONFLICT_INDEX, 10000);

		else if (placement.equalsIgnoreCase("GREEDY-CI-5K"))
			placer = new GreedyFiller(GreedyFiller.MODE_CONFLICT_INDEX, 5000);

		else if (placement.equalsIgnoreCase("GREEDY-CI-2K"))
			placer = new GreedyFiller(GreedyFiller.MODE_CONFLICT_INDEX, 2000);
		
		else if (placement.equalsIgnoreCase("GREEDY-CI-1K"))
			placer = new GreedyFiller(GreedyFiller.MODE_CONFLICT_INDEX, 1000);

		else if (placement.equalsIgnoreCase("GREEDY-CI-100"))
			placer = new GreedyFiller(GreedyFiller.MODE_CONFLICT_INDEX, 100);

		else if (placement.equalsIgnoreCase("GREEDY-CI-50"))
			placer = new GreedyFiller(GreedyFiller.MODE_CONFLICT_INDEX, 50);

		else if (placement.equalsIgnoreCase("REPTX20K-BL"))
			placer = new RowEpitaxial(20000, false);

		else if (placement.equalsIgnoreCase("2D1-SEQ"))
			placer = new TwoDimensionalPartitioning(new SequentialFiller(), 1);

		else if (placement.equalsIgnoreCase("2D4-SEQ"))
			placer = new TwoDimensionalPartitioning(new SequentialFiller(), 4);

		else if (placement.equalsIgnoreCase("2D1-GREEDY"))
			placer = new TwoDimensionalPartitioning(new GreedyFiller(), 1);

		else if (placement.equalsIgnoreCase("2D4-GREEDY"))
			placer = new TwoDimensionalPartitioning(new GreedyFiller(), 4);

		else if (placement.equalsIgnoreCase("2D8-GREEDY"))
			placer = new TwoDimensionalPartitioning(new GreedyFiller(), 8);

		else if (placement.equalsIgnoreCase("2D64-GREEDY"))
			placer = new TwoDimensionalPartitioning(new GreedyFiller(), 64);

		else if (placement.equalsIgnoreCase("2D8-GRASPD"))
			placer = new TwoDimensionalPartitioning(
						new QAPOptimization(
							new GraspDense(), QAPOptimization.MODE_CONFLICT_INDEX),
						8);

		else if (placement.equalsIgnoreCase("LEFT"))
			placer = new ProbeSetEmbeddingWrapper(new LeftMostEmbedding());
			
		else if (placement.equalsIgnoreCase("RIGHT"))
			placer = new ProbeSetEmbeddingWrapper(new RightMostEmbedding());
			
		else if (placement.equalsIgnoreCase("CENTER"))
			placer = new ProbeSetEmbeddingWrapper(new CenteredEmbedding());

		else if (placement.equalsIgnoreCase("EDGE"))
			placer = new ProbeSetEmbeddingWrapper(new EdgeEmbedding());

		else if (placement.equalsIgnoreCase("PIVOT"))
			placer = new ProbeSetEmbeddingWrapper(new PivotEmbedding());

		else if (placement.equalsIgnoreCase("PIVOTPART"))
			placer = new PivotPartitioning(new SequentialFiller());

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
			placer = new PivotPartitioning(new GreedyFiller(
							GreedyFiller.MODE_BORDER_LENGTH, 20000),
							PivotPartitioning.MODE_BORDER_LENGTH, 2);

		else if (placement.equalsIgnoreCase("PIVOTPART4+GREEDY20K-BL"))
			placer = new PivotPartitioning(new GreedyFiller(
							GreedyFiller.MODE_BORDER_LENGTH, 20000),
							PivotPartitioning.MODE_BORDER_LENGTH, 4);

		else if (placement.equalsIgnoreCase("PIVOTPART6+GREEDY20K-BL"))
			placer = new PivotPartitioning(new GreedyFiller(
							GreedyFiller.MODE_BORDER_LENGTH, 20000),
							PivotPartitioning.MODE_BORDER_LENGTH, 6);

		// **************************
		
		else if (placement.equalsIgnoreCase("PIVOTPART2+GREEDY2K-BL"))
			placer = new PivotPartitioning(new GreedyFiller(
							GreedyFiller.MODE_BORDER_LENGTH, 2000),
							PivotPartitioning.MODE_BORDER_LENGTH, 2);

		else if (placement.equalsIgnoreCase("PIVOTPART4+GREEDY2K-BL"))
			placer = new PivotPartitioning(new GreedyFiller(
							GreedyFiller.MODE_BORDER_LENGTH, 2000),
							PivotPartitioning.MODE_BORDER_LENGTH, 4);

		else if (placement.equalsIgnoreCase("PIVOTPART6+GREEDY2K-BL"))
			placer = new PivotPartitioning(new GreedyFiller(
							GreedyFiller.MODE_BORDER_LENGTH, 2000),
							PivotPartitioning.MODE_BORDER_LENGTH, 6);

		else if (placement.equalsIgnoreCase("PIVOTPART8+GREEDY2K-BL"))
			placer = new PivotPartitioning(new GreedyFiller(
							GreedyFiller.MODE_BORDER_LENGTH, 2000),
							PivotPartitioning.MODE_BORDER_LENGTH, 8);

		else if (placement.equalsIgnoreCase("PIVOTPART2+GREEDY2K-CI"))
			placer = new PivotPartitioning(new GreedyFiller(
							GreedyFiller.MODE_CONFLICT_INDEX, 2000),
							PivotPartitioning.MODE_CONFLICT_INDEX, 2);

		else if (placement.equalsIgnoreCase("PIVOTPART4+GREEDY2K-CI"))
			placer = new PivotPartitioning(new GreedyFiller(
							GreedyFiller.MODE_CONFLICT_INDEX, 2000),
							PivotPartitioning.MODE_CONFLICT_INDEX, 4);

		else if (placement.equalsIgnoreCase("PIVOTPART6+GREEDY2K-CI"))
			placer = new PivotPartitioning(new GreedyFiller(
							GreedyFiller.MODE_CONFLICT_INDEX, 2000),
							PivotPartitioning.MODE_CONFLICT_INDEX, 6);

		else if (placement.equalsIgnoreCase("PIVOTPART8+GREEDY2K-CI"))
			placer = new PivotPartitioning(new GreedyFiller(
							GreedyFiller.MODE_CONFLICT_INDEX, 2000),
							PivotPartitioning.MODE_CONFLICT_INDEX, 8);
		
		// **************************

		else if (placement.equalsIgnoreCase("PIVOTPART2+GREEDY1K-BL"))
			placer = new PivotPartitioning(new GreedyFiller(
							GreedyFiller.MODE_BORDER_LENGTH, 1000),
							PivotPartitioning.MODE_BORDER_LENGTH, 2);

		else if (placement.equalsIgnoreCase("PIVOTPART4+GREEDY1K-BL"))
			placer = new PivotPartitioning(new GreedyFiller(
							GreedyFiller.MODE_BORDER_LENGTH, 1000),
							PivotPartitioning.MODE_BORDER_LENGTH, 4);

		else if (placement.equalsIgnoreCase("PIVOTPART6+GREEDY1K-BL"))
			placer = new PivotPartitioning(new GreedyFiller(
							GreedyFiller.MODE_BORDER_LENGTH, 1000),
							PivotPartitioning.MODE_BORDER_LENGTH, 6);

		else if (placement.equalsIgnoreCase("PIVOTPART8+GREEDY1K-BL"))
			placer = new PivotPartitioning(new GreedyFiller(
							GreedyFiller.MODE_BORDER_LENGTH, 1000),
							PivotPartitioning.MODE_BORDER_LENGTH, 8);

		else if (placement.equalsIgnoreCase("PIVOTPART2+GREEDY1K-CI"))
			placer = new PivotPartitioning(new GreedyFiller(
							GreedyFiller.MODE_CONFLICT_INDEX, 1000),
							PivotPartitioning.MODE_CONFLICT_INDEX, 2);

		else if (placement.equalsIgnoreCase("PIVOTPART4+GREEDY1K-CI"))
			placer = new PivotPartitioning(new GreedyFiller(
							GreedyFiller.MODE_CONFLICT_INDEX, 1000),
							PivotPartitioning.MODE_CONFLICT_INDEX, 4);

		else if (placement.equalsIgnoreCase("PIVOTPART6+GREEDY1K-CI"))
			placer = new PivotPartitioning(new GreedyFiller(
							GreedyFiller.MODE_CONFLICT_INDEX, 1000),
							PivotPartitioning.MODE_CONFLICT_INDEX, 6);

		else if (placement.equalsIgnoreCase("PIVOTPART8+GREEDY1K-CI"))
			placer = new PivotPartitioning(new GreedyFiller(
							GreedyFiller.MODE_CONFLICT_INDEX, 1000),
							PivotPartitioning.MODE_CONFLICT_INDEX, 8);
		
		// **************************

		else if (placement.equalsIgnoreCase("PIVOTPART2+GREEDY100-BL"))
			placer = new PivotPartitioning(new GreedyFiller(
							GreedyFiller.MODE_BORDER_LENGTH, 100),
							PivotPartitioning.MODE_BORDER_LENGTH, 2);

		else if (placement.equalsIgnoreCase("PIVOTPART4+GREEDY100-BL"))
			placer = new PivotPartitioning(new GreedyFiller(
							GreedyFiller.MODE_BORDER_LENGTH, 100),
							PivotPartitioning.MODE_BORDER_LENGTH, 4);

		else if (placement.equalsIgnoreCase("PIVOTPART6+GREEDY100-BL"))
			placer = new PivotPartitioning(new GreedyFiller(
							GreedyFiller.MODE_BORDER_LENGTH, 100),
							PivotPartitioning.MODE_BORDER_LENGTH, 6);

		else if (placement.equalsIgnoreCase("PIVOTPART8+GREEDY100-BL"))
			placer = new PivotPartitioning(new GreedyFiller(
							GreedyFiller.MODE_BORDER_LENGTH, 100),
							PivotPartitioning.MODE_BORDER_LENGTH, 8);

		else if (placement.equalsIgnoreCase("PIVOTPART2+GREEDY100-CI"))
			placer = new PivotPartitioning(new GreedyFiller(
							GreedyFiller.MODE_CONFLICT_INDEX, 100),
							PivotPartitioning.MODE_CONFLICT_INDEX, 2);

		else if (placement.equalsIgnoreCase("PIVOTPART4+GREEDY100-CI"))
			placer = new PivotPartitioning(new GreedyFiller(
							GreedyFiller.MODE_CONFLICT_INDEX, 100),
							PivotPartitioning.MODE_CONFLICT_INDEX, 4);

		else if (placement.equalsIgnoreCase("PIVOTPART6+GREEDY100-CI"))
			placer = new PivotPartitioning(new GreedyFiller(
							GreedyFiller.MODE_CONFLICT_INDEX, 100),
							PivotPartitioning.MODE_CONFLICT_INDEX, 6);

		else if (placement.equalsIgnoreCase("PIVOTPART8+GREEDY100-CI"))
			placer = new PivotPartitioning(new GreedyFiller(
							GreedyFiller.MODE_CONFLICT_INDEX, 100),
							PivotPartitioning.MODE_CONFLICT_INDEX, 8);
		
		// **************************

		else if (placement.equalsIgnoreCase("CLUSTER-SEQ-BL"))
			placer = new ClusterPartitioning(new SequentialFiller(),
							PivotPartitioning.MODE_BORDER_LENGTH);

		else if (placement.equalsIgnoreCase("CLUSTER-SEQ-CI"))
			placer = new ClusterPartitioning(new SequentialFiller(),
							PivotPartitioning.MODE_CONFLICT_INDEX);

		// **************************

		else if (placement.equalsIgnoreCase("PIVOTPART8+GREEDY50-CI"))
			placer = new PivotPartitioning(new GreedyFiller(
							GreedyFiller.MODE_CONFLICT_INDEX, 50),
							PivotPartitioning.MODE_CONFLICT_INDEX, 8);

		
		// ****************** GRASP Sparse ******************
		
		else if (placement.equalsIgnoreCase("GRASPS-BL"))
			placer = new QAPOptimization (
							new GraspSparse(), QAPOptimization.MODE_BORDER_LENGTH
						);

		else if (placement.equalsIgnoreCase("GRASPS-INV-BL"))
		{
			QAPSolverAlgorithm qapsolver = new GraspSparse();
			qapsolver.setSwapMatrices(true);
			placer = new QAPOptimization (
							qapsolver, QAPOptimization.MODE_BORDER_LENGTH
						);
		}

		else if (placement.equalsIgnoreCase("GRASPS-CI"))
			placer = new QAPOptimization (
							new GraspSparse(), QAPOptimization.MODE_CONFLICT_INDEX
						);

		else if (placement.equalsIgnoreCase("GRASPS-INV-CI"))
		{
			QAPSolverAlgorithm qapsolver = new GraspSparse();
			qapsolver.setSwapMatrices(true);
			placer = new QAPOptimization (
							qapsolver, QAPOptimization.MODE_CONFLICT_INDEX
						);
		}

		// ****************** GRASP Dense ******************
		
		else if (placement.equalsIgnoreCase("GRASPD-BL"))
			placer = new QAPOptimization (
							new GraspDense(), QAPOptimization.MODE_BORDER_LENGTH
						);

		else if (placement.equalsIgnoreCase("GRASPD-INV-BL"))
		{
			QAPSolverAlgorithm qapsolver = new GraspDense();
			qapsolver.setSwapMatrices(true);
			placer = new QAPOptimization (
							qapsolver, QAPOptimization.MODE_BORDER_LENGTH
						);
		}

		else if (placement.equalsIgnoreCase("GRASPD-CI"))
			placer = new QAPOptimization (
							new GraspDense(), QAPOptimization.MODE_CONFLICT_INDEX
						);

		else if (placement.equalsIgnoreCase("GRASPD-INV-CI"))
		{
			QAPSolverAlgorithm qapsolver = new GraspDense();
			qapsolver.setSwapMatrices(true);
			placer = new QAPOptimization (
							qapsolver, QAPOptimization.MODE_CONFLICT_INDEX
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

		else if (optimization.equalsIgnoreCase("SEQREEMBED-BL-INC"))
			optimizer = new SequentialReembedding(
							OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN, true);
		
		else if (optimization.equalsIgnoreCase("SEQREEMBED-CI-INC"))
			optimizer = new SequentialReembedding(
							OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN, true);

		else if (optimization.equalsIgnoreCase("SEQREEMBED-BL-TOT"))
			optimizer = new SequentialReembedding(
							OptimumSingleProbeEmbedding.BORDER_LENGTH_MIN, false);
		
		else if (optimization.equalsIgnoreCase("SEQREEMBED-CI-TOT"))
			optimizer = new SequentialReembedding(
							OptimumSingleProbeEmbedding.CONFLICT_INDEX_MIN, false);

		// ****************** GRASP Path-relinking ******************

		else if (optimization.equalsIgnoreCase("GRASPPR32-BL"))
		{
			optimizer = new QAPOptimization (
							new GraspPathRelinking(1,32), QAPOptimization.MODE_BORDER_LENGTH
						);
		}

		else if (optimization.equalsIgnoreCase("GRASPPR32-INV-BL"))
		{
			QAPSolverAlgorithm qapsolver = new GraspPathRelinking(1,32);
			qapsolver.setSwapMatrices(true);
			optimizer = new QAPOptimization (
							qapsolver, QAPOptimization.MODE_BORDER_LENGTH
						);
		}

		else if (optimization.equalsIgnoreCase("GRASPPR32-CI"))
		{
			optimizer = new QAPOptimization (
							new GraspPathRelinking(1,32), QAPOptimization.MODE_CONFLICT_INDEX
						);
		}

		else if (optimization.equalsIgnoreCase("GRASPPR32-INV-CI"))
		{
			QAPSolverAlgorithm qapsolver = new GraspPathRelinking(1,32);
			qapsolver.setSwapMatrices(true);
			optimizer = new QAPOptimization (
							qapsolver, QAPOptimization.MODE_CONFLICT_INDEX
						);
		}

		// ****************** GRASP 128 Path-relinking ******************

		else if (optimization.equalsIgnoreCase("GRASPPR128-BL"))
		{
			optimizer = new QAPOptimization (
							new GraspPathRelinking(1, 128), QAPOptimization.MODE_BORDER_LENGTH
						);
		}

		else if (optimization.equalsIgnoreCase("GRASPPR128-INV-BL"))
		{
			QAPSolverAlgorithm qapsolver = new GraspPathRelinking(1, 128);
			qapsolver.setSwapMatrices(true);
			optimizer = new QAPOptimization (
							qapsolver, QAPOptimization.MODE_BORDER_LENGTH
						);
		}

		else if (optimization.equalsIgnoreCase("GRASPPR128-CI"))
		{
			optimizer = new QAPOptimization (
							new GraspPathRelinking(1, 128), QAPOptimization.MODE_CONFLICT_INDEX
						);
		}

		else if (optimization.equalsIgnoreCase("GRASPPR128-INV-CI"))
		{
			QAPSolverAlgorithm qapsolver = new GraspPathRelinking(1, 128);
			qapsolver.setSwapMatrices(true);
			optimizer = new QAPOptimization (
							qapsolver, QAPOptimization.MODE_CONFLICT_INDEX
						);
		}

		else
			throw new IllegalArgumentException
				("Unknown post-placement algorithm '" + optimization + "'");
		
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
		
		if (placer != null)
		{
			System.err.println("Running " + placement + " placement algorithm...");
			
			if (calc_bl)
			{
				System.err.println("Border length before placement: " +
						LayoutEvaluation.borderLength(chip));
			}
			else if (calc_ci)
			{
				System.err.println("Conflict index before placement: " +
						LayoutEvaluation.averageConflictIndex(chip));
			}
			
			start = System.nanoTime();

			// re-place probes on the chip
			unplaced = placer.makeLayout(chip);
			
			end = System.nanoTime();
			
			System.err.println("Elapsed time (placement): " +
					(end - start)/Math.pow(10,9) + " sec");
			
			if (unplaced > 0)
				System.err.println("WARNING: " + unplaced +
					" unplaced probe(s)");
		}

		if (optimizer != null)
		{
			if (calc_bl)
			{
				System.err.println("Border length before optimization: " +
						LayoutEvaluation.borderLength(chip));
			}
			else if (calc_ci)
			{
				System.err.println("Conflict index before optimization: " +
						LayoutEvaluation.averageConflictIndex(chip));
			}
			
			start = System.nanoTime();
			
			optimizer.optimizeLayout (chip);

			end = System.nanoTime();
			
			System.err.println("Elapsed time (optimization): " +
					(end - start)/Math.pow(10,9) + " sec");
		}

		if (calc_bl)
		{
			System.err.println("Final border length: " +
					LayoutEvaluation.borderLength(chip));
		}
		else if (calc_ci)
		{
			System.err.println("Final conflict index: " +
					LayoutEvaluation.averageConflictIndex(chip));
		}

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

	private static void usage ()
	{
		System.err.println (
	"--------------------------\n" +
	"ArrayOpt Microarray Design\n" +
	"--------------------------\n\n" +
	"Usage: ArrayOpt (affy | simple) (fix | nofix) <input> <rows> <columns>\n" +
	"          <probes> <length> <dep-seq> <placer> <optimizer>\n" +
	"          (check | no-check) (calc-bl | calc-ci | no-calc)\n" +
	"          (print | no-print)\n");
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
	"       <placer>    is a placement algorithm or NONE\n" +
	"       <optimizer> is a post-placement optimization or NONE\n" +
	"       'check'     performs validation checks\n" +
	"       'no-check'  does not perform validation checks\n" +
	"       'calc-bl'   computes and prints total border length\n" +
	"       'calc-ci'   computes and prints average conflict indices\n" +
	"       'no-calc'   does not compute/print any quality measure\n" +
	"       'print'     prints the resulting layout on standard output\n" +
	"       'no-print'  does not print the resulting layout\n");
	}
}
