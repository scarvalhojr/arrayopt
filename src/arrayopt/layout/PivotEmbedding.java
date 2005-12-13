/*
 * PivotEmbedding.java
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
 * ???rrayOpt is distributed in the hope that it will be useful, but WITHOUT ANY
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

/**
 * This class represents the pivot embedding method for affymetrix
 * and simple chips. Depending on the number of possible embeddings 
 * of a probe it will arrange probes close to probes which are similar
 * embedded so that the number of conflicts is reduced. Pivot embedding 
 * gives a higher priority to probes with a low number of possible embeddings
 * because probes with a higher number can fit in more places. So that 
 * probes with only one possible embedding, called pivots, are placed first.
 * @author Anna Domanski & Ronny Gaertner 
 */
public class PivotEmbedding implements ProbeSetEmbeddingAlgorithm
{
	/**
	 * Reembedds a whole set of probes in such a way that probes get their optimal 
	 * embedding by ordering probes in regards to their possible embeddings. Probes 
	 * with only one possible embedding (pivots) are the probes on which probes 
	 * with a higher number of possible embeddings are assigned to the pivot that 
	 * minimizes the hamming distance to the probe. After that the probe will be 
	 * embedded optimal in regards to its pivot.
	 */
	public void reembedProbeSet (Chip chip, int probe_id[])
	{
		reembedProbeSet (chip, probe_id, 0, probe_id.length - 1);
	}

	/**
	 * Reembedds a part of a set of probes of a chip starting with first and ending
	 * with last. Reembedds in such a way that probes get their optimal embedding
	 * by ordering probes in regards to their possible embeddings. Probes with only
	 * one possible embedding (pivots) are the probes on which probes with a higher
	 * number of possible embeddings are assigned to the pivot that minimizes the 
	 * hamming distance to the probe. After that the probe will be embedded optimal 
	 * in regards to its pivot.
	 * @param chip chip on which probes reside
	 * @param probe_id[] probe IDs
	 * @param first probe ID of probe to be considered first
	 * @param lasst probe ID of porbe to be considered last
	 */
	public void reembedProbeSet (Chip chip, int probe_id[], int first, int last)
	{
		int temp;
		int border = 0;
		int minhamdistance = Integer.MAX_VALUE;
		int pivot = 0;
		
		// find pivots (probes that have only one possible embedding) and
		// move them to the beginning of the list
		for (int i = first; i <= last;i++)
		{
			if (numberOfEmbeddings(chip, probe_id[i]) == 1)
			{
				temp = probe_id[i];
				probe_id[i] = probe_id[border];
				probe_id[border] = temp;
				border++;
			}
		}
		
		// for each of the remaining probes p, find pivot q which minimizes
		for (int p = border; p <= last;p++)
		{
				for (int q = 0; q < border;q++)
				{
					if(minhamdistance > OptimumEmbedding.minHammingDistance(chip,probe_id[p],probe_id[q]))
							minhamdistance=OptimumEmbedding.minHammingDistance(chip,probe_id[p],probe_id[q]);
							pivot = q;
				}
			
			// reembed p optimally in regards to pivot
			OptimumEmbedding.reembedProbe(chip, p, pivot);
		}
	}

	/**
	 * Returns the number of possible embeddings of a probe in regards to the 
	 * deposition sequence of the chip.
	 * @param chip chip on which the probe resides
	 * @param probe_id probe ID
	 * @return number of possible embeddings 
	 */
	public static int numberOfEmbeddings (Chip chip, int probe_id)
	{
		if (chip instanceof SimpleChip)
			return numberOfEmbeddings ((SimpleChip) chip, probe_id);

		else if (chip instanceof AffymetrixChip)
			return numberOfEmbeddings ((AffymetrixChip) chip, probe_id);

		else
			throw new IllegalArgumentException
				("Unsupported chip type.");
	}

	/**
	 * Returns the number of possible embeddings of a probe residing 
	 * on a Simple Chip.
	 * @param chip simple chip
	 * @param probe_id probe ID
	 * @return number of possible embeddings
	 */
	public static int numberOfEmbeddings (SimpleChip chip, int probe_id)
	{
		
		// try to access "BitSet" directly
		int[] current = new int[(chip.embed_len + 1)];
		int [] previous = new int[(chip.embed_len + 1)];
		
		
		previous[0] = 1;
		current[0] = 1;
		
		int mask = 0x01 << (Integer.SIZE-1);
		
		for (int  j = 0; j < chip.dep_seq.length; j++)
		{
			int pos = 0;
			for(int i = 1; i < chip.embed_len; i++, mask >>= 1)
			{
				if(i % Integer.SIZE == 0)
				{
					mask = 0x01 << (Integer.SIZE-1);
					pos++;
				}
				
				if((chip.embed[probe_id][pos] & mask) != 0)
				{
					//may not be necessary, because previous == current
					current[i] = previous[i];   
					
					if(chip.dep_seq[i-1] == chip.dep_seq[j])
					{
						current[i] += previous[i-1];
					}
					
					if (current[i] == 0)
						break;
				}
			}
			previous = current;
		}

		return current[chip.embed_len];  
	}

	/**
	 * Returns the number of possible embeddings of the Bases of a probe aligned 
	 * the deposition sequence that was used to build the affymetrix chip.
	 * @param chip affymetrix chip
	 * @param probe_id probe ID
	 * @return number of possible embeddings of the probe  
	 */
	public static int numberOfEmbeddings (AffymetrixChip chip, int probe_id)
	{
		//variable for the PM & MM merged embedding
		int[] mergedembedding = new int[chip.embed[probe_id].length];
		//probe_id of the partner of given probe (PM <-> MM) 
		int pair_probe_id;
		
		int[] current = new int[(chip.probe_len + 1)];
		int [] previous = new int[(chip.probe_len + 1)];
		boolean debug = true;
		// int counter = 0;
		
		//setting the probe_id of the partner in regards to the given probe is PM or MM
		if (chip.isPMProbe(probe_id))
				pair_probe_id = probe_id + 1;
		else
				pair_probe_id = probe_id - 1;
		
		if (debug)
			System.out.println("mergedembedding["+2+"] = "+mergedembedding[2] + " PM: " + chip.embed[probe_id][2] + " MM: "+ chip.embed[pair_probe_id][2]);

		// creating embedding of PM and MM merged
		for(int i = 0; i < chip.embed[probe_id].length; i++)
		{
			mergedembedding[i] = chip.embed[probe_id][i] | chip.embed[pair_probe_id][i];
			if (debug)
				System.out.println("mergedembedding["+i+"] = "+mergedembedding[i] + " PM: " + chip.embed[probe_id][i] + " MM: "+ chip.embed[pair_probe_id][i]);
		}
		int basecount = Integer.bitCount(mergedembedding[0]) + Integer.bitCount(mergedembedding[1]) + Integer.bitCount(mergedembedding[2]);
		if (debug)
			System.out.println("Anzahl Basen in merged: "+ basecount);
		
		
		char[] probe = new char[chip.probe_len+1];
		int pos = -1;
		int mask = 0x01 << (Integer.SIZE-1);
		int charpos = 0;
		int counter =0;
		
		//convert BitSet to String  <-- have to fix it, but couldn't figure out failure
		for (int i = 0; i < chip.dep_seq.length; i++)
		{
			if(i % Integer.SIZE == 0)
			{
				mask = 0x01 << (Integer.SIZE-1);
				pos++;
			}
		
			if((mergedembedding[pos] & mask) != 0)
			{
				probe[charpos]=chip.dep_seq[i];
				System.out.println(counter++ +"te mal match!"+" total: "+ i+ " "+chip.dep_seq[i]+ " "+(mergedembedding[pos] & mask) );
				//charpos++;
			}
			mask >>= 1;
		}
		
		//initialize 2 rows of Embeddingsmatrix
		previous[0] = 1;
		current[0] = 1;	
		for(int i =1; i <= chip.probe_len;i++)
		{	previous[i] = 0;
			current[i] = 0;
		}	
		for (int  j = 0; j < chip.dep_seq.length; j++)
		{
			for(int i = 1; i <= chip.probe_len; i++)
			{ 	
				if(probe[i-1] == chip.dep_seq[j])
				{
					// may not be necessary, because previous == current
					// current[i] = previous[i];   
					current[i] += previous[i-1];
				}
					
				if (current[i] == 0)
						break;
				
			}
			previous = current;
		}

		return current[chip.probe_len];	
		
			
			
			
			
			
			
		/*	
			
			// try to access "BitSet" directly
		for(int i =1; i<= chip.embed_len;i++)
		{	previous[i] = 0;
			current[i] = 0;
		}	
		previous[0] = 1;
		current[0] = 1;
		
		int mask = 0x01 << (Integer.SIZE-1);
		int counter2 = 0;
		for (int  j = 0; j < chip.dep_seq.length; j++)
		{
			if (debug)
				System.out.println(counter2 +". Base in DepSeq");
			counter = 0;
			int pos = -1;
			for(int i = 1; i <= chip.embed_len; i++, mask >>= 1)
			{ 	
				if (debug)
					System.out.println(i +". Base in Probe");
				if((i-1) % Integer.SIZE == 0)
				{
					mask = 0x01 << (Integer.SIZE-1);
					pos++;
				}
				if (debug)
					System.out.println("maske: " + (mergedembedding[pos] & mask));
				if((mergedembedding[pos] & mask) != 0)
				{
					if (debug)
					{
						System.out.println(counter + ".mal match "+i+".th Wert: "+ current[i]);
						counter ++;
					}
					// may not be necessary, because previous == current
					// current[i] = previous[i];   
					last = i;
					if(chip.dep_seq[i-1] == chip.dep_seq[j])
					{
						current[i] += previous[i-1];
						if (debug)
						{
							System.out.println("extra");
						}
					}
					
					if (current[i] == 0)
						break;
				}
				else
					current[i]=current[i-1];
			}
			counter2++;
			previous = current;
		}

		return current[last]; */ 
	}	
}
