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
	 * @param last probe ID of probe to be considered last
	 */
	public void reembedProbeSet (Chip chip, int probe_id[], int first, int last)
	{
        OptimumEmbedding embedder = OptimumEmbedding.createEmbedder(chip, OptimumEmbedding.MODE_CONFLICT_INDEX);
		int temp;
		int border = 0;
		double minhamdistance = Double.MAX_VALUE;
        int pivot = 0;
		
		// find pivots (probes that have only one possible embedding) and
		// move them to the beginning of the list --> what if there are no pivots?
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
                    double compare = embedder.minDistanceProbe(probe_id[p],probe_id[q]);
                    if(minhamdistance > compare)
                    {
							minhamdistance = compare;
							pivot = q;
                    }
				}
			
			// reembed p optimally in regards to pivot
			embedder.reembedProbe(p, pivot);
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
		
	    // columns of the embeddingsmatrix - only 2 needed since we are 
        // only interested of the value in the very last cell of this matrix 
        // and one column is computed by the values of the previous column.
        // each column has one additional value -> [0] represents the blank character
		int[] current = new int[(chip.embed_len + 1)];
		int[] previous = new int[(chip.embed_len + 1)];
		
		// at the end of computing the matrix last should contain the number of possible embeddings
        int last = 0;
        
        int mask = 0x01 << (Integer.SIZE-1);
        
        // make sure arrays are initialized with zero!
        for(int i =1; i < chip.embed_len+1 ;i++)
        {   
            previous[i] = 0;
            current[i] = 0;
        }   
        
        // value of first row (blank character) of matrix that will never be changed!
		previous[0] = 1;
		current[0] = 1;
        
        // j indicates the current column of our embeddings matrix where every column is labeled 
        // with the corresponding character of the deposition sequence
		for (int  j = 0; j < chip.dep_seq.length; j++)
		{
		    // counter is used as a pointer that goes over the bitset of our probe to see which bits are set, therefore normally: counter < dep_seq.length
            int counter = 0;
            
            // refers to that integer whose bitset will be accessed
            int pos = -1;
            
            
            // indicates the current row of our embeddings matrix where every row is labeled with the corresponding  
            // character of the probe (i=0 -> blank character) --> j and i point to a cell in the matrix
			for(int i = 0; i < chip.probe_len; counter++,mask >>= 1)
			{
				if(counter % Integer.SIZE == 0)
				{
					mask = 0x01 << (Integer.SIZE-1);
					pos++;
				}
				
                // see if bit is set in the probe --> position of counter indicates which character it is
				if((chip.embed[probe_id][pos] & mask) != 0)
				{
                    i++;
                    last = i;			
					if(chip.dep_seq[counter] == chip.dep_seq[j])
					{
						current[i] += previous[i-1];
					}
					
					if (current[i] == 0)
                    {
                        break;
                    }
				}
			}
            
			// previous = current
            System.arraycopy(current,1,previous,1,last);
		}

		return current[chip.probe_len];  
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
	    int[] merged_embedding = new int[chip.embed[probe_id].length];
	    
	    //probe_id of the partner of given probe (PM <-> MM) 
	    int pair_probe_id;
        
        // merged_embedding has 1 base more than normal probe and one additional value for the first row of the embeddingsmatrix 
        int merged_embedding_length = chip.probe_len + 2;
        
        // columns of the embeddingsmatrix - only 2 needed since we are 
        // only interested of the value in the very last cell of this matrix 
        // and one row is computed by the values of the previous row
	    int[] current = new int[merged_embedding_length];
	    int [] previous = new int[merged_embedding_length];
	    
        // mask for accessing the bitset of the embedding
        int mask = 0x01 << (Integer.SIZE-1);
        
	    
	    
	    //setting the probe_id of the partner in regards to the given probe is PM or MM
	    if (chip.isPMProbe(probe_id))
	    {
            pair_probe_id = probe_id + 1;
        }
	    else
        {
            pair_probe_id = probe_id - 1;
        }
	    
	    // creating embedding of PM and MM merged <-- there should be a better way of creating a correct merdedembedding
	    for(int i = 0; i < chip.embed[probe_id].length; i++)
	    {
	        merged_embedding[i] = chip.embed[probe_id][i] | chip.embed[pair_probe_id][i];
	    }
	      
	    // make sure arrays are initialized with zero!
	    for(int i =1; i< merged_embedding_length ;i++)
	    {	
	        previous[i] = 0;
	        current[i] = 0;
	    }	
	    
        // value of first row of matrix that will never be changed!
	    previous[0] = 1;
	    current[0] = 1;
	    
        // at the end of computing the matrix last should contain the number of possible embeddings
	    int last = 0;
	    
        // j indicates the current column of our embeddings matrix where every column is labeled with the corresponding 
        // character of the deposition sequence
	    for (int  j = 0; j < chip.dep_seq.length; j++)
	    {
            // counter is used as a pointer that goes over the bitset of our merged probe to see which bits are set, therefore normally: counter < dep_seq.length
            int counter = 0;
            
            // refers to that integer whose bitset will be accessed
            int pos = -1;
            
            // indicates the current row of our embeddings matrix where every row is labeled with the corresponding 
            // character of the merged probe
            // since it is a merged probe, merged_probe_len = 26
            for(int i = 0; i <= chip.probe_len; counter++, mask >>>= 1)
	        { 	
	           	            
	            if(((counter) % Integer.SIZE) == 0)
	            {
	                mask = 0x01 << (Integer.SIZE-1);
	                pos++;
	            }
	            
                // see if bit is set for the merged_probe --> position of counter indicates which character it is
	            if(((merged_embedding[pos]) & mask) != 0)
	            {
                    i++;
	          	    last = i;
	                
                    if(chip.dep_seq[counter] == chip.dep_seq[j])
	                {
	                    current[i] += previous[i-1];
	                }
	                
	                if (current[i] == 0 || i == chip.probe_len+2)
	                {
	                    break;
	                }
	                
	            }
             
	            
	        }
	        // previous = current;
            System.arraycopy(current,1,previous,1,last); 
	    }
	    return current[merged_embedding_length-1];  
	}	
}
