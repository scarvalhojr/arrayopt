  /*
 * LeftMostEmbedding.java
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
 * ?rrayOpt is distributed in the hope that it will be useful, but WITHOUT ANY
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
 * This class contains leftmost embedding algorithm for probes of several types of chips. 
 * Leftmost embedding means that each base of a probe is aligned at the leftmost possible spot
 * with the deposition sequence.
 * @author Anna Domanski & Ronny G?rtner
 */
public class LeftMostEmbedding implements SingleProbeEmbeddingAlgorithm,
	ProbeSetEmbeddingAlgorithm
{

	/**
	 * embedd a complete set of probes of a given chip
	 */
	public void reembedProbeSet (Chip chip, int probe_id[])
	{
		reembedProbeSet (chip, probe_id, 0, probe_id.length - 1);
	}

	/**
	 * embedd a set of probes of a chip within a range: from first to last probe
	 */
	public void reembedProbeSet (Chip chip, int probe_id[], int first, int last)
	{
		for (int i = first; i <= last; i++)
			reembedProbe (chip, probe_id[i]);
	}

	/**
	 * reembeds a probe with leftmost embedding and no shift
	 * @param chip either affymetrix or simple chip
	 * @param probe_id probe ID
	 */
	public void reembedProbe (Chip chip, int probe_id)
	{
		reembedProbe (chip, probe_id, 0);
	}

	/**
	 * reembeds a probe with leftmost embedding and a shift
	 * @param chip either affymetrix or simple chip
	 * @param probe_id probe ID
	 * @param shift sets starting position at which point to start with the embedding
	 */
	public void reembedProbe (Chip chip, int probe_id, int shift)
	{
		if (chip instanceof SimpleChip)
			reembedProbe ((SimpleChip) chip, probe_id, shift);

		else if (chip instanceof AffymetrixChip)
			reembedProbe ((AffymetrixChip) chip, probe_id, shift);

		else
			throw new IllegalArgumentException ("Unsupported chip type.");
	}

	/**
	 * reembeds a probe of a simple chip with leftmost embedding and a shift which indicates 
	 * at which step to start
	 *<B> IMPLEMENTATION NOTE:</B> this method has package access to avoid abuse of
	 * the shift argument; it should only be used internally or by the
	 * {@linkplain CenteredEmbedding} class.
	 *
	 * @param chip simple chip
	 * @param probe_id probe ID
	 * @param shift sets starting position at which point to start with the embedding
	 * @throws IllegalArgumentException if desired reembedding is not possible, e.g. because
	 * of a long shift
	 */
	void reembedProbe (SimpleChip chip, int probe_id, int shift)
	{
	int oldmask = 0, newmask, oldint, oldpos, newint, newpos = shift - 1;
	int[] newembedding = new int[chip.embed[probe_id].length];
		
    for (newint = 0; newint < chip.embed[probe_id].length; newint++)
                       newembedding[newint] = 0;

	newmask = 0x01 << ((Integer.SIZE)-(shift % Integer.SIZE));
	newint = (int) Math.floor((double) newpos/Integer.SIZE);
	
	for (oldint = -1, oldpos =0; oldpos < chip.embed_len; oldpos++)
	{	
		if (oldpos % Integer.SIZE == 0)
		{	
			oldint ++;
			oldmask=0x01 << (Integer.SIZE-1);
		}

		if (((chip.embed[probe_id][oldint]) &  oldmask) != 0)
		{
			
			while (++newpos < chip.dep_seq.length)
			{
				if (newpos % Integer.SIZE == 0)
				{
					newint++;
					newmask = 0x01 << (Integer.SIZE-1);
				}
				else
					newmask >>>=1;
			
				if (chip.dep_seq[oldpos] == chip.dep_seq[newpos])
					break;
			}
			
			if (newpos >= chip.dep_seq.length)
			{
				throw new IllegalArgumentException ("Unable to reembed probe.");
			}

			newembedding[newint] |= newmask;
		}

		oldmask >>>=1;
	}

	// overwrite old with new values - newembedding.length == embed.[probe_id].length
	System.arraycopy(newembedding,0,chip.embed[probe_id],0,chip.embed[probe_id].length);
	
		// to do:

		// reset bits up to position (shift - 1)

		// do a left-most embedding from position (shift)

		// reset remaining bits after finished

		// throw IllegalArgumentException ("Unable to reembed probe.")
		// if desired reembedding is not possible, e.g. because of a long shift
	}

	/**
	 * rembeds a probe of an affymetrix chip with leftmost embedding and pair-awareness
	 *
	 * <B>IMPLEMENTATION NOTE:</B> this method has package access to avoid abuse of
	 * the shift argument; it should only be used internally or by the
	 * {@linkplain CenteredEmbedding} class.
	 * @param chip affymetrix Chip
	 * @param probe_id probe ID
	 * @param shift sets starting position at which point to start with the embedding
	 * @throws IllegalArgumentException if desired reembedding is not possible, e.g. because 
	 * of a long shift
	 */
	void reembedProbe (AffymetrixChip chip, int probe_id, int shift)
	{
	int oldmask = 0, newmask = 0x01 << (Integer.SIZE-shift%Integer.SIZE), oldint = -1, oldpos = 0, newint, newpos = shift - 1;
	int basenumber = 0;
	int[] newembedding = new int[chip.embed[probe_id].length];
	boolean middle_synthesized = false, compl_synthesized = false;
	char compl = 'A'; //can be replaced if an appropriate method for returning the complement of a given base in class Chip exists

    //turn all bits off
	for (newint = 0; newint < chip.embed[probe_id].length; newint++)
                       newembedding[newint] = 0;
	
	newmask = 0x01 << ((Integer.SIZE-1)-(newpos % Integer.SIZE));
	newint = (int) Math.floor((double) newpos/Integer.SIZE);
		
	// synthesize up to 12th base
	while (basenumber < AffymetrixChip.AFFY_MIDDLE_BASE - 1)
	{

        if (oldpos % Integer.SIZE == 0)
		{	
			oldint++;
			oldmask=0x01 << (Integer.SIZE-1);
		}
		
		if ((chip.embed[probe_id][oldint] &  oldmask) != 0 )
		{
			// set the right base character		
			while (++newpos < chip.dep_seq.length)
			{
				if (newpos % Integer.SIZE == 0)
				{
					newint++;
					newmask = 0x01 << (Integer.SIZE-1);
				}
				else
					newmask >>>=1;
				if (chip.dep_seq[oldpos] == chip.dep_seq[newpos])
                {
				    break;
                }
			}
			if (newpos >= chip.dep_seq.length)
			{
				throw new IllegalArgumentException ("Unable to reembed probe.");
			}
			//write base to new embedding
			newembedding[newint] |= newmask;
			basenumber++;
		}
		
		oldpos++;
		oldmask >>>= 1;
		
	}

	while (!middle_synthesized)
	{
        if (oldpos % Integer.SIZE == 0)
		{	
			oldint++;
			oldmask=0x01 << (Integer.SIZE-1);
		}

		if ((chip.embed[probe_id][oldint] &  oldmask) != 0 )
		{
			// set the right base character		
			while (++newpos < chip.dep_seq.length)
			{
				if (newpos % Integer.SIZE == 0)
				{
					newint++;
					newmask = 0x01 << (Integer.SIZE-1);
				}
				else
					newmask >>>=1;
				// examine the complement of middle Base -- to do: can be replaced 
				// if an appropriate method for returning the complement of a 
				// given base in class Chip exists
				if (chip.dep_seq[oldpos] == 'A')
					compl = 'T';
				else if (chip.dep_seq[oldpos] == 'T')
					compl = 'A';
				else if	(chip.dep_seq[oldpos] == 'G')
					compl = 'C';
				else
					compl = 'G';
					
				if (compl == chip.dep_seq[newpos])
					compl_synthesized = true;
				
				if (chip.dep_seq[oldpos] == chip.dep_seq[newpos])
                {
                    break;
                } 
                   
			}
			
			if (newpos >= chip.dep_seq.length)
			{
				throw new IllegalArgumentException ("Unable to reembed probe.");
			}
			//write base to new embedding
			newembedding[newint] |= newmask;
			basenumber++;
			middle_synthesized = true;
			}
		oldpos++;
		oldmask >>>= 1;
	}
	int temp = oldpos;
	while (!compl_synthesized)
	{
        
		if (chip.dep_seq[temp] == compl)
			compl_synthesized = true;
		
		if (++newpos >= chip.dep_seq.length)
		{
           			throw new IllegalArgumentException ("Unable to reembed probe.");
		}
		if (newpos % Integer.SIZE == 0)
		{
			newint++;
			newmask = 0x01 << (Integer.SIZE-1);
		}
		else
			newmask >>>=1;
			
		temp++;
	}

	// synthesize rest
	
	while (basenumber < AffymetrixChip.AFFY_PROBE_LENGTH)
	{
	    if (oldpos % Integer.SIZE == 0)
		{	
			oldint++;
			oldmask=0x01 << (Integer.SIZE-1);
		}
		
		if ((chip.embed[probe_id][oldint] &  oldmask) != 0 )
		{
			// set the right base character		
			while (++newpos < chip.dep_seq.length)
			{
				if (newpos % Integer.SIZE == 0)
				{
					newint++;
					newmask = 0x01 << (Integer.SIZE-1);
				}
				else
					newmask >>>=1;
				if (chip.dep_seq[oldpos] == chip.dep_seq[newpos])
                {
                   	break;
                }
			}
			if (newpos >= chip.dep_seq.length)
			{
				throw new IllegalArgumentException ("Unable to reembed probe.");
			}
			//write base to new embedding
			newembedding[newint] |= newmask;
			basenumber++;
		}
		
		oldpos++;
		oldmask >>>= 1;
	}
	System.arraycopy(newembedding,0,chip.embed[probe_id],0,chip.embed[probe_id].length);
	
	

		// to do:

		// same as with SimpleChip but must keep probe pairs synchronized
		// by dividing the task into three steps:
		// 1. do a left-most embedding up to base 12
		// 2. embed middle bases (pair)
		// 3. do a left-most embedding of remaining bases

		// note: from now on, the Affymetrix class can only have probes of length 25,
		// and you should use the constant called AFFY_MIDDLE_BASE to control the limits
		// of the three steps above
		//}
	}
}