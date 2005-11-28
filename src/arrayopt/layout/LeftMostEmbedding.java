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

package arrayopt.layout;

/**
 * This class contains leftmost embedding algorithm for probes of several types of chips
 *
 * @author Anna Domanski & Ronny Gärtner
 */
public class LeftMostEmbedding extends ProbeEmbeddingAlgorithm
{
	private int oldmask, newmask, oldint, oldpos, newint, newpos;
	private int[] newembedding;
		
	/**
	 * reembeds a probe with leftmost embedding and no shift
	 */
	public void reembedProbe (Chip chip, int probe_id)
	{
		reembedProbe (chip, probe_id, 0);
	}

	/**
	 * reembeds a probe with leftmost embedding and a shift
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
	 * reembeds a probe of a simple chip with leftmost embedding and a shift which indicates at which
	 * step to start 
	 * method is not "pair-aware" means no synchronization with other probes
	 *
	 * IMPLEMENTATION NOTE: this method has package access to avoid abuse of
	 * the shift argument; it should only be used internally or by the
	 * {@linkplain CenteredEmbedding} class.
	 */
	void reembedProbe (SimpleChip chip, int probe_id, int shift)
	{
	oldmask = 0; 
	newmask=0;
	newpos = shift;
	newembedding = new int[chip.embed_len];
		

	for (newint = 0; newint < chip.embed[probe_id].length; newint++)
                       newembedding[newint] = 0;
	
	newint= (int) Math.floor((double) newpos/Integer.SIZE);
	
	for (oldint = -1, oldpos =0; oldpos < chip.embed_len; oldpos++)
	{	
		if (oldpos % Integer.SIZE == 0)
		{	
			oldint ++;
			oldmask=0x01 << (Integer.SIZE-1);
		}

		if (newpos % Integer.SIZE == 0)
		{
			newint++;
			newmask = 0x01 << (Integer.SIZE-1);
		}

		if (((chip.embed[probe_id][oldint]) &  oldmask) != 0)
		{
			
			while (chip.dep_seq[oldpos] != chip.dep_seq[newpos])	
			{
				newpos++;
				if (newpos % Integer.SIZE == 0)
				{
					newint++;
					newmask = 0x01 << (Integer.SIZE-1);
				}
				else
					newmask >>>=1;
			}
			

			if (newpos >= chip.dep_seq.length)
			{
				throw new IllegalArgumentException ("Unable to reembed probe.");
			}
			
			newembedding[newint] |= newmask;
			newpos++;
			if (newpos % Integer.SIZE == 0)
			{	
				newint++;
				newmask = 0x01 << (Integer.SIZE-1);
			}
			else
				newmask >>>=1;
		}

		oldmask >>>=1;
	}
	// overwrite old with new values
	for (newint = 0; newint < chip.embed[probe_id].length; newint++)
		chip.embed[probe_id][newint] = newembedding[newint];

	
		// to do:

		// reset bits up to position (shift - 1)

		// do a left-most embedding from position (shift)

		// reset remaining bits after finished

		// throw IllegalArgumentException ("Unable to reembed probe.")
		// if desired reembedding is not possible, e.g. because of a long shift
	}

	/**
	 * rembedds a probe of an Affymetrix Chip with leftmost embedding and pair-awareness
	 * @param shift indicates the start 
	 *
	 * IMPLEMENTATION NOTE: this method has package access to avoid abuse of
	 * the shift argument; it should only be used internally or by the
	 * {@linkplain CenteredEmbedding} class.
	 */
	void reembedProbe (AffymetrixChip chip, int probe_id, int shift)
	{
	
	
	oldmask = 0;
	newmask=0x01 << (Integer.SIZE-1-shift);
	newpos = shift;
	oldpos = 0;
	oldint = -1;
	newembedding = new int[chip.embed_len];
	int basenumber=0, leftborder=0, mpos = 0;
	boolean middle_synthesized = false;

	//turn all bits off
	for (newint = 0; newint < chip.embed[probe_id].length; newint++)
                       newembedding[newint] = 0;
	
	newint = (int) Math.floor((double) newpos/Integer.SIZE);
	
	
	// synthesize up to 12th base
	while (basenumber < chip.AFFY_MIDDLE_BASE)
	{
		if (oldpos % Integer.SIZE == 0)
		{	
			oldint ++;
			oldmask=0x01 << (Integer.SIZE-1);
		}
		
			if ((chip.embed[probe_id][oldint] &  oldmask) != 0 )
			{
				// set the right base character		
				while (chip.dep_seq[oldpos] != chip.dep_seq[newpos])	
				{
					newpos++;
					if (newpos >= chip.dep_seq.length)
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
				}
			//write base to new embedding
			newembedding[newint] |= newmask;
			newpos++;
			if (newpos % Integer.SIZE == 0)
			{
				newint++;
				newmask = 0x01 << (Integer.SIZE-1);
			}
			else
					newmask >>>=1;	
			basenumber++;	
			if (basenumber == chip.AFFY_MIDDLE_BASE)
				leftborder = (newpos-1);
			}
		oldpos++;
		oldmask >>>= 1;
	}
	if (basenumber == chip.AFFY_MIDDLE_BASE)
	{	
		while (!middle_synthesized)
		{
		if (oldpos % Integer.SIZE == 0)
		{	
			oldint ++;
			oldmask=0x01 << (Integer.SIZE-1);
		}
		
			if ((chip.embed[probe_id][oldint] &  oldmask) != 0 )
			{
			
				// set the right base character		
				while (chip.dep_seq[oldpos] != chip.dep_seq[newpos])	
				{
					newpos++;
					if (newpos >= chip.dep_seq.length)
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
				}
			
			
			newembedding[newint] |= newmask;
			middle_synthesized=true;
			mpos = newpos;
			basenumber++;
			}
		oldpos++;
		oldmask >>>= 1;
		}
	if (middle_synthesized)
		{
			if (chip.dep_seq[mpos] == chip.dep_seq[0])
			{
				if (chip.dep_seq[leftborder] == chip.dep_seq[0])
				{
					newpos++;
						if (newpos % Integer.SIZE == 0)
						{
							newint++;
							newmask = 0x01 << (Integer.SIZE-1);
						}
					
						else
							newmask >>>=1;
				}
				if (chip.dep_seq[leftborder] == chip.dep_seq[1])
				{
					newpos++;
						if (newpos % Integer.SIZE == 0)
						{
							newint++;
							newmask = 0x01 << (Integer.SIZE-1);
						}
					
						else
							newmask >>>=1;
				}
				if (chip.dep_seq[leftborder] == chip.dep_seq[2])
				{
					newpos++;
						if (newpos % Integer.SIZE == 0)
						{
							newint++;
							newmask = 0x01 << (Integer.SIZE-1);
						}
					
						else
							newmask >>>=1;
				}
				if (chip.dep_seq[leftborder] == chip.dep_seq[3])	
				{
					for (int i = 0; i < 4; i++)
					{
						newpos++;
						if (newpos % Integer.SIZE == 0)
						{
							newint++;
							newmask = 0x01 << (Integer.SIZE-1);
						}
						else
							newmask >>>=1;
					}
				}
			}
			if (chip.dep_seq[mpos] == chip.dep_seq[1])
			{
				if (chip.dep_seq[leftborder] == chip.dep_seq[0])
				{
						for (int i = 0; i < 2; i++)
						{
							newpos++;
							if (newpos % Integer.SIZE == 0)
							{
								newint++;
								newmask = 0x01 << (Integer.SIZE-1);
							}
							else
								newmask >>>=1;
				
						}
				}
				if (chip.dep_seq[leftborder] == chip.dep_seq[1])
				{
					newpos++;
						if (newpos % Integer.SIZE == 0)
						{
							newint++;
							newmask = 0x01 << (Integer.SIZE-1);
						}
					
						else
							newmask >>>=1;
				}
				if (chip.dep_seq[leftborder] == chip.dep_seq[2])
				{
					for (int i = 0; i < 2; i++)
					{
						newpos++;
						if (newpos % Integer.SIZE == 0)
						{
							newint++;
							newmask = 0x01 << (Integer.SIZE-1);
						}
						else
							newmask >>>=1;
					}
				}
				if (chip.dep_seq[leftborder] == chip.dep_seq[3])
				{
					for (int i = 0; i < 2; i++)
					{
						newpos++;
						if (newpos % Integer.SIZE == 0)
						{
							newint++;
							newmask = 0x01 << (Integer.SIZE-1);
						}
					
						else
							newmask >>>=1;
					}	
				}
			}
			if (chip.dep_seq[mpos] == chip.dep_seq[2])
			{	
				
				if (chip.dep_seq[leftborder] == chip.dep_seq[0])
				{
					newpos++;
						if (newpos % Integer.SIZE == 0)
						{
							newint++;
							newmask = 0x01 << (Integer.SIZE-1);
						}
						else
							newmask >>>=1;
				}
				if (chip.dep_seq[leftborder] == chip.dep_seq[1])
				{
					for (int i = 0; i < 4; i++)
					{
						newpos++;
						if (newpos % Integer.SIZE == 0)
						{
							newint++;
							newmask = 0x01 << (Integer.SIZE-1);
						}
					
						else
							newmask >>>=1;
					}
				}
				if (chip.dep_seq[leftborder] == chip.dep_seq[2])
				{
					newpos++;
						if (newpos % Integer.SIZE == 0)
						{
							newint++;
							newmask = 0x01 << (Integer.SIZE-1);
					
						}
						else
							newmask >>>=1;
				}
				if (chip.dep_seq[leftborder] == chip.dep_seq[3])
				{
					newpos++;
						if (newpos % Integer.SIZE == 0)
						{
							newint++;
							newmask = 0x01 << (Integer.SIZE-1);
						}
						else
							newmask >>>=1;
				}
			}
			if (chip.dep_seq[mpos] == chip.dep_seq[3])
			{
				if (chip.dep_seq[leftborder] == chip.dep_seq[0])
				{
					for (int i = 0; i < 2; i++)
					{
						newpos++;
						if (newpos % Integer.SIZE == 0)
						{
							newint++;
							newmask = 0x01 << (Integer.SIZE-1);
						}
					
						else
							newmask >>>=1;
					}
				}
				if (chip.dep_seq[leftborder] == chip.dep_seq[1])
				{
					for (int i = 0; i < 2; i++)
					{
						newpos++;
						if (newpos % Integer.SIZE == 0)
						{
							newint++;
							newmask = 0x01 << (Integer.SIZE-1);
						}
					
						else
							newmask >>>=1;
					}
				}
				if (chip.dep_seq[leftborder] == chip.dep_seq[2])
				{
					for (int i = 0; i < 2; i++)
					{
						newpos++;
						if (newpos % Integer.SIZE == 0)
						{
							newint++;
							newmask = 0x01 << (Integer.SIZE-1);
						}
					
						else
							newmask >>>=1;
					}
				}
				if (chip.dep_seq[leftborder] == chip.dep_seq[3])
				{
					newpos++;
						if (newpos % Integer.SIZE == 0)
						{
							newint++;
							newmask = 0x01 << (Integer.SIZE-1);
						}
				
						else
							newmask >>>=1;
				}
			}

		
	}
	// synthesize rest
	while (basenumber < chip.AFFY_PROBE_LENGTH)
	{
		if (oldpos % Integer.SIZE == 0)
		{	
			oldint ++;
			oldmask=0x01 << (Integer.SIZE-1);
		}
		
			if ((chip.embed[probe_id][oldint] &  oldmask) != 0 )
			{
				// set the right base character		
				while (chip.dep_seq[oldpos] != chip.dep_seq[newpos])	
				{
					newpos++;
					if (newpos >= chip.dep_seq.length)
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
				}
			
			
			newembedding[newint] |= newmask;
			if (basenumber != chip.AFFY_PROBE_LENGTH-1)
			{
			newpos++;
			if (newpos >= chip.dep_seq.length)
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
			}
			basenumber++;	
			}
			
		oldpos++;
		oldmask >>>= 1;
	}
	
	// write new embedding to probe
	for (newint = 0; newint < chip.embed[probe_id].length; newint++)
		chip.embed[probe_id][newint] = newembedding[newint];

	
	

		// to do:

		// same as with SimpleChip but must keep probe pairs synchronized
		// by dividing the task into three steps:
		// 1. do a left-most embedding up to base 12
		// 2. embed middle bases (pair)
		// 3. do a left-most embedding of remaining bases

		// note: from now on, the Affymetrix class can only have probes of length 25,
		// and you should use the constant called AFFY_MIDDLE_BASE to control the limits
		// of the three steps above
		}
	}
}