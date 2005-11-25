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
 * please document this
 *
 * @author WRITE YOUR NAME HERE
 */
public class LeftMostEmbedding extends ProbeEmbeddingAlgorithm
{
	/**
	 * please document this
	 */
	public void reembedProbe (Chip chip, int probe_id)
	{
		reembedProbe (chip, probe_id, 0);
	}

	/**
	 * please document this
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
	 * please document this
	 *
	 * IMPLEMENTATION NOTE: this method has package access to avoid abuse of
	 * the shift argument; it should only be used internally or by the
	 * {@linkplain CenteredEmbedding} class.
	 */
	void reembedProbe (SimpleChip chip, int probe_id, int shift)
	{
	/*int mask = 0,intpos, oldembed, newembed, deppos = shift,
	
	
	for (intpos = -1, oldembed =0; intpos < chip.embed[prob_id].length; oldembed++, intpos++)
	{	
		if (oldembed % Integer.Size == 0)
		{	
			intpos ++;
			mask=0x01 << (Integer.Size-1);
		}

		if (chip.embed[probe_id][intpos] & mask)
		{
			while (chip.dep_seq[oldembed] != chip.dep_seq[deppos])	
	}*/
		// to do:

		// reset bits up to position (shift - 1)

		// do a left-most embedding from position (shift)

		// reset remaining bits after finished

		// throw IllegalArgumentException ("Unable to reembed probe.")
		// if desired reembedding is not possible, e.g. because of a long shift
	}

	/**
	 * please document this
	 *
	 * IMPLEMENTATION NOTE: this method has package access to avoid abuse of
	 * the shift argument; it should only be used internally or by the
	 * {@linkplain CenteredEmbedding} class.
	 */
	void reembedProbe (AffymetrixChip chip, int probe_id, int shift)
	{
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
