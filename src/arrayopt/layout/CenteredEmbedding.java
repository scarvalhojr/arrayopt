/*
 * CenteredEmbedding.java
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
public class CenteredEmbedding extends ProbeEmbeddingAlgorithm
{
	/**
	 * please document this
	 */
	public void reembedProbe (Chip chip, int probe_id)
	{
		if (chip instanceof SimpleChip)
			reembedProbe ((SimpleChip) chip, probe_id);

		else if (chip instanceof AffymetrixChip)
			reembedProbe ((AffymetrixChip) chip, probe_id);

		else
			throw new IllegalArgumentException ("Unsupported chip type.");
	}

	/**
	 * please document this
	 */
	public void reembedProbe (SimpleChip chip, int probe_id)
	{
		// to do:

		// do a left-most embedding and check how many masked steps
		// are left after the last synthesized probe, call it spaces

		// do a left-most embedding with a shift = spaces / 2
	}

	/**
	 * please document this
	 */
	public void reembedProbe (AffymetrixChip chip, int probe_id)
	{
		// to do:

		// do a left-most embedding and check how many masked steps
		// are left after the last synthesized probe, call it spaces

		// do a left-most embedding with a shift = spaces / 2
	}
}
