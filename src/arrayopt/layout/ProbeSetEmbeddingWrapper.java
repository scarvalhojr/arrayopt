/*
 * ProbeSetEmbeddingWrapper.java
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
 * �rrayOpt is distributed in the hope that it will be useful, but WITHOUT ANY
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
 *
 */
public class ProbeSetEmbeddingWrapper implements PlacementAlgorithm
{
	protected ProbeSetEmbeddingAlgorithm embedder;

	/**
	 *
	 */
	public ProbeSetEmbeddingWrapper (ProbeSetEmbeddingAlgorithm embedder)
	{
		this.embedder = embedder;
	}
	
	/**
	 *
	 */
	public int makeLayout (Chip chip)
	{
		// DO NOT reset current layout
		
		// reembed all probes
		
		/*
		if (embedder instanceof LeftMostEmbedding)
		{
			LeftMostEmbedding lm = (LeftMostEmbedding) embedder;
			int	id[] = chip.getAllProbes();
			
			for (int i = 0; i < id.length; i++)
				lm.reembedProbe (chip, id[i], 4);
		}
		else
		*/
			embedder.reembedProbeSet (chip, chip.getAllProbes());
		
		// since no probe has been moved,
		// no probe may be left out
		return 0;
	}
}