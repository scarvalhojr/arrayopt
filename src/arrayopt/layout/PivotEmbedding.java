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
 */
public class PivotEmbedding implements ProbeSetEmbeddingAlgorithm
{
	/**
	 * please document this
	 */
	public void reembedProbeSet (Chip chip, int probe_id[])
	{
		reembedProbeSet (chip, probe_id, 0, probe_id.length - 1);
	}

	/**
	 * please document this
	 */
	public void reembedProbeSet (Chip chip, int probe_id[], int first, int last)
	{
		// find pivots (probes that have only one possible embedding) and
		// move them to the beginning of the list

		// for each of the remaining probes p, find pivot q which minimizes
		// OptimumEmbedding.minHammingDistance(p, q)

		// reembed p optimally in regards to q:
		// OptimumEmbedding.reembedProbe(p, q)
	}

	/**
	 * please document this
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
	 * please document this
	 */
	public static int numberOfEmbeddings (SimpleChip chip, int probe_id)
	{
		// to do

		return 0;
	}

	/**
	 * please document this
	 */
	public static int numberOfEmbeddings (AffymetrixChip chip, int probe_id)
	{
		// to do

		return 0;
	}
}
