/*
 * ProbeEmbeddingAlgorithm.java
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
 * This class contains different embedding algorithms for probes. 
 * Each probe is a subset of the deposition sequence, i.e. this deposition sequence is a superstring of all
 * probes. They can be embedded in several ways. Each embedding is an alignment with the deposition sequence 
 * so that matches and gaps appear. Embedding algorithms generate such an alignment. 
 * Probes normally appear in tuples but they may also appear in single instances or pairs. So every algorithm 
 * is responsible for keeping pairs (tuples) synchronized if needed.
 *
 * @author Anna Domanski & Ronny Gärtner
 */
public abstract class ProbeEmbeddingAlgorithm
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
	 * embedd a single probe of a chip
	 */
	public abstract void reembedProbe (Chip chip, int probe_id);
}
