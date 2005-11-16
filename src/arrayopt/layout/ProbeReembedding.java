/*
 * ProbeReembedding.java
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
 *
 */
public class ProbeReembedding
{
	public static void leftMostReembedding (Chip chip, int probe_id[], int first, int last)
	{
		for (int i = first; i <= last; i++)
			leftMostReembedding (chip, probe_id[i]);
	}

	public static void leftMostReembedding (Chip chip, int probe_id)
	{
		if (chip instanceof SimpleChip)
			leftMostReembedding ((SimpleChip) chip, probe_id, 0);

		else if (chip instanceof AffymetrixChip)
			leftMostReembedding ((AffymetrixChip) chip, probe_id, 0);

		else
			throw new IllegalArgumentException ("Unsupported chip type.");
	}

	protected static void leftMostReembedding (SimpleChip chip, int probe_id, int start_step)
	{
		// to do
	}

	protected static void leftMostReembedding (AffymetrixChip chip, int probe_id, int start_step)
	{
		// to do
	}

	public static void rightMostReembedding (Chip chip, int probe_id[], int first, int last)
	{
		for (int i = first; i <= last; i++)
			rightMostReembedding (chip, probe_id[i]);
	}

	public static void rightMostReembedding (Chip chip, int probe_id)
	{
		if (chip instanceof SimpleChip)
			rightMostReembedding ((SimpleChip) chip, probe_id, 0);

		else if (chip instanceof AffymetrixChip)
			rightMostReembedding ((AffymetrixChip) chip, probe_id, 0);

		else
			throw new IllegalArgumentException ("Unsupported chip type.");
	}

	protected static void rightMostReembedding (SimpleChip chip, int probe_id, int start_step)
	{
		// to do
	}

	protected static void rightMostReembedding (AffymetrixChip chip, int probe_id, int start_step)
	{
		// to do
	}

	public static void centeredReembedding (Chip chip, int probe_id[], int first, int last)
	{
		for (int i = first; i <= last; i++)
			centeredReembedding (chip, probe_id[i]);
	}

	public static void centeredReembedding (Chip chip, int probe_id)
	{
		if (chip instanceof SimpleChip)
			centeredReembedding ((SimpleChip) chip, probe_id);

		else if (chip instanceof AffymetrixChip)
			centeredReembedding ((AffymetrixChip) chip, probe_id);

		else
			throw new IllegalArgumentException ("Unsupported chip type.");
	}

	protected static void centeredReembedding (SimpleChip chip, int probe_id)
	{
		int shift = 0;

		leftMostReembedding (chip, probe_id, 0);

		// to do

		leftMostReembedding (chip, probe_id, shift);
	}

	protected static void centeredReembedding (AffymetrixChip chip, int probe_id)
	{
		int shift = 0;

		leftMostReembedding (chip, probe_id, 0);

		// to do

		leftMostReembedding (chip, probe_id, shift);
	}
}
