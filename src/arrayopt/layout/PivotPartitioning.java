/*
 * PivotPartitioning.java
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
public class PivotPartitioning implements PlacementAlgorithm
{
	/**
	 * please document this
	 */
	private FillingAlgorithm filler;

	/**
	 * please document this
	 */
	public static final int DEFAULT_STOP_DIMENSION = 4;

	/**
	 * please document this
	 */
	private int stop_dimension;

	/**
	 * please document this
	 */
	public PivotPartitioning (FillingAlgorithm filler)
	{
		this(filler, DEFAULT_STOP_DIMENSION);
	}

	/**
	 * please document this
	 */
	public PivotPartitioning (FillingAlgorithm filler, int stop_dimension)
	{
		this.filler = filler;
		this.stop_dimension = (stop_dimension < 1) ? 1 : stop_dimension;
	}

	/**
	 * please document this
	 */
	public int makeLayout (Chip chip)
	{
		// reset current layout (if any)
		chip.resetLayout();

		// to do: proceed recursively

		// should return the number of unplaced probes (if any)
		return 0;
	}
}
