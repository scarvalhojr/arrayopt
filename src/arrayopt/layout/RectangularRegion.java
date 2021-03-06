/*
 * RectangularRegion.java
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
public class RectangularRegion implements Region
{
	public int first_row;

	public int last_row;

	public int first_col;

	public int last_col;

	public RectangularRegion (int first_row, int last_row, int first_col, int last_col)
	{
		this.first_row = first_row;
		this.last_row  = last_row;
		this.first_col = first_col;
		this.last_col  = last_col;
	}
	
	public int getNumberOfRows ()
	{
		return last_row - first_row + 1;
	}

	public int getNumberOfColumns ()
	{
		return last_col - first_col + 1;
	}

	@Override
	public String toString ()
	{
		return first_row + "-" + last_row + ";" + first_col + "-" + last_col;
	}
}
