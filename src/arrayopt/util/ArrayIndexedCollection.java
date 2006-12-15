/*
 * ArrayIndexedCollection.java
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

package arrayopt.util;

/**
 * This interface defines the operations needed by the {@link QuickSort} sorting
 * algorithm. These methods must be implemented to perform the necessary
 * operations on the data being sorted.
 * 
 * @author Sergio A. de Carvalho Jr.
 */
public interface ArrayIndexedCollection
{
	// TODO remove medianOfThree: QuickSort should implement it using compare
	
	/**
	 * Compares the elements at positions <CODE>a</CODE> and <CODE>b</CODE>. The
	 * result must be negative if <CODE>a</CODE> is less than <CODE>b</CODE>,
	 * zero if <CODE>a</CODE> equals <CODE>b</CODE> or positive if
	 * <CODE>a</CODE> is greater than <CODE>b</CODE>.
	 *  
	 * @param a index of first element
	 * @param b index of second element
	 * @return -1, 0 or +1
	 */
	public abstract int compare (int a, int b);
	
	/**
	 * Swaps the elements at positions <CODE>a</CODE> and <CODE>b</CODE>.
	 * 
	 * @param a index of first element
	 * @param b index of second element
	 */
	public abstract void swap (int a, int b);
	
	/**
	 * Sets one element as the pivot that will be used in later comparisons.
	 * This element will be used in all further comparisons performed by the
	 * {@link #compareToPivot(int)} method.  
	 * 
	 * @param p index of element to be set as pivot
	 */
	public abstract void setPivot (int p);
	
	/**
	 * Compares the element indexed at position <CODE>a</CODE> with the pivot
	 * previously set by the {@link #setPivot(int)} method. As with the
	 * {@link #compare(int, int)} method, the result must be negative if
	 * <CODE>a</CODE> is less than <CODE>b</CODE>, zero if <CODE>a</CODE> equals
	 * <CODE>b</CODE> or positive if <CODE>a</CODE> is greater than
	 * <CODE>b</CODE>.
	 * 
	 * @param a index of element to be compared with the pivot
	 * @return -1, 0 or +1
	 */
	public abstract int compareToPivot (int a);
}
