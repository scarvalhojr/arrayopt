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
	
	/**
	 * Returns the index of the median element. The three elements should be
	 * compared and the one with the median value should have its index
	 * returned. For instance, if <CODE>b < a < c</CODE>, then the method should
	 * return <CODE>a</CODE>.
	 * 
	 * @param a index of first element
	 * @param b index of second element
	 * @param c index of third element
	 * @return index of median element
	 */
	public abstract int medianOfThree (int a, int b, int c);
}
