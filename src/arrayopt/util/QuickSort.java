package arrayopt.util;

/**
 * This class implements the QuickSort sorting algorithm. This implementation
 * is intended to sort collections of objects or any kind of data that is
 * indexed by an array. There is no dependence on the type of objects being
 * sorted nor on the index structure. The sorting is achieved by calling
 * methods of the {@link ArrayIndexedCollection} interface, which defines
 * methods for querying the ordering between the elements and performing swaps.
 * 
 * <P>This implementation is designed to serve as an alternative to the sorting
 * mechanisms provided by the {@link java.util.Arrays} class and Java's
 * Collections framework. The first can only be used to sort simple arrays of
 * numbers and its usefulness is therefore limited. (For instance, it cannot be
 * used to sort one array index by another.) The latter is based on interfaces
 * and generics, and can thus be used to sort anything. However, the Collections
 * framework forces the data to be wrapped as objects, which can be prohibitive
 * in terms of memory if the number of elements is too large.</P>
 * 
 * <P>The solution is a compromise achieved by an interface that defines the
 * necessary sorting operations (comparison, swap, etc.). The result is an
 * implementation that is as general as the Collections framework, but with the
 * drawback of having to provide code that performs the operations on the
 * data.</P> 
 * 
 * <P>In order to use this class to sort a collection of objects or any data
 * indexed by an array, an implementation of the
 * {@link arrayopt.util.ArrayIndexedCollection} interface is required. This
 * implementation is resposible for the comparisons between the elements and for
 * performing the necessary swaps in order to achive the desired ordering.</P>
 * 
 * <P>This implementation is based on the QuickSort code of the
 * {@link java.util.Arrays} class.</P>
 * 
 * @author Sergio A. de Carvalho Jr.
 */
public class QuickSort
{
	private static final long SMALL = 7;

	private static final long BIG = 41;

	/**
	 * Sorts data indexed by an array with the QuickSort algorithm. The data
	 * must be encapsulated by a class that performs the necessary operations
	 * defined by the {@link ArrayIndexedCollection} interface
	 *  
	 * @param col collection to be sorted 
	 * @param off index of first element to be sorted 
	 * @param len number of elements to be sorted
	 */
	public static void sort (ArrayIndexedCollection col, int off, int len)
	{
		int i, j, k, l, m, cmp;
		
		// insertion sort on smallest arrays
		if (len < SMALL)
		{
		    for (i = off + 1; i < off + len; i++)
		    	for (j = i; j > off && col.compare(j - 1, j) > 0; j--)
		    		col.swap(j - 1, j);
		    
		    return;
		}
		
		// choose a partition element, v
		
		// small arrays: middle element
		m = off + (len >> 1);
		
		if (len > SMALL)
		{
		    i = off;
		    j = off + len - 1;
		    
		    // big arrays: pseudomedian of 9
		    if (len >= BIG)
		    {
		    	k = len / 8;
		    	i = col.medianOfThree(i, i + k, i + 2 * k);
				m = col.medianOfThree(m - k, m, m + k);
				j = col.medianOfThree(j - 2 * k, j - k, j);
			}
		    
		    // mid-size: med of 3
		    m = col.medianOfThree(i, m, j);
		}
		
		col.setPivot(m);
		
		// establish invariant: (x = m), (x < m), (x > m), (x = m)
		i = j = off;
		k = l = off + len - 1;
		while(true)
		{
		    while (j <= k && (cmp = col.compareToPivot(j)) <= 0)
		    {
		    	if (cmp == 0)
		    		col.swap(i++, j);
		    	j++;
		    }
		    
		    while (k >= j && (cmp = col.compareToPivot(k)) >= 0)
		    {
		    	if (cmp == 0)
		    		col.swap(k, l--);
		    	k--;
		    }
		    
		    if (j > k) break;
		    
		    col.swap(j++, k--);
		}
		
		// swap partition elements back to middle
		m = Math.min(i - off, j - i);
		vectorSwap (col, off, j - m, m);
		
		m = Math.min(l - k, off + len - l - 1);
		vectorSwap (col, j, off + len - m, m);
		
		// recursively sort non-partition-elements
		if ((m = j - i) > 1)
		    sort(col, off, m);
		if ((m = l - k) > 1)
		    sort(col, off + len - m, m);
	}
	
	/**
	 * Internal method to swap parts of an array. It is based on swaps of
	 * individual elements.
	 * 
	 * @param col collection containing the array
	 * @param a starting position of the first part
	 * @param b starting position of the second part
	 * @param n number of elements
	 */
	private static void vectorSwap (ArrayIndexedCollection col, int a, int b,
			int n)
	{
		for (int i = 0; i < n; i++)
		    col.swap(a++, b++);
	}
}
