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
import java.util.Vector;


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

    private static final double PIVOT_THRESHOLD = 0.2;
    
    private static final double DIV_RATE_PIVOT = .2;
    private static final double DIV_RATE_NON_PIVOT = .2;
    /**
     * Value indicates at which percentage of the total size of the to be sorted array merge sort will stop 
     * the division and insertion sort prepares the pieces of the array for further merging.
     */
    private static final double SORTING_RATIO = 0.1;
    
    private Chip chip;
    private OptimumEmbedding embedder;
    private int[] id;
    private int rows_per_probe;
    
    // each entry contains the sum of 2 pivots, it is most likely that these sum are unique,
    // though it is not a safe, but an easy way!
    private Vector<Integer> visited;
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
	// to do: proceed recursively
    public int makeLayout (Chip chip)
    {
        this.chip = chip;
        id = chip.getMovableProbes();
        final boolean horizontal = true;
        if (chip instanceof AffymetrixChip)
            rows_per_probe = 2;
        else if (chip instanceof SimpleChip)
            rows_per_probe = 1;
        else
            throw new IllegalArgumentException ("Unsupported chip type.");
        
        // reset current layout (if any)
        chip.resetLayout();
        
        // sorting array in order to number of embeddings of probes 
        int pivot_margin = pivotMergeSort(0.05);
        visited = new Vector<Integer>(pivot_margin + 1);
        this.embedder = OptimumEmbedding.createEmbedder(chip, OptimumEmbedding.MODE_CONFLICT_INDEX);
        return partitioning(chip.getChipRegion(),horizontal, 0, pivot_margin, pivot_margin + 1, id.length - 1);
    }
        
    protected int pivotMergeSort(double pivot_treshold)
    {
        int pivot_margin = 0;
        int pivot_property = 0;
        double[] number_of_embeddings = new double[id.length];
        
        for (int i = 0; i < id.length; i++)
        {
            number_of_embeddings[i] = PivotEmbedding.numberOfEmbeddings(chip, id[i]);
        }
        
        synchronousMergeSort(number_of_embeddings, 0, id.length - 1);
        
        while (((pivot_margin + 1) / id.length) < PIVOT_THRESHOLD)
        {
            pivot_property++;
            for (int i = pivot_margin; i < id.length; i++)
            {
                if (number_of_embeddings[i] <= pivot_property)
                {
                    pivot_margin = i;
                }
                else
                {
                    break;
                }
            }
        }
        return pivot_margin;
    }
    
    protected void synchronousMergeSort(double[] property, int start, int stop)
    {
        if (stop - start > SORTING_RATIO * property.length)
        {
            int partition = (int) Math.floor((start + stop)/2);
            synchronousMergeSort(property, start, partition);
            synchronousMergeSort(property, partition + 1 , stop);
            synchronousInsertionSort(property, start, partition);
            synchronousInsertionSort(property, partition + 1, stop);
            synchronousMerge(property, start, partition, stop);
        }
    }
    
    protected void synchronousMerge(double[] property, int start, int partition, int stop)
    {
        double[] left = new double[partition - start + 2];
        double[] right = new double[stop - partition + 1];
        
        System.arraycopy(property, start, left, 0, left.length - 1);
        System.arraycopy(property, partition + 1 , right, 0, right.length - 1);
        left[left.length - 1] = Integer.MAX_VALUE;
        right[right.length - 1] = Integer.MAX_VALUE;
        
        int i = 0;
        int j = 0;
        for (int k = start; k <= stop; k++)
        {
            if (left [i] <= right[j])
            {
                property[k] = left[i];
                switchindex(i, k);
                i++;
            }
            else
            {
                property[k] = right[j];
                switchindex(j + partition + 1, k);
                j++;
            }
        }
    }
    
    protected void synchronousInsertionSort(double[] property, int start, int stop)
    {
    	double property_temp;
    	int id_temp;
    	int i;
    	
    	for(int j = start + 1; j <= stop; j++)
    	{
    		property_temp = property[j];
    		id_temp = id[j];
    		i = j - 1;
    		while(i >= start && property[i] > property_temp)
    		{
    			property[i + 1] = property[i];
    			id[i+1] = id[i];
    			i = i - 1;
    		}
    		property[i + 1] = property_temp;
    		id[i + 1] = id_temp; 
    	}
    }
    
    protected void switchindex(int from_index, int to_index)
    {
        int temp;
        temp = id[from_index];
        id[from_index] = id[to_index];
        id[to_index] = temp;
    }
    // boolean parameter horizontal indicates whether the region is partitioned horizontal or vertical
    protected int partitioning(RectangularRegion region, boolean horizontal, int start_pivot, int stop_pivot, int start_non_pivot, int stop_non_pivot)
    {
        int         cut_pivot;
        int         cut_non_pivot;
        double[]    min_hamming_distance = new double[id.length];
        int         div_rate_non_pivot; 
        int         div_rate_pivot;
        int         rows_first_region;
        int         cols_first_region;
        int         rows_second_region;
        int         cols_second_region;
        int         overflow;
        int         unplaced;
        
        
        // check whether there are enough pivots for further partinioning
        if (stop_pivot - start_pivot < 1)
        {
            return filler.fillRegion(chip, region, id, start_pivot, stop_pivot)
            + filler.fillRegion(chip, region, id, start_non_pivot, stop_non_pivot);    
        }
               
        if (horizontal)
        {
            if (region.last_row - region.first_row + 1 <= stop_dimension * rows_per_probe)
            {
                if (region.last_col - region.first_col + 1 <= stop_dimension)
                {
                    // region cannot be partitioned anymore:
                    // place probes on the specified region and return
                    return filler.fillRegion(chip, region, id, start_pivot, stop_pivot) 
                            + filler.fillRegion(chip, region, id, start_non_pivot, stop_non_pivot);
                }
                else
                {
                    // region can still be vertically partitined
                    return partitioning (region, !horizontal, start_pivot, stop_pivot, 
                                            start_non_pivot, stop_non_pivot);
                }
            }
        }
        
        if (!horizontal)
        {
            if (region.first_col - region.first_col + 1 <= stop_dimension)
            {
                if (region.last_row - region.first_row + 1 <= stop_dimension * rows_per_probe)
                {
                    return filler.fillRegion(chip, region, id, start_pivot, stop_pivot) 
                            + filler.fillRegion(chip, region, id, start_non_pivot, stop_non_pivot);
                }
            }
            else
            {
                return partitioning(region, horizontal, 
                        start_pivot, stop_pivot, start_non_pivot, stop_non_pivot);
            }
        }
            
        
        
        // processing and splitting array
        do
        {
            // find 2 suitable pivots for adjusting the rest of the probes  
            do
            {
                maxMinHammingDistancePivots(start_pivot, stop_pivot, visited);
                // compute hamming distance for non pivots and flag it to which pivot it belongs
                computeDistance(min_hamming_distance, start_pivot, stop_pivot, start_pivot + 1, stop_pivot - 1);
                synchronousMergeSort(min_hamming_distance, start_pivot +1, stop_pivot - 1);
                cut_pivot = getBorder(min_hamming_distance, start_pivot +1, stop_pivot - 1);
                div_rate_pivot = (cut_pivot - start_pivot) / (stop_pivot - start_pivot + 1);
            }
            while (div_rate_pivot < DIV_RATE_PIVOT);
            // compute hamming distance for non pivots and flag it to which pivot it belongs
            computeDistance(min_hamming_distance, start_pivot, stop_pivot, start_non_pivot, stop_non_pivot);
            synchronousMergeSort(min_hamming_distance, start_non_pivot, stop_non_pivot);
            cut_non_pivot = getBorder(min_hamming_distance, start_non_pivot, stop_non_pivot);
            div_rate_non_pivot = (cut_non_pivot - start_non_pivot) / (stop_non_pivot - start_non_pivot + 1);
        }
        while (div_rate_non_pivot < DIV_RATE_NON_PIVOT);
        
        reverseArray(start_pivot + 1, cut_pivot);
        reverseArray(cut_pivot + 1, stop_pivot - 1);
        reverseArray(start_non_pivot, cut_non_pivot);
        reverseArray(cut_non_pivot + 1, stop_non_pivot);
        
        int number_of_probes_first_region = (cut_pivot - start_pivot + cut_non_pivot - start_non_pivot + 2);
        int number_of_probes_second_region = (stop_pivot - cut_pivot + stop_non_pivot - cut_non_pivot);
        
        int div_rate = number_of_probes_first_region / (stop_pivot - start_pivot + stop_non_pivot - start_non_pivot +2);
        // splitting region
        if (horizontal)
        {
           rows_first_region = (int) Math.round(div_rate * (region.last_row - region.first_row + 1));
           cols_first_region = region.last_col - region.first_col + 1;
           rows_second_region = (region.last_row - region.first_row + 1) - rows_first_region;
           cols_second_region = region.last_col - region.first_col + 1;
        }
        else
        {
            rows_first_region = region.last_row - region.first_row + 1;
            cols_first_region =(int) Math.round(div_rate * (region.last_col - region.first_col + 1));
            rows_second_region = region.last_row - region.first_row + 1;
            cols_second_region = (region.last_row - region.first_row +1) - cols_first_region;
        }
        // test if region is suitable to number of probes
        if ((overflow = number_of_probes_first_region - (cols_first_region * rows_first_region) / rows_per_probe) > 0)
        {
            cut_non_pivot -= overflow;
        }
        else if ((overflow = (number_of_probes_second_region) - (cols_second_region * rows_second_region) / rows_per_probe) > 0)
        {
            cut_non_pivot += overflow;
        }
    
        if (horizontal)
        {
            //  create two regions
            int row_div = region.first_row + rows_first_region;
            RectangularRegion first_region = new RectangularRegion(region.first_row, row_div - 1, region.first_col, region.last_col);
            RectangularRegion second_region = new RectangularRegion(row_div, region.last_row, region.first_col, region.last_col);
            //  assign probes to regions 
            unplaced = partitioning(first_region, !horizontal, start_pivot, cut_pivot, start_non_pivot, cut_non_pivot);
            unplaced += partitioning(second_region, !horizontal, cut_pivot + 1, stop_pivot, cut_non_pivot + 1, stop_non_pivot);
        }
        else
        {
            // create two regions
            int col_div = region.first_col + cols_first_region;
            RectangularRegion first_region = new RectangularRegion(region.first_row, region.last_row, region.first_col, col_div - 1);
            RectangularRegion second_region = new RectangularRegion(region.first_row, region.last_row, col_div, region.last_col);
            // assign probes to regions 
            unplaced = partitioning(first_region, horizontal, start_pivot, cut_pivot, start_non_pivot, cut_non_pivot);
            unplaced += partitioning(second_region, horizontal, cut_pivot + 1, stop_pivot, cut_non_pivot + 1, stop_non_pivot);
        }
        
        return unplaced;
    }
    
    private int getBorder(double[] property, int start, int stop)
    {
        int border = start;
        
        for (int i = start; i <= stop; i++)
        {
            if(property[i] < 0)
            {
                border = i;
            }
            else break;
        }
        
        return border;
    }
    
    private void reverseArray(int start, int stop)
    {
        int left;
        int right;
        
        for(int i = start; i <= stop; i++)
        {
            left = i - start;
            right = stop - i;
            if (left > right)
            {
                break;
            }
            else
            {
                switchindex(left, right);
            }
        }
    }
    
    private void computeDistance(double[] min_hamming_distance, int start_pivot, int stop_pivot, int start, int stop)
    {
        double distance_start_pivot;
        double distance_stop_pivot;
        
        
        for (int i = start; i <= stop; i++)
            {
                distance_start_pivot = embedder.minDistanceProbe(id[i], id[start_pivot]);
                distance_stop_pivot = embedder.minDistanceProbe(id[i], id[stop_pivot]);
                if (distance_start_pivot < distance_stop_pivot)
                {
                    min_hamming_distance[i] = -distance_start_pivot;
                }
                else
                {
                    min_hamming_distance[i] = distance_stop_pivot;
                }
            }
        
        
    }
    
    private void maxMinHammingDistancePivots(int start_pivot, int stop_pivot, Vector<Integer> visited)
    {
        double max_distance = -1;
        Integer[] max_distance_pivots = new Integer[2];
        for (int e1 = start_pivot; e1 <= stop_pivot; e1++)
        {
            for (int e2 = start_pivot; e2 <= stop_pivot; e2++)
            {
                if (e1 == e2)
                {
                    continue;
                }
                if (visited.contains(e1 + e2))
                {
                	continue;
                }
                                
                double distance = embedder.minDistanceProbe(id[e1], id[e2]);
                if (max_distance < distance )
                {
                    max_distance = distance;
                    max_distance_pivots[0] = e1;
                    max_distance_pivots[1] = e2;
                }
            }
        }
        if (max_distance_pivots[0] != null && max_distance_pivots[1] != null)
        {
            if(!visited.contains(max_distance_pivots[0] + max_distance_pivots[1]))
            {
            	visited.add(max_distance_pivots[0] + max_distance_pivots[1]);
            }
            
        	switchindex(max_distance_pivots[0], start_pivot);
            switchindex(max_distance_pivots[1], stop_pivot);
        }
        else 
            throw new NullPointerException("No pivots have been found!");
    }
     
}

