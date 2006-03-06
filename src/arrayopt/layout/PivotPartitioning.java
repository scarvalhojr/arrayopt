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
import java.util.Comparator;
import java.util.Arrays;

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
	// to do: proceed recursively
    public int makeLayout (Chip chip)
    {
        Element     elements[];
        int     rows_per_probe;
        int     number_of_embeddings[];
        int     start_pivot = 0;
        int     stop_pivot = 0;
        int     start_non_pivot = 0;
        int     stop_non_pivot = 0;
        int     pivot_margin;
        OptimumEmbedding embedder = OptimumEmbedding.createEmbedder(chip, OptimumEmbedding.MODE_CONFLICT_INDEX);
        
        if (chip instanceof AffymetrixChip)
            rows_per_probe = 2;
        else if (chip instanceof SimpleChip)
            rows_per_probe = 1;
        else
            throw new IllegalArgumentException ("Unsupported chip type.");

        // reset current layout (if any)
        chip.resetLayout();

        // get list of movable probes
        elements = getElements(chip);
        
        Arrays.sort(elements,new NumberOfEmbeddingsComparator());
        pivot_margin = elements[0].number_of_embeddings;
        for (int i = 1; i < elements.length; i++)
        {
            if (elements[i].number_of_embeddings == pivot_margin)
            {
                stop_pivot = i;
            }
            else 
            {
                break;
            }
        }
        start_non_pivot = stop_pivot + 1;
        stop_non_pivot = elements.length -1;
        
        processing(embedder, elements, start_pivot, stop_pivot, start_non_pivot, stop_non_pivot);
       
        
		// should return the number of unplaced probes (if any)
		return 0;
	}
         
    private Element[] getElements(Chip chip)
    {
        int [] ids = chip.getMovableProbes();
        Element[] elements = new Element[ids.length-1];
        for (int i = 0; i < ids.length; i++)
        {
            elements[i] = new Element(chip, ids[i]);
        }
        return elements;
    }
    
    private int processing(OptimumEmbedding embedder, Element[] elements, int start_pivot, int stop_pivot, int start_non_pivot, int stop_non_pivot)
    {
        int cut_pivot;
        int cut_non_pivot;
        
        Integer[] index = getMaxMinHammingDistancePivots(embedder, elements, start_pivot, stop_pivot);
        switchindex(elements, index, start_pivot, stop_pivot);
        computeDistance(embedder, elements, start_pivot, stop_pivot, start_pivot + 1, stop_pivot - 1);
        computeDistance(embedder, elements, start_pivot, stop_pivot, start_non_pivot, stop_non_pivot);
        cut_pivot = reorderToDistance(elements, start_pivot, stop_pivot, start_pivot + 1, stop_pivot - 1);
        cut_non_pivot = reorderToDistance(elements, start_pivot, stop_pivot, start_non_pivot, stop_non_pivot);
        
        processing(embedder, elements, start_pivot, cut_pivot - 1, start_non_pivot, cut_non_pivot - 1);
        processing(embedder, elements, cut_pivot, stop_pivot, cut_non_pivot, stop_non_pivot);
        
        return 0;
    }
    
    
    private int reorderToDistance(Element[] elements, int start_pivot, int stop_pivot, int start, int stop)
    {
        int border =start;
        for (int i = start; i <= stop; i++)
        {
            if (elements[i].min_hamming_distance > 0)
            {
                border = i;
                break;
            }
        }
        
        Arrays.sort(elements, start, border, new HammingDistanceComparator());
        Arrays.sort(elements, border, stop+1, new HammingDistanceComparator());
        return border;
    }
    
    
    private void computeDistance(OptimumEmbedding embedder,Element[] elements, int start_pivot, int stop_pivot, int start, int stop)
    {
        double distance_start;
        double distance_stop;
        for (int i = start; i <= stop; i++)
        {
            distance_start = embedder.minDistanceProbe(elements[i].id, elements[start_pivot].id);
            distance_stop = embedder.minDistanceProbe(elements[i].id, elements[stop_pivot].id);
            if (distance_start < distance_stop)
            {
                elements[i].min_hamming_distance = -distance_start;
            }
            else
            {
                elements[i].min_hamming_distance = distance_stop;
            }
        }
        
    }
    
    private void switchindex(Element[] elements, Integer[] index, int start, int stop)
    {
        switchindex(elements, index[0], start);
        switchindex(elements, index[1], stop);
    }
    
    private void switchindex(Element[] elements, int from_index, int to_index)
    {
        Element temp = elements[from_index];
        elements[from_index] = elements[to_index];
        elements[to_index] = temp;
    }
    
    private Integer[] getMaxMinHammingDistancePivots(OptimumEmbedding embedder, Element[] elements, int start_pivot, int stop_pivot)
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
                double distance = embedder.minDistanceProbe(elements[e1].id, elements[e2].id);
                if (max_distance < distance )
                {
                    max_distance = distance;
                    max_distance_pivots[0] = elements[e1].id;
                    max_distance_pivots[1] = elements[e2].id;
                }
            }
        }
        if (max_distance_pivots[0] != null && max_distance_pivots[1] != null)
        {
            return max_distance_pivots;
        }
        else 
            throw new NullPointerException("No pivots have been found!");
    }
    
    private class Element
    {
        int id;
        int number_of_embeddings;
        // used for holding the hamming distance to the pivots
        double min_hamming_distance = 0;
        
        private Element(Chip chip, int id)
        {
            this.id = id;
            number_of_embeddings = PivotEmbedding.numberOfEmbeddings(chip,id);
        }
    }
    
    private class NumberOfEmbeddingsComparator implements Comparator<Element>
    {
        public int compare(Element e1, Element e2)
        {
            return (((Integer) e2.number_of_embeddings).compareTo((Integer) e1.number_of_embeddings));
        }
    }
    
    private class HammingDistanceComparator implements Comparator<Element>
    {
        public int compare(Element e1, Element e2)
        {
            return (((Double) e1.min_hamming_distance).compareTo( e2.min_hamming_distance));
        }
    }
}

