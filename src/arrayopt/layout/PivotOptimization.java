/*
 * PivotOptimization.java
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
 * ?rrayOpt is distributed in the hope that it will be useful, but WITHOUT ANY
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
import java.util.PriorityQueue;
import java.util.Iterator;
import java.util.Stack;

/**
 * PivotOptimization just reembedds the probes of a chip according to its pivots (normally probes with only 1 possible embedding).
 * If no real pivot can be found, the probes with a minimum of embeddings are handled as pivots.
 * Once the pivots are found, all neighbors around it will be reembedded according to it and after that alle neighbors of the already re
 * embedded probes will be reembedded according to the {@link OptimumEmbedding} and so forth.
 * The probes will not change their spot on the chip.
 * 
 * <P>PivotOptimization is considered for use with simple chips and affymetrix chips.
 * Since it uses the {@link OptimumEmbedding} class for reembedding the probes you can also define which mode should be used for re
 * embedding the probes, i.e.{@link OptimumEmbedding#MODE_CONFLICT_INDEX} or {@link OptimumEmbedding#MODE_BORDER_LENGTH}. If no mode is given the
 * {@link OptimumEmbedding#MODE_CONFLICT_INDEX} is used.</P>
 * 
 * @author Anna Domanski & Ronny Gaertner
 * 
 */
public class PivotOptimization implements PostPlacementAlgorithm
{
    private Chip chip;
    private Chip optimized_chip;
    private Integer embedding_mode;
    private final int default_mode = OptimumEmbedding.MODE_CONFLICT_INDEX;
    private int min_number_of_embeddings = Integer.MAX_VALUE;
    private PriorityQueue<Element> queue;
    
    /**
     * @param chip chip instance 
     * @param mode reembedding mode (either optimization for conflict index or border length)
     */
    public void optimizeLayout (Chip chip, int mode)
    {
           embedding_mode = mode;
           optimizeLayout(chip);
    }
    /**
     * @param
     */
    public void optimizeLayout (Chip chip)
    {
        if (embedding_mode == null)
        {
            embedding_mode = default_mode;
        }
        
        if (chip instanceof SimpleChip);
        else if (chip instanceof AffymetrixChip);
        else
            throw new IllegalArgumentException ("Unsupported chip type.");
    
        this.chip = chip;
        optimized_chip = chip.clone();
        embedding_mode = OptimumEmbedding.MODE_CONFLICT_INDEX;
        CompareProbe comparator = new CompareProbe();
               
        initializeChipOptimization();
        queue = new PriorityQueue<Element>(chip.getNumberOfProbes(), comparator);
        pivotScanning();
        processing();
        missedProbesScanning();
        processing();
        
    }
    
    private void addNeighbors(Element element)
    {
        
        int r;
        int c;
        
        int startrow = -1;
        int stoprow = +1;
        
        if (chip instanceof AffymetrixChip)
        {
            startrow = -2;
            stoprow = +2;
        }
        
        for (int i = startrow; i <= stoprow; i++)
        {
            for (int j = -1; j <= +1; j++)
            {
                if (!(i == 0 && j == 0) && (i == 0 || j == 0))
                {
                    if( ((r = element.row + i) >= 0) && ( r < chip.num_rows) && ((c = element.column + j) >= 0) && c < chip.num_cols)
                    {
                        if (optimized_chip.spot[r][c] != Chip.EMPTY_SPOT || chip.spot[r][c] == Chip.EMPTY_SPOT)
                        {
                            break;
                        }
                        this.queue.add(new Element(r, c));
                    }
                    
                }
            }
            if (chip instanceof AffymetrixChip)
            {
                i++;
            }
        }
    }
      
    private void pivotScanning() 
    {
       int number_of_embeddings;
        
        for (int r = 0; r < chip.num_rows; r++)
        {
            for (int c = 0; c < chip.num_cols; c++)
            {
                number_of_embeddings = PivotEmbedding.numberOfEmbeddings(chip, chip.spot[r][c]);
                if (number_of_embeddings == min_number_of_embeddings)
                {
                    optimized_chip.spot[r][c] = chip.spot[r][c];
                    addNeighbors(new Element(r,c));
                }
            }
            if (chip instanceof AffymetrixChip)
            {
                r++;
            }
                       
        }
    }
    
    // all elements in the collected in the priority queue will be processed. The queue
    // will be updated for each new embedding of a probe
    private void processing()
    {
        Stack<Element> updated_probes = new Stack<Element>();
        OptimumEmbedding embedder = OptimumEmbedding.createEmbedder(chip, embedding_mode);
        
        
        while(queue.size() != 0)
        {
            
            Element current = queue.poll();
            optimized_chip.spot[current.row][current.column] = chip.spot[current.row][current.column];
            Iterator<Element> queue_iterator = queue.iterator();
            
            // 2 while loops are for updating the priority queue considering a probe was
            // reembedded
            while (queue_iterator.hasNext())
            {
                Element element = queue_iterator.next();
                if (element.updateNeighbors())
                {
                    updated_probes.push(element);
                    queue_iterator.remove();
                }
            }
            while (!updated_probes.empty())
            {
                Element element = updated_probes.pop();
                // queue.remove(element);
                // element.updateNeighbors(true);
                queue.add(element);
            }
            
            embedder.reembedSpot(current.row, current.column);
            addNeighbors(current);
            
        }
    }
      
    private void missedProbesScanning()
    {
        for (int r = 0; r < chip.num_rows; r++)
        {
            for (int c = 0; c < chip.num_cols; c++)
            {
                if((chip.spot[r][c] != Chip.EMPTY_SPOT) && (optimized_chip.spot[r][c] == Chip.EMPTY_SPOT))
                {
                    queue.add(new Element(r,c));
                }
            }
            if (chip instanceof AffymetrixChip)
            {
                r++;
            }
        }
    }
//  reset all spots on to be optimized chip to empty spots and fine the minimum number of embeddings
    private void initializeChipOptimization()
    {
        int number_of_embeddings;
        for (int r = 0; r < chip.num_rows; r++)
        {
            for (int c = 0; c < chip.num_cols; c++)
            {
                number_of_embeddings = PivotEmbedding.numberOfEmbeddings(chip, chip.spot[r][c]);
                if (min_number_of_embeddings > number_of_embeddings)
                        {
                            min_number_of_embeddings = number_of_embeddings;
                        }
                optimized_chip.spot[r][c] = Chip.EMPTY_SPOT;
            }
            if (chip instanceof AffymetrixChip)
            {
                r++;
            }
        }
    }
       
    private class CompareProbe implements Comparator<Element>
    {
             
        public int  compare(Element first,Element second)
        {
            int compareproperties = ((Integer) first.num_embed).compareTo(second.num_embed);
            
            if (compareproperties == 0)
            {
                compareproperties = ((Integer) first.immediate_neighbors).compareTo( second.immediate_neighbors);                    
                
                if (compareproperties == 0)
                {
                    compareproperties = ((Integer) first.region_neighbors).compareTo(second.region_neighbors);
                    
                    // if still no decision could be made set it to the second probe! <-- could also be the first one, but a decision must be achieved!
                    if (compareproperties == 0)
                    {
                        compareproperties = 1;
                    }
                    
                }
            }
            return ((int) compareproperties);
        }
    }
    
    private class Element 
    {
        private int row;
        private int column;
        private int id;
        private int num_embed;
        private int immediate_neighbors;
        private int region_neighbors;
        private final int REGION_DIMENSION = 3;
        
        protected Element(int row, int column)
        {
            this.row = row;
            this.column = column;
            id = chip.spot[row][column];
            num_embed = PivotEmbedding.numberOfEmbeddings(chip, id);
        }
       
       
        private boolean updateNeighbors()
        {
            int immediate = getNumberOfImmediateNeighbors();
            int region = getNumberOfRegionNeighbors(REGION_DIMENSION);
  
            if (immediate_neighbors != immediate || region_neighbors != region)
            {
                immediate_neighbors = immediate;
                region_neighbors = region;
                return true;
            }
            else
            {
            	return false;
            }
        }
         
        private int getNumberOfImmediateNeighbors()
        {
            int neighbors = 0;
            int startrow = -1;
            int startcolumn = -1;
            int stoprow = +1;
            int stopcolumn = +1;
            
            if (row == 0)
            {
                startrow = row;
            }
            else if (row == (optimized_chip.num_rows-1))
            {
                stoprow = 0;
            }
            
            if (column == 0)
            {
                startcolumn = column;
            }
            else if (column == optimized_chip.num_cols-1)
            {
                stopcolumn = 0;
            }
            
            for (int i= startrow; i <= stoprow; i++)
            {
                for (int j = startcolumn; j <= stopcolumn; j++)
                {
                    if (!(i == 0 && j == 0) && (i == 0 || j == 0))
                        
                    {
                        if(optimized_chip.spot[row+i][column+j] != Chip.EMPTY_SPOT)
                        {
                            neighbors++;
                        }
                        
                    }
                }
                if ( chip instanceof AffymetrixChip)
                {
                    i++;
                }
            }
            
            return neighbors;
        }
       
        private int getNumberOfRegionNeighbors(int dimension)
        {
            int neighbors = 0;
            int startrow = -dimension; 
            int startcolumn = -dimension;
            int stoprow = +dimension;
            int stopcolumn = +dimension;
            int diff;
            
            if ((diff = row-dimension) < 0)
            {
                startrow =  dimension - diff;
            }
            if ((diff = (optimized_chip.num_rows-1) - row + dimension) < 0)
            {
                stoprow = dimension + diff;
            }
            if ((diff = column-dimension) < 0)
            {
                startcolumn = dimension - diff;
            }
            if ((diff = (optimized_chip.num_cols-1) - column + dimension) < 0)
            {
                stopcolumn = dimension + diff;
            }
            
            for (int i=startrow; i <= stoprow; i++)
            {
                for (int j=startcolumn; j <= stopcolumn; j++)
                {
                    
                    if (!(i == 0 && j == 0))                     
                    {
                        {
                            neighbors++;
                        }
                        
                    }
                }
                if (chip instanceof AffymetrixChip)
                {
                    i++;
                }
            }
            return neighbors;
        }
    }
    
}


