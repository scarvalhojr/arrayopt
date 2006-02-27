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

/**
 * document this 
 */
public class PivotOptimization implements PostPlacementAlgorithm
{
    private Chip chip;
    private Chip optimized_chip;
    private Integer embedding_mode;
    private final int default_mode = OptimumEmbedding.MODE_CONFLICT_INDEX;
    /**
     * document this
     */
    public void optimizeLayout (Chip chip, int mode)
    {
           embedding_mode = mode;
           optimizeLayout(chip);
    }
    
    public void optimizeLayout (Chip chip)
    {
        if (embedding_mode == null)
        {
            embedding_mode = default_mode;
        }
        
        if (chip instanceof SimpleChip)
            optimizeLayout ((SimpleChip) chip);
        
        else if (chip instanceof AffymetrixChip)
            optimizeLayout ((AffymetrixChip) chip);
        
        else
            throw new IllegalArgumentException ("Unsupported chip type.");
    }
    
    
    
    void optimizeLayout(SimpleChip chip)
    {
        this.chip = chip;
        optimized_chip = chip.clone();
        embedding_mode = OptimumEmbedding.MODE_CONFLICT_INDEX;
        OptimumEmbedding embedder = OptimumEmbedding.createEmbedder(chip, embedding_mode);
        CompareProbe comparator = new CompareProbe();
                
        //reset all spots on to be optimized chip to empty spots
        clearChip(optimized_chip);
                               
        PriorityQueue<Element> queue = new PriorityQueue<Element>(chip.getNumberOfProbes(), comparator);
        
              
        for (int r = 0; r < chip.num_rows; r++)
        {
            for (int c = 0; c < chip.num_cols; c++)
            {
                if ((PivotEmbedding.numberOfEmbeddings(chip, chip.spot[r][c])) == 1 )
                {
                    optimized_chip.spot[r][c] = chip.spot[r][c];
                    addNeighbors(queue, new Element(r,c));
                }
            }
            
            
            
        }
    

        
        while(queue.size() != 0)
        {
            Element current = queue.poll();
            optimized_chip.spot[current.row][current.column] = chip.spot[current.row][current.column];
            embedder.reembedSpot(current.row, current.column);
            addNeighbors(queue, current);            
        }
        
        if (!chip.compatible(optimized_chip))
        {
            System.err.println("Warning: optimization of chip failed!");
        }
    }
    
    private void addNeighbors(PriorityQueue<Element> queue, Element element)
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
                        queue.add(new Element(r, c));
                    }
                    
                }
            }
            if (chip instanceof AffymetrixChip)
            {
                i++;
            }
        }
       
    }
    
    private void clearChip(Chip chip)
    {
        for (int r = 0; r < chip.num_rows; r++)
        {
            for (int c = 0; c < chip.num_cols; c++)
            {
                chip.spot[r][c] = Chip.EMPTY_SPOT;
            }
            if (chip instanceof AffymetrixChip)
            {
                r++;
            }
        }
    }
    
    void optimizeLayout(AffymetrixChip chip)
    {
        this.chip = chip;
        optimized_chip = chip.clone();
        OptimumEmbedding embedder = OptimumEmbedding.createEmbedder(chip, OptimumEmbedding.MODE_CONFLICT_INDEX);
        CompareProbe comparator = new CompareProbe();
        
        clearChip(optimized_chip);
        
        PriorityQueue<Element> queue = new PriorityQueue<Element>(chip.getNumberOfProbes(), comparator);
        
        for (int r = 0; r < chip.num_rows; r += 2)
        {
            for (int c = 0; c < chip.num_cols; c++)
            {
                if ((PivotEmbedding.numberOfEmbeddings(chip, chip.spot[r][c]) == 1 ))
                {
                    optimized_chip.spot[r][c] = chip.spot[r][c];
                    addNeighbors(queue, new Element(r,c));
                }
            }
            
            
            
        }
    

        
        while(queue.size() != 0)
        {
            Element current = queue.poll();
            optimized_chip.spot[current.row][current.column] = chip.spot[current.row][current.column];
            embedder.reembedSpot(current.row, current.column);
            addNeighbors(queue, current);            
        }
        
        if (!chip.compatible(optimized_chip))
        {
            System.err.println("Warning: optimization of chip failed!");
        }
        
    }
         
    
    
    private class CompareProbe implements Comparator<Element>
    {
             
        private int  compare(Element first,Element second)
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


