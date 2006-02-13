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
import java.util.PriorityQueue;

/**
 * document this 
 */
public class PivotOptimization implements PostPlacementAlgorithm
{
    private Chip chip;
    private Chip optimized_chip;
    
    /**
     * document this
     */
    public void optimizeLayout (Chip chip)
    {
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
        OptimumEmbedding embedder = OptimumEmbedding.createEmbedder(optimized_chip, OptimumEmbedding.MODE_CONFLICT_INDEX);
        CompareProbe comparator = new CompareProbe();
                
        //reset all spots on to be optimized chip
        for (int r = 0; r < optimized_chip.num_rows; r++)
            for (int c = 0; c < optimized_chip.num_cols; c++)
            {
                optimized_chip.spot[r][c] = Chip.EMPTY_SPOT;
            }
               
        PriorityQueue<Element> queue = new PriorityQueue<Element>(chip.getNumberOfProbes(), comparator);
        
        addNeighbors(queue, new Element(2,2));
        
        for (int r = 0, numberofembeddings = 0; r < chip.num_rows; r++)
        {
            for (int c = 0; c < chip.num_cols; c++)
            {
                /*if (chip.spot[r][c] == Chip.EMPTY_SPOT)
                {
                    optimized_chip.spot[r][c] = chip.spot[r][c];
                }*/
                if ((numberofembeddings = PivotEmbedding.numberOfEmbeddings(chip, chip.spot[r][c])) == 1 )
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
        else
        {
            this.chip =  optimized_chip;
        }
    }
    
    private void addNeighbors(PriorityQueue<Element> queue, Element element)
    {
        int r;
        int c;
        for (int i = -1; i <= +1; i++)
        {
            for (int j = -1; j <= +1; j++)
            {
                if (!(i == 0 && j == 0) && (i == 0 || j == 0))
                {
                    if( ((r = element.row + i) > 0) && ( r < chip.num_rows) && ((c = element.column + j) > 0) && c < chip.num_cols)
                    {
                        if (optimized_chip.spot[r][c] != Chip.EMPTY_SPOT || chip.spot[r][c] == Chip.EMPTY_SPOT)
                        {
                            break;
                        }
                        queue.add(new Element(r, c));
                    }
                    
                }
            }
        }
        
        
        /* Element neighbor;
        if (neighbor = element.getLeftNeighbor != null)
        {
            queue.add(neighbor);
        }
        if (neighbor = element.getRightNeighbor != null)
        {
            queue.add(neighbor);
        }
        if (neighbor = element.getAboveNeighbor != null)
        {
            queue.add(neighbor);
        }
        if (neighbor = element.getBelowNeighbor != null)
        {
            queue.add(neighbor);
        }*/
    }
    
    
    
    void optimizeLayout(AffymetrixChip chip)
    {
        this.chip = chip;
        optimized_chip = chip.clone();
        OptimumEmbedding embedder = OptimumEmbedding.createEmbedder(optimized_chip, OptimumEmbedding.MODE_CONFLICT_INDEX);
        CompareProbe comparator = new CompareProbe();
        
        for (int r = 0; r < optimized_chip.num_rows; r++)
            for (int c = 0; c < optimized_chip.num_cols; c++)
            {
                optimized_chip.spot[r][c] = Chip.EMPTY_SPOT;
            }
    }
    
    private class CompareProbe implements Comparator<Element>
    {
        private final int REGION_DIMENSION = 3;
        
               
        public int  compare(Element first,Element second)
        {
            int compareproperties = ((Integer) (PivotEmbedding.numberOfEmbeddings(chip, first.id))).compareTo(PivotEmbedding.numberOfEmbeddings(chip, second.id));
            
            if (compareproperties == 0)
            {
                compareproperties = ((Integer) first.getNumberOfImmediateNeighbors()).compareTo((Integer) (second.getNumberOfImmediateNeighbors()));                    
                
                if (compareproperties == 0)
                {
                    compareproperties = ((Integer) first.getNumberOfRegionNeighbors(REGION_DIMENSION)).compareTo((Integer) second.getNumberOfRegionNeighbors(REGION_DIMENSION));
                    
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
    
    protected class Element 
    {
        protected  int row;
        protected  int column;
        protected  int id;
        
        protected Element(int row, int column)
        {
            this.row = row;
            this.column = column;
            id = chip.spot[row][column];
        }
       
       
        
        // not necessary anymore is the same as getImmediateNeighbors(finished).length 
        public int getNumberOfImmediateNeighbors()
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
            }
            
            return neighbors;
        }
       
        public int getNumberOfRegionNeighbors(int dimension)
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
                stoprow = dimension - diff;
            }
            if ((diff = column-dimension) < 0)
            {
                startcolumn = dimension - diff;
            }
            if ((diff = (optimized_chip.num_cols-1) - column + dimension) < 0)
            
            for (int i=startrow; i <= stoprow; i++)
            {
                for (int j=startcolumn; j <= stopcolumn; j++)
                {
                    
                    if (!(i == 0 && j == 0))                     
                    {
                        
                        if(optimized_chip.spot[row+i][column+j] != Chip.EMPTY_SPOT)
                        {
                            neighbors++;
                        }
                        
                    }
                }
            }
            return neighbors;
        }
        
      /*  public int[] getImmediateNeighbors(BitSet finished)
        {
            int[] neighbors = new int[3];
            int counter = 0;
            
            for (int i=-1; i <= +1; i++)
            {
                for (int j=-1; j <= +1; j++)
                {
                    if (!(i == 0 && j == 0) && (i == 0 || j == 0))
                        
                    {
                        
                        if(finished.get((row+i)*chip.num_cols+(column+j)))
                        {
                            neighbors[counter] = chip.spot[i][j];
                            counter++;
                        }
                        
                    }
                }
            }
        int[] neighbors_opt = new int[counter];
        System.arraycopy(neighbors,0,neighbors_opt,0,counter);
        return neighbors_opt;
        } */
    }
    
}


