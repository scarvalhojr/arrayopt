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
import java.util.BitSet;
import java.util.PriorityQueue;

/**
 * document this 
 */
public class PivotOptimization implements PostPlacementAlgorithm
{
    protected Chip chip;
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
        CompareProbe comparator = new CompareProbe(chip, chip.getFixedSpots());
        int num_col = chip.getNumberOfColumns();
        int num_row = chip.getNumberOfRows();
        
        PriorityQueue<Element> Queue = new PriorityQueue<Element>(chip.getNumberOfProbes(), comparator);
        
        
        for (int r = 0; r < num_row; r++)
        {
            for (int c = 0; c < num_col; c++)
            {
                if (chip.spot[r][c] == Chip.EMPTY_SPOT)
                {
                    comparator.finished.set(r * num_col + c);
                }
                else if (PivotEmbedding.numberOfEmbeddings(chip, chip.spot[r][c]) == 1 )
                {
                    comparator.finished.set(r * num_col + c);
                    
                    for (int i = r-1; i <= r+1; i++)
                    {
                        for (int j = c - 1; j <= c + 1; j++) 
                        {
                            if (j == c || j < 0 || i < 0 || j > num_col || i > num_row)
                            {
                                break;
                            }
                            
                            Element add = new Element(i,j,chip.spot[i][j]);
                            Queue.add(add);
                        }
                    }
                    
                    
                    
                }
            }
        }
        
        while(comparator.finished.cardinality() != comparator.finished.length())
        {
            Element current = Queue.poll();
            OptimumEmbedding.reembedProbe(chip, current.id, current.getImmediateNeighbors(comparator.finished));
            comparator.finished.set(current.row * num_col + current.column);
        }
        
    }
    
    
    
    void optimizeLayout(AffymetrixChip chip)
    {
        
    }
    
    private class CompareProbe implements Comparator<Element>
    {
        private Chip chip;
        protected BitSet finished;
        private final int region_dimension = 3;
        
        public CompareProbe(Chip chip, BitSet finished)
        {
            this.chip = chip;
            this.finished = finished;
      
        }
        
        public int  compare(Element first,Element second)
        {
            int compareproperties = ((Integer) (PivotEmbedding.numberOfEmbeddings(this.chip, first.id))).compareTo(PivotEmbedding.numberOfEmbeddings(this.chip, second.id));
            
            if (compareproperties == 0)
            {
                compareproperties = ((Integer) first.getNumberOfImmediateNeighbors(finished)).compareTo((Integer) (second.getNumberOfImmediateNeighbors(finished)));                    
                
                if (compareproperties == 0)
                {
                    compareproperties = ((Integer) first.getNumberOfRegionNeighbors(finished, region_dimension)).compareTo((Integer) second.getNumberOfRegionNeighbors(finished, region_dimension));
                    
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
        public  int id;
              
        protected Element(int row, int column, int id)
        {
            this.row = row;
            this.column = column;
            this.id = id;
        }
        
        public int getNumberOfImmediateNeighbors(BitSet finished)
        {
            int neighbors = 0;
            for (int i=-1; i <= +1; i++)
            {
                for (int j=-1; j <= +1; j++)
                {
                    if (!(i == 0 && j == 0) && (i == 0 || j == 0))
                        
                    {
                        
                        if(finished.get((row+i)*chip.num_cols+(column+j)))
                        {
                            neighbors++;
                        }
                        
                    }
                }
            }
            
            return neighbors;
        }
        
        public int getNumberOfRegionNeighbors(BitSet finished, int dimension)
        {
            int neighbors = 0;
            for (int i=-dimension; i <= +dimension; i++)
            {
                for (int j=-dimension; j <= +dimension; j++)
                {
                    
                    if (!(i == 0 && j == 0))                     
                    {
                        
                        if(finished.get((row+i)*chip.num_cols+(column+j)))
                        {
                            neighbors++;
                        }
                        
                    }
                }
            }
            return neighbors;
        }
        
        public int[] getImmediateNeighbors(BitSet finished)
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
        }
    }
    
}


