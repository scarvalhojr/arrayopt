/*
 * ConflictIndex.java
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

/**
 * This class gives the exact definition of the concept of conflict index. Its
 * main purpose is to define the values of the distance-dependent
 * ({@link #distanceWeight(int, int, int, int)}) and postition-dependent
 * ({@link #positionWeight(int, int)}) weighting functions. The latter is
 * usually implemented as a pre-computed two-dimensional array whose dimensions
 * (i.e. its range of values greater than zero) can be queried with the
 * {@link #dimConflictRegion()} method.</P>
 * 
 * <P>Note that this class is abstract and that it also contains several
 * inner-sub-classes that give different definitions of a conflict index. Only
 * one of them is (statically) loaded in this class at any given time. The
 * loaded definition is used by all classes of the layout package.</P>
 * 
 * <P>The default definition ({@link #DEFAULT_DEFINITION}) is a direct
 * implementation of the definition given in the following paper:<BR>
 * 
 * TODO add reference to the paper with the definition of conflict index.<P>
 * 
 * <P>The sub-classes cannot be directly instantiated. Rather, the
 * {@link #loadDefinition} method should be used with the public constants
 * corresponding to the desired definition. See the documentation of each
 * constant for a precise meaning of the available definitions.</P>   
 * 
 * @author Sergio A. de Carvalho Jr.
 */
public abstract class ConflictIndex
{
	/**
	 * This constants represents the default definition of conflict index as
	 * defined in the paper:
	 * 
	 * TODO add reference to the paper with the definition of conflict index.
	 */
	public static final int DEFAULT_DEFINITION = 0;
	
	/**
	 * This constant represents a simplified conflict index defintion. The
	 * weighted distance matrix is balanced and has only three different values:
	 * 0, 0.5 and 1. Only closer spots have non-zero distance weight. The
	 * position-dependent weights are also "simplified" by mapping the default
	 * values to integer values.
	 */
	public static final int SIMPLIFIED_DEFINITION = 1;
	
	/**
	 * This constant represents a conflict index definition that is equivalent
	 * to the definition of border length. That implies that the distance
	 * weight equals 1 for immediate neighbors or zero otherwise.
	 * Position-dependent weights are fixed to 1, independent of the base.
	 */
	public static final int BORDER_LENGTH_EQUIVALENT_DEFINITION = 2;
	
	/**
	 * This constant represents an unbalanced conflict index definition. The
	 * distance weights have no logic or symmetry whatsoever. Position-dependent
	 * weights are also skewed by giving higher values to conflicts in the
	 * beginning of the probe. This definition is mainly used for debugging.
	 */
	public static final int UNBALANCED_DEFINITION = 3;
	
	/**
	 * The currently loaded definition of conflict index.
	 */
	private static ConflictIndex loaded_def = new Default ();

	/**
	 * Loads a pre-defined conflict index definition. Each available defintion
	 * have a corresponding public constant in this class.
	 *   
	 * @param def constant representing the desidered definition
	 */
	public static void loadDefinition (int def)
	{
		switch (def)
		{
			case DEFAULT_DEFINITION:
				loaded_def = new Default ();
				break;
			
			case SIMPLIFIED_DEFINITION:
				loaded_def = new Simplified ();
				break;
			
			case BORDER_LENGTH_EQUIVALENT_DEFINITION:
				loaded_def = new BorderLengthEquivalent ();
				break;

			case UNBALANCED_DEFINITION:
				loaded_def = new Unbalanced ();
				break;

			default:
				throw new IllegalArgumentException
					("Unknown conflict index definition: " + def);
		}
	}
	
	/**
	 * Returns the dimension of the weighting distance matrix. This dimension is
	 * the number of rows and columns around a spot that are considered in the
	 * conflict index computation. If it is 3, for instance, up to three rows
	 * (and columns) up or down (to the left or to the right) of the spot are
	 * considered, resulting in a 7x7 region around the spot.
	 *     
	 * @return dimension of weighted distance matrix
	 */
	public static int dimConflictRegion ()
	{
		return loaded_def.dimConflictRegion_internal ();
	}
	
	/**
	 * Returns the dimension of the weighted distance matrix. This is an
	 * internal method that sub-classes must provide.
	 */
	abstract int dimConflictRegion_internal ();

	/**
	 * Returns the distance-dependent weight of a spot (r2,c2) in regards to a
	 * spot (r1,c1) according to this conflict index definition. Each spots is
	 * given as a coordinate pair (row, column). Note that the caller must
	 * ensure that spot (r2,c2) is inside the conflict region of (r1,c1). This
	 * region's dimension can be queried with the {@link #dimConflictRegion()} 
	 * method.
	 * 
	 * @param r1 row coordinate of the first spot
	 * @param c1 column coordinate of the first spot
	 * @param r2 row coordinate of the second spot
	 * @param c2 column coordinate of the second spot
	 * @return the distance weight of first spot in regards to the second
	 */
	public static double distanceWeight (int r1, int c1, int r2, int c2)
	{
		return loaded_def.distanceWeight_internal(r1, c1, r2, c2);
	}
	
	/**
	 * Returns the distance-dependent weight of a spot (r2,c2) in regards to a
	 * spot (r1,c1) according to this conflict index definition. This is an
	 * internal method that sub-classes must provide.
	 */ 
	abstract double distanceWeight_internal (int r1, int c1, int r2,int c2);
	
	/**
	 * Returns the position-dependent weight of a base according to this
	 * conflict index definition. Usually this weight is determined by the base
	 * number and the probe length. The exact value is, however, defined by
	 * each definition differently (check the the corresponding public constant
	 * representing the definition).
	 *  
	 * @param base number of the base where a conflict would occurr
	 * @param probe_len length of the probe
	 * @return position-dependent wieight
	 */
	public static double positionWeight (int base, int probe_len)
	{
		return loaded_def.positionWeight_internal (base, probe_len);
	}
	
	/**
	 * Returns the position-dependent weight of a base according to this
	 * conflict index definition.  This is an internal method that sub-classes
	 * must provide.
	 */
	abstract double positionWeight_internal (int base, int probe_len);
	
	/**
	 * Default conflict index defintion.
	 */
	private static class Default extends ConflictIndex
	{
		private final int DIM = 3;
		
		private final double THETA_NUM = 5;
		
		private double weight_dist[][];
		
		Default ()
		{
			int size, h, v;
			double d2;
			
			size = 2 * DIM + 1;
			
			weight_dist = new double[size][size];
			
			for (int r = 0; r < size; r++)
			{
				v = DIM - r;
				
				for (int c = 0; c < size; c++)
				{
					h = DIM - c;
					
					d2 = v * v + h * h;
					
					weight_dist[r][c] = d2 > 0 ? 1 / d2 : 0;
				}
			}
		}
		
		@Override
		int dimConflictRegion_internal ()
		{
			return DIM;
		}
		
		@Override
		double positionWeight_internal (int base, int probe_len)
		{
			double theta, c, lambda;
			
			theta = THETA_NUM / probe_len;
			
			c = 1 / Math.exp(theta);
			
			lambda = (base <= probe_len - base) ?
						(base + 1) : (probe_len - base + 1);
			
			return c * Math.exp(theta * lambda);
		}
		
		@Override
		double distanceWeight_internal (int r1, int c1, int r2, int c2)
		{
			return weight_dist[DIM + r2 - r1][DIM + c2 - c1];
		}
	}

	/**
	 * Simplified conflict index defintion.
	 */
	private static class Simplified extends ConflictIndex
	{
		private final int DIM = 3;
		
		private final double THETA_NUM = 5;
		
		private double weight_dist[][] = {
				{	0,	0,		0,		0,		0,		0,		0	},
				{	0,	0,		0,		0.1,	0,		0,		0	},
				{	0,	0,		0.5,	1,		0.5,	0,		0	},
				{	0,	0.1,	1,		0,		1,		0.1,	0	},
				{	0,	0,		0.5,	1,		0.5,	0,		0	},
				{	0,	0,		0,		0.1,	0,		0,		0	},
				{	0,	0,		0,		0,		0,		0,		0	}};
				
		@Override
		int dimConflictRegion_internal ()
		{
			return DIM;
		}
		
		@Override
		double positionWeight_internal (int base, int probe_len)
		{
			double theta, c, lambda;
			
			theta = THETA_NUM / probe_len;
			
			c = 1 / Math.exp(theta);
			
			lambda = (base <= probe_len - base) ?
						(base + 1) : (probe_len - base + 1);
			
			return (int) (c * Math.exp(theta * lambda));
		}
		
		@Override
		double distanceWeight_internal (int r1, int c1, int r2, int c2)
		{
			return weight_dist[DIM + r2 - r1][DIM + c2 - c1];
		}
	}

	/**
	 * Conflict index defintion equivalent to the concept of border length.
	 */
	private static class BorderLengthEquivalent extends ConflictIndex
	{
		private final int DIM = 3;
		
		private double weight_dist[][] = {
				{	0,	0,	0,	0,	0,	0,	0	},
				{	0,	0,	0,	0,	0,	0,	0	},
				{	0,	0,	0,	1,	0,	0,	0	},
				{	0,	0,	1,	0,	1,	0,	0	},
				{	0,	0,	0,	1,	0,	0,	0	},
				{	0,	0,	0,	0,	0,	0,	0	},
				{	0,	0,	0,	0,	0,	0,	0	}};
		
		@Override
		int dimConflictRegion_internal ()
		{
			return DIM;
		}
		
		@Override
		double positionWeight_internal (int base, int probe_len)
		{
			// always return 1
			return 1 + 0 * (base + probe_len);
		}
		
		@Override
		double distanceWeight_internal (int r1, int c1, int r2, int c2)
		{
			return weight_dist[DIM + r2 - r1][DIM + c2 - c1];
		}
	}

	/**
	 * Unbalanced conflict index defintion.
	 */
	private static class Unbalanced extends ConflictIndex
	{
		private final int DIM = 3;
		
		private double weight_dist[][] = {
				{	0.3,	0.83,	0.3,	0.1511,	0.12,	2.1,	0.03	},
				{	0.1,	0.125,	0.2,	0.25,	0.2,	0.325,	0.015	},
				{	0.15,	0.27,	0.5,	3,		0.5,	0.225,	0.3		},
				{	0.5111,	0.25,	2,		0,		1.3,	0.45,	0.6111	},
				{	0.011,	0.24,	0,		1.5003,	0.7,	0.3,	0.3		},
				{	0.18,	0.125,	0.2,	0.35,	0.4,	0,		0.03	},
				{	0,		0,		0.1,	0.3111,	0.2,	0.02,	1.1	}};
				
		@Override
		int dimConflictRegion_internal ()
		{
			return DIM;
		}
		
		@Override
		double positionWeight_internal (int base, int probe_len)
		{
			// conflicts in the beginning of the probe are weighted higher! 
			return probe_len - base;
		}
		
		@Override
		double distanceWeight_internal (int r1, int c1, int r2, int c2)
		{
			return weight_dist[DIM + r2 - r1][DIM + c2 - c1];
		}
	}
}
