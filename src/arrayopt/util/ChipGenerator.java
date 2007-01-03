/*
 * ChipGenerator.java
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

import arrayopt.layout.*;

/**
 * 
 * TODO move the {@link Chip#createRandomLayout} method to this class
 * TODO set the {@link Chip#input_done} flag once a layout is created
 * TODO initialize the {@link Chip#fixed_probe} array once a layout iscreated
 * 
 * @author Sergio A. de Carvalho Jr.
 */
public class ChipGenerator
{
	/**
	 * Constant to indicate the generation of a layout containing all possible
	 * embeddings that result in probes with the appropriate legth.
	 * 
	 * <P>The current implementation generates all possible binary vectors,
	 * checking whether they lead to a probe of the desired length. In case it
	 * does, the new embedding is placed in the next available spot.</P>
	 * 
	 * <P>The number of probes in the chip must be equal to the number of
	 * possible embeddings, which is equal to <CODE>e choose p</CODE>, where
	 * <CODE>e</CODE> is the embedding length ({@link Chip#embed_len}) and
	 * <CODE>p</CODE> is the probe length ({@link Chip#probe_len}). If the
	 * number of probes is less than or greater than that, an
	 * {@link IllegalArgumentException} is thrown.</P> 
	 */
	public static final int ALL_ASYNC_EMBEDDINGS = 0;
	
	public static final int ALL_SYNC_EMBEDDINGS = 1;
	
	private Chip chip;
	
	private int dep_seq_cycle;
	
	private int embed_len;
	
	private int probe_len;
	
	private int num_probes;
	
	private int num_rows;
	
	private int num_cols;
	
	private int id;
	
	private int row;
	
	private int col;
	
	public ChipGenerator (Chip chip)
	{
		this.chip = chip;
		this.dep_seq_cycle = chip.dep_seq_cycle;
		this.embed_len = chip.getEmbeddingLength();
		this.probe_len = chip.getProbeLength();
		this.num_probes = chip.getNumberOfProbes();
		this.num_rows = chip.getNumberOfRows();
		this.num_cols = chip.getNumberOfColumns();
	}
	
	public void generateLayout (int type)
	{
		this.id = 0;
		this.row = 0;
		this.col = 0;
		
		switch (type)
		{
			case ALL_ASYNC_EMBEDDINGS:
				if (chip instanceof SimpleChip)
					allAsyncEmbeddings();
				else
					throw new IllegalArgumentException
						("Generation not supported for this type of chip.");
				break;
			
			case ALL_SYNC_EMBEDDINGS:
				if (chip instanceof SimpleChip)
					allSyncEmbeddings();
				else
					throw new IllegalArgumentException
						("Generation not supported for this type of chip.");
				break;
				
			default:
				throw new IllegalArgumentException
					("Unknown generation .");
		}
	}
	
	private void allSyncEmbeddings ()
	{
		byte probe[];
		
		if (dep_seq_cycle <= 0)
			throw new IllegalArgumentException
				("Deposition sequence must be cyclical.");
		
		if ((chip.dep_seq.length % dep_seq_cycle) != 0)
			throw new IllegalArgumentException
				("Deposition sequence must not contain incomplete cycles.");
		
		probe = new byte[probe_len];
		
		allSyncEmbeddings (probe, 0);

		if (id < num_probes)
			throw new IllegalArgumentException
				("More probes than possible embeddings");
		
		fillEmptySpots ();
	}
	
	private void allSyncEmbeddings (byte probe[], int pos)
	{
		if (pos >= probe_len)
		{
			addProbe (probe);
			return;
		}
		
		for (byte b = 0; b < dep_seq_cycle; b++)
		{
			probe[pos] = b;
			allSyncEmbeddings (probe, pos + 1);
		}
	}
	
	private void allAsyncEmbeddings ()
	{
		byte emb[] = new byte[embed_len];
				
		allAsyncEmbeddings (emb, 0, probe_len);
		
		if (id < num_probes)
			throw new IllegalArgumentException
				("More probes than possible embeddings");
		
		fillEmptySpots ();
	}
	
	private void allAsyncEmbeddings (byte emb[], int pos, int bases)
	{
		if (bases > embed_len - pos) return;
		
		if (pos < embed_len - 1)
		{
			emb[pos] = 0;
			allAsyncEmbeddings (emb, pos + 1, bases);
			if (bases > 0)
			{
				emb[pos] = 1;
				allAsyncEmbeddings (emb, pos + 1, bases - 1);
			}
		}
		else
		{
			if (bases == 0)
			{
				emb[pos] = 0;
				addEmbedding (emb);
			}
			else if (bases == 1)
			{
				emb[pos] = 1;
				addEmbedding (emb);
			}
		}
	}
	
	private void addProbe (byte probe[])
	{
		int w, pos, bitmask = 0;
		
		if (id >= num_probes)
			throw new IllegalArgumentException
				("Insufficient number of probes.");
		
		// turn all bits off
		for (w = 0; w < chip.embed[id].length; w++)
			chip.embed[id][w] = 0;

		for (w = -1, pos = 0; pos < embed_len; pos++)
		{
			if (pos % Integer.SIZE == 0)
			{
				// turn on very first bit of mask only
				bitmask = 0x01 << (Integer.SIZE - 1);
				w++;
			}
			else
				bitmask >>>= 1;
			
			if (probe[pos / dep_seq_cycle] == (pos % dep_seq_cycle))
				chip.embed[id][w] |= bitmask;
		}
		
		chip.spot[row][col] = id;
		if (++row >= num_rows)
		{
			row = 0;
			col++;
		}
		
		id++;
	}
	
	private void addEmbedding (byte emb[])
	{
		int w, pos, bitmask = 0;
		
		if (id >= num_probes)
			throw new IllegalArgumentException
				("Insufficient number of probes.");
		
		// turn all bits off
		for (w = 0; w < chip.embed[id].length; w++)
			chip.embed[id][w] = 0;

		for (w = -1, pos = 0; pos < embed_len; pos++)
		{
			if (pos % Integer.SIZE == 0)
			{
				// turn on very first bit of mask only
				bitmask = 0x01 << (Integer.SIZE - 1);
				w++;
			}
			else
				bitmask >>>= 1;

			if (emb[pos] == 1)
				chip.embed[id][w] |= bitmask;
		}
		
		chip.spot[row][col] = id;
		if (++row >= num_rows)
		{
			row = 0;
			col++;
		}
		
		id++;
	}
	
	private void fillEmptySpots ()
	{
		while (col < num_cols)
		{
			chip.spot[row][col] = Chip.EMPTY_SPOT;
			
			if (++row >= num_rows)
			{
				row = 0;
				col++;
			}
		}
	}
}
