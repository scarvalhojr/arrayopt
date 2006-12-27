/*
 * BMPFile.java
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

import java.io.*;

/**
 * This is a helper class for writing Bitmap files (BMP). The class
 * automatically computes several parameters needed for the BMP headers and
 * contains methods for obtaining the pixel colors in the correct form, and
 * ensuring the proper padding of rows. 
 * 
 * <P>The user of this class should first create an instance of this class
 * giving the dimension of the picture (number of rows and columns of pixels
 * in the BMP image) and a OutputStream where the contents will be written to.
 * The next step consists of writing the file headers with the
 * {@link #writeHeader} method.</P>
 * 
 * <P>Then, the user can write individual pixels directly with the
 * {@link java.io.OutputStream#write(byte[])}, passing an array of bytes
 * representing the color of each pixel. This array of bytes can be easily
 * obtained with the {@link #getRGBColor(int)}. Note that the image is stored
 * from bottom to top, meaning that the first scan line is the last scan line
 * in the image. After each row is completed, the user must call the
 * {@link #finishRow()} method to ensure that each scan line is padded to an
 * even 4-byte boundary.</P>
 * 
 * @author Sergio A. de Carvalho Jr.
 */
public class BMPFile
{
	private final static int BMP_FILEHEADER_SIZE = 14;
	private final static int BMP_INFOHEADER_SIZE = 40;
	
	// bitmap file header
	private final static byte bfh_type [] = {'B', 'M'};
	private int bfh_size;
	private final static int bfh_reserved1 = 0;
	private final static int bfh_reserved2 = 0;
	private final static int bfh_offset = BMP_FILEHEADER_SIZE +
											BMP_INFOHEADER_SIZE;

	// bitmap info header
	private final static int bih_size = BMP_INFOHEADER_SIZE;
	private int bih_width;
	private int bih_height;
	private final static int bih_planes = 1;
	private final static int bih_bitcount = 24;
	private final static int bih_compression = 0;
	private int bih_imagesize;
	private final static int bih_xpixels_per_meter = 0x0;
	private final static int bih_ypixels_per_meter = 0x0;
	private final static int bih_usedcolors = 0;
	private final static int bih_importantcolors = 0;

	private int line_pad;
	
	private OutputStream out;
	
	/**
	 * Creates a new Bitmap (BMP) file with the specified number of rows and
	 * columns on the given OutputStream.
	 * 
	 * @param num_rows number of rows of pixels
	 * @param num_cols number of columns of pixels
	 * @param out output stream where the BMP will be written to
	 */
	public BMPFile (int num_rows, int num_cols, OutputStream out)
	{
		if (num_rows <= 0 || num_cols <= 0)
			throw new IllegalArgumentException
				("Invalid number of rows or columns");
		
		this.bih_width = num_cols;
		this.bih_height = num_rows;
		
		this.line_pad = (4 - ((bih_width * 3) % 4)) % 4;

		this.bih_imagesize = bih_height * (3 * bih_width + line_pad);
		
		this.bfh_size = bih_imagesize + bfh_offset;
		
		this.out = out;
	}
	
	/**
	 * Returns the RGB representation of a color in an array of bytes as it
	 * must be written in the BMP file. Colors should are specificed by an
	 * integer value where each color component (R = red, G = green, B = blue)
	 * is given using eight bits. Examples: <CODE>0xFF0000</CODE> (red),
	 * <CODE>0x0000FF</CODE> (blue),  <CODE>0x000080</CODE> (dark blue), 
	 * <CODE>0xFFFF00</CODE> (yellow), <CODE>0x000000</CODE> (black),
	 * <CODE>0xFFFFFF</CODE> (white).    
	 *  
	 * @param color color in RGB notation
	 * @return byte array to be written in a BMP file
	 */
	public static byte[] getRGBColor (int color)
	{
		byte rgb[] = new byte[3];
		
		rgb [0] = (byte) (color & 0xFF);
		rgb [1] = (byte) ((color >> 8) & 0xFF);
		rgb [2] = (byte) ((color >>  16) & 0xFF);
		
		return rgb;
	}
	
	/**
	 * Writes the Bitmap headers.
	 * 
	 * @throws IOException
	 */
	public void writeHeader () throws IOException
	{
		// write the bitmap file header
		out.write (bfh_type);
		out.write (intToDWord (bfh_size));
		out.write (intToWord (bfh_reserved1));
		out.write (intToWord (bfh_reserved2));
		out.write (intToDWord (bfh_offset));

		// write the bitmap information header
		out.write (intToDWord (bih_size));
		out.write (intToDWord (bih_width));
		out.write (intToDWord (bih_height));
		out.write (intToWord (bih_planes));
		out.write (intToWord (bih_bitcount));
		out.write (intToDWord (bih_compression));
		out.write (intToDWord (bih_imagesize));
		out.write (intToDWord (bih_xpixels_per_meter));
		out.write (intToDWord (bih_ypixels_per_meter));
		out.write (intToDWord (bih_usedcolors));
		out.write (intToDWord (bih_importantcolors));
	}
	
	/**
	 * Fills a row with zeros so that the number of bytes per row is multiple
	 * of four.
	 *  
	 * @throws IOException
	 */
	public void finishRow () throws IOException
	{
		for (int i = 0; i < line_pad; i++)
			out.write (0x00);
	}
	
	/**
	 * Converts an int to a word, returning the value in a 2-byte array.
	 */
	private byte [] intToWord (int value)
	{
		byte word [] = new byte [2];
		word [0] = (byte) (value & 0x00FF);
		word [1] = (byte) ((value >>  8) & 0x00FF);
		return word;
	}

	/**
	 * Converts an int to a double word, returning the value in a 4-byte array.
	 */
	private byte [] intToDWord (int value)
	{
		byte dword [] = new byte [4];
		dword [0] = (byte) (value & 0x00FF);
		dword [1] = (byte) ((value >>  8) & 0x000000FF);
		dword [2] = (byte) ((value >>  16) & 0x000000FF);
		dword [3] = (byte) ((value >>  24) & 0x000000FF);
		return (dword);
	}
}
