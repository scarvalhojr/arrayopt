#!/usr/bin/perl
#
# makemasks.pl
# 
# $Revision$
# 
# $Date$
# 
# Copyright 2005 Sergio Anibal de Carvalho Junior
# 
# This file is part of ArrayOpt.
# 
# --- License ----------------------------------------------------------------
# ArrayOpt is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
# 
# ÂrrayOpt is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
# 
# You should have received a copy of the GNU General Public License along with
# ArrayOpt; if not, write to the Free Software Foundation, Inc., 59 Temple
# Place, Suite 330, Boston, MA 02111-1307, USA.
# ----------------------------------------------------------------------------
# 
# This is the result of a PhD work developed at the Universitaet Bielefeld
# under the supervision of Dr. Sven Rahmann. Proper attribution of the author
# as the source of the software is appreciated.
# 
# Sergio Anibal de Carvalho Jr.  http://www.cebitec.uni-bielefeld.de/~scarvalh
# AG Genominformatik             http://gi.cebitec.uni-bielefeld.de
# Universitaet Bielefeld         http://www.uni-bielefeld.de
#
#
# ==============================================================================
# Creates a visual representation of the set of masks used to create a chip.
#
# Input:  <chip code>_<layout>_<embed>        (layout with probe embedding)
# Output: <chip code>_<layout>_<embed>_masknn (set of mask maps)
#
# The EMBED file contains a list of all empty and non-empty spots in the chip
# with their corresponding embedding in the deposition sequence.
#
# The output is a set of maps depicting each mask of the synthesis steps used to
# produce the chip with the given layout. These files can later be used to
# generate bitmap (BMP) images. At each step, unmasked spots are colored with
# red whereas masked spots are painted black.
#
# Note that only "normal" and control probes have full sequence information,
# while others spots are usually listed as empty. Empty spots can be
# (supposedly) either masked or unmasked (to better integrate with its
# surroundings). The produced mask files will have a distinct color for empty
# spots (as well as those sequence is unknown) so that they can be distinguished.
#
# ------------------------------------------------------------------------------
# To do: - accept range parameters in the command line
# ==============================================================================

use strict;
use warnings;

# ------------------------------------------------------------------------------
# global variables
# ------------------------------------------------------------------------------

# embedding information
my @embed;

# length of the deposition sequence
my $len = 0;

# maximum coordinates
my $xmax = 0;
my $ymax = 0;

# keep track of non-empty spots 
my @spot;

# temp variable
my $tmp;

# ------------------------------------------------------------------------------
# check arguments
# ------------------------------------------------------------------------------

# prefix for input and output files
my $chip_code;
my $layout_code;
my $embed_code;

# $#ARGV contains the index of the last element of @ARGV
if ($#ARGV == 2)
{
	# file name prefixes
	$chip_code = $ARGV[0];
	$layout_code = $ARGV[1];
	$embed_code = $ARGV[2];
}
else
{
	print "Usage: $0 <chip code> <layout> <embedding>\n";
	print "       input file:  <chip code>_<layout>_<embedding>        (layout with embedding)\n";
	print "       output file: <chip code>_<layout>_<embedding>_masknn (set of mask maps)\n";
	exit 1;
}

# file names
my $input_file = $chip_code . "_" . $layout_code . "_" . $embed_code;
my $output_file = $chip_code . "_" . $layout_code . "_" . $embed_code . "_mask";

# ------------------------------------------------------------------------------
# check if output file already exists
# ------------------------------------------------------------------------------

if (-f $output_file . "00")
{
	print "Output file '${output_file}00' already exists. Overwrite all? (Y/N): ";
	chomp($tmp = <STDIN>);
	
	($tmp !~ m/^[Yy]/) and exit 0;
}

# ------------------------------------------------------------------------------
# read input file
# ------------------------------------------------------------------------------

my ($x, $y, $group, $fixed, $pm, $sequence, $embedding);
my ($ln, $count_empty, $count_embed) = (0, 0, 0);

open (INFILE, $input_file) or
	die "Could not open input file '$input_file' for reading: $!.\n";

print "Reading '$input_file'\n";

for (<INFILE>)
{
	$ln++;
	
	($x, $y, $group, $fixed, $pm, $sequence, $embedding) = split /\t|\n|\r\n/;

	# keep track of maximum values
	$xmax = ($x > $xmax) ? $x : $xmax;
	$ymax = ($y > $ymax) ? $y : $ymax;

	# check if spot is already in use
	no warnings 'uninitialized';
	if ($spot[$x][$y] > 0)
	{
		print "Duplicate spot found at position ($x, $y).\n";
		exit 1;
	}
	use warnings 'uninitialized';
	
	if ($embedding eq "-")
	{
		# empty spot: nothing to do
		$count_empty++;
	}
	elsif ($embedding =~ /^[ACGT ]+$/)
	{
		# save embedding
		$embed[$x][$y] = $embedding;
		
		# all embeddings must have the same length
		if ($len == 0)
		{
			$len = length $embedding;
		}
		elsif (length $embedding != $len)
		{
			print "Embeddings of different length found in line $ln.\n";
			exit 1;
		}
		
		# mark position as used
		$spot[$x][$y] = 1;
		$count_embed++;
	}
	else
	{
		print "Invalid embedding in line $ln: [$embedding].\n";
		exit 1;
	}
}

close INFILE;

print "$count_empty empty sequences.\n";
print "$count_embed sequences embedded.\n";
print "Chip has ", $ymax + 1, " rows and ", $xmax + 1, " columns\n";

# ------------------------------------------------------------------------------
# write mask maps
# ------------------------------------------------------------------------------

my ($id, $base, $out_char, $ca, $ra);

my $row_amp = 0;
my $col_amp = 0;

no warnings 'uninitialized';

for (my $m = 0; $m < $len; $m++)
{
	# open output file
	$id = sprintf "%02d", $m;
	open (OUTFILE, '>', ${output_file} . ${id} ) or
		die "Could not open output file '$output_file' for writing: $!.\n";

	print "Writing mask $id to '${output_file}${id}'\n";
	
	$base = "";

	for ($y = 0; $y <= $ymax; $y++)
	{
		for ($ra = 0; $ra <= $row_amp; $ra++)
		{
			for ($x = 0; $x <= $xmax; $x++)
			{
				if ($spot[$x][$y] > 0)
				{
					if (substr($embed[$x][$y], $m, 1) eq " ")
					{
						# spot is masked at this step
						$out_char = "F";
					}
					else
					{
						# spot is unmasked at this step
						$out_char = "3";

						if ($base eq "")
						{
							# base that is synthesized at this step
							$base = substr($embed[$x][$y], $m, 1);
						}
						elsif (substr($embed[$x][$y], $m, 1) ne $base)
						{
							# Invalid embedding!
							die "Invalid embedding for sequence at position ($x, $y), step $m\n";
						}
					}
				}
				else
				{
					# empty spot			
					$out_char = "7";
				}
				
				for ($ca = 0; $ca <= $col_amp; $ca++)
				{
					print OUTFILE $out_char;
				}
			}
			print OUTFILE "\n";
		}
	}
	
	close OUTFILE;
}
