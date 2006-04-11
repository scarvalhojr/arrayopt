#!/usr/bin/perl
#
# analyze-bl.pl
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
# Analyze the border length of a chip layout and embedding.
#
# Input:  <chip code>_<layout>_<embed>          (layout with probe embedding)
# Output: <chip code>_<layout>_<embed>_border   (border length of mask step)
#
# The layout file contains a list of all empty and non-empty spots in the chip
# with their corresponding embedding in the deposition sequence. Note that only
# "normal" and control probes have full sequence information, while others are
# listed as empty spots.
#
# This program output the border length of each synthesis step.
#
# ------------------------------------------------------------------------------
# To do: - documentation
# ==============================================================================

use strict;
use warnings;

# ------------------------------------------------------------------------------
# global variables
# ------------------------------------------------------------------------------

# embedding information
my @embed;

# keep track of non-empty spots 
my @spot;

# length of the deposition sequence
my $ds_len = 0;

# length of probes
my $probe_len = 0;

# maximum coordinates
my $max_x = 0;
my $max_y = 0;

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
	print "       input file:   <chip code>_<layout>_<embed> (layout with probe embedding)\n";
	print "       output file: <chip code>_<embed>_border    (border length of mask steps)\n";
	exit 1;
}

# file names
my $input_file = $chip_code . "_" . $layout_code . "_" . $embed_code;
my $border_file = $chip_code . "_" . $layout_code . "_". $embed_code . "_border";

# ------------------------------------------------------------------------------
# check if output file already exists
# ------------------------------------------------------------------------------

if (-f $border_file)
{
	print "Output file '${border_file}' already exists. Overwrite? (Y/N): ";
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
	$max_x = ($x > $max_x) ? $x : $max_x;
	$max_y = ($y > $max_y) ? $y : $max_y;

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

		# mark position as empty
		$spot[$x][$y] = 0;
	}
	elsif ($embedding =~ /^[ACGT ]+$/)
	{
		# save embedding
		$embed[$x][$y] = $embedding;

		# all sequences must have the same length
		if ($probe_len == 0)
		{
			$probe_len = length $sequence;
		}
		elsif (length $sequence != $probe_len)
		{
			print "Sequence of different length found in line $ln.\n";
			exit 1;
		}
		
		# all embeddings must have the same length
		if ($ds_len == 0)
		{
			$ds_len = length $embedding;
		}
		elsif (length $embedding != $ds_len)
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
		print "Invalid embedding in line $ln: $sequence.\n";
		exit 1;
	}
}

close INFILE;

print "$count_empty empty sequences.\n";
print "$count_embed sequences embedded.\n";
print "Chip has ", $max_y + 1, " rows and ", $max_x + 1, " columns\n";
print "Length of deposition sequence: $ds_len\n";

# ------------------------------------------------------------------------------
# mark empty spots
# ------------------------------------------------------------------------------

no warnings 'uninitialized';

for ($y = 0; $y <= $max_y; $y++)
{
	for ($x = 0; $x <= $max_x; $x++)
	{
		if (not ($spot[$x][$y] > 0))
		{
			$spot[$x][$y] = 0;
		}
	}
}

use warnings 'uninitialized';

# ------------------------------------------------------------------------------
# compute border length
# ------------------------------------------------------------------------------

my ($total_border, $log10) = (0, log(10));
my ($border, $border_norm, $cost, $rx, $ry, $pos_mult, $weight);

# open BORDER output file
open (OUTFILE, '>', ${border_file}) or
	die "Could not open file '$border_file' for writing: $!.\n";

print "Writing border length values to '$border_file'\n";

for (my $m = 0; $m < $ds_len; $m++)
{
	$border = 0;

	for ($y = 0; $y <= $max_y; $y++)
	{
		for ($x = 0; $x <= $max_x; $x++)
		{
			# if spot not empty
			if ($spot[$x][$y] > 0)
			{
				# -> count direct border conflicts:
				
				# top border (horizontal)
				if (($y > 0) and ($spot[$x][$y - 1] > 0) and not
					(substr($embed[$x][$y], $m, 1) eq substr($embed[$x][$y - 1], $m, 1)))
				{
					$border++;
				}

				# left border (vertical)
				if (($x > 0) and ($spot[$x - 1][$y] > 0) and not
					(substr($embed[$x][$y], $m, 1) eq substr($embed[$x - 1][$y], $m, 1)))
				{
					$border++;
				}
			}
		}
	}
	
	# border length normalized by the number of non-empty spots
	$border_norm = $border / $count_embed;
	$total_border += $border_norm;
	
	print OUTFILE "$m\t$border_norm\n";
	print "Step $m => border length: $border (normalized: $border_norm)\n";
}

close OUTFILE;

print "Total border length: $total_border\n";
