#!/usr/bin/perl
#
# analyze.pl
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
# Analyze a chip layout and embedding.
#
# Input:  <chip code>_<layout>_<embed>          (layout with probe embedding)
# Output: <chip code>_<layout>_<embed>_border   (border length of mask step)
#         <chip code>_<layout>_<embed>_conflict (spots' conflict indexes)
#         <chip code>_<layout>_<embed>_confdist (distribution of conflict indexes)
#
# The layout file contains a list of all empty and non-empty spots in the chip
# with their corresponding embedding in the deposition sequence. Note that only
# "normal" and control probes have full sequence information, while others are
# listed as empty spots.
#
# This program output the border length of each synthesis step in the BORDER
# file. It will also compute the conflict index of every non-empty spot of the
# chip, which is written in the CONFLICT file. It will also agreggate these
# indexes to their rounding integers and output their distribution in the
# CONFDIST file.
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

# position pointers
my @pos;

# keep track of non-empty spots 
my @spot;

# spot's conflict index
my @conflict;

# distance-based conflict weight function
my @weight_dist = (	[	0,		0,		0.1,	0.1111,	0.1,	0,		0		],
					[	0,		0.125,	0.2,	0.25,	0.2,	0.125,	0		],
					[	0.1,	0.2,	0.5,	1,		0.5,	0.2,	0.1		],
					[	0.1111,	0.25,	1,		0,		1,		0.25,	0.1111	],
					[	0.1,	0.2,	0.5,	1,		0.5,	0.2,	0.1		],
					[	0,		0.125,	0.2,	0.25,	0.2,	0.125,	0		],
					[	0,		0,		0.1,	0.1111,	0.1,	0,		0		]);

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
	print "       output files: <chip code>_<embed>_border   (border length of mask steps)\n";
	print "       output files: <chip code>_<embed>_conflict (spots' conflict indexes)\n";
	print "       output files: <chip code>_<embed>_confdist (distribution of conflict indexes)\n";
	exit 1;
}

# file names
my $input_file = $chip_code . "_" . $layout_code . "_" . $embed_code;
my $border_file = $chip_code . "_" . $layout_code . "_". $embed_code . "_border";
my $conflict_file = $chip_code . "_" . $layout_code . "_". $embed_code . "_conflict";
my $confdist_file = $chip_code . "_" . $layout_code . "_". $embed_code . "_confdist";

# ------------------------------------------------------------------------------
# check if output files already exist
# ------------------------------------------------------------------------------

if (-f $border_file)
{
	print "Output file '${border_file}' already exists. Overwrite? (Y/N): ";
	chomp($tmp = <STDIN>);
	
	($tmp !~ m/^[Yy]/) and exit 0;
}

if (-f $conflict_file)
{
	print "Output file '${conflict_file}' already exists. Overwrite? (Y/N): ";
	chomp($tmp = <STDIN>);
	
	($tmp !~ m/^[Yy]/) and exit 0;
}

if (-f $confdist_file)
{
	print "Output file '${confdist_file}' already exists. Overwrite? (Y/N): ";
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
# initialize position pointers and conflict indexes
# ------------------------------------------------------------------------------

no warnings 'uninitialized';

for ($y = 0; $y <= $max_y; $y++)
{
	for ($x = 0; $x <= $max_x; $x++)
	{
		$pos[$x][$y] = 0;
		$conflict[$x][$y] = 0;
		
		if (not ($spot[$x][$y] > 0))
		{
			$spot[$x][$y] = 0;
		}
	}
}

use warnings 'uninitialized';

# ------------------------------------------------------------------------------
# compute border length and conflict indexes
# ------------------------------------------------------------------------------

my ($total_border) = (0);
my ($border, $border_norm, $cost, $rx, $ry, $pos_mult);

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
				
				# -> compute spot's conflict index

				if (not (substr($embed[$x][$y], $m, 1) eq " "))
				{
					# spot is in an unmasked step
					# (a nucleotide is synthesized)
					
					# advance position pointer
					$pos[$x][$y]++;
					
					# no conflict in this step (light directed to
					# neighboring spots cannot cause any damage here)
					next;
				}
				
				# compute position multiplier: conflicts
				# are more harmful in the middle of a
				# (a conflict would harm the next nucleotide
				# to be synthesized, if there is any)
				# pos_mult = sqrt ( min (pos + 1, probe_len - pos + 1) )
				$pos_mult = ($pos[$x][$y] + 1 <= $probe_len - $pos[$x][$y] + 1) ?
							 $pos[$x][$y] + 1 :  $probe_len - $pos[$x][$y] + 1;
				
				$pos_mult = sqrt($pos_mult);

				$cost = 0;
				$rx = ($x >= 3) ? $x - 3 : 0;
				for (; $rx <= $x + 3 && $rx <= $max_x; $rx++)
				{
					$ry = ($y >= 3) ? $y - 3 : 0;
					for (; $ry <= $y + 3 && $ry <= $max_y; $ry++)
					{
						# skip if neighbor is not close enough
						# for influencing the conflict index
						# or if (rx, ry) equals (x, y)
						($weight_dist[$rx - $x + 3][$ry - $y + 3] == 0) and next;
						
						# skip if neighbor is empty
						($spot[$rx][$ry] == 0) and next;
						
						# conflict only when neighbor is unmasked:
						# (light directed to neighbor can activate this spot)
						(substr($embed[$rx][$ry], $m, 1) eq " ") and next;
						
						$cost += $pos_mult * $weight_dist[$rx - $x + 3][$ry - $y + 3];
					}
				}
				
				$conflict[$x][$y] += $cost;
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

# ------------------------------------------------------------------------------
# print conflict indexes
# ------------------------------------------------------------------------------

# open CONFLICT output file
open (OUTFILE, '>', ${conflict_file}) or
	die "Could not open file '$conflict_file' for writing: $!.\n";

print "Writing conflict indexes to '$conflict_file'\n";

for ($x = 0; $x <= $max_x; $x++)
{
	for ($y = 0; $y <= $max_y; $y++)
	{
		print OUTFILE "$x\t$y\t$conflict[$x][$y]\n";
	}
}

close OUTFILE;

# ------------------------------------------------------------------------------
# compute distribution of conflict indexes
# ------------------------------------------------------------------------------

my ($avg_conf, $min_conf_x, $min_conf_y, $max_conf_x, $max_conf_y) = (0, -1, -1, -1, -1);
my @conf_bucket;

for ($x = 0; $x <= $max_x; $x++)
{
	for ($y = 0; $y <= $max_y; $y++)
	{
		# skip empty spots
		($spot[$x][$y] == 0) and next;
		
		# average conflict index
		$avg_conf += $conflict[$x][$y] / $count_embed;
		
		if ($min_conf_x == -1 or $conflict[$x][$y] < $conflict[$min_conf_x][$min_conf_y])
		{
			$min_conf_x = $x;
			$min_conf_y = $y;
		}
		
		if ($max_conf_x == -1 or $conflict[$x][$y] > $conflict[$max_conf_x][$max_conf_y])
		{
			$max_conf_x = $x;
			$max_conf_y = $y;
		}
		
		$conf_bucket[int ($conflict[$x][$y] + .5)]++;
	}
}
		
print "Min conflict: $conflict[$min_conf_x][$min_conf_y] ($min_conf_x, $min_conf_y)\n";
print "Avg conflict: $avg_conf\n";
print "Max conflict: $conflict[$max_conf_x][$max_conf_y] ($max_conf_x, $max_conf_y)\n";

# ------------------------------------------------------------------------------
# print distribution of conflict indexes
# ------------------------------------------------------------------------------

# open CONFDIST output file
open (OUTFILE, '>', ${confdist_file}) or
	die "Could not open file '$confdist_file' for writing: $!.\n";

print "Writing distribtion of conflict indexes to '$confdist_file'\n";

no warnings 'uninitialized';

for (my $c = 0; $c <= int($conflict[$max_conf_x][$max_conf_y]); $c++)
{
	if ($conf_bucket[$c] > 0)
	{
		print OUTFILE "$c\t$conf_bucket[$c]\n";
	}
}

close OUTFILE;
