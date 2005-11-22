#!/usr/bin/perl
#
# fixedmap.pl
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
# Make a visual map of a chip, highlighting fixed spots.
#
# Input:  <chip code>_<layout>	        (chip's layout specification)
# Output: <chip code>_<layout>_fixedmap (map highlighting fixed spots)
#
# This program creates a visual map of a chip just like "chipmap.pl". However,
# it tries to highlight the fixed spots, i.e. those which must remain empty or
# whose sequence cannot be moved due to physical constraints of the layout such
# as the arrangement of the control and QC probes.
#
# In order to improve the visibility of such squares, PM and MM probes receive
# the same color (red) when they are not fixed. If they are fixed (which would
# not normally occur, they receive a light blut (aqua) color. Control probes
# (AFFX groups) are colored with a dark red when they are fixed, or with a
# normal red wheny they are not. QC probes receive a dark red color if fixed. If
# they are not fixed (which is also not expected), they get a light blue color.
# Empty fixed spots are yellow colored while those that are not fixed get a
# normal red.
#
# Color table
# --------------------------------------
# 
# 3 (red)         - normal probe (PM or MM), not fixed
# 5 (aqua)        - normal probe (PM or MM), fixed
# 3 (red)         - control probe (AFFX), not fixed
# B (dark red)    - control probe (AFFX), fixed
# 5 (aqua)        - QC probe, not fixed
# B (dark red)    - QC probe, fixed
# 1 (yellow)      - empty spots, fixed
# 3 (red)         - empty spots, not fixed
#
# ------------------------------------------------------------------------------
# To do: - run Pixels2BMP automatically
# ==============================================================================

use strict;
use warnings;

# ------------------------------------------------------------------------------
# global variables
# ------------------------------------------------------------------------------

# map array
my @pixel;

# temp variable
my $tmp;

# ------------------------------------------------------------------------------
# check arguments
# ------------------------------------------------------------------------------

# chip code: prefix for input and output files
my $chip_code;

# layout code: prefix for input and output files
my $layout_code;

# $#ARGV contains the index of the last element of @ARGV
if ($#ARGV == 1)
{
	# get chip code
	$chip_code = $ARGV[0];

	# get layout code
	$layout_code = $ARGV[1];
}
else
{
	print "Usage: $0 <chip code> <layout>\n";
	print "       input file:  <chip code>_<layout>          (chip's layout specification)\n";
	print "       output file: <chip code>_<layout>_fixedmap (map highlighting fixed spots)\n";
	exit 1;
}

# file names
my $input_file = $chip_code . "_" . $layout_code;
my $output_file = $chip_code . "_" . $layout_code . "_fixedmap";

# ------------------------------------------------------------------------------
# check if output file already exists
# ------------------------------------------------------------------------------

if (-f $output_file)
{
	print "Output file '$output_file' already exists. Overwrite? (Y/N): ";
	chomp($tmp = <STDIN>);
	
	($tmp !~ m/^[Yy]/) and exit 0;
}

# ------------------------------------------------------------------------------
# read input file
# ------------------------------------------------------------------------------

my ($xpos, $ypos, $group, $fixed, $pm, $seq, $color);
my ($xmax, $ymax, $line) = (0, 0, 0);

# open input file for reading
open (INFILE, $input_file) or
	die "Could not open CHIP file '$input_file' for reading: $!.\n";

print "Reading '$input_file'\n";

while (<INFILE>)
{
	$line++;
	
	($xpos, $ypos, $group, $fixed, $pm, $seq) = split /\t|\n|\r\n/;
	
	if (!($xpos =~ /[0-9]+/ && $ypos =~ /[0-9]+/) || $xpos < 0 || $ypos < 0)
	{
		print "Invalid position ($xpos, $ypos) at line $line.\n";
		exit 1;
	}
	
	# keep track of maximum values
	$xmax = ($xpos > $xmax) ? $xpos : $xmax;
	$ymax = ($ypos > $ymax) ? $ypos : $ymax;

	# check if spot is already in use
	no warnings 'uninitialized';
	if ($pixel[$xpos][$ypos] > 0)
	{
		print "Duplicate spot found at position ($xpos, $ypos).\n";
		exit 1;
	}
	use warnings 'uninitialized';
	
	if ($fixed eq "Y")
	{
		# spot is fixed
		
		if (($group =~ /^QC-type/) or ($group =~ /^AFFX/))
		{
			# QC and control probes: dark red
			$pixel[$xpos][$ypos] = 0xB;
		}
		elsif ($seq eq "-")
		{
			# empty fixed spots: yellow
			$pixel[$xpos][$ypos] = 0x1;
		}
		else
		{
			# strange... normal probes are not
			# expected to be on fixed spots: aqua
			$pixel[$xpos][$ypos] = 0x5;
		}
	}
	else
	{
		# spot is not fixed
		
		if ($group =~ /^QC-type/)
		{
			# strange... QC probes are expected
			# to be on fixed spots: aqua
			$pixel[$xpos][$ypos] = 0x5;
		}
		else
		{
			# normal probes or empty spots
			# that are not fixed: red
			$pixel[$xpos][$ypos] = 0x3;
		}
	}
}

close INFILE;

print $line, " spots described\n";

# ------------------------------------------------------------------------------
# write output file
# ------------------------------------------------------------------------------

open (OUTFILE, ">", $output_file) or
	die "Could not open '$output_file' for writing: $!.\n";

print "Writing output to '$output_file'\n";

no warnings 'uninitialized';

for (my $y = 0; $y <= $ymax; $y++)
{
	for (my $x = 0; $x <= $xmax; $x++)
	{
		if ($pixel[$x][$y] > 0)
		{
			printf OUTFILE "%X", $pixel[$x][$y];
		}
		else
		{
			# unknown spots are considered
			# empty and not fixed: red
			print OUTFILE "3";
		}
	}
	print OUTFILE "\n";
}

close OUTFILE;

print "Chip has ", $ymax + 1, " rows and ", $xmax + 1, " columns\n";
