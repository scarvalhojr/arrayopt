#!/usr/bin/perl
#
# chipmap.pl
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
# Make a visual map of a chip layout.
#
# Input:  <chip code>_<layout>         (chip's layout specification)
# Output: <chip code>_<layout>_chipmap (visual map of the layout)
#
# This program creates a visual map of all spots in the chip that can later be
# used to generate a bitmap (BMP) image of the chip layout. Spots receive a
# color code according to their types. This code is a hexadecimal number (from
# zero to 9, A, B, C, D, E and F). Each number will be mapped to a color in the
# bitmap image.
#
# The layout file contains a list of all empty and non-empty spots in the chip.
# Non-empty are either normal probes, control probes (those whose group name
# start with AFFX) or QC (quality control) probes. Each type of probe receive a
# distinct color. Furthermore, normal probes receive one of two color codes, one
# for PM (perfect match) probes and one for MM (mismtach probes). QC probes
# receive a color code according to their specific type. Finally, empty spots
# also receive two distinct color codes, one for those that must remain empty
# (blanks) and one for those which  are empty simply because there is not
# enough probes to fill all spots.
#
# Color table
# --------------------------------------
# 3 (red)         - normal probe, PM
# 2 (orange)      - normal probe, MM
# B (dark red)    - control probe (AFFX)
# C (dark green)  - QC probe, type 1
# A (dark orange) - QC probe, type 2
# 9 (dark yellow) - QC probe, type 3
# 1 (yellow)      - QC probe, type 4 or 13
# D (dark blue)   - QC probe, type 5 or 11
# 5 (aqua)        - QC probe, type 6 or 12
# E (violet)      - QC probe, type 9
# 6 (pink)        - QC probe, type 10
# F (black)       - QC probe, type 15
# 8 (dark gray)   - QC probe, type 16
# 4 (green)       - QC probe, other type
# 7 (gray)        - fixed empty spots
# 0 (white)       - empty spots
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
	print "       input file:  <chip code>_<layout>         (chip's layout specification)\n";
	print "       output file: <chip code>_<layout>_chipmap (visual map of the layout)\n";
	exit 1;
}

# file names
my $input_file = $chip_code . "_" . $layout_code;
my $output_file = $chip_code . "_" . $layout_code . "_chipmap";

# ------------------------------------------------------------------------------
# check if output file already exists
# ------------------------------------------------------------------------------

if (-f $output_file)
{
	print "Output file '$output_file' already exists. Overwrite? (Y/N): ";
	chomp($tmp = <STDIN>);
	
	($tmp !~ m/[Yy]$/) and exit 0;
}

# ------------------------------------------------------------------------------
# read input file
# ------------------------------------------------------------------------------

my ($xpos, $ypos, $group, $fixed, $pm, $seq, $color);
my ($xmax, $ymax, $line) = (0, 0, 0);

# open input file for reading
open (INFILE, $input_file) or
	die "Could not open file '$input_file' for reading: $!.\n";

print "Reading '$input_file'\n";

while (<INFILE>)
{
	$line++;
	
	($xpos, $ypos, $group, $fixed, $pm, $seq) = split /\t|\n|\r\n/;
	
	if (!($xpos =~ /[0-9]+/ && $ypos =~ /[0-9]+/) || $xpos < 0 || $ypos < 0)
	{
		print "Invalid spot coordinate ($xpos, $ypos) at line $line.\n";
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

	# does the spot contain a QC probe?
	if ($group =~ /^QC-type/)
	{
		# color distinct QC probes
		if    ($group =~ /^QC-type1$/)  {$color = 0xC;}
		elsif ($group =~ /^QC-type2$/)  {$color = 0xA;}
		elsif ($group =~ /^QC-type3$/)  {$color = 0x9;}
		elsif ($group =~ /^QC-type4$/)  {$color = 0x1;}
		elsif ($group =~ /^QC-type5$/)  {$color = 0xD;}
		elsif ($group =~ /^QC-type6$/)  {$color = 0x5;}
		elsif ($group =~ /^QC-type9$/)	{$color = 0xE;}
		elsif ($group =~ /^QC-type10$/) {$color = 0x6;}
		elsif ($group =~ /^QC-type11$/) {$color = 0xD;}
		elsif ($group =~ /^QC-type12$/) {$color = 0x5;}
		elsif ($group =~ /^QC-type13$/) {$color = 0x1;}
		elsif ($group =~ /^QC-type15$/) {$color = 0xF;}
		elsif ($group =~ /^QC-type16$/) {$color = 0x8;}
		else							{$color = 0x4;}

		$pixel[$xpos][$ypos] = $color;

		next;
	}
	
	# is the spot empty?
	if ($seq eq '-')
	{
		# is it fixed?
		if ($fixed eq 'Y')
		{
			# spot is fixed and empty
			$pixel[$xpos][$ypos] = 0x7;
		}
		else
		{
			# spot is empty but not fixed
			$pixel[$xpos][$ypos] = 0x0;
		}
		
		next;
	}
	
	if ($group =~ /^AFFX/)
	{
		# control probes
		$pixel[$xpos][$ypos] = 0xB;
	}
	else
	{
		# "normal" probes
		if ($pm eq "P")
		{
			$pixel[$xpos][$ypos] = 0x3;
		}
		elsif ($pm eq "M")
		{
			$pixel[$xpos][$ypos] = 0x2;
		}
		else
		{
			print "Unexpected probe type at line $line (must be either P or M).\n";
			exit 1;
		}
	}
}

close INFILE;

print $line, " spots described\n";

# ------------------------------------------------------------------------------
# write output file
# ------------------------------------------------------------------------------

open (OUTFILE, ">", $output_file) or
	die "Could not open file '$output_file' for writing: $!.\n";

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
			# empty spot			
			print OUTFILE "0";
		}
	}
	print OUTFILE "\n";
}

close OUTFILE;

print "Chip has ", $ymax + 1, " rows and ", $xmax + 1, " columns\n";
