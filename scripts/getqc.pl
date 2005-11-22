#!/usr/bin/perl
#
# getqc.pl
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
# Extract QC probes from a CDF file.
#
# Input:  <chip code>_cdf	(CDF file from affymetrix.com)
# Output: <chip code>_qc	(tabular list of QC probes)
#
# The CDF file must be downloaded from affymetrix.com. It contains full layout
# specification of a chip, although there is no probe sequence data. Here we are
# only interested in extracting the information about the QC (quality control)
# probes - those that create, among other things, the checkerboard pattern.
#
# The QC file contains a list of all QC probes in a tabular format. These
# probes are grouped according to their type.
#
# ------------------------------------------------------------------------------
# To do:
# ==============================================================================

use strict;
use warnings;

# ------------------------------------------------------------------------------
# global variables
# ------------------------------------------------------------------------------

# probe information
my @register;

# temp variable
my $tmp;

# ------------------------------------------------------------------------------
# check arguments
# ------------------------------------------------------------------------------

# chip code: prefix for input and output files
my $chip_code;

# $#ARGV contains the index of the last element of @ARGV
if ($#ARGV == 0)
{
	# get chip code
	$chip_code = $ARGV[0];
}
else
{
	print "Usage: $0 <chip code>\n";
	print "       input : <chip code>_cdf (CDF file from affymetrix.com)\n";
	print "       output: <chip code>_qc  (tabular list of QC probes)\n";
	exit 1;
}

my $input_file = $chip_code . "_cdf";
my $output_file = $chip_code . "_qc";

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

my ($group, $xpos, $ypos, $ncells, $type); 
my $state = 0;
my $ln = 0;
my $reg = 0;
my $count = 0;

open (INFILE, $input_file) or
	die "Could not open file '$input_file' for reading: $!.\n";

print "Reading input file '$input_file'\n";

while (my $line = <INFILE>)
{
	$ln++;
	
	if (($state == 0) and ($line =~ /^\[QC[0-9]+\]/))
	{
		$group = substr($line, 1, index($line, "]") - 1);
		
		$line = <INFILE>; $ln++;
		
		if ($line =~ /^Type=[0-9]/)
		{
			($tmp, $type) = split /=|\t|\n|\r\n/, $line;
		}
		else
		{
			print "Expected 'Type' entry at line $ln (group $group), found '$line'.";
			exit 1;
		}
		
		$line = <INFILE>; $ln++;

		if ($line =~ /^NumberCells=/)
		{
			($tmp, $ncells) = split /=|\t|\n|\r\n/, $line;
			
			if ($ncells < 1)
			{
				print "Invalid number of cells for group $group: $ncells.\n";
				exit 1;
			}
		}
		else
		{
			print "Expected 'NumberCells' entry at line $ln (group $group), found '$line'.";
			exit 1;
		}
		
		$line = <INFILE>; $ln++;

		if (!($line =~ /^CellHeader=X\tY\t/))
		{
			print "Expected 'CellHeader=X Y' entry at line $ln (group $group), found '$line'.";
			exit 1;
		}
		
		print "Group $group: type=$type, ncells=$ncells:\n";
		
		$state = 1;
		$count = 0;
		next;
	}
	
	if ($state == 1)
	{
		if ($line =~ /^Cell[0-9]+=[0-9]+\t[0-9]+\t/)
		{
			($tmp, $xpos, $ypos) = split /=|\t|\n|\r\n/, $line;
			$count++;
			
			# save register for later output
			$register[$reg][0] = $xpos;
			$register[$reg][1] = $ypos;
			$register[$reg][2] = $type;
			$reg++;
		}
		else
		{
			# end of group: check number of cells
			if ($count != $ncells)
			{
				print "Wrong number of cells for group $group (type $type): expected $ncells, found $count.\n";
				exit 1;
			}
			
			$state = 0;
		}
	}
}

# make sure last group is also checked for the number of cells
if ($count != $ncells)
{
	print "Wrong number of cells for group $group (type $type): expected $ncells, found $count.\n";
	exit 1;
}

close INFILE;

# ------------------------------------------------------------------------------
# write output file
# ------------------------------------------------------------------------------

open (OUTFILE, '>', $output_file) or
	die "Could not open file '$output_file' for writing: $!.\n";

print "Writing QC file '$output_file'\n";

for (my $r = 0; $r < $reg; $r++)
{
	print OUTFILE "$register[$r][0]\t$register[$r][1]\tQC-type$register[$r][2]\n";
}

print "$reg QC probes\n";

close OUTFILE;

exit 0;
