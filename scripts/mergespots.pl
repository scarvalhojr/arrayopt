#!/usr/bin/perl
#
# mergespots.pl
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
# Merge information from TAB and QC files.
#
# Input:  <chip code>_tab	(tabular list of probes from affymetrix.com)
#         <chip code>_qc	(tabular list of QC probes extracted from CDF file)
# Output: <chip code>_spots	(tabular list of all non-empty spots in the chip)
#         <chip code>_empty	(tabular list of all empty spots in the chip)
#
# The TAB file must be downloaded from affymetrix.com. It contains the list of
# all "normal" probes in the chip in a tabular format. The position of these
# probes in the chip is considered to be NOT fixed (relocation is possible).
# These probes always appear in a PM/MM (perfect match/mismatch) pair, but
# only PM probes are listed in the TAB file.
#
# The QC file is generated from the CDF file. It contains a list of all QC
# (quality control) probes. Their position are considered to be fixed in the
# chip (no relocation is possible). These probes do not appear in PM/MM pairs
# and, their sequence is not publicly available.
#
# The SPOTS file contains a list of all non-empty spots in the chip. These are
# either normal probes (from the TAB file) or QC probes (from the QC file).
# Note that PM/MM probe pairs are created from the PM probe information
# contained in the TAB file.
#
# All spots which are not assigned a probe are regarded as empty and are listed
# in the EMPTY file. These spots are initially regarded as fixed (i.e. no probe
# could be moved to these positions). The EMPTY file, however, is the source for
# the GAPS file. Therefore, some spots can be later turned into non-fixed by
# just removing their corresponding records in the GAPS file.
#
# ------------------------------------------------------------------------------
# To do:
# ==============================================================================

use strict;
use warnings;

# ------------------------------------------------------------------------------
# global variables
# ------------------------------------------------------------------------------

# base-pair complementarity table
my %compl = ("A" => "T", "C" => "G", "G" => "C", "T" => "A");

# number of normal probes
my $normal_probes = 0;

# number of control probes
my $control_probes = 0;

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
	print "       input files are:  <chip code>_tab    (sequences in tabular format)\n";
	print "                         <chip code>_qc     (QC probes extracted from the CDF file)\n";
	print "       output files are: <chip code>_spots  (tabular list of normal and QC probes)\n";
	print "                         <chip code>_empty  (list of empty spots)\n";
	exit 1;
}

# file names
my $tab_file = $chip_code . "_tab";
my $qc_file = $chip_code . "_qc";
my $spots_file = $chip_code . "_spots";
my $empty_file = $chip_code . "_empty";

# ------------------------------------------------------------------------------
# check if output files already exist
# ------------------------------------------------------------------------------

if (-f $spots_file)
{
	print "Output file '$spots_file' already exists. Overwrite? (Y/N): ";
	chomp($tmp = <STDIN>);
	
	($tmp !~ m/^[Yy]/) and exit 0;
}

if (-f $empty_file)
{
	print "Output file '$empty_file' already exists. Overwrite? (Y/N): ";
	chomp($tmp = <STDIN>);
	
	($tmp !~ m/^[Yy]/) and exit 0;
}

# ------------------------------------------------------------------------------
# process TAB file
# ------------------------------------------------------------------------------

# probe information
my ($group, $xpos, $ypos, $sequence);

# chip maximum coordinates
my $xmax = 0;
my $ymax = 0;

# spot flags (mark used
# and unused spots)
my @spot;

# open TAB file for reading
open (INFILE, $tab_file) or
	die "Could not open file '$tab_file' for reading: $!.\n";

# open main output file
open (OUTFILE, ">", $spots_file) or
	die "Could not open file '$spots_file' for writing: $!.";

print "Reading input file '$tab_file'\n";

# skip first line (header)
<INFILE>;

# line number
my $line = 1;

for (; <INFILE>; $line++)
{
	($group, $xpos, $ypos, $tmp, $sequence) = split /\t|\n|\r\n/;
	
	if (!($xpos =~ /[0-9]+/ && $ypos =~ /[0-9]+/) || $xpos < 0 || $ypos < 0)
	{
		print "Invalid position ($xpos, $ypos) at line $line of '$tab_file'.\n";
		exit 1;
	}
	
	# keep track of maximum values
	# probes from this file always occur in pairs
	# (in consecutive lines, with equal x position)
	$xmax = ($xpos > $xmax) ? $xpos : $xmax;
	$ymax = ($ypos + 1 > $ymax) ? $ypos + 1: $ymax;

	# check if spot pair is already in use
	no warnings 'uninitialized';
	if ($spot[$xpos][$ypos] > 0 || $spot[$xpos][$ypos + 1] > 0)
	{
		print "Spot conflict for the pair at x = $xpos, y = $ypos (and y = $ypos + 1).\n";
		exit 1;
	}
	use warnings 'uninitialized';

	# mark spots (PM/MM) as used
	$spot[$xpos][$ypos] = 1;
	$spot[$xpos][$ypos + 1] = 1;
	
	# print probe pair info
	# X, Y, Group, Fixed, P/M, Sequence
	
	# PM probe
	print OUTFILE "$xpos\t$ypos\t$group\tN\tP\t$sequence\n";
	
	$ypos++;
	my $m = int((length $sequence) / 2);
	$sequence = substr($sequence, 0, $m) . $compl{substr($sequence, $m, 1)} . substr($sequence, $m + 1);
	
	# MM probe
	print OUTFILE "$xpos\t$ypos\t$group\tN\tM\t$sequence\n";
	
	# count normal and control probes
	if ($group =~ /^AFFX/)
	{
		$control_probes += 2;
	}
	else
	{
		$normal_probes += 2;
	}
}

close INFILE;

print $line - 1, " probe pairs\n";
print $normal_probes, " normal probes\n";
print $control_probes, " control probes\n";

# ------------------------------------------------------------------------------
# process QC file (quality control probes)
# ------------------------------------------------------------------------------

open (INFILE, $qc_file) or
	die "Could not open file '$qc_file' for reading: $!.\n";

print "Reading input file '$qc_file'\n";

for ($line = 1; <INFILE>; $line++)
{
	($xpos, $ypos, $group) = split /\t|\n|\r\n/;
	
	if (!($xpos =~ /[0-9]+/ && $ypos =~ /[0-9]+/) || $xpos < 0 || $ypos < 0)
	{
		print "Invalid position ($xpos, $ypos) at line $line of '$qc_file'.\n";
		exit 1;
	}
	
	# keep track of maximum values
	$xmax = ($xpos > $xmax) ? $xpos : $xmax;
	$ymax = ($ypos > $ymax) ? $ypos : $ymax;

	# check if spot is already in use
	no warnings 'uninitialized';
	if ($spot[$xpos][$ypos] > 0)
	{
		print "Spot conflict for the QC probe at x = $xpos, y = $ypos.\n";
		exit 1;
	}
	use warnings 'uninitialized';

	# mark spot as used
	$spot[$xpos][$ypos] = 1;
	
	# print probe info
	# X, Y, Group, Fixed, P/M, Sequence
	print OUTFILE "$xpos\t$ypos\t$group\tY\t-\t-\n";
}

close INFILE;

close OUTFILE;

print $line - 1, " QC probes\n";

# ------------------------------------------------------------------------------
# output list of empty spots
# ------------------------------------------------------------------------------

# number of empty spots
my $count = 0;

open (OUTFILE, ">", $empty_file) or
	die "Could not open EMPTY file '$empty_file' for writing: $!.\n";

no warnings 'uninitialized';

for (my $x = 0; $x <= $xmax; $x++)
{
	for (my $y = 0; $y <= $ymax; $y++)
	{
		if (!($spot[$x][$y] > 0))
		{
			# X, Y, Group, Fixed, P/M, Sequence
			printf OUTFILE "$x\t$y\tEMPTY\tY\t-\t-\n";
			$count++;
		}
	}
}

close OUTFILE;

print "$count empty spots\n";
