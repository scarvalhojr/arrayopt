#!/usr/bin/perl
#
# makechip.pl
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
# Create a full specification of the chip's original (Affymetrix) layout.
#
# Input:  <chip code>_spots	 (tabular list of all non-empty spots in the chip)
#         <chip code>_blanks (tabular list of all spots that must remain empty)
#         <chip code>_fixed  (regular expression for groups with fixed spots)
# Output: <chip code>_affy	 (full specification of the original chip layout)
#
# The SPOTS file contains a list of all non-empty spots in the chip. These are
# either normal probes or QC (quality control) probes.
#
# The BLANKS file contains the list of all spots that must be empty in the chip
# (these spots may be empty for some special reason during the Affymetrix
# production process or during the hybridization phase). Note that probe can be
# assigned to these positions (in an attempt to improve the layout). Usually
# this file is manually produced from the EMPTY file since it is difficult to
# programatically separate the spots that must empty from those who are empty
# just because there is not enough probes to fill all spots in the chip.
#
# If the BLANKS file do not exist, the program will prompt the user to use the
# empty files as to generate a first layout specification that can be later 
# changed. This is useful to allow the user an earlier creation of a visual map
# so that the real fixed empty spots can be identified.
#
# Any spot not listed in both input files are also considered empty, however,
# they can be assigned to (relocated) probes in an attempt to improve the chip's
# layout.
#
# The fixed files is optional. If it exists, it must contans a regular
# expression that will match group names whose probes positions are fixed on the
# chip, i.e. these probes can not be moved to other positions in an attempt to
# improve the chip's layout.
#
# The output file CHIP is merely a concatenation of both input files.
#
# ------------------------------------------------------------------------------
# To do:
# ==============================================================================

use strict;
use warnings;

# ------------------------------------------------------------------------------
# global variables
# ------------------------------------------------------------------------------

# line cout
my $lines;

# fixed groups patter
my $fixed_group;

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
	print "       input files are: <chip code>_spots  (tabular list of all non-empty spots in the chip)\n";
	print "                        <chip code>_blanks (tabular list of all spots that must be empty)\n";
	print "                        <chip code>_fixed  (optional regular expression for groups with fixed spots)\n";
	print "       output file is:  <chip code>_affy   (full specification of the original chip layout)\n";
	exit 1;
}

# file names
my $spots_file = $chip_code . "_spots";
my $empty_file = $chip_code . "_empty";
my $blanks_file = $chip_code . "_blanks";
my $fixed_file = $chip_code . "_fixed";
my $chip_file = $chip_code . "_affy";

# ------------------------------------------------------------------------------
# check if output file already exists
# ------------------------------------------------------------------------------

if (-f $chip_file)
{
	print "Output file '$chip_file' already exists. Overwrite? (Y/N): ";
	chomp($tmp = <STDIN>);
	
	($tmp !~ m/^[Yy]/) and exit 0;
}

# ------------------------------------------------------------------------------
# read FIXED regular expression
# ------------------------------------------------------------------------------

if (-f $fixed_file)
{
	# file exists, try to read it
	open (INFILE, $fixed_file) or
		die "Could not open file '$fixed_file' for reading: $!.\n";

	$fixed_group = <INFILE>;
	
	# remove NL and CR chars, if any
	$fixed_group =~ s/\n|\r//g;
	
	close INFILE;
}
else
{
	# warn user
	print "Fixed groups definition not found. Continue anyway? (Y/N): ";
	chomp($tmp = <STDIN>);
	
	($tmp !~ m/^[Yy]/) and exit 0;

	# no group is fixed	
	$fixed_group = qr/^$/;
}

# ------------------------------------------------------------------------------
# open input and output files
# ------------------------------------------------------------------------------

if (-f $blanks_file)
{
	# BLANKS file exists, try to read it
	open (BLANKS_FILE, $blanks_file) or
		die "Could not open file '$blanks_file' for reading: $!.\n";
}
else
{
	# use EMPTY instead?
	print "BLANKS file not found. Use EMPTY instead? (Y/N): ";
	chomp($tmp = <STDIN>);
	
	($tmp !~ m/^[Yy]/) and exit 0;

	# use EMPTY instead of BLANKS
	open (BLANKS_FILE, $empty_file) or
		die "Could not open file '$empty_file' for reading: $!.\n";
}

# open SPOTS file for reading
open (SPOTS_FILE, $spots_file) or
	die "Could not open file '$spots_file' for reading: $!.\n";

# open main output file
open (OUTFILE, ">", $chip_file) or
	die "Could not open file '$chip_file' for writing: $!.";

# ------------------------------------------------------------------------------
# copy SPOTS file
# ------------------------------------------------------------------------------

my ($x, $y, $group, $fixed, $pm, $sequence, $embed);
my ($fixed_count, $nonfixed_count) = (0, 0);

for ($lines = 0; <SPOTS_FILE>; $lines++)
{
	($x, $y, $group, $fixed, $pm, $sequence) = split /\t|\n|\r\n/;
	
	# if probe's group name match a fixed group
	if ($group =~ /$fixed_group/)
	{
		# make spot fixed
		$fixed = "Y";
		$fixed_count ++;
	}
	else
	{
		$nonfixed_count ++;
	}
	
	print OUTFILE "$x\t$y\t$group\t$fixed\t$pm\t$sequence\n";
}

print "$lines spots ($fixed_count fixed, $nonfixed_count not fixed)\n";

close SPOTS_FILE;

# ------------------------------------------------------------------------------
# copy BLANKS file
# ------------------------------------------------------------------------------

for ($lines = 0; <BLANKS_FILE>; $lines++)
{
	print OUTFILE;
}

print "$lines blanks\n";

close SPOTS_FILE;

close OUTFILE;
