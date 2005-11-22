#!/usr/bin/perl
#
# plembed.pl
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
# Compute a "pair-aware" version of the left-most embedding.
#
# Input:  <chip code>_<layout>	       (chip's layout specification)
#         <chip code>_depseq           (deposition sequence)
# Output: <chip code>_<layout>_plembed (layout with probe embedding)
#
# The layout file contains a list of all empty and non-empty spots in the chip.
# Non-empty are either normal probes, control probes (those whose group name
# start with AFFX) or QC (quality control) probes. Empty spots may have special
# functions in the production or hybridization phases.
#
# Only normal and control probes have full sequence information. Therefore, only
# these spots are considered for the embedding, while the others (like QC probes)
# are just regarded as empty.
#
# The output is the PLEMBED file which contains all information contained in the
# layout file plus the corresponding embeddings of all probes whose sequence is
# known.
#
# Probes appear in polymorphic pairs (PM/MM) that are exactly the same except
# for the base in the middle position, where one probe has the complement base
# of the other. For this reason, the first half of both probes are embedded in
# the same way but, starting from the middle position, the embeddings can differ
# a lot. Since PM and MM probes are always located next to each other the
# differences in their second half embeddings can increase the border length and
# the conflict index of the corresponding spots. To avoid this problem, this
# algorithm has been modified to take this into account and embed the sequence
# pairs optimally. With this in mind, the first half of the probes are embedded
# with the normal left-most embedding algorithm. The middle probe is also left-
# most embedded but the second half of the probes are embedded only after a
# point in the deposition sequence where the complement of the middle base can
# also be embedded.
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
	print "       input files: <chip code>_<layout>         (chip's layout specification)\n";
	print "                    <chip code>_depseq           (deposition sequence)\n";
	print "       output file: <chip code>_<layout>_plembed (layout with probe embedding)\n";
	exit 1;
}

# file names
my $chip_file = $chip_code . "_" . $layout_code;
my $depseq_file = $chip_code . "_depseq";
my $output_file = $chip_code . "_" . $layout_code . "_plembed";

# ------------------------------------------------------------------------------
# check if output file already exists
# ------------------------------------------------------------------------------

my $tmp;

if (-f $output_file)
{
	print "Output file '$output_file' already exists. Overwrite? (Y/N): ";
	chomp($tmp = <STDIN>);
	
	($tmp !~ m/^[Yy]/) and exit 0;
}

# ------------------------------------------------------------------------------
# read deposition sequence
# ------------------------------------------------------------------------------

open (INFILE, $depseq_file) or
	die "Could not open file '$depseq_file' for reading: $!.\n";

# deposition sequence
my $dep_seq = <INFILE>;

close INFILE;

# remove NL and CR chars, if any
$dep_seq =~ s/\n|\r//g;

# validates the sequence
($dep_seq !~ /[ACGT]+/) and
	die "Invalid deposition sequence: ($dep_seq).\n";

# length of the deposition sequence
my $len = length $dep_seq;

print "Deposition sequence: $dep_seq (length $len)\n";

# ------------------------------------------------------------------------------
# read input file
# ------------------------------------------------------------------------------

my ($x, $y, $group, $fixed, $pm, $sequence, $embed);
my ($ln, $count_empty, $count_embed) = (0, 0, 0);

open (INFILE, $chip_file) or
	die "Could not open file '$chip_file' for reading: $!.\n";

open (OUTFILE, '>', $output_file) or
	die "Could not open file '$output_file' for writing: $!.\n";

print "Reading '$chip_file'\n";

for (<INFILE>)
{
	$ln++;
	
	($x, $y, $group, $fixed, $pm, $sequence) = split /\t|\n|\r\n/;
	
	if ($sequence eq "-")
	{
		# empty spot: nothing to do
		$embed = "-";
		$count_empty++;
	}
	elsif ($sequence =~ /^[ACGT]+$/)
	{
		# probe sequence: embed
		$embed = embed ($sequence);
		$count_embed++;
	}
	else
	{
		print "Invalid sequence in line $ln: $sequence.\n";
		exit 1;
	}
	
	print OUTFILE "$x\t$y\t$group\t$fixed\t$pm\t$sequence\t$embed\n";
}

close INFILE;
close OUTFILE;

print "$count_empty empty sequences.\n";
print "$count_embed sequences embedded.\n";
print "Output written to '$output_file'\n";

# ------------------------------------------------------------------------------
# compute and return a "pair-aware" left-most embedding of a sequence
# ------------------------------------------------------------------------------
sub embed
{
	my $seq = $_[0];
	my ($s, $p, $mflag, $cflag) = (0, 0, 0, 0);
	my $embed = "";
	
	# position of the middle base
	my $mpos = int((length $seq)/ 2);
	
	# middle base
	my $mbase = substr($seq, $mpos, 1);
	
	# and its complement (which is the
	# midde base of the sequence's mate pair)
	my $cbase = $compl{$mbase};

	# embed the first half of the sequence
	for (; $s < $mpos; $p++)
	{
		if (substr($seq, $s, 1) eq substr($dep_seq, $p, 1))
		{
			$embed = $embed . substr($seq, $s, 1);
			$s++;
		}
		else
		{
			$embed = $embed . " ";
		}
	}
	
	# continue until both the base and its complement can be embedded
	for (; $mflag * $cflag == 0; $p++)
	{
		if (($mflag == 0) and ($mbase eq substr($dep_seq, $p, 1)))
		{
			# embed the middle base
			$embed = $embed . substr($seq, $s, 1);
			$s++;
			
			# flag that middle base has been embedded
			$mflag = 1;
		}
		else
		{
			$embed = $embed . " ";

			if (($cflag == 0) and ($cbase eq substr($dep_seq, $p, 1)))
			{
				# flag that the complement base can already be embedded
				$cflag = 1;
			}
		}
	}
	
	# embed the rest of the sequence
	for (; $p < $len; $p++)
	{
		if (($s < length $seq) and (substr($seq, $s, 1) eq substr($dep_seq, $p, 1)))
		{
			$embed = $embed . substr($seq, $s, 1);
			$s++;
		}
		else
		{
			$embed = $embed . " ";
		}
	}
	
	return $embed;
}