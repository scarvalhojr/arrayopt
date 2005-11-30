#! /bin/awk -f
#
# affy2simple.awk
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
# Creates a simple chip from an Affymetrix chip
#
# This script reads a layout specification of an Affymetrix chip (one that can
# serve as input to the AffymetrixChip class) and outputs a similar chip layout
# that can serve as input to the SimpleChip class. This is done by removing rows
# of MM probes. Moreover, the PM/MM flags are set to '-' since these flags have
# no meaning for a SimpleChip instance.
#
# The removal of MM probes is achieved by removing spots with EVEN row number.
# As a result, the chip's vertical dimension is reduced by nearly a half and the
# row numbers of the remaining spots are mapped to the new grid.
#
# Some chips, however, may have a small rectangular area where PM probes are
# located on spots with EVEN row number. To avoid keeping the MM probes of such
# regions, it is possible to provide its exact borders. Then, inside these
# regions, the script will remove spots with ODD row numbers. If this feature
# is not needed, the no_box variable must be set to 1.
#
# Usage: affy2simple.awk -v no_box=1
#     OR
#        affy2simple.awk -v box_top=<row number> -v box_bottom=<row number>
#        -v box_left=<column number> -v box_right=<column number>
#
# ------------------------------------------------------------------------------
# To do:
# ==============================================================================

BEGIN {
	# input field separator
	FS = "\t"
	
	# output field separator
	OFS = "\t"
	
	if (no_box == 1) {

		# no inverted box
		box_left = -1
		box_right = -1
		box_top = -1
		box_bottom = -1

	} else if (box_left <= 0 || box_right <= 0 || box_top <= 0 || box_bottom <= 0) {
	
		print "ERROR: box borders invalid or not specified."
		print "Usage: affy2simple.awk -v no_box=1"
		print "    OR"
		print "       affy2simple.awk -v box_top=<row number> -v box_bottom=<row number>"
		print "       -v box_left=<column number> -v box_right=<column number>"
		exit
	}
}

# row zero is kept intact
($2 > 0) {

	if ($1 >= box_left && $1 <= box_right && $2 >= box_top && $2 <= box_bottom) {
	
		# inverted spots

		# skip lines with ODD row numbers
		if ($2 % 2 != 0) next
	
		# update lines with EVEN row number
		$2 = $2 / 2
				
	} else {
	
		# normal spots

		# skip lines with EVEN row numbers
		if ($2 % 2 == 0) next

		# update lines with ODD row number
		$2 = 1 + (($2 - 1) / 2)
	}
	
	# reset PM/MM flag
	$5 = "-"
}

{
	# print line with new row number
	print $0
}
