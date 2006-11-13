#!/usr/bin/python
#
# randomchip.py
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
# Create chips with random probes
#
# ------------------------------------------------------------------------------
# To do:
# ==============================================================================

import sys
import random

# alphabet
alpha = 'ACGT'

# base-pair complementarity
base_comp = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

# probe length (number of nucleotides)
probe_len = 25

# Affymetrix deposition sequence
affy_dep_seq = 'TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG'

# middle base position (where PM and MM probes differ)
mid_base = 12

# Deposition sequence for synchronous embeddings (1 base per cycle)
sync_dep_seq = 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT'

# ------------------------------------------------------------------------------
# Create and print on std output a random asynchronous chip (Affymetrix
# deposition sequence) by randomizing a binary embedding vector
# ------------------------------------------------------------------------------
def gen_async_chip_rand_embed (num_rows, num_cols, num_probes):

	# validate arguments
	if (num_rows < 1 or num_cols < 1):
		print "Invalid number of rows/columns:", num_rows, "x", num_cols
		return
	if (num_probes < 1 or num_probes > num_rows * num_cols):
		print "Invalid number of probes:", num_probes
		return

	# number of full and empty spots (if number of
	# probes is less than the total number of spots)
	full = num_probes
	empty = num_rows * num_cols - num_probes

	# use Affy deposition sequence
	dep_seq = affy_dep_seq
	embed_len = len(dep_seq)

	# binary array to generate random embeddings
	embed = [1] * probe_len + [0] * (embed_len - probe_len)

	# fill spots top to bottom, left to right
	for r in range (num_rows):
		for c in range (num_cols):

			# print spot coordinates (column first!)
			sys.stdout.write (str(c) + '\t' + str(r))

			# randomly choose whether the spot is empty or full
			if random.randint (1, empty + full) <= empty:

				# print empty spot
				sys.stdout.write ('\tEMPTY\tN\t-\t-\t-\n')

				empty -= 1

			else:

				# print non-empty spot
				sys.stdout.write ('\tGROUP\tN\t-\t')

				# generate random embedding
				random.shuffle(embed)

				# print probe sequence
				for i in range (embed_len):
					if embed[i] == 1:
						sys.stdout.write (dep_seq[i])

				sys.stdout.write ('\t')

				# print probe sequence
				for i in range (embed_len):
					if embed[i] == 0:
						sys.stdout.write (' ')
					else:
						sys.stdout.write (dep_seq[i])

				sys.stdout.write ('\n')

				full -= 1

# ------------------------------------------------------------------------------
# Create and print on std output a random asynchronous chip (Affymetrix
# deposition sequence) by generating random probe sequences
# ------------------------------------------------------------------------------
def gen_async_chip_rand_probe (num_rows, num_cols, num_probes):

	# validate arguments
	if (num_rows < 1 or num_cols < 1):
		print "Invalid number of rows/columns:", num_rows, "x", num_cols
		return
	if (num_probes < 1 or num_probes > num_rows * num_cols):
		print "Invalid number of probes:", num_probes
		return

	# number of full and empty spots (if number of
	# probes is less than the total number of spots)
	full = num_probes
	empty = num_rows * num_cols - num_probes

	# use Affy deposition sequence
	dep_seq = affy_dep_seq
	embed_len = len(dep_seq)

	# pool of bases for generating random probes
	probe = ['A'] * probe_len + ['C'] * probe_len + ['G'] * probe_len + ['T'] * probe_len

	# binary array to store left-most embedding
	embed = [0] * embed_len

	# fill spots top to bottom, left to right
	for r in range (num_rows):
		for c in range (num_cols):

			# print spot coordinates (column first!)
			sys.stdout.write (str(c) + '\t' + str(r))

			# randomly choose whether the spot is empty or full
			if random.randint (1, empty + full) <= empty:

				# print empty spot
				sys.stdout.write ('\tEMPTY\tN\t-\t-\t-\n')

				empty -= 1

			else:

				# print non-empty spot
				sys.stdout.write ('\tGROUP\tN\t-\t')

				ready = False
				while not ready:
					# generate random probe
					random.shuffle(probe)

					# compute left-most embedding (if possible)
					base = 0
					step = 0
					while base < probe_len and step < embed_len:
						if probe[base] == dep_seq[step]:
							embed[step] = 1
							base += 1
						else:
							embed[step] = 0
						step += 1

					if base == probe_len:
						ready = True
						# turn off remaining steps
						while step < embed_len:
							embed[step] = 0
							step += 1

				# print probe sequence
				for i in range (embed_len):
					if embed[i] == 1:
						sys.stdout.write (dep_seq[i])

				sys.stdout.write ('\t')

				# print probe sequence
				for i in range (embed_len):
					if embed[i] == 0:
						sys.stdout.write (' ')
					else:
						sys.stdout.write (dep_seq[i])

				sys.stdout.write ('\n')

				full -= 1

# ------------------------------------------------------------------------------
# Create and print on std output a random simple chip definition with probes
# synchronously-embedded on a 100-base deposition sequence
# ------------------------------------------------------------------------------
def gen_sync_chip (num_rows, num_cols, num_probes):

	# validate arguments
	if (num_rows < 1 or num_cols < 1):
		print "Invalid number of rows/columns:", num_rows, "x", num_cols
		return
	if (num_probes < 1 or num_probes > num_rows * num_cols):
		print "Invalid number of probes:", num_probes
		return

	# number of full and empty spots (if number of
	# probes is less than the total number of spots)
	full = num_probes
	empty = num_rows * num_cols - num_probes

	# use 100-base deposition sequence
	dep_seq = sync_dep_seq
	embed_len = len(dep_seq)

	# binary array to generate random cycles
	embed = [0] * embed_len
	rnd_cycle = [1, 0, 0, 0]
	len_cycle = len (rnd_cycle)
	num_cycles = embed_len / len_cycle

	# fill spots top to bottom, left to right
	for r in range (num_rows):
		for c in range (num_cols):

			# print spot coordinates (column first!)
			sys.stdout.write (str(c) + '\t' + str(r))

			# randomly choose whether the spot is empty or full
			if random.randint (1, empty + full) <= empty:

				# print empty spot
				sys.stdout.write ('\tEMPTY\tN\t-\t-\t-\n')

				empty -= 1

			else:

				# print non-empty spot
				sys.stdout.write ('\tGROUP\tN\t-\t')

				# generate random embedding
				for c in range (num_cycles):
					random.shuffle(rnd_cycle)
					embed[c * len_cycle:(c + 1) * len_cycle] = rnd_cycle

				# print probe sequence
				for i in range (embed_len):
					if embed[i] == 1:
						sys.stdout.write (dep_seq[i])

				sys.stdout.write ('\t')

				# print probe sequence
				for i in range (embed_len):
					if embed[i] == 0:
						sys.stdout.write (' ')
					else:
						sys.stdout.write (dep_seq[i])

				sys.stdout.write ('\n')

				full -= 1

# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------
# simple_chip (300, 300, 300 * 300)
gen_async_chip_rand_probe (200, 200, 200 * 200)