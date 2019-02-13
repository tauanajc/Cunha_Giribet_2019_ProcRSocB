#!/bin/python
#
# count_species.py
#
# Coded by Christopher E. Laumer, GPL v2
#
# Requires the PyCogent library. Tested with version 1.5.3, python version 2.7.3
#
# This is a simple parser. When run in a directory which contains many multiple alignments named, e.g., OG[N].fa, OG[N+1].fa, it 
# parses all fastas, and outputs a table of species displaying the number of orthogroups in which each species occurs (as well as
# the total number of othrogroups present and the proportion of orthogroups in which that species occurs).  

import glob
from cogent import LoadSeqs, PROTEIN

orthogroups = glob.glob('*.fa')

global_counts = {}
count = 0

for i in orthogroups:
	names = []
	aln = LoadSeqs(i, moltype=PROTEIN)
	names = aln.Names
	for j in names:
		if j not in global_counts.keys():
			global_counts[j] = 1
		elif j in global_counts.keys():
			global_counts[j] += 1
	count += 1

for sp in global_counts:
	print sp, '\t', global_counts[sp], '\t', count, '\t', round(global_counts[sp]/float(count), 2)
		
