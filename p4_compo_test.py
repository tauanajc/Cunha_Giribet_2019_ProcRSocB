#!/usr/bin/env python

"""
Writen by Christopher Laumer.
Edited by Tauana Cunha (https://github.com/tauanajc) on September 2018
	to accommodate more flexible input files.
For each alignment + gene tree input files, runs a p4 simulation test of compositional heterogeneity.
Writes three text files: one with test p-value result,
	and two others with list of alignments that passed or failed test.

Model of evolution has to be manually edited in the script. Currently using LG+gamma4.
p-value/threshold at 0.1
Usage:
p4_compo_test.py -a alignment_file -t tree_file
"""

from p4 import *
from glob import glob
import argparse
import shutil

parser = argparse.ArgumentParser(description='p4 test of compositional heterogeneity based on simulations. \
Takes alignment and gene tree files for one locus. Outputs three files: \
the p-value of the test; two lists of alignments that passed or failed. \
Model of evolution currently has to be edited in the script. Now using LG+gamma4, with a p-value/treshold of 0.1')
parser.add_argument('-a', dest='alignmentfile', required=True, help='input alignment file')
parser.add_argument('-t', dest='treefile', required=True, help='input tree file')
args = parser.parse_args()

read(args.alignmentfile)
a = var.alignments[0]
read(args.treefile)
t = var.trees[0]

d = Data()
t.data = d
t.newComp(free=1, spec='empirical')
t.newRMatrix(free=0, spec='lg')
t.setNGammaCat(nGammaCat=4)
t.newGdasrv(free=1, val=0.5)
t.setPInvar(free=0, val=0.0)
t.optLogLike()

outcome = t.compoTestUsingSimulations(nSims=1000)
handle = open('p4_compotest_outcome.txt','a')
handle.write(args.alignmentfile + '\t' + str(outcome) + '\n')
handle.close()
var.trees,var.alignments=[],[]

#globject = glob(args.alignmentfile + '*')
if outcome < 0.10:
	handle = open('compotest_failed','a')
	handle.write(args.alignmentfile + '\n')
	handle.close()
#	for thing in globject:
#		shutil.move(thing,'compotest_failed')
else:
	handle = open('compotest_passed','a')
	handle.write(args.alignmentfile + '\n')
	handle.close()
#	for thing in globject:
#		shutil.move(thing,'compotest_passed')

print(args.alignmentfile, 'finished')

