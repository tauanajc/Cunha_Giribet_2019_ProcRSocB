#!/usr/bin/python

'''
Author: Tauana Cunha
Date: 04 May 2017
Converts nexus file to phylip format. For amino acid data.
Does not need to be inside folder with input/output.
Usage:
python Nexus2Phylip_AA.py input.nex
'''

import dendropy
import sys

name,extension = sys.argv[1].split('.') # Save the name of the input dataset
name = name + '.phy' # Save string for output in phylip

dataset = dendropy.ProteinCharacterMatrix.get(path=sys.argv[1], schema="nexus", preserve_underscores=True)
dataset.write(path=name, schema="phylip")

