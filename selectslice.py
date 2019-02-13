#!/usr/bin/env python

"""
Created in the Giribet Lab at Harvard University. Modified by Tauana Cunha on Apr 13, 2017. http://github.com/tauanajc
This script creates a folder with all the OMA orthogroups files (or already aligned OG) based on the desired matrix occupancy.
For example, for a 50% occupancy matrix, you should give the X number of taxa that corresponds to half of your sample. It will then
copy all the files for OGs that have X or more taxa present.
It can ask for the number of taxa interactively, or you can give it in the command line of python - change commented line below (j=).

Launch from inside folder with all the fasta files. Need to have .fa in name, but can have other parts after, e.g. .fa.mafft
Usage:
python selectslice.py
python selectslice.py 33
Where 33 is the X number of taxa for the occupancy matrix.
"""

from Bio import SeqIO
import glob
import os
import sys

#j = int(raw_input('Choose a minimum taxon occupancy: ')) # For interactive way to give number of taxa
j = int(sys.argv[1]) # For number of taxa given as argument after name of python script

all_files = glob.glob('*.fa*')

dirname = 'OGslice_' + str(j) + 'taxa'

os.mkdir(dirname)

for i in all_files:
        parsed = SeqIO.to_dict(SeqIO.parse(i, 'fasta'))
        if len(parsed) >= j:
                select = True
#        elif len(parsed) > 4:     # If want to include certain poor taxa, include genes poorly sampled that have that taxa
#                if 'Gnat' in parsed:
#                        select = True
#                else:
#                        select = False
        else:
                select = False
        if select:
                os.chdir(dirname)
                for seq in parsed:
                        out = open(i, 'a')
                        out.write('>' + seq + '\n' + str(parsed[seq].seq) + '\n')
                        out.close()
                os.chdir('..')
