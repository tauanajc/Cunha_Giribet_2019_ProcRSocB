#!/usr/bin/env python

## written by B. Medeiros, starting on 20-Apr-2017
## souzademedeiros@fas.harvard.edu https://github.com/brunoasm
## this script takes as input a bunch of fasta files and outputs an occupancy matrix

import argparse, pandas, os

parser = argparse.ArgumentParser(description='This program takes as input fasta or phylip files, each one representing a gene. It outputs a matrix with taxa on columns and genes on rows, showing whether a taxon is present in the file for each gene. This matrix can be easily loaded on R to print a graphical occupancy matrix.')
parser.add_argument('input', help = 'path to all fasta or phylip files (example: *.fasta)', nargs = '*')
parser.add_argument('-o','--output',default = 'occ.matrix.tsv',help = 'path to output file (default: occ.matrix.tsv)')

args = parser.parse_args()

occupancy = list()

for seqpath in args.input:
    seqfile = open(seqpath,'r')
    taxa_present = set()
    
    #read first line to check if fasta or phylip
    first_line = seqfile.readline()
    if '>' in first_line:
        filetype = 'fasta'
        taxa_present.add(first_line[1:-1])
    else:
        filetype = 'phylip'
    
    #now read line by line and get taxon names according to the file format        
    for line in seqfile:
        if filetype == 'fasta' and line[0] == '>':
            taxa_present.add(line[1:-1])
        elif filetype == 'phylip':
            taxa_present.add(line.split()[0])
    #save taxa present in a dict            
    tempdict = {taxon:int(1) for taxon in taxa_present}
    tempdict.update({'filename':os.path.basename(seqpath)})
    
    #append dict in a list
    occupancy.append(tempdict)
    seqfile.close()

#use pandas to make a dataframe
occupancy_df = pandas.DataFrame(occupancy)
#replace NaNs with 0
occupancy_df.fillna(int(0), inplace = True)

#move filename to the first column
cols = occupancy_df.columns.tolist()
cols.remove('filename')
cols.insert(0, 'filename')
occupancy_df = occupancy_df[cols]

#save resulting matrix
occupancy_df.to_csv(args.output,sep='\t',float_format='%d',index=False)
