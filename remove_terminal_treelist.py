#!/usr/bin/env python

'''
Author: Tauana Cunha | https://github.com/tauanajc
Date: 18 April 2018

Script to remove one species from a list of trees.
Created to mask a species from PhyloBayes treelist files, in case convergence was stuck
based on one terminal difference between chains.
Inputs are treelist file and the name of the terminal to be removed, exactly as it is in the treefile.
Output will be new treelist file with noSPPNAME appended to it.
Usage: python remove_terminal_treelist.py chain1.treelist Crepidula_navicella
'''

import dendropy
import sys

name = sys.argv[1] # Save the name of the input file
spp_delete = sys.argv[2] # Save name of species to delete
new_file = "no" + spp_delete + "_" + name

tree_list = dendropy.TreeList.get(path=sys.argv[1], schema="newick", preserve_underscores=True)

# print(tree_list[0].as_string(schema='newick'))
# tree_list[0].print_plot()

# Prune species
for tree in tree_list:
   tree.prune_taxa_with_labels([spp_delete])
tree_list.write(path=new_file, schema="newick")

