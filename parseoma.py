#!/usr/bin/env python

'''
Created in the Giribet Lab at Harvard University.
Reads OMA output and prints a table with the number of orthologs containing an increasing number of taxa.
Input is OrthologousGroups.txt file from OMA output.

Usage:
python parseoma.py -i OrthologousGroups.txt
'''

import argparse

def parser():
        args = argparse.ArgumentParser()
        args.add_argument('-i', '--input', required=True, help='The path of the OrthologousGroups.txt file output by OMA.')
        args = args.parse_args()
        return args

def main():
        Ortho = open(parser().input, 'r')
        max = 0
        for i in Ortho:
                prelist = i.split('\t')[1:]
                if len(prelist) > max:
                        max = len(prelist)
        Ortho.close()
        D = {}
        while max > 1:
                D[max] = 0
                max -= 1
        Ortho = open(parser().input, 'r')
        for line in Ortho:
                if '#' in line:
                        pass
                else:
                        list = line.split('\t')[1:]
                        D[len(list)] += 1

        print ('---------')
        print 'Number of orthogroups', '\t', 'N taxa in each'
        for k in range(2,len(D)+2):
                print str(D[k]),'\t', str(k)

main()
