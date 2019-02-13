#!/bin/bash
'''
Wrote by Tauana Cunha on August 11, 2017 to work on protein alignments and cut only
positions with more than 80% missing data.
https://github.com/tauanajc
Based on a script for nucleotides by 5heikki in biostars (https://www.biostars.org/p/172496/).

Script to trim ends of alignments. By default, cuts all positions in the ends with more
than 80% missing data.
trimEnds.sh alignment.fa
'''

# Separate the alignment from the names
#cut -f1 -d " " $1 > $1.names # For alignments separated from fasta name by space
#cut -f2 -d " " $1 > $1.ali
grep '>' $1 > $1.names # For alignments separated from fasta name by line break
grep -v '>' $1 > $1.ali

# Set start at 1
START=1

# Set end at alignment length
END=$(head -n1 $1.ali | awk '{print length}')

# Set number of taxa
taxa=$(cat $1.names | wc -l)

# Find out where the good stuff starts
for (( i=$START; i<=$END; i++ )) # For 1 to length of alignment, increasing one by one
    do # Breaks each position and counts how many non Letter in that position
    nonLetterCount=$(cut -c$i $1.ali | tr -d [A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V] | grep . | grep -c .)
    missing=$(bc -l <<< "$nonLetterCount/$taxa*100" | cut -f1 -d.) # Number of non Letter divided by total taxa
    if [ "$missing" -gt 80 ] # If amount of non Letter in that position is more than 80%, will eliminate that position…
        then
        (( START = $i + 1 )) # … and push START of sequence over
    else
        break
    fi
    done

# Find out where the good stuff ends
for (( i=$END; i>=$START; i-- )) # Same thing, but decreasing number
    do
    nonLetterCount=$(cut -c$i $1.ali | tr -d [A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V] | grep . | grep -c .)
    missing=$(bc -l <<< "$nonLetterCount/$taxa*100" | cut -f1 -d.)
    if [ "$missing" -gt 80 ]
        then
        (( END = $i - 1 ))
    else
        break
    fi
    done

# Reunite names with seqs
name=$(basename $1 | cut -f 1 -d .)
paste -d "\n" $1.names <(cut -c$START-$END $1.ali) > $name.trimmed.fa

# Remove unnecessary files
rm $1.names $1.ali
