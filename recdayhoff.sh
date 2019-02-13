#!/bin/bash
'''
################################################
# Recodes an aminoacid alignment to the 6 Dayhoff categories. Phylip or fasta formats only.
# Author: Tauana Cunha | https://github.com/tauanajc/
# Date: February 08, 2018
# Usage: recdayhoff.sh input
################################################
'''

# Exit if no input file
if [ ! $# == 1 ]; then
    echo "Input is one alignment file (phylip or fasta). Usage: recdayhoff.sh inputFile"
    exit
fi

index=$(echo $RANDOM) # To avoid issues if using script for multiple files at the same time


# Detect if fasta or phylip and separate taxa names from alignment itself
if [ $(grep -c ">" $1) != 0 ]; then # Fasta
    grep ">" $1 > names${index}
    grep -v ">" $1 > alignment${index}
else # Phylip
# Separate phylip header:
    head -n 1 $1 >> header${index}
# Separate alignment:
    tail -n +2 $1 > alignment${index}
    sed -i 's:\s\+:\t:' alignment${index} # Converts spaces after taxa names into one tab
    cat alignment${index} | cut -d$'\t' -f1 > names${index} # Separates taxa names from actual alignment
    cat alignment${index} | cut -d$'\t' -f2 > temp${index} && mv temp${index} alignment${index}
# Count lenght of names to add the spaces necessary in each line to get a nicely aligned phylip again:
    spac=$(($(wc -L names${index} | cut -d ' ' -f 1)+2)) # Length of longest taxa name plus 2
    while read line; do
    add=$((${spac}-$(echo $line | wc -L | cut -d ' ' -f 1)))
    perl -E "print ' ' x $add" >> spaces${index}
    echo "" >> spaces${index}
    done < names${index}
    paste -d '' names${index} spaces${index} > temp${index} && mv temp${index} names${index}
    rm spaces${index}
fi


# Dayhoff recoding
sed -i 's:[AGPST]:0:g' alignment${index}
sed -i 's:[DENQ]:1:g' alignment${index}
sed -i 's:[HKR]:2:g' alignment${index}
sed -i 's:[WFY]:3:g' alignment${index}
sed -i 's:[MIVL]:4:g' alignment${index}
sed -i 's:[C]:5:g' alignment${index}


# Put header, names and alignment back in one file
if [ ! -f header${index} ]; then # Fasta
    paste -d '\n' names${index} alignment${index} > dayhoff
else # Phylip
    paste -d '' names${index} alignment${index} >> header${index} && mv header${index} dayhoff
fi

rm alignment${index} names${index} # Remove temporary files


