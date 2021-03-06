#!/usr/bin/env python
import sys
from Bio import SeqIO

inputFile = sys.argv[1]

records = list(SeqIO.parse(inputFile, "fasta"))

with open('LDhat_reformated.fa', 'w') as fileOut:
    # Add required first line for LDhat

    # From Manual: Full sequence data should be aligned and in a modified FASTA format,
    # with the first line detailing the number of sequences/genotypes,
    # the number of sites in the alignment and a flag (1 or 2) that
    # details whether the data is haplotype/phased (1) or genotype/unphased (2).

    numSequences = len(records)
    numSitesInAln = len(str(records[0].seq))
    phaseType = 1

    fileOut.write(f"{numSequences} {numSitesInAln} {phaseType} \n")
    for i in records:
        fileOut.write(">genome_" + i.id + '\n')
        fileOut.write(str(i.seq) + '\n')

    sys.stdout.write(str(numSequences)) # for next steps

