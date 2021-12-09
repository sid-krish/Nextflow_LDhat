#!/usr/bin/env python
import sys

from Bio import SeqIO

inputFasta = sys.argv[1]
prepend_filename = sys.argv[2]

records = list(SeqIO.parse(inputFasta, "fasta"))
# print(">entry_" + records[0].id)  # first record
# print(records[0].seq)

with open(f"to_merge_{prepend_filename}.fa", 'w') as fileOut:
    for i in records:
        fileOut.write(f">{prepend_filename}" + i.id + '\n')
        fileOut.write(str(i.seq) + '\n')
