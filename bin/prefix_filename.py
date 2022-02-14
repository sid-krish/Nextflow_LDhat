#!/usr/bin/env python
import sys

fasta = sys.argv[1]
prepend_file = sys.argv[2]

if prepend_file == "none":
    file_name = fasta.rsplit(".",1)[0]
    print(f"{file_name}_", end = '')

else:
    print(prepend_file, end = '')
