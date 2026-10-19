#!/usr/bin/env python3
"""Assign numeric chromosome aliases in an ADMIXTURE-only copy."""
from pathlib import Path
import sys

source=Path(sys.argv[1]); mapping={}; rows=[]
for line in source.read_text().splitlines():
    fields=line.split()
    if len(fields)!=6:raise ValueError('Invalid BIM row')
    chrom=fields[0]
    if chrom not in mapping:mapping[chrom]=str(len(mapping)+1)
    fields[0]=mapping[chrom];rows.append('\t'.join(fields))
Path(sys.argv[2]).write_text('Original\tNumeric\n'+''.join(k+'\t'+v+'\n' for k,v in mapping.items()))
source.write_text('\n'.join(rows)+'\n')
