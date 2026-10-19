#!/usr/bin/env python3
"""Map assembly accessions to explicit numeric autosomal codes, preserving BED order."""
from pathlib import Path
import sys


def chromosome_map(path):
    rows=[line.split() for line in Path(path).read_text().splitlines() if line.strip() and not line.startswith('#')]
    if not rows or any(len(row)!=2 for row in rows):
        raise ValueError('Invalid chromosome mapping')
    mapping=dict(rows)
    numbers=[int(row[1]) for row in rows]
    if len(mapping)!=len(rows) or len(set(numbers))!=len(rows) or any(n<1 or n>24 for n in numbers):
        raise ValueError('Chromosome map must have unique accessions and numeric codes 1..24')
    return mapping


def main():
    source, map_path, output = sys.argv[1:]
    mapping=chromosome_map(map_path)
    seen,closed=set(),set()
    previous,pos=None,0
    with Path(output).open('w') as target:
        for line in Path(source).read_text().splitlines():
            fields=line.split()
            if len(fields)!=6 or fields[0] not in mapping or fields[1] in seen or fields[1] in {'.','0'}:
                raise ValueError('Invalid, unexpected or duplicate BIM variant: '+line)
            chrom,position=fields[0],int(fields[3])
            if chrom!=previous:
                if chrom in closed:raise ValueError('Chromosomes must be contiguous')
                if previous is not None:closed.add(previous)
                previous,pos=chrom,0
            if position<1 or position<pos:raise ValueError('BIM positions must be sorted and positive')
            pos=position; seen.add(fields[1]); fields[0]=mapping[chrom]
            target.write('\t'.join(fields)+'\n')
    if not seen:raise ValueError('No autosomal SNPs found; check reference assembly and chromosome map')

if __name__=='__main__':main()
