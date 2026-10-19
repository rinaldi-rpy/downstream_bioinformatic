#!/usr/bin/env python3
"""Stream VCF counts and QUAL/INFO-DP/FORMAT-DP descriptive statistics."""
import argparse
import csv
import gzip
import math
from pathlib import Path

class RunningStats:
    def __init__(self):
        self.n=0; self.mean=0.; self.m2=0.; self.minimum=math.inf;self.maximum=-math.inf
    def add(self,value):
        if value in {'.',''}:return
        value=float(value)
        if not math.isfinite(value):raise ValueError('Nonfinite VCF value')
        self.n+=1;delta=value-self.mean;self.mean+=delta/self.n;self.m2+=delta*(value-self.mean)
        self.minimum=min(self.minimum,value);self.maximum=max(self.maximum,value)
    def values(self):
        return [self.n,self.mean if self.n else 'NA',self.minimum if self.n else 'NA',
                self.maximum if self.n else 'NA',math.sqrt(max(0,self.m2)/(self.n-1)) if self.n>1 else 'NA']


def summarize(path):
    opener=gzip.open if str(path).endswith('.gz') else open
    stats={key:RunningStats() for key in ['QUAL','INFO_DP','FORMAT_DP']}
    count,snps=0,0
    with opener(path,'rt') as handle:
        for line in handle:
            if line.startswith('#'):continue
            fields=line.rstrip().split('\t')
            if len(fields)<8:raise ValueError('Invalid VCF row')
            count+=1;snps+=int(len(fields[3])==1 and all(len(a)==1 and a in 'ACGT' for a in fields[4].split(',')))
            stats['QUAL'].add(fields[5])
            info=dict(item.split('=',1) for item in fields[7].split(';') if '=' in item)
            stats['INFO_DP'].add(info.get('DP','.'))
            if len(fields)>9 and 'DP' in fields[8].split(':'):
                index=fields[8].split(':').index('DP')
                for sample in fields[9:]:
                    values=sample.split(':'); stats['FORMAT_DP'].add(values[index] if len(values)>index else '.')
    return count,snps,stats


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw',required=True);parser.add_argument('--qc',required=True)
    parser.add_argument('--bim',required=True);parser.add_argument('--out',required=True)
    args=parser.parse_args();output=Path(args.out)
    with output.open('x',newline='') as handle:
        writer=csv.writer(handle);writer.writerow(['Dataset','Metric','N','Mean','Min','Max','SD'])
        for label,path in [('RAW',args.raw),('QC',args.qc)]:
            count,snps,stats=summarize(path)
            writer.writerow([label,'Variant_records',count,'NA','NA','NA','NA'])
            writer.writerow([label,'SNP_records',snps,'NA','NA','NA','NA'])
            for metric,values in stats.items():writer.writerow([label,metric]+values.values())
        with open(args.bim) as source:n=sum(1 for line in source if line.strip())
        writer.writerow(['LD','Retained_SNPs',n,'NA','NA','NA','NA'])

if __name__=='__main__':main()
