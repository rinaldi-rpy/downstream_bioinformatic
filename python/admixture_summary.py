#!/usr/bin/env python3
"""Keep replicate identity in CV summaries; never average unaligned Q matrices."""
import csv
import math
from pathlib import Path
import re
import shutil
import sys


def main():
    out=Path(sys.argv[1]); kmin,kmax,reps=map(int,sys.argv[2:])
    rows=[];best=[]
    for k in range(kmin,kmax+1):
        current=[]
        for rep in range(1,reps+1):
            directory=out/f'K{k}_rep{rep}'
            matches=re.findall(r'CV error\s*\(K\s*=\s*'+str(k)+r'\)\s*:\s*([0-9.eE+-]+)',(directory/'admixture.log').read_text())
            if len(matches)!=1 or not math.isfinite(float(matches[0])):
                raise ValueError('Expected one finite CV error: '+str(directory))
            row={'K':k,'Replicate':rep,'CV_error':float(matches[0]),'Directory':str(directory)}
            current.append(row);rows.append(row)
        chosen=min(current,key=lambda row:row['CV_error']);best.append(chosen)
        destination=out/'best';destination.mkdir(exist_ok=True)
        for suffix in ['Q','P']:
            shutil.copy2(Path(chosen['Directory'])/f'input.{k}.{suffix}',destination/f'K{k}.{suffix}')
    for name,data in [('cv_all_replicates.csv',rows),('cv_best_per_k.csv',best)]:
        with (out/name).open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=list(data[0]));writer.writeheader();writer.writerows(data)

if __name__=='__main__':main()
