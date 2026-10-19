#!/usr/bin/env python3
"""Record resolved configuration, code hashes, and tool locations."""
import datetime
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys

source,out=map(Path,sys.argv[1:])
keys=['THREADS','MIN_QUAL','MIN_DP','MAX_DP','APPLY_DP','MAF','MISS','LD_WINDOW','LD_STEP',
      'LD_R2','N_PC','K_MIN','K_MAX','N_REPS','CV_FOLDS','SEED','FST_WINDOW','FST_STEP',
      'GROUP_COLUMN','BLOCK_BP','N_BOOT','MIN_CALLED','MAX_M','REPS','BLOCK_SIZE','ROOT']
versions={}
for tool,flag in [('bcftools','--version'),('vcftools','--version'),('plink','--version'),('Rscript','--version')]:
    executable=os.environ.get({'bcftools':'BCFTOOLS','vcftools':'VCFTOOLS','plink':'PLINK'}.get(tool,''),tool)
    path=shutil.which(executable)
    if path:
        result=subprocess.run([path,flag],capture_output=True,text=True)
        versions[tool]=(result.stdout+result.stderr).strip()
record={'utc_time':datetime.datetime.now(datetime.timezone.utc).isoformat(),'python':platform.python_version(),
        'parameters':{k:os.environ.get(k) for k in keys},'versions':versions,
        'executables':{k:shutil.which(os.environ.get(k,k.lower())) for k in ['PLINK','TREEMIX','ADMIXTURE','BCFTOOLS','VCFTOOLS']},
        'code_sha256':{str(p.relative_to(source)):hashlib.sha256(p.read_bytes()).hexdigest()
                       for folder in ['python','r','shell','config'] for p in (source/folder).rglob('*') if p.is_file() and '__pycache__' not in str(p)}}
for p in source.glob('*.sh'):record['code_sha256'][p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
(out/'provenance.json').write_text(json.dumps(record,indent=2)+'\n')
