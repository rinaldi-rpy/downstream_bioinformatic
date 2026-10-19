"""Parse final TreeMix likelihood explicitly; reject ambiguous fallback formats."""
import math
from pathlib import Path
import re

NUMBER = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?'

def read_likelihood(path):
    lines=Path(path).read_text().splitlines()
    labelled=[line for line in lines if re.search(r'likelihood|\bllik\b',line,re.I)]
    if labelled:
        line=labelled[-1]
        match=re.search(r'('+NUMBER+r')\s*$',line)
        if not match:raise ValueError('Invalid final likelihood: '+str(path))
        value=float(match.group(1))
    else:
        numeric=[line.strip() for line in lines if line.strip()]
        if len(numeric)!=1 or not re.fullmatch(NUMBER,numeric[0]):
            raise ValueError('Unrecognized likelihood format: '+str(path))
        value=float(numeric[0])
    if not math.isfinite(value):raise ValueError('Nonfinite likelihood: '+str(path))
    return value
