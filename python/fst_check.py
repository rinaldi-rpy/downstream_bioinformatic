import os, sys
if sys.version_info < (3, 9):
    raise SystemExit('Python >=3.9 is required.')
try:
    import numpy as np
except ImportError:
    raise SystemExit('NumPy is required. Install it in the active environment.')
for key in ['BLOCK_BP', 'N_BOOT', 'CHUNK_SNPS', 'MIN_CALLED']:
    if int(os.environ[key]) < 1:
        raise SystemExit(key + ' must be a positive integer.')
if int(os.environ['MIN_CALLED']) < 2 or int(os.environ['N_BOOT']) < 20:
    raise SystemExit('MIN_CALLED must be >=2 and N_BOOT must be >=20.')
seed = int(os.environ['SEED'])
if not 0 <= seed < 2**32:
    raise SystemExit('SEED must be between 0 and 2^32-1.')
print('Python:', sys.version.split()[0], '; NumPy:', np.__version__)
