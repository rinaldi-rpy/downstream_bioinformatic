import gzip, os, shutil
from pathlib import Path
source, out = Path(os.environ['RUN_DIR']), Path(os.environ['OUT_DIR'])
order = source / 'population_order.txt'
if order.is_file():
    shutil.copy2(str(order), str(out / order.name))
elif (source / 'treemix_input.gz').is_file():
    with gzip.open(str(source / 'treemix_input.gz'), 'rt') as f:
        populations = f.readline().split()
    (out / 'population_order.txt').write_text('\n'.join(populations) + '\n')
else:
    print('Population order unavailable; residual plots will be skipped.')
