import os, re
from pathlib import Path
folder = Path(os.environ['RUN_DIR'])
settings = {}
path = folder / 'settings.txt'
if path.exists():
    for line in path.read_text().splitlines():
        if '=' in line:
            key, value = line.split('=', 1)
            settings[key] = value
pairs = []
for path in (folder / 'runs').glob('*.llik'):
    match = re.fullmatch(r'treemix_m([0-9]+)_rep([0-9]+)\.llik', path.name)
    if match:
        pairs.append(tuple(map(int, match.groups())))
if not pairs:
    raise SystemExit('No treemix_m*_rep*.llik files found.')
max_m = int(settings.get('MAX_M', max(m for m, r in pairs)))
reps = int(settings.get('REPS', max(r for m, r in pairs)))
if max_m < 0 or reps < 1:
    raise SystemExit('Invalid saved model settings.')
expected = {(m, r) for m in range(max_m + 1) for r in range(1, reps + 1)}
missing = sorted(expected - set(pairs))
if missing:
    raise SystemExit('Missing expected likelihood files (m, replicate): ' + str(missing))
print(max_m, reps)
