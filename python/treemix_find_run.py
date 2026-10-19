import os, sys
from pathlib import Path
if sys.version_info < (3, 6):
    raise SystemExit('Python 3.6 or newer is required.')
requested = os.environ['RUN_DIR']
if requested:
    folder = Path(requested).expanduser().resolve()
    if not (folder / 'runs').is_dir():
        raise SystemExit('Expected a run directory containing a runs/ subdirectory: ' + str(folder))
else:
    search = Path(os.environ['SEARCH_DIR'])
    candidates = sorted({p.parent.parent.resolve() for p in search.glob('run_*/runs/*.llik')})
    if (search / 'runs').is_dir():
        candidates.append(search.resolve())
    candidates = sorted(set(candidates))
    if len(candidates) != 1:
        raise SystemExit('Specify one existing run directory as the first argument. Candidates:\n' +
                         '\n'.join(str(p) for p in candidates))
    folder = candidates[0]
print(folder)
