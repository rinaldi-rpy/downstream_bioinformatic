import csv, math, os, re, shutil, statistics
from pathlib import Path
from likelihood import read_likelihood

out = Path(os.environ['OUT_DIR'])
source_dir = Path(os.environ.get('RUN_DIR') or os.environ['OUT_DIR'])
number = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?'
runs, best, summary = [], [], []
for m in range(int(os.environ['MAX_M']) + 1):
    current = []
    for rep in range(1, int(os.environ['REPS']) + 1):
        prefix = source_dir / 'runs' / f'treemix_m{m}_rep{rep}'
        tree = Path(str(prefix) + '.treeout.gz')
        if not tree.is_file() or tree.stat().st_size == 0:
            raise SystemExit('Missing TreeMix output: ' + str(tree))
        final_likelihood = read_likelihood(str(prefix) + '.llik')
        row = dict(M=m, Replicate=rep, LogLikelihood=final_likelihood, Prefix=str(prefix))
        current.append(row)
        runs.append(row)
    chosen = max(current, key=lambda r: r['LogLikelihood'])
    best_prefix = out / 'best' / f'treemix_m{m}'
    for source in (source_dir / 'runs').glob(Path(chosen['Prefix']).name + '.*'):
        suffix = source.name[len(Path(chosen['Prefix']).name):]
        shutil.copy2(source, str(best_prefix) + suffix)
    best.append(dict(chosen, BestPrefix=str(best_prefix)))
    likelihoods = [r['LogLikelihood'] for r in current]
    summary.append(dict(M=m, N=len(current), Mean_LogLikelihood=statistics.mean(likelihoods),
                        SD_LogLikelihood=statistics.stdev(likelihoods) if len(current) > 1 else 0,
                        Best_LogLikelihood=chosen['LogLikelihood']))
for filename, data in [('all_runs.csv', runs), ('best_runs.csv', best), ('model_summary.csv', summary)]:
    with (out / filename).open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(data[0]))
        writer.writeheader()
        writer.writerows(data)
print('Best replicate per model saved. Migration count is not automatically selected.')
