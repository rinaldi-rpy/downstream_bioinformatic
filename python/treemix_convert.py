import csv, gzip, itertools, os
from decimal import Decimal
from pathlib import Path

out = Path(os.environ['OUT_DIR'])
with (out / 'population_sizes.csv').open(newline='') as f:
    populations = [r['Population'] for r in csv.DictReader(f)]
pop_set = set(populations)
stats = dict(total=0, retained=0, missing_population=0, monomorphic=0)

def integer_count(value):
    n = Decimal(value)
    if not n.is_finite() or n < 0 or n != n.to_integral_value():
        raise ValueError(f'Invalid allele count: {value}')
    return int(n)

def rows(handle):
    names = next(handle).split()
    required = {'CHR', 'SNP', 'CLST', 'A1', 'A2', 'MAC', 'NCHROBS'}
    if not required.issubset(names):
        raise SystemExit(f'Unexpected frequency header: {names}')
    for line in handle:
        values = line.split()
        if not values:
            continue
        if len(values) != len(names):
            raise SystemExit('Malformed frequency row.')
        yield dict(zip(names, values))

seen = set()
with (out / 'autosomal.frq.strat').open() as source, \
     gzip.open(out / 'treemix_input.gz', 'wt') as target, \
     (out / 'retained_snps.tsv').open('w') as snp_out:
    target.write(' '.join(populations) + '\n')
    snp_out.write('CHR\tSNP\n')
    for key, group in itertools.groupby(rows(source), key=lambda r: (r['CHR'], r['SNP'])):
        if key in seen:
            raise SystemExit('Frequency rows are not grouped by SNP.')
        seen.add(key)
        stats['total'] += 1
        counts, alleles = {}, None
        for row in group:
            p = row['CLST']
            if p not in pop_set or p in counts:
                raise SystemExit(f'Unexpected or duplicate population at SNP {key}: {p}')
            pair = (row['A1'], row['A2'])
            if alleles is not None and pair != alleles:
                raise SystemExit(f'Inconsistent allele orientation at SNP {key}')
            alleles = pair
            n = integer_count(row['NCHROBS'])
            mac = integer_count(row['MAC']) if n else 0
            if mac > n:
                raise SystemExit(f'Allele count exceeds total at SNP {key}')
            counts[p] = (mac, n - mac)
        if set(counts) != pop_set or any(sum(c) == 0 for c in counts.values()):
            stats['missing_population'] += 1
            continue
        totals = [sum(c[i] for c in counts.values()) for i in (0, 1)]
        if min(totals) == 0:
            stats['monomorphic'] += 1
            continue
        target.write(' '.join(f'{counts[p][0]},{counts[p][1]}' for p in populations) + '\n')
        snp_out.write('\t'.join(key) + '\n')
        stats['retained'] += 1
with (out / 'input_qc.csv').open('w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Metric', 'N_SNPs'])
    writer.writerows(stats.items())
print('Input QC:', stats)
if stats['retained'] < 2 * int(os.environ['BLOCK_SIZE']):
    raise SystemExit('Too few retained SNPs for the selected block size.')
