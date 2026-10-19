import csv, os, re
from collections import Counter
from pathlib import Path

out = Path(os.environ['OUT_DIR'])
bfile = os.environ['BFILE']
from metadata import load_metadata, match_samples
from remap_bim import chromosome_map
metadata, labels = load_metadata(os.environ['META_PATH'], os.environ.get('GROUP_COLUMN','Group Label'))
mapping = match_samples(bfile + '.fam', metadata)
if not mapping:
    raise SystemExit('No samples found.')
counts = Counter(row[2] for row in mapping)
(out / 'population_order.txt').write_text(''.join(p + '\n' for p in sorted(counts)))
if len(counts) < 3:
    raise SystemExit('This pipeline requires at least three populations.')
if min(counts.values()) < 2:
    raise SystemExit('Each population must have at least two individuals for sample-size correction.')
root = os.environ['ROOT']
if root and root not in counts:
    raise SystemExit(f'ROOT not found among sanitized population labels: {root}')
with (out / 'populations.txt').open('w', newline='') as f:
    csv.writer(f, delimiter='\t').writerows(mapping)
for filename, header, rows in [
    ('population_sizes.csv', ['Population', 'N_individuals'], sorted(counts.items())),
    ('population_labels.csv', ['Population', 'Original_label'],
     [(p, labels[p]) for p in sorted(counts)])
]:
    with (out / filename).open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)

# Dataset-specific exclusions: NC_049568.1 = mitochondrion; NC_059181.1 = X.
autosomes = chromosome_map(os.environ['AUTOSOME_MAP'])
accepted = set(autosomes) | set(autosomes.values())
excluded = {'NC_049568.1', 'NC_059181.1', 'X', 'Y', 'MT', '25', '26', '27'}
snps, chromosomes, n_autosomal = set(), Counter(), 0
closed, previous_chr, previous_bp = set(), None, -1
with open(bfile + '.bim') as f, (out / 'autosomal_snps.txt').open('w') as target:
    for line in f:
        fields = line.split()
        if len(fields) != 6:
            raise SystemExit('Invalid BIM row.')
        chrom, snp, _, bp, a1, a2 = fields
        chromosomes[chrom] += 1
        if snp in snps or snp in {'0', '.'}:
            raise SystemExit(f'Duplicate or missing SNP ID: {snp}')
        snps.add(snp)
        if chrom in excluded:
            continue
        if chrom not in accepted:
            raise SystemExit(f'Unexpected chromosome: {chrom}; review autosomal selection.')
        if a1.upper() not in {'A', 'C', 'G', 'T'} or a2.upper() not in {'A', 'C', 'G', 'T'} or a1.upper() == a2.upper():
            continue
        position = int(bp)
        if position < 1:
            raise ValueError("SNP positions must be positive")
        if chrom != previous_chr:
            if chrom in closed:
                raise SystemExit('Chromosomes are not contiguous in BIM; sort the input first.')
            if previous_chr is not None:
                closed.add(previous_chr)
            previous_chr, previous_bp = chrom, -1
        if position < previous_bp:
            raise SystemExit('SNPs are not position-sorted within chromosomes.')
        previous_bp = position
        target.write(snp + '\n')
        n_autosomal += 1
if not n_autosomal:
    raise SystemExit('No autosomal biallelic SNPs retained.')
with (out / 'chromosome_counts.csv').open('w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Chromosome', 'N_variants', 'Excluded'])
    writer.writerows((c, n, c in excluded) for c, n in sorted(chromosomes.items()))
print(f'Matched samples: {len(mapping)}; populations: {len(counts)}')
print(f'Autosomal biallelic SNPs selected: {n_autosomal}')
for population, n in sorted(counts.items()):
    print(f'  {population}: {n}')
if not root:
    print('No explicit outgroup selected; migration direction requires rooting review.')
