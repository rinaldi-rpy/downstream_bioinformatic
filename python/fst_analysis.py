import csv
import itertools
import math
import os
import re
from array import array
from pathlib import Path
import numpy as np

from metadata import load_metadata, match_samples
from remap_bim import chromosome_map

def accepted_autosomes():
    mapping = chromosome_map(os.environ['AUTOSOME_MAP'])
    return set(mapping) | set(mapping.values())

def save_csv(path, header, rows):
    with Path(path).open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)

def read_populations(prefix, meta_path, out):
    metadata, _ = load_metadata(meta_path, os.environ.get('GROUP_COLUMN','Group Label'))
    matched = match_samples(prefix + '.fam', metadata)
    individuals = [(iid,population) for fid,iid,population in matched]
    groups = sorted({g for iid, g in individuals})
    if len(groups) < 4:
        raise ValueError('At least four populations are required for informative unrooted NJ splits.')
    subpops = [np.array([i for i, item in enumerate(individuals) if item[1] == g]) for g in groups]
    if min(len(x) for x in subpops) < int(os.environ['MIN_CALLED']):
        raise ValueError('Population sample counts are smaller than MIN_CALLED.')
    save_csv(out / 'sample_groups.csv', ['ID', 'Population'], individuals)
    save_csv(out / 'population_sizes.csv', ['Population', 'N_individuals'],
             zip(groups, [len(x) for x in subpops]))
    print('Samples: {}; populations: {}; pairs: {}'.format(len(individuals), len(groups), len(groups)*(len(groups)-1)//2), flush=True)
    return groups, subpops, len(individuals)

def read_blocks(prefix, block_bp, out):
    autosomes = accepted_autosomes()
    ids, blocks, lookup = array('i'), [], {}
    previous_chr, previous_pos, closed = None, -1, set()
    excluded = {}
    with open(prefix + '.bim') as handle:
        for line in handle:
            fields = line.split()
            if len(fields) != 6:
                raise ValueError('Invalid BIM row.')
            chrom, snp, cm, pos, a1, a2 = fields
            if chrom not in autosomes or a1.upper() not in 'ACGT' or a2.upper() not in 'ACGT' or len(a1) != 1 or len(a2) != 1 or a1.upper() == a2.upper():
                ids.append(-1)
                excluded[chrom] = excluded.get(chrom, 0) + 1
                continue
            pos = int(pos)
            if pos < 1:
                raise ValueError('Variant positions must be positive.')
            if chrom != previous_chr:
                if chrom in closed:
                    raise ValueError('BIM chromosomes are not contiguous; sort the genotypes first.')
                if previous_chr is not None:
                    closed.add(previous_chr)
                previous_chr, previous_pos = chrom, -1
            if pos < previous_pos:
                raise ValueError('BIM positions are not sorted within chromosomes.')
            previous_pos = pos
            key = (chrom, (pos-1)//block_bp)
            if key not in lookup:
                lookup[key] = len(blocks)
                blocks.append(key)
            ids.append(lookup[key])
    if len(blocks) < 2:
        raise ValueError('At least two SNP-bearing genomic blocks are required.')
    save_csv(out / 'blocks.csv', ['Block_ID', 'CHR', 'START', 'END'],
             ((i, c, b*block_bp+1, (b+1)*block_bp) for i, (c,b) in enumerate(blocks)))
    save_csv(out / 'excluded_variants.csv', ['CHR', 'N_excluded'], sorted(excluded.items()))
    return np.asarray(ids, dtype=np.int32), blocks

def wc_components(n, allele_count, heterozygotes, pair_i, pair_j, min_called):
    n1, n2 = n[:, pair_i], n[:, pair_j]
    nt = n1 + n2
    nb = nt / 2.0
    with np.errstate(divide='ignore', invalid='ignore'):
        p1 = allele_count[:, pair_i] / (2*n1)
        p2 = allele_count[:, pair_j] / (2*n2)
        pb = (allele_count[:, pair_i] + allele_count[:, pair_j]) / (2*nt)
        s2 = (n1*(p1-pb)**2 + n2*(p2-pb)**2) / nb
        hb = (heterozygotes[:, pair_i] + heterozygotes[:, pair_j]) / nt
        nc = nt - (n1*n1 + n2*n2)/nt
        a = nb/nc * (s2 - (pb*(1-pb) - s2/2 - hb/4)/(nb-1))
        b = nb/(nb-1) * (pb*(1-pb) - s2/2 - (2*nb-1)*hb/(4*nb))
        den = a + b + hb/2
    valid = (n1 >= min_called) & (n2 >= min_called) & np.isfinite(a) & np.isfinite(den) & (den > 0)
    return np.where(valid, a, 0), np.where(valid, den, 0), valid

def aggregate_genotypes(prefix, nsamples, subpops, ids, pairs, nblocks, chunk, minimum):
    width = (nsamples+3)//4
    path = Path(prefix + '.bed')
    with path.open('rb') as handle:
        if handle.read(3) != bytes([0x6c, 0x1b, 0x01]):
            raise ValueError('Expected a SNP-major PLINK BED file.')
    if path.stat().st_size != 3 + len(ids)*width:
        raise ValueError('BED size does not match BIM and FAM.')
    bed = np.memmap(str(path), mode='r', dtype=np.uint8, offset=3, shape=(len(ids), width))
    total_a = np.zeros((nblocks, len(pairs)), dtype=float)
    total_d = np.zeros_like(total_a)
    counts = np.zeros(total_a.shape, dtype=np.int64)
    pi, pj = np.array(pairs, dtype=int).T
    sample_index = np.arange(nsamples)
    for start in range(0, len(ids), chunk):
        stop = min(start+chunk, len(ids))
        selected = ids[start:stop] >= 0
        if not selected.any():
            continue
        # PLINK BED: 00=A1/A1, 01=missing, 10=A1/A2, 11=A2/A2.
        codes = (bed[start:stop, sample_index//4] >> (2*(sample_index%4))) & 3
        codes = codes[selected]
        n = np.column_stack([np.sum(codes[:, s] != 1, axis=1) for s in subpops]).astype(float)
        ac = np.column_stack([np.sum(codes[:, s] == 0, axis=1)*2 + np.sum(codes[:, s] == 2, axis=1) for s in subpops]).astype(float)
        het = np.column_stack([np.sum(codes[:, s] == 2, axis=1) for s in subpops]).astype(float)
        a, d, valid = wc_components(n, ac, het, pi, pj, minimum)
        bid = ids[start:stop][selected]
        starts = np.r_[0, np.flatnonzero(bid[1:] != bid[:-1])+1]
        keys = bid[starts]
        total_a[keys] += np.add.reduceat(a, starts, axis=0)
        total_d[keys] += np.add.reduceat(d, starts, axis=0)
        counts[keys] += np.add.reduceat(valid.astype(np.int64), starts, axis=0)
        if start == 0 or stop == len(ids) or start//1000000 != (stop-1)//1000000:
            print('Genotypes processed: {:,}/{:,}'.format(stop, len(ids)), flush=True)
    return total_a, total_d, counts

def bootstrap_weights(rng, strata, nblocks):
    sampled = np.concatenate([indices[rng.randint(0, len(indices), size=len(indices))] for indices in strata])
    return np.bincount(sampled, minlength=nblocks)

def joint_bootstrap(a, d, blocks, nboot, seed):
    strata = [np.array([i for i, key in enumerate(blocks) if key[0] == chrom]) for chrom in sorted({c for c,b in blocks})]
    rng = np.random.RandomState(seed)
    boot = np.full((nboot, a.shape[1]), np.nan)
    for rep in range(nboot):
        weights = bootstrap_weights(rng, strata, len(blocks))
        denominator = weights.dot(d)
        if np.all(denominator > 0):
            estimate = weights.dot(a) / denominator
            if np.isfinite(estimate).all():
                boot[rep] = estimate
        if (rep+1) % 100 == 0:
            print('Bootstrap replicates: {}/{}'.format(rep+1, nboot), flush=True)
    valid = np.isfinite(boot).all(axis=1)
    if valid.sum() < max(20, int(math.ceil(0.9*nboot))):
        raise ValueError('Too many bootstrap replicates lack FST estimates for all pairs. Review missingness and population sizes.')
    return boot[valid], np.flatnonzero(valid)+1, int((~valid).sum())

def compare_old_logs(groups, pairs, estimate, fst_dir, out):
    rows = []
    pattern = r'weighted Fst estimate:\s*([-+0-9.eE]+)'
    for index, (i,j) in enumerate(pairs):
        path = Path(fst_dir) / (groups[i] + '_vs_' + groups[j] + '.log')
        if not path.is_file():
            path = Path(fst_dir) / (groups[j] + '_vs_' + groups[i] + '.log')
        if path.is_file():
            matches = re.findall(pattern, path.read_text(errors='replace'))
            if matches:
                old = float(matches[-1])
                rows.append((groups[i], groups[j], old, estimate[index], estimate[index]-old))
    save_csv(out / 'comparison_with_original_logs.csv',
             ['Population1', 'Population2', 'Original_log_FST', 'Autosomal_recomputed_FST', 'Difference'], rows)

def main():
    out = Path(os.environ['OUT_DIR'])
    import shutil
    shutil.copy2(os.environ['AUTOSOME_MAP'], out / 'autosomes.tsv')
    prefix = os.environ['ANALYSIS_BFILE']
    groups, subpops, nsamples = read_populations(prefix, os.environ['META_PATH'], out)
    block_bp = int(os.environ['BLOCK_BP'])
    ids, blocks = read_blocks(prefix, block_bp, out)
    pairs = list(itertools.combinations(range(len(groups)), 2))
    a, d, counts = aggregate_genotypes(prefix, nsamples, subpops, ids, pairs, len(blocks),
                                      int(os.environ['CHUNK_SNPS']), int(os.environ['MIN_CALLED']))
    denominator = d.sum(axis=0)
    if np.any(denominator <= 0):
        raise ValueError('Some population pairs have no informative sites.')
    estimate = a.sum(axis=0)/denominator
    np.savez_compressed(str(out / 'block_components.npz'), numerator=a, denominator=d, n_sites=counts)
    pair_ids = ['Pair_{:04d}'.format(i+1) for i in range(len(pairs))]
    save_csv(out / 'pair_definitions.csv', ['Pair_ID', 'Population1', 'Population2'],
             ((key, groups[i], groups[j]) for key, (i,j) in zip(pair_ids, pairs)))
    boot, replicate_ids, invalid = joint_bootstrap(a, d, blocks, int(os.environ['N_BOOT']), int(os.environ['SEED']))
    lo, hi = np.percentile(boot, [2.5, 97.5], axis=0)
    save_csv(out / 'bootstrap_fst.csv', ['Replicate'] + pair_ids,
             ((int(rep),)+tuple(values) for rep, values in zip(replicate_ids, boot)))
    save_csv(out / 'pairwise_fst_summary.csv',
             ['Pair_ID', 'Population1', 'Population2', 'FST', 'CI95_Lower', 'CI95_Upper', 'Bootstrap_SE', 'N_sites', 'N_informative_blocks', 'N_bootstrap'],
             ((pair_ids[k], groups[i], groups[j], estimate[k], lo[k], hi[k], np.std(boot[:, k], ddof=1),
               int(counts[:, k].sum()), int((d[:, k]>0).sum()), len(boot)) for k,(i,j) in enumerate(pairs)))
    for name, values in [('fst_matrix.csv', estimate), ('ci95_lower_matrix.csv', lo), ('ci95_upper_matrix.csv', hi)]:
        matrix = np.zeros((len(groups), len(groups)))
        for value, (i,j) in zip(values, pairs):
            matrix[i,j] = matrix[j,i] = value
        save_csv(out / name, ['Population'] + groups, ((g,)+tuple(row) for g,row in zip(groups, matrix)))
    compare_old_logs(groups, pairs, estimate, os.environ['FST_DIR'], out)
    save_csv(out / 'bootstrap_qc.csv', ['Metric', 'Value'],
             [('Requested_replicates', os.environ['N_BOOT']), ('Valid_joint_replicates', len(boot)),
              ('Invalid_joint_replicates', invalid), ('N_blocks', len(blocks)), ('Block_bp', block_bp),
              ('Negative_observed_FST_pairs', int((estimate < 0).sum()))])
    methods = '''Estimator: pairwise Weir-Cockerham theta for diploid biallelic SNPs.
Genome-wide FST = sum(a) / sum(a+b+c), not mean(window FST).
Autosomes: accessions and numeric aliases from the archived autosomes.tsv map; X and mitochondrion excluded.
Partial genotype calls are treated as missing during VCF conversion.
Sites require MIN_CALLED complete diploid individuals in each member of a pair.
Blocks: fixed physical bins with at least one selected SNP; bins never cross chromosomes.
Bootstrap: sample SNP-bearing blocks with replacement within each chromosome,
preserving the original number of blocks per chromosome. The same draws are
used for all pairs. A replicate is retained only when all pair estimates exist.
CI: percentile 95%; negative FST estimates are retained for estimates and CIs.
NJ uses raw FST as a descriptive dissimilarity, with negative distances set to 0.
NJ is unrooted and does not model migration; it is not a sequence phylogeny.
Bootstrap support describes stability of unrooted splits under block resampling.
Genomic resampling conditions on sampled individuals and the QC SNP ascertainment.
Block length is an initial setting; examine sensitivity to LD/block size.
Original windowed FST files and logs are read only, not used to reconstruct components.
Differences from original logs can reflect autosomal/biallelic/site-call filtering.
Estimator reference: https://doi.org/10.1111/j.1558-5646.1984.tb05657.x
Component reference: https://scikit-allel.readthedocs.io/en/stable/stats/fst.html
'''
    (out / 'METHODS.txt').write_text(methods)
    print('Analysis complete: {} pairs; {} valid joint bootstrap replicates; {} invalid.'.format(len(pairs), len(boot), invalid), flush=True)

if __name__ == '__main__':
    main()
