#!/usr/bin/env python3
"""Validated sample-to-population mapping; preserve IDs as strings."""
import argparse
import csv
from pathlib import Path
import re


def load_metadata(path, group='Group Label'):
    path = Path(path)
    if path.suffix.lower() in {'.csv', '.tsv'}:
        with path.open(encoding='utf-8-sig', newline='') as handle:
            reader = csv.DictReader(handle, delimiter='\t' if path.suffix.lower() == '.tsv' else ',')
            if not {'ID', group}.issubset(reader.fieldnames or []):
                raise ValueError('Metadata must contain ID and ' + group)
            rows = list(reader)
    elif path.suffix.lower() in {'.xlsx', '.xls'}:
        import pandas as pd
        frame = pd.read_excel(path, dtype=str, keep_default_na=False)
        if not {'ID', group}.issubset(frame.columns):
            raise ValueError('Metadata must contain ID and ' + group)
        rows = frame.to_dict('records')
    else:
        raise ValueError('Use CSV, TSV, XLSX or XLS metadata.')
    mapping, labels = {}, {}
    for row in rows:
        iid, original = row['ID'], row[group]
        if not iid or iid != iid.strip() or re.search(r'\s', iid) or not original or not original.strip():
            raise ValueError('Empty ID/group or whitespace in sample ID.')
        if iid in mapping:
            raise ValueError('Duplicate metadata ID: ' + iid)
        population = re.sub(r'[^A-Za-z0-9_.-]', '_', original.strip())
        if not re.search(r'[A-Za-z0-9]', population):
            raise ValueError('Invalid population label: ' + original)
        if population in labels and labels[population] != original:
            raise ValueError('Population labels collide after sanitization: ' + original)
        mapping[iid], labels[population] = population, original
    if not mapping:
        raise ValueError('No metadata records.')
    return mapping, labels


def match_samples(path, mapping):
    rows, seen = [], set()
    is_fam = Path(path).suffix == '.fam'
    for line in Path(path).read_text().splitlines():
        fields = line.split()
        if (is_fam and len(fields) != 6) or (not is_fam and len(fields) != 1):
            raise ValueError('Invalid sample/FAM row.')
        fid, iid = fields[:2] if is_fam else (fields[0], fields[0])
        if iid in seen or iid not in mapping:
            raise ValueError('Duplicate or unmatched sample ID: ' + iid)
        seen.add(iid)
        rows.append((fid, iid, mapping[iid]))
    if not rows:
        raise ValueError('No samples in genotype input.')
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['populations', 'cluster', 'popfile', 'summary'])
    parser.add_argument('--meta', required=True)
    parser.add_argument('--group', default='Group Label')
    parser.add_argument('--samples', help='FAM or one VCF sample ID per line; required for population outputs')
    parser.add_argument('--out')
    args = parser.parse_args()
    mapping, labels = load_metadata(args.meta, args.group)
    if args.command == 'summary':
        print('Samples:', len(mapping), '; populations:', len(labels))
        return
    if not args.samples or not args.out:
        parser.error('--samples and --out are required for population outputs')
    rows = match_samples(args.samples, mapping)
    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    if args.command == 'populations':
        output.mkdir(exist_ok=True)
        if list(output.glob('*.txt')):
            raise ValueError('Population output must be fresh (avoid stale groups).')
        for population in sorted({r[2] for r in rows}):
            (output / (population + '.txt')).write_text(''.join(r[1]+'\n' for r in rows if r[2] == population))
        label_path = output/'population_labels.csv'
    else:
        with output.open('w') as handle:
            for fid, iid, population in rows:
                handle.write(('\t'.join((fid, iid, population)) if args.command == 'cluster' else '\t'.join((iid, population)))+'\n')
        label_path = output.parent/'population_labels.csv'
    with label_path.open('w', newline='') as handle:
        writer=csv.writer(handle)
        writer.writerow(['Population', 'Original_label'])
        writer.writerows((p,labels[p]) for p in sorted({r[2] for r in rows}))

if __name__ == '__main__':
    main()
