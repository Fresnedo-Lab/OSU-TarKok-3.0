#!/usr/bin/env python3
"""Predict ncRNAs with Barrnap, tRNAscan-SE and Infernal; resolve overlapping calls."""
import argparse
from collections import defaultdict
import re
from urllib.parse import quote
from _common import Runner, execution_options, path, positive, fasta, write_fasta


def classify(name):
    name = name.lower()
    if 'trna' in name:
        return 'tRNA'
    if any(x in name for x in ('rrna', 'ssu', 'lsu')) or name == '5s' or name.startswith(('5s_', '5_8s')):
        return 'rRNA'
    if 'snor' in name:
        return 'snoRNA'
    if 'snrna' in name or re.match(r'u(1|2|4|5|6|11|12)\b', name):
        return 'snRNA'
    if 'mir' in name:
        return 'miRNA'
    if 'srp' in name:
        return 'SRP_RNA'
    if 'rnase' in name:
        return 'RNase'
    return 'other'


def windows(a):
    if a.overlap >= a.window:
        raise ValueError('--overlap must be smaller than --window')
    a.out.mkdir(parents=True, exist_ok=False)
    batch, size, number, total = [], 0, 0, 0
    def save():
        nonlocal number
        number += 1
        with (a.out / f'chunk_{number:04d}.fa').open('w') as out:
            for header, seq in batch:
                write_fasta(out, header, seq)
    for header, sequence in fasta(a.genome):
        name = header.split()[0]
        total += len(sequence)
        for start in range(0, len(sequence), a.window - a.overlap):
            end = min(start + a.window, len(sequence))
            seq = sequence[start:end]
            if batch and size + len(seq) > a.window:
                save()
                batch, size = [], 0
            batch.append((f'{name}:{start+1}-{end}', seq))
            size += len(seq)
            if end == len(sequence):
                break
    if batch:
        save()
    (a.out / 'search_space_mb.txt').write_text(f'{2 * total / 1e6:.8f}\n')
    print(f'{number} chunks; cmscan -Z {2 * total / 1e6:.8f}')


def combine(a):
    rows, seen = [], set()
    for file in a.inputs:
        with file.open() as src:
            for line in src:
                if line.startswith('#') or not line.strip():
                    continue
                f = line.split()
                if len(f) < 18:
                    raise ValueError(f'Invalid cmscan fmt 2 row in {file}')
                name, acc, query, clan = f[1], f[2], f[3], f[5]
                match = re.fullmatch(r'(.+):(\d+)-(\d+)', query)
                seqid, offset = (match[1], int(match[2]) - 1) if match else (query, 0)
                start, end = sorted((int(f[9]), int(f[10])))
                start, end = start + offset, end + offset
                key = (seqid, start, end, f[11], name)
                if key in seen:
                    continue
                seen.add(key)
                attrs = f'Name={quote(name)};rfam_acc={acc};clan={"" if clan == "-" else clan};evalue={f[17]}'
                rows.append([seqid, 'Infernal', 'ncRNA', start, end, f[16], f[11], '.', attrs])
    with a.out.open('w') as out:
        out.write('##gff-version 3\n')
        for row in sorted(rows, key=lambda r: (r[0], r[3], '\t'.join(map(str, r)))):
            out.write('\t'.join(map(str, row)) + '\n')


def merge(a):
    features = []
    for source, file in [('tRNAscan-SE', a.trnascan), ('barrnap', a.barrnap), ('Infernal', a.infernal)]:
        with file.open() as src:
            for line in src:
                if line.startswith('#') or not line.strip():
                    continue
                f = line.rstrip('\n').split('\t')
                if len(f) != 9:
                    raise ValueError(f'Invalid GFF3 row in {file}')
                f[1] = source
                f[3], f[4] = sorted((int(f[3]), int(f[4])))
                score = float(f[5]) if f[5] != '.' else 0.0
                name = re.search(r'Name=([^;]+)', f[8])
                kind = 'tRNA' if source == 'tRNAscan-SE' else 'rRNA' if source == 'barrnap' else classify(name[1] if name else '')
                features.append((3 if source != 'Infernal' else 2, score, f, kind))
    features.sort(key=lambda r: (-r[0], -r[1]))
    intervals, kept = defaultdict(list), []
    for _, _, f, kind in features:
        key = (f[0], f[6])
        if any(f[3] <= end and start <= f[4] for start, end in intervals[key]):
            continue
        intervals[key].append((f[3], f[4]))
        kept.append((f, kind))
    with a.out.open('w') as out:
        out.write('##gff-version 3\n')
        for n, (f, kind) in enumerate(sorted(kept, key=lambda r: (r[0][0], r[0][3])), 1):
            if not re.search(r'(?:^|;)ID=', f[8]):
                f[8] = f'ID=ncRNA{n:06d};' + f[8]
            f[8] += f';ncRNA_class={kind};predictor={f[1]}'
            f[7] = '.'
            out.write('\t'.join(map(str, f)) + '\n')
    print(f'{len(kept)} calls; {len(features) - len(kept)} overlapping calls removed')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    for name in ('barrnap', 'trnascan', 'cmscan'):
        q = sub.add_parser(name)
        execution_options(q)
        q.add_argument('--genome', type=path, required=True, help='FASTA; cmscan also accepts a window chunk')
        if name == 'trnascan':
            q.add_argument('--config', type=path)
        if name == 'cmscan':
            q.add_argument('--rfam', type=path, required=True, help='cmpress-indexed Rfam.cm')
            q.add_argument('--clans', type=path, required=True)
            q.add_argument('--search-space-mb', type=float, required=True, help='twice the FULL genome size in Mb')
    q = sub.add_parser('windows')
    q.add_argument('--genome', type=path, required=True)
    q.add_argument('--out', type=path, required=True)
    q.add_argument('--window', type=positive, default=5000000)
    q.add_argument('--overlap', type=positive, default=20000)
    q = sub.add_parser('combine', help='remap cmscan fmt 2 windows to genome coordinates')
    q.add_argument('--inputs', type=path, nargs='+', required=True)
    q.add_argument('--out', type=path, required=True)
    q = sub.add_parser('merge')
    for flag in ('trnascan', 'barrnap', 'infernal', 'out'):
        q.add_argument('--' + flag, type=path, required=True)
    a = p.parse_args()
    if a.command in ('windows', 'combine', 'merge'):
        globals()[a.command](a)
        return
    r = Runner(a)
    if a.command == 'barrnap':
        r.run(['barrnap', '--kingdom', 'euk', '--threads', a.threads, a.genome], stdout=a.out/'rrna.gff3')
    elif a.command == 'trnascan':
        cmd = ['tRNAscan-SE', '-E', '--thread', a.threads, '-o', a.out/'trna.out',
               '-f', a.out/'trna.struct', '-m', a.out/'trna.stats', '-b', a.out/'trna.bed',
               '--gff', a.out/'trna.gff3', '-a', a.out/'trna.fa', '-l', a.out/'trna.log']
        if a.config:
            cmd += ['-c', a.config]
        r.run(cmd + [a.genome])
    else:
        if a.search_space_mb <= 0:
            p.error('--search-space-mb must be positive')
        r.run(['cmscan', '-Z', a.search_space_mb, '--cut_ga', '--rfam', '--nohmmonly', '--fmt', '2',
               '--clanin', a.clans, '--oclan', '--oskip', '--cpu', a.threads,
               '--tblout', a.out/'hits.tblout', '-o', a.out/'cmscan.txt', a.rfam, a.genome])


if __name__ == '__main__':
    main()
