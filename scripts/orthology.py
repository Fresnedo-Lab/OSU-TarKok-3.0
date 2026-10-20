#!/usr/bin/env python3
"""Prepare representative proteomes, run OrthoFinder and format CAFE5 counts."""
import argparse
import csv
import re
from _common import Runner, execution_options, path, positive, fasta, text_open, attributes, write_fasta


def representatives(a):
    mapping = {}
    if a.gff3:
        with text_open(a.gff3) as src:
            for line in src:
                if line.startswith('#') or '\tCDS\t' not in line:
                    continue
                fields = attributes(line.rstrip().split('\t')[8])
                pid, gene = fields.get('protein_id'), fields.get('gene') or fields.get('locus_tag')
                if pid and gene:
                    mapping.setdefault(pid, gene)
    best = {}
    for header, sequence in fasta(a.proteins):
        pid = header.split()[0]
        if a.format == 'gwh':
            match = re.search(r'Gene=(\S+)', header)
            gene = match[1] if match else 'UNMAPPED::' + pid
        else:
            gene = mapping.get(pid, 'UNMAPPED::' + pid)
        if gene not in best or len(sequence) > len(best[gene][1]):
            best[gene] = (header, sequence)
    a.out.mkdir(parents=True, exist_ok=False)
    with (a.out/'representatives.fa').open('w') as out, (a.out/'protein_gene.tsv').open('w') as table:
        table.write('protein_id\tgene_id\n')
        for gene, (header, sequence) in best.items():
            write_fasta(out, header, sequence)
            table.write(f'{header.split()[0]}\t{gene}\n')
    print(f'{len(best)} representatives')


def counts(a):
    if a.out.with_suffix('.excluded.tsv') == a.out:
        raise ValueError('Choose an output name other than *.excluded.tsv')
    with a.counts.open() as src, a.out.open('w') as out, a.out.with_suffix('.excluded.tsv').open('w') as large:
        reader = csv.reader(src, delimiter='\t')
        header = next(reader)
        drop_total = header[-1].lower() == 'total'
        names = header[1:-1] if drop_total else header[1:]
        writers = [csv.writer(handle, delimiter='\t', lineterminator='\n') for handle in (out, large)]
        for writer in writers:
            writer.writerow(['Desc', 'Family ID'] + names)
        for row in reader:
            if not row:
                continue
            values = [int(x) for x in (row[1:-1] if drop_total else row[1:])]
            if len(values) != len(names) or any(x < 0 for x in values):
                raise ValueError(f'Invalid gene counts: {row[0]}')
            writers[int(max(values) > a.max_count)].writerow(['(null)', row[0]] + values)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    q = sub.add_parser('representatives')
    q.add_argument('--proteins', type=path, required=True)
    q.add_argument('--format', choices=['ncbi', 'gwh'], required=True)
    q.add_argument('--gff3', type=path, help='required for NCBI protein/gene mapping')
    q.add_argument('--out', type=path, required=True)
    q = sub.add_parser('orthofinder')
    q.add_argument('--proteomes', type=path, required=True, help='directory with one FASTA per taxon')
    q.add_argument('--analysis-threads', type=positive, default=8)
    execution_options(q)
    q = sub.add_parser('counts')
    q.add_argument('--counts', type=path, required=True)
    q.add_argument('--out', type=path, required=True, help='CAFE count TSV')
    q.add_argument('--max-count', type=positive, default=100)
    q = sub.add_parser('ultrametric')
    q.add_argument('--tree', type=path, required=True)
    q.add_argument('--root-age', type=float, default=100, help='tree scaling parameter')
    q.add_argument('--utility', type=path, required=True, help='OrthoFinder make_ultrametric.py')
    q.add_argument('--python', default='python3')
    execution_options(q)
    q = sub.add_parser('cafe')
    q.add_argument('--counts', type=path, required=True)
    q.add_argument('--tree', type=path, required=True, help='ultrametric tree')
    execution_options(q)
    a = p.parse_args()
    if a.command == 'representatives':
        if a.format == 'ncbi' and not a.gff3:
            p.error('NCBI representatives require --gff3')
        representatives(a)
    elif a.command == 'counts':
        counts(a)
    else:
        r = Runner(a)
        if a.command == 'orthofinder':
            r.run(['orthofinder', '-f', a.proteomes, '-t', a.threads, '-a', a.analysis_threads,
                   '-S', 'diamond', '-n', 'analysis', '-o', a.out/'results'])
        elif a.command == 'ultrametric':
            r.run([a.python, a.utility, a.tree, '-r', a.root_age, '-o', a.out/'tree.nwk'])
        else:
            r.run(['cafe5', '-i', a.counts, '-t', a.tree, '-p', '-c', a.threads, '-o', a.out/'results'])


if __name__ == '__main__':
    main()
