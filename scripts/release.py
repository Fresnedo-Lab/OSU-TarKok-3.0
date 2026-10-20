#!/usr/bin/env python3
"""Select representative proteins and apply sequence-ID filters/maps to FASTA and GFF3."""
from pathlib import Path
from collections import defaultdict
import argparse
import csv
from _common import path, fasta, write_fasta


def representatives(a):
    from Bio import SeqIO
    from Bio.Seq import Seq
    for out in (a.output, a.selection_table):
        if out.exists():
            raise ValueError(f'Output already exists: {out}')
    genes = set()
    tx_gene = {}
    starts = set()
    stops = set()
    with a.gff3.open() as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) != 9:
                raise ValueError('GFF3 row must have nine fields')
            attrs = dict((x.split('=', 1) for x in c[8].split(';') if '=' in x))
            if c[2] == 'gene':
                genes.add(attrs['ID'])
            elif c[2] in ('mRNA', 'transcript'):
                tx_gene[attrs['ID']] = attrs['Parent']
            elif c[2] == 'start_codon':
                starts.update(attrs['Parent'].split(','))
            elif c[2] == 'stop_codon':
                stops.update(attrs['Parent'].split(','))
    cds = SeqIO.to_dict(SeqIO.parse(a.cds, 'fasta'))
    proteins = SeqIO.to_dict(SeqIO.parse(a.proteins, 'fasta'))
    by_gene = defaultdict(list)
    rejected = 0
    for (tx, gene) in tx_gene.items():
        if gene not in genes:
            raise ValueError(f'Unknown gene parent: {gene}')
        if tx not in cds or tx not in proteins:
            raise ValueError(f'Missing sequence: {tx}')
        dna = str(cds[tx].seq).upper()
        aa = str(proteins[tx].seq).upper().removesuffix('*')
        if not dna or not aa or len(dna) % 3 or set(dna) - set('ACGTRYSWKMBDHVN') or set(aa) - set('ACDEFGHIKLMNPQRSTVWYBXZJUO') or ('*' in aa):
            rejected += 1
            continue
        translated = str(Seq(dna).translate(table=1)).removesuffix('*')
        if '*' in translated or translated != aa:
            rejected += 1
            continue
        by_gene[gene].append((tx, tx in starts and tx in stops, len(aa)))
    chosen = []
    for gene in sorted(genes):
        choices = by_gene[gene]
        if not choices:
            raise ValueError(f'No translation-valid isoform for {gene}')
        complete = [x for x in choices if x[1]]
        (tx, is_complete, length) = sorted(complete or choices, key=lambda x: (-x[2], x[0]))[0]
        chosen.append((gene, tx, length, 'complete' if is_complete else 'partial'))
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.selection_table.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write([proteins[tx] for (_, tx, _, _) in chosen], a.output, 'fasta')
    with a.selection_table.open('w') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['gene_id', 'protein_id', 'amino_acid_length', 'boundary_status'])
        w.writerows(chosen)
    print(f"{len(chosen)} representatives; {sum((x[3] == 'partial' for x in chosen))} partial fallbacks; {rejected} invalid candidates.")


def subset(a):
    keep = set(a.keep.read_text().split()) if a.keep else None
    rename = {}
    if a.map:
        for line in a.map.read_text().splitlines():
            if not line.strip() or line.startswith('#'):
                continue
            old, new = line.split('\t')
            if old in rename or not new or any(c.isspace() for c in new):
                raise ValueError('Duplicate source ID or invalid destination ID in map')
            rename[old] = new
    a.out.mkdir(parents=True, exist_ok=False)
    retained, lengths, output_ids = set(), {}, set()
    with (a.out/'genome.fa').open('w') as out:
        for header, sequence in fasta(a.genome):
            name = header.split()[0]
            if keep is not None and name not in keep:
                continue
            target = rename.get(name, name)
            if target in output_ids:
                raise ValueError(f'Duplicate destination ID: {target}')
            retained.add(name)
            output_ids.add(target)
            lengths[name] = len(sequence)
            write_fasta(out, target, sequence)
    if keep is not None and keep != retained:
        raise ValueError('Some retained IDs are absent from the FASTA')
    if not retained:
        raise ValueError('No sequences selected')
    if a.gff3:
        with a.gff3.open() as src, (a.out/'genes.gff3').open('w') as out:
            out.write('##gff-version 3\n')
            for name, length in lengths.items():
                out.write(f'##sequence-region {rename.get(name, name)} 1 {length}\n')
            for line in src:
                if line.startswith('##FASTA'):
                    break
                if line.startswith('#') or not line.strip():
                    continue
                f = line.rstrip('\n').split('\t')
                if f[0] not in retained:
                    continue
                if not 1 <= int(f[3]) <= int(f[4]) <= lengths[f[0]]:
                    raise ValueError(f'Invalid coordinates: {f[:5]}')
                f[0] = rename.get(f[0], f[0])
                out.write('\t'.join(f) + '\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(required=True)
    q = sub.add_parser('representatives', help='longest translation-valid complete CDS per gene; partial fallback')
    for flag in ('gff3', 'cds', 'proteins', 'output', 'selection-table'):
        q.add_argument('--' + flag, type=Path, required=True)
    q.set_defaults(func=representatives)
    q = sub.add_parser('subset', help='retain selected sequences and optionally rename FASTA/GFF sequence IDs')
    q.add_argument('--genome', type=path, required=True)
    q.add_argument('--gff3', type=path)
    q.add_argument('--keep', type=path, help='one retained sequence ID per line')
    q.add_argument('--map', type=path, help='headerless TSV: input ID, output ID')
    q.add_argument('--out', type=path, required=True)
    q.set_defaults(func=subset)
    a = p.parse_args()
    a.func(a)


if __name__ == '__main__':
    main()
