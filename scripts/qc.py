#!/usr/bin/env python3
"""Assembly statistics, BUSCO, QUAST and RNA-seq alignment rates."""
import argparse
import os
import re
from _common import Runner, execution_options, path, fasta


def nx(lengths, fraction):
    target, running = sum(lengths) * fraction, 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running >= target:
            return length
    return 0


def stats(a):
    scaffold_lengths, contig_lengths = [], []
    total_n = total_gc = total_acgt = n_runs = anchored = 0
    chromosomes = set(a.chromosomes.read_text().split()) if a.chromosomes else set()
    seen = set()
    for header, sequence in fasta(a.genome):
        name = header.split()[0]
        if name in seen:
            raise ValueError(f'Duplicate sequence ID: {name}')
        seen.add(name)
        sequence = sequence.upper()
        scaffold_lengths.append(len(sequence))
        contig_lengths.extend(len(part) for part in re.split('N+', sequence) if part)
        total_n += sequence.count('N')
        total_gc += sequence.count('G') + sequence.count('C')
        total_acgt += sum(sequence.count(x) for x in 'ACGT')
        n_runs += len(re.findall('N+', sequence))
        if name in chromosomes:
            anchored += len(sequence)
    if not chromosomes <= seen:
        raise ValueError('Chromosome IDs absent from FASTA')
    total = sum(scaffold_lengths)
    if not total:
        raise ValueError('Empty assembly')
    values = {
        'assembly_length_bp': total, 'fasta_sequences_scaffolds': len(scaffold_lengths),
        'contigs_split_at_N_runs': len(contig_lengths),
        'scaffold_N50_bp': nx(scaffold_lengths, .5), 'scaffold_N90_bp': nx(scaffold_lengths, .9),
        'contig_N50_bp': nx(contig_lengths, .5), 'contig_N90_bp': nx(contig_lengths, .9),
        'largest_scaffold_bp': max(scaffold_lengths), 'largest_contig_bp': max(contig_lengths, default=0),
        'N_runs': n_runs, 'N_bases': total_n, 'N_bases_per_100kb': total_n/total*100000,
        'GC_percent_full_sequence': total_gc/total*100,
        'GC_percent_ACGT_only': total_gc/total_acgt*100 if total_acgt else 0,
        'ACGT_bases': total_acgt, 'other_ambiguous_non_N_bases': total-total_acgt-total_n,
    }
    if chromosomes:
        values.update(chromosome_scale_scaffolds=len(chromosomes), chromosome_scale_span_bp=anchored,
                      fraction_anchored_percent=anchored/total*100, unplaced_span_bp=total-anchored,
                      unplaced_fraction_percent=(total-anchored)/total*100)
    a.out.write_text('metric\tvalue\n' + ''.join(f'{key}\t{value}\n' for key, value in values.items()))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    q = sub.add_parser('stats')
    q.add_argument('--genome', type=path, required=True)
    q.add_argument('--chromosomes', type=path, help='optional chromosome IDs, one per line')
    q.add_argument('--out', type=path, required=True, help='metrics TSV')
    q = sub.add_parser('busco')
    q.add_argument('--input', type=path, required=True)
    q.add_argument('--mode', choices=['genome', 'proteins'], required=True)
    q.add_argument('--lineage', type=path, required=True, help='local BUSCO lineage directory')
    execution_options(q)
    q = sub.add_parser('quast')
    q.add_argument('--genomes', type=path, nargs='+', required=True)
    execution_options(q)
    q = sub.add_parser('maprate')
    for flag in ('genome', 'r1', 'r2'):
        q.add_argument('--' + flag, type=path, required=True)
    execution_options(q)
    a = p.parse_args()
    if a.command == 'stats':
        stats(a)
        return
    r = Runner(a)
    if a.command == 'busco':
        cmd = ['busco', '-i', a.input, '-l', a.lineage, '-m', a.mode, '-c', a.threads,
               '-o', 'busco', '--out_path', a.out, '--offline', '--opt-out-run-stats']
        if a.mode == 'genome':
            cmd += ['--miniprot']
        r.run(cmd)
    elif a.command == 'quast':
        r.run(['quast.py', '--large', '--min-contig', '1000', '-t', a.threads, '-o', a.out/'results', *a.genomes])
    else:
        r.run(['hisat2-build', '-p', a.threads, a.genome, a.out/'index'])
        r.run(['hisat2', '-p', a.threads, '-x', a.out/'index', '-1', a.r1, '-2', a.r2,
               '--no-unal', '-S', os.devnull])


if __name__ == '__main__':
    main()
