#!/usr/bin/env python3
"""Reference-guided HiFi plastome consensus and PGA annotation."""
import argparse
from pathlib import Path
import shutil
import sys
from _common import Runner, execution_options, path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    q = sub.add_parser('consensus')
    q.add_argument('--reads', type=path, nargs='+', required=True)
    q.add_argument('--reference', type=path, required=True, help='single plastome FASTA')
    q.add_argument('--fraction', type=float, default=0.043, help='read fraction retained for variant calling')
    q.add_argument('--seed', type=int, default=0)
    execution_options(q)
    q = sub.add_parser('annotate')
    q.add_argument('--genome', type=path, required=True)
    q.add_argument('--references', type=path, required=True, help='directory of reference GenBank files')
    q.add_argument('--pga', type=path, required=True, help='PGA.pl')
    execution_options(q)
    a = p.parse_args()
    if a.command == 'consensus' and not 0 < a.fraction <= 1:
        p.error('--fraction must be in (0, 1]')
    r = Runner(a)
    if a.command == 'annotate':
        target = a.out/'target'
        if not r.dry:
            target.mkdir()
            shutil.copyfile(a.genome, target/'plastome.fa')
        r.run(['perl', a.pga, '-r', a.references, '-t', target, '-o', a.out/'annotation', '-f', 'circular'])
        return
    ref = a.out/'reference.fa'
    if not r.dry:
        shutil.copyfile(a.reference, ref)
    r.run(['samtools', 'faidx', ref])
    recruited = []
    for n, reads in enumerate(a.reads, 1):
        fq = a.out/f'recruited_{n}.fq'
        r.pipe([['minimap2', '-t', a.threads, '-ax', 'map-hifi', ref, reads],
                ['samtools', 'view', '-b', '-F', '2308', '-'], ['samtools', 'fastq', '-']], stdout=fq)
        recruited.append(fq)
    raw = a.out/'recruited.fq'
    if not r.dry:
        with raw.open('wb') as out:
            for file in recruited:
                with file.open('rb') as src:
                    shutil.copyfileobj(src, out)
    reads = a.out/'reads.fq.gz'
    r.pipe([['seqkit', 'rmdup', '-n', raw], ['gzip', '-c']], stdout=reads)
    bam = a.out/'reference.bam'
    r.pipe([['minimap2', '-t', a.threads, '-ax', 'map-hifi', ref, reads],
            ['samtools', 'sort', '-o', bam, '-']])
    r.run(['samtools', 'index', bam])
    subset = a.out/'subsample.bam'
    r.run(['samtools', 'view', '-b', '--subsample', a.fraction, '--subsample-seed', a.seed, bam, '-o', subset])
    r.run(['samtools', 'index', subset])
    calls, filtered = a.out/'calls.vcf.gz', a.out/'calls.filtered.vcf.gz'
    # MAPQ 0 retains reads that map equally well to both inverted repeats.
    r.pipe([['bcftools', 'mpileup', '-d', '600', '-q', '0', '-Q', '13', '-f', ref, subset],
            ['bcftools', 'call', '--ploidy', '1', '-m', '-v', '-Oz', '-o', calls]])
    r.run(['bcftools', 'index', calls])
    r.run(['bcftools', 'view', '-e', 'QUAL<30 || INFO/DP<50', calls, '-Oz', '-o', filtered])
    r.run(['bcftools', 'index', filtered])
    raw_consensus = a.out/'consensus.fa'
    r.run(['bcftools', 'consensus', '-f', ref, filtered], stdout=raw_consensus)
    helper = Path(__file__).resolve().with_name('plastome_analysis.py')
    plastome = a.out/'plastome.fa'
    r.run([sys.executable, helper, 'orient', '--fasta', raw_consensus, '--reference', ref, '--out', plastome])
    r.run([sys.executable, helper, 'structure', '--fasta', plastome, '--out', a.out/'structure.tsv'])
    r.pipe([['minimap2', '-t', a.threads, '-ax', 'map-hifi', plastome, reads],
            ['samtools', 'sort', '-o', a.out/'plastome.bam', '-']])
    r.run(['samtools', 'index', a.out/'plastome.bam'])
    r.run(['samtools', 'depth', '-a', a.out/'plastome.bam'], stdout=a.out/'depth.tsv')


if __name__ == '__main__':
    main()
