#!/usr/bin/env python3
"""RepeatModeler/RepeatMasker, EDTA, repeat classification and GFF3 export."""
import argparse
from collections import defaultdict
import re
import shutil
from urllib.parse import quote, unquote
from _common import Runner, execution_options, path, fasta


def run(a):
    r = Runner(a)
    if a.command == 'model':
        r.run(['BuildDatabase', '-name', 'genome', a.genome])
        r.run(['RepeatModeler', '-database', 'genome', '-threads', a.threads])
    elif a.command == 'mask':
        library = a.out/'repeats.fa'
        if not r.dry:
            with library.open('w') as out:
                for file in a.libraries:
                    with file.open() as src:
                        shutil.copyfileobj(src, out)
                    out.write('\n')
        r.run(['RepeatMasker', '-pa', max(1, a.threads//2), '-e', 'rmblast', '-lib', library,
               '-xsmall', '-gff', '-a', '-dir', a.out, a.genome])
    elif a.command == 'edta':
        if not r.dry:
            shutil.copyfile(a.genome, a.out/'genome.fa')
        r.run(['EDTA.pl', '--genome', 'genome.fa', '--species', 'others', '--step', 'all',
               '--sensitive', '0', '--anno', '1', '--curatedlib', a.library, '--threads', a.threads, '--force', '1'])
        if a.lai:
            r.run(['LAI', '-genome', 'genome.fa.mod', '-intact', 'genome.fa.mod.EDTA.raw/LTR/genome.fa.mod.pass.list',
                   '-all', 'genome.fa.mod.EDTA.anno/genome.fa.mod.out', '-t', a.threads])
    elif a.command == 'classify':
        if not r.dry:
            shutil.copyfile(a.library, a.out/'families.fa')
        r.run(['TEsorter', 'families.fa', '-db', 'rexdb-plant', '-st', 'nucl', '-p', a.threads, '-pre', 'rexdb'])
        if a.repeatclassifier:
            r.run(['RepeatClassifier', '-consensi', 'families.fa', '-pa', a.threads])


def gff_rows(file):
    with file.open() as src:
        for line in src:
            if line.startswith('##FASTA'):
                break
            if line.startswith('#') or not line.strip():
                continue
            f = [x.strip() for x in line.rstrip('\n').split('\t')]
            if len(f) != 9:
                raise ValueError(f'Invalid GFF row in {file}')
            yield f


def encode(value):
    return quote(value, safe=':@+._-!/()|^*$?<>[]{}~')


def export(a):
    lengths = {h.split()[0]: len(s) for h, s in fasta(a.genome)}
    def check(f):
        if f[0] not in lengths or not 1 <= int(f[3]) <= int(f[4]) <= lengths[f[0]]:
            raise ValueError(f'GFF coordinates outside reference: {f[:5]}')
    with a.out.open('w') as out:
        out.write('##gff-version 3\n')
        for name, length in lengths.items():
            out.write(f'##sequence-region {name} 1 {length}\n')
        number = 0
        for part, file in enumerate(a.inputs, 1):
            for f in gff_rows(file):
                check(f)
                number += 1
                if a.format == 'edta':
                    attrs = []
                    for field in f[8].split(';'):
                        if not field:
                            continue
                        key, value = field.split('=', 1)
                        if key in ('ID', 'Parent'):
                            value = ','.join(encode(f'part{part:03d}_{unquote(v)}') for v in value.split(','))
                        else:
                            value = encode(unquote(value))
                        attrs.append(key + '=' + value)
                    f[8] = ';'.join(attrs)
                else:
                    match = re.fullmatch(r'(?:ID=([^;]+);)?Target "([^"]+)" (-?\d+) (-?\d+)', f[8])
                    if not match:
                        raise ValueError(f'Unrecognized RepeatMasker attributes: {f[8]}')
                    group, target, start, end = match.groups()
                    start, end = int(start), int(end)
                    attrs = [f'ID=RM_{number:07d}', 'Name=' + encode(target.removeprefix('Motif:'))]
                    if group:
                        attrs.append(f'repeat_group=part{part:03d}_' + encode(group))
                    if start >= 1 and end >= start:
                        attrs.append(f'Target={encode(target)} {start} {end}')
                    else:
                        attrs += [f'repeat_consensus_start={start}', f'repeat_consensus_end={end}']
                    f[8] = ';'.join(attrs)
                out.write('\t'.join(f) + '\n')
    print(f'{number} repeat features')


def summary(a):
    total = sum(len(s) for _, s in fasta(a.genome))
    if not total:
        raise ValueError('Empty genome')
    classes, masked, interspersed = defaultdict(int), 0, 0
    for file in a.summaries:
        parent = None
        with file.open() as src:
            for line in src:
                if line.startswith('Repeat Stats'):
                    break
                match = re.match(r'^(LTR|TIR|nonTIR)\s+--', line)
                if match:
                    parent = match[1]
                match = re.match(r'\s+(\S+)\s+\d+\s+(\d+)\s+[\d.]+%', line)
                if match and parent:
                    classes[(parent, match[1])] += int(match[2])
                match = re.match(r'^repeat_fragment\s+\d+\s+(\d+)', line)
                if match:
                    classes[('other', 'repeat_fragment')] += int(match[1])
                match = re.match(r'\s*total interspersed\s+\d+\s+(\d+)', line)
                if match:
                    interspersed += int(match[1])
                match = re.match(r'^Total\s+\d+\s+(\d+)\s+[\d.]+%', line)
                if match:
                    masked += int(match[1])
    with a.out.open('w') as out:
        out.write('class\tsuperfamily\tbpMasked\tpercent_genome\n')
        for (parent, family), bp in sorted(classes.items()):
            out.write(f'{parent}\t{family}\t{bp}\t{100*bp/total:.6f}\n')
        out.write(f'total\tinterspersed\t{interspersed}\t{100*interspersed/total:.6f}\n')
        out.write(f'total\tall_repeats\t{masked}\t{100*masked/total:.6f}\n')
    if a.lai_files:
        lengths, weighted = 0, 0
        for file in a.lai_files:
            for line in file.read_text().splitlines():
                if line.startswith('whole_genome'):
                    fields = line.split()
                    length, lai = int(fields[2]), float(fields[6])
                    lengths += length
                    weighted += length * lai
                    break
            else:
                raise ValueError(f'No whole_genome LAI row in {file}')
        if not lengths:
            raise ValueError('Zero total LAI length')
        a.out.with_suffix('.lai.tsv').write_text(f'metric\tvalue\nlength_weighted_LAI\t{weighted/lengths:.6f}\n')


def bucket(cls):
    if not cls:
        return 'Unknown'
    c = cls.lower()
    if 'copia' in c:
        return 'LTR/Copia'
    if 'gypsy' in c:
        return 'LTR/Gypsy'
    if 'caulimovir' in c or 'pararetro' in c:
        return 'LTR/Caulimovirus'
    if 'ltr' in c and 'retro' in c:
        return 'LTR/unknown'
    if c.startswith('ltr') or '/ltr' in c:
        return 'LTR/unknown'
    if 'line' in c or 'l1' in c or 'cre' in c or ('rte' in c) or ('jockey' in c):
        return 'LINE'
    if 'sine' in c:
        return 'SINE'
    if 'helitron' in c or c.startswith('rc') or '/rc' in c:
        return 'Helitron/RC'
    if any((k in c for k in ['mule', 'mudr', 'harbinger', 'pif', 'enspm', 'en-spm', 'cmc', 'hat', 'mariner', 'tc1', 'tir', 'dna', 'mutator', 'cacta', 'p-element', 'transib', 'tag1', 'tip100', 'dtm', 'dtt', 'dth', 'dtc', 'dta', 'dtx'])):
        return 'TIR/DNA'
    if 'rrna' in c or 'trna' in c or 'snrna' in c or ('rdna' in c):
        return 'Other(ncRNA/host)'
    if 'satellite' in c or 'simple' in c or 'low_complex' in c:
        return 'Other(simple/sat)'
    if 'unknown' in c or 'unspecified' in c or 'na' == c:
        return 'Unknown'
    return 'Other(' + cls + ')'

def tesorter_clsstr(t):
    (order, supf, clade) = t
    if supf and supf.lower() not in ('unknown', 'mixture', 'na', ''):
        s = supf
    else:
        s = order
    if clade and clade.lower() not in ('unknown', 'na', ''):
        s = '%s/%s' % (s, clade)
    return s


def families(a):
    import csv
    classified = {}
    for header, sequence in (fasta(a.classified) if a.classified else []):
        name, _, kind = header.split()[0].partition('#')
        classified[name] = kind or 'Unknown'
    tesorter = {}
    for line in a.tesorter.read_text().splitlines():
        if line.startswith('#') or not line.strip():
            continue
        fields = line.split('\t')
        if len(fields) >= 4:
            tesorter[fields[0].split('#')[0]] = tuple(fields[1:4])
    with a.out.open('w', newline='') as out:
        writer = csv.writer(out, delimiter='\t', lineterminator='\n')
        writer.writerow(['family', 'original_class', 'tesorter_rexdb', 'repeatclassifier_dfam',
                         'final_class', 'final_superfamily', 'assigned_by', 'sequence_length'])
        for header, sequence in fasta(a.library):
            name, _, original = header.split()[0].partition('#')
            original = original or 'Unknown'
            ts = tesorter_clsstr(tesorter[name]) if name in tesorter else ''
            rc = classified.get(name, '')
            final, source = 'Unknown', 'none'
            for label, kind in [('RepeatModeler', original), ('TEsorter/REXdb', ts), ('RepeatClassifier/Dfam', rc)]:
                if bucket(kind) != 'Unknown':
                    final, source = kind, label
                    break
            writer.writerow([name, original, ts, rc, final, bucket(final), source, len(sequence)])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    for name in ('model', 'mask', 'edta', 'classify'):
        q = sub.add_parser(name)
        execution_options(q)
        q.set_defaults(func=run)
        if name != 'classify':
            q.add_argument('--genome', type=path, required=True)
        if name in ('edta', 'classify'):
            q.add_argument('--library', type=path, required=True)
        if name == 'mask':
            q.add_argument('--libraries', type=path, nargs='+', required=True, help='de novo and curated repeat FASTAs')
        if name == 'edta':
            q.add_argument('--lai', action='store_true')
        if name == 'classify':
            q.add_argument('--repeatclassifier', action='store_true', help='also run RepeatClassifier')
    q = sub.add_parser('export', help='combine EDTA GFF3 or convert RepeatMasker GFF to GFF3')
    q.add_argument('--format', choices=['edta', 'repeatmasker'], required=True)
    q.add_argument('--genome', type=path, required=True)
    q.add_argument('--inputs', type=path, nargs='+', required=True)
    q.add_argument('--out', type=path, required=True)
    q.set_defaults(func=export)
    q = sub.add_parser('summary', help='sum EDTA masked bases using the supplied genome length')
    q.add_argument('--genome', type=path, required=True)
    q.add_argument('--summaries', type=path, nargs='+', required=True)
    q.add_argument('--lai-files', type=path, nargs='+', help='chromosome LAI files, excluding unplaced sequences')
    q.add_argument('--out', type=path, required=True)
    q.set_defaults(func=summary)
    q = sub.add_parser('families', help='combine repeat-family classifications')
    for flag in ('library', 'tesorter', 'out'):
        q.add_argument('--' + flag, type=path, required=True)
    q.add_argument('--classified', type=path, help='optional RepeatClassifier FASTA')
    q.set_defaults(func=families)
    a = p.parse_args()
    a.func(a)


if __name__ == '__main__':
    main()
