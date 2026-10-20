#!/usr/bin/env python3
"""Curate OH-1004 plastome features and export nucleotide/protein FASTAs."""
from __future__ import annotations
import argparse
import re
from collections import Counter
from pathlib import Path


FEATURE_RE = re.compile('^     (\\S+)\\s+(.+)$')

QUAL_RE = re.compile('^\\s{21}/([^=]+)(?:=(.*))?$')

INTERVAL_RE = re.compile('([<>]?\\d+)\\.\\.([<>]?\\d+)|([<>]?\\d+)')

def read_fasta(path: Path) -> str:
    sequence = []
    headers = 0
    with path.open() as handle:
        for line in handle:
            if line.startswith('>'):
                headers += 1
            else:
                sequence.append(line.strip())
    if headers != 1:
        raise ValueError('Expected a single plastome FASTA record')
    seq = ''.join(sequence).upper()
    if not seq or re.search('[^ACGTN]', seq):
        raise ValueError(f'Unexpected or empty FASTA sequence in {path}')
    return seq

def parse_genbank_features(path: Path):
    features = []
    current = None
    pending_qual = None
    in_features = False

    def finish_qualifier():
        nonlocal pending_qual
        if pending_qual is not None:
            (key, value) = pending_qual
            current['qualifiers'].append((key, value))
            pending_qual = None
    with path.open() as handle:
        for raw in handle:
            line = raw.rstrip('\n')
            if line.startswith('FEATURES'):
                in_features = True
                continue
            if line.startswith('ORIGIN'):
                finish_qualifier()
                if current:
                    features.append(current)
                break
            if not in_features:
                continue
            match = FEATURE_RE.match(line)
            if match:
                finish_qualifier()
                if current:
                    features.append(current)
                current = {'type': match.group(1), 'location': match.group(2).strip(), 'qualifiers': []}
                continue
            if current is None:
                continue
            qualifier = QUAL_RE.match(line)
            if qualifier:
                finish_qualifier()
                (key, value) = qualifier.groups()
                if value is None:
                    current['qualifiers'].append((key, None))
                elif value.startswith('"') and (not value.endswith('"')):
                    pending_qual = [key, value]
                else:
                    current['qualifiers'].append((key, value))
                continue
            continuation = line[21:].strip() if len(line) >= 21 else ''
            if pending_qual is not None:
                separator = '' if pending_qual[0] == 'translation' else ' '
                pending_qual[1] += separator + continuation
                if continuation.endswith('"'):
                    finish_qualifier()
            elif continuation and (not continuation.startswith('/')):
                current['location'] += continuation
    return features

def split_top_level(text: str):
    (parts, depth, start) = ([], 0, 0)
    for (index, char) in enumerate(text):
        if char == '(':
            depth += 1
        elif char == ')':
            depth -= 1
        elif char == ',' and depth == 0:
            parts.append(text[start:index])
            start = index + 1
    parts.append(text[start:])
    return [part.strip() for part in parts]

def parse_location(text: str, strand: int=1):
    text = text.replace(' ', '')
    if text.startswith('complement(') and text.endswith(')'):
        pieces = parse_location(text[11:-1], -strand)
        return list(reversed(pieces))
    for operator in ('join', 'order'):
        prefix = operator + '('
        if text.startswith(prefix) and text.endswith(')'):
            result = []
            for part in split_top_level(text[len(prefix):-1]):
                result.extend(parse_location(part, strand))
            return result
    match = INTERVAL_RE.fullmatch(text)
    if not match:
        raise ValueError(f'Unsupported location: {text}')
    if match.group(3):
        start = end = match.group(3)
    else:
        (start, end) = (match.group(1), match.group(2))
    return [(start, end, strand)]

def clean_value(value: str | None):
    if value is None:
        return None
    if value.startswith('"') and value.endswith('"'):
        value = value[1:-1]
    return value.replace('""', '"').replace('\n', ' ')

def gene_name(feature):
    for (key, value) in feature['qualifiers']:
        if key == 'gene':
            return clean_value(value)
    return None

def adjust_feature(feature):
    if gene_name(feature) == 'ndhF' and feature['type'] in {'gene', 'CDS'}:
        feature = {'type': feature['type'], 'location': feature['location'].replace('108331..110535', '108322..110535'), 'qualifiers': list(feature['qualifiers'])}
    curated_products = {'accD': 'acetyl-CoA carboxylase beta subunit', 'psaA': 'photosystem I chlorophyll a apoprotein A1', 'psaB': 'photosystem I chlorophyll a apoprotein A2', 'psbB': 'photosystem II CP47 reaction center protein', 'ycf2': 'Ycf2 protein'}
    name = gene_name(feature)
    if feature['type'] == 'CDS' and name in curated_products:
        feature = {'type': feature['type'], 'location': feature['location'], 'qualifiers': [(key, f'"{curated_products[name]}"') if key == 'product' else (key, value) for (key, value) in feature['qualifiers']]}
    return feature

def write_tbl(features, path: Path, seqid: str):
    allowed = {'gene', 'product', 'codon_start', 'transl_table', 'trans_splicing', 'note', 'rpt_type', 'exception', 'transl_except'}
    skipped = {'protein_id', 'translation', 'db_xref', 'organism', 'organelle', 'mol_type', 'tissue_type', 'geo_loc_name'}
    duplicate_genes = {name for (name, count) in Counter((gene_name(feature) for feature in features if feature['type'] == 'gene' and gene_name(feature))).items() if count > 1}
    with path.open('w') as out:
        out.write(f'>Feature {seqid}\n')
        for original in features:
            if original['type'] == 'source':
                continue
            feature = adjust_feature(original)
            intervals = parse_location(feature['location'])
            for (index, (left, right, strand)) in enumerate(intervals):
                if strand == -1:
                    (left, right) = (right, left)
                feature_type = feature['type'] if index == 0 else ''
                out.write(f'{left}\t{right}\t{feature_type}\n')
            for (key, raw_value) in feature['qualifiers']:
                if key in skipped:
                    continue
                if key not in allowed:
                    raise ValueError(f"Unreviewed qualifier /{key} on {feature['type']}")
                value = clean_value(raw_value)
                if key == 'exception' and feature['type'] == 'CDS' and (gene_name(feature) == 'psbL'):
                    continue
                if key == 'gene' and feature['type'] != 'gene' and (value in duplicate_genes):
                    continue
                if value is None:
                    out.write(f'\t\t\t{key}\n')
                else:
                    out.write(f'\t\t\t{key}\t{value}\n')
            if feature['type'] == 'CDS' and gene_name(feature) == 'psbL':
                out.write('\t\t\ttransl_except\t(pos:complement(63956..63958),aa:Met)\n')

def write_fasta(sequence: str, path: Path, seqid: str):
    header = f'>{seqid} [organism=Taraxacum kok-saghyz] [cultivar=OH-1004] [location=chloroplast] [topology=circular] [completeness=complete] [gcode=11] [db_xref=taxon:333970]'
    with path.open('w') as out:
        out.write(header + '\n')
        for offset in range(0, len(sequence), 60):
            out.write(sequence[offset:offset + 60] + '\n')

def curate(args):
    sequence = read_fasta(args.target_fasta)
    if len(sequence) != 151316:
        raise ValueError(f'Expected 151316 bp plastome, got {len(sequence)}')
    features = parse_genbank_features(args.reference_gb)
    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    write_fasta(sequence, args.out_prefix.with_suffix('.fsa'), args.seqid)
    write_tbl(features, args.out_prefix.with_suffix('.tbl'), args.seqid)

COMPLEMENT = str.maketrans('ACGTNacgtn', 'TGCANtgcan')

def origin_sequence(path: Path):
    text = path.read_text()
    origin = text.split('ORIGIN', 1)[1].split('//', 1)[0]
    return ''.join(re.findall('[ACGTNacgtn]+', origin)).upper()

def qualifier(feature, requested):
    for (key, value) in feature['qualifiers']:
        if key == requested:
            return clean_value(value)
    return ''

def feature_sequence(sequence, location):
    pieces = []
    for (start, end, strand) in parse_location(location):
        left = int(start.lstrip('<>'))
        right = int(end.lstrip('<>'))
        piece = sequence[left - 1:right]
        if strand == -1:
            piece = piece.translate(COMPLEMENT)[::-1]
        pieces.append(piece)
    return ''.join(pieces)

def write_record(handle, header, sequence, width=60):
    handle.write('>' + header + '\n')
    for offset in range(0, len(sequence), width):
        handle.write(sequence[offset:offset + width] + '\n')

def export(args):
    args.output_dir.mkdir(parents=True, exist_ok=True)
    sequence = origin_sequence(args.genbank_preview)
    features = parse_genbank_features(args.genbank_preview)
    handles = {'gene': (args.output_dir / 'chloroplast.genes_with_introns.fa').open('w'), 'CDS': (args.output_dir / 'chloroplast.cds.fa').open('w'), 'protein': (args.output_dir / 'chloroplast.proteins.fa').open('w'), 'tRNA': (args.output_dir / 'chloroplast.tRNA.fa').open('w'), 'rRNA': (args.output_dir / 'chloroplast.rRNA.fa').open('w')}
    counts = Counter()
    copies = Counter()
    try:
        for feature in features:
            kind = feature['type']
            if kind not in {'gene', 'CDS', 'tRNA', 'rRNA'}:
                continue
            gene = gene_name(feature) or 'unnamed'
            copies[kind, gene] += 1
            copy_number = copies[kind, gene]
            product = qualifier(feature, 'product')
            ident = f'Tks_OH1004_chloroplast|{kind}|{gene}|copy{copy_number}'
            header = f"{ident} location={feature['location']}"
            if product:
                header += f' product={product}'
            nucleotide = feature_sequence(sequence, feature['location'])
            write_record(handles[kind], header, nucleotide)
            counts[kind] += 1
            if kind == 'CDS':
                translation = qualifier(feature, 'translation')
                if not translation:
                    raise ValueError(f'CDS {gene} lacks table2asn translation')
                write_record(handles['protein'], header.replace('|CDS|', '|protein|'), translation)
                counts['protein'] += 1
    finally:
        for handle in handles.values():
            handle.close()
    with (args.output_dir / 'feature_fasta_counts.tsv').open('w') as out:
        out.write('feature_type\tcount\n')
        for kind in ('gene', 'CDS', 'protein', 'tRNA', 'rRNA'):
            out.write(f'{kind}\t{counts[kind]}\n')

def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(required=True)
    q = sub.add_parser('curate-oh1004', help='OR990541.1 features on the coordinate-matched 151316 bp OH-1004 plastome')
    q.add_argument('--reference-gb', type=Path, required=True)
    q.add_argument('--target-fasta', type=Path, required=True)
    q.add_argument('--out-prefix', type=Path, required=True)
    q.add_argument('--seqid', default='Tks_OH1004_chloroplast')
    q.set_defaults(func=curate)
    q = sub.add_parser('export', help='extract genes, CDSs, proteins, tRNAs and rRNAs from GenBank')
    q.add_argument('--genbank-preview', type=Path, required=True)
    q.add_argument('--output-dir', type=Path, required=True)
    q.set_defaults(func=export)
    a = p.parse_args()
    a.func(a)

if __name__ == '__main__':
    main()
