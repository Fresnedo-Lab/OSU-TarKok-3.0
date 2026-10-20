#!/usr/bin/env python3
"""Convert annotation outputs to TSV and join them onto a proteome."""
import argparse
import csv
import glob
import os
import re
from collections import Counter, defaultdict
from pathlib import Path
from _common import path, fasta, attributes


DB_PRIORITY = ['sprot', 'arabidopsis', 'lactuca']

DB_LABEL = {'sprot': 'UniProt_SwissProt', 'arabidopsis': 'Arabidopsis_TAIR', 'lactuca': 'Lactuca_sativa'}

def clean_stitle(stitle, tag):
    """Extract a human-readable protein name from a diamond stitle (col13)."""
    s = stitle.strip()
    if not s:
        return ''
    if tag == 'sprot':
        s = re.sub('^\\S+\\s+', '', s)
        s = re.split('\\s+(?:OS|OX|GN|PE|SV)=', s)[0]
    else:
        m = re.search('\\bdescription:(.+)$', s)
        if m:
            s = m.group(1)
            s = re.split('\\s*\\[Source:', s)[0]
        else:
            s = re.sub('^\\S+\\s+', '', s)
            s = re.split('\\s+\\[[^\\]]*\\]\\s*$', s)[0]
            s = re.split('\\s+(?:gene:|transcript:|gene_biotype:|transcript_biotype:|pep:)', s)[0]
    s = s.strip()
    return s

def parse_blast(blast_dir):
    """Return dict tag -> {qid -> (acc, evalue_float, stitle)} best hit per query."""
    best = {}
    for tag in DB_PRIORITY:
        f13 = os.path.join(blast_dir, f'{tag}.blastp.tsv')
        f = f13 if os.path.isfile(f13) and os.path.getsize(f13) else os.path.join(blast_dir, f'{tag}.blastp.b8.tsv')
        d = {}
        if os.path.isfile(f) and os.path.getsize(f):
            with open(f) as fh:
                for line in fh:
                    c = line.rstrip('\n').split('\t')
                    if len(c) < 12:
                        continue
                    (qid, sid) = (c[0], c[1])
                    try:
                        ev = float(c[10])
                        bs = float(c[11])
                    except ValueError:
                        continue
                    stitle = c[12] if len(c) >= 13 else ''
                    prev = d.get(qid)
                    if prev is None or ev < prev[1] or (ev == prev[1] and bs > prev[3]):
                        d[qid] = (sid, ev, stitle, bs)
        best[tag] = d
    return best

def best_hit_for(qid, best):
    """Return (source_db_label, acc, evalue_str) using DB priority then evalue."""
    for tag in DB_PRIORITY:
        hit = best.get(tag, {}).get(qid)
        if hit:
            return (DB_LABEL[tag], hit[0], f'{hit[1]:.2e}', tag, hit[2])
    return ('', '', '', '', '')

FAMILY_RE = re.compile('^(GH|GT|PL|CE|AA|CBM|SLH|cohesin|dockerin)', re.IGNORECASE)

def norm_family(tok):
    """Return the bare CAZy family from one token, or None.
    'GH32(1-200)' -> 'GH32' ; 'GH32_e123' -> 'GH32' ; 'GH13_31' -> 'GH13'.
    """
    if tok is None:
        return None
    tok = tok.strip()
    if tok in ('', '-', 'N', 'NA', 'None'):
        return None
    tok = re.sub('\\(.*?\\)', '', tok)
    tok = tok.split('+')[0].strip()
    if not tok or tok in ('-',):
        return None
    if not FAMILY_RE.match(tok):
        return None
    m = re.match('^([A-Za-z]+)(\\d+)?', tok)
    if not m:
        return None
    (letters, digits) = (m.group(1), m.group(2))
    return letters.upper() + (digits or '')

def families_from_cell(cell):
    """All distinct bare families mentioned in one tool's cell."""
    fams = []
    if cell is None:
        return fams
    for tok in re.split('[|,;\\s]+', cell.strip()):
        f = norm_family(tok)
        if f and f not in fams:
            fams.append(f)
    return fams

def subfamilies_from_cell(cell):
    """dbCAN_sub subfamily tokens like 'GH32_e123' or 'GH13_31'. Keep full token (minus parens)."""
    subs = []
    if cell is None:
        return subs
    for tok in re.split('[|,;\\s]+', cell.strip()):
        t = re.sub('\\(.*?\\)', '', tok).strip()
        if not t or t in ('-', 'N', 'NA', 'None'):
            continue
        if FAMILY_RE.match(t) and '_' in t:
            if t not in subs:
                subs.append(t)
    return subs

def pick_col(header, *names):
    """Return index of first header cell whose normalised name matches any of *names."""
    norm = [re.sub('[^a-z0-9]', '', h.lower()) for h in header]
    for want in names:
        w = re.sub('[^a-z0-9]', '', want.lower())
        for (i, h) in enumerate(norm):
            if h == w:
                return i
    for want in names:
        w = re.sub('[^a-z0-9]', '', want.lower())
        for (i, h) in enumerate(norm):
            if w and w in h:
                return i
    return None

def gene_from_cds(feature):
    """Return BRAKER gene id (g*) for a CDS feature: prefer gene/locus_tag/protein_id,
    then strip transcript suffix .tN -> gene."""
    q = feature.qualifiers
    val = None
    for key in ('gene', 'locus_tag', 'protein_id'):
        if key in q and q[key]:
            val = q[key][0]
            break
    if val is None:
        return None
    return val.split('.')[0]

def parse_genbank(outdir):
    """Yield (gene_id, cluster_id, cluster_type) from *.final.gbk files."""
    from Bio import SeqIO
    rows = []
    gbks = sorted(glob.glob(os.path.join(outdir, '*.final.gbk')))
    if not gbks:
        return (rows, 0)
    n_clusters = 0
    for gbk in gbks:
        for rec in SeqIO.parse(gbk, 'genbank'):
            clusters = [f for f in rec.features if f.type == 'cluster']
            cds = [f for f in rec.features if f.type == 'CDS']
            for clf in clusters:
                n_clusters += 1
                cnum = None
                for note in clf.qualifiers.get('note', []):
                    if 'Cluster number:' in note:
                        cnum = note.split('Cluster number:')[1].strip()
                        break
                if cnum is None:
                    cnum = str(n_clusters)
                ctype = clf.qualifiers.get('product', ['unknown'])[0]
                rec_id = rec.id.split('.')[0]
                cluster_id = '%s_c%s' % (rec_id, cnum)
                cstart = int(clf.location.start)
                cend = int(clf.location.end)
                for c in cds:
                    cs = int(c.location.start)
                    ce = int(c.location.end)
                    if cs < cend and ce > cstart:
                        gid = gene_from_cds(c)
                        if gid:
                            rows.append((gid, cluster_id, ctype))
    return (rows, n_clusters)


def table(file, header, rows):
    file.parent.mkdir(parents=True, exist_ok=True)
    with file.open('w', newline='') as out:
        writer = csv.writer(out, delimiter='\t', lineterminator='\n')
        writer.writerow(header)
        writer.writerows(rows)


def interpro(a):
    annotations = defaultdict(lambda: [set() for _ in range(4)])
    for file in a.inputs:
        with file.open() as src:
            for line in src:
                if not line.strip() or line.startswith('#'):
                    continue
                f = line.rstrip('\n').split('\t')
                if len(f) < 11:
                    raise ValueError(f'Invalid InterProScan row in {file}')
                fields = annotations[f[0]]
                for index, column in enumerate((11, 12, 13, 14)):
                    if len(f) <= column or f[column] in ('', '-'):
                        continue
                    if column == 13:
                        fields[index].update(re.findall(r'GO:\d+', f[column]))
                    else:
                        fields[index].update(f[column].split('|'))
    ids = [h.split()[0] for h, _ in fasta(a.proteins)]
    rows = ([pid] + [';'.join(sorted(values)) for values in annotations[pid]] for pid in ids)
    table(a.out, ['protein_id', 'ips_IPR', 'ips_IPR_desc', 'ips_GO', 'ips_pathways'], rows)


def eggnog(a):
    columns = [('Preferred_name', 'eggnog_preferred_name'), ('Description', 'eggnog_description'),
               ('GOs', 'eggnog_GOs'), ('KEGG_ko', 'eggnog_KEGG_ko'), ('COG_category', 'eggnog_COG_cat')]
    rows = {}
    for file in a.inputs:
        header = None
        with file.open() as src:
            for line in src:
                if line.startswith('##') or not line.strip():
                    continue
                fields = line.rstrip('\n').split('\t')
                if line.startswith('#'):
                    if fields[0].lstrip('#').lower() in ('query', 'query_name'):
                        header = [x.lstrip('#').lower() for x in fields]
                    continue
                if header is None:
                    raise ValueError(f'No eggNOG header in {file}')
                row = dict(zip(header, fields))
                values = [row.get(name.lower(), '') for name, _ in columns]
                rows[fields[0]] = ['' if value in ('-', 'NA', 'None') else value for value in values]
    table(a.out, ['protein_id'] + [dest for _, dest in columns], ([pid] + vals for pid, vals in rows.items()))


def itak(a):
    rows, seen = [], set()
    for file, kind in ((a.tf, 'TF'), (a.kinases, 'PK')):
        for line in file.read_text().splitlines():
            if not line.strip():
                continue
            f = line.split('\t')
            if kind == 'TF':
                if len(f) < 4:
                    raise ValueError(f'Expected four TF columns: {line}')
                family = f[3].strip() or f[1].strip()
            else:
                family = (f[1].strip() if len(f) > 1 else '') or 'Group-other'
            if (f[0], kind) not in seen:
                rows.append((f[0], kind, family))
                seen.add((f[0], kind))
    table(a.out, ['protein_id', 'type', 'family'], sorted(rows))


def dbcan(a):
    rows = []
    with a.overview.open() as src:
        reader = csv.reader(src, delimiter='\t')
        header = next(reader)
        gene = pick_col(header, 'Gene ID', 'GeneID', 'Gene_ID', 'gene')
        hmm = pick_col(header, 'dbCAN_hmm', 'HMMER', 'dbCAN')
        sub = pick_col(header, 'dbCAN_sub', 'dbCANsub', 'dbCAN-sub')
        dia = pick_col(header, 'DIAMOND', 'CAZy')
        count = pick_col(header, '#ofTools', 'Nofools', 'ofTools', 'NumberOfTools')
        for row in reader:
            if not row:
                continue
            def cell(index):
                return row[index] if index is not None and index < len(row) else ''
            families = [families_from_cell(cell(i)) for i in (hmm, sub, dia)]
            votes = Counter(f for group in families for f in set(group))
            if not votes:
                continue
            best = sorted(votes, key=lambda f: (-votes[f], f))[0]
            subs = subfamilies_from_cell(cell(sub))
            chosen = [s for s in subs if s.upper().startswith(best.upper())] or subs
            n = sum(bool(f) for f in families)
            match = re.search(r'\d+', cell(count))
            if match:
                n = int(match.group())
            rows.append((row[gene if gene is not None else 0], best, '|'.join(chosen) or '-', n))
    table(a.out, ['protein_id', 'CAZy_family', 'CAZy_subfamily', 'n_tools'], rows)


def kofam(a):
    best = {}
    for file in a.inputs:
        with file.open() as src:
            for f in csv.reader(src, delimiter='\t'):
                if len(f) < 7 or f[0] != '*':
                    continue
                if f[1] not in best or float(f[4]) > float(best[f[1]][3]):
                    best[f[1]] = [f[1], f[2], f[3], f[4], f[6]]
    table(a.out, ['protein_id', 'KO', 'KO_threshold', 'KO_score', 'KO_definition'],
          (best[key] for key in sorted(best)))


def targetp(a):
    predictions = {}
    for file in a.inputs:
        for line in file.read_text().splitlines():
            if line.startswith('#') or not line.strip():
                continue
            f = line.split('\t')
            if len(f) < 2:
                raise ValueError(f'Invalid TargetP row: {line}')
            predictions[f[0].split()[0]] = f[1].strip()
    table(a.out, ['protein_id', 'prediction'], sorted(predictions.items()))


def descriptions(a):
    best = parse_blast(a.blast_dir)
    rows = []
    for header, _ in fasta(a.proteins):
        pid = header.split()[0]
        label, acc, value, tag, title = best_hit_for(pid, best)
        description = clean_stitle(title, tag) if title else ''
        rows.append((pid, description or 'Unknown protein', label or 'NA', acc or 'NA', value or 'NA'))
    table(a.out, ['protein_id', 'description', 'source_db', 'best_hit_acc', 'best_hit_evalue'], rows)


def plantismash(a):
    if not list(a.results.glob('*.final.gbk')):
        raise ValueError('No plantiSMASH final GenBank records')
    rows, _ = parse_genbank(str(a.results))
    genes = set()
    with a.gff3.open() as src:
        for line in src:
            if not line.startswith('#') and '\tgene\t' in line:
                genes.add(attributes(line.rstrip().split('\t')[8])['ID'])
    clean = set()
    for gene, cluster, kind in rows:
        # plantiSMASH appends _N to some duplicate CDS labels.
        if gene not in genes:
            candidate = re.sub(r'_\d+$', '', gene)
            if candidate not in genes:
                raise ValueError(f'plantiSMASH gene absent from annotation: {gene}')
            gene = candidate
        clean.add((gene, cluster, kind))
    table(a.out, ['gene_id', 'cluster_id', 'cluster_type'], sorted(clean, key=lambda r: (r[1], r[0])))


def targeted(a):
    rubber = {}
    for line in a.rubber.read_text().splitlines():
        if not line.strip():
            continue
        f = line.split('\t')
        family = f[1].split('__')[0]
        score = float(f[5])
        if f[0] not in rubber or score > rubber[f[0]][0]:
            rubber[f[0]] = (score, family if family in {'CPT', 'CPTL', 'REF', 'SRPP', 'HRT'} else '')
    rubber = {pid: value[1] for pid, value in rubber.items() if value[1]}
    def hits(file):
        return {line.split()[0] for line in file.read_text().splitlines() if line.strip() and not line.startswith('#')}
    p450, tps = hits(a.p450), hits(a.tps)
    rows = [(pid, rubber.get(pid, ''), 'yes' if pid in p450 else '', 'yes' if pid in tps else '')
            for pid in sorted(set(rubber) | p450 | tps)]
    table(a.out, ['protein_id', 'rubber_family', 'is_P450', 'is_TPS'], rows)


def merge(a):
    ids = [h.split()[0] for h, _ in fasta(a.proteins)]
    if len(ids) != len(set(ids)):
        raise ValueError('Duplicate protein IDs')
    idset, gene2tx = set(ids), defaultdict(list)
    if a.gff3:
        with a.gff3.open() as src:
            for line in src:
                if line.startswith('#') or not line.strip():
                    continue
                f = line.rstrip().split('\t')
                if f[2] in ('mRNA', 'transcript'):
                    attrs = attributes(f[8])
                    if attrs['ID'] in idset:
                        for gene in attrs['Parent'].split(','):
                            gene2tx[gene].append(attrs['ID'])
    else:
        for pid in ids:
            gene2tx[re.sub(r'\.t\d+$', '', pid)].append(pid)
    data, columns, labels = defaultdict(dict), [], set()
    for file in sorted(a.inputs):
        label = file.name.removesuffix('.tsv').removesuffix('.perGene')
        if label in labels:
            raise ValueError(f'Duplicate track name: {label}')
        labels.add(label)
        with file.open() as src:
            reader = csv.reader(src, delimiter='\t')
            header = next(reader)
            cols = [label + '.' + c for c in header[1:]]
            columns.extend(cols)
            for row in reader:
                if not row:
                    continue
                key, values = row[0], row[1:]
                targets = [key] if key in idset else gene2tx.get(key, [])
                if not targets:
                    raise ValueError(f'{file.name}: ID absent from proteome/GFF: {key}')
                for tx in targets:
                    for col, value in zip(cols, values):
                        if value and value != '-':
                            if label == 'plantismash':
                                previous = set(data[tx].get(col, '').split(';')) - {''}
                                data[tx][col] = ';'.join(sorted(previous | set(value.split(';'))))
                            else:
                                data[tx][col] = value
    table(a.out, ['protein_id'] + columns, ([pid] + [data[pid].get(c, '') for c in columns] for pid in ids))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    for name in ('interpro', 'eggnog', 'itak', 'dbcan', 'kofam', 'targetp', 'descriptions', 'plantismash', 'targeted', 'merge'):
        q = sub.add_parser(name)
        q.add_argument('--out', type=path, required=True, help='output TSV')
        if name in ('interpro', 'eggnog', 'kofam', 'targetp', 'merge'):
            q.add_argument('--inputs', type=path, nargs='+', required=True)
        if name in ('interpro', 'descriptions', 'merge'):
            q.add_argument('--proteins', type=path, required=True)
        if name == 'itak':
            q.add_argument('--tf', type=path, required=True)
            q.add_argument('--kinases', type=path, required=True)
        if name == 'dbcan':
            q.add_argument('--overview', type=path, required=True)
        if name == 'descriptions':
            q.add_argument('--blast-dir', type=path, required=True)
        if name == 'plantismash':
            q.add_argument('--results', type=path, required=True)
            q.add_argument('--gff3', type=path, required=True)
        if name == 'targeted':
            for flag in ('rubber', 'p450', 'tps'):
                q.add_argument('--' + flag, type=path, required=True)
        if name == 'merge':
            q.add_argument('--gff3', type=path, help='gene/transcript mapping; default IDs use .tN suffix')
        q.set_defaults(func=globals()[name])
    a = p.parse_args()
    a.func(a)


if __name__ == '__main__':
    main()
