#!/usr/bin/env python3
"""Orient a plastome and measure its inverted repeats and reference gene recovery."""
import argparse
import os
import subprocess
import sys
from collections import OrderedDict, defaultdict
from Bio import SeqIO
from Bio.Seq import Seq


def run_paf(target, query, mm2, preset='asm5'):
    out = subprocess.run([mm2, '-cx', preset, '-t', '4', target, query], capture_output=True, text=True, check=True)
    recs = []
    for line in out.stdout.splitlines():
        f = line.split('\t')
        if len(f) < 12:
            continue
        recs.append(dict(q=f[0], qlen=int(f[1]), qs=int(f[2]), qe=int(f[3]), strand=f[4], t=f[5], tlen=int(f[6]), ts=int(f[7]), te=int(f[8]), nmatch=int(f[9]), alen=int(f[10]), mapq=int(f[11])))
    return recs

def write_fasta(name, seq, path, width=70):
    with open(path, 'w') as fh:
        fh.write('>' + name + '\n')
        for i in range(0, len(seq), width):
            fh.write(seq[i:i + width] + '\n')

def canonicalize(contig_fa, ref_fa, out_fa, mm2):
    rec = list(SeqIO.parse(contig_fa, 'fasta'))[0]
    seq = str(rec.seq).upper()
    L = len(seq)
    name = rec.id
    recs = run_paf(ref_fa, contig_fa, mm2)
    if not recs:
        raise ValueError('No alignment to plastome reference')
    plus = sum((r['alen'] for r in recs if r['strand'] == '+'))
    minus = sum((r['alen'] for r in recs if r['strand'] == '-'))
    flip = minus > plus
    if flip:
        seq = str(Seq(seq).reverse_complement())
        tmp = out_fa + '.flip.tmp.fa'
        write_fasta('tmp', seq, tmp)
        recs = run_paf(ref_fa, tmp, mm2)
        os.remove(tmp)
    recs_p = [r for r in recs if r['strand'] == '+' and r['alen'] >= 1000]
    if not recs_p:
        raise ValueError('No forward alignment block of at least 1 kb')
    r = min(recs_p, key=lambda x: x['ts'])
    q0 = (r['qs'] - r['ts']) % L
    seq2 = seq[q0:] + seq[:q0]
    write_fasta(name, seq2, out_fa)
    print(f"canonicalize: L={L} strand_flip={flip} rotate_by={q0} anchor(ref_ts={r['ts']} contig_qs={r['qs']} alen={r['alen']})")

def gc(s):
    s = s.upper()
    n = len(s)
    return 100.0 * (s.count('G') + s.count('C')) / n if n else 0.0

def structure(fa, makeblastdb, blastn, out_tsv):
    rec = list(SeqIO.parse(fa, 'fasta'))[0]
    seq = str(rec.seq).upper()
    L = len(seq)
    db = out_tsv + '.selfdb'
    subprocess.run([makeblastdb, '-in', fa, '-dbtype', 'nucl', '-out', db], capture_output=True, text=True, check=True)
    bl = subprocess.run([blastn, '-query', fa, '-db', db, '-evalue', '1e-20', '-word_size', '50', '-outfmt', '6 qstart qend sstart send length pident sstrand'], capture_output=True, text=True, check=True)
    for ext in ('.selfdb.ndb', '.selfdb.nhr', '.selfdb.nin', '.selfdb.njs', '.selfdb.not', '.selfdb.nsq', '.selfdb.ntf', '.selfdb.nto'):
        try:
            os.remove(out_tsv + ext)
        except OSError:
            pass
    ir_hsps = []
    for line in bl.stdout.splitlines():
        f = line.split('\t')
        if len(f) < 7:
            continue
        (qs, qe, ss, se, length) = map(int, f[:5])
        pident = float(f[5])
        strand = f[6]
        if strand == 'minus' and length >= 1000:
            ir_hsps.append((length, pident, qs, qe, ss, se))
    lines = []
    lines.append(f'metric\tvalue')
    lines.append(f'total_length_bp\t{L}')
    lines.append(f'GC_percent\t{gc(seq):.2f}')
    if not ir_hsps:
        lines.append('IR_detected\tNO')
        lines.append('note\tno_minus_strand_repeat>=1kb;_possible_IR_loss_or_incomplete')
        with open(out_tsv, 'w') as fh:
            fh.write('\n'.join(lines) + '\n')
        print('WARN: no IR detected')
        return
    ir_hsps.sort(reverse=True)
    (length, pident, qs, qe, ss, se) = ir_hsps[0]
    (a1, a2) = (min(qs, qe), max(qs, qe))
    (b1, b2) = (min(ss, se), max(ss, se))
    (ira, irb) = ((a1, a2), (b1, b2))
    if a1 > b1:
        (ira, irb) = ((b1, b2), (a1, a2))
    (ia1, ia2) = ira
    (ib1, ib2) = irb
    ir_len = ia2 - ia1 + 1
    ir_len2 = ib2 - ib1 + 1
    sc_mid = ib1 - ia2 - 1
    sc_wrap = L - ib2 + (ia1 - 1)
    lsc = max(sc_mid, sc_wrap)
    ssc = min(sc_mid, sc_wrap)
    mid_seq = seq[ia2:ib1 - 1] if ib1 - 1 > ia2 else ''
    ir_seq = seq[ia1 - 1:ia2]
    lines.append(f'IR_detected\tYES')
    lines.append(f'IR_identity_percent\t{pident:.2f}')
    lines.append(f'IRa_start\t{ia1}')
    lines.append(f'IRa_end\t{ia2}')
    lines.append(f'IRb_start\t{ib1}')
    lines.append(f'IRb_end\t{ib2}')
    lines.append(f'IR_length_bp\t{ir_len}')
    lines.append(f'IR_length_bp_copy2\t{ir_len2}')
    lines.append(f'LSC_length_bp\t{lsc}')
    lines.append(f'SSC_length_bp\t{ssc}')
    lines.append(f'IR_GC_percent\t{gc(ir_seq):.2f}')
    lines.append(f'midSC_GC_percent\t{gc(mid_seq):.2f}')
    lines.append(f'sum_LSC+SSC+2xIR\t{lsc + ssc + ir_len + ir_len2}')
    lines.append(f"quadripartite_closes\t{('YES' if abs(lsc + ssc + ir_len + ir_len2 - L) <= 200 else 'CHECK')}")
    with open(out_tsv, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    print('\n'.join(lines))

def merged_cov(intervals, qlen):
    """Fraction of the query [1..qlen] covered by the union of intervals."""
    if not qlen:
        return 0.0
    ivs = sorted(((min(a, b), max(a, b)) for (a, b) in intervals))
    cov = 0
    cs = ce = None
    for (a, b) in ivs:
        if cs is None:
            (cs, ce) = (a, b)
        elif a <= ce + 1:
            ce = max(ce, b)
        else:
            cov += ce - cs + 1
            (cs, ce) = (a, b)
    if cs is not None:
        cov += ce - cs + 1
    return cov / qlen

def inventory(plastome_fa, gb_file, makeblastdb, blastn, out_prefix):
    rec = SeqIO.read(gb_file, 'genbank')
    full = rec.seq
    genes = OrderedDict()
    for feat in rec.features:
        if feat.type not in ('CDS', 'tRNA', 'rRNA'):
            continue
        q = feat.qualifiers
        name = (q.get('gene') or q.get('locus_tag') or q.get('product') or ['?'])[0]
        name = name.strip().replace(' ', '_')
        try:
            s = str(feat.extract(full)).upper()
        except Exception:
            continue
        if len(s) < 20:
            continue
        if name not in genes:
            genes[name] = [feat.type, s, 1]
        else:
            genes[name][2] += 1
            if feat.type == 'CDS':
                genes[name][0] = 'CDS'
            if len(s) > len(genes[name][1]):
                genes[name][1] = s
    genefa = out_prefix + '.ref_genes.fa'
    with open(genefa, 'w') as fh:
        for (name, (typ, s, cnt)) in genes.items():
            fh.write(f'>{name}|{typ}\n{s}\n')
    db = out_prefix + '_blastdb'
    r1 = subprocess.run([makeblastdb, '-in', plastome_fa, '-dbtype', 'nucl', '-out', db], capture_output=True, text=True, check=True)
    bl = subprocess.run([blastn, '-task', 'blastn', '-query', genefa, '-db', db, '-evalue', '1e-5', '-max_target_seqs', '10', '-outfmt', '6 qseqid qlen qstart qend length pident sstart send bitscore'], capture_output=True, text=True, check=True)
    if r1.returncode != 0 or bl.returncode != 0:
        sys.stderr.write('BLAST/makeblastdb error:\n' + r1.stderr + '\n' + bl.stderr + '\n')
    hsps = defaultdict(list)
    qlen_map = {}
    for line in bl.stdout.splitlines():
        f = line.split('\t')
        if len(f) < 9:
            continue
        qid = f[0]
        qlen = int(f[1])
        qs = int(f[2])
        qe = int(f[3])
        alen = int(f[4])
        pid = float(f[5])
        ss = int(f[6])
        se = int(f[7])
        bits = float(f[8])
        qlen_map[qid] = qlen
        hsps[qid].append((qs, qe, pid, ss, se, bits, alen))
    best = {}
    for (qid, lst) in hsps.items():
        qlen = qlen_map[qid]
        cov = merged_cov([(x[0], x[1]) for x in lst], qlen)
        top = max(lst, key=lambda x: x[6])
        bestpid = max((x[2] for x in lst))
        best[qid] = (bestpid, cov, min(top[3], top[4]), max(top[3], top[4]))
    for e in ('.ndb', '.nhr', '.nin', '.njs', '.not', '.nsq', '.ntf', '.nto'):
        try:
            os.remove(db + e)
        except OSError:
            pass
    rows = []
    cls_total = {'CDS': 0, 'tRNA': 0, 'rRNA': 0}
    cls_rec = {'CDS': 0, 'tRNA': 0, 'rRNA': 0}
    missing = []
    for (name, (typ, s, cnt)) in genes.items():
        qid = f'{name}|{typ}'
        cls_total[typ] = cls_total.get(typ, 0) + 1
        if qid in best:
            (pid, qcov, hs, he) = best[qid]
            recovered = pid >= 70.0 and qcov >= 0.6
        else:
            (pid, qcov, hs, he) = (0.0, 0.0, 0, 0)
            recovered = False
        if recovered:
            cls_rec[typ] = cls_rec.get(typ, 0) + 1
        else:
            missing.append(f'{name}({typ})')
        rows.append((name, typ, len(s), cnt, f'{pid:.1f}', f'{qcov:.2f}', hs, he, 'Y' if recovered else 'N'))
    with open(out_prefix + '.tsv', 'w') as fh:
        fh.write('gene\tclass\tref_len\tref_copies\tbest_pident\tmerged_qcov\thit_start\thit_end\trecovered\n')
        for r in rows:
            fh.write('\t'.join((str(x) for x in r)) + '\n')
    total = sum(cls_total.values())
    rec_total = sum(cls_rec.values())
    with open(out_prefix + '.summary.txt', 'w') as fh:
        fh.write(f'reference\t{os.path.basename(gb_file)}\n')
        fh.write(f'unique_genes_reference\t{total}\n')
        fh.write(f'genes_recovered\t{rec_total}\n')
        for c in ('CDS', 'tRNA', 'rRNA'):
            fh.write(f'{c}_reference\t{cls_total.get(c, 0)}\n')
            fh.write(f'{c}_recovered\t{cls_rec.get(c, 0)}\n')
        fh.write(f"missing_genes\t{(';'.join(missing) if missing else 'NONE')}\n")
    print(open(out_prefix + '.summary.txt').read())

def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    for name in ('orient', 'structure', 'inventory'):
        q = sub.add_parser(name)
        q.add_argument('--fasta', required=True)
        q.add_argument('--out', required=True)
        if name == 'orient':
            q.add_argument('--reference', required=True)
            q.add_argument('--minimap2', default='minimap2')
        else:
            q.add_argument('--makeblastdb', default='makeblastdb')
            q.add_argument('--blastn', default='blastn')
        if name == 'inventory':
            q.add_argument('--genbank', required=True)
    a = p.parse_args()
    if a.command == 'orient':
        canonicalize(a.fasta, a.reference, a.out, a.minimap2)
    elif a.command == 'structure':
        structure(a.fasta, a.makeblastdb, a.blastn, a.out)
    else:
        inventory(a.fasta, a.genbank, a.makeblastdb, a.blastn, a.out)

if __name__ == '__main__':
    main()
