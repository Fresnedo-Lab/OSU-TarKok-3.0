#!/usr/bin/env python3
"""RNA-seq preparation, BRAKER3 prediction and UTR annotation."""
import argparse
import re
import subprocess
from _common import Runner, execution_options, path, fasta


def sample_pairs(file):
    """Read a headerless TSV: sample, R1, R2; paths relative to the TSV."""
    pairs, seen = [], set()
    for line in file.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        name, r1, r2 = line.split("\t")
        if name in seen or not re.fullmatch(r"[A-Za-z0-9_.-]+", name) or name in (".", ".."):
            raise ValueError(f"Invalid or duplicate sample: {name}")
        seen.add(name)
        pairs.append((name, (file.parent / r1).resolve(), (file.parent / r2).resolve()))
    if not pairs:
        raise ValueError("Empty sample table")
    return pairs


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    for name in ("trim", "map", "braker", "utr"):
        q = sub.add_parser(name)
        execution_options(q)
        if name in ("trim", "map"):
            q.add_argument("--samples", type=path, required=True, help="headerless TSV: sample, R1, R2")
        if name in ("map", "braker"):
            q.add_argument("--genome", type=path, required=True)
        if name in ("braker", "utr"):
            q.add_argument("--bam", type=path, required=True)
        if name == "braker":
            q.add_argument("--proteins", type=path, required=True, help="protein evidence FASTA")
            q.add_argument("--species", required=True, help="new AUGUSTUS species name")
        if name == "utr":
            q.add_argument("--gtf", type=path, required=True)
            q.add_argument("--stringtie2utr", type=path, required=True, help="BRAKER stringtie2utr.py")
            q.add_argument("--python", default="python3")
    a = p.parse_args()
    pairs = sample_pairs(a.samples) if a.command in ("trim", "map") else []
    r = Runner(a)
    if a.command == "trim":
        for name, r1, r2 in pairs:
            r.run(["trim_galore", "--paired", "--cores", a.threads, "-o", a.out / name, r1, r2])
    elif a.command == "map":
        index = a.out / "genome"
        r.run(["hisat2-build", "-p", a.threads, a.genome, index])
        bams = []
        for name, r1, r2 in pairs:
            bam = a.out / f"sample.{name}.bam"
            r.pipe([["hisat2", "-p", a.threads, "--dta", "-x", index, "-1", r1, "-2", r2],
                    ["samtools", "sort", "-o", bam, "-"]])
            bams.append(bam)
        r.run(["samtools", "merge", "-@", a.threads, a.out / "rnaseq.bam", *bams])
        r.run(["samtools", "index", a.out / "rnaseq.bam"])
    elif a.command == "braker":
        if not r.dry:
            genome = {h.split()[0]: len(s) for h, s in fasta(a.genome)}
            header = subprocess.check_output(["samtools", "view", "-H", a.bam], text=True)
            bam = {}
            for line in header.splitlines():
                if line.startswith("@SQ\t"):
                    fields = dict(x.split(":", 1) for x in line.split("\t")[1:])
                    bam[fields["SN"]] = int(fields["LN"])
            if genome != bam:
                raise ValueError("FASTA and BAM sequence names/lengths differ")
        r.run(["braker.pl", f"--species={a.species}", f"--genome={a.genome}",
               f"--bam={a.bam}", f"--prot_seq={a.proteins}", "--softmasking", "--gff3",
               f"--threads={a.threads}", f"--workingdir={a.out}"])
    else:
        r.run(["stringtie", a.bam, "-p", a.threads, "-o", a.out / "stringtie.gtf"])
        r.run([a.python, a.stringtie2utr, "-g", a.gtf, "-s", a.out / "stringtie.gtf",
               "-o", a.out / "genes.utr.gtf"])
        r.run(["agat_convert_sp_gxf2gxf.pl", "-g", a.out / "genes.utr.gtf", "-o", a.out / "genes.utr.gff3"])


if __name__ == "__main__":
    main()
