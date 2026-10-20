#!/usr/bin/env python3
"""Run protein annotation and plant biosynthetic cluster detection."""
import argparse
from _common import Runner, execution_options, path, positive, fasta, write_fasta


def prepare(a):
    proteins, seen = [], set()
    for header, sequence in fasta(a.proteins):
        name = header.split()[0]
        sequence = sequence.removesuffix("*")
        if name in seen or not sequence or "*" in sequence:
            raise ValueError(f"Duplicate ID, empty sequence or internal stop: {name}")
        seen.add(name)
        proteins.append((name, sequence))
    if not proteins:
        raise ValueError("Empty protein FASTA")
    a.out.mkdir(parents=True, exist_ok=False)
    with (a.out / "proteins.fa").open("w") as out:
        for name, sequence in proteins:
            write_fasta(out, name, sequence)
    size = (len(proteins) + a.chunks - 1) // a.chunks
    if a.chunks > 1:
        for n, start in enumerate(range(0, len(proteins), size), 1):
            with (a.out / f"part_{n:03d}.fa").open("w") as out:
                for name, sequence in proteins[start:start + size]:
                    write_fasta(out, name, sequence)
    print(f"{len(proteins)} proteins")


def annotate(a):
    r = Runner(a)
    if not r.dry and a.command in ("interpro", "eggnog", "kofam", "targetp"):
        (a.out / "tmp").mkdir()
    prefix = a.out / "proteins"
    if a.command == "interpro":
        r.run(["interproscan.sh", "-i", a.proteins, "-f", "TSV,GFF3", "--goterms", "--pathways",
               "--disable-precalc", "--cpu", a.threads, "--tempdir", a.out / "tmp", "-b", prefix])
    elif a.command == "eggnog":
        r.run(["emapper.py", "-i", a.proteins, "--itype", "proteins", "-m", "diamond",
               "--data_dir", a.db, "--cpu", a.threads, "--temp_dir", a.out / "tmp",
               "-o", "proteins", "--output_dir", a.out])
    elif a.command == "itak":
        r.run(["itak", "-p", a.threads, "-o", a.out / "results", a.proteins], env={"ITAK_DB_DIR": str(a.db)})
    elif a.command == "dbcan":
        r.run(["run_dbcan", "CAZyme_annotation", "--input_raw_data", a.proteins, "--mode", "protein",
               "--output_dir", a.out, "--db_dir", a.db, "--threads", a.threads])
    elif a.command == "kofam":
        r.run(["exec_annotation", "--cpu", a.threads, "-p", a.profiles, "-k", a.ko_list,
               "--tmp-dir", a.out / "tmp", "-f", "detail-tsv", "-o", a.out / "kofam.tsv", a.proteins])
    elif a.command == "targetp":
        r.run([a.install / "bin/targetp", "-fasta", a.proteins, "-org", "pl", "-format", "short",
               "-tmp", a.out / "tmp", "-prefix", prefix], cwd=a.install)
    elif a.command == "plantismash":
        gff = a.out / "genes.gff3"
        if not r.dry:
            seqids = {h.split()[0] for h, _ in fasta(a.genome)}
            with a.gff3.open() as src, gff.open("w") as dst:
                for line in src:
                    if line.startswith("##FASTA"):
                        break
                    if line.startswith("#") or not line.strip():
                        continue
                    f = line.rstrip().split("\t")
                    if f[2] in {"gene", "mRNA", "CDS"}:
                        if f[0] not in seqids:
                            raise ValueError(f"GFF sequence missing from FASTA: {f[0]}")
                        dst.write(line)
        r.run([a.python, a.install / "run_antismash.py", "--cpus", a.threads,
               "--taxon", "plants", "--input-type", "nucl", "--eukaryotic", "--genefinding", "none",
               "--gff3", gff, "--cdhit", "--limit", "-1", "--outputfolder", a.out / "results", a.genome], cwd=a.install)
    elif a.command == "descriptions":
        for label in ("sprot", "arabidopsis", "lactuca"):
            r.run(["diamond", "blastp", "--query", a.proteins, "--db", getattr(a, label),
                   "--out", a.out / f"{label}.blastp.tsv", "--max-target-seqs", "10", "--evalue", "1e-5",
                   "--threads", a.threads, "--outfmt", "6", "qseqid", "sseqid", "pident", "length",
                   "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle"])
    elif a.command == "targeted":
        r.run(["diamond", "makedb", "--in", a.rubber_refs, "-d", a.out / "rubber", "--threads", a.threads])
        r.run(["diamond", "blastp", "--query", a.proteins, "--db", a.out / "rubber", "--out", a.out / "rubber.tsv",
               "--evalue", "1e-10", "--max-target-seqs", "5", "--threads", a.threads,
               "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "stitle"])
        for label, models in (("p450", [a.p450_model]), ("tps", [a.tps_n_model, a.tps_c_model])):
            pieces = []
            for i, model in enumerate(models):
                piece = a.out / f"{label}.{i}.hmm"
                r.run(["hmmfetch", a.pfam, model], stdout=piece)
                pieces.append(piece)
            hmm = a.out / f"{label}.hmm"
            if not r.dry:
                hmm.write_text("".join(piece.read_text() for piece in pieces))
            r.run(["hmmsearch", "--cpu", a.threads, "--cut_ga", "--domtblout", a.out / f"{label}.domtbl",
                   "--tblout", a.out / f"{label}.tbl", "-o", a.out / f"{label}.txt", hmm, a.proteins])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    q = sub.add_parser("prepare", help="remove terminal stops; optionally split the proteome")
    q.add_argument("--proteins", type=path, required=True)
    q.add_argument("--out", type=path, required=True)
    q.add_argument("--chunks", type=positive, default=1)
    q.set_defaults(func=prepare)
    for name in ("interpro", "eggnog", "itak", "dbcan", "kofam", "targetp", "plantismash", "descriptions", "targeted"):
        q = sub.add_parser(name)
        execution_options(q)
        q.set_defaults(func=annotate)
        if name != "plantismash":
            q.add_argument("--proteins", type=path, required=True)
        if name in ("eggnog", "itak", "dbcan"):
            q.add_argument("--db", type=path, required=True)
        if name == "kofam":
            q.add_argument("--profiles", type=path, required=True)
            q.add_argument("--ko-list", type=path, required=True)
        if name in ("targetp", "plantismash"):
            q.add_argument("--install", type=path, required=True)
        if name == "plantismash":
            q.add_argument("--genome", type=path, required=True)
            q.add_argument("--gff3", type=path, required=True)
            q.add_argument("--python", default="python", help="Python executable for plantiSMASH")
        if name == "descriptions":
            for label in ("sprot", "arabidopsis", "lactuca"):
                q.add_argument(f"--{label}", type=path, required=True, help="DIAMOND database")
        if name == "targeted":
            q.add_argument("--rubber-refs", type=path, required=True, help="FASTA IDs prefixed FAMILY__ (CPT/CPTL/REF/SRPP/HRT)")
            q.add_argument("--pfam", type=path, required=True)
            q.add_argument("--p450-model", default="PF00067.28")
            q.add_argument("--tps-n-model", default="PF01397.28")
            q.add_argument("--tps-c-model", default="PF03936.22")
    a = p.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
