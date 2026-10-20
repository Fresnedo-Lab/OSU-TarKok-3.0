#!/usr/bin/env python3
"""HiFi assembly, reference-guided scaffolding and Racon polishing."""
import argparse
from pathlib import Path
import shutil
from _common import Runner, execution_options, path, positive


def assemble(a):
    r = Runner(a)
    if a.command == "hifiasm":
        if bool(a.hic_r1) != bool(a.hic_r2):
            raise ValueError("Supply both Hi-C read files")
        cmd = ["hifiasm", "-o", a.out / "assembly", "-t", a.threads]
        if a.hic_r1:
            cmd += ["--h1", a.hic_r1, "--h2", a.hic_r2]
        r.run(cmd + [a.reads])
        suffixes = ["hic.hap1.p_ctg", "hic.hap2.p_ctg"] if a.hic_r1 else ["bp.p_ctg"]
        if not r.dry:
            for suffix in suffixes:
                with (a.out / f"assembly.{suffix}.gfa").open() as src, (a.out / f"{suffix}.fa").open("w") as dst:
                    for line in src:
                        f = line.rstrip().split("\t")
                        if f[0] == "S":
                            if f[2] == "*":
                                raise ValueError("GFA segment has no sequence")
                            dst.write(f">{f[1]}\n{f[2]}\n")
    elif a.command == "verkko":
        r.run(["verkko", "-d", a.out / "assembly", "--hifi", a.reads, "--no-nano",
               "--local-cpus", a.threads, "--local-memory", a.memory_gb])
    else:
        r.run(["flye", "--pacbio-hifi", a.reads, "--out-dir", a.out / "assembly", "--threads", a.threads])


def ragtag(a):
    r = Runner(a)
    command = ['ragtag.py', a.action, '-t', a.threads, '-w', '-o', a.out / 'results']
    if a.aligner:
        command += ['--aligner', a.aligner]
    r.run(command + [a.reference, a.query])


def scaffold(a):
    r = Runner(a)
    agps = []
    for label, reference, preset in [("framework", a.framework, "asm5"), ("reference", a.reference, "asm10")]:
        out = a.out / label
        r.run(["ragtag.py", "scaffold", reference, a.backbone, "--aligner", "minimap2",
               "--mm2-params", f"-x {preset} -t {a.threads}", "-u", "-w", "-o", out])
        agps.append(out / "ragtag.scaffold.agp")
    weights = a.out / "weights.csv"
    r.write(weights, f"{agps[0]},2\n{agps[1]},1\n")
    r.run(["ragtag.py", "merge", a.backbone, *agps, "-f", weights,
           "-l", "1000", "-u", "-w", "-o", a.out / "merged"])
    current = a.out / "merged/ragtag.merge.fasta"
    for label, donor in [("hifiasm_patch", a.donor), ("framework_patch", a.framework)]:
        r.run(["ragtag.py", "patch", current, donor, "--aligner", "minimap2",
               "--mm2-params", f"-x asm5 -t {a.threads}", "-u", "-w", "-o", a.out / label])
        current = a.out / label / "ragtag.patch.fasta"
    if not r.dry:
        shutil.copyfile(current, a.out / "scaffolded.fa")


def polish(a):
    r = Runner(a)
    current = a.genome
    for i in range(1, a.rounds + 1):
        paf, result = a.out / f"round{i}.paf", a.out / f"round{i}.fa"
        r.run(["minimap2", "-x", "map-hifi", "-t", a.threads, current, a.reads], stdout=paf)
        r.run(["racon", "-t", a.threads, a.reads, paf, current], stdout=result)
        current = result
    if not r.dry:
        shutil.copyfile(current, a.out / "polished.fa")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    for name in ("hifiasm", "verkko", "flye"):
        q = sub.add_parser(name)
        q.add_argument("--reads", type=path, required=True, help="combined HiFi FASTQ")
        execution_options(q)
        q.set_defaults(func=assemble)
        if name == "verkko":
            q.add_argument("--memory-gb", type=positive, default=1400)
        elif name == "hifiasm":
            q.add_argument("--hic-r1", type=path)
            q.add_argument("--hic-r2", type=path)
    q = sub.add_parser("ragtag", help="individual correction, scaffolding or patching step")
    q.add_argument("--action", choices=["correct", "scaffold", "patch"], required=True)
    q.add_argument("--reference", type=path, required=True)
    q.add_argument("--query", type=path, required=True)
    q.add_argument("--aligner", choices=["minimap2", "nucmer", "unimap"])
    execution_options(q)
    q.set_defaults(func=ragtag)
    q = sub.add_parser("scaffold", help="weighted RagTag merge followed by two donor patches")
    for name in ("backbone", "framework", "reference", "donor"):
        q.add_argument(f"--{name}", type=path, required=True)
    execution_options(q)
    q.set_defaults(func=scaffold)
    q = sub.add_parser("polish")
    q.add_argument("--genome", type=path, required=True)
    q.add_argument("--reads", type=path, required=True)
    q.add_argument("--rounds", type=positive, default=2)
    execution_options(q)
    q.set_defaults(func=polish)
    a = p.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
