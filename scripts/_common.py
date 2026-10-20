"""Command execution and sequence I/O shared by the workflow scripts."""
import argparse
import gzip
import os
from pathlib import Path
import shlex
import subprocess
import sys


def path(value):
    return Path(value).expanduser().resolve()


def positive(value):
    value = int(value)
    if value < 1:
        raise argparse.ArgumentTypeError("must be positive")
    return value


def execution_options(parser):
    parser.add_argument("--out", type=path, required=True, help="new output directory")
    parser.add_argument("--threads", type=positive, default=8)
    parser.add_argument("--dry-run", action="store_true", help="print commands without running them")


class Runner:
    def __init__(self, args):
        self.dry = args.dry_run
        self.out = args.out
        if not self.dry:
            self.out.mkdir(parents=True, exist_ok=False)

    def run(self, command, *, stdout=None, cwd=None, env=None):
        self.pipe([command], stdout=stdout, cwd=cwd, env=env)

    def pipe(self, commands, *, stdout=None, cwd=None, env=None):
        commands = [[str(x) for x in cmd] for cmd in commands]
        line = " | ".join(shlex.join(cmd) for cmd in commands)
        if stdout:
            line += " > " + shlex.quote(str(stdout))
        print(line, file=sys.stderr)
        if self.dry:
            return
        processes = []
        destination = open(stdout, "wb") if stdout else None
        previous = None
        try:
            for i, command in enumerate(commands):
                process = subprocess.Popen(
                    command, cwd=cwd or self.out,
                    env=None if env is None else {**os.environ, **env},
                    stdin=previous,
                    stdout=subprocess.PIPE if i < len(commands) - 1 else destination,
                )
                if previous is not None:
                    previous.close()
                previous = process.stdout
                processes.append(process)
            codes = [process.wait() for process in processes]
            for command, code in zip(commands, codes):
                if code:
                    raise subprocess.CalledProcessError(code, command)
        finally:
            if previous is not None:
                previous.close()
            for process in processes:
                if process.poll() is None:
                    process.terminate()
                    process.wait()
            if destination:
                destination.close()

    def write(self, file, text):
        if not self.dry:
            Path(file).write_text(text)


def text_open(file):
    return gzip.open(file, "rt") if str(file).endswith(".gz") else open(file)


def fasta(file):
    """Yield the full header and sequence, one record at a time."""
    header, sequence = None, []
    with text_open(file) as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header, sequence = line[1:].strip(), []
            elif line.strip():
                if header is None:
                    raise ValueError("Sequence before FASTA header")
                sequence.append(line.strip())
    if header is not None:
        yield header, "".join(sequence)


def write_fasta(handle, header, sequence):
    handle.write(f">{header}\n")
    for start in range(0, len(sequence), 60):
        handle.write(sequence[start:start + 60] + "\n")


def attributes(text):
    return dict(item.split("=", 1) for item in text.split(";") if "=" in item)
