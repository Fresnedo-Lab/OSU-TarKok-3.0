# OSU TarKok 3.0

Scripts for assembly, scaffolding and annotation of the *Taraxacum kok-saghyz* OH-1004 genome.

Use Python 3.10+ and install `requirements.txt`. Each analysis requires its tools and databases to be installed separately. Programs are called from `PATH`; inputs, databases and output locations are command-line arguments.

| Script | Contents |
|---|---|
| `assembly.py` | hifiasm, Verkko, Flye, RagTag and Racon |
| `structural.py` | Trim Galore, HISAT2, BRAKER3 and StringTie UTRs |
| `functional.py` | InterProScan, eggNOG, iTAK, dbCAN, Kofam, TargetP, plantiSMASH and targeted searches |
| `ncrna.py` | Barrnap, tRNAscan-SE, Infernal windows and merged ncRNAs |
| `repeats.py` | RepeatModeler, RepeatMasker, EDTA, LAI, TEsorter and GFF3 export |
| `chloroplast.py` | Reference-guided consensus and PGA annotation |
| `plastome_analysis.py` | Orientation, inverted repeats and gene recovery |
| `plastome_features.py` | OH-1004 feature curation and sequence export |
| `release.py` | Sequence filtering, renaming and representative proteins |
| `orthology.py` | Representative proteomes, OrthoFinder and CAFE5 |
| `qc.py` | Assembly statistics, BUSCO, QUAST and RNA-seq mapping rates |

Scripts are in `scripts/`; `_common.py` contains shared command and FASTA functions. Run `python scripts/SCRIPT.py --help`, then `python scripts/SCRIPT.py COMMAND --help` for arguments. Analysis commands accept `--threads` and `--dry-run`. Use a new output directory for each run.

```bash
python scripts/assembly.py verkko --reads hifi.fq.gz --out verkko --threads 32
python scripts/assembly.py ragtag --action scaffold --reference framework.fa --query assembly.fa --aligner minimap2 --out scaffold
python scripts/assembly.py polish --genome scaffold/results/ragtag.scaffold.fasta --reads hifi.fq.gz --out polished
python scripts/assembly.py scaffold --backbone backbone.fa --framework framework.fa --reference tk1151.fa --donor hifiasm.fa --out combined
python scripts/assembly.py polish --genome combined/scaffolded.fa --reads hifi.fq.gz --out final
```

Backbone preparation uses RagTag scaffolding, two Racon rounds, RagTag correction and scaffolding, then two more Racon rounds. The final merge weights the chromosome framework 2:1 against TK1151 and patches with hifiasm and the framework. The `ragtag` and `polish` subcommands cover these individual steps.

```bash
python scripts/repeats.py mask --genome genome.fa --libraries denovo_repeats.fa plant_repeats.fa --out masked
python scripts/structural.py map --genome genome.fa --samples trimmed_reads.tsv --out rnaseq
python scripts/structural.py braker --genome masked/genome.fa.masked --bam rnaseq/rnaseq.bam --proteins evidence.fa --species taraxacum --out braker
python scripts/functional.py prepare --proteins braker/braker.aa --out proteome
python scripts/functional.py eggnog --proteins proteome/proteins.fa --db eggnog_db --out eggnog
python scripts/functional_tables.py eggnog --inputs eggnog/proteins.emapper.annotations --out eggnog.tsv
python scripts/functional_tables.py merge --proteins proteome/proteins.fa --gff3 braker/braker.gff3 --inputs eggnog.tsv itak.tsv plantismash.tsv --out annotation.tsv
```

RNA sample tables have three tab-separated columns without a header: sample, R1, R2. Read paths are relative to the table. BRAKER needs a writable `AUGUSTUS_CONFIG_PATH`. For UTRs, supply BRAKER's `stringtie2utr.py` and its Python dependencies.

Functional tables use the TSV filename as the column prefix. Name the cluster table `plantismash.tsv` to retain multiple cluster memberships. Protein IDs must match the proteome; gene-level rows expand through the GFF3. TargetP uses plant mode. Rubber-reference IDs use prefixes such as `CPT__` or `SRPP__`.

The plastome consensus uses HiFi reads aligned to OR990541.1. The default sampling fraction, 0.043, belongs to this dataset and should be adjusted for other read depths. PGA uses a directory of reference GenBank files, including OR990541.1 and NC_007578.1. `curate-oh1004` applies the OH-1004 corrections to a coordinate-matched 151,316 bp sequence; its feature table can be processed with table2asn before FASTA export.

`release.py representatives` selects the longest translation-valid, boundary-complete isoform per gene, falling back to a valid partial isoform when necessary. `orthology.py representatives` supports NCBI and GWH. OrthoFinder uses DIAMOND; CAFE excludes families with over 100 copies in any taxon. The ultrametric root age of 100 is a scaling parameter.

Tool/database versions used include InterProScan 5.78-109.0, eggNOG-mapper 2.1.12 with database 5.0.2, iTAK 2.0.7, dbCAN 5, TargetP 2, Rfam 15.1, EDTA 2.2.2 and BUSCO 6 with `eudicotyledons_odb12`. Keep software and database releases compatible. Reference sequences and licensed databases are supplied separately.

Validation covers command construction and comparisons with existing sequence and annotation outputs. Full assembly and annotation programs were not rerun during script preparation.
