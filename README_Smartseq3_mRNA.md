# Drosophila melanogaster RNA-seq Pipeline (dm6 / SS3)

A paired-end RNA-seq analysis pipeline for *Drosophila melanogaster* (reference genome dm6), with a batch-processing wrapper for running multiple samples sequentially or in parallel.

---

> [!IMPORTANT]
> **Before running, you must update all system-specific paths.** The table below lists every path that needs to be changed. They are also marked with `# <<< EDIT THIS PATH` in the scripts.
>
> | Variable / Location | What it points to |
> |---------------------|-------------------|
> | `PERLCODE` | Directory containing the RPKM Perl script |
> | `SCRIPTDIR` | Directory containing helper scripts |
> | `GENOME_FAI` | dm6 genome `.fa.fai` index file |
> | `GTF_FILE` | dm6 RefSeq GTF annotation file |
> | `STAR_INDICES` | STAR genome index directory |
> | `DM6_BEDFILE` | dm6 BED file for RSeQC |
> | `DM6_CDS_LENGTH` | CDS lengths TSV for TPM calculation |
> | `PICARD_JAR` | Full path to `picard.jar` |
> | `htseq_to_tpm_CDS.py` (in script body) | Python TPM calculation script |
> | `read_distribution.py` (in script body) | RSeQC script path |
> | `PIPELINE_SCRIPT` in `dm6_SS3_batch_process.sh` | Full path to the single-sample pipeline script |

---

## Overview

This repository contains two shell scripts:

| Script | Purpose |
|--------|---------|
| `dm6_SS3_annotated_2026.sh` | Processes a **single sample** through 9 steps: QC trimming → alignment → deduplication → quantification → QC reporting |
| `dm6_SS3_batch_process.sh` | Wraps the single-sample pipeline to process **multiple samples** from a directory, with optional parallelism and resume support |

---

## Pipeline Steps

The single-sample pipeline (`dm6_SS3_annotated_2026.sh`) runs the following steps:

1. **Quality trimming** — `cutadapt` removes low-quality bases and short reads
2. **Alignment** — `STAR` aligns paired-end reads to the dm6 genome
3. **BAM conversion & sorting** — `samtools` converts, sorts, and indexes the alignment
4. **Duplicate removal** — `Picard MarkDuplicates` removes PCR duplicates
5. **Sorted SAM generation** — produces a sorted SAM for downstream tools
6. **BigWig coverage** — `bamCoverage` (deepTools) generates a CPM-normalized `.bw` track
7. **Gene quantification** — `HTSeq` counts reads over whole genes and CDS regions
8. **Normalized expression** — calculates RPKM (via Perl) and TPM (via Python) values
9. **QC & transcript assembly** — `RSeQC` read distribution analysis; `StringTie` transcript assembly

---

## Requirements

### Software

| Tool | Minimum Version |
|------|----------------|
| cutadapt | 3.0 |
| STAR | 2.7 |
| samtools | 1.10 |
| Java | (for Picard) |
| Picard | 2.25.0 |
| HTSeq | 0.11 |
| StringTie | 2.0 |
| deepTools (`bamCoverage`) | any recent |
| RSeQC (`read_distribution.py`) | 4.0+ |
| Python | 3.x |
| Perl | 5.x |

### Reference Files

Before running the pipeline, update the configuration paths at the top of `dm6_SS3_annotated_2026.sh` to point to your local copies of:

- `genome.fa.fai` — dm6 genome FASTA index
- `dm6.ncbiRefSeq.gtf` — RefSeq GTF annotation
- STAR genome index directory
- dm6 BED file for RSeQC
- CDS length TSV (`gene_cds_lengths_by_symbol.tsv`)
- `picard.jar`
- `gene_expressions_from_raw_reads_dm6.pl` — Perl RPKM normalization script
- `htseq_to_tpm_CDS.py` — Python TPM calculation script

---

## Installation

```bash
git clone https://github.com/<your-username>/<your-repo>.git
cd <your-repo>
chmod +x dm6_SS3_annotated_2026.sh dm6_SS3_batch_process.sh
```

Edit the **Configuration Variables** section near the top of `dm6_SS3_annotated_2026.sh` to match your system paths and `PIPELINE_SCRIPT` in `dm6_SS3_batch_process.sh` to point to the same script.

---

## Usage

### Single Sample

```bash
./dm6_SS3_annotated_2026.sh <read1.fastq.gz> <read2.fastq.gz> <sample_prefix>
```

**Example:**
```bash
./dm6_SS3_annotated_2026.sh CS_248_R1_0001.fastq.gz CS_248_R2_0001.fastq.gz CS_248
```

### Batch Processing

```bash
./dm6_SS3_batch_process.sh [OPTIONS] <input_directory> <output_directory>
```

**Options:**

| Flag | Description | Default |
|------|-------------|---------|
| `-p, --parallel N` | Process N samples simultaneously | 1 |
| `-s, --suffix STR` | Input file suffix | `_0001.fastq.gz` |
| `-r, --resume` | Skip samples that already have output | off |
| `-h, --help` | Show help message | — |

**Examples:**
```bash
# Process all samples one at a time
./dm6_SS3_batch_process.sh raw_data/ results/

# Process 4 samples in parallel
./dm6_SS3_batch_process.sh -p 4 raw_data/ results/

# Resume a previously interrupted run
./dm6_SS3_batch_process.sh --resume raw_data/ results/
```

---

## Input File Naming

FASTQ files must follow one of these naming conventions:

```
{SAMPLE}_R1_0001.fastq.gz  /  {SAMPLE}_R2_0001.fastq.gz
{SAMPLE}_1_0001.fastq.gz   /  {SAMPLE}_2_0001.fastq.gz
```

The `--suffix` flag can be used to accommodate other suffixes.

---

## Output Structure

Each sample produces an output directory `{SAMPLE_PREFIX}_analysis/` with the following layout:

```
{SAMPLE_PREFIX}_analysis/
├── bam/
│   ├── {prefix}.sort.bam          # Final deduplicated, sorted BAM
│   ├── {prefix}.sort.bam.bai      # BAM index
│   └── {prefix}_sorted.sam        # Sorted SAM
├── counts/
│   ├── {prefix}_unique_raw.txt          # HTSeq whole-gene raw counts
│   ├── {prefix}_unique_CDS.raw.txt      # HTSeq CDS raw counts
│   ├── {prefix}_expressions.xls         # RPKM (whole gene)
│   ├── {prefix}_expressions_CDS.xls     # RPKM (CDS)
│   ├── {prefix}_mRNA_CDS_TPM.tsv        # TPM values (CDS)
│   ├── {prefix}_stringtie_assembled_transcript.txt.gz
│   └── {prefix}_gene_abundance.txt
├── qc/
│   ├── {prefix}_CPM.bw                          # BigWig coverage track
│   ├── {prefix}_dup_metrics.txt                 # Picard duplicate metrics
│   ├── {prefix}_read_counts.txt                 # Mapping statistics
│   └── {prefix}_read_distribution_RSeQC.txt     # RSeQC read distribution
└── logs/
    ├── {prefix}_pipeline.log        # Full pipeline log
    └── {prefix}_pipeline_summary.txt
```

For batch runs, an additional `batch_processing.log` and `batch_summary.txt` are written to the output directory.

---

## Configuration

Key parameters in `dm6_SS3_annotated_2026.sh`:

| Variable | Default | Description |
|----------|---------|-------------|
| `THREADS` | 8 | Number of CPU threads |
| `QUALITY_CUTOFF` | 20 | Cutadapt quality cutoff (Phred) |
| `MIN_LENGTH` | 10 | Minimum read length after trimming |
| `MISMATCH_RATE` | 0.06 | STAR max mismatch rate |
| `MAX_MULTIMAP` | 1 | STAR max multimapping locations |
| `MAX_MATE_GAP` | 25000 | STAR max gap between mates |
| `JAVA_MEM` | 4g | Java heap memory for Picard |
| `CREATE_SUBDIRS` | true | Organise outputs into subdirectories |

---

## Citation / Acknowledgements

If you use this pipeline in published work, please cite the underlying tools:

- [STAR](https://doi.org/10.1093/bioinformatics/bts635) — Dobin et al., 2013
- [cutadapt](https://doi.org/10.14806/ej.17.1.200) — Martin, 2011
- [samtools](https://doi.org/10.1093/bioinformatics/btp352) — Li et al., 2009
- [Picard](https://broadinstitute.github.io/picard/)
- [HTSeq](https://doi.org/10.1093/bioinformatics/btu638) — Anders et al., 2015
- [StringTie](https://doi.org/10.1038/nbt.3122) — Pertea et al., 2015
- [deepTools](https://doi.org/10.1093/nar/gkw257) — Ramírez et al., 2016
- [RSeQC](https://doi.org/10.1093/bioinformatics/bts356) — Wang et al., 2012
