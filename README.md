# Ribo-seq Preprocessing Pipeline

A comprehensive bioinformatics pipeline for processing ribosome profiling (Ribo-seq) sequencing data from *Drosophila melanogaster* (dm6 genome assembly). This pipeline transforms raw FASTQ files into quantified gene expression data with rigorous quality control and UMI-based deduplication.

## Overview

Ribosome profiling (Ribo-seq) provides genome-wide information about translation by sequencing ribosome-protected mRNA fragments. This pipeline processes raw sequencing reads through adapter trimming, small RNA removal, genome alignment, deduplication, and quantification to produce high-quality translation profiles.

## Features

- **Robust error handling** with automatic validation of inputs and dependencies
- **Organized output structure** with separate directories for different file types
- **Comprehensive logging** with timestamps and software version tracking
- **Parallel processing** with configurable thread allocation
- **UMI-based deduplication** to remove PCR duplicates
- **Multi-stage small RNA filtering** to ensure clean mRNA-specific reads
- **Multiple quantification outputs** including raw counts, normalized expressions, and TPM values
- **Quality control reports** at multiple pipeline stages

## Pipeline Workflow

```
Raw FASTQ
    ↓
1. Quality Control (FastQC)
    ↓
2. Adapter Trimming & Length Filtering (cutadapt)
    ↓
3. UMI Extraction (umi_tools)
    ↓
4. Small RNA Removal (Bowtie)
    ↓
5. Genome Mapping (Bowtie)
    ↓
6. UMI Deduplication (umi_tools)
    ↓
7. Final Small RNA Removal (samtools)
    ↓
8. Read Distribution Analysis (RSeQC) & Coverage Tracks (bamCoverage)
    ↓
9. Gene Quantification (HTSeq)
    ↓
Final Outputs: BAM, bigWig, Count Tables, TPM Values
```

## Requirements

### Software Dependencies

- **cutadapt** (≥3.0) - Adapter trimming
- **FastQC** - Quality control
- **umi_tools** - UMI extraction and deduplication
- **Bowtie** (v1) - Read alignment
- **samtools** (≥1.10) - BAM file manipulation
- **bamCoverage** (deepTools) - Coverage track generation
- **HTSeq** - Read counting
- **RSeQC** - Read distribution analysis
- **Python 3** with required libraries
- **Perl** - Gene expression calculations
- **bedtools** - BED file manipulation

### Reference Files Required

Update the following paths in the script for your system:

```bash
GENOME_FAI_FILE="/path/to/dm6/genome.fa.fai"
gtf_file="/path/to/dm6.ncbiRefSeq.gtf"
small_RNA="/path/to/small_RNAs_plus_tRNAs.bed"
BOWTIE_INDEXES="/path/to/bowtie_dm6_index/"
all_sRNA_BOWTIE_INDEXES="/path/to/small_RNA_bowtie_index"
dm6_bedfile="/path/to/dm6_RefSeq_All.bed"
dm6_CDS_length="/path/to/gene_cds_lengths_by_symbol.tsv"
```

## Installation

1. Clone this repository:
```bash
git clone https://github.com/yourusername/riboseq-pipeline.git
cd riboseq-pipeline
```

2. Install dependencies (example using conda):
```bash
conda create -n riboseq python=3.9
conda activate riboseq
conda install -c bioconda cutadapt fastqc umi_tools bowtie samtools deeptools htseq rseqc bedtools perl
```

3. Download or prepare reference files for dm6 genome (see Reference Files section)

4. Update file paths in the script to match your system

5. Make the script executable:
```bash
chmod +x riboseq_pipeline_improved.sh
```

## Usage

### Basic Usage

```bash
./riboseq_pipeline_improved.sh <input.fastq.gz> <output_prefix>
```

### Example

```bash
./riboseq_pipeline_improved.sh sample1_R1.fastq.gz sample1
```

### Advanced Usage

Customize thread count for parallel processing:
```bash
THREADS=16 ./riboseq_pipeline_improved.sh sample1_R1.fastq.gz sample1
```

### Processing Multiple Samples

Create a simple batch script:
```bash
#!/bin/bash
for fastq in *.fastq.gz; do
    prefix=$(basename $fastq .fastq.gz)
    ./riboseq_pipeline_improved.sh $fastq $prefix
done
```

## Output Files

The pipeline creates an organized output directory structure:

```
<output_prefix>_output/
├── bam/
│   ├── <prefix>.dedup.sRNA.bam          # Final filtered alignment
│   └── <prefix>.dedup.sRNA.bam.bai      # BAM index
├── counts/
│   ├── <prefix>_unique_raw.txt          # HTSeq raw counts (all genes)
│   ├── <prefix>_unique_CDS_raw.txt      # HTSeq raw counts (CDS only)
│   ├── <prefix>_expressions.xls         # Normalized expression (all genes)
│   ├── <prefix>_expressions_CDS.xls     # Normalized expression (CDS)
│   └── <prefix>_CDS_TPM.tsv             # TPM values
├── tracks/
│   └── <prefix>_CPM.bw                  # Normalized coverage track (CPM)
├── fastqc/
│   ├── <input>_fastqc.html              # QC report (raw data)
│   └── <prefix>_adaptor_trimmed_fastqc.html  # QC report (trimmed)
├── logs/
│   └── <prefix>_pipeline.log            # Complete pipeline log
├── <prefix>_summary.txt                 # Pipeline summary with read counts
├── <prefix>_read_distribution.txt       # Genomic feature distribution
├── <prefix>_reads.bed                   # Read coordinates
└── <prefix>_midpoints.bed               # P-site positions
```

### Key Output Descriptions

- **BAM file**: Final deduplicated, small RNA-filtered alignment
- **bigWig file**: Genome browser track for visualization (CPM-normalized)
- **Count files**: Raw read counts per gene for downstream analysis
- **TPM file**: Transcripts per million for cross-sample comparisons
- **Summary file**: Read counts at each processing step + metadata
- **Read distribution**: Percentage of reads in CDS, UTRs, intergenic regions, etc.

## Pipeline Parameters

Default parameters can be modified in the configuration section:

| Parameter | Default | Description |
|-----------|---------|-------------|
| MIN_READ_LENGTH | 30 | Minimum read length after trimming (nt) |
| MAX_READ_LENGTH | 50 | Maximum read length after trimming (nt) |
| UMI_LENGTH | 10 | Length of UMI barcode (bp) |
| MISMATCHES_SRNA | 1 | Allowed mismatches for small RNA mapping |
| MISMATCHES_GENOME | 1 | Allowed mismatches for genome mapping |
| ADAPTER_SEQ | TGGAATTCTCGGGTGCCAAGG | Illumina 3' adapter sequence |
| THREADS | 8 | Number of CPU threads to use |

## Quality Control Checkpoints

The pipeline includes multiple QC steps:

1. **FastQC reports** on raw and trimmed reads
2. **Read count tracking** at each filtering step (logged to summary file)
3. **Mapping statistics** (uniquely mapped, multi-mappers, unmapped)
4. **Deduplication metrics** (reads before/after UMI deduplication)
5. **Read distribution analysis** (genomic feature assignment)

Review the `*_summary.txt` file to ensure:
- Adapter trimming is efficient (>80% reads retained)
- Small RNA contamination is low (<10% of reads)
- Genome mapping rate is acceptable (>70%)
- Deduplication rate is reasonable (typically 20-50% depending on library)

## Troubleshooting

### Common Issues

**Error: Input file not found**
- Verify the input FASTQ file path is correct
- Ensure the file has read permissions

**Error: Required tool not found in PATH**
- Install missing dependencies
- Activate the conda environment if using one

**Error: Required file not found**
- Update file paths in the configuration section
- Verify all reference files are downloaded and uncompressed

**Low mapping rate**
- Check adapter sequence is correct for your library prep
- Verify genome index matches your organism
- Review FastQC reports for quality issues

**High duplication rate**
- May be normal for highly expressed genes
- Check library complexity and sequencing depth
- Verify UMI extraction pattern matches your protocol

## Citation

If you use this pipeline, please cite:

```
Abruzzi et al. (2024). Ribo-seq Preprocessing Pipeline for Drosophila melanogaster.
```

And cite the tools used:
- **cutadapt**: Martin, M. (2011). DOI:10.14806/ej.17.1.200
- **Bowtie**: Langmead, B. et al. (2009). DOI:10.1186/gb-2009-10-3-r25
- **samtools**: Li, H. et al. (2009). DOI:10.1093/bioinformatics/btp352
- **HTSeq**: Anders, S. et al. (2015). DOI:10.1093/bioinformatics/btu638
- **UMI-tools**: Smith, T. et al. (2017). DOI:10.1101/gr.209601.116

## License

This pipeline is provided as-is for academic and research use.

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Commit your changes
4. Submit a pull request

## Support

For questions or issues:
- Open an issue on GitHub
- Contact: [your contact information]

## Version History

- **v2.0** (2024) - Improved version with error handling, parallel processing, and organized outputs
- **v1.0** (2024) - Initial release (Abruzzi et al.)

## Acknowledgments

Developed by the Abruzzi Laboratory for circadian rhythm and ribosome profiling studies in *Drosophila melanogaster*.
