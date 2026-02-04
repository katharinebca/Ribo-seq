#!/bin/bash

################################################################################
# Ribo-seq Preprocessing Pipeline (Abruzzi 2024 - Improved Version)
################################################################################
#
# Description: Processes ribosome profiling sequencing data from raw reads to
#              gene expression quantification for Drosophila melanogaster (dm6)
#
# Usage: ./riboseq_pipeline_improved.sh <input.fastq.gz> <output_prefix>
#
# Requirements:
#   - cutadapt, fastqc, umi_tools
#   - bowtie, samtools, bamCoverage
#   - HTSeq, RSeQC
#   - Python 3 with required libraries
#
# Author: Abruzzi Lab
# Modified: 2024
################################################################################

# Exit on any error, undefined variable, or pipe failure
set -euo pipefail

# Optional: Enable command tracing for debugging (uncomment if needed)
# set -x

################################################################################
# CONFIGURATION
################################################################################

# Input parameters
if [ $# -ne 2 ]; then
    echo "ERROR: Incorrect number of arguments"
    echo "Usage: $0 <input_fastq.gz> <output_prefix>"
    echo ""
    echo "Example: $0 sample1.fastq.gz sample1"
    exit 1
fi

infile1=$1
prefix1=$2

# Processing parameters
MIN_READ_LENGTH=30
MAX_READ_LENGTH=50
UMI_LENGTH=10
MISMATCHES_SRNA=1
MISMATCHES_GENOME=1
ADAPTER_SEQ="TGGAATTCTCGGGTGCCAAGG"

# Resource allocation
THREADS=${THREADS:-8}  # Use environment variable if set, otherwise default to 8

# File paths - UPDATE THESE FOR YOUR SYSTEM
PERLCODE="/home/analysis/perl_code"
SCRIPTDIR="/home/analysis/scripts"
GENOME_FAI_FILE="/home/analysis/genome/dm6/genome.fa.fai"
gtf_file="/home/analysis/genome/dm6/refseq_2021/dm6.ncbiRefSeq.gtf"
small_RNA="/home/katea/small_RNAs_plus_tRNAs.bed"
BOWTIE_INDEXES="/home/katea/bowtie_dm6_2023/"
all_sRNA_BOWTIE_INDEXES="/home/katea/dm6_bowtie_index_RNAcentral_2026/trsnRNA_2026"
dm6_bedfile="/home/katea/dm6_RefSeq_All.Dec_14_2021_RSeQC.bed"
dm6_CDS_length="/home/analysis/genome/dm6/refseq_2021/gene_cds_lengths_by_symbol.tsv"

################################################################################
# VALIDATION
################################################################################

echo "================================================================================"
echo "Ribo-seq Preprocessing Pipeline - Starting"
echo "================================================================================"
echo "Input file: $infile1"
echo "Output prefix: $prefix1"
echo "Threads: $THREADS"
echo "Start time: $(date)"
echo ""

# Check input file exists
if [ ! -f "$infile1" ]; then
    echo "ERROR: Input file '$infile1' not found"
    exit 1
fi

# Check required directories and files exist
echo "Validating required files and directories..."
required_files=(
    "$GENOME_FAI_FILE"
    "$gtf_file"
    "$small_RNA"
    "$dm6_bedfile"
    "$dm6_CDS_length"
)

required_dirs=(
    "$PERLCODE"
    "$SCRIPTDIR"
    "$BOWTIE_INDEXES"
)

for file in "${required_files[@]}"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: Required file not found: $file"
        exit 1
    fi
done

for dir in "${required_dirs[@]}"; do
    if [ ! -d "$dir" ]; then
        echo "ERROR: Required directory not found: $dir"
        exit 1
    fi
done

# Check required software
echo "Validating required software..."
required_tools=(cutadapt fastqc umi_tools bowtie samtools bamCoverage python)
for tool in "${required_tools[@]}"; do
    if ! command -v "$tool" &> /dev/null; then
        echo "ERROR: Required tool '$tool' not found in PATH"
        exit 1
    fi
done

echo "Validation complete!"
echo ""

################################################################################
# SETUP
################################################################################

# Create organized output directories
mkdir -p "${prefix1}_output"/{fastqc,bam,counts,tracks,logs}

# Initialize summary and log files
output_summary_file="${prefix1}_output/${prefix1}_summary.txt"
log_file="${prefix1}_output/logs/${prefix1}_pipeline.log"

# Redirect all output to log file while displaying on screen
exec > >(tee -a "$log_file")
exec 2>&1

# Initialize summary file with header
cat > "$output_summary_file" << EOF
# Ribo-seq Preprocessing Summary
# Sample prefix: $prefix1
# Input file: $infile1
# Pipeline started: $(date)
# 
# Software versions:
# cutadapt: $(cutadapt --version)
# bowtie: $(bowtie --version 2>&1 | head -1)
# samtools: $(samtools --version | head -1)
# umi_tools: $(umi_tools --version)
#
# Parameters:
# Min read length: $MIN_READ_LENGTH
# Max read length: $MAX_READ_LENGTH
# UMI length: $UMI_LENGTH
# Adapter sequence: $ADAPTER_SEQ
# Threads: $THREADS
#
EOF

################################################################################
# HELPER FUNCTIONS
################################################################################

# Function to count reads in FASTQ file
count_reads() {
    local file=$1
    awk 'END {print NR/4}' "$file"
}

# Function to log read counts to summary
log_count() {
    local description=$1
    local count=$2
    echo "$description: $count" >> "$output_summary_file"
}

# Function to print step headers
print_step() {
    local step_num=$1
    local step_desc=$2
    echo ""
    echo ">>> Step $step_num: $step_desc"
    echo ">>> Time: $(date)"
    echo ""
}

################################################################################
# PIPELINE STEPS
################################################################################

# Step 1: Quality Control on Original Data
print_step "1/9" "Running FastQC on raw data..."
fastqc -t $THREADS -o "${prefix1}_output/fastqc" "$infile1"

# Step 2: Adapter Trimming and Read Length Filtering
print_step "2/9" "Trimming adapters and filtering by read length..."

cutadapt -a "$ADAPTER_SEQ" "$infile1" \
  --untrimmed-output "${prefix1}_no_adaptor.fastq" \
  -m $MIN_READ_LENGTH --too-short-output "${prefix1}_too_short.fastq" \
  -M $MAX_READ_LENGTH --too-long-output "${prefix1}_too_long.fastq" \
  -o "${prefix1}_adaptor_trimmed.fastq" \
  -j $THREADS

# QC on trimmed reads
fastqc -t $THREADS -o "${prefix1}_output/fastqc" "${prefix1}_adaptor_trimmed.fastq"

# Count and log reads
log_count "Reads with no adapter detected" "$(count_reads ${prefix1}_no_adaptor.fastq)"
log_count "Reads too short (<${MIN_READ_LENGTH} nt)" "$(count_reads ${prefix1}_too_short.fastq)"
log_count "Reads too long (>${MAX_READ_LENGTH} nt)" "$(count_reads ${prefix1}_too_long.fastq)"
log_count "Reads retained after adapter trimming (${MIN_READ_LENGTH}-${MAX_READ_LENGTH} nt)" "$(count_reads ${prefix1}_adaptor_trimmed.fastq)"

# Clean up intermediate files
rm "${prefix1}_no_adaptor.fastq" "${prefix1}_too_short.fastq" "${prefix1}_too_long.fastq"

# Step 3: UMI Extraction
print_step "3/9" "Extracting UMIs from 5' end..."

umi_tools extract \
  --stdin="${prefix1}_adaptor_trimmed.fastq" \
  --extract-method=string \
  --bc-pattern=$(printf 'N%.0s' $(seq 1 $UMI_LENGTH)) \
  --stdout="${prefix1}_adaptor_trimmed_UMI.fastq"

rm "${prefix1}_adaptor_trimmed.fastq"

# Step 4: Remove Small RNAs
print_step "4/9" "Mapping and removing small RNAs..."

bowtie --sam -p $THREADS -v $MISMATCHES_SRNA \
  --un "${prefix1}_reads_not_mapped_sRNA.fastq" \
  --al "${prefix1}_sRNA_aligned.fastq" \
  "$all_sRNA_BOWTIE_INDEXES" \
  "${prefix1}_adaptor_trimmed_UMI.fastq" \
  "${prefix1}_map_sRNA.sam"

log_count "Reads mapped to small RNAs" "$(count_reads ${prefix1}_sRNA_aligned.fastq)"
log_count "Reads not mapped to small RNAs (proceeding to genome mapping)" "$(count_reads ${prefix1}_reads_not_mapped_sRNA.fastq)"

# Clean up
rm "${prefix1}_adaptor_trimmed_UMI.fastq" "${prefix1}_map_sRNA.sam"

# Step 5: Genome Mapping
print_step "5/9" "Mapping reads to dm6 genome..."

bowtie --sam -p $THREADS -v $MISMATCHES_GENOME --best \
  --un "${prefix1}_not_mapped_genome.fastq" \
  "$BOWTIE_INDEXES" \
  "${prefix1}_reads_not_mapped_sRNA.fastq" \
  "${prefix1}.sam"

rm "${prefix1}_reads_not_mapped_sRNA.fastq"

# Convert to BAM
samtools view -@ $THREADS -bhS "${prefix1}.sam" > "${prefix1}.bam"
rm "${prefix1}.sam"

# Extract only mapped reads and count
samtools view -h -F 4 "${prefix1}.bam" > "${prefix1}_mapped.sam"
log_count "Reads mapped to genome (before UMI deduplication)" "$(samtools view -c ${prefix1}_mapped.sam)"
rm "${prefix1}_mapped.sam"

# Sort and index BAM
samtools sort -@ $THREADS "${prefix1}.bam" -o "${prefix1}.sort.bam"
samtools index "${prefix1}.sort.bam"
rm "${prefix1}.bam"

# Step 6: UMI Deduplication
print_step "6/9" "Deduplicating reads using UMIs..."

umi_tools dedup \
  -I "${prefix1}.sort.bam" \
  --method=unique \
  --read-length \
  -S "${prefix1}.dedup.sort.bam"

log_count "Reads after UMI deduplication" "$(samtools view -c ${prefix1}.dedup.sort.bam)"

rm "${prefix1}.sort.bam" "${prefix1}.sort.bam.bai"

# Step 7: Final Small RNA Removal
print_step "7/9" "Performing second small RNA removal..."

samtools view -@ $THREADS -L "$small_RNA" \
  -U "${prefix1}.dedup.sRNA.bam" \
  -o "${prefix1}_sRNA_contaminants.bam" \
  "${prefix1}.dedup.sort.bam"

samtools index "${prefix1}.dedup.sRNA.bam"

log_count "Final mRNA reads after second sRNA removal" "$(samtools view -c ${prefix1}.dedup.sRNA.bam)"

rm "${prefix1}.dedup.sort.bam"

# Move final BAM to organized directory
mv "${prefix1}.dedup.sRNA.bam"* "${prefix1}_output/bam/"
final_bam="${prefix1}_output/bam/${prefix1}.dedup.sRNA.bam"

# Step 8: Read Distribution Analysis and Coverage Tracks
print_step "8/9" "Analyzing read distribution and generating coverage tracks..."

# Read distribution
python /home/katea/RSeQC-4.0.0/scripts/read_distribution.py \
  -i "$final_bam" \
  -r "$dm6_bedfile" > "${prefix1}_output/${prefix1}_read_distribution.txt"

# Generate normalized bigWig file
bamCoverage --bam "$final_bam" \
  --normalizeUsing CPM \
  --ignoreForNormalization chrM \
  --numberOfProcessors $THREADS \
  -o "${prefix1}_output/tracks/${prefix1}_CPM.bw"

# Create BED file with read coordinates
convert2bed -i bam < "$final_bam" > "${prefix1}_reads.bed"

# Calculate midpoints for P-site analysis
awk 'BEGIN{OFS="\t"}{
  mid = int(($2 + $3) / 2);
  print $1, mid, mid+1, $4, $5, $6
}' "${prefix1}_reads.bed" > "${prefix1}_midpoints.bed"

# Determine location of midpoints and reads
if [ -f "/home/katea/scripts/dm6_location_midpoints_2026.sh" ]; then
    /home/katea/scripts/dm6_location_midpoints_2026.sh "${prefix1}_midpoints.bed"
fi

if [ -f "/home/katea/scripts/dm6_location_2026.sh" ]; then
    /home/katea/scripts/dm6_location_2026.sh "${prefix1}_reads.bed"
fi

# Step 9: Gene Quantification
print_step "9/9" "Quantifying gene expression..."

# Count reads per gene
python -m HTSeq.scripts.count \
  -q -f bam \
  --stranded=no \
  --minaqual=10 \
  "$final_bam" \
  "$gtf_file" > "${prefix1}_output/counts/${prefix1}_unique_raw.txt"

# Count reads in CDS regions
python -m HTSeq.scripts.count \
  -q -t CDS -f bam \
  --stranded=no \
  --minaqual=10 \
  "$final_bam" \
  "$gtf_file" > "${prefix1}_output/counts/${prefix1}_unique_CDS_raw.txt"

# Calculate gene expressions
readnum=$(samtools view -c "$final_bam")

perl "${PERLCODE}/gene_expressions_from_raw_reads_dm6.pl" \
  -r "$readnum" \
  "${prefix1}_output/counts/${prefix1}_unique_raw.txt" > \
  "${prefix1}_output/counts/${prefix1}_expressions.xls"

perl "${PERLCODE}/gene_expressions_from_raw_reads_dm6.pl" \
  -r "$readnum" \
  "${prefix1}_output/counts/${prefix1}_unique_CDS_raw.txt" > \
  "${prefix1}_output/counts/${prefix1}_expressions_CDS.xls"

# Calculate TPM values
python /home/katea/scripts/htseq_to_tpm_CDS.py \
  "${prefix1}_output/counts/${prefix1}_unique_CDS_raw.txt" \
  "$dm6_CDS_length" \
  "${prefix1}_output/counts/${prefix1}_CDS_TPM.tsv"

################################################################################
# CLEANUP AND FINALIZATION
################################################################################

print_step "Cleanup" "Organizing output files..."

# Move intermediate BED files to output directory
mv "${prefix1}_reads.bed" "${prefix1}_midpoints.bed" "${prefix1}_output/" 2>/dev/null || true

# Compress large intermediate files if they exist
if [ -f "${prefix1}_sRNA_aligned.fastq" ]; then
    gzip "${prefix1}_sRNA_aligned.fastq" &
fi

if [ -f "${prefix1}_not_mapped_genome.fastq" ]; then
    gzip "${prefix1}_not_mapped_genome.fastq" &
fi

# Wait for background compression jobs
wait

# Add completion timestamp
echo "" >> "$output_summary_file"
echo "# Pipeline completed: $(date)" >> "$output_summary_file"

################################################################################
# SUMMARY
################################################################################

echo ""
echo "================================================================================"
echo "Pipeline Completed Successfully!"
echo "================================================================================"
echo "Output directory: ${prefix1}_output/"
echo ""
echo "Key output files:"
echo "  - Summary: $output_summary_file"
echo "  - Final BAM: $final_bam"
echo "  - Coverage track: ${prefix1}_output/tracks/${prefix1}_CPM.bw"
echo "  - Gene counts: ${prefix1}_output/counts/${prefix1}_unique_raw.txt"
echo "  - TPM values: ${prefix1}_output/counts/${prefix1}_CDS_TPM.tsv"
echo "  - Read distribution: ${prefix1}_output/${prefix1}_read_distribution.txt"
echo ""
echo "Log file: $log_file"
echo "End time: $(date)"
echo "================================================================================"

exit 0
