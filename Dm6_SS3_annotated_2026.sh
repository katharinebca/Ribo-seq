!/bin/bash

################################################################################
# RNA-seq Analysis Pipeline - Improved Version
################################################################################
# Description: Processes paired-end RNA-seq data through quality control,
#              alignment, deduplication, and quantification
#
# Author: Improved pipeline with error handling and logging
# Date: 2026-02-04
# Reference genome: dm6 (Drosophila melanogaster)
#
# Requirements:
#   - cutadapt >= 3.0
#   - STAR >= 2.7
#   - samtools >= 1.10
#   - Java (for Picard)
#   - HTSeq >= 0.11
#   - stringtie >= 2.0
#   - deeptools (bamCoverage)
#   - RSeQC
#
# Usage:
#   ./rnaseq_pipeline_improved.sh <read1.fastq> <read2.fastq> <sample_prefix>
#
# Example:
#   ./rnaseq_pipeline_improved.sh sample_R1.fastq.gz sample_R2.fastq.gz MySample
#
################################################################################

# Exit on error, undefined variables, and pipe failures
set -euo pipefail

################################################################################
# Configuration Variables
################################################################################

# Paths - Update these for your system
PERLCODE="/home/analysis/perl_code"
SCRIPTDIR="/home/analysis/scripts"
GENOME_FAI="/home/analysis/genome/dm6/genome.fa.fai"
GTF_FILE="/home/analysis/genome/dm6/refseq_2021/dm6.ncbiRefSeq.gtf"
STAR_INDICES="/home/analysis/genome/dm6/star_indices_2.7.8"
DM6_BEDFILE="/home/katea/dm6_RefSeq_All.Dec_14_2021_RSeQC.bed"
DM6_CDS_LENGTH="/home/analysis/genome/dm6/refseq_2021/gene_cds_lengths_by_symbol.tsv"
PICARD_JAR="/opt/PICARD/2.25.0/picard.jar"

# Processing parameters
THREADS=8
QUALITY_CUTOFF=20
MIN_LENGTH=10
MISMATCH_RATE=0.06
MIN_MATCH=15
MAX_MULTIMAP=1
MAX_MATE_GAP=25000
JAVA_MEM="4g"

# Output organization
CREATE_SUBDIRS=true  # Set to false to keep all files in current directory

################################################################################
# Color codes for output
################################################################################
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

################################################################################
# Functions
################################################################################

# Print colored messages
log_info() {
    echo -e "${BLUE}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} ℹ $1"
}

log_success() {
    echo -e "${GREEN}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} ✓ $1"
}

log_warning() {
    echo -e "${YELLOW}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} ⚠ $1"
}

log_error() {#
    echo -e "${RED}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} ✗ ERROR: $1"
}

# Print section headers
print_header() {
    echo ""
    echo "=========================================================================="
    echo "  $1"
    echo "=========================================================================="
}

# Run a command with error checking
run_step() {
    local step_name="$1"
    shift
    log_info "Running: $step_name"
    
    if "$@"; then
        log_success "$step_name completed"
        return 0
    else
        log_error "$step_name failed with exit code $?"
        return 1
    fi
}

# Check if a file exists and is not empty
validate_file() {
    local file="$1"
    local description="$2"
    
    if [ ! -f "$file" ]; then
        log_error "$description not found: $file"
        exit 1
    fi
    
    if [ ! -s "$file" ]; then
        log_error "$description is empty: $file"
        exit 1
    fi
}

# Check disk space (requires at least 50GB free)
check_disk_space() {
    local required_gb=50
    local available_kb=$(df -k . | awk 'NR==2 {print $4}')
    local available_gb=$((available_kb / 1024 / 1024))
    
    if [ "$available_gb" -lt "$required_gb" ]; then
        log_warning "Low disk space: ${available_gb}GB available (recommended: ${required_gb}GB)"
    else
        log_info "Disk space available: ${available_gb}GB"
    fi
}

# Cleanup function (called on exit)
cleanup() {
    local exit_code=$?
    
    if [ $exit_code -ne 0 ]; then
        log_error "Pipeline failed with exit code $exit_code"
        log_info "Check log file for details: ${LOG_FILE}"
    fi
    
    # Remove temporary directory if it exists
    if [ -n "${TMPDIR:-}" ] && [ -d "$TMPDIR" ]; then
        log_info "Cleaning up temporary directory: $TMPDIR"
        rm -rf "$TMPDIR"
    fi
    
    log_info "Pipeline ended: $(date)"
    log_info "Total runtime: $SECONDS seconds"
}

trap cleanup EXIT

################################################################################
# Input Validation
################################################################################

print_header "RNA-seq Pipeline - Starting"

# Check command line arguments
if [ $# -ne 3 ]; then
    echo "Usage: $0 <read1.fastq> <read2.fastq> <sample_prefix>"
    echo ""
    echo "Arguments:"
    echo "  read1.fastq     : Forward reads (R1) in FASTQ format"
    echo "  read2.fastq     : Reverse reads (R2) in FASTQ format"
    echo "  sample_prefix   : Prefix for output files"
    echo ""
    echo "Example:"
    echo "  $0 sample_R1.fastq.gz sample_R2.fastq.gz MySample"
    exit 1
fi

INPUT_FILE1="$1"
INPUT_FILE2="$2"
SAMPLE_PREFIX="$3"

# Convert to absolute paths to avoid issues when changing directories
INPUT_FILE1=$(realpath "$INPUT_FILE1")
INPUT_FILE2=$(realpath "$INPUT_FILE2")

log_info "Input file 1: $INPUT_FILE1"
log_info "Input file 2: $INPUT_FILE2"
log_info "Sample prefix: $SAMPLE_PREFIX"

# Validate input files
validate_file "$INPUT_FILE1" "Input file 1 (R1)"
validate_file "$INPUT_FILE2" "Input file 2 (R2)"

# Check disk space
check_disk_space

################################################################################
# Create Output Directory Structure
################################################################################

if [ "$CREATE_SUBDIRS" = true ]; then
    OUTDIR="${SAMPLE_PREFIX}_analysis"
    mkdir -p "$OUTDIR"/{bam,counts,qc,logs}
    log_info "Created output directory: $OUTDIR"
    cd "$OUTDIR"
    
    # Set subdirectory paths
    BAM_DIR="bam"
    COUNT_DIR="counts"
    QC_DIR="qc"
    LOG_DIR="logs"
else
    OUTDIR="."
    BAM_DIR="."
    COUNT_DIR="."
    QC_DIR="."
    LOG_DIR="."
fi

# Create temporary directory
TMPDIR=$(mktemp -d -t "${SAMPLE_PREFIX}_XXXXXX")
log_info "Temporary directory: $TMPDIR"

################################################################################
# Setup Logging
################################################################################

LOG_FILE="${LOG_DIR}/${SAMPLE_PREFIX}_pipeline.log"
SUMMARY_FILE="${LOG_DIR}/${SAMPLE_PREFIX}_pipeline_summary.txt"

# Redirect all output to both console and log file
exec 1> >(tee -a "$LOG_FILE")
exec 2>&1

log_info "Log file: $LOG_FILE"

################################################################################
# Check Required Tools
################################################################################

print_header "Checking Required Software"

REQUIRED_TOOLS=(cutadapt STAR samtools java python stringtie bamCoverage perl)
MISSING_TOOLS=()

for tool in "${REQUIRED_TOOLS[@]}"; do
    if command -v "$tool" &> /dev/null; then
        version=$($tool --version 2>&1 | head -n1 || echo "version unknown")
        log_success "$tool found: $version"
    else
        log_error "$tool not found in PATH"
        MISSING_TOOLS+=("$tool")
    fi
done

if [ ${#MISSING_TOOLS[@]} -gt 0 ]; then
    log_error "Missing required tools: ${MISSING_TOOLS[*]}"
    exit 1
fi

# Check required files
print_header "Validating Reference Files"
validate_file "$GENOME_FAI" "Genome FAI file"
validate_file "$GTF_FILE" "GTF annotation file"
validate_file "$DM6_BEDFILE" "BED file for RSeQC"
validate_file "$DM6_CDS_LENGTH" "CDS length file"
validate_file "$PICARD_JAR" "Picard JAR file"

if [ ! -d "$STAR_INDICES" ]; then
    log_error "STAR indices directory not found: $STAR_INDICES"
    exit 1
fi
log_success "All reference files validated"

################################################################################
# Step 1: Quality Trimming with Cutadapt
################################################################################

print_header "Step 1/9: Quality Trimming"

TRIMMED_R1="${SAMPLE_PREFIX}.1.fastq"
TRIMMED_R2="${SAMPLE_PREFIX}.2.fastq"

run_step "Cutadapt trimming" \
    cutadapt \
    -q "$QUALITY_CUTOFF" \
    --trim-n \
    --minimum-length="$MIN_LENGTH" \
    -o "$TRIMMED_R1" \
    -p "$TRIMMED_R2" \
    "$INPUT_FILE1" \
    "$INPUT_FILE2"

# Count total read pairs
TOTAL_READ_PAIRS=$(( $(wc -l < "$TRIMMED_R1") / 4 ))
log_info "Total read pairs after trimming: $TOTAL_READ_PAIRS"

################################################################################
# Step 2: Alignment with STAR
################################################################################

print_header "Step 2/9: STAR Alignment"

run_step "STAR alignment" \
    STAR \
    --runThreadN "$THREADS" \
    --outFilterMismatchNoverLmax "$MISMATCH_RATE" \
    --outFileNamePrefix "${SAMPLE_PREFIX}_" \
    --outFilterMatchNmin "$MIN_MATCH" \
    --outFilterMultimapNmax "$MAX_MULTIMAP" \
    --genomeDir "$STAR_INDICES" \
    --outSJfilterReads Unique \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonical \
    --alignMatesGapMax "$MAX_MATE_GAP" \
    --readFilesIn "$TRIMMED_R1" "$TRIMMED_R2"

# Rename SAM file
SAM_FILE="${BAM_DIR}/${SAMPLE_PREFIX}.sam"
mv "${SAMPLE_PREFIX}_Aligned.out.sam" "$SAM_FILE"
log_success "SAM file created: $SAM_FILE"

# Delete trimmed files to save space
rm "$TRIMMED_R1" "$TRIMMED_R2"
log_info "Removed trimmed FASTQ files"

################################################################################
# Step 3: SAM to BAM Conversion and Sorting
################################################################################

print_header "Step 3/9: BAM Conversion and Sorting"

BAM_FILE="${BAM_DIR}/${SAMPLE_PREFIX}.bam"
SORTED_BAM="${BAM_DIR}/${SAMPLE_PREFIX}.sort.bam"

# Convert to BAM
run_step "SAM to BAM conversion" \
    samtools view -@ "$THREADS" -bhS "$SAM_FILE" -o "$BAM_FILE"

# Sort BAM
run_step "BAM sorting" \
    samtools sort -@ "$THREADS" "$BAM_FILE" -o "$SORTED_BAM"

run_step "BAM indexing" \
    samtools index "$SORTED_BAM"

rm "$BAM_FILE"

# Count mapped reads (pre-deduplication)
MAPPED_READS=$(samtools view -c "$SORTED_BAM")
MAPPED_READ_PAIRS=$((MAPPED_READS / 2))
log_info "Mapped read pairs (pre-dedup): $MAPPED_READ_PAIRS"

################################################################################
# Step 4: Remove Duplicates with Picard
################################################################################

print_header "Step 4/9: Duplicate Removal"

DEDUP_BAM="${BAM_DIR}/${SAMPLE_PREFIX}_nodup.bam"
DUP_METRICS="${QC_DIR}/${SAMPLE_PREFIX}_dup_metrics.txt"

run_step "Picard MarkDuplicates" \
    java -Xmx"$JAVA_MEM" -jar "$PICARD_JAR" MarkDuplicates \
    --INPUT "$SORTED_BAM" \
    --OUTPUT "$DEDUP_BAM" \
    --METRICS_FILE "$DUP_METRICS" \
    --VALIDATION_STRINGENCY LENIENT \
    --REMOVE_DUPLICATES true \
    --TMP_DIR "$TMPDIR" \
    --ASSUME_SORTED true

rm "$SORTED_BAM" "${SORTED_BAM}.bai"

# Re-sort and index deduplicated BAM
FINAL_BAM="${BAM_DIR}/${SAMPLE_PREFIX}.sort.bam"
run_step "Re-sorting deduplicated BAM" \
    samtools sort -@ "$THREADS" "$DEDUP_BAM" -o "$FINAL_BAM"

run_step "Indexing final BAM" \
    samtools index "$FINAL_BAM"

rm "$DEDUP_BAM"

# Count deduplicated reads
DEDUP_READS=$(samtools view -c "$FINAL_BAM")
DEDUP_READ_PAIRS=$((DEDUP_READS / 2))
log_info "Deduplicated read pairs: $DEDUP_READ_PAIRS"

################################################################################
# Step 5: Generate Sorted SAM for Downstream Analysis
################################################################################

print_header "Step 5/9: Generate Sorted SAM"

SORTED_SAM="${BAM_DIR}/${SAMPLE_PREFIX}_sorted.sam"
run_step "Generate sorted SAM" \
    samtools view -h "$FINAL_BAM" -o "$SORTED_SAM"

################################################################################
# Step 6: Generate BigWig Coverage Track
################################################################################

print_header "Step 6/9: Generate BigWig Coverage"

BIGWIG_FILE="${QC_DIR}/${SAMPLE_PREFIX}_CPM.bw"
run_step "bamCoverage (CPM normalization)" \
    bamCoverage \
    --bam "$FINAL_BAM" \
    --normalizeUsing CPM \
    --ignoreForNormalization chrM \
    --numberOfProcessors "$THREADS" \
    -o "$BIGWIG_FILE"

################################################################################
# Step 7: Gene Quantification with HTSeq
################################################################################

print_header "Step 7/9: Gene Quantification"

RAW_COUNTS="${COUNT_DIR}/${SAMPLE_PREFIX}_unique_raw.txt"
CDS_COUNTS="${COUNT_DIR}/${SAMPLE_PREFIX}_unique_CDS.raw.txt"

run_step "HTSeq count (whole gene)" \
    python -m HTSeq.scripts.count \
    -q -r pos -f bam \
    --stranded=no \
    "$FINAL_BAM" \
    "$GTF_FILE" > "$RAW_COUNTS"

run_step "HTSeq count (CDS only)" \
    python -m HTSeq.scripts.count \
    -q -r pos -t CDS -f bam \
    --stranded=no \
    "$FINAL_BAM" \
    "$GTF_FILE" > "$CDS_COUNTS"

################################################################################
# Step 8: Calculate Normalized Expression Values
################################################################################

print_header "Step 8/9: Normalized Expression Calculation"

# Calculate mapping percentages
PCT_MAPPED=$(awk -v m="$MAPPED_READ_PAIRS" -v t="$TOTAL_READ_PAIRS" \
    'BEGIN { if (t>0) printf "%.2f", (m/t)*100; else printf "0.00" }')

PCT_AFTER_DEDUP=$(awk -v d="$DEDUP_READ_PAIRS" -v m="$MAPPED_READ_PAIRS" \
    'BEGIN { if (m>0) printf "%.2f", (d/m)*100; else printf "0.00" }')

# Write read count summary
READ_COUNTS_FILE="${QC_DIR}/${SAMPLE_PREFIX}_read_counts.txt"
{
    printf "sample\t%s\n" "$SAMPLE_PREFIX"
    printf "total_read_pairs\t%s\n" "$TOTAL_READ_PAIRS"
    printf "mapped_read_pairs\t%s\n" "$MAPPED_READ_PAIRS"
    printf "pct_mapped\t%s\n" "$PCT_MAPPED"
    printf "deduplicated_read_pairs\t%s\n" "$DEDUP_READ_PAIRS"
    printf "pct_remaining_after_dedup\t%s\n" "$PCT_AFTER_DEDUP"
} > "$READ_COUNTS_FILE"

log_success "Read count summary written to $READ_COUNTS_FILE"

# Generate normalized read counts (RPKM/FPKM)
EXPR_FILE="${COUNT_DIR}/${SAMPLE_PREFIX}_expressions.xls"
CDS_EXPR_FILE="${COUNT_DIR}/${SAMPLE_PREFIX}_expressions_CDS.xls"

if [ -f "$PERLCODE/gene_expressions_from_raw_reads_dm6.pl" ]; then
    run_step "Calculate normalized expressions (whole gene)" \
        perl "$PERLCODE/gene_expressions_from_raw_reads_dm6.pl" \
        -r "$DEDUP_READ_PAIRS" "$RAW_COUNTS" > "$EXPR_FILE"
    
    run_step "Calculate normalized expressions (CDS)" \
        perl "$PERLCODE/gene_expressions_from_raw_reads_dm6.pl" \
        -r "$DEDUP_READ_PAIRS" "$CDS_COUNTS" > "$CDS_EXPR_FILE"
else
    log_warning "Perl script not found, skipping RPKM calculation"
fi

# Generate TPM normalized counts
TPM_FILE="${COUNT_DIR}/${SAMPLE_PREFIX}_mRNA_CDS_TPM.tsv"
if [ -f "/home/katea/scripts/htseq_to_tpm_CDS.py" ]; then
    run_step "Calculate TPM (CDS)" \
        python /home/katea/scripts/htseq_to_tpm_CDS.py \
        "$CDS_COUNTS" \
        "$DM6_CDS_LENGTH" \
        "$TPM_FILE"
else
    log_warning "TPM calculation script not found, skipping"
fi

################################################################################
# Step 9: Read Distribution Analysis and Transcript Assembly
################################################################################

print_header "Step 9/9: QC and Transcript Assembly"

# Read distribution with RSeQC
READ_DIST_FILE="${QC_DIR}/${SAMPLE_PREFIX}_read_distribution_RSeQC.txt"
if command -v read_distribution.py &> /dev/null; then
    run_step "RSeQC read distribution" \
        python /home/katea/RSeQC-4.0.0/scripts/read_distribution.py \
        -i "$FINAL_BAM" \
        -r "$DM6_BEDFILE" > "$READ_DIST_FILE"
else
    log_warning "RSeQC not found, skipping read distribution analysis"
fi

# StringTie transcript assembly
ASSEMBLED_TRANSCRIPT="${COUNT_DIR}/${SAMPLE_PREFIX}_stringtie_assembled_transcript.txt"
GENE_ABUNDANCE="${COUNT_DIR}/${SAMPLE_PREFIX}_gene_abundance.txt"

run_step "StringTie transcript assembly" \
    stringtie \
    -A "$GENE_ABUNDANCE" \
    -o "$ASSEMBLED_TRANSCRIPT" \
    -p "$THREADS" \
    -e \
    -G "$GTF_FILE" \
    "$FINAL_BAM"

# Compress large output files
run_step "Compress assembled transcripts" \
    gzip -f "$ASSEMBLED_TRANSCRIPT"

################################################################################
# Generate Final Summary Report
################################################################################

print_header "Generating Summary Report"

{
    echo "=================================="
    echo "RNA-seq Pipeline Summary Report"
    echo "=================================="
    echo ""
    echo "Sample: $SAMPLE_PREFIX"
    echo "Date: $(date)"
    echo "Pipeline version: 2.0 (improved)"
    echo ""
    echo "Input Files:"
    echo "  R1: $INPUT_FILE1"
    echo "  R2: $INPUT_FILE2"
    echo ""
    echo "Processing Statistics:"
    echo "  Total read pairs: $TOTAL_READ_PAIRS"
    echo "  Mapped read pairs: $MAPPED_READ_PAIRS ($PCT_MAPPED%)"
    echo "  Deduplicated read pairs: $DEDUP_READ_PAIRS ($PCT_AFTER_DEDUP% of mapped)"
    echo ""
    echo "Output Files:"
    if [ "$CREATE_SUBDIRS" = true ]; then
        echo "  Structure: $OUTDIR/"
        echo "    ├── bam/          (BAM/SAM files)"
        echo "    ├── counts/       (gene counts & expression)"
        echo "    ├── qc/           (QC metrics & reports)"
        echo "    └── logs/         (pipeline logs)"
        echo ""
    fi
    echo "  Final BAM: $FINAL_BAM"
    echo "  BigWig coverage: $BIGWIG_FILE"
    echo "  Raw counts: $RAW_COUNTS"
    echo "  CDS counts: $CDS_COUNTS"
    echo "  Gene abundance: $GENE_ABUNDANCE"
    if [ -f "$TPM_FILE" ]; then
        echo "  TPM values: $TPM_FILE"
    fi
    echo ""
    echo "QC Files:"
    echo "  Read counts: $READ_COUNTS_FILE"
    echo "  Duplicate metrics: $DUP_METRICS"
    if [ -f "$READ_DIST_FILE" ]; then
        echo "  Read distribution: $READ_DIST_FILE"
    fi
    echo ""
    echo "Parameters Used:"
    echo "  Threads: $THREADS"
    echo "  Quality cutoff: $QUALITY_CUTOFF"
    echo "  Min read length: $MIN_LENGTH"
    echo "  Mismatch rate: $MISMATCH_RATE"
    echo ""
} > "$SUMMARY_FILE"

log_success "Summary report written to $SUMMARY_FILE"
cat "$SUMMARY_FILE"

################################################################################
# Final Steps
################################################################################

print_header "Pipeline Completed Successfully"

log_success "All steps completed without errors"
log_info "Total runtime: $SECONDS seconds ($((SECONDS/60)) minutes)"
log_info "Output directory: $OUTDIR"
log_info "Log file: $LOG_FILE"
log_info "Summary: $SUMMARY_FILE"

echo ""
echo "Next steps:"
echo "  1. Review summary report: $SUMMARY_FILE"
echo "  2. Check quality metrics: $READ_COUNTS_FILE"
echo "  3. Examine gene counts: $RAW_COUNTS"
echo "  4. View coverage in IGV: $BIGWIG_FILE"
echo ""

exit 0

