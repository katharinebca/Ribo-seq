#!/bin/bash

################################################################################
# Batch RNA-seq Pipeline Processing Script
################################################################################
# Description: Process multiple RNA-seq samples in batch
# Usage: ./batch_process.sh <input_directory> <output_directory>
#
# Expected input structure:
#   input_dir/
#     ├── CS_248_R1_0001.fastq.gz
#     ├── CS_248_R2_0001.fastq.gz
#     ├── Sample2_R1_0001.fastq.gz
#     ├── Sample2_R2_0001.fastq.gz
#     └── ...
#
################################################################################

set -euo pipefail

# Color codes
GREEN='\033[0;32m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m'

# Configuration
PIPELINE_SCRIPT="/home/katea/scripts/dm6_SS3_annotated_2026.sh"
MAX_PARALLEL=1  # Number of samples to process simultaneously

################################################################################
# Functions
################################################################################

usage() {
    cat << EOF
Usage: $0 [OPTIONS] <input_directory> <output_directory>

Process multiple paired-end RNA-seq samples in batch.

Arguments:
  input_directory    Directory containing FASTQ files
  output_directory   Directory for output files

Options:
  -p, --parallel N   Process N samples in parallel (default: 1)
  -s, --suffix STR   Input file suffix (default: _0001.fastq.gz)
  -r, --resume       Skip samples that already have output
  -h, --help         Show this help message

Examples:
  # Process all samples one at a time
  $0 raw_data/ results/

  # Process 4 samples in parallel
  $0 -p 4 raw_data/ results/

  # Resume previous run
  $0 --resume raw_data/ results/

Input file naming convention:
  Files should be named: {SAMPLE}_R1_0001.fastq.gz and {SAMPLE}_R2_0001.fastq.gz
  Or: {SAMPLE}_1_0001.fastq.gz and {SAMPLE}_2_0001.fastq.gz

EOF
    exit 1
}

log_info() {
    echo -e "${BLUE}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} ✓ $1"
}

log_error() {
    echo -e "${RED}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} ✗ ERROR: $1"
}

################################################################################
# Parse Arguments
################################################################################

SUFFIX="_0001.fastq.gz"
RESUME=false

while [[ $# -gt 0 ]]; do
    case $1 in
        -p|--parallel)
            MAX_PARALLEL="$2"
            shift 2
            ;;
        -s|--suffix)
            SUFFIX="$2"
            shift 2
            ;;
        -r|--resume)
            RESUME=true
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            break
            ;;
    esac
done

if [ $# -ne 2 ]; then
    usage
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"

# Validate inputs
if [ ! -d "$INPUT_DIR" ]; then
    log_error "Input directory not found: $INPUT_DIR"
    exit 1
fi

if [ ! -f "$PIPELINE_SCRIPT" ]; then
    log_error "Pipeline script not found: $PIPELINE_SCRIPT"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Convert paths to absolute so redirects still work after cd
INPUT_DIR="$(realpath "$INPUT_DIR")"
OUTPUT_DIR="$(realpath "$OUTPUT_DIR")"

################################################################################
# Find Sample Pairs
################################################################################

log_info "Scanning for paired-end samples in: $INPUT_DIR"
log_info "Using suffix: $SUFFIX"

# Find all R1 files
declare -a SAMPLES
for R1 in "$INPUT_DIR"/*_R1${SUFFIX} "$INPUT_DIR"/*_1${SUFFIX}; do
    [ -f "$R1" ] || continue

    # Determine R2 filename
    if [[ "$R1" =~ _R1${SUFFIX}$ ]]; then
        R2="${R1/_R1${SUFFIX}/_R2${SUFFIX}}"
    else
        R2="${R1/_1${SUFFIX}/_2${SUFFIX}}"
    fi

    # Check if R2 exists
    if [ ! -f "$R2" ]; then
        log_error "Paired file not found for $R1"
        continue
    fi

    # Extract sample name
    basename=$(basename "$R1")
    if [[ "$basename" =~ _R1${SUFFIX}$ ]]; then
        sample_name="${basename/_R1${SUFFIX}/}"
    else
        sample_name="${basename/_1${SUFFIX}/}"
    fi

    SAMPLES+=("$sample_name|$R1|$R2")
done

if [ ${#SAMPLES[@]} -eq 0 ]; then
    log_error "No paired-end samples found in $INPUT_DIR"
    log_info "Expected naming: {SAMPLE}_R1${SUFFIX} and {SAMPLE}_R2${SUFFIX}"
    exit 1
fi

log_success "Found ${#SAMPLES[@]} sample pairs"

################################################################################
# Process Samples
################################################################################

BATCH_LOG="${OUTPUT_DIR}/batch_processing.log"
BATCH_SUMMARY="${OUTPUT_DIR}/batch_summary.txt"

# Initialize summary file
{
    echo "Batch RNA-seq Processing Summary"
    echo "================================="
    echo "Started: $(date)"
    echo "Input directory: $INPUT_DIR"
    echo "Output directory: $OUTPUT_DIR"
    echo "Parallel jobs: $MAX_PARALLEL"
    echo "Total samples: ${#SAMPLES[@]}"
    echo ""
} > "$BATCH_SUMMARY"

log_info "Processing ${#SAMPLES[@]} samples (max $MAX_PARALLEL in parallel)"
log_info "Batch log: $BATCH_LOG"

PROCESSED=0
FAILED=0
SKIPPED=0

# Process samples
for sample_info in "${SAMPLES[@]}"; do
    IFS='|' read -r sample_name R1 R2 <<< "$sample_info"

    # Check if sample already processed (resume mode)
    if [ "$RESUME" = true ]; then
        if [ -f "$OUTPUT_DIR/${sample_name}.sort.bam" ]; then
            log_info "Skipping $sample_name (already processed)"
            # Avoid ((SKIPPED++)) under set -e
            SKIPPED=$((SKIPPED+1))
            continue
        fi
    fi

    log_info "Processing sample: $sample_name"

    # Wait if we've reached max parallel jobs
    while [ "$(jobs -r | wc -l)" -ge "$MAX_PARALLEL" ]; do
        sleep 10
    done

    # Run pipeline in background
    (
        cd "$OUTPUT_DIR"
        if "$PIPELINE_SCRIPT" "$R1" "$R2" "$sample_name" >> "$BATCH_LOG" 2>&1; then
            log_success "Completed: $sample_name"
            echo "SUCCESS: $sample_name" >> "$BATCH_SUMMARY"
        else
            log_error "Failed: $sample_name"
            echo "FAILED: $sample_name" >> "$BATCH_SUMMARY"
        fi
    ) &

    # Avoid ((PROCESSED++)) under set -e
    PROCESSED=$((PROCESSED+1))
done

# Wait for all background jobs to finish
log_info "Waiting for all jobs to complete..."
wait

################################################################################
# Generate Final Summary
################################################################################

# Count successes and failures
SUCCEEDED=$(grep -c "^SUCCESS:" "$BATCH_SUMMARY" || echo 0)
FAILED=$(grep -c "^FAILED:" "$BATCH_SUMMARY" || echo 0)

{
    echo ""
    echo "Batch Processing Complete"
    echo "========================="
    echo "Completed: $(date)"
    echo "Total samples: ${#SAMPLES[@]}"
    echo "Successfully processed: $SUCCEEDED"
    echo "Failed: $FAILED"
    echo "Skipped (already done): $SKIPPED"
    echo ""
    echo "Output files in: $OUTPUT_DIR"
    echo ""

    if [ "$FAILED" -gt 0 ]; then
        echo "Failed samples:"
        grep "^FAILED:" "$BATCH_SUMMARY" | sed 's/^FAILED: /  - /'
        echo ""
    fi
} >> "$BATCH_SUMMARY"

log_info "============================================"
log_success "Batch processing complete!"
log_info "Total samples: ${#SAMPLES[@]}"
log_info "Succeeded: $SUCCEEDED"
if [ "$FAILED" -gt 0 ]; then
    log_error "Failed: $FAILED"
fi
if [ "$SKIPPED" -gt 0 ]; then
    log_info "Skipped: $SKIPPED"
fi
log_info "Summary: $BATCH_SUMMARY"
log_info "Detailed log: $BATCH_LOG"
log_info "============================================"

cat "$BATCH_SUMMARY"

exit 0
