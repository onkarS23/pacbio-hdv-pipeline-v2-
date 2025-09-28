#!/bin/bash

#SBATCH --job-name=hdv_pipeline
#SBATCH --output=hdv_pipeline_%j.out
#SBATCH --error=hdv_pipeline_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=cpu

# Script to run PacBio HDV Analysis Pipeline
# Usage: bash run_hdv_pipeline.sh [input_dir] [output_dir]

set -e  # Exit on any error

# Default parameters
INPUT_DIR="${1:-test_input}"
OUTPUT_DIR="${2:-hdv_results}"
THREADS="${3:-16}"

# Activate conda environment
echo "Activating conda environment..."
source ~/miniconda3/etc/profile.d/conda.sh
conda activate pacbio-hdv-v2

# Verify environment
echo "Checking environment..."
which nextflow
which medaka
which minimap2
which samtools

# Change to pipeline directory
cd /home/osingh/pacbio-hdv-pipeline/pipeline_v2

# Print pipeline info
echo "========================================="
echo "PacBio HDV Pipeline Execution"
echo "========================================="
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Threads: $THREADS"
echo "Started at: $(date)"
echo "========================================="

# Run the pipeline
echo "Starting Nextflow pipeline..."
nextflow run main_v2.nf \
    --input_dir "$INPUT_DIR" \
    --outdir "$OUTPUT_DIR" \
    --threads "$THREADS" \
    -with-report "${OUTPUT_DIR}/execution_report.html" \
    -with-timeline "${OUTPUT_DIR}/timeline.html" \
    -with-dag "${OUTPUT_DIR}/flowchart.html" \
    -resume

# Check completion status
if [ $? -eq 0 ]; then
    echo "========================================="
    echo "Pipeline completed successfully!"
    echo "Completed at: $(date)"
    echo "Results in: $OUTPUT_DIR"
    echo "========================================="
else
    echo "========================================="
    echo "Pipeline failed!"
    echo "Failed at: $(date)"
    echo "Check logs for details"
    echo "========================================="
    exit 1
fi

# Generate summary
echo "Generating summary..."
find "$OUTPUT_DIR" -name "*.txt" -exec echo "=== {} ===" \; -exec head -10 {} \;

echo "Pipeline execution completed. Check $OUTPUT_DIR for results."