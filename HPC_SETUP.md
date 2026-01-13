# HPC Setup Guide

Quick guide for running the PacBio HDV pipeline on HPC clusters.

## Quick Start
```bash
# Connect to HPC
ssh user@hpc-address

# Load conda
module load miniconda3

# Clone repo
git clone https://github.com/onkarS23/pacbio-hdv-pipeline-v2-.git
cd pacbio-hdv-pipeline-v2-

# Setup
bash setup_v2.sh
conda activate pacbio-hdv-v2

# Download reference
cd data
wget -O hg38.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
samtools faidx hg38.fa
cd ..

# Prepare your data
cp /path/to/your/*.fastq test_input/

# Run pipeline
nextflow run main_v2.nf --input_dir test_input --outdir results --threads 16
```

## SLURM Job Script

Create `submit_pipeline.sh`:
```bash
#!/bin/bash
#SBATCH --job-name=hdv_pipeline
#SBATCH --output=hdv_%j.log
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00

module load miniconda3
conda activate pacbio-hdv-v2

nextflow run main_v2.nf \
  --input_dir test_input \
  --outdir results \
  --threads 16 \
  -with-report execution_report.html \
  -with-timeline timeline.html
```

Submit: `sbatch submit_pipeline.sh`

Check status: `squeue -u $USER`

## Troubleshooting

- **Java error**: `conda install -y openjdk=17`
- **Memory error**: Reduce `--threads` or request more memory
- **hg38 exists?**: Check with admin: `ln -s /path/to/hg38.fa data/hg38.fa`

See [README.md](README.md) for complete documentation.
