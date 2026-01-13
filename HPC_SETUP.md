# HPC Setup Guide

Quick setup guide for running the pipeline on HPC clusters via MobaXterm.

## Prerequisites

- SSH access to HPC cluster
- MobaXterm or terminal

## Setup Steps

### 1. Connect to HPC
```bash
ssh username@your-hpc-address
```

### 2. Load or Install Conda

**Check for existing conda:**
```bash
module avail conda
```

**If available, load it:**
```bash
module load miniconda3
```

**If NOT available, install:**
```bash
cd ~
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

### 3. Clone and Setup
```bash
git clone https://github.com/onkarS23/pacbio-hdv-pipeline-v2-.git
cd pacbio-hdv-pipeline-v2-
bash setup_v2.sh
conda activate pacbio-hdv-v2
```

### 4. Download hg38 Reference
```bash
cd data
wget -O hg38.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
samtools faidx hg38.fa
cd ..
```

**Or use existing reference (ask admin):**
```bash
ln -s /path/to/hg38.fa data/hg38.fa
```

### 5. Add Your Data
```bash
cp /path/to/your/*.fastq test_input/
```

### 6. Run Pipeline

**Create job script `submit.sh`:**
```bash
#!/bin/bash
#SBATCH --job-name=hdv
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00

module load miniconda3  # or: source ~/miniconda3/etc/profile.d/conda.sh
conda activate pacbio-hdv-v2

nextflow run main_v2.nf \
  --input_dir test_input \
  --outdir results \
  --threads 16
```

**Submit:**
```bash
sbatch submit.sh
squeue -u $USER
```

## Troubleshooting

- **Java error**: `conda install -y openjdk=17`
- **Memory error**: Reduce `--threads` or increase `#SBATCH --mem`
- **Conda not found**: `module load miniconda3` or `source ~/.bashrc`

## Complete Documentation

See [README.md](README.md) for:
- All pipeline parameters
- Output file descriptions
- Advanced usage examples
