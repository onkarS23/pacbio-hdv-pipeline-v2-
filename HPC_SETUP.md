# HPC Setup Guide

Quick guide for running the PacBio HDV pipeline on HPC clusters.

## Quick Start

### Step 1: Connect to HPC
```bash
ssh user@hpc-address
```

### Step 2: Setup Conda

**Option A: If conda module is available (check with your HPC admin)**
```bash
# List available conda modules
module avail conda
# OR
module avail miniconda

# Load the conda module (adjust name as shown above)
module load miniconda3
# OR
module load anaconda3
```

**Option B: If conda is not available, install Miniconda**
```bash
cd ~
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
# Follow prompts: yes to accept, yes to initialize
source ~/.bashrc

# Accept conda terms
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
```

### Step 3: Clone Repository
```bash
git clone https://github.com/onkarS23/pacbio-hdv-pipeline-v2-.git
cd pacbio-hdv-pipeline-v2-
```

### Step 4: Run Setup
```bash
bash setup_v2.sh
conda activate pacbio-hdv-v2
```

### Step 5: Download hg38 Reference
```bash
cd data
wget -O hg38.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
samtools faidx hg38.fa
cd ..
```

**Note:** Check with your HPC admin if hg38 is already available:
```bash
# Common locations - ask your admin
ls /data/references/hg38*
ls /shared/genomes/hg38*

# If found, create a symlink instead:
ln -s /path/to/existing/hg38.fa data/hg38.fa
ln -s /path/to/existing/hg38.fa.fai data/hg38.fa.fai
```

### Step 6: Prepare Your Data
```bash
# Copy your FASTQ files
cp /path/to/your/data/*.fastq test_input/

# Or create symlinks to save space
ln -s /path/to/your/data/*.fastq test_input/
```

### Step 7: Run Pipeline
```bash
nextflow run main_v2.nf \
  --input_dir test_input \
  --outdir results \
  --threads 16
```

---

## Running as Batch Job (Recommended for HPC)

Create `submit_pipeline.sh`:
```bash
#!/bin/bash
#SBATCH --job-name=hdv_pipeline
#SBATCH --output=hdv_%j.log
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00

# Load conda (if using module system)
module load miniconda3

# OR if installed in home directory, activate it
# source ~/miniconda3/etc/profile.d/conda.sh

conda activate pacbio-hdv-v2

nextflow run main_v2.nf \
  --input_dir test_input \
  --outdir results \
  --threads 16 \
  -with-report execution_report.html \
  -with-timeline timeline.html
```

Submit the job:
```bash
sbatch submit_pipeline.sh
```

Check job status:
```bash
squeue -u $USER
```

---

## Troubleshooting

**Conda command not found:**
- If using modules: `module load miniconda3`
- If installed locally: `source ~/.bashrc`

**Java not found:**
```bash
conda activate pacbio-hdv-v2
conda install -y openjdk=17
```

**Out of memory:**
- Reduce threads: `--threads 8`
- Request more memory in SLURM script: `#SBATCH --mem=128G`

**hg38 download too slow:**
- Ask your HPC admin if it's already available
- Use symlink instead of downloading

---

## Pipeline Parameters
```bash
nextflow run main_v2.nf \
  --input_dir test_input \      # Directory with FASTQ files
  --outdir results \             # Output directory
  --human_ref data/hg38.fa \     # Human reference
  --hdv_ref data/hdv_ref.fa \    # HDV reference (included)
  --threads 16                   # Number of threads
```

See [README.md](README.md) for complete documentation and all parameters.
