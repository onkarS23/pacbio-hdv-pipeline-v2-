# PacBio HDV Analysis Pipeline v2.0

> **🆕 Version 2.0 Updates:**
> - ✅ Fixed medaka dependency (v1.11.0+ - resolves numpy conflicts)
> - ✅ Added OpenJDK 17 for Nextflow
> - ✅ HPC-ready with SLURM/PBS support
> - ✅ Tested and optimized for cluster environments
>
> **📘 HPC Users (Recommended):** See [HPC_SETUP.md](HPC_SETUP.md) for streamlined setup guide
>
> **💻 WSL/Local Users:** If running on Windows WSL or systems with limited RAM (<16GB), use `wsl.config` for reduced memory requirements

---

# PacBio HDV Analysis Pipeline v2.0

A comprehensive Nextflow pipeline for analyzing Hepatitis Delta Virus (HDV) from PacBio HiFi sequencing data, including haplotype reconstruction using RVHaplo.

## Features

- **Human read filtering**: Remove human contamination from PacBio reads
- **Read processing**: Quality filtering, adapter trimming, and primer removal  
- **HDV analysis**: Mapping to HDV reference and consensus generation
- **Haplotype reconstruction**: RVHaplo-based viral haplotype analysis (raw + polished)
- **Quality control**: Comprehensive QC analysis with NanoPlot and MultiQC
- **Coverage visualization**: Detailed coverage plots and statistics

## Quick Start Guide

### For new users cloning this repository:

```bash
# 1. Clone the repository
git clone <repository-url>
cd pacbio-hdv-pipeline-clean

# 2. Set up the environment (automated)
bash setup_v2.sh
conda activate pacbio-hdv-v2

# 3. Download large reference files (not included due to size)
# Create data directory if needed
mkdir -p data
# Download hg38 reference (required but not included due to 3GB size)
wget -O data/hg38.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip data/hg38.fa.gz
samtools faidx data/hg38.fa

# 4. Prepare your test data
# Place your PacBio FASTQ files in test_input/ directory
# (Original large test files removed - use your own data)

# 5. Optional: Run Quality Control first
nextflow run qc_only.nf --reads "test_input/*.fastq"

# 6. Run complete HDV analysis pipeline
nextflow run main_v2.nf \
  --input_dir test_input \
  --outdir hdv_results \
  --threads 16

# 7. Review results
# Check hdv_results/ for individual sample analysis
# Look for haplotype files in 04_rvhaplo_analysis/
```

### Alternative: Use the submission script

```bash
# Use the provided submission script
bash run_hdv_pipeline.sh
# or with custom parameters:
bash run_hdv_pipeline.sh test_input my_results 8
```

## Requirements

- **System**: Linux/macOS environment
- **Software**: 
  - Nextflow ≥ 21.0
  - Conda/Miniconda
  - Git LFS (if using large file storage)
- **Resources**: 
  - 16+ CPUs recommended
  - 32+ GB RAM
  - 100+ GB disk space for intermediate files
- **Reference genomes**: See [Reference Data Setup](#reference-data-setup)

## Installation

### Quick Setup (Recommended)
```bash
# Navigate to the pipeline directory
cd pacbio-hdv-pipeline-clean

# Run the automated setup script for v2
bash setup_v2.sh

# Activate the environment
conda activate pacbio-hdv-v2

# Verify installation
which nextflow
which medaka
which minimap2
which racon
```

### Manual Installation
```bash
# Create conda environment from the provided file
conda env create -f environment_v2.yml
conda activate pacbio-hdv-v2

# Install Nextflow if not available
curl -s https://get.nextflow.io | bash
mkdir -p $HOME/bin
mv nextflow $HOME/bin/
export PATH="$HOME/bin:$PATH"
```

## Reference Data Setup

### Included References
- `data/hdv_ref.fa` - HDV reference sequence (NC_001653, ~1.8KB)
- `data/hdv_refs/` - Additional HDV reference variants

### Missing Large References (Must Download)
Due to GitHub file size limits, you need to download:

```bash
# Download human reference genome (required, ~3.1GB)
cd data
wget -O hg38.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
samtools faidx hg38.fa

# Verify references are ready
ls -la data/
# Should see: hdv_ref.fa, hdv_ref.fa.fai, hg38.fa, hg38.fa.fai
```

## Test Data

### Provided Test Data
- `test_small.fastq` - Small test file (681KB) for quick testing

### Large Test Data (Not Included)
Original large test files were removed due to size. To test with realistic data:

1. **Use your own PacBio HiFi data**:
   ```bash
   # Place your FASTQ files in test_input/
   cp your_data/*.fastq test_input/
   ```

2. **Or download public datasets**:
   ```bash
   # Example: Download from SRA or other public repositories
   # Place in test_input/ directory
   ```

## Usage

### Main Pipeline (Recommended)

Run the complete HDV analysis pipeline:

```bash
# Basic usage
nextflow run main_v2.nf \
  --input_dir test_input \
  --outdir hdv_results

# Advanced usage with custom parameters
nextflow run main_v2.nf \
  --input_dir test_input \
  --outdir hdv_results \
  --human_ref data/hg38.fa \
  --hdv_ref data/hdv_ref.fa \
  --threads 16 \
  --min_length 500 \
  --min_coverage 5 \
  --alt_threshold 0.5
```

### Quality Control Only

For standalone QC analysis:

```bash
# Run QC analysis only
nextflow run qc_only.nf \
  --reads "test_input/*.fastq" \
  --outdir qc_results
```

### Using the Submission Script

The pipeline includes a submission script for easier execution:

```bash
# Default run
bash run_hdv_pipeline.sh

# Custom parameters
bash run_hdv_pipeline.sh [input_dir] [output_dir] [threads]

# Example
bash run_hdv_pipeline.sh test_input my_hdv_analysis 24
```

## Pipeline Parameters

### Required Parameters
| Parameter | Description | Example |
|-----------|-------------|---------|
| `--input_dir` | Directory containing FASTQ files | `test_input` |
| `--outdir` | Output directory | `hdv_results` |

### Optional Parameters
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--human_ref` | Path to human reference genome | `data/hg38.fa` |
| `--hdv_ref` | Path to HDV reference genome | `data/hdv_ref.fa` |
| `--threads` | Number of threads | `16` |
| `--min_length` | Minimum read length | `500` |
| `--min_coverage` | Minimum coverage for consensus | `5` |
| `--alt_threshold` | Alternative allele threshold | `0.5` |

## Output Structure

```
hdv_results/                    # Main output directory
├── sample_bc1002/             # Results for each sample
│   ├── 01_human_mapping/
│   │   ├── bc1002.unmapped.fastq
│   │   └── bc1002_human_mapping_stats.txt
│   ├── 02_processed_reads/
│   │   ├── cleaned_bc1002.fastq
│   │   └── bc1002_processing_stats.txt
│   ├── 03_hdv_analysis/
│   │   ├── bc1002.hdv.sam                    # Original SAM for RVHaplo
│   │   ├── bc1002.hdv.sorted.bam             # Sorted alignments
│   │   ├── bc1002.hdv.consensus.fa           # Pipeline consensus
│   │   ├── bc1002.hdv.coverage.txt           # Coverage data
│   │   ├── bc1002.original.coverage.png      # Coverage plots
│   │   ├── bc1002.consensus.coverage.png
│   │   └── bc1002_hdv_results.txt
│   └── 04_rvhaplo_analysis/
│       ├── bc1002_rvhaplo_results/           # Complete RVHaplo output
│       ├── bc1002_raw_haplotypes.fasta       # Raw haplotypes (unpolished)
│       ├── bc1002_polished_consensus.fasta   # Polished haplotypes (racon)
│       └── bc1002_rvhaplo_stats.txt
├── execution_report.html       # Nextflow execution report
├── timeline.html              # Pipeline timeline
└── pipeline_summary_report.html # Combined summary
```

## Key Output Files

### RVHaplo Haplotype Files
- **Raw haplotypes**: `*_raw_haplotypes.fasta` - Unpolished haplotype sequences
- **Polished consensus**: `*_polished_consensus.fasta` - Racon-polished final sequences

### Analysis Files
- **HDV consensus**: `*.hdv.consensus.fa` - Pipeline-generated consensus
- **Coverage plots**: `*.coverage.png` - Visual coverage analysis
- **Statistics**: `*_stats.txt` files - Detailed metrics for each step

## Expected Results

### Successful Analysis Indicators
1. **Human mapping**: 70-95% of reads typically map to human genome
2. **HDV detection**: Presence of HDV-mapped reads and consensus sequences
3. **Haplotype reconstruction**: Both raw and polished haplotype files generated
4. **Coverage analysis**: Adequate depth across HDV genome

### Typical Metrics
- **High-quality PacBio HiFi**: >95% reads above Q20, mean length 1-20kb
- **HDV coverage**: Variable, depends on viral load in sample
- **Processing efficiency**: ~90%+ reads retained after quality filtering

## Troubleshooting

### Common Issues

1. **Missing reference files**
   ```bash
   # Check if references exist and are indexed
   ls -la data/hg38.fa*
   ls -la data/hdv_ref.fa*
   ```

2. **Environment activation problems**
   ```bash
   # Ensure correct environment is active
   conda activate pacbio-hdv-v2
   which medaka  # Should show path in pacbio-hdv-v2 environment
   ```

3. **RVHaplo fails to produce haplotypes**
   - Check input SAM files have mapped reads
   - Verify racon is available in environment
   - Review RVHaplo output directory for error logs

4. **Memory issues**
   - Reduce thread count: `--threads 8`
   - Process samples individually
   - Monitor system resources during execution

5. **Large file handling**
   ```bash
   # Clean up intermediate files if disk space is low
   rm -rf work/  # Nextflow work directory
   rm -rf hdv_results/*/03_hdv_analysis/*.sam  # Large SAM files
   ```

### Getting Help
```bash
# Show pipeline help
nextflow run main_v2.nf --help

# Check environment setup
conda list -n pacbio-hdv-v2

# Test with small dataset
nextflow run main_v2.nf --input_dir . --outdir test_run
```

## Tools and Dependencies

### Core Analysis Tools
- **minimap2**: Long-read mapping to references
- **samtools**: SAM/BAM file processing and indexing
- **bedtools**: Genomic interval operations and coverage analysis
- **seqkit**: FASTA/FASTQ sequence processing
- **cutadapt**: Adapter and primer trimming

### Specialized Tools
- **medaka**: Consensus sequence polishing (for ONT support)
- **racon**: Consensus sequence polishing (PacBio optimized)
- **RVHaplo**: Viral haplotype reconstruction
- **MCL**: Graph clustering (used by RVHaplo)

### Quality Control
- **NanoPlot**: Long-read sequencing QC and visualization
- **MultiQC**: Report aggregation and summary

## Configuration

### Cluster Execution
For HPC environments, use the provided configuration:

```bash
# Slurm cluster
nextflow run main_v2.nf \
  -c nextflow.config \
  --input_dir test_input \
  --outdir hdv_results

# Custom resource allocation
nextflow run main_v2.nf \
  --input_dir test_input \
  --outdir hdv_results \
  --threads 24 \
  -with-report execution_report.html \
  -with-timeline timeline.html
```

### Resource Requirements
- **CPU**: 16+ cores recommended (adjustable with `--threads`)
- **Memory**: 32+ GB RAM for large datasets
- **Storage**: 10-20x input size for intermediate files
- **Runtime**: 2-6 hours depending on data size and complexity

## File Organization

```
pacbio-hdv-pipeline-clean/
├── main_v2.nf                 # Main pipeline script
├── qc_only.nf                 # QC-only pipeline
├── run_hdv_pipeline.sh        # Submission script
├── setup_v2.sh                # Environment setup
├── environment_v2.yml         # Conda environment specification
├── nextflow.config            # Pipeline configuration
├── data/                      # Reference genomes
│   ├── hdv_ref.fa             # HDV reference (included)
│   ├── hdv_refs/              # Additional HDV references
│   └── README.md              # Reference setup instructions
├── test_input/                # Test data directory
│   └── README.md              # Test data instructions
├── tools/                     # External tools
│   └── RVHaplo/               # RVHaplo installation
├── docs/                      # Documentation
└── examples/                  # Example configurations
```

## Citation

If you use this pipeline in your research, please cite:

- **RVHaplo**: Cai, D. et al. RVHaplo: Reconstructing viral haplotypes from long reads. *Bioinformatics* (2022).
- **minimap2**: Li, H. Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics* 34, 3094-3100 (2018).
- **medaka**: Oxford Nanopore Technologies. medaka: Sequence correction provided by ONT Research.
- **Nextflow**: Di Tommaso, P. et al. Nextflow enables reproducible computational workflows. *Nature Biotechnology* 35, 316-319 (2017).

## License

This pipeline is released under the MIT License. See LICENSE file for details.

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request with detailed description

## Support

For issues and questions:
1. Check this README and troubleshooting section
2. Review the pipeline logs and error messages
3. Submit an issue on the GitHub repository with:
   - Error messages
   - System information
   - Input data characteristics
   - Steps to reproduce the problem
