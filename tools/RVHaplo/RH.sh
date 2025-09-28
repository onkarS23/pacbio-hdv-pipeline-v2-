#!/bin/bash
#SBATCH --job-name=rvhaplo07  # Name of job
#SBATCH --output=/home/osingh/rvhaplo/RVHaplo/sample_007.out    # stdout
#SBATCH --error=/home/osingh/rvhaplo/RVHaplo/sample_007.err     # stderr
#SBATCH --partition=cpu      # partition to use (check with sinfo)
#SBATCH --nodes=1         # Number of nodes
#SBATCH --threads-per-core=1   # Ensure we only get one logical CPU per core
#SBATCH --cpus-per-task=32     # Number of cores per task
#SBATCH --mem=50G         # Memory per node | Alternative: --mem-per-cpu
#SBATCH --time=24:00:00      # wall time limit (HH:MM:SS)
#SBATCH --qos=long



# Run the script
./rvhaplo.sh -i /home/osingh/HDAV_ALL_DATA/New_1_Feb_All_data_files/JC_mm2_sam/JC_Sample_007_output.sam -r /home/osingh/rvhaplo/RVHaplo/HDV_JC-reference-genome.fa -o result_007 -p Sample-007_GT1_WR_WC_09 -t 32 -sg 5 -wr 0.9 -wc 0.9 -m 3 -l 20

