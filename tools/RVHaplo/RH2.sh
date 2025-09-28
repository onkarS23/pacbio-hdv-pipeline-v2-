#!/bin/bash
#SBATCH --job-name=rvhaplo  # Name of job
#SBATCH --output=/home/osingh/rvhaplo/RVHaplo/All_sample.out    # stdout
#SBATCH --error=/home/osingh/rvhaplo/RVHaplo/All_sample.err     # stderr
#SBATCH --partition=cpu      # partition to use (check with sinfo)
#SBATCH --nodes=1         # Number of nodes
#SBATCH --threads-per-core=1   # Ensure we only get one logical CPU per core
#SBATCH --cpus-per-task=32     # Number of cores per task
#SBATCH --mem=50G         # Memory per node | Alternative: --mem-per-cpu
#SBATCH --time=240:00:00      # wall time limit (HH:MM:SS)
#SBATCH --qos=long

# Array of sample input files
samples=(
    "JC_Sample_002_output.sam"
    "JC_Sample_003_output.sam"
    "JC_Sample_004_output.sam"
    "JC_Sample_005_output.sam"
    "JC_Sample_006_output.sam"
    "JC_Sample_007_output.sam"
    "JC_Sample_008_output.sam"
    "JC_Sample_009_output.sam"
    "JC_Sample_010_output.sam"
    "JC_Sample_011_output.sam"
    "JC_Sample_012_output.sam"
    "JC_Sample_013_output.sam"
    "JC_Sample_014_output.sam"
    "JC_Sample_015_output.sam"
    "JC_Sample_016_output.sam"

    # Add all other sample files here
)

# Reference genome path
reference_genome="/home/osingh/rvhaplo/RVHaplo/LT604957_GT5-reference_reoriented.fa"

# Loop through each sample
for sample in "${samples[@]}"; do
    # Extract the sample name without the extension for naming
    sample_name=$(basename "$sample" .sam)
    
    # Create unique output directories and filenames
    output_dir="/home/osingh/rvhaplo/RVHaplo/${sample_name}_result"
    stdout_file="${output_dir}/${sample_name}.out"
    stderr_file="${output_dir}/${sample_name}.err"
    
    # Create the output directory if it doesn't exist
    mkdir -p "$output_dir"
    
    # Set the SLURM output and error file paths
    sbatch --output="$stdout_file" --error="$stderr_file" <<EOT
#!/bin/bash
#SBATCH --job-name=rvhaplo_${sample_name}
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=50G
#SBATCH --time=240:00:00
#SBATCH --qos=long

# Run the script
./rvhaplo.sh -i /home/osingh/HDAV_ALL_DATA/New_1_Feb_All_data_files/JC_mm2_sam/${sample} -r $reference_genome -o ${output_dir} -p ${sample_name} -t 32 -sg 1 -wr 0.60 -wc 0.60
EOT

done


