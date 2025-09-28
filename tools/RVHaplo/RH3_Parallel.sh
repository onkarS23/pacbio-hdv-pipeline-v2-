#!/bin/bash

# List of input files
input_files=(
    "/home/osingh/HDAV_ALL_DATA/New_1_Feb_All_data_files/JC_mm2_sam/JC_Sample_002_output.sam"
    "/home/osingh/HDAV_ALL_DATA/New_1_Feb_All_data_files/JC_mm2_sam/JC_Sample_003_output.sam"
    "/home/osingh/HDAV_ALL_DATA/New_1_Feb_All_data_files/JC_mm2_sam/JC_Sample_004_output.sam"
    "/home/osingh/HDAV_ALL_DATA/New_1_Feb_All_data_files/JC_mm2_sam/JC_Sample_005_output.sam"
    "/home/osingh/HDAV_ALL_DATA/New_1_Feb_All_data_files/JC_mm2_sam/JC_Sample_006_output.sam"
    "/home/osingh/HDAV_ALL_DATA/New_1_Feb_All_data_files/JC_mm2_sam/JC_Sample_007_output.sam"
    "/home/osingh/HDAV_ALL_DATA/New_1_Feb_All_data_files/JC_mm2_sam/JC_Sample_008_output.sam"
    "/home/osingh/HDAV_ALL_DATA/New_1_Feb_All_data_files/JC_mm2_sam/JC_Sample_009_output.sam"
    "/home/osingh/HDAV_ALL_DATA/New_1_Feb_All_data_files/JC_mm2_sam/JC_Sample_010_output.sam"
    "/home/osingh/HDAV_ALL_DATA/New_1_Feb_All_data_files/JC_mm2_sam/JC_Sample_011_output.sam"
    "/home/osingh/HDAV_ALL_DATA/New_1_Feb_All_data_files/JC_mm2_sam/JC_Sample_012_output.sam"
    "/home/osingh/HDAV_ALL_DATA/New_1_Feb_All_data_files/JC_mm2_sam/JC_Sample_013_output.sam"
    "/home/osingh/HDAV_ALL_DATA/New_1_Feb_All_data_files/JC_mm2_sam/JC_Sample_014_output.sam"
    "/home/osingh/HDAV_ALL_DATA/New_1_Feb_All_data_files/JC_mm2_sam/JC_Sample_015_output.sam"
    "/home/osingh/HDAV_ALL_DATA/New_1_Feb_All_data_files/JC_mm2_sam/JC_Sample_016_output.sam"
    
)

# Common parameters
reference="/home/osingh/rvhaplo/RVHaplo/LT604957_GT5-reference_reoriented.fa"
base_output_dir="/home/osingh/rvhaplo/RVHaplo/results"
threads=32

# Loop over each input file and submit a job
for input_file in "${input_files[@]}"; do
    sample_name=$(basename "$input_file" .sam) # Extract sample name from input file name
    job_name="rvhaplo_${sample_name}"          # Create a unique job name
    sample_output_dir="${base_output_dir}/${sample_name}" # Create a unique output directory for each sample

    # Create the sample output directory if it doesn't exist
    mkdir -p "$sample_output_dir"

    sbatch --job-name="$job_name" \
           --output="${sample_output_dir}/${sample_name}.out" \
           --error="${sample_output_dir}/${sample_name}.err" \
           --partition=cpu \
           --nodes=1 \
           --threads-per-core=1 \
           --cpus-per-task=$threads \
           --mem=50G \
           --time=240:00:00 \
           --qos=long \
           --wrap="./rvhaplo.sh -i $input_file -r $reference -o $sample_output_dir -p ${sample_name}_NP -t $threads -sg 5 -wr 0.5 -wc 0.5 -m 3 -l 20"
done
