#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Show help message
def helpMessage() {
    log.info"""
    =========================================
     PacBio HDV Analysis Pipeline v2.0
    =========================================
    
    Usage:
      nextflow run main_v2.nf --input_dir /path/to/fastq --outdir results
    
    Required arguments:
      --input_dir       Directory containing PacBio FASTQ files
      --outdir          Output directory
    
    Optional arguments:
      --human_ref       Path to human reference genome (default: data/hg38.fa)
      --hdv_ref         Path to HDV reference genome (default: data/hdv_ref.fa)
      --threads         Number of threads (default: 16)
      --min_length      Minimum read length (default: 500)
      --min_coverage    Minimum coverage for consensus (default: 5)
      --alt_threshold   Alternative allele threshold (default: 0.5)
      --help            Show this help message
    """.stripIndent()
}

// Initialize parameters with defaults
params.help = false
params.human_ref = params.human_ref ?: "${projectDir}/data/hg38.fa"
params.hdv_ref = params.hdv_ref ?: "${projectDir}/data/hdv_ref.fa"
params.threads = params.threads ?: 16
params.min_length = params.min_length ?: 500
params.min_coverage = params.min_coverage ?: 5
params.alt_threshold = params.alt_threshold ?: 0.5
params.sam_filter = params.sam_filter ?: "-F 2304"
params.mp_min_mapq = params.mp_min_mapq ?: 20
params.mp_min_baseq = params.mp_min_baseq ?: 20

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Parameter validation
if (!params.input_dir) {
    log.error "Please specify --input_dir"
    exit 1
}

if (!params.outdir) {
    log.error "Please specify --outdir"
    exit 1
}

log.info """
=========================================
 PacBio HDV Analysis Pipeline v2.0
=========================================
Input directory    : ${params.input_dir}
Output directory   : ${params.outdir}
Human reference    : ${params.human_ref}
HDV reference      : ${params.hdv_ref}
Threads            : ${params.threads}
Min read length    : ${params.min_length}
Min coverage       : ${params.min_coverage}
Alt threshold      : ${params.alt_threshold}
=========================================
"""

// Create reference file channels
human_ref_ch = Channel.value(file(params.human_ref))
hdv_ref_ch = Channel.value(file(params.hdv_ref))

/*
 * STEP 1: Map to human genome and extract unmapped reads
 */
process HUMAN_MAPPING {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/01_human_mapping", mode: 'copy'
    
    input:
    tuple val(sample_id), path(fastq)
    path human_ref
    
    output:
    tuple val(sample_id), path("${sample_id}.unmapped.fastq"), emit: unmapped_reads
    path "${sample_id}_human_mapping_stats.txt", emit: stats
    
    script:
    """
    # Count total reads
    TOTAL_READS=\$(zcat -f ${fastq} | grep -c "^@")
    
    # Map to human genome
    minimap2 -a -x map-hifi -t ${params.threads} -k 19 -w 19 -O 4,24 -E 2,1 -A 1 -B 4 -z 400 -r 2000 -g 5000 \\
        --eqx ${human_ref} ${fastq} > ${sample_id}.human.sam 2> ${sample_id}.human.mapping.log
    
    # Filter and extract unmapped reads
    samtools view ${params.sam_filter} -b -f 4 ${sample_id}.human.sam > ${sample_id}.unmapped.bam
    samtools fastq ${sample_id}.unmapped.bam > ${sample_id}.unmapped.fastq
    
    # Calculate statistics using awk instead of bc
    HUMAN_MAPPED=\$(samtools view -c -F 260 ${sample_id}.human.sam)
    UNMAPPED_READS=\$(grep -c "^@" ${sample_id}.unmapped.fastq)
    HUMAN_PERCENT=\$(awk -v a=\$HUMAN_MAPPED -v b=\$TOTAL_READS 'BEGIN{ if(b==0) print 0; else printf "%.2f", a*100/b }')
    UNMAPPED_PERCENT=\$(awk -v a=\$UNMAPPED_READS -v b=\$TOTAL_READS 'BEGIN{ if(b==0) print 0; else printf "%.2f", a*100/b }')
    
    # Write stats
    cat > ${sample_id}_human_mapping_stats.txt << EOS
Sample: ${sample_id}
Total reads: \$TOTAL_READS
Human mapped reads: \$HUMAN_MAPPED (\$HUMAN_PERCENT%)
Unmapped reads: \$UNMAPPED_READS (\$UNMAPPED_PERCENT%)
EOS
    
    # Cleanup
    rm -f ${sample_id}.human.sam ${sample_id}.unmapped.bam
    """
}

/*
 * STEP 2: Filter and process reads
 */
process PROCESS_READS {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/02_processed_reads", mode: 'copy'
    
    input:
    tuple val(sample_id), path(unmapped_fastq)
    
    output:
    tuple val(sample_id), path("cleaned_${sample_id}.fastq"), emit: cleaned_reads
    path "${sample_id}_processing_stats.txt", emit: stats
    
    script:
    """
    # Step 1: Filter reads ≥ min_length
    seqkit seq -g -m ${params.min_length} ${unmapped_fastq} -o filtered_${sample_id}.fastq
    
    # Step 2: Primer trimming using Cutadapt
    cutadapt \\
        -g TGCCATGCCGACCCGAAGAGG...GATGCCCAGGTCGGACCGCGA \\
        -g TCGCGGTCCGACCTGGGCATC...CCTCTTCGGGTCGGCATGGCA \\
        --revcomp \\
        --times=1 \\
        --error-rate=0.05 \\
        --minimum-length ${params.min_length} \\
        --discard-untrimmed \\
        -o trimmed_${sample_id}.fastq \\
        filtered_${sample_id}.fastq > cutadapt_log.txt
    
    # Step 3: Detect residual primers
    grep -B1 -A2 -E 'TGCCATGCCGACCCGAAGAGG|GATGCCCAGGTCGGACCGCGA|TCGCGGTCCGACCTGGGCATC|CCTCTTCGGGTCGGCATGGCA' \\
        trimmed_${sample_id}.fastq > reads_with_primers.fastq || true
    
    # Step 4: Remove contaminated reads if needed
    if [ -s reads_with_primers.fastq ]; then
        cp reads_with_primers.fastq residual_primers_${sample_id}.fastq
        grep "^@" reads_with_primers.fastq | cut -d' ' -f1 | sed 's/^@//' > contaminated_ids.txt
        seqkit grep -v -f contaminated_ids.txt trimmed_${sample_id}.fastq -o cleaned_${sample_id}.fastq
        rm contaminated_ids.txt
        PRIMER_CONTAMINATED=\$(grep -c "^@" residual_primers_${sample_id}.fastq || echo 0)
        CLEANED_COUNT=\$(grep -c "^@" cleaned_${sample_id}.fastq)
    else
        cp trimmed_${sample_id}.fastq cleaned_${sample_id}.fastq
        PRIMER_CONTAMINATED=0
        CLEANED_COUNT=\$(grep -c "^@" cleaned_${sample_id}.fastq)
    fi
    
    # Step 5: Summary Stats
    RAW_COUNT=\$(grep -c "^@" ${unmapped_fastq})
    FILTERED_COUNT=\$(grep -c "^@" filtered_${sample_id}.fastq)
    TRIMMED_COUNT=\$(grep -c "^@" trimmed_${sample_id}.fastq)
    
    cat > ${sample_id}_processing_stats.txt << EOS
Sample: ${sample_id}
Total raw reads: \$RAW_COUNT
Reads ≥${params.min_length}bp: \$FILTERED_COUNT
Reads <${params.min_length}bp (filtered out): \$((RAW_COUNT - FILTERED_COUNT))
Reads matched and trimmed by cutadapt: \$TRIMMED_COUNT
Reads discarded by cutadapt (no primer): \$((FILTERED_COUNT - TRIMMED_COUNT))
Reads with residual primer (removed): \$PRIMER_CONTAMINATED
Final cleaned reads: \$CLEANED_COUNT
EOS
    
    # Cleanup intermediate files
    rm -f filtered_${sample_id}.fastq trimmed_${sample_id}.fastq reads_with_primers.fastq || true
    """
}

/*
 * STEP 3: Map to HDV and generate consensus
 */
process HDV_ANALYSIS {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/03_hdv_analysis", mode: 'copy'
    
    input:
    tuple val(sample_id), path(cleaned_fastq)
    path hdv_ref
    
    output:
    tuple val(sample_id), path("${sample_id}.hdv.sam"), emit: sam_files
    path "${sample_id}.hdv.sorted.bam*", emit: bam_files
    path "${sample_id}_hdv_results.txt", emit: results
    path "${sample_id}.hdv.consensus.fa", emit: consensus
    path "${sample_id}*.coverage.*", emit: coverage_files
    
    script:
    """
    # Create coverage plotting script
    cat > plot_coverage_simple.py << 'PYEOF'
#!/usr/bin/env python3
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

if len(sys.argv) != 3:
    sys.exit(1)

cov_file, out_png = sys.argv[1], sys.argv[2]
pos, cov = [], []
with open(cov_file) as f:
    for line in f:
      if not line.strip(): continue
      p = line.rstrip().split('\\t')
      if len(p) < 3: continue
      pos.append(int(p[1])); cov.append(int(p[2]))

plt.figure(figsize=(12,3))
plt.plot(pos, cov, linewidth=1)
plt.xlabel("Position (bp)"); plt.ylabel("Coverage")
plt.tight_layout(); plt.savefig(out_png, dpi=200, bbox_inches='tight')
PYEOF
    chmod +x plot_coverage_simple.py
    
    # Create indel-aware consensus fallback script
    cat > indel_consensus_fallback.py << 'PYEOF'
#!/usr/bin/env python3
import sys
from collections import Counter

if len(sys.argv) != 6:
    sys.exit(1)

ref_fa, pileup, out_fa, MIN_COV, MIN_FREQ = sys.argv[1:]
MIN_COV = int(MIN_COV); MIN_FREQ = float(MIN_FREQ)

ref_name = None
with open(ref_fa) as f:
  for line in f:
    if line.startswith(">") and ref_name is None:
      ref_name = line[1:].strip().split()[0]
      break

def parse_bases(bases, ref_base):
  subs, ins, dels = [], [], []
  i = 0; L = len(bases)
  while i < L:
    c = bases[i]
    if c in '.,': subs.append(ref_base)
    elif c in 'ACGTN': subs.append(c)
    elif c == '^': i += 1
    elif ord(c) in (36,42): pass  # \$ and *
    elif c == '+':
      i += 1; num=[]
      while i < L and bases[i].isdigit(): num.append(bases[i]); i+=1
      ln = int("".join(num)) if num else 0
      ins.append(bases[i:i+ln].upper()); i += ln-1
    elif c == '-':
      i += 1; num=[]
      while i < L and bases[i].isdigit(): num.append(bases[i]); i+=1
      ln = int("".join(num)) if num else 0
      dels.append(ln); i += ln-1
    i += 1
  return subs, ins, dels

cons = []
with open(pileup) as f:
  for line in f:
    if not line.strip(): continue
    parts = line.rstrip().split('\\t')
    if len(parts) < 5: continue
    _, _, ref_base, depth, bases = parts[:5]
    depth = int(depth); ref_base = ref_base.upper()
    if depth < MIN_COV: cons.append('N'); continue
    subs, ins_list, dels = parse_bases(bases, ref_base)
    if len(subs) < MIN_COV: cons.append('N'); continue

    base, n = Counter(subs).most_common(1)[0]
    chosen = base if (n/len(subs)) >= MIN_FREQ else ref_base

    del_len = 0
    if dels:
      dl, dn = Counter(dels).most_common(1)[0]
      if dl > 0 and dn/len(subs) >= MIN_FREQ: del_len = dl

    ins_seq = ""
    if ins_list:
      ins_choice, inum = Counter(ins_list).most_common(1)[0]
      if inum/len(subs) >= MIN_FREQ: ins_seq = ins_choice

    if del_len > 0:
      if ins_seq: cons.append(ins_seq)
    else:
      cons.append(chosen)
      if ins_seq: cons.append(ins_seq)

s = "".join(cons)
with open(out_fa, "w") as w:
  w.write(f">{ref_name}_consensus\\n")
  for i in range(0, len(s), 80): w.write(s[i:i+80] + "\\n")
PYEOF
    chmod +x indel_consensus_fallback.py
    
    # Ensure HDV reference has faidx index
    if [ ! -f "${hdv_ref}.fai" ]; then
        samtools faidx ${hdv_ref}
    fi
    
    # Step 1: Map to HDV reference and KEEP the SAM file
    minimap2 -a -x map-hifi -t ${params.threads} -k 19 -w 19 -O 4,24 -E 2,1 -A 1 -B 4 -z 400 -r 2000 -g 5000 \\
        ${hdv_ref} ${cleaned_fastq} > ${sample_id}.hdv.sam 2> ${sample_id}.hdv.mapping.log
    
    # Step 2: BAM processing (create from SAM copy)
    samtools view ${params.sam_filter} -b ${sample_id}.hdv.sam > ${sample_id}.hdv.bam
    samtools sort ${sample_id}.hdv.bam -o ${sample_id}.hdv.sorted.bam
    samtools index ${sample_id}.hdv.sorted.bam
    
    # Step 3: Coverage analysis
    bedtools genomecov -d -ibam ${sample_id}.hdv.sorted.bam > ${sample_id}.hdv.coverage.txt
    python3 plot_coverage_simple.py ${sample_id}.hdv.coverage.txt ${sample_id}.original.coverage.png || true
    
    # Step 4: Consensus generation
    samtools mpileup -B -f ${hdv_ref} -d 8000 -A -Q ${params.mp_min_baseq} -q ${params.mp_min_mapq} -a ${sample_id}.hdv.sorted.bam \\
        > ${sample_id}.hdv.pileup
    python3 indel_consensus_fallback.py \\
        ${hdv_ref} ${sample_id}.hdv.pileup ${sample_id}.hdv.consensus.fa ${params.min_coverage} ${params.alt_threshold}
    
    # Validate consensus file
    if [ ! -s ${sample_id}.hdv.consensus.fa ]; then
        echo "Consensus FASTA is empty"; exit 1
    fi
    
    # Update header
    sed -i "1s/^>.*/>$sample_id\\_consensus/" ${sample_id}.hdv.consensus.fa
    samtools faidx ${sample_id}.hdv.consensus.fa
    
    # Step 5: Remap to consensus
    minimap2 -a -x map-hifi -t ${params.threads} -k 19 -w 19 -O 4,24 -E 2,1 -A 1 -B 4 -z 400 -r 2000 -g 5000 \\
        ${sample_id}.hdv.consensus.fa ${cleaned_fastq} > ${sample_id}.hdv.consensus.sam 2>> ${sample_id}.hdv.mapping.log
    
    samtools view ${params.sam_filter} -b ${sample_id}.hdv.consensus.sam > ${sample_id}.hdv.consensus.bam
    samtools sort ${sample_id}.hdv.consensus.bam -o ${sample_id}.hdv.consensus.sorted.bam
    samtools index ${sample_id}.hdv.consensus.sorted.bam
    
    # Step 6: Consensus coverage
    bedtools genomecov -d -ibam ${sample_id}.hdv.consensus.sorted.bam > ${sample_id}.hdv.consensus.coverage.txt
    python3 plot_coverage_simple.py ${sample_id}.hdv.consensus.coverage.txt ${sample_id}.consensus.coverage.png || true
    
    # Step 7: Final statistics using awk instead of bc
    TOTAL_READS=\$(grep -c "^@" ${cleaned_fastq})
    HDV_MAPPED=\$(samtools view -c -F 260 ${sample_id}.hdv.sorted.bam)
    HDV_CONS_MAPPED=\$(samtools view -c -F 260 ${sample_id}.hdv.consensus.sorted.bam)
    HDV_PERCENT=\$(awk -v a=\$HDV_MAPPED -v b=\$TOTAL_READS 'BEGIN{ if(b==0) print 0; else printf "%.2f", a*100/b }')
    HDV_CONS_PERCENT=\$(awk -v a=\$HDV_CONS_MAPPED -v b=\$TOTAL_READS 'BEGIN{ if(b==0) print 0; else printf "%.2f", a*100/b }')
    
    cat > ${sample_id}_hdv_results.txt << EOS
Sample: ${sample_id}
Total cleaned reads: \$TOTAL_READS
HDV mapped reads (original): \$HDV_MAPPED (\$HDV_PERCENT%)
HDV mapped reads (consensus): \$HDV_CONS_MAPPED (\$HDV_CONS_PERCENT%)
Original SAM file: ${sample_id}.hdv.sam
Original coverage: ${sample_id}.hdv.coverage.txt
Consensus coverage: ${sample_id}.hdv.consensus.coverage.txt
Consensus FASTA: ${sample_id}.hdv.consensus.fa
EOS
    
    # Cleanup intermediate files but KEEP the original SAM file
    rm -f ${sample_id}.hdv.bam ${sample_id}.hdv.consensus.sam ${sample_id}.hdv.consensus.bam ${sample_id}.hdv.pileup || true
    """
}

/*
 * STEP 4: RVHaplo Analysis
 */
process RVHAPLO {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/04_rvhaplo_analysis", mode: 'copy'
    
    input:
    tuple val(sample_id), path(sam_file)
    path hdv_ref

    output:
    path "${sample_id}_rvhaplo_results/", emit: rvhaplo_dir, optional: true
    path "${sample_id}_raw_haplotypes.fasta", emit: raw_haplotypes
    path "${sample_id}_polished_consensus.fasta", emit: polished_consensus
    path "${sample_id}_rvhaplo_stats.txt", emit: stats

    script:
    """
    # Check if SAM file exists and has content
    if [ ! -s "${sam_file}" ]; then
        echo "SAM file is empty - no reads for haplotype reconstruction" > ${sample_id}_rvhaplo_stats.txt
        touch ${sample_id}_raw_haplotypes.fasta
        touch ${sample_id}_polished_consensus.fasta
        exit 0
    fi
    
    MAPPED_COUNT=\$(samtools view -c -F 4 ${sam_file})
    if [ \$MAPPED_COUNT -eq 0 ]; then
        echo "No mapped reads in SAM file - skipping haplotype reconstruction" > ${sample_id}_rvhaplo_stats.txt
        touch ${sample_id}_raw_haplotypes.fasta
        touch ${sample_id}_polished_consensus.fasta
        exit 0
    fi

    # Use the SAM file directly (it should already be named correctly)
    SAM_FILE="${sam_file}"
    
    # Set up absolute paths for RVHaplo
    ABS_IN=\$(readlink -f "\$SAM_FILE")
    ABS_REF=\$(readlink -f "${hdv_ref}")
    ABS_OUT=\$(readlink -f "${sample_id}_rvhaplo_results")
    mkdir -p "\$ABS_OUT"

    # Run RVHaplo from its directory
    cd ${projectDir}/tools/RVHaplo
    bash ./rvhaplo.sh \\
      -i "\$ABS_IN" \\
      -r "\$ABS_REF" \\
      -o "\$ABS_OUT" \\
      -p "${sample_id}_GT1_WR_WC_95" \\
      -t 8 \\
      -sg 5 \\
      -wr 0.95 \\
      -wc 0.95 \\
      -m 3 \\
      -l 20 \\
      --polisher racon || echo "RVHaplo execution completed with warnings"

    # Return to work directory
    cd - > /dev/null

    # Look for both raw haplotypes and polished consensus files
    RAW_HAPLOTYPES="${sample_id}_rvhaplo_results/${sample_id}_GT1_WR_WC_95_haplotypes.fasta"
    POLISHED_CONSENSUS="${sample_id}_rvhaplo_results/${sample_id}_GT1_WR_WC_95_consensus.fasta"
    
    # Alternative naming patterns RVHaplo might use
    if [ ! -f "\$RAW_HAPLOTYPES" ]; then
        RAW_HAPLOTYPES=\$(find ${sample_id}_rvhaplo_results/ -name "*haplotypes*.fasta" -o -name "*haplotype*.fasta" | head -1)
    fi
    
    if [ ! -f "\$POLISHED_CONSENSUS" ]; then
        POLISHED_CONSENSUS=\$(find ${sample_id}_rvhaplo_results/ -name "*consensus*.fasta" -o -name "*polished*.fasta" | head -1)
    fi
    
    # Copy raw haplotypes if found
    if [ -f "\$RAW_HAPLOTYPES" ] && [ -s "\$RAW_HAPLOTYPES" ]; then
        cp "\$RAW_HAPLOTYPES" ${sample_id}_raw_haplotypes.fasta
        RAW_COUNT=\$(grep -c "^>" ${sample_id}_raw_haplotypes.fasta || echo 0)
    else
        touch ${sample_id}_raw_haplotypes.fasta
        RAW_COUNT=0
    fi
    
    # Copy polished consensus if found
    if [ -f "\$POLISHED_CONSENSUS" ] && [ -s "\$POLISHED_CONSENSUS" ]; then
        cp "\$POLISHED_CONSENSUS" ${sample_id}_polished_consensus.fasta
        POLISHED_COUNT=\$(grep -c "^>" ${sample_id}_polished_consensus.fasta || echo 0)
    else
        touch ${sample_id}_polished_consensus.fasta
        POLISHED_COUNT=0
    fi
    
    # Generate comprehensive stats
    cat > ${sample_id}_rvhaplo_stats.txt << EOF
RVHaplo Analysis Results
Sample: ${sample_id}
Input SAM: \$SAM_FILE
Input mapped reads: \$MAPPED_COUNT

Raw haplotypes file: \$RAW_HAPLOTYPES
Raw haplotypes count: \$RAW_COUNT

Polished consensus file: \$POLISHED_CONSENSUS
Polished haplotypes count: \$POLISHED_COUNT

Status: \$([ \$RAW_COUNT -gt 0 ] || [ \$POLISHED_COUNT -gt 0 ] && echo "Success" || echo "No haplotypes reconstructed")
EOF

    # Add haplotype details if available
    if [ \$RAW_COUNT -gt 0 ]; then
        echo "" >> ${sample_id}_rvhaplo_stats.txt
        echo "Raw haplotype headers:" >> ${sample_id}_rvhaplo_stats.txt
        grep "^>" ${sample_id}_raw_haplotypes.fasta >> ${sample_id}_rvhaplo_stats.txt
    fi
    
    if [ \$POLISHED_COUNT -gt 0 ]; then
        echo "" >> ${sample_id}_rvhaplo_stats.txt
        echo "Polished consensus headers:" >> ${sample_id}_rvhaplo_stats.txt
        grep "^>" ${sample_id}_polished_consensus.fasta >> ${sample_id}_rvhaplo_stats.txt
    fi
    """
}

/*
 * STEP 5: Generate comprehensive summary report
 */
process GENERATE_REPORT {
    publishDir "${params.outdir}/", mode: 'copy'
    
    input:
    path human_stats
    path processing_stats
    path hdv_results
    path rvhaplo_stats
    
    output:
    path "pipeline_summary_report.html", emit: report
    path "pipeline_stats.txt", emit: summary
    
    script:
    """
    # Combine all statistics
    echo "======================================" > pipeline_stats.txt
    echo " PacBio HDV Pipeline Summary Report" >> pipeline_stats.txt
    echo "======================================" >> pipeline_stats.txt
    echo "Generated: \$(date)" >> pipeline_stats.txt
    echo "" >> pipeline_stats.txt
    
    for file in ${human_stats}; do
        echo "--- Human Mapping ---" >> pipeline_stats.txt
        cat "\$file" >> pipeline_stats.txt
        echo "" >> pipeline_stats.txt
    done
    
    for file in ${processing_stats}; do
        echo "--- Read Processing ---" >> pipeline_stats.txt
        cat "\$file" >> pipeline_stats.txt
        echo "" >> pipeline_stats.txt
    done
    
    for file in ${hdv_results}; do
        echo "--- HDV Analysis ---" >> pipeline_stats.txt
        cat "\$file" >> pipeline_stats.txt
        echo "" >> pipeline_stats.txt
    done
    
    for file in ${rvhaplo_stats}; do
        echo "--- RVHaplo Analysis ---" >> pipeline_stats.txt
        cat "\$file" >> pipeline_stats.txt
        echo "" >> pipeline_stats.txt
    done
    
    # Create comprehensive HTML report
    cat > pipeline_summary_report.html << EOF
<!DOCTYPE html>
<html>
<head>
    <title>PacBio HDV Pipeline Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        h1 { color: #2E8B57; }
        h2 { color: #4682B4; border-bottom: 1px solid #ccc; }
        pre { background-color: #f5f5f5; padding: 10px; border-radius: 5px; }
        .section { margin-bottom: 30px; }
        .highlight { background-color: #ffffcc; padding: 5px; }
    </style>
</head>
<body>
    <h1>PacBio HDV Analysis Pipeline Report</h1>
    <p><strong>Generated:</strong> \$(date)</p>
    
    <div class="section">
        <h2>Pipeline Summary</h2>
        <div class="highlight">
            <p><strong>Key Outputs:</strong></p>
            <ul>
                <li>Original HDV alignments (SAM files) for RVHaplo</li>
                <li>HDV consensus sequences</li>
                <li>Coverage plots and statistics</li>
                <li>Haplotype reconstruction results</li>
            </ul>
        </div>
        <pre>
\$(cat pipeline_stats.txt)
        </pre>
    </div>
</body>
</html>
EOF
    """
}

/*
 * Main workflow
 */
workflow {
    // Collect FASTQ files with flexible naming
    fastq_ch = Channel
        .fromPath("${params.input_dir}/*.{fastq,fastq.gz}")
        .map { file -> 
            def sample_id = file.baseName.replaceAll(/\\.fastq.*\$/, '')
            // Try to extract barcode if present (bc1002, bc1003, etc.)
            def barcode = file.name.find(/bc\\d+/)
            if (barcode) {
                sample_id = barcode
            }
            tuple(sample_id, file)
        }
    
    // Run pipeline steps
    human_results = HUMAN_MAPPING(fastq_ch, human_ref_ch)
    processed_results = PROCESS_READS(human_results.unmapped_reads)
    hdv_results = HDV_ANALYSIS(processed_results.cleaned_reads, hdv_ref_ch)
    
    // Run RVHaplo with the original SAM files
    rvhaplo_results = RVHAPLO(hdv_results.sam_files, hdv_ref_ch)
    
    // Generate comprehensive report
    GENERATE_REPORT(
        human_results.stats.collect(),
        processed_results.stats.collect(),
        hdv_results.results.collect(),
        rvhaplo_results.stats.collect()
    )
}

/*
 * Workflow completion message
 */
workflow.onComplete {
    log.info """
    =========================================
    Pipeline completed!
    =========================================
    Results directory: ${params.outdir}
    Success: ${workflow.success}
    Duration: ${workflow.duration}
    =========================================
    
    Key Outputs:
    - Human mapping statistics
    - Processed read statistics  
    - HDV consensus sequences
    - Coverage plots (PNG format)
    - Original SAM files (preserved for RVHaplo)
    - RVHaplo haplotype reconstruction
    - Comprehensive HTML summary report
    =========================================
    """
}