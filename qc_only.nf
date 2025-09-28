nextflow.enable.dsl=2

// ---------- Defaults (override via CLI) ----------
params.qc_env    = params.qc_env    ?: "${projectDir}/.envs/nanoplot-env"
params.reads     = params.reads     ?: 'data/*.{fastq,fq,fastq.gz,fq.gz}'
params.outdir    = params.outdir    ?: 'production_results'
params.threads   = params.threads   ?: 4

// ---------- Per-sample NanoPlot QC ----------
process NANOPLOT {
    tag "${file(reads).getBaseName()}"
    publishDir "${params.outdir}/nanoplot", mode: 'copy'
    conda "${params.qc_env}"
    
    input:
    path reads
    
    output:
    path "nanoplot_*", emit: report_dir
    
    script:
    def bn = reads.getBaseName()
    def outdir = "nanoplot_${bn}"
    """
    set -euo pipefail
    
    echo "=== NanoPlot QC Analysis ===" >&2
    echo "Input: ${reads}" >&2
    echo "Output dir: ${outdir}" >&2
    
    # Run NanoPlot with comprehensive options
    NanoPlot \\
        --fastq ${reads} \\
        --outdir ${outdir} \\
        --prefix ${bn}_ \\
        --threads ${params.threads} \\
        --N50 \\
        --tsv_stats \\
        --info_in_report \\
        --plots hex dot kde \\
        --format png pdf svg
    
    echo "NanoPlot analysis completed for ${reads}" >&2
    ls -la ${outdir}/ >&2
    """
}

// ---------- Aggregate with MultiQC ----------
process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    conda "${params.qc_env}"
    
    input:
    path nanoplot_dirs
    
    output:
    path "multiqc_report.html", optional: true
    path "multiqc_data", optional: true
    path "summary_report.txt"
    
    script:
    """
    set -euo pipefail
    
    echo "=== MultiQC Report Generation ===" >&2
    
    # Create comprehensive summary report
    echo "PacBio HiFi QC Analysis Summary" > summary_report.txt
    echo "===============================" >> summary_report.txt
    echo "Date: \$(date)" >> summary_report.txt
    echo "Analysis tool: NanoPlot" >> summary_report.txt
    echo "" >> summary_report.txt
    
    # Process each NanoPlot output
    for dir in ${nanoplot_dirs.join(' ')}; do
        if [ -d "\$dir" ]; then
            SAMPLE=\$(basename "\$dir" | sed 's/^nanoplot_//')
            echo "Sample: \$SAMPLE" >> summary_report.txt
            
            # Extract key statistics from NanoPlot output
            if [ -f "\$dir"/*_NanoStats.txt ]; then
                echo "  Status: Complete" >> summary_report.txt
                
                # Extract key metrics
                TOTAL_READS=\$(grep "Number of reads" "\$dir"/*_NanoStats.txt | cut -f2 || echo "N/A")
                TOTAL_BASES=\$(grep "Total bases" "\$dir"/*_NanoStats.txt | cut -f2 || echo "N/A")
                MEAN_LENGTH=\$(grep "Mean read length" "\$dir"/*_NanoStats.txt | cut -f2 || echo "N/A")
                N50=\$(grep "Read length N50" "\$dir"/*_NanoStats.txt | cut -f2 || echo "N/A")
                MEAN_QUALITY=\$(grep "Mean read quality" "\$dir"/*_NanoStats.txt | cut -f2 || echo "N/A")
                
                echo "  Total reads: \$TOTAL_READS" >> summary_report.txt
                echo "  Total bases: \$TOTAL_BASES" >> summary_report.txt
                echo "  Mean read length: \$MEAN_LENGTH bp" >> summary_report.txt
                echo "  N50: \$N50 bp" >> summary_report.txt
                echo "  Mean quality: \$MEAN_QUALITY" >> summary_report.txt
            else
                echo "  Status: No statistics file found" >> summary_report.txt
            fi
            echo "" >> summary_report.txt
        fi
    done
    
    # Try to run MultiQC (it may recognize NanoPlot outputs)
    if command -v multiqc >/dev/null 2>&1; then
        if multiqc ${nanoplot_dirs.join(' ')} -o . --force --ignore "**/tmp" 2>/dev/null; then
            echo "MultiQC: Successfully integrated NanoPlot results" >> summary_report.txt
        else
            echo "MultiQC: Could not process NanoPlot outputs, summary report generated instead" >> summary_report.txt
        fi
    else
        echo "MultiQC: Not available, summary report generated" >> summary_report.txt
    fi
    
    # Create basic HTML report if MultiQC didn't work
    if [ ! -f "multiqc_report.html" ]; then
        cat > multiqc_report.html << 'EOF'
<!DOCTYPE html>
<html>
<head>
    <title>PacBio HiFi QC Results</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        h1, h2 { color: #2c3e50; }
        .summary { background-color: #f8f9fa; padding: 15px; border-radius: 5px; }
        .metric { margin: 5px 0; }
    </style>
</head>
<body>
    <h1>PacBio HiFi QC Analysis Results</h1>
    <div class="summary">
        <h2>Analysis Summary</h2>
        <p>Quality control analysis completed using NanoPlot.</p>
        <p>Individual sample reports are available in the nanoplot output directories.</p>
        <p>See summary_report.txt for detailed statistics.</p>
    </div>
    
    <h2>Individual Reports</h2>
    <p>Each sample has a dedicated NanoPlot report with:</p>
    <ul>
        <li>Read length distribution plots</li>
        <li>Quality score distributions</li>
        <li>Yield over time analysis</li>
        <li>Comprehensive statistics</li>
        <li>Publication-quality figures</li>
    </ul>
</body>
</html>
EOF
        mkdir -p multiqc_data
        echo "NanoPlot QC analysis completed" > multiqc_data/summary.txt
    fi
    
    echo "Report generation completed" >&2
    """
}

// ---------- Workflow ----------
workflow {
    log.info """
    =====================================
    PacBio HiFi QC Pipeline (NanoPlot)
    =====================================
    reads      : ${params.reads}
    outdir     : ${params.outdir}
    qc_env     : ${params.qc_env}
    threads    : ${params.threads}
    =====================================
    """
    
    reads_ch = Channel.fromPath(params.reads, checkIfExists: true)
                     .ifEmpty { error "No FASTQ files matched: ${params.reads}" }
    
    reads_ch.view { "Processing: $it" }
    
    nanoplot_out = NANOPLOT(reads_ch)
    MULTIQC(nanoplot_out.report_dir.collect())
}