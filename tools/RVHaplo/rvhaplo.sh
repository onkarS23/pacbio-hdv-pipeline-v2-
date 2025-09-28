#! /bin/bash
set -euo pipefail

### required arguments
file_sam=""
file_ref=""

### optional arguments
file_path='./result'
prefix="rvhaplo"
mq=0
thread=32            # fixed stray backticks
error_rate=0.1
signi_level=0.05
cond_pro=0.65
fre_snv=0.8
num_read_1=10
num_read_2=5
gap=15
smallest_snv=20
only_snv=0
ovlap_read=5
weight_read=0.85
mcl_inflation=2
lar_cluster=50
ovlap_cluster=10
depth=5
weight_cluster=0.8
abundance=0.005
s_pos=1
e_pos=10000000000
sub_graph=1
polisher="medaka"    # NEW: 'medaka' (legacy default) or 'racon'

# ---- helper: Racon polish per-cluster (2 rounds), writes *_haplotypes.fasta
racon_polish_clusters () {
  local base_path="$1"    # e.g. $file_path
  local pfx="$2"          # e.g. $prefix
  local clusters_dir="${base_path}/clusters"
  local polished_all="${base_path}/${pfx}_haplotypes.fasta"
  : > "$polished_all"

  shopt -s nullglob
  for fq in "${clusters_dir}"/cluster_*.fastq; do
    local base
    base=$(basename "$fq" .fastq)

    # Try per-cluster draft first; fall back to global consensus
    local draft="${clusters_dir}/${base}_consensus.fasta"
    [[ -s "$draft" ]] || draft="${base_path}/${pfx}_consensus.fasta"
    if [[ ! -s "$draft" ]]; then
      echo "WARN: no draft consensus for ${base}; skipping racon polish." >&2
      continue
    fi

    # Round 1
    minimap2 -x map-hifi -a "$draft" "$fq" | samtools sort -o "${clusters_dir}/${base}.r1.bam"
    samtools index "${clusters_dir}/${base}.r1.bam"
    racon -u "$fq" "${clusters_dir}/${base}.r1.bam" "$draft" > "${clusters_dir}/${base}.p1.fa"

    # Round 2
    minimap2 -x map-hifi -a "${clusters_dir}/${base}.p1.fa" "$fq" | samtools sort -o "${clusters_dir}/${base}.r2.bam"
    samtools index "${clusters_dir}/${base}.r2.bam"
    racon -u "$fq" "${clusters_dir}/${base}.r2.bam" "${clusters_dir}/${base}.p1.fa" > "${clusters_dir}/${base}.polished.fa"

    cat "${clusters_dir}/${base}.polished.fa" >> "$polished_all"
  done
  shopt -u nullglob
  [[ -s "$polished_all" ]] || : > "$polished_all"
}

function help_info() {
	echo "Usage: $0 -i alignment.sam -r ref_genome.fasta [options]"
	echo ""
	echo "RVHaplo: Reconstructing viral haplotypes using long reads"
	echo ""
	echo "Author: Dehan CAI"
	echo "Date:   May 2022"
	echo "Version 2: Support multi-thread processing; Use a C package of MCL; Cost less memory"
	echo ""
	echo "    -i  | --input                    alignment file (SAM)"
	echo "    -r  | --ref_genome               reference genome (FASTA)"
	echo ""
	echo "    Options:"
	echo "    -o  | --out                      Path where to output results (default: ./result)"
	echo "    -p  | --prefix STR               Prefix of output files (default: rvhaplo)"
	echo "    -t  | --thread INT               CPU threads (default: 32)"
	echo "    -e  | --error_rate FLOAT         Sequencing error rate (default: 0.1)"
	echo "    -mq | --map_qual INT             Min mapping quality (default: 0)"
	echo "    -s  | --signi_level FLOAT        Significance level for binomial tests (default: 0.05)"
	echo "    -c  | --cond_pro FLOAT           Max conditional prob threshold for SNV (default: 0.65)"
	echo "    -f  | --fre_snv FLOAT            Dominant base frequency threshold (default: 0.80)"
	echo "    -n1 | --num_read_1 INT           Min reads for p(ai|aj) (default: 10)"
	echo "    -n2 | --num_read_2 INT           Min reads for p(ai|aj1..ajp) (default: 5)"
	echo "    -g  | --gap INT                  Min gap between SNVs (default: 15)"
	echo "    -ss | --smallest_snv INT         Min #SNVs for reconstruction (default: 20)"
	echo "    -os | --only_snv 0|1             Only output SNVs, skip reconstruction (default: 0)"
	echo "    -or | --ovlap_read INT           Min overlap between reads (default: 5)"
	echo "    -wr | --weight_read FLOAT        Min edge weight in read graph (default: 0.85)"
	echo "    -sg | --sub_graph INT            #subgraphs for MCL (default: 1)"
	echo "    -m  | --mcl_inflation FLOAT      MCL inflation (default: 2)"
	echo "    -l  | --lar_cluster INT          Large-cluster size threshold (default: 50)"
	echo "    -oc | --ovlap_cluster INT        Min overlap between cluster consensuses (default: 10)"
	echo "    -d  | --depth INT                Depth limit for cluster consensus (default: 5)"
	echo "    -wc | --weight_cluster FLOAT     Min cluster weight (default: 0.8)"
	echo "    -sp | --start_pos INT            Start position (default: 1)"
	echo "    -ep | --end_pos INT              End position (default: 1e10)"
	echo "    -a  | --abundance FLOAT          Min haplotype abundance (default: 0.005)"
	echo "    -P  | --polisher STR             'medaka' (ONT) or 'racon' (PacBio); default: medaka"
	echo "    -h  | --help                     Print help message"
	echo ""
	echo "    For details, see https://github.com/dhcai21/RVHaplo"
	echo ""
	exit 0
}

[[ "${1:-}" == "" ]] && { help_info; exit 1; }

while [[ "$#" -gt 0 ]]; do
	case "$1" in
		-h|--help) help_info; exit 1 ;;
		-i|--input)       file_sam="${2:-}"; [[ -z "${file_sam}" || "${file_sam:0:1}" == "-" ]] && { echo "Error: no sam file as input"; exit 1; }; shift 2 ;;
		-r|--ref_genome)  file_ref="${2:-}"; [[ -z "${file_ref}" || "${file_ref:0:1}" == "-" ]] && { echo "Error: no fasta file as input"; exit 1; }; shift 2 ;;
		-o|--out)         file_path="${2:-}"; [[ -z "${file_path}" || "${file_path:0:1}" == "-" ]] && { echo "Error: no output path"; exit 1; }; shift 2 ;;  # fixed variable
		-p|--prefix)      prefix="${2:-}"; [[ -z "${prefix}" ]] && { echo "Error: no input for $1"; exit 1; }; shift 2 ;;
		-mq|--map_qual)   mq="${2:-}"; shift 2 ;;
		-t|--thread)      thread="${2:-}"; shift 2 ;;
		-e|--error_rate)  error_rate="${2:-}"; shift 2 ;;
		-s|--signi_level) signi_level="${2:-}"; shift 2 ;;
		-c|--cond_pro)    cond_pro="${2:-}"; shift 2 ;;
		-f|--fre_snv)     fre_snv="${2:-}"; shift 2 ;;
		-n1|--num_read_1) num_read_1="${2:-}"; shift 2 ;;
		-n2|--num_read_2) num_read_2="${2:-}"; shift 2 ;;
		-g|--gap)         gap="${2:-}"; shift 2 ;;
		-ss|--smallest_snv) smallest_snv="${2:-}"; shift 2 ;;
		-os|--only_snv)   only_snv="${2:-}"; shift 2 ;;
		-or|--ovlap_read) ovlap_read="${2:-}"; shift 2 ;;
		-wr|--weight_read) weight_read="${2:-}"; shift 2 ;;
		-sg|--sub_graph)  sub_graph="${2:-}"; shift 2 ;;
		-m|--mcl_inflation) mcl_inflation="${2:-}"; shift 2 ;;
		-oc|--ovlap_cluster) ovlap_cluster="${2:-}"; shift 2 ;;
		-wc|--weight_cluster) weight_cluster="${2:-}"; shift 2 ;;
		-d|--depth)       depth="${2:-}"; shift 2 ;;
		-l|--lar_cluster) lar_cluster="${2:-}"; shift 2 ;;
		-sp|--start_pos)  s_pos="${2:-}"; shift 2 ;;
		-ep|--end_pos)    e_pos="${2:-}"; shift 2 ;;
		-a|--abundance)   abundance="${2:-}"; shift 2 ;;
		-P|--polisher)    polisher="${2:-}"; shift 2 ;;  # NEW
		*) echo "Error: unknown parameter $1"; exit 1 ;;
	esac
done

[[ -z "$file_sam" ]] && { echo "Error: no sam file input"; exit 1; }
[[ -z "$file_ref" ]] && { echo "Error: no reference genome input"; exit 1; }

# normalize output path/prefix with or without trailing slash
if [[ "${file_path: -1}" == "/" ]]; then
  file_prefix="${file_path}${prefix}"
  file_path="${file_path%/}"
else
  file_prefix="${file_path}/${prefix}"
fi

##########  count nucleotide occurrence  ##########
echo "count nucleotide occurrence"
if [[ "$file_path" != "." ]]; then
	rm -rf "$file_path"
	mkdir -p "$file_path"
fi
rm -rf "$file_path/alignment"
mkdir -p "$file_path/alignment"

unique_sam="$file_path/alignment/${prefix}.sam"
samtools view -h -F 0x900 -q "$mq" "$file_sam" > "$unique_sam"
file_bam="$file_path/alignment/${prefix}.bam"
samtools view -b "$unique_sam" > "$file_bam"
rm "$unique_sam"
file_bam_sorted="$file_path/alignment/${prefix}_sorted.bam"
samtools sort "$file_bam" -o "$file_bam_sorted"
samtools index "$file_bam_sorted"

file_acgt="${file_prefix}_acgt.txt"
python ./src/count_frequency.py "$file_bam_sorted" "$file_acgt"

########## two binomial tests  ##########
echo "SNV detection"
file_snv="${file_prefix}_snv.txt"
python ./src/two_binomial.py "$error_rate" "$signi_level" "$file_acgt" "$file_snv" "$thread" "$s_pos" "$e_pos"

## judge number of detected SNV sites
size="$(wc -l < "$file_snv" || echo 0)"
if (( size != 0 )); then
	python ./src/out_haplotypes.py "${file_prefix}_clusters.pickle" "$file_bam_sorted" "$file_path" "$file_acgt" 1 "${file_prefix}_consensus.fasta" "$s_pos" "$e_pos"
	python ./src/extract_reads.py "$file_path" "$prefix" 1
	if [[ "$polisher" == "medaka" ]]; then
		python ./src/run_medaka.py "$file_path" "$prefix" 1
	else
		racon_polish_clusters "$file_path" "$prefix"
	fi
	exit 0
fi

## maximum conditional probability and construct reads graph
python ./src/mcp_read_graph.py "$file_bam_sorted" "$file_snv" "$cond_pro" "$smallest_snv" "$num_read_1" "$num_read_2" "$gap" \
	"$weight_read" "$ovlap_read" "$file_prefix" "$fre_snv" "$thread" "$only_snv" "$sub_graph"

## judge number of detected SNV sites (again)
size="$(wc -l < "$file_snv" || echo 0)"
if (( size != 0 )); then
	python ./src/out_haplotypes.py "${file_prefix}_clusters.pickle" "$file_bam_sorted" "$file_path" "$file_acgt" 1 "${file_prefix}_consensus.fasta" "$s_pos" "$e_pos"
	python ./src/extract_reads.py "$file_path" "$prefix" 1
	if [[ "$polisher" == "medaka" ]]; then
		python ./src/run_medaka.py "$file_path" "$prefix" 1
	else
		racon_polish_clusters "$file_path" "$prefix"
	fi
	exit 0
fi

(( only_snv != 0 )) && exit 0

## check the number of reads with overlaps
rg_file="${file_prefix}_reads_graph.txt"
size="$(wc -l < "$rg_file" || echo 0)"
if (( size == 0 )); then
	echo "No enough reads with overlaps"
	exit 0
fi

# MCL clustering
echo "MCL clustering"
python ./src/run_mcl.py "$file_prefix" "$thread" "$mcl_inflation" "$sub_graph"

## hierarchical clustering
echo "hierarchical clustering"
python ./src/hierarchical_cluster.py "${file_prefix}_matrix.pickle" "$lar_cluster" "$depth" \
	"$ovlap_cluster" "$weight_cluster" "$abundance" "$file_prefix"

## reconstruct haplotypes
rm -rf "$file_path/clusters"
mkdir -p "$file_path/clusters"

echo "haplotypes reconstruction"
python ./src/out_haplotypes.py "${file_prefix}_clusters.pickle" "$file_bam_sorted" "$file_path" "$file_acgt" x "${file_prefix}_consensus.fasta" "$s_pos" "$e_pos"

echo "haplotypes polishing ($polisher)"
python ./src/extract_reads.py "$file_path" "$prefix" x
if [[ "$polisher" == "medaka" ]]; then
	python ./src/run_medaka.py "$file_path" "$prefix" x
else
	racon_polish_clusters "$file_path" "$prefix"
fi

rm -f "${file_prefix}_matrix.pickle" "${file_prefix}_reads_cluster.txt" "${file_prefix}_clusters.pickle" || true
rm -rf "$file_path/medaka" || true

echo "complete reconstructing haplotypes"
exit 0
