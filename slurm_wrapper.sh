#!/bin/bash

set -euo pipefail

#### READ CONFIG (YAML) ####

read_config() {
  local key="$1"
  python3 -c "
import yaml
with open('config.yaml') as f:
    cfg = yaml.safe_load(f)
keys = '$key'.split('.')
val = cfg
for k in keys:
    val = val.get(k, None)
    if val is None:
        break
print(val if val else '')
" 2>/dev/null
}

#### CONFIG VALUES ####

RESOURCE_PROFILE=$(read_config resource_profile || echo "low")
THREADS=$(read_config threads || echo "4")
USE_SLURM=false
[[ "$RESOURCE_PROFILE" == "high" ]] && USE_SLURM=true

# Tool settings
ALIGNER=$(read_config genome.aligner || echo "hisat2")
RRNA_FILTER=$(read_config genome.rrna_filter || echo "ribodetector")

# Paths
RAW_DIR=$(read_config paths.raw_dir || echo "raw_data")
RESULTS_DIR=$(read_config paths.results_dir || echo "results")
GENOME_DIR=$(read_config paths.genome_dir || echo "genome")
MERGED_DIR="$RESULTS_DIR/merged"
QC_DIR="$RESULTS_DIR/qc"
TRIM_DIR="$RESULTS_DIR/trimmed"
ALIGN_DIR="$RESULTS_DIR/aligned"
COUNT_DIR="$RESULTS_DIR/counts"
DESEQ_DIR="$RESULTS_DIR/deseq2"
LOG_DIR="logs"

# Genome files
GENOME_FASTA="$GENOME_DIR/genome.fa"
GTF_FILE="$GENOME_DIR/annotation.gtf"
RRNA_LSU_FASTA=$(read_config genome.silva_db_lsu)
RRNA_SSU_FASTA=$(read_config genome.silva_db_ssu)
GENOME_INDEX="$GENOME_DIR/genome"
GENOME_URL=$(read_config genome.fasta_url)
GTF_URL=$(read_config genome.gtf_url)
RRNA_LSU_URL=$(read_config genome.silva_db_lsu)
RRNA_SSU_URL=$(read_config genome.silva_db_ssu)

mkdir -p "$RESULTS_DIR" "$LOG_DIR" "$MERGED_DIR" "$QC_DIR" "$TRIM_DIR" "$ALIGN_DIR" "$COUNT_DIR" "$DESEQ_DIR" "$GENOME_DIR"

#### CONDA ENV SETUP ####

ENV_NAME="rnaseq_env"

source "$(dirname "$0")/pipeline_functions.sh"

step="$1"
shift
"$@"