#!/bin/bash

# RNA-seq pipeline for bacterial data
# Features: Conda support, SLURM or local, YAML config, logging, genome auto-download

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

#check_conda_env() {
#  echo "[Setup] Checking Conda environment..."
#  if ! command -v conda &> /dev/null; then
#    echo "❌ Conda not found. Please install Miniconda or Anaconda first."
#    exit 1
#  fi

#  if ! conda info --envs | grep -q "$ENV_NAME"; then
#    echo "[Setup] Creating Conda environment '$ENV_NAME'..."
#    conda env create -y --name "$ENV_NAME" \
#      -f environment.yml
#  fi

#  echo "[Setup] Activating environment '$ENV_NAME'..."
#  source "$(conda info --base)/etc/profile.d/conda.sh"
#  conda activate "$ENV_NAME"
#  export PATH="$CONDA_PREFIX/bin:$PATH"
#} 

check_mamba_env() {
  echo "[Setup] Checking Mamba environment..."
  if ! command -v mamba &> /dev/null; then
    echo "❌ Mamba not found. Please install it !"
    exit 1
  fi

  if ! conda info --envs | grep -q "$ENV_NAME"; then
    echo "[Setup] Creating Mamba environment '$ENV_NAME'..."
    mamba env create -y --name "$ENV_NAME" -f environment.yml
  fi

  echo "[Setup] Activating environment '$ENV_NAME'..."
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "$ENV_NAME"
  export PATH="$CONDA_PREFIX/bin:$PATH"
}

submit_or_run() {
  local step="$1"
  shift
  local cmd="$*"

  if $USE_SLURM; then
    local extra_slurm_opts=""

#    if [[ "$step" == "download" ]]; then
#      extra_slurm_opts="--wait"
    if [[ "$step" == "rrna_filter" && "$cmd" == *"ribodetector"* ]]; then
      extra_slurm_opts="--threads-per-core=1"
    fi

    sbatch --job-name="$step" --cpus-per-task=$THREADS $extra_slurm_opts \
           --output="$LOG_DIR/${step}.log" --wait \
           --wrap="./slurm_wrapper.sh \"$step\" $cmd"
  else
    run_logged "$step" "$@"
  fi
}

source "$(dirname "$0")/pipeline_functions.sh"

#### RUN PIPELINE ####

check_mamba_env

submit_or_run "download" download_genome_data
submit_or_run "merge" merge_lanes
submit_or_run "qc" quality_control
submit_or_run "trim" trim_reads
submit_or_run "rrna_filter" rrna_filter
submit_or_run "align" align_reads
submit_or_run "count" count_features
submit_or_run "deseq2" run_deseq2

echo "✅ Pipeline complete. Logs are in $LOG_DIR"
