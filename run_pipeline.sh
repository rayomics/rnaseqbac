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
RRNA_LSU_FASTA="$GENOME_DIR/rrna_lsu.fa"
RRNA_SSU_FASTA="$GENOME_DIR/rrna_ssu.fa"
GENOME_INDEX="$GENOME_DIR/genome"
GENOME_URL=$(read_config genome.fasta_url)
GTF_URL=$(read_config genome.gtf_url)
RRNA_LSU_URL=$(read_config genome.silva_db_lsu)
RRNA_SSU_URL=$(read_config genome.silva_db_ssu)

mkdir -p "$RESULTS_DIR" "$LOG_DIR" "$MERGED_DIR" "$QC_DIR" "$TRIM_DIR" "$ALIGN_DIR" "$COUNT_DIR" "$DESEQ_DIR" "$GENOME_DIR"

#### CONDA ENV SETUP ####

ENV_NAME="rnaseq_env"

check_conda_env() {
  echo "[Setup] Checking Conda environment..."
  if ! command -v conda &> /dev/null; then
    echo "❌ Conda not found. Please install Miniconda or Anaconda first."
    exit 1
  fi

  if ! conda info --envs | grep -q "$ENV_NAME"; then
    echo "[Setup] Creating Conda environment '$ENV_NAME'..."
    conda env create -y --name "$ENV_NAME" \
      -f environment.yml
  fi

  echo "[Setup] Activating environment '$ENV_NAME'..."
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "$ENV_NAME"
}

#### LOGGING ####

run_logged() {
  local step="$1"
  shift
  echo "[Running] $step"
  {
    echo ">>> Step: $step"
    date
    "$@"
    echo "<<< Completed: $step"
  } > "$LOG_DIR/${step}.log" 2>&1
}

submit_or_run() {
  local step="$1"
  shift
  local cmd="$*"
  if $USE_SLURM; then
    sbatch --job-name="$step" --cpus-per-task=$THREADS \
           --output="$LOG_DIR/${step}.log" --wrap="$cmd"
  else
    run_logged "$step" "$cmd"
  fi
}

#### GENOME DOWNLOAD & INDEXING ####

download_genome_data() {
  if [ ! -f "$GENOME_FASTA" ]; then
    echo "[Download] Genome FASTA..."
    curl -L "$GENOME_URL" | gunzip -c > "$GENOME_FASTA"
  fi

  if [ ! -f "$GTF_FILE" ]; then
    echo "[Download] GTF annotation..."
    curl -L "$GTF_URL" | gunzip -c > "$GTF_FILE"
  fi

  if [ ! -f "$RRNA_LSU_FASTA" ]; then
    echo "[Download] rRNA LSU data..."
    curl -L "$RRNA_LSU_URL" | gunzip -c > "$RRNA_LSU_FASTA"
  fi

  if [ ! -f "$RRNA_SSU_FASTA" ]; then
    echo "[Download] rRNA SSU data..."
    curl -L "$RRNA_SSU_URL" | gunzip -c > "$RRNA_SSU_FASTA"
  fi

  echo "[Indexing] Building $ALIGNER index..."
  case "$ALIGNER" in
    hisat2)
      if [ ! -f "${GENOME_INDEX}.1.ht2" ]; then
        hisat2-build -p "$THREADS" "$GENOME_FASTA" "$GENOME_INDEX"
      fi
      ;;
    bowtie2)
      if [ ! -f "${GENOME_INDEX}.1.bt2" ]; then
        bowtie2-build --threads "$THREADS" "$GENOME_FASTA" "$GENOME_INDEX"
      fi
      ;;
    *)
      echo "ERROR: Unknown aligner '$ALIGNER'" >&2
      exit 1
      ;;
  esac 
}


#### PIPELINE STEPS ####

merge_lanes() {
  for R1 in "$RAW_DIR"/*_R1_001.fastq.gz; do
    sample=$(basename "$R1" | cut -d_ -f1)
    R1_out="$MERGED_DIR/${sample}_R1.fastq.gz"
    R2_out="$MERGED_DIR/${sample}_R2.fastq.gz"
    cat "$RAW_DIR/${sample}"*_R1_001.fastq.gz > "$R1_out"
    cat "$RAW_DIR/${sample}"*_R2_001.fastq.gz > "$R2_out"
  done
}

quality_control() {
  fastqc -t "$THREADS" -o "$QC_DIR" "$MERGED_DIR"/*.fastq.gz
  multiqc -o "$QC_DIR" "$QC_DIR"
}

trim_reads() {
  for R1 in "$MERGED_DIR"/*_R1.fastq.gz; do
    sample=$(basename "$R1" _R1.fastq.gz)
    R2="$MERGED_DIR/${sample}_R2.fastq.gz"
    fastp -i "$R1" -I "$R2" \
          -o "$TRIM_DIR/${sample}_R1.trimmed.fastq.gz" \
          -O "$TRIM_DIR/${sample}_R2.trimmed.fastq.gz" \
          -w "$THREADS"
  done
}

rrna_filter() {
  case "$RRNA_FILTER" in
    sortmerna)
      for R1 in "$TRIM_DIR"/*_R1.trimmed.fastq.gz; do
        sample=$(basename "$R1" _R1.trimmed.fastq.gz)
        R2="$TRIM_DIR/${sample}_R2.trimmed.fastq.gz"
        sortmerna --reads "$R1" --reads "$R2" --ref "$RRNA_FULL_FASTA" \
        --out2 --aligned "$ALIGN_DIR/${sample}.trimmed.aligned" --other "$ALIGN_DIR/${sample}.trimmed.nonrrna" \
        --fastx --workdir "$ALIGN_DIR" --idx-dir "$ALIGN_DIR/idx" --threads $THREADS
        rm -r "$ALIGN_DIR/kvdb"
      done
      ;;

    ribodetector)
      for R1 in "$TRIM_DIR"/*_R1.trimmed.fastq.gz; do
        sample=$(basename "$R1" _R1.trimmed.fastq.gz)
        R2="$TRIM_DIR/${sample}_R2.trimmed.fastq.gz"
        ribodetector_cpu -t "$THREADS" -l "$(zcat $R1 | head | awk 'NR==2 {print length($0)}')" \
        -i "$R1" "$R2" -e none --chunk_size 256 \
        -o "$ALIGN_DIR/${sample}.trimmed.nonrrna.1.fq.gz" "$ALIGN_DIR/${sample}.trimmed.nonrrna.2.fq.gz"            
      done
      
      ;;
    *)
      echo "ERROR: Unknown rRNA filtering tool '$RRNA_FILTER'" >&2
      exit 1
      ;;
  esac 
}

align_reads() {
  for R1 in "$ALIGN_DIR"/*.trimmed.nonrrna.1.fq.gz; do
    sample=$(basename "$R1" .trimmed.nonrrna.1.fq.gz)
    R2="$ALIGN_DIR/${sample}.trimmed.nonrrna.2.fq.gz"
    SAM="$ALIGN_DIR/${sample}.sam"

    case "$ALIGNER" in
      hisat2)
        hisat2 -p "$THREADS" -x "$GENOME_INDEX" -1 "$R1" -2 "$R2" -S "$SAM"
        ;;
      bowtie2)
        bowtie2 -p "$THREADS" -x "$GENOME_INDEX" -1 "$R1" -2 "$R2" -S "$SAM"
        ;;
      *)
        echo "ERROR: Unknown aligner '$ALIGNER'" >&2
        exit 1
        ;;
    esac
  done
}

count_features() {
  BAM_LIST=()
  for SAM in "$ALIGN_DIR"/*.sam; do
    BAM="${SAM%.sam}.bam"
    samtools view -@ "$THREADS" -Sb "$SAM" | samtools sort -@ "$THREADS" -o "$BAM"
    rm "$SAM"
    BAM_LIST+=("$BAM")
  done
  featureCounts -p -T "$THREADS" -a "$GTF_FILE" -o "$COUNT_DIR/counts.txt" "${BAM_LIST[@]}"
}

run_deseq2() {
  Rscript scripts/run_deseq2.R "$COUNT_DIR/counts.txt" "$DESEQ_DIR"
}

#### RUN PIPELINE ####

check_conda_env
download_genome_data

submit_or_run "merge" merge_lanes
submit_or_run "qc" quality_control
submit_or_run "trim" trim_reads
submit_or_run "rrna_filter" rrna_filter
submit_or_run "align" align_reads
submit_or_run "count" count_features
submit_or_run "deseq2" run_deseq2

echo "✅ Pipeline complete. Logs are in $LOG_DIR"
