#!/bin/bash

set -euo pipefail





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

  #if [ ! -f "$RRNA_LSU_FASTA" ]; then
  #  echo "[Download] rRNA LSU data..."
  #  curl -L "$RRNA_LSU_URL" | gunzip -c > "$RRNA_LSU_FASTA"
  #fi

  #if [ ! -f "$RRNA_SSU_FASTA" ]; then
  #  echo "[Download] rRNA SSU data..."
  #  curl -L "$RRNA_SSU_URL" | gunzip -c > "$RRNA_SSU_FASTA"
  #fi

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
    sample=$(basename "$R1" | cut -d_ -f1-2)
    R1_out="$MERGED_DIR/${sample}_R1.fastq.gz"
    R2_out="$MERGED_DIR/${sample}_R2.fastq.gz"
    cat "$RAW_DIR/${sample}"*_R1_001.fastq.gz > "$R1_out" &
    cat "$RAW_DIR/${sample}"*_R2_001.fastq.gz > "$R2_out" &
  done
  wait
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
      --detect_adapter_for_pe \
      --disable_adapter_trimming \
      --disable_trim_poly_g \
      --disable_quality_filtering \
      --disable_length_filtering \
      --qualified_quality_phred 20 \
      --unqualified_percent_limit 0 \
      --length_required 75 \
      --thread "$THREADS" &
  done
  wait
}

rrna_filter() {
  case "$RRNA_FILTER" in
    sortmerna)
      for R1 in "$TRIM_DIR"/*_R1.trimmed.fastq.gz; do
        sample=$(basename "$R1" _R1.trimmed.fastq.gz)
        R2="$TRIM_DIR/${sample}_R2.trimmed.fastq.gz"
        sortmerna --reads "$R1" --reads "$R2" --ref "$RRNA_LSU_FASTA" --ref "$RRNA_SSU_FASTA" \
        --out2 --other "$ALIGN_DIR/${sample}.trimmed.nonrrna" \
        --fastx --workdir "$ALIGN_DIR/$sample" --idx-dir "$(dirname "$RRNA_LSU_FASTA")/idx" --threads 16 --paired_in --num_alignments 1 -v
#       (( count++ ))
#        echo "$count"
#        if (( count % MAX_JOBS == 0 )); then
#          wait
#        fi
        rm -rf "$ALIGN_DIR/$sample/kvdb"
        rm -rf "$ALIGN_DIR/$sample/readb"
        rm -rf "$ALIGN_DIR/*.sam"
        rm -rf "$ALIGN_DIR/$sample/out"
      done
      #wait
      ;;

    ribodetector)
      for R1 in "$TRIM_DIR"/*_R1.trimmed.fastq.gz; do
        sample=$(basename "$R1" _R1.trimmed.fastq.gz)
        R2="$TRIM_DIR/${sample}_R2.trimmed.fastq.gz"
        chunk_opt=""
        if ! $USE_SLURM ; then
          chunk_opt="--chunk_size 5000"
        fi
        ribodetector_cpu -t "$THREADS" -l "$(zcat $R1 | head | awk 'NR==2 {print length($0)}')" \
        -i "$R1" "$R2" -e norrna $chunk_opt \
        -o "$ALIGN_DIR/${sample}.trimmed.nonrrna_fwd.fq.gz" "$ALIGN_DIR/${sample}.trimmed.nonrrna_rev.fq.gz"
      done

      ;;
    *)
      echo "ERROR: Unknown rRNA filtering tool '$RRNA_FILTER'" >&2
      exit 1
      ;;
  esac
}

align_reads() {
  for R1 in "$ALIGN_DIR"/*.trimmed.nonrrna_fwd.fq.gz; do
    sample=$(basename "$R1" .trimmed.nonrrna_fwd.fq.gz)
    R2="$ALIGN_DIR/${sample}.trimmed.nonrrna_rev.fq.gz"
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
  for SAM in "$ALIGN_DIR"/*.sam; do
    BAM="${SAM%.sam}.bam"
    samtools view -@ "$THREADS" -Sb "$SAM" | samtools sort -@ "$THREADS" -o "$BAM"
    rm "$SAM"
  done
}

count_features() {
  BAM_LIST=()
  for SAM in "$ALIGN_DIR"/*.bam; do
    BAM="${SAM%.sam}"
    BAM_LIST+=("$BAM")
  done
  featureCounts -p -B -C -t CDS -s 2 -g "$GENE_IDENTIFIER" -T "$THREADS" -a "$GTF_FILE" -o "$COUNT_DIR/counts.txt" "${BAM_LIST[@]}"
}

run_deseq2() {
  Rscript scripts/run_deseq2.R "$COUNT_DIR/counts.txt" "sample_sheet.tsv" "$DESEQ_DIR"
}
