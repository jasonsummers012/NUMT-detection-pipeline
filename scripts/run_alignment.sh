#!/bin/bash
set -euo pipefail

# Config
PROJECT_DIR="."
REF_DIR="./data"
REF_NAME="combined_reference.fa"
READS1="$REF_DIR/SRR13269374_1.fastq.gz"
READS2="$REF_DIR/SRR13269374_2.fastq.gz"
THREADS=4

# Output
OUT_DIR="$PROJECT_DIR/results"
PREFIX="sample_alignment"
BWA_INDEX_PREFIX="$REF_DIR/$REF_NAME"

# Setup
mkdir -p "$OUT_DIR"

# Indexing Step
if [ ! -f "${BWA_INDEX_PREFIX}.bwt" ]; then
    echo "Indexing reference genome..."
    bwa index "$BWA_INDEX_PREFIX"
else
    echo "Reference already indexed."
fi

# Alignment Step
echo "Running BWA MEM alignment..."
bwa mem -t "$THREADS" "$BWA_INDEX_PREFIX" "$READS1" "$READS2" > "$OUT_DIR/${PREFIX}.sam"

# SAM to BAM
echo "Converting SAM to BAM..."
samtools view -bS "$OUT_DIR/${PREFIX}.sam" > "$OUT_DIR/${PREFIX}.bam"

# Sorting BAM
echo "Sorting BAM file..."
samtools sort "$OUT_DIR/${PREFIX}.bam" -o "$OUT_DIR/${PREFIX}_sorted.bam"

# Index BAM
echo "Indexing sorted BAM..."
samtools index "$OUT_DIR/${PREFIX}_sorted.bam"

echo "Pipeline completed successfully."