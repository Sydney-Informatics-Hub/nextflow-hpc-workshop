#!/bin/bash
# Script to build a STAR genome index

# Set paths to genome FASTA and annotation GTF files
GENOME_FASTA="genome.fa"
ANNOTATION_GTF="annotation.gtf"
INDEX_DIR="star_index"
OVERHANG=99

# Number of threads to use
THREADS=8

# Create index directory if it doesn't exist
mkdir -p "$INDEX_DIR"

# Run STAR to generate the genome index
STAR --runThreadN "$THREADS" \
    --runMode genomeGenerate \
    --genomeDir "$INDEX_DIR" \
    --genomeFastaFiles "$GENOME_FASTA" \
    --sjdbGTFfile "$ANNOTATION_GTF" \
    --sjdbOverhang "$OVERHANG"

echo "STAR genome index created in $INDEX_DIR"