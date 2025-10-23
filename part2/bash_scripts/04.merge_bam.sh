#!/bin/bash
#PBS -q normal
#PBS -P er01
#PBS -l ncpus=1
#PBS -l mem=1GB
#PBS -l walltime=00:05:00
#PBS -l storage=scratch/er01
#PBS -l wd

set -euo pipefail

module load samtools

mkdir -p results/merge_${SAMPLE_ID}
samtools cat results/align_${SAMPLE_ID}.split_*/*.bam | samtools sort -O bam -o results/merge_${SAMPLE_ID}/${SAMPLE_ID}.bam
samtools index results/merge_${SAMPLE_ID}/${SAMPLE_ID}.bam

echo DONE
