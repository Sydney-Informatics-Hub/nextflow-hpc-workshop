#!/bin/bash
#PBS -P ab01
#PBS -N concat
#PBS -q normalbw
#PBS -l ncpus=1
#PBS -l mem=1GB
#PBS -l walltime=00:01:00
#PBS -l storage=scratch/ab01
#PBS -l wd

module load samtools

set -euo pipefail

mkdir -p results/align_${SAMPLE_ID}
samtools cat -o results/align_${SAMPLE_ID}/${SAMPLE_ID}.bam results/align_${SAMPLE_ID}/${SAMPLE_ID}.split_*.bam
