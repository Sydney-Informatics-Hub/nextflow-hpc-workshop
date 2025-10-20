#!/bin/bash
#PBS -q normal
#PBS -P er01
#PBS -l ncpus=2
#PBS -l mem=1GB
#PBS -l walltime=00:02:00
#PBS -l storage=scratch/er01
#PBS -l wd

set -euo pipefail

module load fastqc

mkdir -p results/fastqc_${SAMPLE_ID}_logs
fastqc -t 2 --outdir results/fastqc_${SAMPLE_ID}_logs --format fastq data/${SAMPLE_ID}.R1.fq data/${SAMPLE_ID}.R2.fq

echo DONE
