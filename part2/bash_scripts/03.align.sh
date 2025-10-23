#!/bin/bash
#PBS -q normal
#PBS -P er01
#PBS -l ncpus=1
#PBS -l mem=1GB
#PBS -l walltime=00:05:00
#PBS -l storage=scratch/er01
#PBS -l wd

set -euo pipefail

module load bwa
module load samtools

FQDIR=$(dirname ${FQ1})
SAMPLE_SPLIT=$(basename ${FQ1} .R1.fq)
FQ2=${FQDIR}/${SAMPLE_SPLIT}.R2.fq

mkdir -p results/align_${SAMPLE_SPLIT}
bwa mem -t 1 -R "@RG\tID:${SAMPLE_ID}\tPL:ILLUMINA\tPU:${SAMPLE_ID}\tSM:${SAMPLE_ID}\tLB:${SAMPLE_ID}\tCN:SEQ_CENTRE" ref/BWAIndex/Homo_sapiens_assembly38.20-22.fasta ${FQ1} ${FQ2} | samtools view -b -o results/align_${SAMPLE_SPLIT}/${SAMPLE_SPLIT}.bam

echo DONE
