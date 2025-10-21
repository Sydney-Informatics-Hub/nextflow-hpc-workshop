#!/bin/bash
#PBS -P ab01
#PBS -N align
#PBS -q normalbw
#PBS -l ncpus=1
#PBS -l mem=1GB
#PBS -l walltime=00:01:00
#PBS -l storage=scratch/ab01
#PBS -l wd

module load bwa
module load samtools

set -euo pipefail

mkdir -p results/align_${SAMPLE_ID}
bwa mem -t 1 ref/BWAIndex/Homo_sapiens_assembly38.20-22.fasta ${READS_1} ${READS_2} | samtools view -b -o results/align_${SAMPLE_ID}/$(basename ${READS_1} .R1.fq).bam
