#!/bin/bash
#PBS -P ab01
#PBS -N split_fastq
#PBS -q normalbw
#PBS -l ncpus=1
#PBS -l mem=1GB
#PBS -l walltime=00:01:00
#PBS -l storage=scratch/ab01
#PBS -l wd

SAMPLE_ID="tiny"
READS_1="data/${SAMPLE_ID}.R1.fq"
READS_2="data/${SAMPLE_ID}.R2.fq"

mkdir -p results/fastq_split
split -l 12 -d --additional-suffix=".R1.fq" ${READS_1} results/fastq_split/${SAMPLE_ID}.split_
split -l 12 -d --additional-suffix=".R2.fq" ${READS_2} results/fastq_split/${SAMPLE_ID}.split_

for f in results/fastq_split/${SAMPLE_ID}.split_*.R1.fq
do
    qsub -v SAMPLE_ID="${SAMPLE_ID}",READS_1="${f}",READS_2="${f:0:(-4)}2.fq" align.gadi.sh
done