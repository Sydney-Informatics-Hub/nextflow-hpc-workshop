#!/bin/bash
#SBATCH --account=pawsey1234
#SBATCH --job-name=split_fastq
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH --time=00:01:00

SAMPLE_ID="tiny"
READS_1="data/${SAMPLE_ID}.R1.fq"
READS_2="data/${SAMPLE_ID}.R2.fq"

mkdir -p results/fastq_split
split -l 12 -d --additional-suffix=".R1.fq" ${READS_1} results/fastq_split/${SAMPLE_ID}.split_
split -l 12 -d --additional-suffix=".R2.fq" ${READS_2} results/fastq_split/${SAMPLE_ID}.split_

for f in results/fastq_split/${SAMPLE_ID}.split_*.R1.fq
do
    sbatch --export SAMPLE_ID="${SAMPLE_ID}",READS_1="${f}",READS_2="${f:0:(-4)}2.fq" align.setonix.sh
done