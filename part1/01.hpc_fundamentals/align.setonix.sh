#!/bin/bash
#SBATCH --account=pawsey1234
#SBATCH --job-name=align
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH --time=00:01:00

module load bwa/0.7.17--h7132678_9
module load samtools/1.15--h3843a85_0

set -euo pipefail

mkdir -p results/align_${SAMPLE_ID}
bwa mem -t 1 ref/BWAIndex/Homo_sapiens_assembly38.20-22.fasta ${READS_1} ${READS_2} | samtools view -b -o results/align_${SAMPLE_ID}/$(basename ${READS_1} .R1.fq).bam
