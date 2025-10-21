#!/bin/bash
#SBATCH --account=pawsey1234
#SBATCH --job-name=concat
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH --time=00:01:00

module load samtools/1.15--h3843a85_0

set -euo pipefail

mkdir -p results/align_${SAMPLE_ID}
samtools cat -o results/align_${SAMPLE_ID}/${SAMPLE_ID}.bam results/align_${SAMPLE_ID}/${SAMPLE_ID}.split_*.bam
