#!/bin/bash
#PBS -q normal
#PBS -P er01
#PBS -l ncpus=4
#PBS -l mem=4GB
#PBS -l walltime=00:10:00
#PBS -l storage=scratch/er01
#PBS -l wd

set -euo pipefail

module load gatk

mkdir -p results/geno_${SAMPLE_ID}
gatk --java-options "-Xmx4g" HaplotypeCaller -R ref/Homo_sapiens_assembly38.20-22.fasta -I results/merge_${SAMPLE_ID}/${SAMPLE_ID}.bam -O results/geno_${SAMPLE_ID}/${SAMPLE_ID}.g.vcf.gz -ERC GVCF

echo DONE
