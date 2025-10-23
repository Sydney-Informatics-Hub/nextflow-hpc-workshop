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

mkdir -p results/joint_geno
gatk --java-options "-Xmx4g" CombineGVCFs -R ref/Homo_sapiens_assembly38.20-22.fasta $(for f in results/geno_*/*.g.vcf.gz; do printf -- "--variant $f "; done) -O results/joint_geno/cohort.g.vcf.gz
gatk --java-options "-Xmx4g" GenotypeGVCFs -R ref/Homo_sapiens_assembly38.20-22.fasta -V results/joint_geno/cohort.g.vcf.gz -O results/joint_geno/cohort.vcf.gz

echo DONE
