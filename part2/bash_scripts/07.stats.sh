#!/bin/bash
#PBS -q normal
#PBS -P er01
#PBS -l ncpus=1
#PBS -l mem=1GB
#PBS -l walltime=00:05:00
#PBS -l storage=scratch/er01
#PBS -l wd

set -euo pipefail

module load bcftools

mkdir -p results/vcf_stats
bcftools stats results/joint_geno/cohort.vcf.gz > results/vcf_stats/bcftools_stats.txt

echo DONE
