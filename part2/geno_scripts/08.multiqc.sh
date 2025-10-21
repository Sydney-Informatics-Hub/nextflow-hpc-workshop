#!/bin/bash
#PBS -q normal
#PBS -P er01
#PBS -l ncpus=1
#PBS -l mem=1GB
#PBS -l walltime=00:05:00
#PBS -l storage=scratch/er01
#PBS -l wd

set -euo pipefail

module load multiqc

mkdir -p results/multiqc
for f in ${PWD}/results/fastqc_*/*; do ln -s $f results/multiqc; done
ln -s ${PWD}/results/vcf_stats/bcftools_stats.txt results/multiqc
cd results/multiqc

multiqc .

echo DONE
