#!/bin/bash -l
#PBS -N nf_genotype
#PBS -P za08
#PBS -q copyq
#PBS -l ncpus=1
#PBS -l mem=8GB
#PBS -l walltime=00:20:00
#PBS -j oe
#PBS -o nf_genotype.log
#PBS -m ae
#PBS -M your@email.com
#PBS -l storage=scratch/za08
#PBS -l wd

#load the nextflow module
module load nextflow/24.04.5

# Run the nextflow workflow with custom config
nextflow run main.nf -c gadi.config -resume
