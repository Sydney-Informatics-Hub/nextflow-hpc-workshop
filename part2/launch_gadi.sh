#!/bin/bash -l
#PBS -N hello_nextflow
#PBS -P za08
#PBS -q copyq
#PBS -l ncpus=1
#PBS -l mem=8GB
#PBS -l walltime=01:00:00
#PBS -j oe
#PBS -o hello_nextflow.log
#PBS -m ae
#PBS -M your@email.com
#PBS -l storage=scratch/za08
#PBS -l wd

#load the nextflow module
module load nextflow/24.04.4

#Run the nextflow workflow with custom config
nextflow run main.nf --project za08 -c conf/gadi.config 
