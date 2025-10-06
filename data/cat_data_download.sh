#!/bin/bash -l
#PBS -N nffetchngs
#PBS -l select=1:ncpus=2:mem=4gb
#PBS -l walltime=24:00:00

## Data origin https://pubmed.ncbi.nlm.nih.gov/33780570/

ids=/scratch/er01/gs5517/nxf-hpc/nextflow-hpc-workshop/data/samplesheet_cat.csv
outdir=/scratch/er01/gs5517/nxf-hpc/nextflow-hpc-workshop/data 
config=/scratch/er01/gs5517/nxf-hpc/nextflow-hpc-workshop/data/fetchngs_gadi.config

module load nextflow 
module load singularity 

# configure singularity 
SINGULARITY_CACHEDIR=/scratch/er01/gs5517/singularity

# nextflow run 
nextflow run /scratch/er01/gs5517/nxf-hpc/fetchngs/main.nf \
  --input $ids \
  --outdir $outdir \
  -c $config -resume
