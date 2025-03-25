#!/bin/bash
#PBS -N multiqc
#PBS -l wd
#PBS -l storage=gdata/za08+scratch/za08
#PBS -l walltime=1:00:00
#PBS -l mem=1GB
#PBS -l ncpus=1

module load singularity

singularity exec /g/data/za08/container_cache/multiqc_1.25.1--pyhdfd78af_0.sif multiqc --outdir results/ results/
