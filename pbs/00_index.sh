#!/usr/bin/bash
#PBS -N index
#PBS -l wd
#PBS -l storage=gdata/za08+scratch/za08
#PBS -l walltime=1:00:00
#PBS -l mem=1GB
#PBS -l ncpus=1

module load singularity 

mkdir "results"
singularity exec /g/data/za08/container_cache/salmon_1.10.3--h45fbf2d_4.sif \
	salmon index \
	--transcripts data/ggal/transcriptome.fa \
	--index results/salmon_index
