#!/usr/bin/bash
#PBS -N fastqc
#PBS -q normalbw
#PBS -l wd
#PBS -l storage=gdata/za08+scratch/za08
#PBS -l walltime=00:00:30
#PBS -l mem=1GB
#PBS -l ncpus=1

module load singularity

FASTQDIR="data/${SAMPLE_ID}"

mkdir -p "results/fastqc_${SAMPLE_ID}_logs"
singularity exec /g/data/za08/container_cache/fastqc_0.12.1--hdfd78af_0.sif fastqc \
	--outdir "results/fastqc_${SAMPLE_ID}_logs" \
	--format fastq ${FASTQDIR}/*.fastq.gz
