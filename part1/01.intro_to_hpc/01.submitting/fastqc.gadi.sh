#!/usr/bin/bash
#PBS -N fastqc
#PBS -l wd
#PBS -l storage=gdata/za08+scratch/za08
#PBS -l walltime=1:00:00
#PBS -l mem=1GB
#PBS -l ncpus=1

module load singularity

SAMPLE_ID="gut"
READS_1="data/ggal/${SAMPLE_ID}_1.fq"
READS_2="data/ggal/${SAMPLE_ID}_2.fq"

mkdir -p "results/fastqc_${SAMPLE_ID}_logs"
singularity exec /g/data/za08/container_cache/fastqc_0.12.1--hdfd78af_0.sif fastqc \
	--outdir "results/fastqc_${SAMPLE_ID}_logs" \
	--format fastq ${READS_1} ${READS_2}
