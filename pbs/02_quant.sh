#!/usr/bin/bash
#PBS -N quant
#PBS -l wd
#PBS -l storage=gdata/za08+scratch/za08
#PBS -l walltime=1:00:00
#PBS -l mem=1GB
#PBS -l ncpus=1

module load singularity

SAMPLE_ID="gut"
READS_1="data/ggal/${SAMPLE_ID}_1.fq"
READS_2="data/ggal/${SAMPLE_ID}_2.fq"

singularity exec /g/data/za08/container_cache/salmon_1.10.3--h45fbf2d_4.sif \
	salmon quant \
	--libType=U \
	-i results/salmon_index \
	-1 ${READS_1} \
	-2 ${READS_2} \
	-o results/${SAMPLE_ID}
