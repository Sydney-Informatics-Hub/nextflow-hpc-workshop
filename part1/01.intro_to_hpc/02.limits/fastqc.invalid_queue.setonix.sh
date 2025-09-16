#!/bin/bash -l
#SBATCH --account=$PAWSEY_PROJECT
#SBATCH --job-name=fastqc
#SBATCH --output=fastqc_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8GB
#SBATCH --partition=work
#SBATCH --time=36:00:00
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=END,FAIL

module load singularity

FASTQDIR="data/${SAMPLE_ID}"

mkdir -p "results/fastqc_${SAMPLE_ID}_logs"
singularity exec /g/data/za08/container_cache/fastqc_0.12.1--hdfd78af_0.sif fastqc \
	--outdir "results/fastqc_${SAMPLE_ID}_logs" \
	--format fastq ${FASTQDIR}/*.fastq.gz
