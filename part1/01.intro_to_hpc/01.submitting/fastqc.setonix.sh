#!/bin/bash -l
#SBATCH --account=$PAWSEY_PROJECT
#SBATCH --job-name=fastqc
#SBATCH --output=fastqc_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH --partition=work
#SBATCH --time=00:10:00
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=END,FAIL

module load pawseyenv/2023.08
module load singularity/3.11.4-nompi

singularity pull docker://biocontainers/fastqc:0.12.1--hdfd78af_0

SAMPLE_ID="gut"
READS_1="data/ggal/${SAMPLE_ID}_1.fq"
READS_2="data/ggal/${SAMPLE_ID}_2.fq"

mkdir -p "results/fastqc_${SAMPLE_ID}_logs"
singularity exec fastqc:0.12.1--hdfd78af_0 fastqc \
	--outdir "results/fastqc_${SAMPLE_ID}_logs" \
	--format fastq ${READS_1} ${READS_2}
