#!/bin/bash -l
#SBATCH --account=$PAWSEY_PROJECT
#SBATCH --job-name=index
#SBATCH --output=index_%j.log
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

singularity pull https://depot.galaxyproject.org/singularity/salmon:1.10.3--h6dccd9a_2

mkdir "results"
singularity exec salmon:1.10.3--h6dccd9a_2 salmon index \
	--transcripts data/ggal/transcriptome.fa \
	--index results/salmon_index
