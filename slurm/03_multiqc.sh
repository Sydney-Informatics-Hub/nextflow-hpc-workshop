#!/bin/bash -l
#SBATCH --account=$PAWSEY_PROJECT
#SBATCH --job-name=multiqc
#SBATCH --output=multiqc_%j.log
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

singularity pull https://depot.galaxyproject.org/singularity/multiqc:1.25.1--pyhdfd78af_0

singularity exec multiqc:1.25.1--pyhdfd78af_0 --outdir results/ results/
