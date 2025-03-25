#!/bin/bash -l
#SBATCH --account=$PAWSEY_PROJECT
#SBATCH --job-name=hello_nextflow
#SBATCH --output=hello_nextflow_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH --partition=work
#SBATCH --time=01:00:00
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=END,FAIL

#load the nextflow module
module load nextflow/24.04.3

#Run the nextflow workflow with custom config
nextflow run main.nf -c conf/setonix.config 
