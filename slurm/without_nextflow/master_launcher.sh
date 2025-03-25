#!/bin/bash -l
#SBATCH --account=$PAWSEY_PROJECT
#SBATCH --job-name=sequential_pipeline
#SBATCH --output=pipeline_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4GB
#SBATCH --partition=debug
#SBATCH --time=1:00:00
#SBATCH --mail-user=you@email.com
#SBATCH --mail-type=END,FAIL

# Path to your scripts
SCRIPT1="./00_index.sh"
SCRIPT2="./01_fastqc.sh"
SCRIPT3="./02_quant.sh"
SCRIPT4="./03_multiqc.sh"

# Make sure scripts are executable
chmod +x $SCRIPT1 $SCRIPT2 $SCRIPT3 $SCRIPT4

# Submit first job
echo "Submitting job 1..."
JOB1=$(sbatch --parsable $SCRIPT1)
echo "Job 1 submitted with ID: $JOB1"

# Submit second job with dependency on first job
echo "Submitting job 2 (depends on job 1)..."
JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 $SCRIPT2)
echo "Job 2 submitted with ID: $JOB2"

# Submit third job with dependency on second job
echo "Submitting job 3 (depends on job 2)..."
JOB3=$(sbatch --parsable --dependency=afterok:$JOB2 $SCRIPT3)
echo "Job 3 submitted with ID: $JOB3"

# Submit fourth job with dependency on third job
echo "Submitting job 4 (depends on job 3)..."
JOB4=$(sbatch --parsable --dependency=afterok:$JOB3 $SCRIPT4)
echo "Job 4 submitted with ID: $JOB4"

echo "All jobs submitted successfully!"
echo "Pipeline: $JOB1 -> $JOB2 -> $JOB3 -> $JOB4"
