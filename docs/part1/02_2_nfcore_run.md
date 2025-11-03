# 2.2 Running nf-core pipelines

!!! info "Learning objectives"

    - Build a basic run script for a nf-core pipeline
    - Run an nf-core pipeline
    - Understand the need for HPC-specific configurations
    
## 2.2.1 Prepraring the environment

!!! example "Exercise: review the nf-core/sarek environment"

    TODO:

    - Move into the nf-core directory of the materials repo
    - Review the input files
        - BAM/BAI files
        - Reference FASTA
        - samplesheet.csv

## 2.2.2 Write a simple run script

!!! example "Exercise: Create a run script for nf-core/sarek"

    TODO:

    - Create a new blank file called `run.sh`
    - Build up the basic run command
        - `--input samplesheet.csv`: Defines the input samples and the paths to their BAMs
        - `--step markduplicates`: Defines the starting point of the pipeline
        - `--skip_tools baserecalibrator,mosdepth,samtools`: Tells the pipeline to skip certain tools - means we only run `markduplicates`
        - `--fasta/--fasta_fai`: Define the location of the reference FASTA and its index
        - `--outdir results`: Define the output directory
        - `-resume`: Let's us resume from previously finished runs
    - Add a few extra options to simplify the run for our example
        - `--no_intervals true`: Prevents a few initial processes from running that aren't needed for our example
        - `--igenomes_ignore true`: Prevents downloading external resources that we don't need
    - Make the script executable

## 2.2.3 Running the pipeline

!!! example "Exercise: Run the workflow"

    TODO: Have participants run the script in the terminal, see what happens. It will quickly fail due to missing software, as we haven't defined the use of singularity or any executors. Step through the error messages and explain what happened and what we need to define in the next section.

    - Singularity: Lets us use containers for running the tools - don't need the software pre-installed on the system
    - Executor: Lets us specify that we want to use an HPC scheduler like PBSPro or SLURM. Without this, Nextflow will try to run on the head node (BAD!). The executor definition will let us handle resource requirements by letting Nextflow talk to the scheduler and submit each task as a job to the cluster.