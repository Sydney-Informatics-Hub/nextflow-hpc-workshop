# Configuring custom Nextflow pipelines

!!! info "Learning objectives"

    - Identify best practices for benchmarking and configuring pipelines,
    step-by-step
    - Troubleshoot common HPC and scheduling errors by inspecting task logs
    and interpretting exit codes and error messages
    - Know how to find and specify system-dependent HPC values such as queues
    - Recall how Nextflow interacts with HPC components

We start with running things on a single sample. Then, optimise and configure the
pipeline with a single sample so that you can conserve SUs as you benchmark.

!!! example "Exercises"

    TODO: Run the pipeline with bare bones configuration. Specify the correct
    `profile` for your system.

    === "Gadi (PBS)"

        ```bash
        nextflow run main-nf -profile pbspro --pbspro_account vp91
        ```

    === "Setonix (Slurm)"

        ```bash
        nextflow run main-nf -profile slurm --slurm_account courses01
        ```

What are the errors?

!!! example "Exercises"
    
    TODO troubleshoot, fix, and rerun

!!! example "Exercises"

    TODO inspect a `.command.run` see how resources are allocated


Relate this back and compare why Nextflow is powerful for HPC vs. serially
running pbs/slurm scripts.

Rest of the content should be knowing how to set up system-specific config,
how to ensure the resourcing aligns well with the setup of the infrastructure.

Adapted from CW suggestion: Reiterate that HPC architecture differs across
platforms and that the queue/partition names and resources on that queue affect
the config files that needs to be created for that platform.

Tie back in that nextflow code can run on any platform, but when using HPC, the
config needs to be correct for that specific infrastructure.

Demo this with things like environmental variables, queue/partition names

!!! example "Exercises"

    TODO add the QoL additions from nf-core config 1.2, such as 

    === "Gadi (PBS)"

        ```groovy
        process {
            clusterOptions = "-P ${System.getenv('PROJECT')}"
            cache = 'lenient'
            stageInMode = 'symlink'
        }
        ```

    === "Setonix (Slurm)"

        ```bash
        process {
            clusterOptions = "--account=${System.getenv('PAWSEY_PROJECT')}"
            cache = 'lenient'
            stageInMode = 'symlink'
        }
        ```

    Build up to the config used in part 1.2 for
    [Gadi](https://github.com/Sydney-Informatics-Hub/nextflow-on-hpc-materials/blob/main/part1/config/gadi.config) and
    [Setonix](https://github.com/Sydney-Informatics-Hub/nextflow-on-hpc-materials/blob/main/part1/config/setonix.config)
