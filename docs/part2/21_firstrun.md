# Configuring custom Nextflow pipelines

!!! info "Learning objectives"

    - Identify best practices for benchmarking and configuring pipelines,
    step-by-step
    - Troubleshoot common HPC and scheduling errors by inspecting task logs
    and interpretting exit codes and error messages
    - Recall how Nextflow interacts with HPC components
    - Recall how executors queues, and work directories control task execution
    on HPC

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

!!! example "Exercises"

    TODO add the QoL additions from nf-core config, such as 

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

    TODO modularise pipeline? or from the very beginning?
