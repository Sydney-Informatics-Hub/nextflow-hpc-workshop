# Assigning process resources

!!! info "Learning objectives"

    - Identify how to size resource requests appropriately for each process
    - Apply resource-aware design principles to improve job efficiency
    - Optimise processes for time, noting the trade-offs with cost (SU usage)
    - Understand how to efficiently configure jobs to fit system queues

## Efficiency and scalability

Optimising for time and throughput - refer to concepts like more nodes, less walltime is more efficient than long walltime with less nodes (on Gadi, but ensure to generalise to any HPC/scheduler).

## Configuring processes

!!! example "Exercise"

    Create a file `conf/custom.config` and copy the following code chunk.
    Save the file.

    ```groovy title="custom.config"
    
    ```
    View trace, configure resources for each process in nextflow.config

Consider 

withName

withLabel

!!! example "Exercise"

    TODO 

Important: How to get "free" resources by correctly configuring to queue 
- Instructions to look at Gadi and Setonix queues/partitions

https://sydney-informatics-hub.github.io/training.gadi.intro/07-Optimisation/index.html

## Configuring to queues and partitions

Different queues/partitions are intended for different types of jobs

Provide a guided example


TODO excalidraw showing tradeoffs between 

!!! example "Exercise"

    TODO Review Gadi queues, Setonix partitions

    === "Gadi (PBS)"

        Review [Cascade Lake queues](https://opus.nci.org.au/spaces/Help/pages/236881198/Queue+Limits...#QueueLimits...-CascadeLakequeuelimits) (for simplicity, ignore bw, sl, gpus)
        
    === "Pawsey (Slurm)"


    TODO Then configure conf/<sched>.nf. Provide specific queues e.g. express
    vs. normal based on cpu/mem requirements, instead of configuring for all.

https://sydney-informatics-hub.github.io/usyd-gadi-onboarding-guide/notebooks/08_job_script.html#queue-selection-examples
    
Know how to set up system-specific config, how to ensure the resourcing aligns
well with the setup of the infrastructure.

HPC architecture differs across platforms and that the queue/partition
names and resources on that queue affect
the config files that needs to be created for that platform.

Tie back in that nextflow code can run on any platform, but when using HPC, the
config needs to be correct for that specific infrastructure.

## Configuring java heap sizes

!!! example "Exercises"

    TODO Update GENOTYPE and JOINT_GENOTYPE processes with -Xmx${tasks.memory}

Other things to consider - when writing custom R or Python scripts, writing
them efficiently. Utilising things like vectorisation, libraries such as numpy
etc., OpenMPI - link out.
