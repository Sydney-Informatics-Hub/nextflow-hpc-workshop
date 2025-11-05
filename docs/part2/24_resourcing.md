# Assigning process resources

!!! info "Learning objectives"

    - Identify how to size resource requests appropriately for each process
    - Apply resource-aware design principles to improve job efficiency
    - Optimise processes for time, noting the trade-offs with cost (SU usage)
    - Understand how to efficiently configure jobs to fit system queues

## Efficiency and scalability

Optimising for time and throughput - refer to concepts like more nodes, less walltime is more efficient than long walltime with less nodes (on Gadi, but ensure to generalise to any HPC/scheduler).

## Configuring to queues and partitions

!!! example "Exercise"

    TODO Review Gadi queues, Setonix partitions

    TODO Then configure conf/<sched>.nf. Provide specific queues e.g. express
    vs. normal based on cpu/mem requirements, instead of configuring for all.


Know how to set up system-specific config, how to ensure the resourcing aligns
well with the setup of the infrastructure.

HPC architecture differs across platforms and that the queue/partition
names and resources on that queue affect
the config files that needs to be created for that platform.

Tie back in that nextflow code can run on any platform, but when using HPC, the
config needs to be correct for that specific infrastructure.

## Configuring processes

!!! example "Exercise"

    TODO View trace, configure resources for each process in nextflow.config

Consider 

withName

withLabel

!!! example "Exercise"

    TODO 

Important: How to get "free" resources by correctly configuring to queue 
- Instructions to look at Gadi and Setonix queues/partitions


https://sydney-informatics-hub.github.io/training.gadi.intro/07-Optimisation/index.html

## Configuring java heap sizes

!!! example "Exercises"

    TODO Update GENOTYPE and JOINT_GENOTYPE processes with -Xmx${tasks.memory}

Other things to consider - when writing custom R or Python scripts, writing
them efficiently. Utilising things like vectorisation, libraries such as numpy
etc., OpenMPI - link out.
