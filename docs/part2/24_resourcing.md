# Assigning process resources

!!! info "Learning objectives"

    - Identify how to size resource requests appropriately for each process
    - Apply resource-aware design principles to improve job efficiency
    - Optimise processes for time, noting the trade-offs with cost (SU usage)
    - Understand how to efficiently configure jobs to fit system queues
    - Apply infrastructure requirements of memory:CPU ratios to match the node
    architecture of the target HPC queue/partition

## Efficiency and scalability

Optimising for time and throughput - refer to concepts like more nodes, less walltime is more efficient than long walltime with less nodes (on Gadi, but ensure to generalise to any HPC/scheduler).

Tetris explanation - requesting the number of CPUs based on the average available
memory based on that node, and time

## Configuring processes

nf-core configurations tend to be too general and do not fit architecture. This
can slow down scheduling by requesting more resources than required (block other
users - HPC is a shared system, and underutilise node allocations.

Example trace files

=== "Gadi (PBS)"

    | name                       | status    | exit | duration | realtime | cpus | %cpu   | memory | %mem | rss      |
    | -------------------------- | --------- | ---- | -------- | -------- | ---- | ------ | ------ | ---- | -------- |
    | FASTQC (fastqc on NA12878) | COMPLETED | 0    | 24.3s    | 2s       | 4    | 213.9% | 2 GB   | 0.3% | 329 MB   |
    | FASTQC (fastqc on NA12889) | COMPLETED | 0    | 24.3s    | 2s       | 4    | 224.7% | 2 GB   | 0.3% | 337.9 MB |
    | ALIGN (1)                  | COMPLETED | 0    | 24.4s    | 0ms      | 4    | 262.6% | 2 GB   | 0.0% | 26.6 MB  |
    | FASTQC (fastqc on NA12877) | COMPLETED | 0    | 29.4s    | 3s       | 4    | 226.3% | 2 GB   | 0.2% | 322.6 MB |
    | GENOTYPE (1)               | COMPLETED | 0    | 1m 20s   | 50s      | 4    | 145.6% | 2 GB   | 1.0% | 1.2 GB   |
    | JOINT_GENOTYPE (1)         | COMPLETED | 0    | 39.9s    | 9s       | 4    | 225.6% | 2 GB   | 0.4% | 511.3 MB |
    | STATS (1)                  | COMPLETED | 0    | 29.9s    | 0ms      | 4    | 161.3% | 2 GB   | 0.0% | 4 MB     |
    | MULTIQC                    | COMPLETED | 0    | 34.9s    | 3.9s     | 4    | 100.7% | 2 GB   | 0.1% | 98.6 MB  |

=== "Setonix (Slurm)"
    | name                       | status    | exit | duration | realtime | cpus | %cpu   | memory | %mem | rss      |
    | -------------------------- | --------- | ---- | -------- | -------- | ---- | ------ | ------ | ---- | -------- |
    | FASTQC (fastqc on NA12878) | COMPLETED | 0    | 14.6s    | 4s       | 2    | 154.2% | 2 GB   | 0.1% | 250.5 MB |
    | FASTQC (fastqc on NA12889) | COMPLETED | 0    | 14.5s    | 4s       | 2    | 164.9% | 2 GB   | 0.1% | 261.5 MB |
    | ALIGN (1)                  | COMPLETED | 0    | 14.6s    | 1s       | 2    | 110.7% | 2 GB   | 0.0% | 96.4 MB  |
    | FASTQC (fastqc on NA12877) | COMPLETED | 0    | 14.6s    | 3s       | 2    | 142.4% | 2 GB   | 0.1% | 251.2 MB |
    | GENOTYPE (1)               | COMPLETED | 0    | 44.9s    | 32s      | 2    | 156.7% | 2 GB   | 0.7% | 1.7 GB   |
    | JOINT_GENOTYPE (1)         | COMPLETED | 0    | 20s      | 8s       | 2    | 227.5% | 2 GB   | 0.2% | 458.8 MB |
    | STATS (1)                  | COMPLETED | 0    | 10.3s    | 0ms      | 2    | 123.9% | 2 GB   | 0.0% | 2 MB     |
    | MULTIQC                    | COMPLETED | 0    | 19s      | 5.1s     | 2    | 72.5%  | 2 GB   | 0.0% | 86.7 MB  | 

!!! example "Exercise"

    Create a file `conf/custom.config` and copy the following code chunk.
    Save the file.

    ```groovy title="custom.config"
    process {
        withName: "FASTQC" {
            cpus = 2
            memory = 1.GB
            time = 2.min
        }

        withName: "ALIGN" {
            cpus = 1
            memory = 1.GB
            time = 5.min
        }
        
        withName: "GENOTYPE" {
            cpus = 
            memory = 
            time = 
        }
        
        withName: "JOINT_GENOTYPE" {
            cpus = 
            memory = 
            time = 
        }
        
        withName: "STATS" {
            cpus = 1
            memory = 1.GB
            time = 5.min
        }
        
        withName: "MULTIQC" {
            cpus = 1
            memory = 1.GB
            time = 5.min
        }
    }
    ```

    View trace, configure resources for GENOTYPE and JOINT_GENOTYPE

!!! example "Exercise"

    TODO Replace similar ones withLabel

Important: How to get "free" resources by correctly configuring to queue 
- Instructions to look at Gadi and Setonix queues/partitions

https://sydney-informatics-hub.github.io/training.gadi.intro/07-Optimisation/index.html

## Configuring to queues and partitions

Different queues/partitions are intended for different types of jobs

Provide a guided example

TODO excalidraw showing properties of different job types (e.g. long running
highmem vs. short and quick 

!!! example "Exercise"

    TODO Review Gadi queues, Setonix partitions

    === "Gadi (PBS)"

        Review [Cascade Lake queues](https://opus.nci.org.au/spaces/Help/pages/236881198/Queue+Limits...#QueueLimits...-CascadeLakequeuelimits) 
        
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

## Writing good custom scripts

Other things to consider - when writing custom R or Python scripts, writing
them efficiently. Utilising things like vectorisation, libraries such as numpy
etc., OpenMPI - link out.
