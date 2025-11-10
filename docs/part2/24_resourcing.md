# Assigning process resources

!!! info "Learning objectives"

    - Identify how to size resource requests appropriately for each process
    - Apply resource-aware design principles to improve job efficiency
    - Optimise processes for time, noting the trade-offs with cost (SU usage)
    - Understand how to efficiently configure jobs to fit system queues
    - Apply infrastructure requirements of memory/CPU proportions to match the node
    architecture of the target HPC queue/partition

## Efficiency

Optimising for time and throughput - refer to concepts like more nodes, less
walltime is more efficient than long walltime with less nodes (on Gadi, but
ensure to generalise to any HPC/scheduler). This can also be more
cost-effective using more resources for less time.

Requesting only the resources you need for the process can help ensure
yours (and others'!) jobs can be scheduled, and is scheduled appropriately
in to the correct queue or partition.

Like seating patrons in a restaurant or cafe example - smaller your group,
the likely you will get a table?

## A note on queue selection

In the workshop all the jobs require minimal resources and run quickly.
Hence, we can run them all on the same low-SU queues. When you configure
your own processes

Different queues/partitions are intended for different types of jobs

https://sydney-informatics-hub.github.io/usyd-gadi-onboarding-guide/notebooks/08_job_script.html#queue-selection-examples
    
Know how to set up system-specific config, how to ensure the resourcing aligns
well with the setup of the infrastructure.

HPC architecture differs across platforms and that the queue/partition
names and resources on that queue affect
the config files that needs to be created for that platform.

Tie back in that nextflow code can run on any platform, but when using HPC, the
config needs to be correct for that specific infrastructure.

Links to Setonix and Gadi docs, nf-core instutional configs.

- https://pawsey.atlassian.net/wiki/spaces/US/pages/51929058/Running+Jobs+on+Setonix
- https://opus.nci.org.au/spaces/Help/pages/236881198/Queue+Limits...#QueueLimits...-Broadwellqueuelimits

## Configuring processes

nf-core configurations that come with the pipelines tend to be too general and do not
fit architecture. This can slow down scheduling by requesting more resources than
required (block other users - HPC is a shared system, and can underutilise node allocations.

Example trace files we use to configure our resources.

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


We will now configure our scheduler-specific configs so it fits the node
infrastructure you are using. Recall that we specified extra CPUs to get
the workflow running:

=== "Gadi (PBS)"

    ```groovy title='custom.config
    process {
        cpu = 4 // 'queue' normalbw = 128 GB / 28 CPU ~ 4.6
        memory = 2.GB
    }
    ```

=== "Pawsey (Slurm)"

    ```groovy title='custom.config
    process {
        cpu = 2 // 'work' partition = 230 GB / 128 CPU ~ 1.8
        memory = 2.GB
    }
    ```

- Go to explain pages
- https://pawsey.atlassian.net/wiki/spaces/US/pages/51929058/Running+Jobs+on+Setonix
- https://opus.nci.org.au/spaces/Help/pages/236881198/Queue+Limits...#QueueLimits...-Broadwellqueuelimits

Identify the limiting resource. In the above case, it is memory -> find the
average CPUs that need to be allocated to this node.

!!! example "Exercise"

    Add the following to your process block.

    === "Gadi (PBS)"

        ```groovy title="conf/pbspro.config"
        process {
            withName: "FASTQC" {
                cpus = 
                memory = 2.GB
                time = 2.min
            }

            withName: "ALIGN" {
                cpus = 4
                memory = 2.GB
                time = 2.min
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
                cpus = 4
                memory = 2.GB
                time = 2.min
            }
            
            withName: "MULTIQC" {
                cpus = 4
                memory = 2.GB
                time = 2.min
            }
        }
        ```

    === "Setonix (Slurm)"

        ```groovy title="conf/slurm.config"
        process {
            withName: "FASTQC" {
                cpus = 2
                memory = 2.GB
                time = 2.min
            }

            withName: "ALIGN" {
                cpus = 2
                memory = 2.GB
                time = 2.min
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
                cpus = 2
                memory = 2.GB
                time = 2.min
            }
            
            withName: "MULTIQC" {
                cpus = 2
                memory = 2.GB
                time = 2.min
            }
        }
        ```

We will give `FASTQC` two CPUs to process each of the paired-end reads.
According to the trace, it does not require much memory, so the limiting
resource here is CPU.

Finding the effective usable RAM/core.

!!! example "Exercise:

    Refering to the HPC queue/partition documentation, how much memory
    should you allocate given the `cpus = 2`?

    ??? tip

        - Divide the usuable RAM on that queue/partition, by the highest
        number of CPUs.
        - Multiply that value with 2 CPUs
        - Round up/down.

    === "Gadi (PBS)"

        Review [Queue Limits](https://opus.nci.org.au/spaces/Help/pages/236881198/Queue+Limits...#QueueLimits...)
        for `normalbw`.

        A: 128GB/28CPU ~ 4.6GB per CPU
        A: 4.6GB x 2 CPU required = 9.2
        A: 9 GB memory

    === "Pawsey (Slurm)"

        Review [partitions](https://pawsey.atlassian.net/wiki/spaces/US/pages/51929058/Running+Jobs+on+Setonix)
        for `work`

        A: 230GB/128CPU ~ 1.8GB per CPU
        A: 1.8GB x 2 CPU required = 3.6
        A: 4 GB memory


!!! example "Exercise"

    === "Gadi (PBS)"

        ```groovy title="conf/pbspro.config"
        process {
            withLabel "small" {
                cpus = 2
                memory = 4.GB
            }
            ...
        }
        ```

    === "Pawsey (Slurm)"

        ```groovy title="conf/slurm.config"
        process {
            withLabel "small" {
                cpus = 2
                memory = 4.GB
            }
            ...
        }
        ```

Takeaway: Specifying the number of resources is the first step of
ensuring you don't ask for resources you don't need. On systems with
a lot of freedom (cloud instances, workstations) this is sufficient.

However on shared HPC systems, we need to be more explicit with what
we can use. Providing the extra resources can provide extra processing
power in comparison to being stringent.

## Configuring java heap sizes

!!! example "Exercises"

    Update GENOTYPE and JOINT_GENOTYPE processes with -Xmx${tasks.memory}

## Writing good custom scripts

Other things to consider - when writing custom R or Python scripts, writing
them efficiently. Utilising things like vectorisation, libraries such as numpy
etc., OpenMPI - link out.
