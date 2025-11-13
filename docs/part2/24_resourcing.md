# Assigning process resources

!!! info "Learning objectives"

    - Recall how to size resource requests appropriately for each process
    - Apply resource-aware design principles to improve job efficiency
    - Understand how to efficiently configure jobs to fit system queues
    - Apply infrastructure requirements of memory/CPU proportions to match the node
    architecture of the target HPC queue/partition

## Requesting resources efficiently

When running workflows on an HPC, the way you configure your resources (number of CPUs, the amount of memory, and the expected walltime) can make a big difference to how efficiently your jobs run, how long they sit in the queue, and even how much service units (SUs) the job costs.

A common assumption is that requesting fewer resources is safer or more considerate. In practice,
this may not be the best approach. Jobs that are under-resourced often run for much longer than they need to, or fail partway through. On the other hand, using **more CPUs and memory for a shorter time can often be more efficient and more cost-effective**. This principle applies broadly across HPC systems.

The scheduler also plays a role here. HPC schedulers look for jobs that “fit” into available space. If your job is small and short, it might slide into a gap easily, but if it’s too large, or asks for more time than necessary, it might be left waiting. In general, jobs that request just enough resources to run efficiently and cleanly tend to move through the system faster and cause fewer problems for other users.

!!! tip 

    You can think of it like finding a table at a busy café. If you’re alone or in a small group, it’s easier to get seated quickly, but if you camp out all day at a large table and only order a coffee, you’re not making good use of the space. The HPC scheduler works in a similar way: it tries to fill in available “seats” with jobs that fit neatly and will clear out efficiently.

In this section, we’ll explore how to make smart decisions about resource configuration such as:

* How many resources your process actually needs
* How to estimate the average number of usable memory per CPU on a node to access additional resources with minimal drawbacks

These adjustments don’t just benefit your own workflow, they make better use of shared infrastructure.

## Configuring processes

To begin tuning our workflow, we first need to understand how many resources each process actually used. We will use these example trace summaries, generated in the previous lesson, as a baseline for how long each process ran, and how many CPUs and memory it used.

=== "Gadi (PBS)"

    | name                       | status    | exit | duration | realtime | cpus | %cpu  | memory | %mem | peak_rss |
    | -------------------------- | --------- | ---- | -------- | -------- | ---- | ----- | ------ | ---- | -------- |
    | ALIGN (1)                  | COMPLETED | 0    | 29.6s    | 1s       | 1    | 93.7% | 4 GB   | 0.0% | 95.8 MB  |
    | FASTQC (fastqc on NA12877) | COMPLETED | 0    | 34.5s    | 5s       | 1    | 76.2% | 4 GB   | 0.1% | 286.6 MB |
    | GENOTYPE (1)               | COMPLETED | 0    | 59.9s    | 45s      | 1    | 97.6% | 4 GB   | 0.5% | 950.3 MB |
    | JOINT_GENOTYPE (1)         | COMPLETED | 0    | 34.8s    | 16s      | 1    | 93.3% | 4 GB   | 0.3% | 508.8 MB |
    | STATS (1)                  | COMPLETED | 0    | 19.9s    | 0ms      | 1    | 73.4% | 4 GB   | 0.0% | 3.1 MB   |
    | MULTIQC                    | COMPLETED | 0    | 29.9s    | 4.7s     | 1    | 79.5% | 4 GB   | 0.0% | 97.2 MB  |

=== "Setonix (Slurm)"

    | name                       | status    | exit | duration | realtime | cpus | %cpu   | memory | %mem | peak_rss |
    | -------------------------- | --------- | ---- | -------- | -------- | ---- | ------ | ------ | ---- | -------- |
    | FASTQC (fastqc on NA12877) | COMPLETED | 0    | 13.8s    | 4s       | 1    | 135.0% | 2 GB   | 0.1% | 240.6 MB |
    | ALIGN (1)                  | COMPLETED | 0    | 13.8s    | 2s       | 1    | 100.1% | 2 GB   | 0.0% | 98.2 MB  |
    | GENOTYPE (1)               | COMPLETED | 0    | 39.9s    | 28s      | 1    | 164.8% | 2 GB   | 0.5% | 1.4 GB   |
    | JOINT_GENOTYPE (1)         | COMPLETED | 0    | 19.2s    | 8s       | 1    | 204.2% | 2 GB   | 0.2% | 466 MB   |
    | STATS (1)                  | COMPLETED | 0    | 14.9s    | 1s       | 1    | 45.2%  | 2 GB   | 0.0% | 2 MB     |
    | MULTIQC                    | COMPLETED | 0    | 19.9s    | 5.3s     | 1    | 62.4%  | 2 GB   | 0.0% | 78.6 MB  |

While we could configure each process to match these values, we’re instead going to take a broader view. We'll explore how HPC systems allocate memory per CPU, and how to align our process requests to match this architecture more effectively.

This approach lets us make the most of the available resources - sometimes even getting "extra" memory at no additional cost and still remain scheduled quickly.

## Effective memory per CPU

Recall that this is the configuration we used in Part 2.1 to get the pipeline running on the respective HPC systems:

=== "Gadi (PBS)"

    ```groovy title='custom.config'
    process {
        cpu = 1 // 'normalbw' queue = 128 GB / 28 CPU ~ 4.6
        memory = 4.GB
    }
    ```

=== "Pawsey (Slurm)"

    ```groovy title='custom.config'
    process {
        cpu = 1 // 'work' partition = 230 GB / 128 CPU ~ 1.8
        memory = 2.GB
    }
    ```
    
The reason for these values come down to requesting a balanced amount of memory, relative to your CPU.

![](../part1/figs/00_smarter_proportions.png)

Most HPC systems allocate jobs to nodes based on both CPU and memory requests, where each queue or partition is associated with a specific type of node with: 

- A fixed number of CPUs
- A fixed amount of memory

This means there is an **average amount of memory per CPU** - this becomes an important consideration for optimising resource requests for your Nextflow processes.

## Exploring resource options for FASTQC

Based on the Part 1 discussion on the optimal number of `FASTQC` threads, we will give `FASTQC` two CPUs to process each of the paired-end reads. According to the trace file, it does not require much memory, so the limiting resource here is CPU.

Let's find the effective usable memory per CPU for the `normalbw` queue on Gadi, and the `work` partition on Setonix.

Note: Some queues and partitions may have different nodes with different hardware specifications, such as nodes with more or less memory. Using the effective memory per cpu to fit the smaller node has the benefit of running on both, whereas the larger node may take longer, but provides extra resources.

!!! example "Exercise"

    Refering to the HPC queue/partition documentation, how much memory
    should you allocate given a process that requires `cpus = 2`?

    Steps:

    - Divide the usuable RAM on that queue/partition, by the highest
    number of CPUs.
    - Multiply that value with 2 CPUs
    - Round up/down.

    === "Gadi (PBS)"

        Review [Queue Limits](https://opus.nci.org.au/spaces/Help/pages/236881198/Queue+Limits...#QueueLimits...)
        for `normalbw`.

        ??? note "Answer"

            - 128GB/28CPU ~ 4.6GB per CPU
            - 4.6GB x 2 CPU required = 9.2
            - 9 GB memory

            OR

            - 18 GB to fit 256GB nodes

    === "Setonix (Slurm)"

        Review [partitions](https://pawsey.atlassian.net/wiki/spaces/US/pages/51929058/Running+Jobs+on+Setonix)
        for `work`

        ??? note "Answer"

            - 230GB/128CPU ~ 1.8GB per CPU
            - 1.8GB x 2 CPU required = 3.6
            - 4 GB memory

Now that we’ve calculated the effective memory per CPU on each system, the next consideration is - how much should we actually request for a small, fast process like `FASTQC`?

Let’s look at the trade-offs.

If we request:

    === "Gadi (PBS Pro)"

    - 2 CPUs and 8 GB memory (on the Gadi `normalbw` partition), this takes advantage of all the memory you’re entitled to, but `FASTQC` won't actually use that memory. So you're not getting any extra performance and may lengthen the time in queue.
    - 2 CPUs and 1 GB memory, on the other hand, still gives `FASTQC` enough to run, and because you're requesting less RAM, your job may be scheduled fasterm as it can fit into more available nodes. This is more memory efficient too.

    === "Setonix (Slurm)"

    - 2 CPUs and 4 GB memory (on the Setonix `work` partition), this takes advantage of all the memory you’re entitled to, but `FASTQC` won't actually use that memory. So you're not getting any extra performance and may lengthen the time in queue.
    - 2 CPUs and 1 GB memory, on the other hand, still gives `FASTQC` enough to run, and because you're requesting less RAM, your job may be scheduled fasterm as it can fit into more available nodes. This is more memory efficient too.
        
We will proceed with the 2 CPUs 1 GB memory option for `FASTQC` as the job won't benefit from the extra memory.

### Configuring with process names (`withName`)

Recall that configuration can be workflow-specific to run across different systems, but needs to be system-specific to use the infrastructure effectively. Here, we will apply resource configuration settings in the `config/custom.config` as we are tuning according to the HPC we are running it on.

`withName` is a powerful tool to:
    - Specifically targets individual modules
    - Specify multiple module names using wildcards (`*` or `|`)
    - Avoid editing the module.nf file to add a process label (remember: separation of workflow logic and system-specific tuning)
    - Has a higher priority than withLabel

!!! question "Configuring `withLabel` and configuration priorities"

    Processes that require similar resources can be configured using the `withLabel` process selector. Processes can be tagged with multiple labels to flexibly categorise different runtime needs (e.g. high memory). However, `withLabel` may be overwritten by settings defined by `withName`.

    For more information, see [custom resource configuration using process labels](https://sydney-informatics-hub.github.io/customising-nfcore-workshop/notebooks/2.3_configEnv.html#custom-resource-configuration-using-process-labels).

To summarise and group the resource usage from the trace file:  

| Process                         | Resources | Rationale                                        |
| ------------------------------- | --------- | ------------------------------------------------ |
| FASTQC                          | 2CPU, 1GB | Fix 2CPU to process R1 and R2, memory sufficient |
| ALIGN, GENOTYPE, JOINT_GENOTYPE | 2CPU, 1GB | High CPU utilisation >90%                        |
| GENOTYPE                        | 2CPU, 2GB |                                                  |
| STATS, MULTIQC                  | 1CPU, 2GB |                                                  |

Let's record these in our configs, with a bit of buffer so things don't fail.

Note that we add them to the config file, and not the modules. This keeps the workflow logic and system-specific configuration separate.

!!! example "Exercise"

    === "Gadi (PBS)"

        ```groovy title="conf/custom.config"
        process {
            cpu = 1 // 'normalbw' queue = 128 GB / 28 CPU ~ 4.6
            memory = 4.GB
        
            withName: /FASTQC|ALIGN|JOINT_GENOTYPE/ {
                cpus = 2
                memory = 1.GB
                time = 2.minutes
            }
    
            withName: GENOTYPE {
                cpus = 2
                memory = 2.GB
                time = 2.minutes
            }

            withName: /STATS|MULTIQC/ {
                cpus = 1
                memory = 2.GB
                time = 2.minutes
            }

        }
        ```

    === "Setonix (Slurm)"

        ```groovy title="conf/custom.config"
        process {
            cpu = 1 // 'work' partition = 230 GB / 128 CPU ~ 1.8
            memory = 2.GB

            withName: /FASTQC|ALIGN|JOINT_GENOTYPE/ {
                cpus = 2
                memory = 1.GB
                time = 2.minutes
            }
    
            withName: GENOTYPE {
                cpus = 2
                memory = 2.GB
                time = 2.minutes
            }

            withName: /STATS|MULTIQC/ {
                cpus = 1
                memory = 2.GB
                time = 2.minutes
            }
        }
        ```

    ```bash
    ./run.sh
    ```

Review the new trace file. What has changed? What happened to our FASTQC process? Faster or slower?

=== "Gadi (PBS)"

    | name                       | status    | exit | duration | realtime | cpus | %cpu   | memory | %mem | peak_rss |
    | -------------------------- | --------- | ---- | -------- | -------- | ---- | ------ | ------ | ---- | -------- |
    | FASTQC (fastqc on NA12877) | COMPLETED | 0    | 49.6s    | 3s       | 2    | 136.4% | 1 GB   | 0.1% | 213.8 MB |
    | ALIGN (1)                  | COMPLETED | 0    | 49.5s    | 1s       | 2    | 122.4% | 1 GB   | 0.1% | 102.5 MB |
    | GENOTYPE (1)               | COMPLETED | 0    | 1m 20s   | 36s      | 2    | 136.3% | 2 GB   | 0.8% | 1.5 GB   |
    | JOINT_GENOTYPE (1)         | COMPLETED | 0    | 49.8s    | 9s       | 2    | 163.8% | 1 GB   | 0.2% | 391.6 MB |
    | STATS (1)                  | COMPLETED | 0    | 34.9s    | 0ms      | 1    | 62.0%  | 2 GB   | 0.0% | 3.1 MB   |
    | MULTIQC                    | COMPLETED | 0    | 49.9s    | 5.6s     | 1    | 67.3%  | 2 GB   | 0.0% | 93.5 MB  |

=== "Setonix (Slurm)"

    | name                       | status    | exit | duration | realtime | cpus | %cpu   | memory | %mem | peak_rss | 
    | -------------------------- | --------- | ---- | -------- | -------- | ---- | ------ | ------ | ---- | -------- | 
    | FASTQC (fastqc on NA12877) | COMPLETED | 0    | 14s      | 4s       | 2    | 145.9% | 1 GB   | 0.1% | 259 MB   | 
    | ALIGN (1)                  | COMPLETED | 0    | 14.1s    | 1s       | 2    | 139.3% | 1 GB   | 0.0% | 2 MB     | 
    | GENOTYPE (1)               | COMPLETED | 0    | 44.8s    | 28s      | 2    | 164.0% | 2 GB   | 0.5% | 1.3 GB   | 
    | JOINT_GENOTYPE (1)         | COMPLETED | 0    | 19.4s    | 9s       | 2    | 211.3% | 1 GB   | 0.1% | 343.3 MB | 
    | STATS (1)                  | COMPLETED | 0    | 14.5s    | 0ms      | 1    | 132.7% | 2 GB   | 0.0% | 2 MB     | 
    | MULTIQC                    | COMPLETED | 0    | 15.3s    | 4.3s     | 1    | 77.6%  | 2 GB   | 0.0% | 76.4 MB  | 


## Passing allocated resources into process scripts

Nextflow allows you to request specific resources for each process (like CPU and memory), but this doesn’t automatically tell the tool inside the process to use those resources. Some bioinformatics tools require you to explicitly specify how much memory or how many threads to use inside the script block, or uses a suboptimal default.

If this is not done, **the scheduler may assign the resources, but the process itself might underutilise them** - reducing efficiency and potentially wasting queue time.

If you do not update this, the resources will be allocated but not all
completely utilised. In our current pipeline the memory is hardcoded for
processes `GENOTYPE` and `JOINT_GENOTYPE` with 4 GB of memory.

In our current workflow, the `GENOTYPE` and `JOINT_GENOTYPE` processes hardcode the Java heap size (-Xmx4g). This is not ideal, because:

- It ignores any values set in the config file (custom.config)
- It prevents automatic adjustment if more memory is allocated
- It can lead to unnecessary limits or, worse, unexpected failures

Instead, we can make these values dynamic using `task.memory.toGiga()`. This allows the process to always match the memory allocated by the scheduler.

In this case the size of the data will not impact the memory usage drastically
as the previous trace files show that both processes use less than 4 GB. In
your own pipelines, a task with high memory utilisation may benefit from the
added memory allocated, and reduce walltime.

!!! example "Exercises"

    For both systems, replace the hardcoded `"-Xmx4g"` setting with the memory
    allocated to that task using `"-Xmx${task.memory.toGiga()}g"` for:

    1. Process `GENOTYPE` in `modules/genotype.nf`.

    ```groovy title='modules/genotype.nf'
    process GENOTYPE {
    // truncated

        output:
        tuple val(sample_id), path("${sample_id}.g.vcf.gz"), path("${sample_id}.g.vcf.gz.tbi"), emit: gvcf

        script:
        """
        gatk --java-options "-Xmx${task.memory.toGiga()}g" HaplotypeCaller -R $ref_fasta -I $bam -O ${sample_id}.g.vcf.gz -ERC GVCF
        """

    }
    ```

    2. Process `JOINT_GENOTYPE` in `modules/joint_genotype.nf`

    ```groovy title='modules/joint_genotype.nf'
    process JOINT_GENOTYPE {
    // truncated

        output:
        tuple val(cohort_id), path("${cohort_id}.vcf.gz"), path("${cohort_id}.vcf.gz.tbi"), emit: vcf

        script:
        variant_params = gvcfs.collect { f -> "--variant ${f}" }.join(" ")
        """
        gatk --java-options "-Xmx${task.memory.toGiga()}g" CombineGVCFs -R $ref_fasta $variant_params -O cohort.g.vcf.gz
        gatk --java-options "-Xmx${task.memory.toGiga()}g" GenotypeGVCFs -R $ref_fasta -V cohort.g.vcf.gz -O cohort.vcf.gz
        """
    }
    ```

!!! example "Exercises"

    Once updated, run your newly configured pipeline using.
    
    ```bash
    ./run.sh
    ```

    Once complete, inspect the new trace file.

You should see no major changes in memory usage or efficiency - but now your script matches the config, and will adjust based on any changes to the memory allocation.

!!! info "Writing efficient custom scripts"

    Beyond resource flags, it's also important to write efficient code inside your processes. If you're writing custom scripts (e.g. in Python or R):

    - Prefer vectorised operations over loops
    - Use optimised libraries like `numpy` for Python scripts
    - Consider parallelisation strategies (e.g. OpenMP)
    - Avoid holding large objects in memory unnecessarily

    While these are outside the scope of this workshop, they’re essential if you want to scale up workflows on HPC.

## Summary

Specifying the number of resources is the first step of
ensuring you don't ask for resources you don't need. On systems with
a lot of freedom (cloud instances, workstations) this is sufficient.

However on shared HPC systems, we need to be more explicit with what
we can use. Providing the extra resources can provide extra processing
power for supported tools, in comparison to being stringent.

Whilst the way the data is processed stays the same, it is important to review how your tools work (read their documentation!) and whether they can utilise extra resources.