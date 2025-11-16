# Parallelisation

!!! info "Learning objectives"

    - Consider the limitations of parallelisation and scenarios where it should
    not be applied
    - Differentiate between multithreading and scatter-gather approaches to
    parallelisation
    - Implement parallelisation approaches (multithreading, multiprocessing)
    while preserving biological correctness
    - Recall that faster jobs that use more resources are more efficient than
    long-running jobs with less resources

Part 1.4 introduced parallelisation approaches with the goal of speeding up your jobs by utilising more resources. As your data gets larger, or more samples are required to be processed, it needs to run efficiently. In this lesson we will explore how Nextflow supports different forms of parallelisation to help you scale your workflows.

![](figs/00_benchmark_at_scale_theme.png)

Recall that splitting your data up across too many cores can lead to diminishing returns, such as increased SU usage and walltime. Parallelisation requires benchmarking to find the right balance between the walltime, CPU efficiency, and service unit consumption. We want to avoid over-parallelising our workflows.

## Multithreading `bwa mem`

In this section we will look at implementing another multithreading example with `bwa mem`, used in the `ALIGN` process. These are the example benchmarking results from Part 1, with the CPU efficiency calculated for you. 

| Cores | Walltime (s) | CPU time (s) | CPU efficiency |
|-------|--------------|--------------|----------------|
| 2     | 0.744        | 1.381        | 93%            |
| 4     | 0.456        | 1.537        | 84%            |
| 6     | 0.355        | 1.618        | 76%            |
| 8     | 0.291        | 1.628        | 70%            |

These values were taken from the `time` command output:

- Walltime = `real`
- CPU time = `user` + `sys`
- CPU efficiency = `cpu time / (cpu time * cores)`

!!! info "CPU efficiency"

    Recall that CPU efficiency is a measure of how many cpus were actually used, in comparison to how many cpus were requested. A high CPU efficiency (100%) means that all of the CPUs were utilised, while a low efficiency suggests that too many were requested.

In this example, we have to consider the trade offs between each run and what we would like to optimise for. 

Providing 2 cores has the slowest walltime but utilises the 2 CPUs efficiently (93%).

On the other hand, providing 8 cores provides ~40% speed up in walltime with reduced CPU efficiency.

As responsible users of shared systems, we will select the option that maintains high CPU efficiency. While this is not prescriptive, aim for >80% CPU efficiency ensures we are not reserving resources in excess, that others' can use for their own jobs.

!!! question "Poll"

    TODO: Add first question to poll

    1. How many cores would you choose to provide `ALIGN` to ensure that it still uses the CPUs efficiently, but with a speed up in walltime? 
    2. Which `.config` file would you want to use? (`nextflow.config` - workflow-specific, system-agnostic; `custom.config` - workflow-specific, system-specific)
    3. How much extra memory can you utilise if required? (Consider the effective RAM/CPUs
    proportion of the queue or partition)

    === "Gadi (PBS)"

        ??? note "Answers"

            - 4 CPUs has > 80% CPU efficiency  
            - `custom.config` to ensure it is tuned for the `normalbw` queue.
            - 4 CPUs with 16-18 GB memory

    === "Setonix (Slurm)"

        ??? note "Answers"

            - 4 CPUs has > 80% CPU efficiency  
            - `custom.config` to ensure it fits the `work` partition.
            - 4 CPUs with 7-8 GB memory

!!! example "Exercise"

    TODO untruncate

    In `conf/custom.config`, add the following inside the `process` scope:

    === "Gadi (PBS)"

        ```groovy title="conf/custom.config" hl_lines="3-8"
        process {
        // truncated

        	withName: 'ALIGN' {
        		cpus = 4 // Efficient number of cores
        		memory = 18.GB // 4 CPUs * 4.6 GB
        		time = 2.min
        	}
        }
        ```

    === "Setonix (Slurm)"

        ```groovy title="conf/custom.config" hl_lines="3-8"
        process {
        // truncated

        	withName: 'ALIGN' {
        		cpus = 4 // Efficient number of cores
        		memory = 7.GB // 4 CPUs * 1.8 GB
        		time = 2.min
        	}
        }
        ```

    Note: our `ALIGN` process (`modules/align.nf`) has `-t $task.cpus` already defined, so you do not need to amend it.

    Save your file and run with:
    
    ```
    ./run.sh
    ```

!!! info "Remember to read the tool documentation!"

    All software and bioinformatics tools are all built differently. Some support multi-threading, some can only run things with a single thread. Overlooking these details may not be crucial when running on systems where you have autonomy and access to all resources (personal compute, cloud instances), however, these are important parts of configuring your workflow on HPC shared systems to set reasonable limits and requests.

## Scatter-gathering alignment

![](figs/00_Scatter_gather_fig.png)

One of the core benefits of running bioinformatics workflows on HPC is access to increased processing power and hardware. For jobs that can be conducted indepdenently of each other, if configured correctly, we can run many jobs simultaneously and reduce the overall walltime required to run the workflow. One strategy to implement this is by:

1. Splitting/scattering the data
2. Processing each of the data chunks separately
3. Combining/gathering the processed outputs back into a single file

!!! note "Not everything can or should be split"

    Recall from Part 1 that we can't split everything - it should only be done if the particular processing step can be conducted independently of each other. Scattering taks does not make sense when results depend on comparing all data together, such as detecting structural variants across multiple chromosomes.

We will scatter-gather the alignment step. This is a widely approach for mapping reads, as whole genome data is large, can be time-consuming, and mapping can be conducted independently of each other. To do so, we will leverage Nextflow's built-in [`splitFastq`](https://www.nextflow.io/docs/latest/reference/operator.html#splitfastq) operator.

!!! example "Exercise"

    Add the following chunk to your `main.nf` file:

    ```groovy title="main.nf"
    // Split FASTQs for each sample
    split_fqs = reads
        .splitFastq(limit: 3, pe: true, file: true)
        .view()
    ```

- The `reads` channel is taken as input. It contains the `[ sample_name, fastq_r1, fastq_r2 ]`
- `.splitFastq` splits each paired `.fastq` file (`pe: true`) into three files (`limit: 3`)
- `file: true` stores each split `.fastq` file in the work directory and avoids out-of-memory issues
- We include `.view()` to inspect the contents of the `split_fqs` channel we just created

Next, we need to update the inputs to `ALIGN`, so it takes the split `.fastq` files.

!!! example "Exercise"

    In `main.nf`, in the `workflow` scope, replace the input argument to `ALIGN` from `ALIGN(reads, bwa_index)` to `ALIGN(split_fqs, bwa_index)`.

    ```groovy title="main.nf"
    // Split FASTQs for each sample
    split_fqs = reads
        .splitFastq(limit: 3, pe: true, file: true)
        .view()

    ALIGN(split_fqs, bwa_index)
    ```

    Save the file, and run:
    ```bash
    ./run.sh
    ```

    ??? abstract Show output

        ```console title="Output"
        [NA12877, /scratch/pawsey1227/fjaya/nextflow-on-hpc-materials/part2/work/86/7ca90d09bee8f0f952b82e0b791d66/NA12877_chr20-22.R1.1.fq, /scratch/pawsey1227/fjaya/nextflow-on-hpc-materials/part2/work/c8/f6d055eae92f2336aa98b2d2790623/NA12877_chr20-22.R2.1.fq]
        [NA12877, /scratch/pawsey1227/fjaya/nextflow-on-hpc-materials/part2/work/86/7ca90d09bee8f0f952b82e0b791d66/NA12877_chr20-22.R1.2.fq, /scratch/pawsey1227/fjaya/nextflow-on-hpc-materials/part2/work/c8/f6d055eae92f2336aa98b2d2790623/NA12877_chr20-22.R2.2.fq]
        [NA12877, /scratch/pawsey1227/fjaya/nextflow-on-hpc-materials/part2/work/86/7ca90d09bee8f0f952b82e0b791d66/NA12877_chr20-22.R1.3.fq, /scratch/pawsey1227/fjaya/nextflow-on-hpc-materials/part2/work/c8/f6d055eae92f2336aa98b2d2790623/NA12877_chr20-22.R2.3.fq]
        executor >  slurm (8)
        [f1/dea318] process > FASTQC (fastqc on NA12877) [100%] 1 of 1 ✔
        [79/d85144] process > ALIGN (2)                  [100%] 3 of 3 ✔
        [4d/c80183] process > GENOTYPE (1)               [100%] 1 of 1 ✔
        [63/8a26ab] process > JOINT_GENOTYPE (1)         [100%] 1 of 1 ✔
        [61/5157d0] process > STATS (1)                  [100%] 1 of 1 ✔
        [1e/b0a7c0] process > MULTIQC                    [100%] 1 of 1 ✔
        Completed at: 16-Nov-2025 21:24:28
        Duration    : 1m 47s
        CPU hours   : (a few seconds)
        Succeeded   : 8
        ```

Let's take a look at the stdout printed.

!!! tip "Scatter-gather patterns"

    How you implement the scatter-gather pattern in Nextflow will be highly dependent on your workflow structure, and input and output files. [Nextflow patterns](https://nextflow-io.github.io/patterns/) provides examples of commonly used patterns that support a range of different needs, such as splitting text and CSV files, and collecting outputs multiple outputs into a single file or groups.

!!! note

    Recall that modules are useful to keep things modular and avoid cluttering our main.nf! It is far easier to swap out the module imports, in comparison to deleting the process definitions or commenting them out. We may want to reuse our revert back to our original processes too. Modules, combined with the workflow definition and groovy operators are what allow applying scatter-gather patterns (multiprocessing) to your Nextflow with ease.

    NOTE: The alignment module was updated to use a scatter-gather approach. Instead of aligning the entire FASTQ in one go with the ALIGN module, the workflow now:

    1. Splits the FASTQ into 3 chunks (SPLIT_FASTQ).
    2. Aligns each chunk in parallel (ALIGN_CHUNK, see the 3 of 3 completed).
    3. Merges the aligned chunks into a single BAM (MERGE_BAMS).

    This change optimises performance for large datasets by leveraging parallel processing.

!!! What about multiple samples?

    You have now applied a multi-processing approach on a single-sample. As processes are run independently of each other this does not always need to apply to single sample that is split. Running multiple samples is also a form of multi-processing and comes shipped with Nextflow's dataflow model. Once your pipeline is configured to run well with a single sample, [queue channels](https://sydney-informatics-hub.github.io/hello-nextflow-2025/part1/05_inputs/#queue-channels) make adding additional samples relatively easy.

    We will revisit this in the next section.

## A note on dynamic resourcing

Since our data is small and similar-sized, we can apply the same resource configurations within the same process and it will still run successfully. However, it is common that we need to **run the same process with input data of widely variying sizes**. For example, if we were to run variant calling with reads from the whole genome, human chromosome 1 is nearly 4x larger than chromosome 20.

One option may be to configure the resource usage so it runs on the largest data (chr. 1). This will ensure all processes run sucessfully at the cost of **vastly underutilising the resources you have requested** on the smaller data:

- Excess resources will be reserved other users could access
- Your jobs stay in queue for longer
- You can be charged more SUs than required

The alternative would be to take a dynamic approach when you need to process data that require unequal CPUs, memory, or walltime.

If a job fails, you can tell Nextflow to automatically re-run it with additional resources. For example:

| Directive | Closure Example             | Attempt 1 (Initial Run) | Attempt 2 (First Retry) |
| --------- | --------------------------- | ----------------------- | ----------------------- |
| `memory`  | `{ 2.GB * task.attempt }`   | 2 GB                    | 4 GB                    |
| `time`    | `{ 1.hour * task.attempt }` | 1 hour                  | 2 hours                 |

Another approach is to ynamically assign a resource based on properties of the input data. For example, by the size of the file:

```groovy
process {
  withName: 'ALIGN' {
    memory = { reads.size() * 2 }
  }
}
```

For more information, see Nextflow's training on:

- [Retry strategies](https://training.nextflow.io/2.1.1/advanced/configuration/#retry-strategies)
- [Dynamic directives](https://training.nextflow.io/2.1.1/advanced/configuration/#dynamic-directives)

The same concepts of configuring resources will apply here, aim to fit the
appropriate queues/partitions, but will require addtional benchmarking. This is worthwhile if developing and running high-throughput workflows.

## Summary

Recall that we do not always want to parallelise everything. There are reasons
to avoid this due to the data set or analysis requiring data be analysed as a
whole, or computational reasons where processes become less efficient.

These are some strategies to consider for your own pipelines. Run benchmarks,
identify long-running or inefficent processes and consider which one of these
approaches are supported by the tools (multi-threading), or can be split and
processed in parallel.

TODO add final code changes