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

As responsible users of shared systems, we will select the option that maintains high CPU efficiency. While this is not prescriptive, aiming for >80% CPU efficiency ensures we are not reserving resources in excess, that others' can use for their own jobs.

!!! question "Poll"

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

### Scatter: splitting our reads

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
         N E X T F L O W   ~  version 24.10.0

        Launching `main.nf` [suspicious_wiles] DSL2 - revision: de5d65b946

        [47/e84f4d] FASTQC (fastqc on NA12877) [100%] 1 of 1 ✔
        [c5/6d69ea] ALIGN (3)                  [100%] 3 of 3 ✔
        [a7/7a6424] GENOTYPE (2)               [100%] 3 of 3 ✔
        [47/e84f4d] FASTQC (fastqc on NA12877) [100%] 1 of 1 ✔
        [c5/6d69ea] ALIGN (3)                  [100%] 3 of 3 ✔
        [a7/7a6424] GENOTYPE (2)               [100%] 3 of 3 ✔
        [-        ] JOINT_GENOTYPE             -
        [-        ] STATS                      -
        [-        ] MULTIQC                    -
        [NA12877, /scratch/pawsey1227/fjaya/nextflow-on-hpc-materials/part2/work/a0/569c8d068367a6f922be0841dce142/NA12877_chr20-22.R1.1.fq, /scratch/pawsey1227/fjaya/nextflow-on-hpc-materials/part2/work/9e/c4cab021b4ecfa15e1f9a059ffd8e7/NA12877_chr20-22.R2.1.fq]
        [NA12877, /scratch/pawsey1227/fjaya/nextflow-on-hpc-materials/part2/work/a0/569c8d068367a6f922be0841dce142/NA12877_chr20-22.R1.2.fq, /scratch/pawsey1227/fjaya/nextflow-on-hpc-materials/part2/work/9e/c4cab021b4ecfa15e1f9a059ffd8e7/NA12877_chr20-22.R2.2.fq]
        [NA12877, /scratch/pawsey1227/fjaya/nextflow-on-hpc-materials/part2/work/a0/569c8d068367a6f922be0841dce142/NA12877_chr20-22.R1.3.fq, /scratch/pawsey1227/fjaya/nextflow-on-hpc-materials/part2/work/9e/c4cab021b4ecfa15e1f9a059ffd8e7/NA12877_chr20-22.R2.3.fq]
        ERROR ~ Error executing process > 'JOINT_GENOTYPE (1)'

        Caused by:
          Process `JOINT_GENOTYPE` input file name collision -- There are multiple input files for each of the following file names: NA12877.g.vcf.gz, NA12877.g.vcf.gz.tbi



        Container:
          /scratch/pawsey1227/fjaya/nextflow-on-hpc-materials/singularity/quay.io-biocontainers-gatk4-4.6.2.0--py310hdfd78af_1.img

        Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

         -- Check '.nextflow.log' file for details
        ```

Let's take a look at the `stdout` printed.

The output of `reads.splitFastq()` include three separate arrays that contain:

- The name of the sample
- The path to the R1 `.fastq` file
- The path to the R2 `.fastq` file

Note that each `.fastq` file is now identified with a chunk number (e.g. `.../NA12877_chr20-22.R2.1.fq`) - we have successfully split the reads into three.

The following line indicates that the `ALIGN` and `GENOTYPE` processes now run three times, successfully:

```console
        [c5/6d69ea] ALIGN (3)                  [100%] 3 of 3 ✔
        [a7/7a6424] GENOTYPE (2)               [100%] 3 of 3 ✔
```

!!! note

    Three separate tasks were automatically generated and scheduled by Nextflow without any extra instructions from us. Because we have already configured `ALIGN` and `GENOTYPE` in our config files, all scattered tasks used those resource settings automatically.

    Setting up your channels and using groovy operators can be a bit tricky at first, however, once these are set up correctly, Nextflow will take care of the scatter–gather orchestration for you. This makes it straightforward to parallelise work at scale with minimal additional code.

However, you should have received an error before `JOINT_GENOTYPE` was run:

```console title="Output"
 Caused by:
          Process `JOINT_GENOTYPE` input file name collision -- There are multiple input files for each of the following file names: NA12877.g.vcf.gz, NA12877.g.vcf.gz.tbi
```

Let's troubleshoot by inspecting the output of the `GENOTYPE` process

!!! example "Exercise"
    
    First, inspect the outputs of the `ALIGN` process by adding `.view()`.

    ```groovy title="main.nf" hl_lines="3"
    // Run genotyping with aligned bam and genome reference
    GENOTYPE(ALIGN.out.aligned_bam, ref)
    GENOTYPE.out.view()
    ```

    Save the file, and run with `-resume`

    ```
    ./run.sh 
    ```

    ??? abstract Output

        ```console title="Output"
        # workdirs have been truncated with '...' for readability
        [NA12877, .../NA12877_chr20-22.R1.1.fq, .../NA12877_chr20-22.R2.1.fq]
        [NA12877, .../NA12877_chr20-22.R1.2.fq, .../NA12877_chr20-22.R2.2.fq]
        [NA12877, .../NA12877_chr20-22.R1.3.fq, .../NA12877_chr20-22.R2.3.fq]
        [NA12877, .../NA12877.g.vcf.gz, .../NA12877.g.vcf.gz.tbi]
        [NA12877, .../NA12877.g.vcf.gz, .../NA12877.g.vcf.gz.tbi]
        [NA12877, .../NA12877.g.vcf.gz, .../NA12877.g.vcf.gz.tbi]
        ERROR ~ Error executing process > 'JOINT_GENOTYPE (1)'

        Caused by:
          Process `JOINT_GENOTYPE` input file name collision -- There are multiple input files for each of the following file names: NA12877.g.vcf.gz, NA12877.g.vcf.gz.tbi


        Container:
          /scratch/pawsey1227/fjaya/nextflow-on-hpc-materials/singularity/quay.io-biocontainers-gatk4-4.6.2.0--py310hdfd78af_1.img

        Tip: view the complete command output by changing to the process work dir and entering the command `cat .command.out`

         -- Check '.nextflow.log' file for details

        ```

        - The first three lines are the outputs of our `.splitFastq()` operation, this has not changed since the last time the workflow was run
        - The last three lines are the outputs emitted from the `GENOTYPE` process. One output reflects the run for one of the chunks processed. However, all the **output names are the same**. This is the cause of the error.

There are several ways to resolve this, such as adding the chunk name to the outputs, similar to how `.splitFastq()` manages this. In the next section, we will implement the "gather" pattern in our pipeline to fix this.

### Gather: combining our alignments

Now that we have sucessfully split our reads, we will implement a gather pattern to bring our alignments into a single file again.

!!! tip "Different patterns for different needs"

    There is no one-size-fits-all approach for scattering and gathering. How this is implemented in Nextflow will be highly dependent on your workflow structure, and input and output files. [Nextflow patterns](https://nextflow-io.github.io/patterns/) provides examples of commonly used patterns that support a range of different needs, such as splitting text and CSV files, to re-grouping and organising results for downstream processing.

    Whole‑genome alignment is an excellent use case for scatter–gather as the alignment of each chunk can run independently, dramatically reducing walltime for this heavy step. Once alignment is complete, the BAM files are merged and the workflow proceeds as normal.

    Where you choose to re‑gather your data will depend on where your bottlenecks are and at which points you need to process the dataset as a whole again.

!!! example "Exercise"

    Import the `MERGE_BAMS` module in your `main.nf` file.

    ```groovy title="main.nf" hl_lines="4"
    include { FASTQC } from './modules/fastqc'
    include { SPLIT_FASTQ } from './modules/split_fastq'
    include { ALIGN } from './modules/align'
    include { MERGE_BAMS } from './modules/merge_bams'
    include { GENOTYPE } from './modules/genotype'
    include { JOINT_GENOTYPE } from './modules/joint_genotype'
    include { STATS } from './modules/stats'
    include { MULTIQC } from './modules/multiqc'
    ```

    Update the workflow script:

    ```groovy title ="main.nf"
        // Split FASTQs for each sample
    split_fqs = reads
        .splitFastq(limit: 3, pe: true, file: true)
        .map { sample, r1, r2 ->
            def chunk_id = r1.toString().tokenize('.')[2]
            return [ sample, r1, r2, chunk_id ]
        }
    
    ALIGN(split_fqs, bwa_index)

    gathered_bams = ALIGN.out.aligned_bam
        .groupTuple()

    MERGE_BAMS(gathered_bams)

    // Run genotyping with aligned bam and genome reference
    GENOTYPE(MERGE_BAMS.out.aligned_bam, ref)
    ```

    And update the `ALIGN` module to take the `chunk_id`:

    ```
    process ALIGN {

        container "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:1bd8542a8a0b42e0981337910954371d0230828e-0"
        publishDir "${params.outdir}/alignment"

        input:
        tuple val(sample_id), path(reads_1), path(reads_2), val(chunk_id)
        tuple val(ref_name), path(bwa_index)

        output:
        tuple val(sample_id), path("${sample_id}.${chunk_id}.bam"), path("${sample_id}.${chunk_id}.bam.bai"), emit: aligned_bam

        script:
        """
        bwa mem -t $task.cpus -R "@RG\\tID:${sample_id}\\tPL:ILLUMINA\\tPU:${sample_id}\\tSM:${sample_id}\\tLB:${sample_id}\\tCN:SEQ_CENTRE" ${bwa_index}/${ref_name} $reads_1 $reads_2 | samtools sort -O bam -o ${sample_id}.${chunk_id}.bam
        samtools index ${sample_id}.${chunk_id}.bam
        """

    }
    ```

    Re-run the pipeline
    ```
    ./run.sh
    ```

Now, let's re-inspect that the merge worked as intended.

!!! example "Exercise"

    ```bash
    tree -a
    ```
    ```console title="Output"
    .
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .command.trace
    ├── .exitcode
    ├── NA12877.1.bam -> /scratch/pawsey1227/fjaya/nextflow-on-hpc-materials/part2/work/fc/bebb980d3d81cba7dacb6d052faf08/NA12877.1.bam
    ├── NA12877.1.bam.bai -> /scratch/pawsey1227/fjaya/nextflow-on-hpc-materials/part2/work/fc/bebb980d3d81cba7dacb6d052faf08/NA12877.1.bam.bai
    ├── NA12877.2.bam -> /scratch/pawsey1227/fjaya/nextflow-on-hpc-materials/part2/work/ec/b16bfd4ccc0c627553eb0e55337b21/NA12877.2.bam
    ├── NA12877.2.bam.bai -> /scratch/pawsey1227/fjaya/nextflow-on-hpc-materials/part2/work/ec/b16bfd4ccc0c627553eb0e55337b21/NA12877.2.bam.bai
    ├── NA12877.3.bam -> /scratch/pawsey1227/fjaya/nextflow-on-hpc-materials/part2/work/c7/b93c40ec2e8c9eebe4cce19ac96ef7/NA12877.3.bam
    ├── NA12877.3.bam.bai -> /scratch/pawsey1227/fjaya/nextflow-on-hpc-materials/part2/work/c7/b93c40ec2e8c9eebe4cce19ac96ef7/NA12877.3.bam.bai
    ├── NA12877.bam
    └── NA12877.bam.bai
    ```

    ```bash
    cat .command.sh
    ```
    ```console title="Output"
    #!/bin/bash -ue
    samtools cat NA12877.3.bam NA12877.2.bam NA12877.1.bam | samtools sort -O bam -o NA12877.bam
    samtools index NA12877.bam
    ```

    We can now see that each bam was merged into a single bam file for the sample, and the rest of the workflow progressed as normal.

To sum up this step, you successfully:

    1. Split the paried FASTQ reads into 3 chunks using `.splitFastq`
    2. Aligned each chunk in parallel
    3. Merged the aligned chunks into a single BAM (`MERGE_BAMS()`)
    4. Ran the remainder of the workflow as usual

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