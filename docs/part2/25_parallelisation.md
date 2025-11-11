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

![](figs/00_benchmark_at_scale_theme.png)

- Trade-offs between time x SUs x efficiency
- recall we do not always want to provide more threads/split

## Multithreading

Benchmarking results from Part 1 suggest that 4 threads is the most optimal
setting for `bwa mem`, with diminishing returns observed at 6 and 8 threads.

The `script` block for the `ALIGN` process is already configured to
dynamically use the number of threads specified via `task.cpus`.

Let's update our configuration so that the `ALIGN` process requests 4 CPUs
and the memory proportional the CPU required.

!!! question "Discussion"

    Given we want the `ALIGN` process to consistently use 4 CPUs, what approaches
    to configuration would you take on the system you are using?

    Things to consider include:

    - Which `.config` file would you want to use? (Consider whether this is
    something that needs to be portable across systems vs. system specific)
    - How much memory would you provide? (Consider the effective RAM/CPUs
    proportion of the queue or partition).
    - Based on the CPU and memory requirements, which directive would be
    more suitable to use - `withLabel` or `withName`? (Consider whether
    it matches an existing configuration)

    === "Gadi (PBS)"

        ??? note "Answers"

            - `custom.config` to ensure it is tuned for the `normalbw` queue.
            - 4 CPUs with 9 GB memory
            - This is the same configuration as `FASTQC`, therefore `withLabel`
            will be used
            - TODO: confirm 4 CPUs with 16 GB vs 4CPU9GB

    === "Setonix (Slurm)"

        ??? note "Answers"

            - `custom.config` to ensure it fits the `work` partition.
            - 4 CPUs with 7 GB memory
            - No other configuration are similar, therefore a new configuration
            is needed using `withName`.

!!! example "Exercise"

    In `conf/custom.config`, add the following inside the `process` scope:

    === "Gadi (PBS)"

        ```groovy title="conf/custom.config" hl_lines="3-8"
        process {
        // truncated

        	withName: 'ALIGN' {
        		cpus = 4 // Optimal threads
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
        		cpus = 4 // Optimal threads
        		memory = 7.GB // 4 CPUs * 1.8 GB
        		time = 2.min
        	}
        }
        ```

    Save your file and run with `./run.sh`.

## Multi-processing with scatter-gather

![](figs/00_Scatter_gather_fig.png)

We will parallelise the alignment step. Scatter-gathering in this step is widely
used as whole genome data is large, time-consuming, and the reads can be mapped
to a reference independently of each other.

This requires including the pre-defined modules:

- `SPLIT_FASTQ`
- `ALIGN_CHUNK`
- `MERGE_BAMS`

!!! example "Exercise"

    Replace the module imports:

    _Before:_

    ```groovy title="main.nf" hl_lines="2"
    include { FASTQC } from './modules/fastqc'
    include { ALIGN } from './modules/align'
    include { GENOTYPE } from './modules/genotype'
    include { JOINT_GENOTYPE } from './modules/joint_genotype'
    include { STATS } from './modules/stats'
    include { MULTIQC } from './modules/multiqc'
    ```

    _After:_

    ```groovy title="main.nf" hl_lines="2-4"
    include { FASTQC } from './modules/fastqc'
    include { SPLIT_FASTQ } from './modules/split_fastq'
    include { ALIGN_CHUNK } from './modules/align_chunk'
    include { MERGE_BAMS } from './modules/merge_bams'
    include { GENOTYPE } from './modules/genotype'
    include { JOINT_GENOTYPE } from './modules/joint_genotype'
    include { STATS } from './modules/stats'
    ```

Next, we need to update our workflow scope to facilitate the scatter-gather.

TODO: Note some key parts. e.g. getting chunk_id, .groupTuple().map{}.MERGE_BAMS()
to merge again correctly by sample.

!!! example "Exercise"

    _Before:_

    ```groovy title="main.nf" hl_lines="7-11"
    workflow {
    // truncated

        // Run the fastqc step with the reads_in channel
        FASTQC(reads)

        // Run the align step with the reads_in channel and the genome reference
        ALIGN(reads, bwa_index)

        // Run genotyping with aligned bam and genome reference
        GENOTYPE(ALIGN.out.aligned_bam, ref)

        // Gather gvcfs and run joint genotyping
        all_gvcfs = GENOTYPE.out.gvcf
            .map { _sample_id, gvcf, gvcf_idx -> [ params.cohort_name, gvcf, gvcf_idx ] }
            .groupTuple()
        JOINT_GENOTYPE(all_gvcfs, ref)

    // truncated
    }
    ```

    _After:_

    ```groovy title="main.nf" hl_lines="7-42"
    workflow {
    // truncated

        // Run the fastqc step with the reads_in channel
        FASTQC(reads)

        // Perform a scatter/gather workflow by splitting the FASTQs into multiple chunks,
        // aligning each chunk separately, then merging the BAMs together at the end

        // Split FASTQs for each sample
        reads_to_split = reads
            .map { it + [ params.split_n ] }
        SPLIT_FASTQ(reads_to_split)

        // Extract the chunk ID from the split FASTQs and run ALIGN_CHUNK
        // TODO: simplify
        split_fqs_r1 = SPLIT_FASTQ.out.split_fq
            .map { sample_id, fqs_1, _fqs_2 -> [ sample_id, fqs_1 ] }
            .transpose()
            .map { sample_id, fq1 -> {
                def chunk_id = fq1.baseName.tokenize(".")[0]
                [ sample_id, chunk_id, fq1 ]
            } }
        split_fqs_r2 = SPLIT_FASTQ.out.split_fq
            .map { sample_id, _fqs_1, fqs_2 -> [ sample_id, fqs_2 ] }
            .transpose()
            .map { sample_id, fq2 -> {
                def chunk_id = fq2.baseName.tokenize(".")[0]
                [ sample_id, chunk_id, fq2 ]
            } }
        split_fqs = split_fqs_r1.join(split_fqs_r2, by: [0, 1])

        ALIGN_CHUNK(split_fqs, bwa_index)

        // Gather all BAM chunks for a sample and run MERGE_BAMS
        gathered_bams = ALIGN_CHUNK.out.aligned_bam
            .groupTuple()
            .map { sample_id, chunk_id, bams, bais -> [ sample_id, bams, bais ] }
        MERGE_BAMS(gathered_bams)

        // Run genotyping with aligned bam and genome reference
        GENOTYPE(MERGE_BAMS.out.aligned_bam, ref)

        // Gather gvcfs and run joint genotyping
        all_gvcfs = GENOTYPE.out.gvcf
            .map { _sample_id, gvcf, gvcf_idx -> [ params.cohort_name, gvcf, gvcf_idx ] }
            .groupTuple()
        JOINT_GENOTYPE(all_gvcfs, ref)

    // truncated
    }
    ```

 <!--- If you run the workflow now, it will fail because params.split_n is
 not defined. Could be a troubleshooting exercise --->

Lastly, a new parameter, `params.split_n` was defined. This needs to be
added.

!!! example "Exercise"

    - In `nextflow.config` within the `params` scope, define `split_n`.
    - Provide a default value of `split_n = 3`. This will split the FASTQ file
    into three chunks.

    ??? note "Answer"

        ```groovy title="nextflow.config" hl_lines="12"
        // Define params
        params {
            samplesheet = "$projectDir/samplesheet_single.csv"
            ref_prefix = "$projectDir/../data/ref/Hg38.subsetchr20-22"
            ref_fasta = "${params.ref_prefix}.fasta"
            ref_fai = "${params.ref_prefix}.fasta.fai"
            ref_dict = "${params.ref_prefix}.dict"
            bwa_index = "$projectDir/../data/ref"
            bwa_index_name = "Hg38.subsetchr20-22.fasta"
            cohort_name = "cohort"
            outdir = "results"
            split_n = 3
        }
        ```

Lastly, add the `process_small` labels to each of the modules:

- `SPLIT_FASTQ`
- `ALIGN_CHUNK`
- `MERGE_BAMS`

!!! example "Exercise"

    For `SPLIT_FASTQ`:

    ```groovy title="modules/split_fastq.nf"
    process SPLIT_FASTQ {

        tag "split fastqs for ${sample_id}"
        container "quay.io/biocontainers/fastp:1.0.1--heae3180_0"
        label "process_small"

        input:
        tuple val(sample_id), path(reads_1), path(reads_2), val(n)

        output:
        tuple val(sample_id), path("*.${sample_id}.R1.fq"), path("*.${sample_id}.R2.fq"), emit: split_fq

        script:
        """
        fastp -Q -L -A -i $reads_1 -I $reads_2 -o ${sample_id}.R1.fq -O ${sample_id}.R2.fq -s $n
        """
    }
    ```
    For `MERGE_BAMS`:

    ```groovy title="modules/split_fastq.nf"
    process MERGE_BAMS {

        container "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:1bd8542a8a0b42e0981337910954371d0230828e-0"
        publishDir "${params.outdir}/alignment"
        label "process_small"

        input:
        tuple val(sample_id), path(bams), path(bais)

        output:
        tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: aligned_bam

        script:
        """
        samtools cat ${bams} | samtools sort -O bam -o ${sample_id}.bam
        samtools index ${sample_id}.bam
        """

    }
    ```
    For `ALIGN_CHUNK`:

    ```groovy title="modules/align_chunk.nf.nf"
    process ALIGN_CHUNK {

        container "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:1bd8542a8a0b42e0981337910954371d0230828e-0"
        label "process_small"

        input:
        tuple val(sample_id), val(chunk), path(reads_1), path(reads_2)
        tuple val(ref_name), path(bwa_index)

        output:
        tuple val(sample_id), val(chunk), path("${sample_id}.${chunk}.bam"), path("${sample_id}.${chunk}.bam.bai"), emit: aligned_bam

        script:
        """
        bwa mem -t $task.cpus -R "@RG\\tID:${sample_id}\\tPL:ILLUMINA\\tPU:${sample_id}\\tSM:${sample_id}\\tLB:${sample_id}\\tCN:SEQ_CENTRE" ${bwa_index}/${ref_name} $reads_1 $reads_2 | samtools sort -O bam -o ${sample_id}.${chunk}.bam
        samtools index ${sample_id}.${chunk}.bam
        """

    }
    ```
    now run the parallelised pipeline
    ```bash
    ./run.sh
    ```

    Your output should look similar to this
    ```groovy
    [82/d867f8] FASTQC (fastqc on NA12877)             [100%] 1 of 1 ✔
    [43/40bd7e] SPLIT_FASTQ (split fastqs for NA12877) [100%] 1 of 1 ✔
    [0a/0ef899] ALIGN_CHUNK (2)                        [100%] 3 of 3 ✔
    [c1/79d4cf] MERGE_BAMS (1)                         [100%] 1 of 1 ✔
    [f8/c1663c] GENOTYPE (1)                           [100%] 1 of 1 ✔
    [a7/a67522] JOINT_GENOTYPE (1)                     [100%] 1 of 1 ✔
    [49/418d69] STATS (1)                              [100%] 1 of 1 ✔
    [2e/1b5029] MULTIQC                                [100%] 1 of 1 ✔
    Completed at: 11-Nov-2025 11:59:08
    Duration    : 4m 31s
    CPU hours   : (a few seconds)
    Succeeded   : 10
    ```

    NOTE: The alignment module was updated to use a scatter-gather approach. Instead of aligning the entire FASTQ in one go with the ALIGN module, the workflow now:

    1. Splits the FASTQ into 3 chunks (SPLIT_FASTQ).
    2. Aligns each chunk in parallel (ALIGN_CHUNK, see the 3 of 3 completed).
    3. Merges the aligned chunks into a single BAM (MERGE_BAMS).

    This change optimises performance for large datasets by leveraging parallel processing.

## A note on dynamic resourcing

Since our data is small, this will not make an impact. These are nice strategies
when you need to process data that require unequal CPUs, memory, or walltime.
For example, if we were to process all the chromosomes which considerably vary
in size.

Two approaches:

| Directive | Closure Example             | Attempt 1 (Initial Run) | Attempt 2 (First Retry) |
| --------- | --------------------------- | ----------------------- | ----------------------- |
| `memory`  | `{ 2.GB * task.attempt }`   | 2 GB                    | 4 GB                    |
| `time`    | `{ 1.hour * task.attempt }` | 1 hour                  | 2 hours                 |

```groovy
process {
  withName: 'ALIGN' {
    memory = { reads.size() * 2}
  }
}
```

The same concepts of configuring resources will apply here, aim to fit the
appropriate queues/partitions.

## Summary

Recall that we do not always want to parallelise everything. There are reasons
to avoid this due to the data set or analysis requiring data be analysed as a
whole, or computational reasons where processes become less efficient.

These are some strategies to consider for your own pipelines. Run benchmarks,
identify long-running or inefficent processes and consider which one of these
approaches are supported by the tools (multi-threading), or can be split and
processed in parallel.
