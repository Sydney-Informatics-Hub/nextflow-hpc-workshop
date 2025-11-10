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

bwa mem benchmarking results from part 1 suggests 2N threads is ideal

!!! example "Exercise"

    TODO Update bwa mem script with ${task.cpus}

## Multi-processing with scatter-gather

![](figs/00_Scatter_gather_fig.png)

!!! example "Exercise"

    _Before:_

    ```groovy title="main.nf"
    include { FASTQC } from './modules/fastqc'
    include { ALIGN } from './modules/align'
    include { GENOTYPE } from './modules/genotype'
    include { JOINT_GENOTYPE } from './modules/joint_genotype'
    include { STATS } from './modules/stats'
    include { MULTIQC } from './modules/multiqc'

    // Define the workflow
    workflow {

        // Define the fastqc input channel
        reads = Channel.fromPath(params.samplesheet)
            .splitCsv(header: true)
            .map { row -> {
                // def strandedness = row.strandedness ? row.strandedness : 'auto'
                [ row.sample, file(row.fastq_1), file(row.fastq_2) ]
            }}

        bwa_index = Channel.fromPath(params.bwa_index)
            .map { idx -> [ params.bwa_index_name, idx ] }
        ref = Channel.of( [ file(params.ref_fasta), file(params.ref_fai), file(params.ref_dict) ] )

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

        // Get VCF stats
        STATS(JOINT_GENOTYPE.out.vcf)

        // Collect summary data for MultiQC
        multiqc_in = FASTQC.out.qc_out
            .mix(STATS.out.stats_out)
            .collect()

        /*
        * Generate the analysis report with the
        * outputs from fastqc and bcftools stats
        */
        MULTIQC(multiqc_in)

    }
    ```

    _After:_
    ```groovy title="main.nf"

    ```

!!! example "Exercise"

    TODO Update the channels in main.nf

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
