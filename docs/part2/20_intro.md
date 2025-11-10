# Intro to part 2

!!! info "Learning objectives"

    - Recognise how Nextflow can be run across different systems with shared
    and system-specific configuration
    - Compare the trade-offs between running pipelines with pbs/sbatch scripts
    vs. Nextflow, such as ease-of-use, flexibility, and performance
    - Understand why we need custom pipelines and configuration

## Overview

Introduce our workflow use case and lesson structure:

- Getting "unoptimised" workflow running
- Benchmarking, testing and optimisation on a single sample.
- Review HPC optimisation strategies 
- Implement them
- Scale up to multiple samples

TODO: Figures comparing serially run bash scripts vs. unoptimised vs.
optimised workflow structure. Figure to include walltime comparison of each
approach, perhaps SU usage/cost and resource usage too.

- Outline what are we optimising for? Time/cost/throughput? Focusing on time an
throughput:

Throughput may be important for clinicians who need a fast turnaround
of results by expending more SUs.

What part 2 is NOT:

- How to write and build nextflow pipelines and logic
- Understanding the contents input and output files in detail

## Why do you need custom pipelines?

Reasons for using custom pipeline:

- Fit for your own purposes
- Existing pipelines (nf-core) may not do everything required, pipelines
doesn't exist
- nf-core can be expensive and misconfigured
- It doesn't exist!

## The pipeline file anatomy

TODO: `tree` of part2 structure - showcasing locations of /conf, /modules,
main.nf, nextflow.config

### `main.nf` and `modules/`

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

Set up as modules - a recurring theme is that there are Nextflow files that can
be left as is, and some that need tweaking. This makes your pipelines portable,
organised, reproducible, and easy to set up across different systems.

Here, we will not edit the modules, but configure how and where to run them on
HPCs.

- Mention modules imported.
- Why modules?
- During scatter-gather later, easy to swap out modules, rather than coding
in processes in main.nf. Reinforce reproducibility, not all files need to be
touched
- https://sydney-informatics-hub.github.io/template-nf-guide/notebooks/modules.html

### `nextflow.config` and `conf/`

![](figs/schema.png)

We start with basic preconfigured - executors, the default queues, and
singularity enabled.

As well as the default parameters for the workflow to run.

[Config profiles](https://www.nextflow.io/docs/latest/config.html#config-profiles)
as one of NXFs powerful features that make porting code easy
while maintaining reproducibility.

We have defined one profile for each slurm and pbspro - this makes
running this workshop across two HPCs possible!!

```groovy title="nextflow.config"
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
}

// Define HPC profiles to run with job scheduler
profiles {
    // Use this profile to interact with the scheduler on setonix
    slurm { includeConfig "conf/slurm.config" }

    // Use this profile to interact with the scheduler on gadi
    pbspro { includeConfig "conf/pbspro.config" }
}
```

=== "Gadi (PBS)"

    ```groovy

    ```

=== "Pawsey (Slurm"

TODO: snippet of starting `nextflow.config` and `conf/pbs-or-slurm.config`
- This is what we will be starting with, building up towards the optimised
pipeline by end of part 2
- Show for pbspro and slurm

