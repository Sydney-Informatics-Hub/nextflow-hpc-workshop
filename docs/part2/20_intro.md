# 2.0 Introduction 

!!! info "Learning objectives"

    - Recall the overall structure of a modular Nextflow pipeline, including the roles of `main.nf`, `modules/`, and configuration files
    - Recognise the importance of separating workflow logic from system-specific configuration for portability and reproducibility
    - Identify which parts of the pipeline need to be adapted when moving between systems (e.g. from local to HPC)

## 2.0.1 Log back in to your assigned HPC

Log in to your assigned HPC with the user account and password provided to you on day 1:

=== "Gadi"

    ```bash
    ssh username@gadi.nci.org.au
    ```

=== "Setonix"

    ```bash
    ssh username@setonix.pawsey.org.au
    ```

!!! note

    Be sure substitute your assigned user name for `username` in the above code example.

Navigate to the scratch space for the workshop project, then open your cloned part2 repository:

=== "Gadi"

    ```bash
    cd /scratch/vp91/$USER/nextflow-on-hpc-materials/part2
    ```

=== "Setonix"

    ```bash
    cd /scratch/courses01/$USER/nextflow-on-hpc-materials/part2
    ```

## 2.0.2 Configuring a custom pipeline

Part 2 of this workshop builds on the foundational HPC and Nextflow configuration concepts introduced in Part 1. We will now apply these concepts to configure a variant calling pipeline for efficient execution on HPC systems.

To keep the focus on configuration, the pipeline code and logic are provided for you and will not be reviewing the contents of input and output files in detail, beyond configuration needs. We’ll begin by getting the pipeline running on the HPC, then progressively explore how to benchmark performance, understand HPC-specific constraints, and implement optimisations to improve efficiency and resource use.

Throughout this section, we’ll continue using the variant calling example to deepen your understanding of the key decisions involved in tuning pipelines for the specific HPC infrastructure you work on.

!!! warning "We will be using small data sets!"

    In this section, we will continue using the small paired-end reads from Part 1 to demonstrate key concepts related to running and configuring workflows on HPC systems.
    
    Because the data is deliberately small, the processes will run quickly and with minimal resource requirements. This helps illustrate how the workflow behaves, but it also means that:
    
    - Some resource-related errors (e.g. out-of-memory, walltime exceeded) **will not occur**, even if your configuration is suboptimal.
    - Performance trade-offs (e.g. with different `cpus` or `memory` settings) may be **less noticeable** than with real-world data.
    - HPC scheduling behaviour may differ slightly - smaller jobs are often scheduled more quickly and occupy less system space.
    
    When running real data sets on HPC systems, you may encounter different behaviours, longer runtimes, and additional errors. 

## 2.0.3 Why do you need custom pipelines?

There are several reasons why you might need to develop or adapt your own
Nextflow pipeline:

- **Tailored to your specific needs** - Custom pipelines give you full control over input/output formats, tool parameters, workflow logic, and configuration options.
- **Gaps in available pipelienes** - Existing pipelines (e.g. nf-core) may not cover your use case, or a relevant pipeline may not exist at all.
- **Resource optimisation** - nf-core pipelines are generalised by design and may be over-provisioned or misconfigured for your HPC environment. Although easier to get running out-of-the box, this could lead to inefficient use of HPC resources or being charged excess service units (SUs)!

## 2.0.3 The scenario: variant calling on HPC

TODO: Describe starting variant calling pipeline briefly

![](figs/00_workflow_illustration.png)

TODO: scenario for benchmarking and optimising

## 2.0.4 The pipeline file anatomy

This pipeline builds on the structure introduced in our introductory
[Nextflow for the life sciences](https://sydney-informatics-hub.github.io/hello-nextflow-2025/part2/00_intro/#205-nextflowing-the-workflow) workshop, where the workflow is separated into the 
**data processing logic**, and **system-specific configuration**. This layout helps keep things reproducible, easy
to maintain, and simple to adapt across different environments - like moving from your laptop to an HPC!

At a glance:

```bash
tree
```
```console
# (some folders are truncated)
.
├── conf                    # System-specific configuration files
│   ├── pbspro.config
│   └── slurm.config
├── main.nf                 # Main workflow structure
├── modules                 # Individual process modules
│   ├── align_chunk.nf
│   ├── align.nf
│   ├── fastqc.nf
│   ├── genotype.nf
│   └── ...
├── nextflow.config         # Base parameters and profile defs
├── samplesheet_full.csv
└── samplesheet_single.csv
```

### 2.0.4.1 `main.nf` and `modules/`: What to run?

The core workflow logic lives in `main.nf` and the `modules/` directory. This is where you define what to do, like running `FASTQC`, `ALIGN`, or `GENOTYPE`. Our pipelines are structured using modules, allowing you to **separate workflow components from system-specific configuration**. They’re written once and can be used across different systems without needing to change the underlying code.

In this part, we **don’t edit the module files directly**. Instead, we focus on configuring _how and where_ these steps run on HPCs, using separate config files. This separation makes it easy to test a pipeline locally and later scale it up on a system like Gadi or Setonix, without rewriting processes.

[TODO]: add figure explaining below concept

!!! tip "Why modules?"

    Nextflow pipelines are often set up using modules, which help keep things clean and organised.
    The key idea is that some files can stay the same no matter where you run them, 
    while others are easily tweaked to match the system you're using. 
    This is especially useful when you're testing a pipeline locally first, 
    then moving to run it on an HPC - something many researchers do. 
    
    If everything were hardcoded, you'd end up changing lots of lines across multiple 
    files, which can quickly get messy and error-prone. By using modules, you only
    need to configure how and where things run, without rewriting the pipeline a
    itself. It also makes it easier to reuse parts of a pipeline or swap out tools 
    later on, without starting from scratch.

    For more information, see ["What's in `modules/`"](https://sydney-informatics-hub.github.io/template-nf-guide/notebooks/modules.html).
    
```bash
cat main.nf
```
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

This structure also makes it easier to swap in alternative tools and processes, especially later when working with scatter-gather patterns - all without cluttering `main.nf` or compromising reproducibility.

### 2.0.4.2 `nextflow.config` and `conf/`: How and where to run it

![](figs/00_workflow_illustration.png)

Nextflow’s configuration files define how and where each process runs - including what resources to request, which job scheduler to use, and whether to run using containers.

In the context of HPCs, this means specifying:

- How many CPUs, memory, and time a process should use
- The appropriate **executor** (e.g. PBS on Gadi or SLURM on Setonix)
- The default **queue/partition** and optional account/project codes
- Whether and how to use Singularity containers
- Plus many other useful features

These settings are defined in the main `nextflow.config`, and extended using config profiles - one for each target system. This separation allows you to **run the same pipeline across different HPCs just by switching profiles, without modifying the core workflow.**

Here's what the `nextflow.config` file looks like:

```bash
cat nextflow.config
```
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

It contains the default parameters required for `main.nf`, and the profiles for Gadi (PBS) and Setonix (Slurm).

Each profile brings in a system-specific config file from the `conf/` folder:

=== "Gadi (PBS)"

    ```bash
    cat conf/pbspro.config
    ```
    ```groovy title="conf/pbspro.config"
    params.pbspro_account = ""

    process {
      executor = 'pbspro'
      queue = 'normalbw'
      clusterOptions = "-P ${params.pbspro_account}"
      module = 'singularity'
    }

    singularity {
      enabled = true
      autoMounts = true
      cacheDir = "${projectDir}/singularity"
    }
    ```

=== "Setonix (Slurm)"

    ```bash
    cat conf/slurm.config
    ```
    ```groovy title="conf/slurm.config"
    params.slurm_account = ""

    process {
      executor = 'slurm'
      queue = 'work'
      clusterOptions = "--account=${params.slurm_account}"
      module = 'singularity/4.1.0-slurm'
    }

    singularity {
      enabled = true
      autoMounts = true
      cacheDir = "${projectDir}/singularity"
    }
    ```

This setup makes it easy to test and run the same workflow in different environments. By the end of Part 2, we’ll have extended these configs to better reflect the characteristics of each system - improving efficiency without touching the pipeline logic itself.

## 2.0.5 Summary

This section introduced the basic structure of the custom variant calling Nextflow pipeline for the remainder of Part 2, emphasising the separation between workflow logic (`main.nf`, `modules/`) and system-specific configuration (`nextflow.config`, `conf/`). We reviewed how this separation supports portability, reproducibility, and ease of adaptation across environments, such as when transitioning from local testing to running on HPC systems like Gadi (PBS) and Setonix (Slurm).