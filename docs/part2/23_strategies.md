# Optimising Nextflow for HPC

!!! info "Learning objectives"

    - Identify the key factors that impact workflow performance on HPC systems
    - Describe common workflow optimisation strategies in the context of Nextflow

## Overview

Our workflow is now functional, it runs successfully with scheduler-specific settings and outputs useful trace files showing resource usage. Recall the principles covered in Lesson 1.4, we will apply them to our custom workflow.

From here, we will explore ways of optimising our workflow for more efficient execution on the HPC. 

!!! note "Why do we care about optimisation?"

    Bioinformatics data processing on HPCs has a significant carbon footprint due to the energy required to run our computations. We optimise our workflows to not only reduce their runtime but also adopt more sustaintable computing practices. 

    *[This short editorial](https://www.nature.com/articles/s43588-023-00506-2) from Nature Computational Science highlights the challenges research computing faces in the context of the climate crisis*

Workflow optimisation involves fine-tuning your pipeline to make it more efficient. It can be used to reduce runtime, idling hardware, and cost (e.g. when on a system that charges use based on service units). 

## Why optimise?

A small workflow run a handful of times might not benefit dramatically from optimisation. Many Nextflow workflows that employ good practices (e.g. nf-core) will run with default configuration, but defaults might not always fit your data and therefore the behaviour of your processes, or the constraints of your cluster. Think back to Part 1 and the configuration customisations we implemented for our nf-core workflow. 

Optimising workflows on HPC becomes especially important when:

- You are analysing large datasets or many samples 
- You will execute your pipeline repeatedly 
- You are operating on a system that uses a service unit or time-limited allocation model
- Your processes have data-dependent resource usage 

By optimising, you are making your workflow resilient and scalable. 

## What affects performance?

Efficiency of any workflow on HPC dependens on the interaction of three factors: 

### 1. Your HPC system 

We have already witnessed many differences between Gadi and Setonix in previous lessons. A workflow that performs well on one cluster may perform poorly on another simply because the underlying architecture and scheduler rules differ. 

Good optimisation respects the boundaries of the system you're working on. When planning an optimisation approach, consider:

| HPC characteristic | What it means | Why it matters for optimisation                     |
| --------- | ---------------------- | ------------- |
| **Default scheduler behaviour**         | Policies set by administrators: fair-share, job priorities, backfill rules, default limits  | Affects queue wait time, job placement efficiency, and how many tasks can run in parallel                            |
| **Queue limits**                        | Maximum walltime, CPU count, and memory allowed per queue or partition                      | Determines which queues you can use, how large each job can be, and whether your workflow gets delayed               |
| **Node architecture**                   | Hardware layout: cores per node, memory per node, CPU type (Intel/AMD), GPUs, local scratch | Ensures you request resources that “fit” the node, avoid resource fragmentation, and maximise throughput             |
| **Charging model** | How HPC usage is accounted (CPU proportion, memory proportion, or the maximum of both)      | Guides you to request only what you need over requesting directly increases SU consumption without improving runtime |

### 2. The characteristics of your data 

Data shapes the computational behaviour of bioinformatics workflows. Even two workflows with identical code can perform very differently depending on the file sizes, sample numbers, and data complexity. Understanding these factors can help you anticipate bottlenecks and assign resources more accurately. When planning an optimisation approach, consider: 

| Data characteristic     | What it means         | Why it matters for optimisation      |
| ------------------------------------ | -------------------------- | -------------------------------- |
| **File size**                        | The total size of FASTQ, BAM/CRAM, reference genomes, annotation files, or intermediate outputs                      | Larger files increase memory requirements, disk I/O, runtime, and queue time; they also influence whether single-threading or multithreading is more efficient            |
| **Sample number**                    | Total number of samples in the analysis, including replicates or cohorts                                             | More samples → more processes → heavier scheduler load; the workflow may require scatter–gather to parallelise effectively and avoid bottlenecks                          |
| **Data heterogeneity**               | Variability in file sizes, read depth, sample complexity, or quality across inputs                                   | Highly variable samples produce uneven resource usage; some processes may require per-sample resource overrides to prevent memory kills or slowdowns                      |
| **Data type**                        | Whether data are short reads, long reads, single-cell, imaging derivatives, matrices, VCFs, etc.                     | Different data modalities have different computational profiles (I/O-heavy, CPU-heavy, memory-heavy); optimisation strategies should account for the modality’s behaviour |
| **I/O intensity**                    | Frequency and volume of read/write operations (large temporary files, sort steps, indexing, BAM ↔ FASTQ conversions) | I/O-heavy processes benefit from local SSD or node-local scratch; misconfigured I/O can add hours to runtime on shared filesystems                                        |
| **Parallelisability**                | Whether samples or chunks of data can be processed independently                                                     | Determines when scatter–gather is useful, how many jobs can run concurrently, and how well the workflow scales on HPC                                                     |
| **Compression and indexing formats** | gzip vs bgzip, BAM vs CRAM, presence of .bai/.crai/.fai, CCS vs raw reads                                            | Impacts CPU time, memory, and I/O behaviour; inefficient formats slow down the entire workflow                                                                            |


### 3. The structure of your workflow 

Even with the same tools and data, two workflows can behave differently depending on their structure: 

* Number of processes 
* Order of dependencies 
* Opportunities for parallelism
* Whether steps are CPU-bound, memory-bound, or I/O boind 
* Incorporated tool's ability to multithread 

## What will we optimise?

For the remainder of Part 2 we will apply the strategies introduced from Part 1 to optimise our custom workflow. In particular, we will:

1. Assign appropriate resources for each per process
    - Use trace files to fine-tune `cpus`, `memory`, and `time`
    - Tune these to fit with the architecutre of your HPC

2. Enable multi-threading for BWA MEM
    - Increase speed of alignment by using multiple threads

3. Parallelise alignment using scatter-gather
    - Split our FASTQ files, align "simultaneously", and combine again

These same principles can be applied to any HPC environment, but in this workshop, we’ll focus on NCI Gadi and Pawsey Setonix as practical examples.