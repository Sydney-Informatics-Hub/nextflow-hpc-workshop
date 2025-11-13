# Optimising for HPC

!!! info "Learning objectives"

    - Identify the key factors that impact workflow performance on HPC systems
    - Describe common workflow optimisation strategies

## Overview

Our workflow is now functional - it runs sucessfully with scheduler-specific settings and outputs useful trace files showing resource usage.

Now, we shift focus from making it run, to getting to run efficiently.

TODO: figure benchmarking progress e.g. it runs, it runs well, it scales

Workflow optimisation involves fine-tuning your pipeline to reduce runtime, resource waste, and cost (e.g. service units). This is particularly important on share HPC infrastructure where jobs are charged based on compute and memory usage.

## Why optimise?

Small inefficiencies for small, or infreuqently-run workflows may not make a large impact. Many users run nf-core pipelines on HPCs successfully, using the default configuration. Optimising workflows on HPC becomes increasingly valuable when:

- You are running many samples or large data sets (high-throughput)
- Your workflow will be run repeatedly
- You are approaching allocation limits or want it to be more efficient 

TODO: figures would be useful to demo these

## What affects performance?

Awareness of these factors can help you make decisions about the optimisation approach to use.

- The specific HPC **system** (e.g. queue limits, node sizes, charging model, default settings)
- The characteristics of your **data** (e.g. file size, sample number, processing steps)
- The **workflow** structure (e.g. the number of steps/processes, memory or CPU-intensive tasks)

!!! note

    Each HPC has different hardware and node architectures and it is up to the user to understand how their workflow fits within the limits

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