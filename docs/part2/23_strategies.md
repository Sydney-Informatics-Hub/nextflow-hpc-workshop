# Optimising Nextflow for HPC

!!! info "Learning objectives"

    - Identify the key factors that impact workflow performance on HPC systems
    - Describe common workflow optimisation strategies

## Overview

Our workflow is now functional, it runs successfully with scheduler-specific settings and outputs useful trace files showing resource usage. From here, we will explore ways of optimising our workflow for more efficient execution on the HPC. 

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