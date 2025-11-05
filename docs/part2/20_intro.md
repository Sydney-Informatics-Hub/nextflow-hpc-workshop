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
throughput.

What part 2 is NOT:
- How to write and build nextflow pipelines and logic
- Understanding the contents input and output files in detail

## Why do you need custom pipelines?

Reasons for using custom pipeline:
- Fit for your own purposes
- Existing pipelines (nf-core) may not do everything required, pipelines
doesn't exist
- it doesn't exist!

## The pipeline file anatomy

TODO: `tree` of part2 structure - showcasing locations of /conf, /modules,
main.nf, nextflow.config

### `main.nf` and `modules/`

TODO: snippet of `main.nf`

Set up as modules - a recurring theme is that there are Nextflow files that can
be left as is, and some that need tweaking. This makes your pipelines portable,
organised, reproducible, and easy to set up across different systems.

Here, we will not edit the modules, but configure how and where to run them on
HPCs.

- mentioned modules imported.
- Why modules?
- During scatter-gather later, easy to swap out modules, rather than coding
in processes in main.nf. Reinforce reproducibility, not all files need to be
touched
- https://sydney-informatics-hub.github.io/template-nf-guide/notebooks/modules.html

### `nextflow.config` and `conf/`

We start with basic preconfigured - executors, the default queues, and
singularity enabled.

As well as the default parameters for the workflow to run.

TODO: snippet of starting `nextflow.config` and `conf/pbs-or-slurm.config`
- This is what we will be starting with, building up towards the optimised
pipeline by end of part 2
- Show for pbspro and slurm

