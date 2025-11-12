# Optimising for HPC

!!! info "Learning objectives"

    - Understand how HPC resource requests affect job scheduling and performance
    - Identify process and configuration bottlenecks in trace reports

## Overview

Our configuraiton is now functional! The workflow runs sucessfully with scheduler-specific settings and outputs useful trace files showing resource usage.

Now, we shift focus from making it run, to making it efficient.

TODO: figure benchmarking progress e.g. it runs, it runs well, it scales

Optimisation depends on several key factors, such as:

- The specific HPC system (e.g. queue limits, node sizes, charging model)
- The characteristics of your data (e.g. file size, sample number, processing steps)
- The type of processes your workflow runs (e.g. high-memory)

The same principles can be applied to any HPC environment, but in this workshop, we’ll focus on NCI Gadi and Pawsey Setonix as practical examples.

Each HPC has different hardware and node architectures and it is up to the user to understand how their workflow fits within the limits. 

Small inefficiencies for small, or infreuqently-run workflows may not make a large impact. Many users run nf-core pipelines on HPCs successfully, using the default configuration.
     
However, if you are running high-throughput workflows that consume a considerable amount of service units, it is worthwhile to benchmark and optimise your workflow for the HPC you will run it on.

## 

For the remainder of Part 2 we will apply the various optimisation strategies introduced in Part 1. We'll build off the benchmarks to configure our custom pipeline for HPC.

1. Right resourcing for each process, to fit your assigned HPC
2. multithreading bwa mem
3. Multiprocessing alignment using scatter-gather patterns in Nextflow