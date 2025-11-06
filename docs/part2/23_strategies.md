# Overcoming performance bottlenecks with HPC optimisation strategies

!!! info "Learning objectives"

    - Understand how HPC resource requests affect job scheduling and performance
    - Identify process and configuration bottlenecks in trace reports
    
## Overview on optimisation strategies

Recap on Part 1.4 Work smarter, not harder

Key strategies:

1. Right sizing: match resources to what the process actually uses (CPU, memory, time)
2. Multithreading
3. Multi-processing

## Identifying process and configuration bottlenecks

!!! example "Exercise"

    TODO Review trace files, identify processes that are hitting close to 100% CPU/memory usage.

    TODO Any tools that support multithreading? Is it always good to apply?

    TODO Review workflow structure - can anything be broken down and run in parallel?


!!! question "Poll"

    Have a blank table, ask attendees to put "dots" on processes that can be
    optimised? PollEv interactive image with regions

    | Process        | Resourcing       | Multi-threading | Multi-processing (scatter-gather) |
    | -------------- | ---------------- | --------------- | --------------------------------- |
    | FASTQC         | Yes              | Yes - but..     | Yes                               |
    | ALIGN          | Yes              | No              | Yes                               |
    | GENOTYPE       | Yes + directives | Yes             | NO                                |
    | JOINT_GENOTYPE | Yes + directives | Yes             | NO                                |
    | STATS          | Yes              | No              | NO                                |
    | MULTIQC        | Yes              | No              | NO                                |

TODO: provide image for Polls

### The false positive and negative answers

Discuss!

### The No's

- JOINT_GENOTYPE, STATS, MULTIQC cannot be split
    - 1. Biology: Variants need to be called across all samples in the cohort in the scaling up section
    - 2. Technical: MULTIQC needs all the files to compile the report
- Multi-threading is not supported by all tools.

### The Yes'

- Resourcing is extremely important when working on HPC:
    1. it fits the requirements for the system infrastructure(queues or partitions)
    2. so you jobs are scheduled quickly
    3. so your jobs are run efficiently and request the right number of resources

- For the ones that support multi-threading... should we..? how many threads?

- FASTQC-ALIGN can be scattered, done in parallel (leveraging HPC resources, lowering walltime), and brought back together

These three points will be revisited in the next sections sequentially.
