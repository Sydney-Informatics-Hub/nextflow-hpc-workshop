# Optimising Nextflow for HPC

!!! info "Learning objectives"

    - Identify the key factors that impact workflow performance on HPC systems
    - Describe common workflow optimisation strategies in the context of Nextflow



Our workflow is now functional: it runs successfully with scheduler-specific settings on HPC and outputs useful trace files showing resource usage. In [Lesson 1.4](../part1/01_4_smarter.md), we learnt about the importance of appropriate resource requests on HPC and the value of running an efficient and optimised workflow. We will now begin to apply these principles to our custom workflow. 



## 2.3.1 Why optimise?

Workflow optimisation involves fine-tuning your workflow to make it more efficient. It can be used to reduce runtime, facilitate higher throughput and larger more powerful studies, avoid idle hardware that could be freed up for other researchers, and decrease computing cost. Many HPCs charge use per core hour, yet even for those where the cost is not billed to the researcher, the cost of the resources must be covered by someone, and this may be the Faculty/School, institution, funding body, or federal budget. Demonstrating efficient use of compute funds can be important in grant success for projects requiring substantial compute. 

Another important reason optimisation is worthwhile is the impact on the environment. Bioinformatics data processing on HPCs has a significant carbon footprint due to the energy required to run our computations. Making our workflows complete in a faster time using less hardware contributes to sustaintable computing practices. 

*[This short editorial](https://www.nature.com/articles/s43588-023-00506-2) from Nature Computational Science highlights the challenges research computing faces in the context of the climate crisis*.

Today we will apply workflow optimisation from two perspectives:

* Resource efficiency: benchmarking and adjusting resources for efficient process execution on HPC
* Speed-up: introducing code changes to make the workflow run faster without costing more compute hours

## 2.3.2 When to optimise 

A small workflow run a handful of times might not benefit dramatically from optimisation. Many Nextflow workflows that employ good practices (e.g. nf-core) will often run "out of the box" on HPC with default resources, but defaults might not always fit your data and therefore the behaviour of your processes, or the constraints of your cluster. Think back to Part 1 and the configuration customisations we implemented for our nf-core workflow. 

Optimising workflows on HPC becomes especially important when:

- You are analysing large datasets or many samples 
- You will execute your pipeline repeatedly 
- You are operating on a system that uses a service unit or time-limited allocation model
- Your processes have data-dependent resource usage 

The ideal time to optimise a workflow is while it is being developed. This is often more simple to achieve than back-tracking and adding improvements to an existing and potentially large codebase, and will also ensure that the benefits of an efficient workflow are enjoyed as early as possible. Notably, this may add development time to producing a finished workflow, but the efforts are rewarded with a resilient and scalable workflow that will reduce the time interval between data acquisition and final results.

If you have an existing workflow that is reliable but inefficient, there is always value in taking the time to optimise your regularly used workflow. This endeavour also provides the opportunity to update tool versions and introduce other enhancements such as the use of Singularity containers.

When optimising a workflow, resource efficiency should always be considered. Optimisation through parallelisation however is not always possible or recommended. It is vital to always consider what components of work can be parallelised or split up in a ***biologically valid*** way. If you are unsure, one means of testing is to run the analysis with and without parallelism on a subset of data to observe any impact on results. Due to heuristics of tool algorithms, in many cases a small amount of difference in the final result is valid and tolerated. In other cases, it may be expected that a tool produces identical results irrespective of multi-threading or parallelisation. It is important to check the tool documentation, consider the nature of the data and underlying biology, test the effects of parallelisation, and only apply parallelisation when it makes biological sense to do so. 

Consider the example of sequenced DNA fragments in whole genome sequencing. Each fragment is completely independent on every other fragment in the library, so it can be aligned to the reference sequence independently without affecting the mapping results. This is a perfect case of "embarrassingly parallel" processing, i.e. running the same analysis numerous times on slightly different input data. In reality, we would not align one read at a time as the overhead of submitting millions of tiny jobs and loading the same tools and reference over and over would overload the scheduler, impact performance, and would also likely result in stern reproach from the HPC system administrators! So in this case, we have biological validity to apply scatter-gather parallelism, but we need to also apply sound HPC practices to ensure the level of parallelism is compatible with the HPC and produces an actual speed-up. 

Now consider another example of identifying large-scale changes to an organism's genome. These are known as *structural variants* and can include millions of nucleotides and span multiple chromosomes. Do you think it is biologically valid to split this up into chunks for parallel processing? The answer is no, because the tool needs to "see" all of the reads and all of the reference genome at once to be able to describe such a variant. In cases like this, the best level of parallelism we can achieve is **parallel by sample**. Parallel by sample is logical for numerous common bioinforamtics processing tasks, and is typically only invalid when all samples must be analysed together for example when collating final results.  

## 2.3.3 What affects performance?

Efficiency of any workflow on HPC depends on the interaction of three factors: 

### 1. Your HPC system 

We have already witnessed many differences between Gadi and Setonix in previous lessons. A workflow that performs well on one cluster may perform poorly on another simply because the underlying architecture and scheduler rules differ. 

Good optimisation respects the boundaries of the system you're working on. When planning an optimisation approach, consider:

| HPC characteristic | What it means | Why it matters for optimisation                     |
| --------- | ---------------------- | ------------- |
| **Default scheduler behaviour**         | Policies set by administrators: fair-share, job priorities, backfill rules, default limits  | Affects queue wait time, job placement efficiency, and how many tasks can run in parallel                            |
| **Queue limits**                        | Maximum walltime, cores, and memory allowed per queue or partition                      | Determines which queues you can use, how large each job can be, and whether your workflow gets delayed               |
| **Node architecture**                   | Hardware layout: cores per node, memory per node, CPU type (Intel/AMD), GPUs, local scratch | Ensures you request resources that “fit” the node, avoid resource fragmentation, and maximise throughput             |
| **Charging model** | How HPC usage is accounted (CPU proportion, memory proportion, or the maximum of both)      | Guides you to request only what you need: over requesting directly increases SU consumption without improving runtime |

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
* Whether steps are CPU-bound, memory-bound, or I/O bound 
* Incorporated tool's ability to multithread 

## 2.3.4 What will we optimise today?

For the remainder of Part 2 we will apply the strategies introduced from Part 1 to optimise then scale our custom workflow. 

In the next section, we will **assign appropriate resources for each process** by using trace files to fine-tune `cpus`, `memory`, and `time` and align these to the resources on the compute nodes of our HPCs. 

In Lesson 2.5, we will introduce **parallel processing** firstly by enabling **multi-threading** in a thread-aware tool (BWA), and then by coding **scatter-gather parallelism** into the workflow. 

In today's example workflow, we will be applying scatter-gather to run alignment with BWA. Note that the `GENOTYPE` process can also be parallelised in a biologically valid way using a parallisation strategy known as "interval chunking", but for simplicity we will not be optimising that process today. 

Finally, in Lesson 2.6, we will **scale to multiple samples**. This will consolidate all of the resource optimsiations and parallelisation strategies (multi-threading; scatter-gather; parallel by sample) that we have built up today into one efficient, end-to-end run optimised for our respective HPCs. 

