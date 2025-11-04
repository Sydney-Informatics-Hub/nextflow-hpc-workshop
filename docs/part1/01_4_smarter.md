# 1.4 Work smarter, not harder  

!!! info "Learning objectives"

    - Understand how HPC resources affect job scheduling and performance  
    - Describe how HPC scheduling and resource limitations shape pipeline configuration.
    - Differentiate between multithreading and scatter–gather approaches to parallelisation  
    - Identify how to size resource requests appropriately for each process in a workflow  
    - Apply resource-aware design principles to improve job efficiency and reproducibility  

HPC systems give us access to large amounts of compute, but that doesn’t mean we should use resources carelessly. Misusing compute leads to long queue times, wasted allocation, unstable workflows and unhappy HPC administrators. Designing resource-aware workflows is essential for performance and fair use. 

HPC systems are constantly measuring your resource usage. You can use their built in tools to measure actual use. The tools available to you will depend on the job scheduler and the administrator's implementation of the scheduler. 

!!! example "Inspect a previous job"

    TODO design an exercise for checking the previous lessons' jobs.
    === "Gadi"
        ```bash
        qstat -xf <jobid> | grep -E "Resource_List|Used"
        ```
    === "Setonix"
        ```bash
        sacct -j <jobid> --format=JobID,JobName,Elapsed,State,AllocCPUS,TotalCPU,MaxRSS
        seff <jobid>
        ```

## 1.4.1 Resource awareness: right sizing 

At its core, HPC efficiency is about matching the structure of your workflow to the available compute. It is therefore beneficial to be "resource aware" in your approach to running jobs. This involves understanding how much time, CPU, and memory each tool in your workflow actually needs and requesting enough. 

Here we focus on 3 metrics: 

- Elapsed time: the total walltime for the job 
- CPU time: how much CPU time was actually used 
- Maximum memory (RAM): the peak memory footprint of the job 

In the context of running a workflow made up of multiple steps, each running a different tool, we consider each step separately in optimising for efficiency. 

Let's look at our workflow:  

![](./figs/00_workflow01.png){width=80%}

TODO fix this table to make it make sense/honest/relevant. 

| Step | Dominant resource | How to tell |
|------|-------------------|-------------|
| **Quality control** | I/O-bound | Reads many small files; CPU idle time high |
| **Read alignment** | CPU-bound | CPU pegged near 100%; memory stable |
| **Variant callingr** | CPU + memory | CPU ~90%, high steady memory usage |

We will observe these constraints in subsequent lessons when we run and optimise our workflows.

!!! note "Where can I find this information?" 

    Most mature bioinformatics tools document their approximate resource usage. Keep in mind, documentation will not reflect the reality of your specific dataset and environment. Look at your own resource usage and scaling behaviour to test this.  

    The more you practice tuning processes for efficiency, the faster you'll develop an intuition for scaling. 

### CPU efficiency 

This refers to how well the CPU requested for the job are being utilised. It is commonly expressed as a percentage or range from 0-1. It is calculated as: 

```
CPU time / walltime / number of CPUs 
```

Values near 1 mean your job used all the CPUs efficiently. Values much lower than 1 suggest your job was waiting on I/O, memory, or over-allocated resources.

### Memory (RAM) efficiency 

Memory efficiency describes how much of your allocated memory was truly needed. If your job used close to the requested memory, you're right sized. If maximum memory is much lower than requested, then you over-allocated, and your job will fail if maximum memory exceeds the requested allocation.    

### Walltime awareness 

Walltime is the time it takes for your job to run from start to finish. Schedulers use this value to plan future jobs being run by yourself and others. Overestimating walltime will keep you job in the queue for longer. If you underestimate walltime, the job will be killed mid-run.

## 1.4.2 Optimising your jobs for efficiency

Once you understand how your workflow uses resources, you can start to optimise it. Optimisation is about balancing speed, efficiency, and fair use of the system. 

In bioinformatics workflows, there are 2 main strategies used for increasing efficiency and throughput: 

| Level           | Definition       | Where it runs  | Memory model                | Example                |
| ------------------------------- | --------------------------------- | -------------------------- | --------------------------- | --------------------- |
| **Multi-threading**              | Multiple threads share memory space inside a single process. Threads cooperate to perform a single task faster.                        | One node                   | Shared memory               | Using the `-t` or `--threads` flag on a tool        |
| **Multi-processing**             | Multiple independent processes run in parallel, each with its own memory space. Often managed by a workflow engine or batch scheduler. | One node or multiple nodes | Separate memory per process | Running `fastqc` on multiple samples simultaneously        |

We can explore parallelisation methods of multi-threading and multi-processing in the context of our variant calling workflow and some small dummy data. 

!!! warning "Beware of the overheads" 

    While parallelism can be a great way to speed up your data processing, it doesn't always make things run faster. Splitting your jobs should only be done where it makes practical sense to do so. Keep in mind, the more you split your work up, the more issues you will have to contend with: 

    - More job scheduling, more file-handling, more merging 
    - Increased memory footprint and I/O for each sub-job 
    - Diminishing returns when the time saved is offset by coordination cost 

    See [this great explainer](https://gatk.broadinstitute.org/hc/en-us/articles/360035532012-Parallelism-Multithreading-Scatter-Gather) of parallelism from the GATK team.

## Parallelisation: multi-threading 

!!! note "What is multi-threading?"

    Multithreading means one program is using multiple cores on the same node to complete a single task faster. Our ability to do this is dependent on the tool being used. Some bioinformatics tools support multithreading via the `--threads` or `-t` flag.

In our variant calling workflow, some tools work best when given multiple cores on the same node, e.g. alignment with `bwa mem`. When you run: 

```
bwa mem -t 4 ref.fasta sample_R1.fq.gz sample_R2.fq.gz > alignment.sam
```

You're telling BWA to use 4 worker threads for parts of the alignment process that can be parallelised. Not all parts of the process can be parallelised. 

!!! example "Does more threads always speed up a job?" 

    ```
    TODO time bwa mem
    ```

 Not all tools can make use of threads effectively, and some don't implement the `--threads` flag effectively at all. 

!!! example "Does more threads always speed up a job?" 
    
    While many tools benefit from multi-threading, [FastQC is not one of them](https://www.biostars.org/p/9598335/). It is designed to process one file per thread rather than splitting a single file across multiple threads.

    ```
    TODO fastqc example
    ```

    So here, increasing threads does not accelerate processing of a single fastq file. 

!!! warning "Threads ≠ cores " 

    - **A core** is a physical processing unit on a computer node, it can execute one thread at a time. 
    - **A thread** is a software-level unit of execution, a single stream of instructions within a process

    In practice, when you run a command like: 

    ```bash
    bwa mem -t 8
    ```

    You're telling bwa to spawn 8 threads. For those threads to be executed in parallel, you must request 8 cores to run those threads on. If you were to request only 4 cores, those 8 threads will be competing for 4 cores which causes inefficiency and slower run times. If you ask for more CPUs than threads, those extra CPUs will sit idle, wasting allocation and increasing queue time. 

## Parallelisation: multiprocessing 

!!! note "What is multi-processing?"

    Multi-processing means running independent processes at the same time, each with its own memory space and input data. Rather than speeding up a single task, it increases overall throughput by running many tasks in parallel. In bioinformatics workflows, this usually looks like processing multiple samples or genomic regions simultaneously. 
    
    Multi-processing is typically implemented by the user in workflow design, rather than by a tool itself.

![](./figs/00_workflow02.png){width=80%}

TODO maybe change this diagram to focus on one process and show different ways it can be implemented? 

!!! example "TODO an exercise" 

    ```
    An exercise in workflow design
    ```

!!! warning "Does it make sense biologically to process in parallel?" 

    Not all parallelisation makes sense, it depends on what you're analysing and how the tool interprets the data. Parallelisation makes sense when the data are independent, e.g. running FastQC on multiple fastq files, aligning multiple samples with bwa-mem, or calling SNPs and indels per chromosome. 

    Paralellisation does not make sense when results depend on comparing all data together e.g. joint genotyping, genome assembly, or detecting structural variants across multiple chromosomes. 