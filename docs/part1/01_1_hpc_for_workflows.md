# 1.1 HPC for bioinformatics workflows

!!! info "Learning objectives"

    - Know what an HPC is
    - Identify when a workflow requires an HPC
    - Understand some of the constraints that come with an HPC

## 1.1.1 What is an HPC?

High Performance Computing (HPC) systems are large **clusters** of computers with lots of CPUs, memory (RAM) and storage space. They are specially-built to run large numbers of computational jobs efficiently and concurrently. Bioinformatics analysis often involves many steps, many tools, many samples, and large datasets, which can quickly overwhelm your average laptop or desktop computer. When this becomes the case, HPCs can be the perfect solution to scaling up your workflows and running them efficiently and quickly. However, HPCs expect work to be submitted in a particular way, following specific rules. This means our workflows often need to be designed for HPC, not just moved to HPC.

## 1.1.2 When does a workflow need HPC?

In bioinformatics, a workflow is simply a defined series of steps that take data as input and transform that data into processed data and/or analytical results. This is true whether you are doing whole genome variant calling, proteomics quantification, single-cell transcriptomics, or metagenomics assembly. Each step in the pipeline performs one job, and each job depends on some form of computation and storage.

![](figs/00_workflow.png)

### Signs your workflow is ready for HPC

Not every workflow needs a supercomputer. Many analyses start on a laptop and stay there — especially during method development, testing small datasets, or when turnaround is more important than throughput. HPC becomes necessary when your workflow starts to hit practical limits of time, memory, storage, reliability, or governance.

A workflow is usually ready for HPC when scale becomes a problem. This might be scale in data size (more gigabytes than your laptop can hold), compute time (weeks of serial runs), memory usage (jobs crash due to insufficient RAM), or workflow complexity (tens of jobs become too painful to run manually).


| Challenge | Example scenario |
|-----------|----------|
| Runtime is too long | A single sample takes >12 hours to process |
| Data size is too big | Multiple large FASTQs need to be processed |
| Memory requirements are too large | R or Python crashes loading matrices |
| Scaling samples manually is painful | Manually running multiple scripts across multiple samples |
| Storage is a bottleneck | Local disk is constantly full due to raw and processed data size |
| Serial execution is too slow | Workflow is too slow to run one sample after another - multi-sample analysis must run faster |
| Data governance, ethics, and security constraints | Legal and/or ethical requirements mean highly-protected data must stay on institutional, secure systems |

If any of the above scenarios sound familiar, then your workflow is likely ready to be moved to and configured for running on an HPC.

## 1.1.3 From your laptop to HPC 

Before running a workflow, it is important to understand the system we are running it on. Running workloads on HPC is very different from running them on your laptop or a local workstation. HPCs are not just bigger, they are also: 

- Shared
- Scheduled
- Resource constrained. 

This introduces an important trade-off. HPCs give you access to massive computational power but at the cost of flexibility. On your laptop or a local workstation you can run whatever you like, whenever you like so long as it fits within the resource limitations of the system. On HPC, you gain scale and speed but you must work within system policies and limits. 

![](figs/00_hpc_use.png){width=70%}

### Shared systems

HPCs are large-scale institutional computing clusters that are intended to be used by many users at once. Indeed, their size and available resources mean than dozens or even hundreds of users can be using them at the same time and still manage to run large scale workflows concurrently and in a timely manner. However, this shared nature puts a significant constraint on how they can be used.

The primary constraint is that you don't have the freedom to install whatever software you want on the system. This means that you need other solutions to running the tools that you want. We will explore this issue in the [next section](./01_2_hpc_software.md).

Another constraint is the file system. While HPCs typically have huge shared file systems, they are neither infinite in size nor speed. Running workflows that generate lots of files, or read and write to the file system too frequently, will degrade the performance of the system for all users. Therefore, we need to be conscious of what our workflow is doing and make sure we design it to use the system fairly and efficiently. We will discuss storage limitations in [1.3 HPC architecture](./01_3_hpc_architecture.md).

### Schedulers

On your local laptop, you will be used to running things whenever you like, but on shared systems like HPCs, this is not the case. Instead, HPCs require you to **submit** jobs to a **scheduler**, which decides where and when to run your job based on its resource requirements and the requirements of all other jobs in the **queue**. This makes HPCs **asynchronous** and **non-interactive**: job execution doesn't happen immediately and jobs won't necessarily execute in the order that they were submitted. As such, an HPC workflow needs to be designed to handle this delayed and potentially out-of-order execution style. As we will see later today, Nextflow is perfectly suited to writing workflows that work in this way.


### Resource constraints

Finally, HPCs may have large amounts of computing resources, but they aren't infinite, and they also need to be shared between many users. Therefore, it is vital when running jobs on an HPC to define exactly how many resources you require, including the number of CPUs you need, the amount of memory/RAM, and how much time your jobs needs. As you will see later in this workshop, it is very important to optimise these requests as best as you can, as under- and over-requesting resources can negatively impact your jobs.

## 1.1.4 Introducing our workshop scenario: WGS short variant calling
 
!!! warning "Don't worry if you don't have prior knowledge of this workflow"

    The focus of this workflow is on learning Nextflow; the experimental context we are using (WGS short variant calling) is just a practical example to help you understand workflow design principles for HPC and how Nextflow works. You are not expected to have prior knowledge of variant calling workflows or best practices.

For this workshop, we will be focussing on a common bioinformatics analysis workflow used in genomics to identify genetic variants (SNPs and indels) from short-read whole genome sequencing data. This workflow involves multiple processes and tools and is computationally intensive. At a high level, the general procedure is: 

![TODO add fig](../img/wgs_variant_calling.png)

1. Quality control of raw sequences, e.g. filtering & trimming reads
2. Alignment of reads to a reference genome
3. Post alignment processing: sorting, marking duplicates, indexing
4. Variant calling: call SNVs and indels for each sample against reference
5. Joint genotyping: combining samples from a cohort into a single callset
6. Quality filtering and annotation

Running this workflow end-to-end captures many challenges that running on HPC using Nextflow can solve: 

- Many independent jobs: each sample can be processed separately for many steps 
- Resource diversity: tools used at each step require different amounts of CPU, memory, and walltime  
- Large IO demands: reading and writing of multi-gigabyte files benefits from parallel filesystems 

Throughout the workshop we will implement and explore different parts of this workflow in slightly different ways in order to highlight the lessons being taught. 

!!! example "Discussion: why does this workflow need HPC?"

    Consider the workflow described above: 

    1. Which parts of the workflow are the most computationally expensive? 
    2. What would happen if we tried to run this workflow on a personal computer?  

## Conclusion

TODO