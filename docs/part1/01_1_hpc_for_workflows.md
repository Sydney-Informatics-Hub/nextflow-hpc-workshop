# 1.1 HPC for bioinformatics workflows

!!! info "Learning objectives"
    TODO revise these

    - Describe the HPC system components that workflows interact with.
    - Identify why containerisation and resource-aware design are essential for scalable workflows.
    - Describe how HPC scheduling and resource limitations shape pipeline configuration.
    - Connect HPC principles to the Nextflow workflow management systems.

High Performance Computing (HPC) systems are built to run large numbers of computational jobs efficiently. Bioinformatics analysis often involves many steps, many tools, and many samples, making it a perfect match for HPC. However, HPCs expect work to be submitted in a particular way, following specific rules. This means our workflows often need to be designed for HPC, not just moved to HPC.

## 1.1.1 When does a workflow need HPC?

In bioinformatics, a workflow is simply a defined series of steps that take data as input and transform that data into processed data and/or analytical results. This is true whether you are doing whole genome variant calling, proteomics quantification, single-cell transcriptomics, or metagenomics assembly. Each step in the pipeline performs one job, and each job depends on some form of computation and storage.

![](figs/00_workflow.png)

### Signs your workflow is ready for HPC

TODO this is not very good, can come up with some better examples here, can be communicated better. 

Not every workflow needs a supercomputer. Many analyses start on a laptop and stay there—especially during method development, testing small datasets, or when turnaround is more important than throughput. HPC becomes necessary when your workflow starts to hit practical limits of time, memory, storage, reliability, or governance.

A workflow is usually ready for HPC when scale becomes a problem. This might be scale in data size (more gigabytes than your laptop can hold), compute time (weeks of serial runs), memory usage (jobs crash due to insufficient RAM), or workflow complexity (tens of jobs become too painful to run manually).


| Challenge | Scenario |
|-----------|----------|
| Runtime is too long | A single sample takes >12 hours to process |
| Data size is too big | Multiple large FASTQs to be processed |
| Memory limits hit | R or Python crashes loading matrices |
| Scaling samples manually is painful | Running multiple scripts across multiple samples |
| Storage is a bottleneck | Local disk constantly full |
| You need parallel execution | Multi-sample analysis must run faster |
| Workflow reliability matters | Need checkpointing and recovery |
| Data must stay on institutional systems | Governance, ethics, security |

## 1.1.2 From your laptop to HPC 

Before running a workflow, it is important to understand the system we are running it on. Running workloads on HPC is very different from running them on your laptop or a local workstation. HPCs are not just bigger, they are also: 

- Shared
- Scheduled
- Resource constrained. 

This introduces an important trade-off. HPCs give you access to massive computational power but at the cost of flexibility. On your laptop or a local workstation you can run whatever you like, whenever you like so long as it fits within the resource limitations of the system. On HPC, you gain scale and speed but you must work within system policies and limits. 

![](figs/00_hpc_use.png){width=70%}

## 1.1.3 Our experiment: WGS short variant calling
 
!!! warning "Don't worry if you don't have prior knowledge of this workflow"
    The focus of this workflow is on learning Nextflow, the experimental context we are using is just a practical example to help you understand workflow design principles for HPC and how Nextflow works. 

We are going to implement a common pipeline used in genomics to identify genetic variants (SNPs and indels) from short-read whole genome sequencing data. Through the workshop we will implement this workflow in slightly different ways. It involves multiple processes and tools and is computationally intensive. At a high level, it does the following: 

TODO add fig 

1. Quality control of raw sequences 
2. Alignment of reads to ref genome 
3. Post alignment processing: sorting, marking duplicates, indexing 
4. Variant calling: call SNVs and indels for each sample against reference 
5. Joint genotyping: combining samples for a cohort into a single callset 
6. Quality filtering and annotation 

Running this workflow end-to-end captures many challenges that running on HPC using Nextflow can solve: 

- Many independent jobs: each sample can be processed separately for many steps 
- Resource diversity: tools used at each step require different amounts of CPU, memory, and walltime  
- Large IO demands: reading and writing of multi-gigabyte files benefits from parallel filesystems 

!!! exercise "Discussion: why does this workflow need HPC?"
    Consider the workflow described above: 

    1. Which parts of the workflow are the most computationally expensive? 
    2. What would happen if we tried to run this workflow on a personal computer?  

## Conclusion 