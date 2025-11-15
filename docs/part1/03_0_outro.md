# 1.10 Summary

In today's session we have learned about what HPCs are, what sets them apart from your standard computer and how they require you to work in a particular way, by explicitly requesting the resources you need and submitting jobs to a scheduler to be executed whenever those resources become available. We have also learned how Nextflow uses executors to talk to the scheduler and manage submitting the individual workflow tasks to be executed on the HPC compute nodes.

## HPCs are shared systems

Remember, you are working on a shared system with hundreds of other users. Because of this, you need to be explicit about how many resources your jobs require. The scheduler's job is to take those requests and find a suitable time and place to run the jobs. It is important to optimise your pipelines; if you under-request resources, you job will fail; but if you over-request, your job may wait in the queue for longer and cost more than it needed to.

## Containers are your friend

Controlling the verions of the software you use is vital for reproducible research. Out of the various different methods for running software on HPCs, containers provide the most freedom, allowing you to use any version of a tool you would like, without having to worry about dependency issues, version conflicts, or building the software yourself. Their self-contained nature and portability makes containers the recommended method for using software within your Nextflow workflows.

## Nextflow is portable

The goal of Nextflow is to allow for the development of portable and reproducible workflows that are written once and run anywhere. To achieve this, they separate out the actual workflow logic from the configuration. When properly written, the Nextflow code itself should never need to be altered when moving the workflow between systems, whether that be a local laptop, an institutional HPC, or the cloud. Instead, one or more configuration files can be used to handle system-specific details, such as communicating with an HPC scheduler. Remember that there are different levels of configuration that a pipeline will typically require, and each should have their own configuration files:

- Pipeline-specific configuration details should be kept in `nextflow.config` and other `.config` files bundled with the workflow code. These should never need to be altered from system to system.
- Institutional-specific configuration details, such as HPC executor configuration, should be kept in a separate file (e.g. `gadi.config` or `setonix.config`) that can be used by multiple different pipelines. Ideally, this shouldn't need to be altered between pipelines on the same system.
- Run-level configuration details, such as fine-tuning memory and CPU requriements for a particular dataset, should be kept in yet another configuration file that is specific to that instance of that pipeline.

## Fine-tuning workflows for your data

While a well-written pipeline should be able to handle inputs of different sizes and dynamically request resources accordingly, it is impossible to forsee every possibility. As such, it is often necessary to have a custom configuration file to fine-tune a pipeline to your specific dataset. We saw today that for very small datasets, a pipeline may by default request too many resources than necessary. In other cases, for a particularly large dataset, you may need to increase the resources you require to prevent out-of-memory failures and jobs exceeding their requested walltime.

## Next steps: configuring a custom pipeline

While nf-core is an amazing resource for community-built bioinformatics tools, it still has its limitations. There isn't a pipeline for every purpose you might have, and when there is, you might find it is still under development, or contains bugs that break your analyses, or simply doesn't quite do what you need in just the right way. As such, it is often necessary to write your own custom pipelines that handle your data in just the way you like. In tomorrow's session, we will apply the concepts we learned today to configuring such a custom workflow for HPC. We will be continuing along the same theme as today: WGS short variant calling. In doing so, we will explore how we can update a pipeline's logic to better take advantage of the parallel nature of HPCs and implement our own scatter-gather pattern for parallel processing of sequencing data.