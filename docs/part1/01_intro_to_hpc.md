# 1.1 Introduction to HPC

Before we discuss running Nexflow on HPCs, we must first understand the basic concepts of these systems. Running jobs on an HPC is a bit different to running them on a laptop or desktop, and careful consideration needs to be made with respect to the resources we need.

## 1.1.1 Anatomy of an HPC

There are several different HPC systems that exist, but they all follow the same basic structure and organisation. They typically consist of a **cluster** of separate computers, called **nodes**, which are networked together and can talk to each other. There are two main types of nodes: **login nodes** and **compute nodes**. The login nodes are relatively low-powered computers that users connect to directly, while the compute nodes typically have large amounts of memory and high numbers of CPUs available for use, and are where the jobs are actually run. From the login nodes, useres can **submit** jobs to the HPC **scheduler**, which orchestrates which jobs are run on which compute nodes, and when. The HPC system will also typically have a very large networked data storage system available, which all of the nodes can read from and write to.

![Anatomy of an HPC](/docs/assets/hpc_anatomy.png)

The structure of HPCs has some huge advantages. First, you can run compute-heavy tasks that you couldn't otherwise run on a local machine. But importantly, jobs that are independent of one another - for example, processing separate samples - can be run in parallel rather than sequentially. This can significantly speed up run times for these jobs. However, if you are more familiar with running tools interactively via the command line, the HPC architecture can seem a little off-putting, as it requires constructing **scripts** for each job and **submitting** them in a non-interactive manner. This also means that it is often difficult to monitor jobs in real-time, and instead you will often need to read through log files at the end of the job.

## 1.1.2 HPC resources

A major difference between running a job on an HPC and running it locally on your laptop or desktop is that you need to specify the resources that your job requires. That is, you need to tell the HPC's scheduler how many CPUs you want, how much memory you need, and how long your job needs to run. The scheduler needs this information so that it can efficiently orchestrate where and when your job should run amongst all the other jobs from all the other users.

It is important to optimise the resources for your jobs as best as you can. It is not a good idea to simply request as much as possible, for a number of reasons.

1. First, it can be quite expensive! HPCs typically charge you for the amount of resources you use. Requesting too much can be unnecessarily expensive!
2. Second, requesting too much can actually slow down your work! Suppose, as an extreme example, you have a job that requires 1 CPU and will run in approximately 10 minutes; but instead, you request 8 CPUs and 12 hours! Now, the scheduler needs to find a node that has 8 CPUs available for a 12 hour stretch. With lots of users using the system simultaneously, you might be waiting several hours for a suitable slot to open up. On the other hand, if you had requested just 1 CPU and 10 minutes, the scheduler may have been able to find a spot for your job immediately!
3. Third, each node has a finite amount of resources available. By requesting too many resources, you reduce the amounts available to other users on the system.

It is also important to submit your job to an appropriate **queue**. HPCs typically have different queues for tasks with different resource requirements. Some queues have lots of resources available for large, long-running tasks, while others are optimised for smaller, regular tasks.

![HPC scheduling](/docs/assets/scheduler.png)

Of course, it is also problematic if you don't request enough resources for your job. If you request too few CPUs, too little memory, or not enought time, your job may fail or stop prematurely. As such, benchmarking of your jobs is typically required to tune the resources needed carefully. This is typically done by first running small-scale tests of your workflow and extrapolating the resources required for larger, full-scale datasets.

Remember, HPCs are a **shared resource**, and operate on a **fair share** policy. We as users request a certain amount of resources for our jobs and the scheduler attempts to fairly distribute those resources to all users. By requesting appropriate amounts of resources for your tasks, and by submitting your jobs to the appropriate queues, you can ensure that they move as quickly as possible through the queue, while maintaining fair use of the system as a whole by all users.

## 1.1.3 The do's and do not do's of HPCs

HPCs are designed to be used by many users simultaneously. This means that we have to be mindful of how we use these systems.

First and foremost, you should never run compute-heavy tasks on the login nodes. It may sometimes be tempting to run a tool or a script interactively in the login nodes, but keep in mind that: 1. the login nodes are relatively low-powered compared to the compute nodes; and 2. many people are likely logged in at the same time and trying to use the same login node for their work. This means that compute-heavy tasks can quickly use up all the login node resources, and will either significantly slow down the node, or fail altogether.

Another important consideration when working on HPCs is intermediate file storage. Some tasks can create dozens or even hundreds of intermediate files. Constant reading and writing from the file system can also slow things down for both yourself and other users. As such, there are usually different levels of storage systems available for use by your jobs.

- Main network storage. Most HPCs will have a primary network storage, which is intended for storing large input and output files for your project. This storage class is typically the slowest, and should be used for files that aren't frequently written to. For example, if you have DNA sequencing data in FASTQ format and intend to map the reads to the genome and conduct variant calling, the primary network storage is a good place to store your input FASTQs, as well as the output BAM and VCF files from your analyses.
- Scratch storage. This is also network storage available to all nodes, and is typically faster than the primary storage, but is also **temporary**; many HPCs will have a deletion policy in place for scratch storage to ensure it doesn't fill up. Scratch is intended for temporary files, not for long-term storage. It can be a good place to run your jobs, before moving the outputs to the primary storage space.
- Node storage. This is the storage of the compute node itself, and is often much faster than either primary or scratch storage. Most HPCs will let you access a temporary directory on the node, which is extremely useful for writing lots of intermediate files during your run, that don't need to be kept afterwards. For example, many bioinformatics tools create lots of small intermediate files during the run, which then get assembled into a single large file at the end. The intermediate files don't need to be kept, so there is no sense in writing them to network storage, which would unnecessarily slow down the system due.

| HPC component | Used for | Should not be used for |
| --------- | -------- | ---------------------- |
| Login node | Submitting jobs. Short-running, low-intensive tasks (e.g. text editing, creating scripts) | Running compute-intensive tasks (e.g. bioinformatics tools). |
| Primary network storage | Large input/output files. | Lots of small intermediate files. Heavy I/O operations. Temporary data storage. |
| Scratch storage | Temporary data storage. Medium I/O operations. Intermediate files that you might want to keep short-term. | Heavy I/O operations. |
| Node storage | Heavy I/O operations. Intermediate files that don't need to be kept after the job is complete. | Temporary or long-term storage. |