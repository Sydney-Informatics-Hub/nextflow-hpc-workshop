# 1.3 Parallelisation on an HPC

HPCs are clusters of numerous nodes, each with a high number of CPU cores available for use. This makes them excellent environemnts for parallel workflows; that is, workflows with many independent tasks that can be run simultaneously, rather than sequentially. This can considerably speed up large workflows that may otherwise take hours or even days to complete on a standard computer.

## 1.3.1 Parallel workflows

The most obvious use case for paralellisation is to process multiple samples through the same workflow at the same time. If you have a workflow that requires 4 CPUs and 8GB of memory, then this will create a bottleneck on a standard computer with limited resources, and likely you will have to process each sample one-by-one, sequentially. But, on an HPC, with hundreds or even thousands of CPUs and many TB of memory available, suddenly you have an opportunity to run each sample in parallel to one another.

Another very common example in bioinformatics is analysing sequencing reads that map to different chromosomes in the genome. For example, suppose you want to perform single nucleotide variant (SNV) calling on genomic sequencing data. This involves looking at all the reads that pile up at each position in the genome and determining whether or not the base at that position differs from the reference genome. Since chromosomes are physically separate from one another, you can **parallelise** this task by initially splitting your reads into one file for each chromosome, then performing the variant calling task on each chromosome separately and simultaneously. Once all the chromosomes have been processed, you can later concatenate the results back together again. This is a common pattern known as **scatter-gather**: you **scatter** your data into independent tasks that can be run simultaneously, then later **gather** the results back together again.

![The scatter-gather pattern](/docs/assets/scatter_gather.png)

## 1.3.2 The limits of parallelisation

Keep in mind, however, that not everything can take advantage of the parallel nature of HPCs, and some tasks will still need to be sequential in nature. In the scatter-gather pattern, for example, the gathering step naturally requires waiting for all of the parallel tasks to complete, so the workflow can only be as fast as the slowest parallel task.

It is also important to **always consider biology** when attempting to parallelise a task. In our variant calling example above, we could safely split our data up by chromosome because the chromosomes are physically distinct from one another, and sequencing reads on one chromosome don't need to be considered when looking for SNVs on another chromosome. But now suppose you want to perform structural variant calling. Structural variants include things like translocation events, where a sequence from one chromosome has moved to another chromosome. In this case, you can't simply analyse the chromosomes independently of one another, since information from one chromosome may be linked to information from another. Instead, we need to consider the entire genome as a whole, and as such, the task is much more difficult to parallelise.

![Parallel vs non-parallel tasks](/docs/assets/parallel_non_parallel.png)