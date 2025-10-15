# 1.3 Parallelisation on an HPC

HPCs are clusters of numerous nodes, each with a high number of CPU cores available for use. This makes them excellent environemnts for parallel workflows; that is, workflows with many independent tasks that can be run simultaneously, rather than sequentially. This can considerably speed up large workflows that may otherwise take hours or even days to complete on a standard computer.

A very common example in bioinformatics is analysing sequencing reads that map to different chromosomes in the genome. For example, suppose you want to count the number of RNA sequencing reads that map to every gene in the genome. Since chromosomes are physically separate from one another and there is no overlap between genes on separate chromosomes, you can **parallelise** this task by initially splitting your reads into one file for each chromosome, then performing the counting task on each file separately and simultaneously. Once all the chromosomes have been counted, you can later concatenate the results back together again. This is a common pattern known as **scatter-gather**: you **scatter** your data into independent tasks that can be run simultaneously, then later **gather** the results back together again.

![The scatter-gather pattern](/docs/assets/scatter_gather.png)

