# 1.0 Nextflow on HPC - Part 1

## 1.0.1 Welcome!

Welcome to our Nextflow on HPC workshop! This workshop extends upon our [Nextflow for the Life Sciences](https://sydney-informatics-hub.github.io/hello-nextflow-2025/) workshop and explores running Nextflow pipelines on high performance computing (HPC) systems to enable running intensive bioinformatics workflows in a reproducible and scalable manner.

## 1.0.2 Workshop overview

This workshop is split into two parts. In Part 1 of this workshop, we start by introducing HPC systems, how they work, and how they are used to run resource-intensive jobs in an efficient, parallel manner. We will then explore how Nextflow pipelines can be configured, optimised, and run on these systems. For this purpose, we will also introduce `nf-core`, a community-driven collection of Nextflow pipelines for common bioinformatics tasks.

In [Part 2](/docs/part2/00_intro.md), we will build up a new custom Nextflow pipeline and explore in more depth how we can configure it for running on an HPC, as well as some common workflow patterns that help us process data in parallel on HPCs.

This workshop introduces two specific Austraian HPC systems: Gadi, which is housed at the National Computational Infrastructure (NCI); and Setonix, which is hosted by Pawsey. Through these two HPC systems we will explore the similarities and differences between two of the major HPC schedulers: PBSPro and SLURM, respectively.

## 1.0.3 Why run Nextflow on an HPC?

We saw in [Nextflow for the Life Sciences](https://sydney-informatics-hub.github.io/hello-nextflow-2025/) the benefits of Nextflow for creating reproducible bioinformatics workflows. What we didn't touch on, however, was **where** we should ideally run these pipelines. In that workshop, we created a very simple RNA sequencing pipeline that ran very quickly on a tiny dataset. But in the real world, we often deal with huge datasets that can quickly overwhelm laptop and desktop computers. In these cases, one common solution is to use high performance computing (HPC) systems: clusters of specialised computers that have lots of memory, storage, and CPU power. These systems are particularly useful in cases where the dataset can be split into multiple jobs that can be run in parallel - as is often the case with biological data.

As we will see in this workshop, Nextflow has a built-in ability to run on several common HPC systems, and due to its flexibility and configurability, can be finely tuned to scale up to large biological datasets.
