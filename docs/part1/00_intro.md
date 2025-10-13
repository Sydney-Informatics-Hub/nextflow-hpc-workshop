# 1.0 Nextflow on HPC - Part 1

## 1.0.1 Welcome!

Welcome to Part 1 of the Nextflow on HPC workshop. This workshop extends upon our [Nextflow for the Life Sciences](https://sydney-informatics-hub.github.io/hello-nextflow-2025/) workshop and explores running Nextflow pipelines on high performance computing (HPC) systems to enable running intensive bioinformatics workflows in a reproducible and scalable manner.

## 1.0.1 Session overview

This session starts with a recap of the RNA sequencing Nextflow pipeline we built in [Nextflow for the Life Sciences](https://sydney-informatics-hub.github.io/hello-nextflow-2025/) and explores the limitations of running Nextflow pipelines locally when working with larger datasets. We then introduce HPC systems, how they work, and how they are used to run resource-intensive jobs in an efficient, parallel manner. We will then explore how Nextflow pipelines can be configured, optimised, and run on these systems. For this purpose, we will also introduce `nf-core`, a community-driven collection of Nextflow pipelines for common bioinformatics tasks. In [Part 2](/docs/part2/00_intro.md), we will return to our custom RNA sequencing Nextflow pipeline, configure it for running on an HPC, and build upon it to explore some common workflow patterns that help us process data in parallel on HPCs.

This workshop introduces two specific Austraian HPC systems: Gadi, which is housed at the National Computational Infrastructure (NCI); and Setonix, which is hosted by Pawsey. Through these two HPC systems we will explore the similarities and differences between two of the major HPC schedulers: PBSPro and SLURM, respectively.
