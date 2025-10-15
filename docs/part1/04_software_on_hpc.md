# 1.4 HPC software

Back in [Section 1.2 - Submitting jobs to an HPC](./02_submitting_jobs.md), we touched on an important issue with HPCs: that installing and using software isn't as straight-forward as it can be on a local computer. There are a few reasons for this. First of all, many software packages require elevated permissions to install. For shared systems like HPCs, it would be a pretty major security risk to allow every user to have such permissions! Additionally, these software packages are usually installed globally, that is in places accessible to all users. If everyone was allowed to install these packages, we could quickly run into compatibility issues and the whole system may break.

Luckily, there are a few ways to address this issue.

First, some - but not all! - programs can be install **locally** - that is, in your home directory or some other location that is only accessible to you. These are typically smaller, less complex programs. The downside to this approach is that you will usually need to build the software from source code yourself; very few packages on Linux - the operating system of choice for HPCs - come as ready-to-use executable files. If you aren't comfortable with building software, this is likely not a good choice for you.

## 1.4.1 The `module` system

Thankfully, most HPCs have lots of pre-installed packages for the most common tasks. To avoid conflicts between packages, these are typically not loaded by default, but are availabe as **modules** which you can load at any time. We have already touch upon this module system in [Section 1.2](./02_submitting_jobs.md) and saw that FASTQC was available on both Gadi and Setonix. In fact, both systems have lots of common bioinformatics tools, such as `bcftools` and `samtools`, `bwa-mem`, and the GATK suite, to name a few. You can load modules with the `module load ...` command, unload them with `module unload ...`, and clean your environment completely with `module purge`. You can also list all available modules with `module avail`:

```bash
module avail
```

```console title="Output"
------------------------------------------------------------------------------- /opt/Modules/modulefiles --------------------------------------------------------------------------------
pbs  singularity  

---------------------------------------------------------------------------- /opt/Modules/v4.3.0/modulefiles ----------------------------------------------------------------------------
dot  module-git  module-info  modules  null  use.own  

------------------------------------------------------------------------------- /apps/Modules/modulefiles -------------------------------------------------------------------------------
abaqus/2020                gamess/2020-R2                 intel-dnnl/2024.2.0              julia/1.9.1               openmpi/4.1.7                                 visrtx/0.5.0     
abaqus/2021                gamess/2021-R1                 intel-dnnl/2024.2.1              julia/1.10.4              openmpi/4.1.7-debug                           vmd/1.9.3        
abyss/2.2.3                gamess/2021-R2-p1              intel-dnnl/2025.0.1              julia/1.11.0              openmpi/5.0.5                                 vtk/8.1.0        

...
```

The module system is great for one-off, common tasks, but you might already see a problem with it: you are dependent upon whatever software the system administrators have pre-installed, which may vary from system to system. A package that is availabe on Gadi might not be available on Setonix, and vice versa. Even if it is, different versions may be installed, which can often result in different results at best, or errors at worst. Our entire goal with using tools like Nextflow is to write reproducible, consistent workflows that work on any system we take them to: that is, "write once, run anywhere". So, how can we accomplish this?

## 1.4.2 Containers

Most HPC systems - including Gadi and Setonix - also support **containers**. A container is an extremely useful concept, where a piece of software is packaged up into a self-contained environment, including all the dependencies it requires. This container can then be run on any system that supports them. This allows us to run virtually any software package we want, without it needing to be pre-installed on our system.

There are vast libraries of pre-made containers available for all sorts of tasks, especially for bioinformatics, which means you usually won't have to write your own containers. For bioinformatics, we generally recommend to use [Biocontainers](https://biocontainers.pro/registry) wherever possible. Biocontainers are pre-built and tested containers specifically for bioinformatics tools, with a huge library and great community support. 

You can find Biocontainers at the following repositories:  

- [Biocontainers registry](https://biocontainers.pro/registry)
- [Quay.io](https://quay.io/organization/biocontainers)
- [DockerHub](https://hub.docker.com/r/biocontainers/biocontainers)
- [Seqera containers](https://seqera.io/containers/)

So, how can we use these containers? There are a few different "flavours" of containers out there. You may have heard of Docker before, which is a very popular container system. However, HPCs typically don't use Docker, as it has some security concerns for shared systems. Instead, another quite popular container system is used on HPCs like Gadi and Setonix, called **Singularity**.

### 1.4.2.1 Pulling a Singularity image

