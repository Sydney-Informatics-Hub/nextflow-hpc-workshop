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

Most HPC systems - including Gadi and Setonix - also support **containers**. A container is an extremely useful concept, where a piece of software is packaged up into a self-contained environment, including all the dependencies it requires. This container can then be run on any system that supports them. This allows us to run virtually any software package we want, without it needing to be pre-installed on our system. Because containers let us ensure that we are running the exact same version of a package on every system we run it on, they are the **recommended way to deal with software when writing reproducible workflows in Nextflow**. Later on in today's workshop we will see how we configure Nextflow to use containers.

There are vast libraries of pre-made containers available for all sorts of tasks, especially for bioinformatics, which means you usually won't have to write your own containers. For bioinformatics, we generally recommend to use [Biocontainers](https://biocontainers.pro/registry) wherever possible. Biocontainers are pre-built and tested containers specifically for bioinformatics tools, with a huge library and great community support. 

You can find Biocontainers at the following repositories:  

- [Biocontainers registry](https://biocontainers.pro/registry)
- [Quay.io](https://quay.io/organization/biocontainers)
- [DockerHub](https://hub.docker.com/r/biocontainers/biocontainers)
- [Seqera containers](https://seqera.io/containers/)

So, how can we use these containers? There are a few different "flavours" of containers out there. You may have heard of Docker before, which is a very popular container system. However, HPCs typically don't use Docker, as it has some security concerns for shared systems. Instead, another quite popular container system is used on HPCs like Gadi and Setonix, called **Singularity**.

### 1.4.2.1 Pulling a Singularity image

To start working with Singularity, we first need to load the Singularity module. This will differ slightly between different HPCs:

=== "Gadi"

    ```bash
    module load singularity
    ```

=== "Setonix"

    ```bash
    module load singularity/4.1.0-slurm
    ```

Next, we need a Singularity **image** - a file containing the entire container environment for the software we are interested in. Let's continue with our FASTQC example and **pull** a FASTQC image from the Biocontainers repository on Quay.io:

```bash
singularity pull fastqc.sif docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0
```

```console title="Output"
INFO:    Converting OCI blobs to SIF format
INFO:    Starting build...
INFO:    Fetching OCI image...
5.4MiB / 5.4MiB [===================================================================================================================================================] 100 % 5.5 MiB/s 0s
276.6MiB / 276.6MiB [===============================================================================================================================================] 100 % 5.5 MiB/s 0s
INFO:    Extracting OCI image...
INFO:    Inserting Singularity configuration...
INFO:    Creating SIF file...
```

Now, you might notice that the command contains the protocol identifier `docker://`. Without confusing things too much, we are actually pulling a Docker image from Quay.io, not a native Singularity image. One of the great things about Singularity is that it can (most of the time) convert Docker images to its own format, meaning that images that are only available as Docker images can usually still be used with Singularity!

When you run the above command, there will be a new file created in your directory called `fastqc.sif`. This is the image file that we will use to run FASTQC.

### 1.4.2.2 Running `fastqc` within a container

So then, how do we actually use this image file and run the `fastqc` command on our data? First, let's try pulling up the FASTQC help information:

```bash
singularity exec fastqc.sif fastqc --help
```

As you can see, using Singularity is fairly simple. We start by calling `singularity exec`, which means we want to **execute** some command within a container. We then provide the path to the singularity image (in this case `fastqc.sif` in the current directory). Finally, we write the actual command we want to run, in this case, `fastqc --help`. Running this, we should see the help page on the terminal:

```console title="Output"

            FastQC - A high throughput sequence QC analysis tool

SYNOPSIS

        fastqc seqfile1 seqfile2 .. seqfileN

    fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam] 
           [-c contaminant file] seqfile1 .. seqfileN

DESCRIPTION

    FastQC reads a set of sequence files and produces from each one a quality
    control report consisting of a number of different modules, each one of 
    which will help to identify a different potential type of problem in your
    data.

...
```

Great! Being able to print out the help page is a good sign that our image works correctly and has a properly-installed version of `fastqc` available.

So then, let's try adding the FASTQC container to our `fastqc.sh` script.

=== "Gadi"

    ```bash title="fastqc.sh" hl_lines="11,19"
    #!/usr/bin/bash
    #PBS -P ab01
    #PBS -N fastqc
    #PBS -q normal
    #PBS -l ncpus=1
    #PBS -l mem=1GB
    #PBS -l walltime=00:10:00
    #PBS -l storage=scratch/ab01
    #PBS -l wd

    module load singularity

    SAMPLE_ID="gut"
    READS_1="data/ggal/${SAMPLE_ID}_1.fq"
    READS_2="data/ggal/${SAMPLE_ID}_2.fq"

    mkdir -p "results/fastqc_${SAMPLE_ID}_logs"

    singularity exec fastqc.sif \
    fastqc \
        --outdir "results/fastqc_${SAMPLE_ID}_logs" \
        --format fastq ${READS_1} ${READS_2}
    ```

=== "Setonix"

    ```bash title="fastqc.sh" hl_lines="11,19"
    #!/usr/bin/bash
    #SBATCH --account=pawsey1234
    #SBATCH --job-name=fastqc
    #SBATCH --partition=work
    #SBATCH --nodes=1
    #SBATCH --ntasks=1
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=1GB
    #SBATCH --time=00:10:00

    module load singularity/4.1.0-slurm

    SAMPLE_ID="gut"
    READS_1="data/ggal/${SAMPLE_ID}_1.fq"
    READS_2="data/ggal/${SAMPLE_ID}_2.fq"

    mkdir -p "results/fastqc_${SAMPLE_ID}_logs"

    singularity exec fastqc.sif \
    fastqc \
        --outdir "results/fastqc_${SAMPLE_ID}_logs" \
        --format fastq ${READS_1} ${READS_2}
    ```

All we have done is remove the `fastqc` module and instead loaded the `singularity` module, then added `singularity exec fastqc.sif` to the start of our `fastqc` command (on a preceding line followed by a backslash `\`, indicating that the command continues on the next line).

And that's it! Go ahead and try submitting the job again, and after a few minute it should finish successfully!

=== "Gadi"

    ```bash
    qsub fastqc.sh
    ```

=== "Setonix"

    ```bash
    sbatch fastqc.sh
    ```
