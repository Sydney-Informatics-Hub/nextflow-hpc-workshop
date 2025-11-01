# 1.2 Software installation is different on HPC

As soon as we start using HPC, we’re working in a shared environment where we can’t control everything and that includes software installation.  

!!! warning "No sudo for you!" 
    Unlike your laptop, you do not have administrative (`sudo`) privileges on HPC systems. On a laptop, you can install software however you like. On HPC, thousands of users share the same system, so unrestricted installs would break environments, cause version conflicts, and introduce security risks. That’s why HPC systems block `sudo`.

Bioinformatics workflows need software, and often many versions of it. Consider the tools we want to use in our mapping and variant calling workflow:

- fastqc/TODO version for quality control
- bwa/TODO version for alignment 
- samtools/TODO version for alignment post-processing
- gatk/TODO version for variant calling 
- multiqc/TODO version for report aggregation 

These tools don’t always play nicely together. You can install and manage software for your own workflows, you just need to use HPC-approved methods. There are three main approaches:

| Approach | Description | Pros | Cons |
|----------|-------------|------|------|
| **Environment modules** | Pre-installed software provided by HPC admins, loaded with `module load` | Fast, easy to use, no setup required | Limited versions; may conflict with workflow needs |
| **Conda/Mamba environments** | User-managed Python/R environments | Flexible, easy for development | Slow installs; dependency conflicts common; not fully reproducible |
| **Containers** (Apptainer/Singularity) | Portable, isolated software environments | Reproducible, portable, avoids dependency issues | Requires container knowledge |

!!! note "Why we love containers" 

    Containers bundle all the software a workflow needs, including tools, dependencies, libraries, OS layers, into a single portable image. On HPC, that means:

    - No dependency conflicts: every container runs in its own isolated environment, unaffected by other users or system modules  
    - Reproducibility: the same container image can be used across clusters, clouds, or laptops, ensuring identical software behaviour everywhere  
    - Reduced maintenace: you don't need to worry about installing, updating, and debugging complex software stacks

## 1.2.1 A simple script 

Let's see this in practice by working with a simple script that runs fastqc: 

```bash title="fastqc.sh"
#!/bin/bash
SAMPLE_ID="tiny"
READS_1="data/${SAMPLE_ID}.R1.fq"
READS_2="data/${SAMPLE_ID}.R2.fq"

mkdir -p "results/fastqc_${SAMPLE_ID}_logs"
fastqc \
    --outdir "results/fastqc_${SAMPLE_ID}_logs" \
    --format fastq ${READS_1} ${READS_2}
```

A small fastq dataset has also been provided:

```bash
TODO update the path according to updated directory structure
ls data
```

```console
TODO update this to work with updated dataset
tiny.R1.fq  tiny.R2.fq
```

!!! example "Try to run fastqc" 
    Try running the script:

    ```bash
    chmod +x fastqc.sh
    ./fastqc.sh
    ```

??? question "Result..."
    You’ll see:
    ```console
    ./fastqc.sh: line 9: fastqc: command not found
    ```
    The command fails because `fastqc` is **not installed by default**.  
    On HPC, you’ll need to either load a pre-installed module or use a container.

## 1.2.2 Using modules

Modules are the standard way to access centrally installed software on HPC systems. There are lots of modules pre-installed on HPCs like Gadi and Setonix that allow you to use many common tools, including lots of bioinformatics tools. `fastqc` is one such tool that is pre-installed on both of these HPC systems.

!!! example "Exercise 1.2.2.1: Finding `fastqc` module"

    You can list available modules on your system with:

    ```bash
    module avail
    ```

    The terminal will fill up with a long list of available modules:
   
    === "Gadi"

        ```console title="Available modules"
        -------------------------------------------------------------------- /opt/Modules/modulefiles ---------------------------------------------------------------------
        pbs  singularity  

        ----------------------------------------------------------------- /opt/Modules/v4.3.0/modulefiles -----------------------------------------------------------------
        dot  module-git  module-info  modules  null  use.own  

        -------------------------------------------------------------------- /apps/Modules/modulefiles --------------------------------------------------------------------
        abaqus/2020                gamess/2021-R2-p1              intel-dnnl/2025.2.0              lammps/3Aug2022           openquake/3.10.1                              
        abaqus/2021                gamess/2022-R2                 intel-dpct/2021.1.1              lammps/3Mar2020           openquake/3.11.2                              
        ```

    === "Setonix"

        ```console title="Available modules"
        ----------------------------------------------- /opt/cray/pe/lmod/modulefiles/mpi/gnu/12.0/ofi/1.0/cray-mpich/8.0 ------------------------------------------------
        cray-hdf5-parallel/1.14.3.1        cray-mpixlate/1.0.5        cray-parallel-netcdf/1.12.3.13 (D)
        cray-hdf5-parallel/1.14.3.3        cray-mpixlate/1.0.6        cray-parallel-netcdf/1.12.3.15
        cray-hdf5-parallel/1.14.3.5 (D)    cray-mpixlate/1.0.7 (D)    cray-parallel-netcdf/1.12.3.17

        ---------------------------------------------- /software/setonix/2025.08/modules/zen3/gcc/14.2.0/astro-applications ----------------------------------------------
        apr-util/1.6.3                  casacore/3.4.0-adios2        giant-squid/2.3.0           mwalib/1.8.7                wcstools/3.9.7
        apr/1.7.5                       casacore/3.4.0-openmp        hyperbeam/0.10.2-cpu        pgplot/5.2.2                wsclean/2.9-idg
        ```

Note that on Setonix some modules are marked with a `(D)`: this indicates that it is the default version of that module. However, Setonix disables automatic loading of default modules; instead, you must explicitly specify the version you want when loading a module.

This list can be navigated with the `j` and `k` keys to move down and up, respectively. To exit the list, press `q`.

Note that the modules appear in a format of `<TOOL NAME>/<VERSION>`. Often, several versions of a tools will be pre-installed on the system for you to choose between.

!!! example "Exercise 1.2.2.2: Loading the `fastqc` module"

    Find the available `fastqc` versions on your system with:


    ```bash
    module avail fastqc
    ```

    Load the module on your system and confirm the version: 
    === "Gadi"

        ```bash
        module load fastqc
        fastqc --version
        ```
        ```console title="Output"
        FastQC v0.12.1
        ```

    === "Setonix"

        ```bash
        module load fastqc/0.11.9--hdfd78af_1
        fastqc --version
        ```
        ```console title="Output: fastqc version"
        FastQC v0.11.9
        ```

    Once loaded, you can rerun the `fastqc.sh` script and verify that it now runs successfully.

    ```bash
    ./fastqc.sh
    ```

    ```console
    TODO update this for updated materials 
    Started analysis of tiny.R1.fq
    Analysis complete for tiny.R1.fq
    Started analysis of tiny.R2.fq
    Analysis complete for tiny.R2.fq
    ```

Modules are quick and convenient, but they depend on what your HPC administrators provide. As we saw above, Gadi and Setonix provided different versions of FastQC. Relying on administrators to install modules for you can be a barrier to running your workflows. 

!!! note "Running your own modules" 
    If you like, you can make your own software "module loadable" for all members of your project space. This is a nice way to share software within a project but can be tedious to maintain as the number of tools and dependencies grow. 

    If you're interested in trying this our for yourself, take a look at [this tutorial](https://hpc.ncsu.edu/Documents/user_modules.php) from NCSU on custom modules.

When you need specific versions or multiple conflicting tools, containers are a better option.

## 1.2.3 Using containers 

Containers are portable software environments, they package everything your tool needs to run. 

!!! example "Exercise 1.2.3.1: Load the Singularity module"

    Unload any existing module first:

    ```bash
    module unload fastqc
    ```

    Singularity, like other software, is not loaded by default. On Gadi, there is just one default Singularity module, and can be simply loaded with `module load`. On Setonix, you can find the versions of Singularity available with `module avail`:
    
    === "Gadi"
  
        ```bash
        module load singularity
        ```

    === "Setonix"

        ```bash
        module avail singularity
        module load singularity/4.1.0-slurm
        ```
        ```console
        -------------------------------------------------- /software/setonix/2025.08/modules/zen3/gcc/14.2.0/utilities ---------------------------------------------------
        singularityce/3.11.4    singularityce/4.1.0 (D)

        ------------------------------------------------------------ /software/setonix/2025.08/pawsey/modules ------------------------------------------------------------
        singularity/3.11.4-askap-gpu    singularity/3.11.4-mpi       singularity/3.11.4-slurm       singularity/4.1.0-mpi-gpu    singularity/4.1.0-nompi
        singularity/3.11.4-askap        singularity/3.11.4-nohost    singularity/4.1.0-askap-gpu    singularity/4.1.0-mpi        singularity/4.1.0-slurm (D)
        singularity/3.11.4-mpi-gpu      singularity/3.11.4-nompi     singularity/4.1.0-askap        singularity/4.1.0-nohost
        ```

For the purposes of today's workshop, we will be using Singularity's default parameters, so we don't need to provide any options. You will also find a pre-loaded Singularity image for running the `fastqc` command in your current directory under the `singularity/fastqc.sif`.

!!! example "Exercise 1.2.3.2: Run FastQC in a Singularity container"

    TODO confirm this - A pre-built FastQC container image has been provided for you at `singularity/fastqc.sif`.  
    
    You can execute your script inside it with:

    ```bash
    singularity exec singularity/fastqc.sif ./fastqc.sh
    ```

    The output should look familiar:
    ```console
    TODO update this with updated materials
    Started analysis of tiny.R1.fq
    Analysis complete for tiny.R1.fq
    Started analysis of tiny.R2.fq
    Analysis complete for tiny.R2.fq
    ```