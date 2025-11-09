# 1.2 Software on HPC

!!! info "Learning objectives"

    - Explain software constraints of HPC systems 
    - Demonstrate how to load and use software modules on HPC
    - Describe how containerisation works and understand how containers encapsulate software dependencies 
    - Justify why constainerisation is essential for building portable and reproducible workflows 
    
As soon as we start using HPC, we’re working in a shared environment where we can’t control everything and that includes software installation.  

!!! warning "No sudo for you!"

    Unlike your laptop, you do not have administrative (`sudo`) privileges on HPC systems. This puts some restrictions on what software you can install and where. On a laptop, you can install software however you like. On HPC, thousands of users share the same system, so unrestricted installs would break environments, cause version conflicts, and introduce security risks. That’s why HPC systems block `sudo`.

## 1.2.1 Software installation is different on HPCs

There are a few ways in which we can install and use software on HPCs. Some of these are more complex than others, and they each have their benefits and drawbacks:

| Approach | Description | Pros | Cons |
|----------|-------------|------|------|
| **Pre-built executables** | Pre-built, ready-to-run executable files that you simply download and run | Easy to use - just download and run | Few tools provide pre-built executables on Linux; they must be built for the same operating system and CPU architecture that you are using |
| **DIY installation** | Compiling and installing software yourself from open-source code | Full control over where the software gets installed | Not beginner-friendly; build times can be long; often requires lots of configuration and specific dependencies |
| **Conda/Mamba environments** | User-managed Python/R environments | Flexible; easy for development; well-documented; lots of tutorials available | Slow installs; dependency conflicts common; not fully reproducible; can use up lots of storage space |
| **Environment modules** | Pre-installed software provided by HPC admins, loaded with `module load` | Fast; easy to use; no setup required | Limited versions; may conflict with workflow needs; not every tool is available on every system |
| **Containers** (Apptainer/Singularity) | Portable, isolated software environments | Reproducible; portable; avoids dependency issues; huge repositories of pre-built containers available | Requires container knowledge; introduces slightly more complexity into scripts and workflows |

!!! note "Why we love containers" 

    Of the several options described above, we highly recommend using containers for workflow development, especially when using workflow languages like Nextflow. Containers bundle up everything a tool needs, including, dependencies, libraries, OS layers, into a single portable image. This image is self-contained (hence the name!) and doesn't interact or conflict with any of the software installed on your computer. On HPC, that means:

    - No dependency conflicts: every container runs in its own isolated environment, unaffected by other users or system modules
    - Reproducibility: the same container image can be used across clusters, clouds, or laptops, ensuring identical software behaviour everywhere
    - Reduced maintenace: you don't need to worry about installing, updating, and debugging complex software stacks.

    We source prebuilt containers from the [BioContainers](https://biocontainers.pro/) on the [quay.io registry](https://quay.io/) and [Seqera containers](https://seqera.io/containers/). 

!!! note "One container per tool"

    Containers are an ideal way to package up a tool with all of its dependencies, but are still susceptible to version and dependency conflicts if you try to package up too many tools into the one image. As such, the general best practice is that **one container is built around one tool**. As a consequence, when writing workflows, we also typically want to aim for **one tool and one container per process**. In fact, nexflow only lets us use one container for a single process. This forces us to break up our workflow into smaller, modular chunks, which helps improve readability and maintainability of our pipelines. 

    Note that there are cases where we may package up two very closely related tools into a single container and Nextflow process, either because they are part of a larger suite of software, or they are known to work well together, or it would introduce unnecessary complexity into our pipeline to separate them. However, this is the exception, not the rule.

## 1.2.2 A simple command 

Let's start exploring HPC software with a simple command: accessing the help page for `fastqc`.

!!! example "Exercise: Try to run fastqc"

    Try running the following command:

    ```bash
    fastqc --help
    ```

??? question "Result..."

    You will see:

    ```console
    bash: fastqc: command not found
    ```

    The command fails because `fastqc` is **not available by default**.
    On HPC, you’ll need to either load a pre-installed module or use a container.

## 1.2.3 Using modules

Modules are the standard way to access centrally installed software on HPC systems. There are lots of modules pre-installed on HPCs like Gadi and Setonix that allow you to use many common tools, including lots of bioinformatics tools. `fastqc` is one such tool that is pre-installed on both of these HPC systems.

!!! example "Exercise: Finding the `fastqc` module"

    You can list available modules on your system with:

    ```bash
    module avail
    ```

    The terminal will fill up with a long list of available modules:
   
    === "Gadi (PBS)"

        ```console title="Available modules"
        -------------------------------------------------------------------- /opt/Modules/modulefiles ---------------------------------------------------------------------
        pbs  singularity  

        ----------------------------------------------------------------- /opt/Modules/v4.3.0/modulefiles -----------------------------------------------------------------
        dot  module-git  module-info  modules  null  use.own  

        -------------------------------------------------------------------- /apps/Modules/modulefiles --------------------------------------------------------------------
        abaqus/2020                gamess/2021-R2-p1              intel-dnnl/2025.2.0              lammps/3Aug2022           openquake/3.10.1                              
        abaqus/2021                gamess/2022-R2                 intel-dpct/2021.1.1              lammps/3Mar2020           openquake/3.11.2                              
        ```

    === "Setonix (Slurm)"

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

!!! example "Exercise: Loading the `fastqc` module"

    Find the available `fastqc` versions on your system with:

    ```bash
    module avail fastqc
    ```

    Load the module on your system and confirm the version: 
    === "Gadi (PBS)"

        ```bash
        module load fastqc
        fastqc --version
        ```

        ```console title="Output"
        FastQC v0.12.1
        ```

    === "Setonix (Slurm)"

        ```bash
        module load fastqc/0.11.9--hdfd78af_1
        fastqc --version
        ```

        ```console title="Output: fastqc version"
        FastQC v0.11.9
        ```

    Once loaded, you can rerun the `fastqc --help` command and verify that it now runs successfully.

    ```bash
    fastqc --help
    ```

    ```console
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

Modules are quick and convenient, but they depend on what your HPC administrators provide. As we saw above, Gadi and Setonix provided different versions of FastQC. Relying on administrators to install modules for you can be a barrier to running your workflows. 

!!! note "Running your own modules"

    If you like, you can make your own software "module loadable" for all members of your project space. This is a nice way to share software within a project but can be tedious to maintain as the number of tools and dependencies grow. 

    If you're interested in trying this our for yourself, take a look at [this tutorial](https://hpc.ncsu.edu/Documents/user_modules.php) from NCSU on custom modules.

When you need specific versions or multiple conflicting tools, **containers are a far better option**.

## 1.2.4 Using containers 

Containers are portable software environments: they package everything your tool needs to run. 

!!! example "Exercise: Load the Singularity module"

    Unload any existing module first:

    ```bash
    module unload fastqc
    ```

    Singularity, like other software, is not loaded by default. On Gadi, there is just one default Singularity module, and can be simply loaded with `module load`. On Setonix, you can find the versions of Singularity available with `module avail`:
    
    === "Gadi (PBS)"
  
        ```bash
        module load singularity
        ```

    === "Setonix (Slurm)"

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

To run a command or script inside a singularity container, you simply run `singularity exec /path/to/image.sif <your command>`, where `/path/to/image.sif` is the path to the container image that you wish to use. There are also numerous options and flags that you can provide to the `singularity exec` command itself. For the purposes of today's workshop, we will be using Singularity's default parameters, so we don't need to provide any options. As part of these defaults, the command will be run within your current working directory, as if you simply ran the command directly.

!!! example "Exercise: Run FastQC in a Singularity container"

    You should find a pre-built container image for running the `fastqc` command at `singularity/fastqc.sif`.  
    
    You can execute the `fastqc --help` inside it with:

    ```bash
    singularity exec singularity/fastqc.sif fastqc --help
    ```

    The output should look familiar:

    ```console
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

Later, when we set up Nextflow to run on the HPC, we will configure it to use singularity containers. Behind the scenes, this process of running a command within a container is essentially what Nextflow does for us.
