# 1.0 Part 1 Introduction 

In the first part of this workshop, we will familarise ourselves with some foundational concepts required to effectively run bioinformatics workflows on HPC clusters. We will first apply these concepts to the configuration of a popular nf-core pipeline, Sarek. In part 2, we will further apply these concepts to a custom Nextflow workflow. 

!!! warning "Note the code blocks!"

    You’ll notice that many code examples in this workshop are presented in **tabs** for the different HPC systems we're using in this workshop.

    The commands provided are largely the same, but each system has its own scheduler, queues, and job submission syntax.

    Select the tab that matches the system you’re using — the content will stay synced across the page.  
    
    If you switch to **Setonix** once, all other tabs on this page will automatically follow.


    === "Gadi"
        If you've been assigned a training account on Gadi, your commands and configurations will be specific to the PBS Pro job scheduler running on Gadi.

        ```bash
        An example command to run on Gadi
        ```

    === "Setonix"
        If you've been assigned a training account on Setonix, your commands and configurations will be specific to the SLURM job scheduler running on Setonix.

        ```bash
        An example command to run on Setonix
        ```

## 1.0.1 Log in to your assigned HPC

Log in to your assigned HPC with the user account and password provided to you:

=== "Gadi"

    ```bash
    ssh username@gadi.nci.org.au
    ```

=== "Setonix"

    ```bash
    ssh username@setonix.pawsey.org.au
    ```

For this workshop, we will be working within the scratch storage system of the HPCs. Navigate to the scratch space for the workshop project. The project ID will be provided to you on the day of the workshop.

=== "Gadi"

    NCI projects have randomly-assigned IDs with a two letter and two digit pattern, e.g. `ab01`. The scratch space for the project can be found at `/scratch/<PROJECT_ID>`:

    ```bash
    cd /scratch/ab01
    ```

=== "Setonix"

    Pawsey projects have IDs of the form `pawsey1234`. The scratch space for the project can be found at `/scratch/<PROJECT_ID>`:

    ```bash
    cd /scratch/pawsey1234
    ```

Within the scratch space, you will find a folder with your user name. Navigate into this folder:

=== "Gadi"

    ```bash
    cd /scratch/ab01/username
    ```

=== "Setonix"

    ```bash
    cd /scratch/pawsey1234/username
    ```

In here, you will find several folders containing pre-loaded materials for this workshop:

TODO change this to match new structure that accommodates materials and demo workflow repos. 

```bash
ls
```

```console title="Output"
01.hpc_fundamentals
02.nf-core
03.diy_workflow
```

## 1.0.2 Setup the project space 

For this section, we will be working in `01.hpc_fundamentals/`. Navigate to this directory and inspect its contents; you will find some pre-loaded scripts and data for this section:

```bash
cd 01.hpc_fundamentals
ls
```

```console title="Output"
TODO outputs that reflect directory set up
```