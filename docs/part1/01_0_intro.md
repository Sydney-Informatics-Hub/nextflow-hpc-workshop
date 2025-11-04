# 1.0 Part 1 Introduction 

In the first part of this workshop, we will familarise ourselves with some foundational concepts required to effectively run bioinformatics workflows on HPC clusters. We will then apply these concepts to the configuration of a popular nf-core pipeline, Sarek. In part 2, we will further apply these concepts to a custom Nextflow workflow. 

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

!!! note

    Be sure substitute your assigned user name for `username` in the above code example.

For this workshop, we will be working within the scratch storage system of the HPCs. Navigate to the scratch space for the workshop project:

=== "Gadi"

    NCI projects have randomly-assigned IDs with a two letter and two digit pattern. For the purposes of this training session, everyone using Gadi will be a member of the specially-created training project `vp91`. The scratch space for the project can be found at `/scratch/vp91`:

    ```bash
    cd /scratch/vp91
    ```

=== "Setonix"

    For the purposes of this training session, everyone using Setonix will be a member of the specially-created Pawsey training project `courses01`. The scratch space for the project can be found at `/scratch/courses01`:

    ```bash
    cd /scratch/courses01
    ```

Within the scratch space, you will find a folder with your user name. Navigate into this folder:

=== "Gadi"

    ```bash
    cd /scratch/vp91/$USER
    ```

=== "Setonix"

    ```bash
    cd /scratch/courses01/$USER
    ```

## 1.0.2 Setup the project space 

When you first log in, your directory in the scratch space will be an empty folder. The first job for the day will be to clone the workshop materials into this space. To do this, run the following commands:

```bash
git clone https://github.com/Sydney-Informatics-Hub/nextflow-on-hpc-materials.git
cd nextflow-on-hpc-materials
ls
```

You should see a few folders and files inside here:

```console title="Output"
data/      part1/     part2/     README.md
```

Next, look at the contents of the `data/` folder. You will see a few sub-folders containing the data we will be using for this workshop:

```bash
ls data/
```

```console title="Output"
bams/ fqs/  ref/
```

If you navigate into the `ref/` subfolder, you will see a file called `Hg38.subsetchr20-22.tar.gz`. This is a compressed archive of our reference genome data that we will need for this workshop:

```bash
cd data/ref/
ls
```

```console title="Output"
Hg38.subsetchr20-22.tar.gz
```

Run the following command to extract the reference data:

```bash
tar -xzf Hg38.subsetchr20-22.tar.gz
```

After a few seconds, the command should complete, and if you list the directory contents again, you will see the full reference data present in the folder:

```bash
ls
```

```console title="Output"
Hg38.subsetchr20-22.dict      Hg38.subsetchr20-22.fasta.ann Hg38.subsetchr20-22.fasta.pac
Hg38.subsetchr20-22.fasta     Hg38.subsetchr20-22.fasta.bwt Hg38.subsetchr20-22.fasta.sa
Hg38.subsetchr20-22.fasta.amb Hg38.subsetchr20-22.fasta.fai Hg38.subsetchr20-22.tar.gz
```

For the part of this setup procedure, navigate back up and into the `part1/` directory:

```bash
cd ../../part1/
ls
```

```console title="Output"
config/    README.md scripts/
```

=== "Gadi"

    Within the `scripts/` directory is an executable file called `pull_sarek.pbs.sh`. Go ahead and run this from the current directory:

    ```bash
    ./scripts/pull_sarek.pbs.sh
    ```

=== "Setonix"

    Within the `scripts/` directory is an executable file called `pull_sarek.slurm.sh`. Go ahead and run this from the current directory:

    ```bash
    ./scripts/pull_sarek.slurm.sh
    ```

This script will pull the `nf-core/sarek` code from GitHub that we will use in the second half of today's session, and will also pull all of the relevant files required to run it.

Once the script completes, we are fully set up and ready to get started with the workshop!
