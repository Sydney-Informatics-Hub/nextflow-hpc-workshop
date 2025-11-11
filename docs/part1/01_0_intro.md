# 1.0 Part 1 Introduction 

In the first part of this workshop, we will familarise ourselves with some foundational concepts required to effectively run bioinformatics workflows on HPC clusters. We will then apply these concepts to the configuration of a popular nf-core pipeline, Sarek. In part 2, we will further apply these concepts to a custom Nextflow workflow. 

!!! warning "Note the code blocks!"

    You’ll notice that many code examples in this workshop are presented in **tabs** for the different HPC systems we're using in this workshop.

    The commands provided are largely the same, but each system has its own scheduler, queues, and job submission syntax.

    Select the tab that matches the system you’re using — the content will stay synced across the page.  
    
    If you switch to **Setonix** once, all other tabs on this page will automatically follow.


    === "Gadi (PBS)"
        If you've been assigned a training account on Gadi, your commands and configurations will be specific to the PBS Pro job scheduler running on Gadi.

        ```bash
        An example command to run on Gadi
        ```

    === "Setonix (Slurm)"
        If you've been assigned a training account on Setonix, your commands and configurations will be specific to the SLURM job scheduler running on Setonix.

        ```bash
        An example command to run on Setonix
        ```

## 1.0.1 Log in to your assigned HPC

If you haven't already done so, follow the [setup instructions](../setup.md) to set up VSCode and log in to your assigned HPC with the user account and password provided to you.

For this workshop, we will be working within the scratch storage system of the HPCs. Navigate to the scratch space for the workshop project.

1. In the left-hand side bar, click on the "Explorer" tab (an icon that looks like two sheets of paper).

    ![](../img/vscode_explorer.png)

2. Click on "Open Folder"

3. In the text box that appears, enter the path of your assigned directory in the HPC scratch space:

    === "Gadi (PBS)"

        On Gadi, we will be working entirely within the scratch space for the `vp91` project: `/scratch/vp91/`. In this directory, everyone will have their own folder, labelled with the same name as their user ID. This is the path you should enter into the text box when prompted. For example, for the username `usr123`, you would enter:

        ```
        /scratch/vp91/usr123
        ```

    === "Setonix (Slurm)"

        On Setonix, we will be working entirely within the scratch space for the `courses01` project: `/scratch/courses01/`. In this directory, everyone will have their own folder, labelled with the same name as their user ID. This is the path you should enter into the text box when prompted. For example, for the username `usr123`, you would enter:

        ```
        /scratch/courses01/usr123
        ```

## 1.0.2 Setup the project space 

When you first log in, your directory in the scratch space will be an empty folder. The first job for the day will be to clone the workshop materials into this space. To do this, open the VSCode terminal (`Ctrl + J`) and run the following commands:

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

Navigate back up and into the `part1/` directory:

```bash
cd ../../part1/
ls
```

```console title="Output"
config/    README.md scripts/
```

=== "Gadi (PBS)"

    Within the `scripts/` directory is an executable file called `pull_sarek.pbs.sh`. Go ahead and run this from the current directory:

    ```bash
    ./scripts/pull_sarek.pbs.sh
    ```

=== "Setonix (Slurm)"

    Within the `scripts/` directory is an executable file called `pull_sarek.slurm.sh`. Go ahead and run this from the current directory:

    ```bash
    ./scripts/pull_sarek.slurm.sh
    ```

This script will pull the `nf-core/sarek` code from GitHub that we will use in the second half of today's session.

Once the `sarek` pipeline has been pulled from GitHub, there is one more script to run.

=== "Gadi (PBS)"

    The `scripts/` directory contains another file called `setup_singularity.pbs.sh`. Go ahead and run this as well:

    ```bash
    ./scripts/setup_singularity.pbs.sh
    ```

=== "Setonix (Slurm)"

    The `scripts/` directory contains another file called `setup_singularity.slurm.sh`. Go ahead and run this as well:

    ```bash
    ./scripts/setup_singularity.slurm.sh
    ```

This script will copy all of the singularity files required by our exercises today and tomorrow.

Once the scripts have completed, our working directory is fully set up.

As a final step, go to VSCode's "File" menu and select "Open Folder...". Enter the full path to the `part1` directory in the text box:

=== "Gadi (PBS)"

    `/scratch/vp91/<username>/nextflow-on-hpc-materials/part1`

=== "Setonix (Slurm)"

    `/scratch/courses01/<username>/nextflow-on-hpc-materials/part1`

Press Enter and the window will refresh. Now the window is loaded at today's working directory, since everything we will be doing in part 1 of this workshop will be within this folder.

We are now ready to get started with the workshop!
