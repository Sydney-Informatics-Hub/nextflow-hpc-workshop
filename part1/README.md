# Part 1 scripts and exercises

Part 1 is split into two sections. Section 1 introduces participants to HPC environments, in particular NCI's Gadi and Pawsey's Setonix. Gadi uses the PBSPro scheduler, while Setonix uses SLURM. In Section 2, we introduce participants to nf-core, how their pipelines are strucutred, and how to configure them for running on an HPC using layered configuration files and the profiles available at [nf-co.re/configs/](https://nf-co.re/configs/).

## Section 1: Introduction to HPC

Folder: `01.intro_to_hpc/`

### Exercise set 1: Submitting jobs

Folder: `01.submitting`

There are two sets of scripts in this directory. The first is the script `fastqc.sh` along with the HPC-specific submission scripts, `submit.{gadi,setonix}.sh`. The submission scripts are intended to be used within the training materials as examples of code that participants can run directly from the terminal.

The second set of scripts are the self-contained `fastqc.{gadi,setonix}.sh` scripts which participants can submit with a simple `qsub` or `sbatch` command, with all scheduler parameters contained within the script headers. The intention is to first demonstrate the verbose method of submitting scripts, then to introduce the more convenient approach of using the script headers.

### Exercise set 2: Queue limits

Folder: `02.limits`

These scripts are intended to demonstrate the importance of optimising the job parameters. There are three scripts for each scheduler in this directory.

- `fastqc.short_walltime.{gadi,setonix}.sh`: Has a 30 second walltime - demonstrate the need to have a reasonable walltime for the job to complete.
- `fastqc.invalid_queue.{gadi,setonix}.sh`: Has a 32 hour (Setonix) or 72 hour (Gadi) walltime - demonstrate an invalid selection for the queue. Both will fail to submit.
- `fastqc.many_cpus.{gadi,setonix}.sh`: Requests 12 CPUs and specifies 12 threads for the FastQC command - demonstrate a well-provisioned job.

### Exercise set 3: Parallelisation

Folder: `03.parallel`

The intention for this section is to demonstrate the power of HPCs to parallelise jobs. We again run FastQC on multiple samples, but now we submit each sample as a separate job. The scripts `submit.{gadi,setonix}.sh` use a loop to iterate through the samples and submit a job for each one. Note that the Gadi and Setonix jobs are slightly different due to differences in how parameters are passed to the jobs.

## Section 2: Introduction to nf-core

Folder: `02.nf-core`

In this section, we use the `nf-core/rnaseq` pipeline to demonstrate nf-core and how to run these pipelines on an HPC.

### Exercise set 1: Running nf-core pipelines

Folder: `01.run_nf-core`

The script `run.sh` simply runs the `rnaseq` pipeline using singularity but no scheduler.

### Exercise set 2: Configuring nf-core

Folder: `02.config`

Here, `run.sh` has been modified to include a `custom.config` configuration file. This still needs filling out, but we can demonstrate how we can layer a configuration file on top of the existing ones to override some parameters and resource requirements to suit our needs. This could be used to demonstrate poorly-optimised jobs that intentionally fail.

### Exercise set 3: Configuring nf-core for HPC

Folder: `03.hpc_configs`

In this exercise, we have Gadi- and Setonix-specific run scripts (`run.{gadi,setonix}.sh`) that include `gadi` and `setonix` profiles in their `custom.config` files. These include the `gadi.config` and `setonix.config` configuration files for submitting to the clusters.

We also have nf-core run scripts (`run.{gadi,setonix}.nfcore.sh`) that use the public `nci_gadi` and `pawsey_setonix` configurations.