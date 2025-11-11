# Configuring custom Nextflow pipelines

!!! info "Learning objectives"

    - Identify best practices for benchmarking and configuring pipelines,
    step-by-step
    - Troubleshoot common HPC and scheduling errors by inspecting task logs
    and interpretting exit codes and error messages
    - Know how to find and specify system-dependent HPC values such as queues/partitions
    - Recall how Nextflow interacts with HPC components

!!! warning Running the head job on the correct node

    For the workshop, we have pre-pulled containers and use them from the cache,
    and as data is small and the workflow runs quickly, we will be running them
    on login nodes.

    **However, this should not be done when developing and running your own
    pipelines on HPC systems.**

    When developing or running workflows with real data, you should follow your HPC's
    recommended approach for long-running, low-resource jobs. This may include being
    used in conjunction with terminal multiplexers such as `tmux` or `screen` within
    these sessions to keep the job running when you log out or connection drops out.
    
    Alternatively, starting your Nextflow run within an **interactive job** is useful for testing and debugging 
    workflows directly on compute nodes. If you instead schedule your Nextflow
    job, ensure that your `run.sh` script includes the appropriate scheduler options.

    See the recommendations on running the head job for [Gadi](https://opus.nci.org.au/spaces/Help/pages/241926895/Persistent+Sessions) and [Setonix](https://pawsey.atlassian.net/wiki/spaces/US/pages/286097469/How+to+Run+Workflows+on+the+Workflow+Nodes)

    Finally, be aware of **network access on compute nodes**. If Nextflow cannot locate
    the container, it will attempt to pull from an online repository. This occurs during run time and will fail if the compute node does not have internet access. In these cases, you can pre-pull containers, or schedule the head job on a queue/partition with network access.

Before launching a full-scale analysis, it is important to optimise your pipeline using a small, representative subset of your data. This helps you:

- Estimate the computational requirements of each step
- Configure your pipeline
- Avoid wasting service units (SUs) during development

In this case, we’re starting with the reads from a single individual (NA1287) from the cohort. This is a good proxy for the other two samples as they all contain the same subset of chromosomes (20, 21, and 22).

Once we’re confident everything works as intended, we can scale up to run on the full dataset.

In real-world pipelines, throughput becomes a key consideration of how you choose to configure pipelines. It’s not just about how fast one sample runs, but how many samples can be processed concurrently. A well-optimised configuration ensures that your pipeline makes efficient use of the HPC resources available, reducing queue times, avoiding bottlenecks, and increasing the number of samples processed per unit time.

The first step is to get our starting pipeline working on the HPC. This means running the head job to schedule the processes and executed on compute nodes.

!!! example "Exercises"

    First load the Nextflow and singularity modules, following the same method we learnt yesterday:

    === "Gadi (PBS)"
        ```bash
        module load nextflow/24.04.5 singularity
        ```

    === "Setonix (Slurm)"
        ```bash
        module load nextflow/24.10.0 singularity/4.1.0-slurm
        ```

    Then execute your Nextflow command. **(Note: Your run may fail here - that's ok for this step!)**

    === "Gadi (PBS)"

        ```bash
        nextflow run main.nf -profile pbspro --pbspro_account vp91
        ```

        ??? abstract "Output"

            ```
            executor >  pbspro (5)
            [e4/dc7852] FASTQC (fastqc on NA12889) | 3 of 3 ✔
            [3f/c0e6ae] ALIGN (1)                  | 1 of 1 ✔
            [d0/6677c4] GENOTYPE (1)               | 1 of 1, failed: 1 ✘
            [-        ] JOINT_GENOTYPE             -
            [-        ] STATS                      -
            [-        ] MULTIQC                    | 0 of 1
            ERROR ~ Error executing process > 'GENOTYPE (1)'

            Caused by:
              Process `GENOTYPE (1)` terminated with an error exit status (247)


            Command executed:

              gatk --java-options "-Xmx4g" HaplotypeCaller -R Hg38.subsetchr20-22.fasta -I NA12877.bam -O NA12877.g.vcf.gz -ERC GVCF

            Command exit status:
              247

            Command output:
              (empty)

            Command error:
            ```

    === "Setonix (Slurm)"

        ```bash
        nextflow run main.nf -profile slurm --slurm_account courses01
        ```

        ??? abstract "Output"

            ```console
            N E X T F L O W   ~  version 24.10.0                 14:13:57

            Launching `main.nf` [dreamy_cuvier] DSL2 - revision: 5e5c4f57e0

            executor >  slurm (8)
            [c4/babfde] FASTQC (fastqc on NA12877) | 3 of 3 ✔
            [6f/7a523b] ALIGN (1)                  | 1 of 1 ✔
            [bf/525fc9] GENOTYPE (1)               | 1 of 1 ✔
            [5c/ea7cb0] JOINT_GENOTYPE (1)         | 1 of 1 ✔
            [7f/4e64dc] STATS (1)                  | 1 of 1 ✔
            [70/254c9e] MULTIQC                    | 1 of 1 ✔
            Completed at: 05-Nov-2025 13:59:02
            Duration    : 2m 21s
            CPU hours   : (a few seconds)
            Succeeded   : 8
            ```

The jobs were scheduled to run on the compute nodes successfully, indicated by the `executor > [name]`. However, the pipeline may have failed!

!!! question "Zoom react!"

    1. If your job has finished succesfully, react "Yes" on Zoom, and "No" if it returned an error
    2. Similarly, react "Yes" if you are running it on Gadi, and "No" for Setonix

A common reason that pipelines fail on HPC is due to improper configuration. Here, we have yet to configure the resources properly, to the default allocation from the system's queue or partition were used.

Let's explore what resources were actually used and compare them to what was allocated, by inspecting the logs and system job information using the methods from Part 1.

!!! example "Exercises"

    Find the Job ID of the failed or completed `GENOTYPE` process in your `.nextflow.log`:

    ```bash
    grep GENOTYPE .nextflow.log
    ```

    Then, inspect the job resource usage:

    === "Gadi (PBS)"

        Use `qstat -xf <job_id>` to query the job info and display the top 12 lines (there is no relevant information beyond the 12th line for this exercise).

        ```bash
        qstat -xf <job_id> | head -12
        ```
        ```console
        Job Id: 154038474.gadi-pbs
        Job_Name = nf-GENOTYPE_1
        Job_Owner = fj9712@gadi-login-05.gadi.nci.org.au
        resources_used.cpupercent = 94
        resources_used.cput = 00:02:50
        resources_used.jobfs = 0b
        resources_used.mem = 512124kb
        resources_used.ncpus = 1
        resources_used.vmem = 512124kb
        resources_used.walltime = 00:03:03
        job_state = F
        queue = normal-exec
        ```

    === "Setonix (Slurm)"

        Use `sacct` with formatting to view key stats:

        ```bash
        sacct -j <job_id> --format=JobID,JobName,User,CPUTime,TotalCPU,NCPUS,Elapsed,State,MaxRSS,MaxVMSize,Partition
        ```
        ```console
        JobID           JobName      User    CPUTime   TotalCPU      NCPUS    Elapsed      State     MaxRSS  MaxVMSize  Partition
        ------------ ---------- --------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------
        34325657     nf-GENOTY+     fjaya   00:01:50  01:07.026          2   00:00:55  COMPLETED                             work
        34325657.ba+      batch             00:01:50  01:07.023          2   00:00:55  COMPLETED   1566932K          0
        34325657.ex+     extern             00:01:50  00:00.003          2   00:00:55  COMPLETED          0          0
        ```

In both cases, we can observe that the jobs were assigned the following resources:

|         | CPU   | Memory   |
| ------- | ----- | -------- |
| Gadi    | 1     | 512 MB   |
| Setonix | 2     | 1.6 GB   |

This is why **explicit resource configuration** is important. Even though the pipeline technically ran (or failed), these defaults are unsuitable for real data.

In this case, it shows that the `GENOTYPE` process needs at least 2 GB of memory. Let's explicitly configure that in the next step.

## 2.1.2 Building a custom config

![](figs/00_custom_configs.png)

[TODO] update figure. Add explanation for why you're building out the pipeline like this

!!! tip

    Note that different values are provided based on the specific, low-cost
    queue and partition.

    Here we assign the average number of cores available based on memory
    requirements. This will be revisited in the resourcing section.

!!! example "Exercises"

    1. Create a new file `conf/custom.config`

        ```bash
        touch conf/custom.config
        ```

    2. Add the following contents based on your HPC

    === "Gadi (PBS)"

        ```groovy title='custom.config'
        process {
            cpu = 4 // 'normalbw' queue = 128 GB / 28 CPU ~ 4.6
            memory = 2.GB
        }
        ```

    === "Setonix (Slurm)"

        ```groovy title='custom.config'
        process {
            cpu = 2 // 'work' partition = 230 GB / 128 CPU ~ 1.8
            memory = 2.GB
        }
        ```

Next, add the additional options from part 1.

- Institutional config

!!! example "Exercises"

    === "Gadi (PBS)"

        ```groovy title="conf/pbspro.config"
        singularity {
            enabled = true
            autoMounts = true
            cacheDir = "$projectDir/singularity"
        }

        executor {
            queueSize = 300
            pollInterval = '30 sec'
            queueStatInterval = '30 sec'
        }

        process {
            executor = 'pbspro'
            storage = "scratch/${System.getenv('PROJECT')}"
            module = 'singularity'
            cache = 'lenient'
            stageInMode = 'symlink'
            queue = 'normalbw'
            clusterOptions = "-P ${System.getenv('PROJECT')}"
        }
        ```

    === "Setonix (Slurm)"

        ```groovy title="conf/slurm.config"
        singularity {
            enabled = true
            autoMounts = true
            cacheDir = "$projectDir/singularity"
        }

        executor {
            queueSize = 300
            pollInterval = '30 sec'
            queueStatInterval = '30 sec'
        }

        process {
            executor = 'slurm'
            module = 'singularity/4.1.0-slurm'
            cache = 'lenient'
            stageInMode = 'symlink'
            queue = 'work'
            clusterOptions = "--account=${System.getenv('PAWSEY_PROJECT')}"
        }
        ```

Nextflow is powerful for HPC vs. serially running pbs/slurm scripts.

Next we will wrap this up in a run script.

We will add the new custom configs using `-c`. This will be revisited in
resourcing.

!!! example "Exercises"

    1. Create a new file called `run.sh`

        ```bash
        touch run.sh
        ```

    2. Copy and paste the following code based on your HPC:

    === "Gadi (PBS)"

        ```groovy title="run.sh"
        #!/bin/bash

        module load nextflow/24.04.5
        module load singularity

        nextflow run main.nf -profile pbspro --pbspro_account vp91 -c conf/custom.config
        ```

    === "Setonix (Slurm)"

        ```groovy title="run.sh"
        #!/bin/bash

        module load nextflow/24.10.0
        module load singularity/4.1.0-slurm

        nextflow run main.nf -profile slurm --slurm_account courses01 -c conf/custom.config
        ```

    3. Save the run.sh file (Windows: Ctrl+S, macOS: Cmd+S).
    4. Provide execute permission by running

        ```bash
        chmod +x run.sh
        ```

!!! example "Exercise"

    Run your newly configured pipeline using by executing `./run.sh` in the terminal.
    ??? note "Results"

        On both Gadi and Setonix, both runs should now be successful and
        executed on the respective scheduler.

        === "Gadi (PBS)"

            ```bash
            Loading nextflow/24.04.5
            Loading requirement: java/jdk-17.0.2

             N E X T F L O W   ~  version 24.04.5

            Launching `main.nf` [determined_picasso] DSL2 - revision: 5e5c4f57e0

            executor >  <scheduler> (6)
            [7b/16f7c2] FASTQC (fastqc on NA12877) | 1 of 1 ✔
            [ec/5fd924] ALIGN (1)                  | 1 of 1 ✔
            [17/803fe5] GENOTYPE (1)               | 1 of 1 ✔
            [0f/e21d58] JOINT_GENOTYPE (1)         | 1 of 1 ✔
            [40/0c6e28] STATS (1)                  | 1 of 1 ✔
            [6a/94910b] MULTIQC                    | 1 of 1 ✔
            Completed at: 10-Nov-2025 12:19:38
            Duration    : 4m 1s
            CPU hours   : (a few seconds)
            Succeeded   : 6
            ```

        === "Pawsey (Slurm)"


            ```bash
             N E X T F L O W   ~  version 24.10.0

            Launching `main.nf` [nostalgic_bell] DSL2 - revision: 5e5c4f57e0

            executor >  slurm (6)
            [a4/1eae6a] FASTQC (fastqc on NA12877) [100%] 1 of 1 ✔
            [76/9e6fca] ALIGN (1)                  [100%] 1 of 1 ✔
            [96/f57a2e] GENOTYPE (1)               [100%] 1 of 1 ✔
            [ec/547a9c] JOINT_GENOTYPE (1)         [100%] 1 of 1 ✔
            [88/b49a40] STATS (1)                  [100%] 1 of 1 ✔
            [8d/b5b351] MULTIQC                    [100%] 1 of 1 ✔
            Completed at: 10-Nov-2025 10:28:08
            Duration    : 3m 31s
            CPU hours   : (a few seconds)
            Succeeded   : 6
            ```

qstat/sacct each job is inefficient, especially with pipelines with more
processes, and running on more than one sample. The next section will introduce
how this can be automated using Nextflow's in-built monitoring features.
