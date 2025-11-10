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
    pipelines.**

    The main Nextflow job is a low-resource and long running job that schedules
    the individual tasks to be run on a compute node. Running this on the login
    node will result in your pipeline being killed prematurely.

    TODO: WHY they get killed (e.g. max 30 mins, 4 GB). Throttling is standard
    on any HPC

    When testing and developing your own pipeline, it is recommended to refer
    back to the HPC's recommended way of running long-running, low-memory head
    jobs. These could include persistent sessions with tools like `tmux` or
    `screen` or dedicated workflow nodes.

    TODO: include links to these

    Alternatively, interactive jobs are a useful option to run and debug your
    workflows interactively. If you need to schedule your Nextflow job, **the
    scheduling options should be included in the `run.sh` script**.

    The last consideration should be whether the compute nodes have network
    access or not. When using containers, it is like that these need to be pulled
    when the process runs on a compute node. If compute nodes do not have network
    access, this will fail to pull the container, and consequently the pipeline.

We start with running things on a single sample. This should be representative
of all the data we run i.e. we will be processing chromosomes 20, 21, 22.

Then, continue to optimise and configure the pipeline with a single sample so
that you can conserve SUs as you develop and benchmark.

Run the pipeline with bare bones configuration. Specify the correct `profile`
for your system.

!!! example "Exercises"

    First load the nextflow and singularity modules, following the same method we learnt yesterday:

    === "Gadi"
        ```bash
        module load nextflow/24.04.5 singularity
        ```
    === "Setonix"
        ```bash
        module load nextflow/24.10.0 singularity/4.1.0-slurm
        ```

    Then execute your nextflow command:

    === "Gadi (PBS)"

        ```bash
        nextflow run main.nf -profile pbspro
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
        nextflow run main.nf -profile slurm
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

!!! example "Exercises"

    Find the `GENOTYPE` job id from `.nextflow.log`

    ```bash
    grep GENOTYPE .nextflow.log
    ```

    Then, view the resources used vs. allocated:

    === "Gadi (PBS)"

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

Note that we have not configured the resources, and default were used for that
queue, or partition:

|         | CPU   | MEM    |
| ------- | ----- | ------ |
| Gadi    | 1     | 512 MB |
| Setonix | 2 x 1 | 1.6 GB |

There are pros (you don't get allocated too many resources) and cons (the job
fails). This both indicates that we need to be intentional how we need to be
explicit with configuration, and being aware of differences to note when
running on different systems.

- With "real" larger data sets, this will likely fail

TODO brief explanation on "2 x 1" for Setonix - virtual cores for the partition

We know that 2GB memory will be sufficient for this process. Let's explicitly
configure that.

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

        nextflow run main.nf -profile pbspro -c conf/custom.config
        ```

    === "Setonix (Slurm)"

        ```groovy title="run.sh"
        #!/bin/bash

        module load nextflow/24.10.0
        module load singularity/4.1.0-slurm

        nextflow run main.nf -profile slurm -c conf/custom.config
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
