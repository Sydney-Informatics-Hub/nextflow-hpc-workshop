# Running a custom pipeline on HPC

!!! info "Learning objectives"

    - Execute a custom Nextflow pipeline on HPC, applying system-specific configurations
    - Diagnose and troubleshoot workflow failures by inspecting logs and files in work directories 
    - Construct a reproducible configuration file for an HPC environment.

!!! warning "Testing and developing in the right environment"

    **Containers**
    For the workshop, we have pre-pulled containers and use them from the cache. In real-world scenarios, pulling containers may take time. Consideration to whether your HPC system allows internet access from compute nodes is also important.

    **Running the Nextflow head job**
    As our test data is small and the workflow runs quickly, we will be running the Nextflow `run` command directly on the login nodes. On HPC, this is neither permitted nor feasible for real data analyses, as HPC login nodes typically kill long-running commands. 
    
    We recommend using **interactive jobs** for testing and debugging workflows directly on compute nodes instead of the login nodes. This has the added benefit of testing in the same hardware environment as the final workflow will run. 
    
    Some HPCs have dedicated workflow queues/partitions to which the Nextflow head job can be submitted. These nodes are set up to support long-running jobs with low resource requirements, as workflow processes that require higher resources are each submitted to different queues/partitions based on the resources specified within the custom infrastructure config. 

    See the recommendations on running the head job for [Gadi](https://opus.nci.org.au/spaces/Help/pages/241926895/Persistent+Sessions) and [Setonix](https://pawsey.atlassian.net/wiki/spaces/US/pages/286097469/How+to+Run+Workflows+on+the+Workflow+Nodes).

It is useful to develop your pipelines using a small, representative subset of your data. This allows you:

- Rapidly iterate your workflow design
- Validate environment and software set up
- Tune pipeline configuration without burning service units
- Reduce queue wait times
- Estimate basic performance characteristics of each process

## 2.1.1 Run without configuration

We will start configuring our custom pipeline with using a subset of raw reads from a single individual (NA1287) in our sample cohort. This is a good proxy for the other two samples as they are all of the same input data type and about the same size. Once we’re confident everything works as intended, we will scale up to run on the full set of samples.

Let's start by running the pipeline out of the box to identify what we need to configure: 


!!! example "Exercise: Run out of the box"

    1. Load the Nextflow module, following the same method we learnt yesterday:

    === "Gadi (PBS pro)"
        ```bash
        module load nextflow/24.04.5 
        ```

    === "Setonix (Slurm)"
        ```bash
        module load nextflow/24.10.0
        ```

    2. Run your Nextflow command out of the box:

    === "Gadi (PBS pro)"

        ```bash
        nextflow run main.nf
        ```

        ??? abstract "Output"

            ```
            N E X T F L O W   ~  version 24.04.5

            Launching `main.nf` [suspicious_almeida] DSL2 - revision: 5e5c4f57e0

            executor >  local (2)
            [a3/a6da06] FASTQC (fastqc on NA12877) [100%] 1 of 1, failed: 1 ✘
            [16/860647] ALIGN (1)                  [100%] 1 of 1, failed: 1 ✘
            [-        ] GENOTYPE                   -
            [-        ] JOINT_GENOTYPE             -
            [-        ] STATS                      -
            [-        ] MULTIQC                    -
            ERROR ~ Error executing process > 'FASTQC (fastqc on NA12877)'

            Caused by:
            Process `FASTQC (fastqc on NA12877)` terminated with an error exit status (127)


            Command executed:

            mkdir -p "fastqc_NA12877"
            fastqc -t 1 --outdir "fastqc_NA12877" --format fastq NA12877_chr20-22.R1.fq.gz NA12877_chr20-22.R2.fq.gz

            Command exit status:
            127

            Command output:
            (empty)

            Command error:
            .command.sh: line 3: fastqc: command not found

            Work dir:
            /scratch/project/username/nextflow-on-hpc-materials/part2/work/a3/a6da061dd65d454add3f000923235d

            Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

            -- Check '.nextflow.log' file for details
            ```

    === "Setonix (Slurm)"

        ```bash
        nextflow run main.nf
        ```

        ??? abstract "Output"

            ```console
             N E X T F L O W   ~  version 24.10.0

            Launching `main.nf` [trusting_brown] DSL2 - revision: 5e5c4f57e0

            executor >  local (2)
            executor >  local (2)
            [e6/21a49e] FASTQC (fastqc on NA12877) [100%] 1 of 1, failed: 1 ✘
            [68/bda517] ALIGN (1)                  [100%] 1 of 1, failed: 1 ✘
            [-        ] GENOTYPE                   -
            [-        ] JOINT_GENOTYPE             -
            [-        ] STATS                      -
            [-        ] MULTIQC                    -
            ERROR ~ Error executing process > 'ALIGN (1)'

            Caused by:
            Process `ALIGN (1)` terminated with an error exit status (127)


            Command executed:

            bwa mem -t 1 -R "@RG\tID:NA12877\tPL:ILLUMINA\tPU:NA12877\tSM:NA12877\tLB:NA12877\tCN:SEQ_CENTRE" ref/Hg38.subsetchr20-22.fasta NA12877_chr20-22.R1.fq.gz NA12877_chr20-22.R2.fq.gz | samtools sort -O bam -o NA12877.bam
            samtools index NA12877.bam

            Command exit status:
            127

            Command output:
            (empty)

            Command error:
            .command.sh: line 2: samtools: command not found
            .command.sh: line 2: bwa: command not found

            Work dir:
            /scratch/project/username/nextflow-on-hpc-materials/part2/work/68/bda517eeaa47f5ba24a9c401bdeb0e

            Tip: view the complete command output by changing to the process work dir and entering the command `cat .command.out`

            -- Check '.nextflow.log' file for details
            ```

!!! question "Why did it fail?"
    Use the output printed to your screen, `.nextflow.log` file in your current directory, as well as the `.command.log`, `.command.err`, and `.exitcode ` files in the work directory of the failed task to identify what caused your workflow run to fail. 

    ??? abstract "Answer"
        Our run stopped because one or more processes failed. [Exit status `127` in Linux environments](https://linuxconfig.org/how-to-fix-bash-127-error-return-code) means your system was not able to find the command referenced in the process. This suggests the software is not available in our environment. 

## 2.1.2 Enabling containers

All our process modules specify a container to run inside. This can only happen if Singularity is explicitly enabled in our configuration. Let's enable this in our system-specific configuration files and attempt to run again:  

!!! example "Exercise: Enable Singularity"

    1. Load the Singularity module, following the same method we learnt yesterday:

    === "Gadi (PBS pro)"
        ```bash
        module load singularity
        ```

    === "Setonix (Slurm)"
        ```bash
        module load singularity/4.1.0-slurm
        ```

    2. Add the following to your system-specific config file that you can find in `config/`. Remember, we have already enabled profiles in our `nextflow.config`, so no need to edit that file. 

    === "Gadi (PBS pro)"
        ```groovy title="config/pbspro.sh"
        process {
            // Load the globally installed singularity module before running any process
            module = 'singularity'
        }

        singularity {
            // Explicitly turns on container execution
            enabled = true
            // Automatically bind-mount working directory on scratch and common system paths
            autoMounts = true
            // Define location of stored container images 
            cacheDir = "/scratch/${System.getenv('PROJECT')}/${System.getenv('USER')}/nextflow-on-hpc-materials/singularity"
        }
        ```

    === "Setonix (Slurm)"
        ```groovy title="config/slurm.config"
        process {
            // Load the globally installed singularity/4.1.0-slurm module before running any process
            module = 'singularity/4.1.0-slurm'
        }

        singularity {
            // Explicitly turns on container execution
            enabled = true
            // Automatically bind-mount working directory on scratch and common system paths
            autoMounts = true
            // Define location of stored container images 
            cacheDir = "/scratch/${System.getenv('PAWSEY_PROJECT')}/${System.getenv('USER')}/nextflow-on-hpc-materials/singularity"
        }
        ```

    3. Run your updated Nextflow command:

    === "Gadi (PBS pro)"

        ```bash
        nextflow run main.nf -profile pbspro
        ```

        ??? abstract "Output"

            ```console title="Output"
             N E X T F L O W   ~  version 24.04.5

            Launching `main.nf` [big_picasso] DSL2 - revision: e34a5e5f9d

            executor >  local (6)
            [c3/dda9c5] FASTQC (fastqc on NA12877) | 1 of 1 ✔
            [b4/210425] ALIGN (1)                  | 1 of 1 ✔
            [89/5b15d9] GENOTYPE (1)               | 1 of 1 ✔
            [7d/a18133] JOINT_GENOTYPE (1)         | 1 of 1 ✔
            [3f/1db7ee] STATS (1)                  | 1 of 1 ✔
            [be/fd841e] MULTIQC                    | 1 of 1 ✔
            ```

    === "Setonix (Slurm)"

        ```bash
        nextflow run main.nf -profile slurm
        ```

        ??? abstract "Output"

            ```console
             N E X T F L O W   ~  version 24.10.0

            Launching `main.nf` [distracted_angela] DSL2 - revision: e34a5e5f9d

            executor >  local (6)
            [bb/ddf5ff] FASTQC (fastqc on NA12877) | 1 of 1 ✔
            [1f/b0906a] ALIGN (1)                  | 1 of 1 ✔
            [34/f0234b] GENOTYPE (1)               | 1 of 1 ✔
            [07/08074e] JOINT_GENOTYPE (1)         | 1 of 1 ✔
            [e0/b78049] STATS (1)                  | 1 of 1 ✔
            [27/35c181] MULTIQC                    | 1 of 1 ✔
            ```

Your workflow should have run successfully, however, there is one grave mistake when running on the HPC - **all processes were run on the login node!** This is Nextflow's default behaviour when no executor is specified. 

## 2.1.3 Scheduling jobs

Recall from Part 1 that Nextflow's executor is the part of the workflow engine that talks to the computing environment (whether it's a laptop or HPC). When running on a shared HPC system, these settings are important to include so they are **queued properly on a compute node**.

To have processes run on the compute nodes, `executor` needs to be set to the appropriate job scheduler in our system-specific config files to avoid the default executor (`local`) being used. 'Local' executions means to run in the same environment as the head job, whether that be your local laptop or the login node of an HPC. 

Although Gadi and Setonix both run Nextflow with Singularity, each has different environment variables, filesystem layouts, job schedulers, queue structures, module names/versions, and container cache behaviour. These differences affect how Nextflow executes each process. To have Nextflow submit our processes as separate compute jobs, we need to instruct which **job scheduler** and **queue/partition** should be used. 

The job scheduler is specified with the Nextflow configuration option `executor` while the queue/partition to submit jobs to is specified by the `queue` option. 

!!! note

    In Part 2 we will use the same queues and partitions:

    - The `normalbw` queue on Gadi
    - The `work` partition on Pawsey

    These are low-cost queues suitable for general compute jobs.

Let's add both the `executor` and `queue` configuration options to our system-specific configs to tell Nextflow **where** to run the processes. 

!!! example "Exercise: Configure executor and queue"

    === "Gadi (PBS pro)"

        ```groovy title="config/pbspro.config" hl_lines="4-6"
        process {
            // Load the globally installed singularity module before running any process
            module = 'singularity'
            // Run using the pbspro scheduler on the 'normalbw' queue
            executor = 'pbspro'
            queue = 'normalbw'
        }

        singularity {
            // Explicitly turns on container execution
            enabled = true
            // Automatically bind-mount working directory on scratch and common system paths
            autoMounts = true
            // Define location of stored container images 
            cacheDir = "/scratch/${System.getenv('PROJECT')}/${System.getenv('USER')}/nextflow-on-hpc-materials/singularity"
        }
        ```

    === "Setonix (Slurm)"

        ```groovy title="config/slurm.config" hl_lines="4-6"
        process {
            // Load the globally installed singularity/4.1.0-slurm module before running any process
            module = 'singularity/4.1.0-slurm'
            // Run using the pbspro scheduler on the 'normalbw' queue
            executor = 'slurm'
            queue = 'work'
        }

        singularity {
            // Explicitly turns on container execution
            enabled = true
            // Automatically bind-mount working directory on scratch and common system paths
            autoMounts = true
            // Define location of stored container images 
            cacheDir = "/scratch/${System.getenv('PAWSEY_PROJECT')}/${System.getenv('USER')}/nextflow-on-hpc-materials/singularity"
        }
        ```



Before we re-run our pipeline, we want to add a few more settings to ensure we are using the HPC responsibly. These include options so you don't overwhelm the system by submitting too many jobs or job queries at once (`queueSize`, `pollInterval`, `queueStatInterval`), and options to avoid generating duplicate files (`cache`, `stageInMode`). 

Note that for Setonix, we have specified a `reservation`. This enables us to use reserved resources on Setonix for the duration of the workshop. If you are running Nextflow on Setonix outside of the workshop, this option should be omitted.

!!! example "Exercise: Further config options"

        1. Add the additional config options to your system-specific config file: 

    === "Gadi (PBS pro)"
    
        ```groovy title="config/pbspro.config" hl_lines="1 2 3 4 5 6 14 15 16 17 18"
        executor {
            // For high-throughput jobs, these values should be higher
            queueSize = 30
            pollInterval = '5 sec'
            queueStatInterval = '5 sec'
        }

        process {
            // Load the globally installed singularity module before running any process
            module = 'singularity'
            // Run using the pbspro scheduler on the 'normalbw' queue
            executor = 'pbspro'
            queue = 'normalbw'
            clusterOptions = "-P ${System.getenv('PROJECT')}"
            storage = "scratch/${System.getenv('PROJECT')}"
            cache = 'lenient'
            stageInMode = 'symlink'
        }
    
        singularity {
            // Explicitly turns on container execution
            enabled = true
            // Automatically bind-mount working directory on scratch and common system paths
            autoMounts = true
            // Define location of stored container images 
            cacheDir = "/scratch/${System.getenv('PROJECT')}/${System.getenv('USER')}/nextflow-on-hpc-materials/singularity"
        }
        ```

        2. Save the file and run the pipeline:

        ```bash
        nextflow run main.nf -profile pbspro
        ```


        ??? abstract "Output"

            ```console title="Output"
            executor >  pbspro (3)
            [51/e95140] FASTQC (fastqc on NA12877) | 1 of 1 ✔
            [bb/c54519] ALIGN (1)                  | 1 of 1 ✔
            [93/a8ecf4] GENOTYPE (1)               | 1 of 1, failed: 1 ✘
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
                ...

            Work dir:
             /scratch/vp91/xyz777/nextflow-on-hpc-materials/part2/work/93/a8ecf466a33f67f23eb0ca12b0b753

            Tip: view the complete command output by changing to the process work dir and entering the command `cat .command.out`

            -- Check '.nextflow.log' file for details
            ```
    
    === "Setonix (Slurm)"
    
        ```groovy title="config/slurm.config" hl_lines="1 2 3 4 5 6 14 15 16"
        executor {
            // For high-throughput jobs, these values should be higher
            queueSize = 30
            pollInterval = '5 sec'
            queueStatInterval = '5 sec'
        }

        process {
            // Load the globally installed singularity/4.1.0-slurm module before running any process
            module = 'singularity/4.1.0-slurm'
            // Run using the pbspro scheduler on the 'normalbw' queue
            executor = 'slurm'
            queue = 'work'
            clusterOptions = "--account=${System.getenv('PAWSEY_PROJECT')} --reservation=NextflowHPC"
            cache = 'lenient'
            stageInMode = 'symlink'
        }
    
        singularity {
            // Explicitly turns on container execution
            enabled = true
            // Automatically bind-mount working directory on scratch and common system paths
            autoMounts = true
            // Define location of stored container images 
            cacheDir = "/scratch/${System.getenv('PAWSEY_PROJECT')}/${System.getenv('USER')}/nextflow-on-hpc-materials/singularity"
        }
        ```

        Save your file and run the pipeline:
        ```bash
        nextflow run main.nf -profile slurm
        ```

        ??? abstract "Output"

            ```
             N E X T F L O W   ~  version 24.10.0

            Launching `main.nf` [wise_venter] DSL2 - revision: e34a5e5f9d

            executor >  slurm (6)
            [ec/2e1a61] FASTQ | 1 of 1 ✔
            [6a/e63199] ALIGN | 1 of 1 ✔
            [07/5dea67] GENOT | 1 of 1 ✔
            [76/bac7dd] JOINT | 1 of 1 ✔
            [1e/4c035e] STATS | 1 of 1 ✔
            [f3/c78e59] MULTI | 1 of 1 ✔
            Completed at: 17-Nov-2025 12:41:33
            Duration    : 2m 46s
            CPU hours   : (a few seconds)
            Succeeded   : 6
            ```

!!! question "Zoom react!"

    If your job has finished succesfully, react "Yes" on Zoom, and "No" if it returned an error

Let's explore what resources were actually used and compare them to what was allocated, by inspecting the logs and system job information using the methods from Part 1.

We will use the respective job scheduler introspection tools to observe the resources used on both both systems.

!!! example "Exercise: Inspect resource usage"

    1. Find the Job ID of the failed or completed `GENOTYPE` process in your `.nextflow.log`: 

    ```bash
    grep -w GENOTYPE .nextflow.log | grep jobId
    ```

    The output should look something like:
    ```console
    Nov-11 12:26:17.249 [Task submitter] DEBUG nextflow.executor.GridTaskHandler - [PBSPRO] submitted process GENOTYPE (1) > jobId: 154278383.gadi-pbs; workDir: /scratch/ad78/fj9712/nextflow-on-hpc-materials/part2/work/0b/465f7eacfa500698687ff12df65060
    Nov-11 12:27:47.144 [Task monitor] DEBUG n.processor.TaskPollingMonitor - Task completed > TaskHandler[jobI d: 154278383.gadi-pbs; id: 5; name: GENOTYPE (1); status: COMPLETED; exit: 0; error: -; workDir: /scratch/a d78/fj9712/nextflow-on-hpc-materials/part2/work/0b/465f7eacfa500698687ff12df65060 started: 1762828007140; exited: 2025-11-11T02:27:35Z; ]
    ```

    The value we need in this example is `154278383` - this corresponds to the job id that was scheduled.

    2. Then, use this job ID and the specific job introspection command for your HPC system to view the resource usage of the `GENOTYPE` process:

    === "Gadi (PBS pro)"

        Use `qstat -xf <job_id>` for a comprehensive summary of job allocation, environment variables, and resource usage:

        ```bash
        qstat -xf <job_id>
        ```
        ```console
        Job Id: 154038474.gadi-pbs
        Job_Name = nf-GENOTYPE_1
        Job_Owner = user@gadi-login-05.gadi.nci.org.au
        resources_used.cpupercent = 94
        resources_used.cput = 00:02:50
        resources_used.jobfs = 0b
        resources_used.mem = 512124kb
        resources_used.ncpus = 1
        resources_used.vmem = 512124kb
        resources_used.walltime = 00:03:03
        job_state = F
        queue = normal-exec
        ...
        ```

    === "Setonix (Slurm)"

        Use `seff` for a simple summary of job resource usage:

        ```bash
        seff <job_id>
        ```
        ```console
        Job ID: 34903205
        Cluster: setonix
        User/Group: cou057/cou057
        State: COMPLETED (exit code 0)
        Nodes: 1
        Cores per node: 2
        CPU Utilized: 00:00:19
        CPU Efficiency: 36.54% of 00:00:52 core-walltime
        Job Wall-clock time: 00:00:26
        Memory Utilized: 682.74 MB
        Memory Efficiency: 37.11% of 1.80 GB (920.00 MB/core)
        ```

These reporting tools show that the jobs were assigned the following resources:

|         | CPU   | Memory   |
| ------- | ----- | -------- |
| Gadi    | 1     | 512 MB   |
| Setonix | 2     | 1.8 GB   |


We have observed that the default of 1.8 GB RAM allocated to the Setonix `GENOTYPE` job was sufficient, but the default of 512 MB on Gadi was not.

This is why **explicit resource configuration** is important. Even though the pipeline technically ran (or failed), these defaults are unsuitable for real data.

Many Nextflow workflows provide test data that can be used to help get the pipleline running on your HPC. Researchers can be left confused when a successful test run is followed by a failed run on their own data. The test data is often small and requires minimal resources, while real data can be orders of magnitude larger and more complex. This discrepancy can lead to unexpected failures if the workflow is not properly configured for the actual data being processed. 

We will now go on to look at how we can create a custom configuration file to specify appropriate resources for our data and HPC environment.


## 2.1.4 Custom workflow configuration files

So far we have applied the default `nextflow.config` file that contains details required to run the workflow on any platform. We have then gone on to create and apply an infrastructure-specific configuration file (`config/pbspro.config` or `config/slurm.config`) to tell Nextflow how and where to run on a particular HPC system. Next we will create a third configuration file (`config/custom.config`) to specify resource requirements that are appropriate for our data and HPC environment.

![](figs/00_custom_configs.png)

At this point, you may be wondering "**Why do we have so many configuration files?** We use three different configuration files to keep our Nextflow workflows reproducible, modular, and portable across different systems. 

This setup:

- Ensures consistency when running the same pipeline in different environments
- Allows re-use of configuration components across multiple pipelines
- Simplifies maintenance by isolating system-specific settings from workflow logic

Each configuration file serves a distinct purpose:

- `nextflow.config` is the main configuration file that defines the core behaviour of the workflow itself (e.g. `main.nf`). It includes parameters (params), and references to profiles. To maintain reproducibility, **this file should not be modified during system-specific tuning**. It should only change if the underlying workflow logic changes - that is, what gets run.

- `config/pbspro.config` and `config/slurm.config` define how the pipeline should run on a particular type of HPC system. These files specify details such as which executor to use (e.g. PBS Pro or Slurm), whether to use Singularity or Docker, and other runtime behaviour. They do not control the internal logic of the pipeline. These files should be tailored to match the requirements and setup of the HPC infrastructure you are targeting. Working on a new HPC? You'll need to make a new config file for it! But the good news is you can still ***use the same `nextflow.config` file***.

- `config/custom.config` is an additional system-specific customisation layer that defines process settings such as CPU and memory requests. When developing or adapting a custom pipeline for an HPC environment, this is typically where most tuning happens to fit the specific node architecture, queue constraints, and resource requirements of the data being processed. While these settings *could* be included within the same config that defines the executor and other system-specific settings, separating them into a distinct file allows for easier management and modification of resource allocations without altering the core system configuration. This modularity is especially useful when experimenting with different resource configurations during pipeline development and testing or when adapting the pipeline for different datasets with varying resource needs.

While this structure is a useful starting point, it is not the only way to structure your configuration. The nf-core community have their own set of standards with some presets for some instutitions (including Gadi and Setonix!). However, it is important to double check that these configs are suitable and optimal for your purposes. For more information see [nf-core/configs](https://nf-co.re/configs/).

## 2.1.5 Assigning default resources

Now that we have an appreciation of the layers of configuration for running Nextflow workflows on HPC, we will continue to get the pipeline running by specifying what resources should be run, instead of relying on the system-specific default. We want to ensure that:

- Jobs are being scheduled correctly
- All process tasks and the pipeline complete successfully

Note that the resources assigned differ across sytems - these values are based on the average memory available per core and will be revisited in the upcoming lesson on resourcing.

!!! example "Exercise: Create custom resource config"

    1. Create a new file `config/custom.config`

        ```bash
        touch config/custom.config
        ```

    2. Open the new empty file, and add the following contents based on your HPC

    === "Gadi (PBS pro)"

        ```groovy title='custom.config'
        process {
            cpu = 1 // 'normalbw' queue = 128 GB / 28 CPU ~ 4.6 OR 9.1
            memory = 4.GB
        }
        ```

    === "Setonix (Slurm)"

        ```groovy title='custom.config'
        process {
            cpu = 1 // 'work' partition = 230 GB / 128 CPU ~ 1.8
            memory = 2.GB
        }
        ```

Since we have been repeatedly running the same Nextflow run command, it makes sense to save this into a run script. This can reduce human error (for example missing a flag, typos) and is easier to re-run. This is especially useful in the testing and benchmarking stages and can later be adapted into a job submission script if your HPC has a queue where you can run Nextflow head jobs

!!! example "Exercise: Create and use a run script"

    1. Create a new file called `run.sh`

        ```bash
        touch run.sh
        ```

    2. Open the new empty file, and copy-paste the following code based on your HPC:

    === "Gadi (PBS pro)"

        ```groovy title="run.sh"
        #!/bin/bash

        module load nextflow/24.04.5
        module load singularity

        nextflow run main.nf -profile pbspro -c config/custom.config
        ```

    === "Setonix (Slurm)"

        ```groovy title="run.sh"
        #!/bin/bash

        module load nextflow/24.10.0
        module load singularity/4.1.0-slurm

        nextflow run main.nf -profile slurm -c config/custom.config
        ```

    3. Save the file (Windows: Ctrl+S, macOS: Cmd+S)

    4. Provide execute permission by running:

        ```bash
        chmod +x run.sh
        ```

    5. Run the workflow using the new run script: 

    ```bash
    ./run.sh
    ```

    ??? note "Results"

        On both Gadi and Setonix, the workflow should now be successful and executed end-to-end on the respective scheduler.

        === "Gadi (PBS pro)"

            ```bash
            Loading nextflow/24.04.5
            Loading requirement: java/jdk-17.0.2

             N E X T F L O W   ~  version 24.04.5

            Launching `main.nf` [determined_picasso] DSL2 - revision: 5e5c4f57e0

            executor >  pbspro (6)
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

        === "Setonix (Slurm)"


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

## 2.1.6 Summary

You’ve now built the scaffolding needed to begin fine-tuning your resource requests and exploring monitoring and optimisation techniques. In the next section, we'll start measuring actual resource usage and configuring processes more precisely for efficient implementation on the underlying HPC system.

## 2.1.7 Code checkpoint


??? abstract "Show complete code at the end of Section 2.1"

    === "Gadi (PBS pro)"

        ```groovy title="config/pbspro.config"
        executor {
            // For high-throughput jobs, these values should be higher
            queueSize = 30
            pollInterval = '5 sec'
            queueStatInterval = '5 sec'
        }

        process {
            // Load the globally installed singularity module before running any process
            module = 'singularity'
            // Run using the pbspro scheduler on the 'normalbw' queue
            executor = 'pbspro'
            queue = 'normalbw'
            clusterOptions = "-P ${System.getenv('PROJECT')}"
            storage = "scratch/${System.getenv('PROJECT')}"
            cache = 'lenient'
            stageInMode = 'symlink'
        }

        singularity {
            // Explicitly turns on container execution
            enabled = true
            // Automatically bind-mount working directory on scratch and common system paths
            autoMounts = true
            // Define location of stored container images 
            cacheDir = "/scratch/${System.getenv('PROJECT')}/${System.getenv('USER')}/nextflow-on-hpc-materials/singularity"
        }
        ```

        ```groovy title="run.sh"
        #!/bin/bash

        module load nextflow/24.04.5
        module load singularity

        nextflow run main.nf -profile pbspro -c config/custom.config
        ```

    === "Setonix (Slurm)"  

        ```groovy title="slurm.config"
        executor {
            // For high-throughput jobs, these values should be higher
            queueSize = 30
            pollInterval = '5 sec'
            queueStatInterval = '5 sec'
        }

        process {
            // Load the globally installed singularity/4.1.0-slurm module before running any process
            module = 'singularity/4.1.0-slurm'
            // Run using the pbspro scheduler on the 'normalbw' queue
            executor = 'slurm'
            queue = 'work'
            clusterOptions = "--account=${System.getenv('PAWSEY_PROJECT')}"
            cache = 'lenient'
            stageInMode = 'symlink'
        }

        singularity {
            // Explicitly turns on container execution
            enabled = true
            // Automatically bind-mount working directory on scratch and common system paths
            autoMounts = true
            // Define location of stored container images 
            cacheDir = "/scratch/${System.getenv('PAWSEY_PROJECT')}/${System.getenv('USER')}/nextflow-on-hpc-materials/singularity"
        }
        ```

        ```groovy title="run.sh"
        #!/bin/bash

        module load nextflow/24.10.0
        module load singularity/4.1.0-slurm

        nextflow run main.nf -profile slurm -c config/custom.config
        ```
