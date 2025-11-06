# Configuring custom Nextflow pipelines

!!! info "Learning objectives"

    - Identify best practices for benchmarking and configuring pipelines,
    step-by-step
    - Troubleshoot common HPC and scheduling errors by inspecting task logs
    and interpretting exit codes and error messages
    - Know how to find and specify system-dependent HPC values such as queues/partitions
    - Recall how Nextflow interacts with HPC components
We start with running things on a single sample. This should be representative
of all the data we run i.e. we will be processing chromosomes 20, 21, 22.

Then, continue to optimise and configure the pipeline with a single sample so
that you can conserve SUs as you develop and benchmark.

Run the pipeline with bare bones configuration. Specify the correct `profile`
for your system.

!!! example "Exercises"

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

!!! example "Exercises"

    TODO add the QoL additions from nf-core config 1.2.

    === "Gadi (PBS)"

        TODO: Build up to the config used in part 1.2 for
        [Gadi](https://github.com/Sydney-Informatics-Hub/nextflow-on-hpc-materials/blob/main/part1/config/gadi.config) and

    === "Setonix (Slurm)"

        TODO: Build up to the config used in part 1.2 for
        [Setonix](https://github.com/Sydney-Informatics-Hub/nextflow-on-hpc-materials/blob/main/part1/config/setonix.config)

TODO: add why we use each scope and directive. Mainly for reference and demo-ing
QoL things for running pipelines, rather that teaching HPC concepts. e.g. this
is how we usually do things.

<!-- Attendees to copy and paste things into their configs, lead trainer to add
in code in scopes, and talk through what each thing does. -->

TODO Relate this back and compare why Nextflow is powerful for HPC vs. serially
running pbs/slurm scripts.

!!! example "Exercises"

    For both Pawsey (Slurm) and Gadi (PBS) create a new file
    `conf/custom.config` and add the following:

    ```groovy title='custom.config
    process {
        cpu = 1
        memory = 2.GB
    }
    ```

TODO: Explain the `-c` flag. When would you use it?

!!! example "Exercises"

    TODO add run script. Goal: so we don't have to specify the same flags each time
    (e.g. -profile, --<sched>_account). Tie back to reproducibility and flexibility
    of running NXF across different envs - run.sh is one of the things we adapt
    for the system (i.e. project name, scheduler) so everything else can stay
    the same.

    Run script should look something like:

    ```bash
    #load modules ...

    
    nextflow run main.nf -profile <sched> --<sched>_account <project> -c custom.config
    ```

!!! example "Exercise"

    TODO Maybe, qstat and sacct/seff the GENOTYPE process again. What changed?

qstat/sacct-ing each job is inefficient, especially with pipelines with more
processes, and running on more than one sample. Segue into the next section
where this can be automated using Nextflow's in-built monitoring features.
