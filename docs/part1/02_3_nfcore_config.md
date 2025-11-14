# 2.3 Configuring nf-core

!!! info "Learning objectives"

    - Understand Nextflow executors
    - Understand how Nextflow uses executor definitions to talk to the HPC schedulere
    - Understand the singularity configuration scope
    - Know how to dynamically define the queue based on required resources

## 2.3.1 Executors

Executors are the back-end that Nextflow talks to to run the tasks in your workflow. By default, Nextflow will assume that you want to run everything on the same computer that you ran the `nextflow run` command on. But, as we have learned, that is definitely **not** how we want to run things on a shared HPC system: we want each task to be submitted to the scheduler to run on a compute node with all the appropriate resources.

We saw earlier today that there are a fair few parameters that need to be configured when submitting jobs to an HPC, and these differ between systems. Luckily, Nextflow includes native support for the two HPC schedulers that we are working with in this workshop: PBSPro and Slurm.

To set up Nextflow to use an HPC executor, we simply define the `process.executor` configuration option in the Nextflow configuration. We can also configure a few other parameters to control how many jobs get submitted to the HPC and how frequently; this is useful with large pipelines to avoid overwhelming the system (and angering the admins!). To keep things clean, we will create a new, blank configuration file called either `gadi.config` or `setonix.config`, depending on your system, and specify it in our run command. As we will see, Nextflow lets you **layer** configurations on top of one another and combines them in a predictable way to allow fine control of how each process runs.

!!! example "Exercise: Define the HPC executor"

    Start by creating a new blank file in the `config` directory called `hpc.config`. You can do this via the VSCode explorer window (right-click the `config` folder and select "New File...") or via the terminal:

    === "Gadi (PBS)"

        ```bash
        touch config/gadi.config
        ```

    === "Setonix (Slurm)"

        ```bash
        touch config/setonix.config
        ```

    Open the new file in the VSCode editor and add a `process {}` scope:

    ```groovy
    process {

    }
    ```

    Within the process scope, define the `executor` option and set it to the relevant executor for your system:

    === "Gadi (PBS)"

        ```groovy hl_lines="2"
        process {
            executor = 'pbspro'
        }
        ```

    === "Setonix (Slurm)"

        ```groovy hl_lines="2"
        process {
            executor = 'slurm'
        }
        ```

    We also want to set some limits to how many jobs can be submitted at once and how frequently they get submitted. These settings are important, because many large pipelines can create potentially hundreds of jobs that may overwhelm the system. Most HPCs will set a limit for how many jobs a user can submit at once, and your pipeline may fail if it tries to submit more than this limit.

    For our purposes, we will keep our queued job limit to 30, and limit the number of jobs we can submit at once to 20 per minute. We will also tell Nextflow to request for status updates on our jobs once every 15 seconds.

    === "Gadi (PBS)"

        ```groovy hl_lines="5-10"
        process {
            executor = 'pbspro'
        }

        executor {
            queueSize = 30
            submitRateLimit = '20 min'
            pollInterval = '15 sec'
            queueStatInterval = '15 sec'
        }
        ```

    === "Setonix (Slurm)"

        ```groovy hl_lines="5-10"
        process {
            executor = 'slurm'
        }

        executor {
            queueSize = 30
            submitRateLimit = '20 min'
            pollInterval = '15 sec'
            queueStatInterval = '15 sec'
        }
        ```

    Now we have defined our executor and some relevant settings for it, we will need to tell Nextflow to actually use this new configuration file; by default, Nextflow will only use the `nextflow.config` file in the project directory, and will only load other configuration files when explicitly told to do so.

    In the `run.sh` script, update the following highlighted lines by adding a ` \` to the end of the old command and adding the new configuration file with the `-c` option:

    === "Gadi (PBS)"

        ```bash title="run.sh" linenums="1" hl_lines="16-17"
        #!/bin/bash

        module load nextflow/24.04.5
        module load singularity

        nextflow run sarek/main.nf \
            --input ../data/fqs/samplesheet.single.csv \
            --fasta ../data/ref/Hg38.subsetchr20-22.fasta \
            --fasta_fai ../data/ref/Hg38.subsetchr20-22.fasta.fai \
            --dict ../data/ref/Hg38.subsetchr20-22.dict \
            --bwa ../data/ref \
            --step mapping \
            --skip_tools markduplicates,baserecalibrator,mosdepth,samtools \
            --outdir results \
            --no_intervals true \
            --igenomes_ignore true \
            -c config/gadi.config
        ```

    === "Setonix (Slurm)"

        ```bash title="run.sh" linenums="1" hl_lines="17"
        #!/bin/bash

        module load nextflow/24.10.0
        module load singularity/4.1.0-slurm

        nextflow run sarek/main.nf \
            --input ../data/fqs/samplesheet.single.csv \
            --fasta ../data/ref/Hg38.subsetchr20-22.fasta \
            --fasta_fai ../data/ref/Hg38.subsetchr20-22.fasta.fai \
            --dict ../data/ref/Hg38.subsetchr20-22.dict \
            --bwa ../data/ref \
            --step mapping \
            --skip_tools markduplicates,baserecalibrator,mosdepth,samtools \
            --outdir results \
            --no_intervals true \
            --igenomes_ignore true \
            -c config/setonix.config
        ```

    The new line adds the `-c <config>` option to the `nextflow run` command. We provide the path to our new configuration file which tells Nextflow to load it and merge it with the existing configuration set up by the `nextflow.config` file. Note that the settings in configuration files provided by the `-c` command will take precedence over those set in the `nextflow.config` file, so if any options are specified in both files, the setting in `config/gadi.config` or `config/setonix.config` will be used. We will explore layering configurations further in the next section of the workshop.

    Go ahead and try running the script:

    ```bash
    ./run.sh
    ```

    What happens?

    ??? question "Result..."

        The workflow will begin to run and submit several jobs to the executor:

        ```console title="Initial output"
        N E X T F L O W   ~  version 24.04.5

        Launching `sarek/main.nf` [exotic_cray] DSL2 - revision: 3954909713

        ...

        executor >  pbspro (1)
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:BWAMEM1_INDEX                                   -
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:BWAMEM2_INDEX                                   -
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:DRAGMAP_HASHTABLE                               -
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:GATK4_CREATESEQUENCEDICTIONARY                  -
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:MSISENSORPRO_SCAN                               -
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:SAMTOOLS_FAIDX                                  -
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:TABIX_BCFTOOLS_ANNOTATIONS                      -
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:TABIX_DBSNP                                     -
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:TABIX_GERMLINE_RESOURCE                         -
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:TABIX_KNOWN_SNPS                                -
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:TABIX_KNOWN_INDELS                              -
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:TABIX_PON                                       -
        [-        ] process > NFCORE_SAREK:PREPARE_INTERVALS:TABIX_BGZIPTABIX_INTERVAL_COMBINED           [  0%] 0 of 1
        [-        ] process > NFCORE_SAREK:SAREK:SPRING_DECOMPRESS_TO_FQ_PAIR                             -
        [-        ] process > NFCORE_SAREK:SAREK:SPRING_DECOMPRESS_TO_R1_FQ                               -
        [-        ] process > NFCORE_SAREK:SAREK:SPRING_DECOMPRESS_TO_R2_FQ                               -
        [-        ] process > NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:SAMTOOLS_VIEW_MAP_MAP                -
        [-        ] process > NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:SAMTOOLS_VIEW_UNMAP_UNMAP            -
        [-        ] process > NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:SAMTOOLS_VIEW_UNMAP_MAP              -
        [-        ] process > NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:SAMTOOLS_VIEW_MAP_UNMAP              -
        [-        ] process > NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:SAMTOOLS_MERGE_UNMAP                 -
        [-        ] process > NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:COLLATE_FASTQ_UNMAP                  -
        [-        ] process > NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:COLLATE_FASTQ_MAP                    -
        [-        ] process > NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:CAT_FASTQ                            -
        [5f/63da05] process > NFCORE_SAREK:SAREK:FASTQC (test_sample1-all)                                [  0%] 0 of 1
        [-        ] process > NFCORE_SAREK:SAREK:FASTP                                                    [  0%] 0 of 1
        [-        ] process > NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:BWAMEM1_MEM     -
        [-        ] process > NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:BWAMEM2_MEM     -
        [-        ] process > NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:DRAGMAP_ALIGN   -
        [-        ] process > NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:SENTIEON_BWAMEM -
        [-        ] process > NFCORE_SAREK:SAREK:BAM_MERGE_INDEX_SAMTOOLS:MERGE_BAM                       -
        [-        ] process > NFCORE_SAREK:SAREK:BAM_MERGE_INDEX_SAMTOOLS:INDEX_MERGE_BAM                 -
        [-        ] process > NFCORE_SAREK:SAREK:BAM_TO_CRAM_MAPPING                                      -
        [-        ] process > NFCORE_SAREK:SAREK:CRAM_QC_NO_MD:SAMTOOLS_STATS                             -
        [-        ] process > NFCORE_SAREK:SAREK:CRAM_QC_NO_MD:MOSDEPTH                                   -
        [-        ] process > NFCORE_SAREK:SAREK:CRAM_TO_BAM                                              -
        [-        ] process > NFCORE_SAREK:SAREK:CRAM_SAMPLEQC:BAM_NGSCHECKMATE:BCFTOOLS_MPILEUP          -
        [-        ] process > NFCORE_SAREK:SAREK:CRAM_SAMPLEQC:BAM_NGSCHECKMATE:NGSCHECKMATE_NCM          -
        [-        ] process > NFCORE_SAREK:SAREK:MULTIQC
        ```

        However, soon after the jobs start running on the compute nodes, they fail:

        ```console title="Final output"
        -[nf-core/sarek] Pipeline completed with errors-
        ERROR ~ Error executing process > 'NFCORE_SAREK:SAREK:FASTP (test_sample1-all)'

        Caused by:
        Process `NFCORE_SAREK:SAREK:FASTP (test_sample1-all)` terminated with an error exit status (127)

        ...

        Command exit status:
        127

        Command output:
        (empty)

        Command error:
        .command.sh: line 11: fastp: command not found

        Work dir:
        /scratch/er01/PIPE-5866-NF4HPC/nextflow-on-hpc-materials/part1/work/92/74509fb605260a345f7c3f29419f1c

        Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

        -- Check '.nextflow.log' file for details
        ERROR ~ Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting

        -- Check '.nextflow.log' file for details
        ```

        The above example shows exactly what went wrong: `fastp: command not found`. We haven't yet configured Nextflow to use Singularity, so it is assuming the software is installed on the compute node. In the next section, we will set up Nextflow to use Singularity to run each tool.

## 2.3.2 Containers in nf-core

!!! example "Exercise: Define the singularity configuration"

    At the bottom of our configuration file, we will now add a `singularity` scope and enable the containerisation software. At the same time, we will also define the Singularity cache directory. This is the directory where Singularity should store all downloaded containers so that it doesn't need to download them over and over again whenever the same tool is required. As part of the setup work we did earlier today, we have already created this cache directory within the parent folder, at `../singularity/`. We can define this in the `singularity` configuration scope by setting the `cacheDir` option. We will provide the full path to this file, which is at `/scratch/<PROJECT>/<USER>/nextflow-on-hpc-materials/singularity`. The current project ID and username are both accessible as environment variables on both Gadi and Setonix; we access these with the groovy function `System.getenv()`:

    === "Gadi (PBS)"

        ```groovy hl_lines="11-15"
        process {
            executor = 'pbspro'
        }

        executor {
            queueSize = 30
            submitRateLimit = '20 min'
            pollInterval = '15 sec'
            queueStatInterval = '15 sec'
        }

        singularity {
            enabled = true
            cacheDir = "/scratch/${System.getenv('PROJECT')}/${System.getenv('USER')}/nextflow-on-hpc-materials/singularity"
        }
        ```

    === "Setonix (Slurm)"

        ```groovy hl_lines="11-15"
        process {
            executor = 'slurm'
        }

        executor {
            queueSize = 30
            submitRateLimit = '20 min'
            pollInterval = '15 sec'
            queueStatInterval = '15 sec'
        }

        singularity {
            enabled = true
            cacheDir = "/scratch/${System.getenv('PAWSEY_PROJECT')}/${System.getenv('USER')}/nextflow-on-hpc-materials/singularity"
        }
        ```

    Enabling the use of Singularity will tell Nextflow to run our tasks using a `singularity exec` command, similar to what we used earlier today. However, you may remember that the `singularity` command isn't available to use by default on the HPC systems: we needed to run `module load` first. If we tried to run the workflow now, we would get an error that `singularity` couldn't be found. Luckily, Nextflow has us covered here once again: the `process.module` configuration option lets us define modules that we want to load when running a process. Go ahead and update the `process` scope to define the `singularity` module:

    === "Gadi (PBS)"

        ```groovy hl_lines="3"
        process {
            executor = 'pbspro'
            module = 'singularity'
        }

        executor {
            queueSize = 30
            submitRateLimit = '20 min'
            pollInterval = '15 sec'
            queueStatInterval = '15 sec'
        }

        singularity {
            enabled = true
            cacheDir = "/scratch/${System.getenv('PROJECT')}/${System.getenv('USER')}/nextflow-on-hpc-materials/singularity"
        }
        ```

    === "Setonix (Slurm)"

        ```groovy hl_lines="3"
        process {
            executor = 'slurm'
            module = 'singularity/4.1.0-slurm'
        }

        executor {
            queueSize = 30
            submitRateLimit = '20 min'
            pollInterval = '15 sec'
            queueStatInterval = '15 sec'
        }

        singularity {
            enabled = true
            cacheDir = "/scratch/${System.getenv('PAWSEY_PROJECT')}/${System.getenv('USER')}/nextflow-on-hpc-materials/singularity"
        }
        ```

We now have a configuration file with both our executor defined and singularity enabled. There are just a few finishing touches we need to make, primarily around defining the HPC queue that we want to use and some additional options for how to handle files.

## 2.3.3 Configuring HPC resources

!!! example "Exercise: Finalise the config with resource requirements"

    To finish off our HPC configuration file, we are going to set the following options:

    - Specifically tell the HPC scheduler which project we want to use
    - Define the HPC queue that we want our jobs submitted to
    - Define how files should be staged in the working directories
    - Define how Nextflow should determine when to use the outputs of previous runs when using the `-resume` flag
    - Enable the trace file so we can track our resource usage and optimise our jobs

    Let's start with defining the project we want to use. Most HPCs will assign each user a default project to use when submitting jobs, but it is good practice to be explicit about it, especially if you are part of several HPC projects and switch between the often. We can do this with the `process.clusterOptions` setting, which lets us pass arbitrary parameters to the scheduler:

    === "Gadi (PBS)"

        On Gadi, we set the project via the `-P` option. We will again use the groovy function `System.getenv()` to grab the value of the `$PROJECT` environment variable, which holds our default project ID, and pass that to the `-P` option:

        ```groovy hl_lines="4"
        process {
            executor = 'pbspro'
            module = 'singularity'
            clusterOptions = "-P ${System.getenv('PROJECT')}"
        }

        executor {
            queueSize = 30
            submitRateLimit = '20 min'
            pollInterval = '15 sec'
            queueStatInterval = '15 sec'
        }

        singularity {
            enabled = true
            cacheDir = "/scratch/${System.getenv('PROJECT')}/${System.getenv('USER')}/nextflow-on-hpc-materials/singularity"
        }
        ```

        On Gadi, we also need to explicitly tell the HPC to mount the scratch space for our project. Normally, we would do this with the `-l storage=scratch/<PROJECT ID>` option to the `qsub` command; in this case, we could use the `clusterOptions` setting to specify the scratch space parameter as well. However, the Nextflow version that is available on Gadi has a special configuration option built in that lets us simplify this a bit: instead, we can use the `storage` option within the `process` scope to specify the scratch space (and any other storage locations we might like to mount). **Note** that this is **NOT** a standard Nextflow feature, and is only **specific to Gadi**:

        ```groovy hl_lines="5"
        process {
            executor = 'pbspro'
            module = 'singularity'
            clusterOptions = "-P ${System.getenv('PROJECT')}"
            storage = "scratch/${System.getenv('PROJECT')}"
        }

        executor {
            queueSize = 30
            submitRateLimit = '20 min'
            pollInterval = '15 sec'
            queueStatInterval = '15 sec'
        }

        singularity {
            enabled = true
            cacheDir = "/scratch/${System.getenv('PROJECT')}/${System.getenv('USER')}/nextflow-on-hpc-materials/singularity"
        }
        ```

    === "Setonix (Slurm)"

        On Setonix, we set the project via the `--account` option. We will again use the groovy function `System.getenv()` to grab the value of the `$PAWSEY_PROJECT` environment variable, which holds our default project ID, and pass that to the `--account` option:

        ```groovy hl_lines="4"
        process {
            executor = 'slurm'
            module = 'singularity/4.1.0-slurm'
            clusterOptions = "--account=${System.getenv('PAWSEY_PROJECT')}"
        }

        executor {
            queueSize = 30
            submitRateLimit = '20 min'
            pollInterval = '15 sec'
            queueStatInterval = '15 sec'
        }

        singularity {
            enabled = true
            cacheDir = "/scratch/${System.getenv('PAWSEY_PROJECT')}/${System.getenv('USER')}/nextflow-on-hpc-materials/singularity"
        }
        ```

    The next thing that we will define is the HPC queue that we wish to submit our jobs to. For the purposes of our example today, we only need the basic queue on each system. However, it is good practice to **dynamically specify** the queue based on resource requirements; that way, large jobs won't fail or be rejected entirely by the HPC scheduler due to invalid resource requests.

    Nextflow lets us define the queue that we want via the `queue` option in the `process` scope. We can dynamically specify the queue by using curly braces to wrap around a conditional statement that tests how much memory each job needs:

    === "Gadi (PBS)"

        On Gadi, the `normalbw` queue supports tasks with up to 128GB of memory. If we need more than that, we want to use the `hugemembw` queue. We can achieve this using a short-hand if-else statement in groovy: `<condition> ? <value if true> : <value if false>`. We can ask whether the memory required by the current task (`task.memory`) is less than 128GB; if so, we set `queue` to `normalbw`, otherwise we set it to `hugemembw`. By wrapping the whole statement in curly braces, we ensure that it is evaluated when the task runs:

        ```groovy hl_lines="6"
        process {
            executor = 'pbspro'
            module = 'singularity'
            clusterOptions = "-P ${System.getenv('PROJECT')}"
            storage = "scratch/${System.getenv('PROJECT')}"
            queue = { task.memory < 128.GB ? 'normalbw' : 'hugemembw' }
        }

        executor {
            queueSize = 30
            submitRateLimit = '20 min'
            pollInterval = '15 sec'
            queueStatInterval = '15 sec'
        }

        singularity {
            enabled = true
            cacheDir = "/scratch/${System.getenv('PROJECT')}/${System.getenv('USER')}/nextflow-on-hpc-materials/singularity"
        }
        ```

        !!! note "Read the docs!"

            Remember to check [NCI's "Queue Limits" page](https://opus.nci.org.au/spaces/Help/pages/236881198/Queue+Limits) when configuring the `queue` for your pipelines on Gadi, as it contains important information about how and when to select each particular queue. We are using a fairly naive selection method today for simplicity, but more complex queue selection methods are possible and advisable for larger pipelines.

    === "Setonix (Slurm)"

        On Setonix, the `work` queue supports tasks with up to 230GB of memory. If we need more than that, we want to use the `highmem` queue. We can achieve this using a short-hand if-else statement in groovy: `<condition> ? <value if true> : <value if false>`. We can ask whether the memory required by the current task (`task.memory`) is less than 230GB; if so, we set `queue` to `work`, otherwise we set it to `highmem`. By wrapping the whole statement in curly braces, we ensure that it is evaluated when the task runs:

        ```groovy hl_lines="5"
        process {
            executor = 'slurm'
            module = 'singularity/4.1.0-slurm'
            clusterOptions = "--account=${System.getenv('PAWSEY_PROJECT')}"
            queue = { task.memory < 230.GB ? 'work' : 'highmem' }
        }

        executor {
            queueSize = 30
            submitRateLimit = '20 min'
            pollInterval = '15 sec'
            queueStatInterval = '15 sec'
        }

        singularity {
            enabled = true
            cacheDir = "/scratch/${System.getenv('PAWSEY_PROJECT')}/${System.getenv('USER')}/nextflow-on-hpc-materials/singularity"
        }
        ```

        !!! note "Read the docs!"

            Remember to check [Pawsey's "Running Jobs on Setonix" page](https://pawsey.atlassian.net/wiki/spaces/US/pages/51929058/Running+Jobs+on+Setonix) when configuring the `queue` for your pipelines on Setonix, as it contains important information about how and when to select each particular queue. We are using a fairly naive selection method today for simplicity, but more complex queue selection methods are possible and advisable for larger pipelines.

    We want to add just a couple of extra options to the process definition. The first option is the `stageInMode` option. We will explicitly tell Nextflow that we want to use **symbolic links**. These are essentially shortcuts that point to another file on the system, and let us refer to inputs within our working directory without physically copying them in, which would use up lots of additional storage space. To set this, we define `stageInMode = 'symlink'` in the `process` scope.

    The second option we want to set is the `cache` mode. Nextflow lets us use the outputs of previous runs when re-running a pipeline, by specifying the `-resume` flag on the command line. This is very useful for avoiding re-running jobs when we don't have to. By default, Nextflow uses various features of a file, including its timestamp, to determine if it has changed and whether a job needs to be re-run or not. However, the shared filesystem on HPCs can interfere with the timestamps and cause jobs to re-run when they don't need to. Nextflow provides a workaround for this, by letting us use a `lenient` mode that ignores the timestamp.

    Let's set these two options now:

    === "Gadi (PBS)"

        ```groovy hl_lines="7-8"
        process {
            executor = 'pbspro'
            module = 'singularity'
            clusterOptions = "-P ${System.getenv('PROJECT')}"
            storage = "scratch/${System.getenv('PROJECT')}"
            queue = { task.memory < 128.GB ? 'normalbw' : 'hugemembw' }
            stageInMode = 'symlink'
            cache = 'lenient'
        }

        executor {
            queueSize = 30
            submitRateLimit = '20 min'
            pollInterval = '15 sec'
            queueStatInterval = '15 sec'
        }

        singularity {
            enabled = true
            cacheDir = "/scratch/${System.getenv('PROJECT')}/${System.getenv('USER')}/nextflow-on-hpc-materials/singularity"
        }
        ```

    === "Setonix (Slurm)"

        ```groovy hl_lines="6-7"
        process {
            executor = 'slurm'
            module = 'singularity/4.1.0-slurm'
            clusterOptions = "--account=${System.getenv('PAWSEY_PROJECT')}"
            queue = { task.memory < 230.GB ? 'work' : 'highmem' }
            stageInMode = 'symlink'
            cache = 'lenient'
        }

        executor {
            queueSize = 30
            submitRateLimit = '20 min'
            pollInterval = '15 sec'
            queueStatInterval = '15 sec'
        }

        singularity {
            enabled = true
            cacheDir = "/scratch/${System.getenv('PAWSEY_PROJECT')}/${System.getenv('USER')}/nextflow-on-hpc-materials/singularity"
        }
        ```

    We now have a fully-functioning HPC configuration file! We will, however, add just one more feature that will help us monitor the resources we are using and optimise our workflow. This is the `trace` file, and the next few lines that we add will enable it, set its file name (including a time stamp), and set the values that we want to keep track of:

    === "Gadi (PBS)"

        ```groovy hl_lines="22-30"
        process {
            executor = 'pbspro'
            module = 'singularity'
            clusterOptions = "-P ${System.getenv('PROJECT')}"
            storage = "scratch/${System.getenv('PROJECT')}"
            queue = { task.memory < 128.GB ? 'normalbw' : 'hugemembw' }
            stageInMode = 'symlink'
            cache = 'lenient'
        }

        executor {
            queueSize = 30
            submitRateLimit = '20 min'
            pollInterval = '15 sec'
            queueStatInterval = '15 sec'
        }

        singularity {
            enabled = true
            cacheDir = "/scratch/${System.getenv('PROJECT')}/${System.getenv('USER')}/nextflow-on-hpc-materials/singularity"
        }

        params.trace_timestamp = new java.util.Date().format('yyyy-MM-dd_HH-mm-ss')

        trace {
            enabled = true
            overwrite = false
            file = "./runInfo/trace-${params.trace_timestamp}.txt"
            fields = 'name,status,exit,duration,realtime,cpus,%cpu,memory,%mem,peak_rss'
        }
        ```

    === "Setonix (Slurm)"

        ```groovy hl_lines="21-30"
        process {
            executor = 'slurm'
            module = 'singularity/4.1.0-slurm'
            clusterOptions = "--account=${System.getenv('PAWSEY_PROJECT')}"
            queue = { task.memory < 230.GB ? 'work' : 'highmem' }
            stageInMode = 'symlink'
            cache = 'lenient'
        }

        executor {
            queueSize = 30
            submitRateLimit = '20 min'
            pollInterval = '15 sec'
            queueStatInterval = '15 sec'
        }

        singularity {
            enabled = true
            cacheDir = "/scratch/${System.getenv('PAWSEY_PROJECT')}/${System.getenv('USER')}/nextflow-on-hpc-materials/singularity"
        }

        params.trace_timestamp = new java.util.Date().format('yyyy-MM-dd_HH-mm-ss')

        trace {
            enabled = true
            overwrite = false
            file = "./runInfo/trace-${params.trace_timestamp}.txt"
            fields = 'name,status,exit,duration,realtime,cpus,%cpu,memory,%mem,peak_rss'
        }
        ```

    Let's quickly break down this new code:

    - `params.trace_timestamp = new java.util.Date().format('yyyy-MM-dd_HH-mm-ss')`: This sets a custom parameter called `trace_timestamp` and sets it to the current date and time. This will let us create a unique file for every run.
    - `trace { ... }`: This defines the trace file scope, and all options within are specific to defining that file.
        - `enabled = true`: This simply enables the use of the trace file
        - `overwrite = false`: This prevents a trace file from being overwritten
        - `file = "./runInfo/trace-${params.trace_timestamp}.txt"`: This sets the file path for the trace file, using the `trace_timestamp` parameter we set just above
        - `fields = 'name,status,exit,duration,realtime,cpus,%cpu,memory,%mem,peak_rss'`: This defines the fields that we want to capture within the trace file, including the task name, its status, how long it ran for, and how efficiently it used the CPUs and memory provided to it.

    And that's it! You are now ready to re-run the workflow and Nextflow will now know how to submit the jobs to your assigned HPC and how to use Singularity to run each job.

    ```bash
    ./run.sh
    ```

    ??? question "Result..."

        After a few moments as the pipeline starts up, you should notice the tasks getting submitted to the HPC:

        ```console title="Output"
        N E X T F L O W   ~  version 24.04.5

        Launching `sarek/main.nf` [exotic_cray] DSL2 - revision: 3954909713

        ...

        executor >  pbspro (1)
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:BWAMEM1_INDEX                                   -
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:BWAMEM2_INDEX                                   -
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:DRAGMAP_HASHTABLE                               -
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:GATK4_CREATESEQUENCEDICTIONARY                  -
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:MSISENSORPRO_SCAN                               -
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:SAMTOOLS_FAIDX                                  -
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:TABIX_BCFTOOLS_ANNOTATIONS                      -
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:TABIX_DBSNP                                     -
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:TABIX_GERMLINE_RESOURCE                         -
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:TABIX_KNOWN_SNPS                                -
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:TABIX_KNOWN_INDELS                              -
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:TABIX_PON                                       -
        [-        ] process > NFCORE_SAREK:PREPARE_INTERVALS:TABIX_BGZIPTABIX_INTERVAL_COMBINED           [  0%] 0 of 1
        [-        ] process > NFCORE_SAREK:SAREK:SPRING_DECOMPRESS_TO_FQ_PAIR                             -
        [-        ] process > NFCORE_SAREK:SAREK:SPRING_DECOMPRESS_TO_R1_FQ                               -
        [-        ] process > NFCORE_SAREK:SAREK:SPRING_DECOMPRESS_TO_R2_FQ                               -
        [-        ] process > NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:SAMTOOLS_VIEW_MAP_MAP                -
        [-        ] process > NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:SAMTOOLS_VIEW_UNMAP_UNMAP            -
        [-        ] process > NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:SAMTOOLS_VIEW_UNMAP_MAP              -
        [-        ] process > NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:SAMTOOLS_VIEW_MAP_UNMAP              -
        [-        ] process > NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:SAMTOOLS_MERGE_UNMAP                 -
        [-        ] process > NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:COLLATE_FASTQ_UNMAP                  -
        [-        ] process > NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:COLLATE_FASTQ_MAP                    -
        [-        ] process > NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:CAT_FASTQ                            -
        [5f/63da05] process > NFCORE_SAREK:SAREK:FASTQC (test_sample1-all)                                [  0%] 0 of 1
        [-        ] process > NFCORE_SAREK:SAREK:FASTP                                                    [  0%] 0 of 1
        [-        ] process > NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:BWAMEM1_MEM     -
        [-        ] process > NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:BWAMEM2_MEM     -
        [-        ] process > NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:DRAGMAP_ALIGN   -
        [-        ] process > NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:SENTIEON_BWAMEM -
        [-        ] process > NFCORE_SAREK:SAREK:BAM_MERGE_INDEX_SAMTOOLS:MERGE_BAM                       -
        [-        ] process > NFCORE_SAREK:SAREK:BAM_MERGE_INDEX_SAMTOOLS:INDEX_MERGE_BAM                 -
        [-        ] process > NFCORE_SAREK:SAREK:BAM_TO_CRAM_MAPPING                                      -
        [-        ] process > NFCORE_SAREK:SAREK:CRAM_QC_NO_MD:SAMTOOLS_STATS                             -
        [-        ] process > NFCORE_SAREK:SAREK:CRAM_QC_NO_MD:MOSDEPTH                                   -
        [-        ] process > NFCORE_SAREK:SAREK:CRAM_TO_BAM                                              -
        [-        ] process > NFCORE_SAREK:SAREK:CRAM_SAMPLEQC:BAM_NGSCHECKMATE:BCFTOOLS_MPILEUP          -
        [-        ] process > NFCORE_SAREK:SAREK:CRAM_SAMPLEQC:BAM_NGSCHECKMATE:NGSCHECKMATE_NCM          -
        [-        ] process > NFCORE_SAREK:SAREK:MULTIQC
        ```

        You might notice that the jobs take a while to start. The reason for this is because we're using the `sarek` default resource requirements for each process. `sarek` has been designed around whole genome sequencing data, which is usually quite large and requires a lot of computing power to run. However, for our workshop today, we are working with a very small test dataset that only requires 1 CPU for each job and at most 5GB of memory. If your pipeline has finished, you can inspect the trace file within the `runInfo` folder and see that the CPU and memory percentage is quite low for many of the tasks:

        ```bash
        # Your trace file will have a unique name based on the time it was run
        cat runInfo/trace-2025-11-18_13-14-15.txt
        ```

        ```console title="Output"
        name    status  exit    duration        realtime        cpus    %cpu    memory  %mem    peak_rss
        NFCORE_SAREK:PREPARE_INTERVALS:TABIX_BGZIPTABIX_INTERVAL_COMBINED (no_intervals)        COMPLETED       0       12.9s   0ms     1       78.4%   1 GB    0.0%    2 MB
        NFCORE_SAREK:SAREK:FASTQC (test_sample1-all)    COMPLETED       0       13.3s   2s      4       201.4%  4 GB    0.1%    328.9 MB
        NFCORE_SAREK:SAREK:FASTP (test_sample1-all)     COMPLETED       0       16.3s   1s      12      85.9%   4 GB    0.0%    2 MB
        NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:BWAMEM1_MEM (test_sample1)  COMPLETED       0       9.6s    0ms     24      102.3%  30 GB   0.0%    2 MB
        NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:BWAMEM1_MEM (test_sample1)  COMPLETED       0       12.8s   1s      24      96.4%   30 GB   0.0%    2 MB
        NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:BWAMEM1_MEM (test_sample1)  COMPLETED       0       9.8s    0ms     24      105.9%  30 GB   0.0%    2 MB
        NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:BWAMEM1_MEM (test_sample1)  COMPLETED       0       11.5s   0ms     24      104.2%  30 GB   0.0%    2 MB
        NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:BWAMEM1_MEM (test_sample1)  COMPLETED       0       10.4s   0ms     24      27.1%   30 GB   0.0%    2 MB
        NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:BWAMEM1_MEM (test_sample1)  COMPLETED       0       13.4s   1s      24      44.5%   30 GB   0.0%    2 MB
        NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:BWAMEM1_MEM (test_sample1)  COMPLETED       0       12.5s   0ms     24      230.8%  30 GB   0.0%    2 MB
        NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:BWAMEM1_MEM (test_sample1)  COMPLETED       0       14.5s   0ms     24      275.7%  30 GB   0.0%    2 MB
        NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:BWAMEM1_MEM (test_sample1)  COMPLETED       0       11.5s   0ms     24      34.9%   30 GB   0.0%    2 MB
        NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:BWAMEM1_MEM (test_sample1)  COMPLETED       0       13.4s   0ms     24      106.8%  30 GB   0.0%    2 MB
        NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:BWAMEM1_MEM (test_sample1)  COMPLETED       0       15.4s   0ms     24      90.1%   30 GB   0.0%    2 MB
        NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:BWAMEM1_MEM (test_sample1)  COMPLETED       0       12.4s   1s      24      53.0%   30 GB   0.0%    2 MB
        NFCORE_SAREK:SAREK:BAM_MERGE_INDEX_SAMTOOLS:MERGE_BAM (test_sample1)    COMPLETED       0       14.3s   0ms     2       78.0%   12 GB   0.0%    2 MB
        NFCORE_SAREK:SAREK:BAM_MERGE_INDEX_SAMTOOLS:INDEX_MERGE_BAM (test_sample1)      COMPLETED       0       9.9s    0ms     1       87.4%   1 GB    0.0%    4 MB
        NFCORE_SAREK:SAREK:BAM_TO_CRAM_MAPPING (test_sample1)   COMPLETED       0       14.9s   1s      2       65.0%   4 GB    0.0%    4 MB
        NFCORE_SAREK:SAREK:MULTIQC      COMPLETED       0       29.9s   20.3s   4       72.6%   12 GB   0.2%    661.9 MB
        ```

        The amount of memory each job requested is in the `mem` column of this file, while the actual amount used is in the final `peak_rss` column. Note how most jobs are requesting several GB of memory, but are actually only using a few MB:
        
        ```bash
        # Pull out just the task name, memory requested, and memory used columns
        cut -f 1,8,10 runInfo/trace-2025-11-18_13-14-15.txt
        ```

        ```console title="Output"
        name    memory  peak_rss
        NFCORE_SAREK:PREPARE_INTERVALS:TABIX_BGZIPTABIX_INTERVAL_COMBINED (no_intervals)        1 GB    2 MB
        NFCORE_SAREK:SAREK:FASTQC (test_sample1-all)    4 GB    328.9 MB
        NFCORE_SAREK:SAREK:FASTP (test_sample1-all)     4 GB    2 MB
        NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:BWAMEM1_MEM (test_sample1)  30 GB   2 MB
        NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:BWAMEM1_MEM (test_sample1)  30 GB   2 MB
        NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:BWAMEM1_MEM (test_sample1)  30 GB   2 MB
        NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:BWAMEM1_MEM (test_sample1)  30 GB   2 MB
        NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:BWAMEM1_MEM (test_sample1)  30 GB   2 MB
        NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:BWAMEM1_MEM (test_sample1)  30 GB   2 MB
        NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:BWAMEM1_MEM (test_sample1)  30 GB   2 MB
        NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:BWAMEM1_MEM (test_sample1)  30 GB   2 MB
        NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:BWAMEM1_MEM (test_sample1)  30 GB   2 MB
        NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:BWAMEM1_MEM (test_sample1)  30 GB   2 MB
        NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:BWAMEM1_MEM (test_sample1)  30 GB   2 MB
        NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:BWAMEM1_MEM (test_sample1)  30 GB   2 MB
        NFCORE_SAREK:SAREK:BAM_MERGE_INDEX_SAMTOOLS:MERGE_BAM (test_sample1)    12 GB   2 MB
        NFCORE_SAREK:SAREK:BAM_MERGE_INDEX_SAMTOOLS:INDEX_MERGE_BAM (test_sample1)      1 GB    4 MB
        NFCORE_SAREK:SAREK:BAM_TO_CRAM_MAPPING (test_sample1)   4 GB    4 MB
        NFCORE_SAREK:SAREK:MULTIQC      12 GB   661.9 MB
        ```

        Also note how many instances of `BWAMEM1_MEM` ran: 12. This is a scatter-gather pattern at work: the `FASTP` job breaks up the input FASTQs into multiple smaller files (12 in this case), each of which gets independently processed by `bwa mem`. The BAMs generated by `bwa mem` then get merged back together by the `MERGE_BAM` process. However, for our dataset, 12 parallel `bwa mem` processes is a bit overkill.
        
        These are signs that we could significantly optimise the pipeline for our dataset, which is what we will do in the next lesson.

        !!! note
        
            If your pipeline hasn't finished after a few minutes, you can cancel the run with a `Ctrl + C` keyboard combination. In the final section for today, we will create another configuration file to layer on top of our exisitng configuration and fine-tune our tasks to run more efficiently.

!!! question "How are you going?"

    If you're following along so far, let us know by reacting on zoom with a **":material-check:{ .check } Yes"**.
    
    If you're running into any issues, please react with a **":material-close:{ .close } No"** and we can help out before we move on to the next section.