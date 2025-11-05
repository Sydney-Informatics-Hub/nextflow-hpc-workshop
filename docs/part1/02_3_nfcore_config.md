# 2.3 Configuring nf-core

!!! info "Learning objectives"

    - Understand Nextflow executors
    - Understand how Nextflow uses executor definitions to talk to the HPC schedulere
    - Understand the singularity configuration scope
    - Know how to dynamically define the queue based on required resources

## 2.3.1 Executors

Executors are the back-end that Nextflow talks to to run the tasks in your workflow. By default, Nextflow will assume that you want to run everything on the same computer that you ran the `nextflow run` command on. But, as we have learned, that is definitely **not** how we want to run things on a shared HPC system: we want each task to be submitted to the scheduler to run on a compute node with all the appropriate resources.

We saw earlier today that there are a fair few parameters that need to be configured when submitting jobs to an HPC, and these differ between systems. Luckily, Nextflow includes native support for the two HPC schedulers that we are working with in this workshop: PBSPro and Slurm.

To set up Nextflow to use an HPC executor, we simply define the `process.executor` configuration option in the Nextflow configuration. We can also configure a few other parameters to control how many jobs get submitted to the HPC and how frequently; this is useful with large pipelines to avoid overwhelming the system (and angering the admins!). To keep things clean, we will create a new, blank configuration file called `hpc.config` and specify it in our run command. As we will see, Nextflow lets you **layer** configurations on top of one another and combines them in a predictable way to allow fine control of how each process runs.

!!! example "Exercise: Define the HPC executor"

    Start by creating a new blank file in the `config` directory called `hpc.config`. You can do this via the VSCode explorer window (right-click the `config` folder and select "New File...") or via the terminal:

    ```bash
    touch config/hpc.config
    ```

    Open the new file in the VSCode editor and add a `process {}` scope:

    ```groovy
    process {

    }
    ```

    Within the process scope, define the `executor` option and set it to the relevant executor for your system:

    === "Gadi"

        ```groovy hl_lines="2"
        process {
            executor = 'pbspro'
        }
        ```

    === "Setonix"

        ```groovy hl_lines="2"
        process {
            executor = 'slurm'
        }
        ```

    We also want to set some limits to how many jobs can be submitted at once and how frequently they get submitted. These settings are important, because many large pipelines can create potentially hundreds of jobs that may overwhelm the system. Most HPCs will set a limit for how many jobs a user can submit at once, and your pipeline may fail if it tries to submit more than this limit.

    For our purposes, we will keep our queued job limit to 30, and limit the number of jobs we can submit at once to 20 per minute. We will also tell Nextflow to request for status updates on our jobs once every 30 seconds.

    === "Gadi"

        ```groovy hl_lines="5-10"
        process {
            executor = 'pbspro'
        }

        executor {
            queueSize = 30
            submitRateLimit = '20 min'
            pollInterval = '30 sec'
            queueStatInterval = '30 sec'
        }
        ```

    === "Setonix"

        ```groovy hl_lines="5-10"
        process {
            executor = 'slurm'
        }

        executor {
            queueSize = 30
            submitRateLimit = '20 min'
            pollInterval = '30 sec'
            queueStatInterval = '30 sec'
        }
        ```

    Now we have defined our executor and some relevant settings for it, we will need to tell Nextflow to actually use this new configuration file; by default, Nextflow will only use the `nextflow.config` file in the project directory, and will only load other configuration files when explicitly told to do so.

    In the `run.sh` script, add the following highlighted line to the `nextflow run` command:

    === "Gadi"

        ```bash title="run.sh" linenums="1" hl_lines="15"
        #!/bin/bash

        module load nextflow/24.04.5
        module load singularity

        nextflow run sarek/main.nf \
            --input ../data/bams/samplesheet.csv \
            --fasta ../data/ref/Hg38.subsetchr20-22.fasta \
            --fasta_fai ../data/ref/Hg38.subsetchr20-22.fasta.fai \
            --step markduplicates \
            --skip_tools baserecalibrator,mosdepth,samtools \
            --outdir results \
            --no_intervals true \
            --igenomes_ignore true \
            -c config/hpc.config \
            -resume
        ```

    === "Setonix"

        ```bash title="run.sh" linenums="1" hl_lines="15"
        #!/bin/bash

        module load nextflow/24.10.0
        module load singularity/4.1.0-slurm

        nextflow run sarek/main.nf \
            --input ../data/bams/samplesheet.csv \
            --fasta ../data/ref/Hg38.subsetchr20-22.fasta \
            --fasta_fai ../data/ref/Hg38.subsetchr20-22.fasta.fai \
            --step markduplicates \
            --skip_tools baserecalibrator,mosdepth,samtools \
            --outdir results \
            --no_intervals true \
            --igenomes_ignore true \
            -c config/hpc.config \
            -resume
        ```

    The new line, `-c config/hpc.config \`, tells Nextflow to load the new configuration file and merge it with the existing configuration set up by the `nextflow.config` file. Note that the settings in configuration files provided by the `-c` command will take precedence over those set in the `nextflow.config` file, so if any options are specified in both files, the setting in `config/hpc.config` will be used. We will explore layering configurations further in the next section of the workshop.

We've now set up a basic configuration file to specify our HPC executor, but we aren't done yet: we still haven't set up Nextflow to use Singularity. If we were to run the script now, it would still fail, just slower than before, as the tasks would get submitted to the default queue, wait for a little bit before starting to run on a compute node, then quickly fail when they couldn't find the appropriate software installed. In the next section, we will set up Nextflow to use Singularity to run each tool.

## 2.3.2 Containers in nf-core

!!! example "Exercise: Define the singularity configuration"

    At the bottom of our `config/hpc.config` configuration file, we will now add a `singularity` scope and enable the containerisation software. At the same time, we will also define the Singularity cache directory. This is the directory where Singularity should store all downloaded containers so that it doesn't need to download them over and over again whenever the same tool is required. As part of the setup work we did earlier today, we have already created this cache directory within the `sarek` folder, at `sarek/singularity/`. We can define this in the `singularity` configuration scope by setting `cacheDir = "$projectDir/singularity"`, where `"$projectDir"` is a Nextflow variable that refers to the directory in which the `main.nf` script is located (in our case, the `sarek/` directory).

    === "Gadi"

        ```groovy hl_lines="11-15"
        process {
            executor = 'pbspro'
        }

        executor {
            queueSize = 30
            submitRateLimit = '20 min'
            pollInterval = '30 sec'
            queueStatInterval = '30 sec'
        }

        singularity {
            enabled = true
            cacheDir = "$projectDir/singularity"
        }
        ```

    === "Setonix"

        ```groovy hl_lines="11-15"
        process {
            executor = 'slurm'
        }

        executor {
            queueSize = 30
            submitRateLimit = '20 min'
            pollInterval = '30 sec'
            queueStatInterval = '30 sec'
        }

        singularity {
            enabled = true
            cacheDir = "$projectDir/singularity"
        }
        ```

    Enabling the use of Singularity will tell Nextflow to run our tasks using a `singularity exec` command, similar to what we used earlier today. However, you may remember that the `singularity` command isn't available to use by default on the HPC systems: we needed to run `module load` first. If we tried to run the workflow now, we would get an error that `singularity` couldn't be found. Luckily, Nextflow has us covered here once again: the `process.module` configuration option lets us define modules that we want to load when running a process. Go ahead and update the `process` scope to define the `singularity` module:

    === "Gadi"

        ```groovy hl_lines="3"
        process {
            executor = 'pbspro'
            module = 'singularity'
        }

        executor {
            queueSize = 30
            submitRateLimit = '20 min'
            pollInterval = '30 sec'
            queueStatInterval = '30 sec'
        }

        singularity {
            enabled = true
            cacheDir = "$projectDir/singularity"
        }
        ```

    === "Setonix"

        ```groovy hl_lines="3"
        process {
            executor = 'slurm'
            module = 'singularity/4.1.0-slurm'
        }

        executor {
            queueSize = 30
            submitRateLimit = '20 min'
            pollInterval = '30 sec'
            queueStatInterval = '30 sec'
        }

        singularity {
            enabled = true
            cacheDir = "$projectDir/singularity"
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

    Let's start with defining the project we want to use. Most HPCs will assign each user a default project to use when submitting jobs, but it is good practice to be explicit about it, especially if you are part of several HPC projects and switch between the often.

    Nextflow doesn't have a way to natively set the project parameter for our HPCs, but it does let us define arbitrary parameters to pass to the scheduler via the `process.clusterOptions` setting. We will use that now:

    === "Gadi"

        On Gadi, we set the project via the `-P` option. We will use the groovy function `System.getenv()` to grab the value of the `$PROJECT` environment variable, which holds our default project ID, and pass that to the `-P` option:

        ```groovy hl_lines="4"
        process {
            executor = 'pbspro'
            module = 'singularity'
            clusterOptions = "-P ${System.getenv('PROJECT')}"
        }

        executor {
            queueSize = 30
            submitRateLimit = '20 min'
            pollInterval = '30 sec'
            queueStatInterval = '30 sec'
        }

        singularity {
            enabled = true
            cacheDir = "$projectDir/singularity"
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
            pollInterval = '30 sec'
            queueStatInterval = '30 sec'
        }

        singularity {
            enabled = true
            cacheDir = "$projectDir/singularity"
        }
        ```

    === "Setonix"

        On Setonix, we set the project via the `--account` option. We will use the groovy function `System.getenv()` to grab the value of the `$PAWSEY_PROJECT` environment variable, which holds our default project ID, and pass that to the `--account` option:

        ```groovy hl_lines="4"
        process {
            executor = 'slurm'
            module = 'singularity/4.1.0-slurm'
            clusterOptions = "--account=${System.getenv('PAWSEY_PROJECT')}"
        }

        executor {
            queueSize = 30
            submitRateLimit = '20 min'
            pollInterval = '30 sec'
            queueStatInterval = '30 sec'
        }

        singularity {
            enabled = true
            cacheDir = "$projectDir/singularity"
        }
        ```

    The next thing that we will define is the HPC queue that we wish to submit our jobs to. For the purposes of our example today, we only need the basic queue on each system. However, it is good practice to **dynamically specify** the queue based on resource requirements; that way, large jobs won't fail or be rejected entirely by the HPC scheduler due to invalid resource requests.

    Nextflow lets us defin the queue that we want via the `queue` option in the `process` scope. We can dynamically specify the queue by using curly braces to wrap around a conditional statement that tests how much memory each job needs:

    === "Gadi"

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
            pollInterval = '30 sec'
            queueStatInterval = '30 sec'
        }

        singularity {
            enabled = true
            cacheDir = "$projectDir/singularity"
        }
        ```

    === "Setonix"

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
            pollInterval = '30 sec'
            queueStatInterval = '30 sec'
        }

        singularity {
            enabled = true
            cacheDir = "$projectDir/singularity"
        }
        ```

    We want to add just a couple of extra options to the process definition. The first option is the `stageInMode` option. We will explicitly tell Nextflow that we want to use **symbolic links**. These are essentially shortcuts that point to another file on the system, and let us refer to inputs within our working directory without physically copying them in, which would use up lots of additional storage space. To set this, we define `stageInMode = 'symlink'` in the `process` scope.

    The second option we want to set is the `cache` mode. Nextflow lets us use the outputs of previous runs when re-running a pipeline, by specifying the `-resume` flag on the command line. This is very useful for avoiding re-running jobs when we don't have to. By default, Nextflow uses various features of a file, including its timestamp, to determine if it has changed and whether a job needs to be re-run or not. However, the shared filesystem on HPCs can interfere with the timestamps and cause jobs to re-run when they don't need to. Nextflow provides a workaround for this, by letting us use a `lenient` mode that ignores the timestamp.

    Let's set these two options now:

    === "Gadi"

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
            pollInterval = '30 sec'
            queueStatInterval = '30 sec'
        }

        singularity {
            enabled = true
            cacheDir = "$projectDir/singularity"
        }
        ```

    === "Setonix"

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
            pollInterval = '30 sec'
            queueStatInterval = '30 sec'
        }

        singularity {
            enabled = true
            cacheDir = "$projectDir/singularity"
        }
        ```

    We now have a fully-functioning HPC configuration file! We will, however, add just one more feature that will help us monitor the resources we are using and optimise our workflow. This is the `trace` file, and the next few lines that we add will enable it, set its file name (including a time stamp), and set the values that we want to keep track of:

    === "Gadi"

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
            pollInterval = '30 sec'
            queueStatInterval = '30 sec'
        }

        singularity {
            enabled = true
            cacheDir = "$projectDir/singularity"
        }

        params.trace_timestamp = new java.util.Date().format('yyyy-MM-dd_HH-mm-ss')

        trace {
            enabled = true
            overwrite = false
            file = "./runInfo/trace-${params.trace_timestamp}.txt"
            fields = 'name,status,exit,duration,realtime,cpus,%cpu,memory,%mem,rss'
        }
        ```

    === "Setonix"

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
            pollInterval = '30 sec'
            queueStatInterval = '30 sec'
        }

        singularity {
            enabled = true
            cacheDir = "$projectDir/singularity"
        }

        params.trace_timestamp = new java.util.Date().format('yyyy-MM-dd_HH-mm-ss')

        trace {
            enabled = true
            overwrite = false
            file = "./runInfo/trace-${params.trace_timestamp}.txt"
            fields = 'name,status,exit,duration,realtime,cpus,%cpu,memory,%mem,rss'
        }
        ```

    Let's quickly break down this new code:

    - `params.trace_timestamp = new java.util.Date().format('yyyy-MM-dd_HH-mm-ss')`: This sets a custom parameter called `trace_timestamp` and sets it to the current date and time. This will let us create a unique file for every run.
    - `trace { ... }`: This defines the trace file scope, and all options within are specific to defining that file.
        - `enabled = true`: This simply enables the use of the trace file
        - `overwrite = false`: This prevents a trace file from being overwritten
        - `file = "./runInfo/trace-${params.trace_timestamp}.txt"`: This sets the file path for the trace file, using the `trace_timestamp` parameter we set just above
        - `fields = 'name,status,exit,duration,realtime,cpus,%cpu,memory,%mem,rss'`: This defines the fields that we want to capture within the trace file, including the task name, its status, how long it ran for, and how efficiently it used the CPUs and memory provided to it.

    And that's it! You are now ready to re-run the workflow and Nextflow will now know how to submit the jobs to your assigned HPC and how to use Singularity to run each job.

    ```bash
    ./run.sh
    ```

    ??? question "Result..."

        After a few moments as the pipeline starts up, you should notice the tasks getting submitted to the HPC:

        ```console title="Output"
        executor >  pbspro (5)
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:BWAMEM1_INDEX                                              -
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:BWAMEM2_INDEX                                              -
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:DRAGMAP_HASHTABLE                                          -
        [f3/b449ba] process > NFCORE_SAREK:PREPARE_GENOME:GATK4_CREATESEQUENCEDICTIONARY (Hg38.subsetchr20-22.fasta) [  0%] 0 of 1
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:MSISENSORPRO_SCAN                                          -
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:SAMTOOLS_FAIDX                                             -
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:TABIX_BCFTOOLS_ANNOTATIONS                                 -
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:TABIX_DBSNP                                                -
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:TABIX_GERMLINE_RESOURCE                                    -
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:TABIX_KNOWN_SNPS                                           -
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:TABIX_KNOWN_INDELS                                         -
        [-        ] process > NFCORE_SAREK:PREPARE_GENOME:TABIX_PON                                                  -
        [39/60ebaa] process > NFCORE_SAREK:PREPARE_INTERVALS:TABIX_BGZIPTABIX_INTERVAL_COMBINED (no_intervals)       [  0%] 0 of 1
        [ab/088f09] process > NFCORE_SAREK:SAREK:BAM_MARKDUPLICATES:GATK4_MARKDUPLICATES (test_sample1)              [  0%] 0 of 3
        [-        ] process > NFCORE_SAREK:SAREK:BAM_MARKDUPLICATES:CRAM_QC_MOSDEPTH_SAMTOOLS:SAMTOOLS_STATS         -
        [-        ] process > NFCORE_SAREK:SAREK:BAM_MARKDUPLICATES:CRAM_QC_MOSDEPTH_SAMTOOLS:MOSDEPTH               -
        [-        ] process > NFCORE_SAREK:SAREK:CRAM_TO_BAM                                                         -
        [-        ] process > NFCORE_SAREK:SAREK:CRAM_SAMPLEQC:BAM_NGSCHECKMATE:BCFTOOLS_MPILEUP                     -
        [-        ] process > NFCORE_SAREK:SAREK:CRAM_SAMPLEQC:BAM_NGSCHECKMATE:NGSCHECKMATE_NCM                     -
        [-        ] process > NFCORE_SAREK:SAREK:MULTIQC
        ```

        You might notice that the jobs take a while to start. The reason for this is because we're using the `sarek` default resource requirements for each process. `sarek` has been designed around whole genome sequencing data, which is usually quite large and requires a lot of computing power to run. However, for our workshop today, we are working with a very small test dataset that only requires 1 CPU for each job and at most 5GB of memory. If your pipeline has finished, you can inspect the trace file within the `runInfo` folder and see that the CPU and memory percentage is quite low for many of the tasks:

        ```bash
        # Your trace file will have a unique name based on the time it was run
        cat runInfo/trace-2025-11-18_13-14-15.txt
        ```

        ```console title="Output"
        name	status	exit	duration	realtime	cpus	%cpu	memory	%mem	rss
        NFCORE_SAREK:PREPARE_GENOME:GATK4_CREATESEQUENCEDICTIONARY (Hg38.subsetchr20-22.fasta)	COMPLETED	0	27.1s	4s	6	216.9%	36 GB	0.1%	307 MB
        NFCORE_SAREK:SAREK:BAM_MARKDUPLICATES:GATK4_MARKDUPLICATES (test_sample3)	COMPLETED	0	29.2s	7s	6	205.9%	30 GB	2.8%	7.1 GB
        NFCORE_SAREK:PREPARE_INTERVALS:TABIX_BGZIPTABIX_INTERVAL_COMBINED (no_intervals)	COMPLETED	0	15.6s	0ms	1	77.8%	1 GB	0.0%	2 MB
        NFCORE_SAREK:SAREK:BAM_MARKDUPLICATES:GATK4_MARKDUPLICATES (test_sample1)	COMPLETED	0	21.2s	7s	6	204.6%	30 GB	2.3%	5.8 GB
        NFCORE_SAREK:SAREK:BAM_MARKDUPLICATES:GATK4_MARKDUPLICATES (test_sample2)	COMPLETED	0	24.2s	6s	6	197.9%	30 GB	2.5%	6.2 GB
        NFCORE_SAREK:SAREK:MULTIQC	COMPLETED	0	28.8s	9.2s	4	83.5%	12 GB	0.2%	455.6 MB
        ```

        If your pipeline hasn't finished after a few minutes, you can cancel the run with a `Ctrl + C` keyboard combination. In the final section for today, we will create another configuration file to layer on top of our exisitng configuration and fine-tune our tasks to run more efficiently.
