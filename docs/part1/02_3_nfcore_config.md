# 2.3 Configuring nf-core

!!! info "Learning objectives"

    - Understand Nextflow executors
    - Understand how Nextflow uses executor definitions to talk to the HPC schedulere
    - Understand the singularity configuration scope
    - Know how to dynamically define the queue based on required resources
    
## 2.3.1 Executors

TODO: Give an overview of the two executors we will use in this workshop - PBSPro and SLURM - and which task resources they request under the hood (cpus, memory, walltime, queue)

!!! example "Exercise: Define the HPC executor"

    TODO:
    
    - Create a new blank file called `gadi.config` or `setonix.config`
    - Define a `process {}` scope in the config file and add `executor = 'pbspro'` or `executor = 'slurm'`
    - Define the `executor {}` scope in the config file and add some basic limits:
        - `queueSize = 30`: How many jobs can be submitted at once
        - `pollInterval = '30 sec'`: How often to check for termination
        - `queueStatInterval = '30 sec'`: How often to check for queue status
        - `submitRateLimit = '20 min'`: Limit rate of job submission: 20/minute
    - Update the run script with `-c <gadi|setonix>.config`

TODO: Either explain that this won't run yet, or have participants try to run again and see for themselves. Won't run because:

- No queue specified, so there's no way to tell the HPC what queue to use. We will configure this at the end of this section.
- Singularity still needs to be defined, otherwise the software would need to be pre-installed on the compute nodes.
- Some additional HPC-specific resources need to be defined:
    - Gadi: `process.storage = "scratch/${System.getenv('PROJECT')}"`
        - Need to note that this is Gadi-specific and won't work on generic PBSPro systems
    - Setonix: `process.clusterOptions = "--account=${System.getenv('PAWSEY_PROJECT')}"`

## 2.3.2 Containers in nf-core

!!! example "Exercise: Define the singularity configuration"

    TODO:

    - Add the `singularity {}` scope with `enabled = true`
        - Optional? Add `autoMounts = true` and `autoCleanUp = true`
            - Need to test if these are necessary. autoCleanUp isn't documented.
    - Add `module = 'singularity/<version>'` to process scope

## 2.3.3 Configuring HPC resources

!!! example "Exercise: Finalise the config with resource requirements"

    TODO:

    - Add common process options:
        - `cache = 'lenient'`: Use file path and size alone to determine when to use the cached process outputs
        - `stageInMode = 'symlink'`: Use symbolic links to stage input files. Required for HPC systems where input files might be stored on different physical drives.
    - Add HPC-specific process options:
        - Gadi: `storage = "scratch/${System.getenv('PROJECT')}"`: Tell gadi to mount the scratch space for your project
            - Need to note that this is Gadi-specific and won't work on generic PBSPro systems
        - Setonix: `process.clusterOptions = "--account=${System.getenv('PAWSEY_PROJECT')}"`: Tell Setonix to submit the job using your project
    - Add dynamic queue specification:
        - Gadi: `queue = { task.memory < 128.GB ? 'normalbw' : 'hugemembw' }`: Use the `normalbw` queue unless you have large memory requirements
            - Simplified from institutional config for readibility, simplicity, and to keep it similar to setonix config
        - Setonix: `queue = { task.memory <= 230.GB ? 'work' : 'highmem' }`: Use the `work` queue unless you have large memory requirements
    - Add `params.trace_timestamp = new java.util.Date().format('yyyy-MM-dd_HH-mm-ss)`
    - Add the trace config for benchmarking purposes
        - `trace.enabled = true`: Enable the trace file
        - `trace.overwrite = false`: Don't overwrite an older trace file: important for keeping track of multiple runs
        - `trace.fields = 'name,status,exit,duration,realtime,cpus,%cpu,memory,%mem,rss'`: Specify the fields you want to track
        - `trace.file = "results/runInfo/trace-${params.trace_timestamp}.txt"`
    - Run the workflow!
        - Note that jobs might be slow to start. After a few minutes we will get everyone to cancel their runs if they haven't completed yet, and explain that we need to fine-tune in the next section to improve our overall run time.