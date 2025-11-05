# 1.5 Running Nextflow on HPC

!!! info "Learning objectives"

    - Understand how Nextflow interacts with HPC components
    - Learn how executors, queues, and work directories control task execution on HPC
    - Run and inspect a simple workflow using a provided HPC configuration
    - Identify which parts of a Nextflow config correspond to scheduler options

!!! note "Nextflow Refresher"

    **Core Concepts**

    A nextflow pipeline consists of three primary components:

    - **Processes** define what to run. Each process can use any Linux-compatible language (e.g., Bash, Python, R, Perl).
    - **Channels** define how data flows between processes. Channels asynchronously carry data between processes and can fan-out (parallel tasks) or fan-in (merge results).
    - **Workflows** Define the order in which processes connect. They orchestrate execution, specifying dependencies and the overall structure of the pipeline.

    Each process runs independently. When a channel contains multiple inputs, Nextflow automatically creates parallel tasks, each running in isolation, connected only by data passed through channels.

    When running locally, these tasks all execute on your own computer using the local executor. This is great for development and small test runs.

    But as datasets grow, your laptop quickly runs out of CPU and memory. That’s where the HPC scheduler takes over.

## 1.5.1 From your laptop to the cluster

In earlier lessons, we saw that HPCs are shared, scheduled, and resource-limited. Nextflow acts as an intermediary, it:

- Submits your workflow's processes to the scheduler
- Handls the movement of data between filesystems
- Monitors job completion.

Each process in a Nextflow pipeline becomes a separate batch job or task array on the cluster if and only if you configure the workflow to interact with the cluster's job scheduler. Nextflow then:

- Prepares a `work/` directory on shared storage
- Submits the process commands to the scheduler
- Moves data between the filesystem and compute nodes as needed
- Checks for completion and retrieves logs, outputs, and exit codes
- Publishes the output data to the shared filesystem

TODO diagram demonstrating this, something like: https://sydney-informatics-hub.github.io/template-nf-guide/notebooks/demo.html

## 1.5.2 Our first HPC workflow

We'll use a demo workflow, [config-demo-nf](https://github.com/Sydney-Informatics-Hub/config-demo-nf) to see this in action. This workflow contains a single process that splits a sequence into multiple files.

![](figs/00_confid_demo_nf.png)

[TODO] get oversight on figure and update (fix file). Fix demo to include appropriate channels

!!! example "Download and run the workflow"

    Use git to clone the workflow code base to your working directory:

    ```bash
    git clone https://github.com/Sydney-Informatics-Hub/config-demo-nf.git
    ```

    Then load the nextflow module, following the same method we used in lesson 1.2:

    === "Gadi"
        ```bash
        module load nextflow
        ```
    === "Setonix"
        ```bash
        module load nextflow/version
        ```

    Now run the workflow:

    ```bash
    nextflow run config-demo-nf/main.nf
    ```
    The output of your command should look something like this
    === "Gadi"
        ```console
        N E X T F L O W   ~  version 24.04.5
        Launching `main.nf` [jovial_lalande] DSL2 - revision: 4e4c91df36
        executor >  local (1)
        [c1/104707] splitSequences | 1 of 1 ✔
        ```

    === "Setonix"
        ```console
        N E X T F L O W   ~  version 24.10.0
        Launching `main.nf` [prickly_lovelace] DSL2 - revision: 66bfcf5bb9
        executor >  local (1)
        [6b/e8feb6] splitSequences [100%] 1 of 1 ✔
        ```

    What does each line mean?

    1. The version of Nextflow that was executed
    2. The script and version names
    3. The executor used (in the above case: local)
    4. The process that was executed once, which means there is one task. The line starts with a unique hexadecimal value, and ends with the task completion information

### Task directories and the `work/` folder

When you run a Nextflow pipeline, it automatically creates a `work/` directory. This is where all computation happens behind the scenes.
Inside this directory, each process execution (or task) runs in its own isolated subdirectory, identified by a unique hash, in the above example, `work/6b/e8feb6` (NOTE: your unique hash will be different).

!!! note

    You can execute `tree work/` to view the work directory structure.

    ```bash
    tree work/
    ```

    ```console title="Output"
    work
    └── 6b
    └── e8feb6a83bb78a7a6661ccc1211857
        ├── seq_1
        ├── seq_2
        ├── seq_3
        └── sequence.fa -> /scratch/PROJECT/USER/nextflow-on-hpc-materials/part1/config-demo-nf/sequence.fa
    ```

Here’s what happens inside each task directory:

1. Setup: Nextflow stages (copies or links) the input files, plus a small script (.command.sh) that defines what to run.
2. Execution: The process runs inside that folder, writing its results there.
3. Cleanup: Nextflow collects the output files and makes them available for downstream processes or publishing.

Each directory is independent so tasks don’t share writable space. If one process needs data from another, it’s passed through Nextflow channels, not shared files. This isolation is especially important on HPC systems, where tasks may run on different compute nodes.

### Where did my task actually run?

Our first run used the local executor, which means all computation happened directly on the login node, the same machine where we typed the command.
This is perfectly fine for quick tests or debugging, but not suitable for real workloads on HPC systems.

On HPC, heavy computation should be handled by compute nodes managed by a job scheduler.

To make that happen, we’ll use Nextflow configuration files.

## 1.5.3 Configuring for the scheduler

A Nextflow configuration file (`nextflow.config` or files inside a `config/` directory) defines how and where your workflow runs without changing the workflow code itself.

It can specify:

- Executor: Which system to use (e.g., local, slurm, pbspro).
- Queue: Defines where jobs run within the scheduler (e.g., normal, highmem, gpu).
- Work Directory: Defines where intermediate files are stored so compute nodes can access them.
- Resources: CPU, memory, and time per process.
- Environment: Modules, containers, or conda environments to load.
- Storage: Where to store temporary and output files.

Configs are powerful on HPC systems, because they are what connect your workflow to the scheduler.
They translate Nextflow’s processes into properly submitted batch jobs.

Because configs are separate from the workflow logic, you can:

- Run the same pipeline on laptop, HPC, or cloud by changing only the config.
- Tune performance by adjusting resources or queues per process
- Adapt workflows to site-specific environments (modules, scratch paths, queue names)
- Share portable workflows that others can run on their HPC without code changes

In short, configs are what make Nextflow workflows portable, scalable, and cluster-aware.

!!! example "Running the workflow on the compute nodes"

    === "Gadi"

    We have pre-made a very simple configuration file, `pbspro.config`, that will allow the example Nextflow pipeline to run on Gadi. Go ahead and re-run the workflow, adding the new configuration file with the `-c pbspro.config` option. You will also need to define a new parameter: `pbspro_account`:

    ```bash
    nextflow run config-demo-nf/main.nf -c pbspro.config --pbspro_account vp91
    ```

    === "Setonix"

    We have pre-made a very simple configuration file, `slurm.config`, that will allow the example Nextflow pipeline to run on Setonix. Go ahead and re-run the workflow, adding the new configuration file with the `-c slurm.config` option. You will also need to define a new parameter: `slurm_account`:

    ```bash
    nextflow run config-demo-nf/main.nf -c slurm.config --slurm_account courses01
    ```

    The output of your command should now look something like this:

    === "Gadi"
        ```bash
        N E X T F L O W   ~  version 24.04.5

        Launching `main.nf` [lethal_gilbert] DSL2 - revision: 4e4c91df36

        executor >  pbspro (1)
        [a8/5345da] splitSequences | 1 of 1 ✔
        ```
    === "Setonix"
        ```bash
        N E X T F L O W   ~  version 24.10.0

        Launching `main.nf` [nice_boltzmann] DSL2 - revision: 4e4c91df36

        executor >  slurm (1)
        [67/d497fa] splitSequences [100%] 1 of 1 ✔
        ```
    Notice that the executor now matches your HPC’s system, slurm on Setonix or pbspro on Gadi.

## 1.5.4 Profiles

Another very useful feature of Nextflow is the ability to bundle up configuration options into **profiles**. This can help to simplify the command line arguments to Nextflow by using the `-profile <profile name>` syntax, rather than having to provide the path to the relevant configuration file. We have already set up the `nextflow.config` file to define two profiles, `pbspro` and `slurm`, which import the relevant configuraiton files when they are used:

```groovy title="nextflow.config" linenums="4"
// Define HPC profiles to run with job scheduler
profiles {
  // Use this profile to interact with the scheduler on setonix 
  slurm { includeConfig "slurm.config" }

  // Use this profile to interact with the scheduler on gadi   
  pbspro { includeConfig "pbspro.config" }
}
```

!!! example "Running the workflow on the compute nodes with profiles"

    Run the workflow once more, this time using the executor profiles:

    === "Gadi"
        ```bash
         nextflow run config-demo-nf/main.nf -profile pbspro --pbspro_account vp91
        ```
    === "Setonix"
        ```bash
        nextflow run config-demo-nf/main.nf -profile slurm --slurm_account courses01
        ```

    The output of your command should be the same as before:

    === "Gadi"
        ```bash
        N E X T F L O W   ~  version 24.04.5

        Launching `main.nf` [lethal_gilbert] DSL2 - revision: 4e4c91df36

        executor >  pbspro (1)
        [a8/5345da] splitSequences | 1 of 1 ✔
        ```
    === "Setonix"
        ```bash
        N E X T F L O W   ~  version 24.10.0

        Launching `main.nf` [nice_boltzmann] DSL2 - revision: 4e4c91df36

        executor >  slurm (1)
        [67/d497fa] splitSequences [100%] 1 of 1 ✔
        ```

Now that we've recapped the basics of Nextflow and seen how a simple Nextflow pipeline can be run on an HPC, in the next section we will look at how we can start running some more complex workflows on HPCs, and how we can configure them to more efficiently utilise the resources available to them.