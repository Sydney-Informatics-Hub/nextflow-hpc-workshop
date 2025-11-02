# 1.5 Running Nextflow on HPC

!!! info "Learning objectives"

    - Understand how Nextflow interacts with HPC components
    - Learn how executors, queues, and work directories control task execution on HPC
    - Run and inspect a simple workflow using a provided HPC configuration
    - Identify which parts of a Nextflow config correspond to scheduler options

## 1.5.1 From your laptop to the cluster

In earlier lessons, we saw that HPCs are shared, scheduled, and resource-limited. Nextflow acts as an intermediary, it:

- Submits your workflow's processes to the scheduler
- Handls the movement of data between filesystems
- Monitors job completion.

Each process in a Nextflow pipeline becomes a separate batch job or task array on the cluster if and only if you configure the workflow to interact with the cluster's job scheduler. Nextflow then: 

- Prepares a `work/` directory on shared storage
- Submits the process commands to the scheduler 
- Waits for the scheduler to run it on a compute node
- Checks for completion and retrieves logs, outputs, and exit codes 
- Publishes the output data to the shared filesystem

TODO diagram demonstrating this, something like: https://sydney-informatics-hub.github.io/template-nf-guide/notebooks/demo.html 

## 1.5.2 Our first HPC workflow 

We'll use a demo workflow, [config-demo-nf](https://github.com/Sydney-Informatics-Hub/config-demo-nf) to see this in action. This workflow contains a single process that splits a sequence into multiple files.

TODO workflow diagram, code base components, connected to nextflow run command like: https://sydney-informatics-hub.github.io/template-nf-guide/notebooks/components.html

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
    ```console
    TODO output
    ```

TODO explore the output printed to the screen, like https://sydney-informatics-hub.github.io/template-nf-guide/notebooks/demo.html#run-the-demo

### 1.5.3 Configuring for the scheduler

