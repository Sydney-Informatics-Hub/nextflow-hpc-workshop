# Pipeline monitoring and reporting

!!! info "Learning objectives"

    - Apply Nextflow's monitoring features to understand process-level resource usage on HPC
    - Analyse key resource metrics to assess workflow efficiency and identify processes that require tuning 
    - Use nextflow log and custom trace fields to locate, track, and debug workflow execution

In bioinformatics workflows, resource requirements are often not fixed for a given analysis type or tool. They can vary greatly for many reasons, for exampleinout data size, genomic complexity, tissue expression profile, and species. This means we can't assume CPU, memory, or time values will work the same for every dataset or even every sample within a dataset. On HPC systems - where resources are shared and allocations may be finite and/or charged - these differences matter. 

For these reasons, it is not efficient or reliable to just "set and forget" resource values for an entire workflow: you need to build in flexibility. This requires you to have visiblity over pipeline behaviour at the process level. Nextflow provides several monitoring and reporting tools that help you understand this behaviour including the `nextflow log` command, trace files, and reports. 

Now that our workflow is running without error on the HPC compute nodes, we will enable [Nextflow's reports](https://nextflow.io/docs/latest/reports.html). This allows us to view the resource usage of each process for our representative sample.

## 2.2.1 Execution reports and timelines

Nextflow can produce an [execution report](https://nextflow.io/docs/latest/reports.html#execution-report) at the end of your workflow run that summarises all process execution metrics. Similarly, it can create an [execution timeline](https://nextflow.io/docs/latest/reports.html#execution-timeline). These reports summarise how your workflow ran, which processes were executed, and how long they took. They are very helpful during development, troubleshooting, and performance optimisation. 

These reports can be created when running the pipeline using the `-with-report` and `-with-timeline` flags, or by adding `timeline{}` and `report{}` directives to our configuration file. Inclusion within the configuration file has many advantages, providing a means of ensuring the reports are always generated when the pipeline is run, regardless of the execution command used, and allowing for further customisation.

Both `-with-report` and `-with-timeline` allow us to specify a custom file name, and choose whether or not their output files can be overwritten each run. Rather than overwriting the same report file every time we run the pipeline, we will add a timestamp parameter that automatically labels each report with the exact date and time the workflow was launched to easily ensure over-write never occurs. This makes it easier to track multiple runs, especially when you are iterating quickly and comparing resource usage. 

!!! example "Exercise: Enable execution and timeline reports"

    We will enable execution and timeline reports in our `custom.config` file. 
    
    1. Copy the following code and paste at the end of your `custom.config`:

    === "Gadi (PBS pro)"

        ```groovy title='config/custom.config'
        // Name the reports according to when they were run
        params.timestamp = new java.util.Date().format('yyyy-MM-dd_HH-mm-ss')

        // Generate timeline-timestamp.html timeline report 
        timeline {
            enabled = true
            overwrite = false
            file = "./runInfo/timeline-${params.timestamp}.html"
        }

        // Generate report-timestamp.html execution report 
        report {
            enabled = true
            overwrite = false
            file = "./runInfo/report-${params.timestamp}.html"
        }
        ```

    === "Setonix (Slurm)"

        ```groovy title='config/custom.config'
        // Name the reports according to when they were run
        params.timestamp = new java.util.Date().format('yyyy-MM-dd_HH-mm-ss')

        // Generate timeline-timestamp.html timeline report 
        timeline {
            enabled = true
            overwrite = false
            file = "./runInfo/timeline-${params.timestamp}.html"
        }

        // Generate report-timestamp.html execution report 
        report {
            enabled = true
            overwrite = false
            file = "./runInfo/report-${params.timestamp}.html"
        }
        ```
  
    2. Save the confug, and use your run script to run the workflow: 
    ```bash
    ./run.sh
    ```

Once the job completes, you should have a new folder called `runInfo`, where your timeline and report files are saved.



!!! question "Questions"
    
    1. In VScode, right click on the `report.html` file and download to your local computer
    2. Open the file in your local browser

    Which process had required the highest memory usage?

    ??? abstract "Answer"
        `GENOTYPE` used 0.8-1.6 GB memory

## 2.2.2 Resource reporting with `nextflow log` 

Recall we generated a trace file in lesson [1.8.3](../part1/01_8_nfcore_config.md). For the runs we have already executed, which lacked the trace config option or command line flag, we can generate a trace report using the `nextflow log` command. This command is not identical to the trace file, but provides similar information about resource usage for each process in a tabular format. Since we don't have a trace file yet, we will first perform resource tracing using the `nextflow log` utility.

 
!!! example "Exercise: Review resource usage with nextflow log"

    1. Run the `nextflow log` command with no arguments to view a summary of previous runs:

    ```bash
    nextflow log
    ```

    ??? abstract "Output"

        | TIMESTAMP           | DURATION | RUN NAME          | STATUS | REVISION ID | SESSION ID                           | COMMAND |
        | ------------------- | -------- | ----------------- | ------ | ----------- | ------------------------------------ | ------- |
        | 2025-11-17 12:59:51 | 2m 24s   | naughty_bartik    | OK     | e34a5e5f9d  | 0ad50a9d-4e39-401e-ae46-65a1ee4e7933 | ...     |
        | 2025-11-17 13:07:26 | 2m 13s   | jolly_stonebraker | OK     | e34a5e5f9d  | fecd0d8c-0b9d-4cd6-9409-874ab4b2f976 | ...     |                   

    This information is extracted from file saved in the work directory called `.nextflow.log`, `.nextflow.log.1` etc. These files are renamed each run from the same directory, so that latest run is always `.nextflow.log` and the highest-numbered log is the oldest. In order to perform resource tracing using this method, it is vital to not delete these logs!
    
    2. Include the `-list-fields` flag to view all of the available fields for this utility:

    ```bash
    nextflow log -list-fields
    ```

    That's quite a few! Just a handful of fields are shown by default, but there are >40 that can be optionally displayed.

    3. Extract some specific fields for a recent run. Choose any fields you like, but we suggest focusing on those that provide relevant run info. Ensure to substitute `<run name>` in the command below with the actual run name from one of your runs. 

    ```groovy
    nextflow log <run name> -f name,status,exit,realtime,cpus,pcpu,memory,pmem,rss
    ```

!!! question "Question"
    
    Which process had the highest `realtime`?

    ??? abstract "Answer"

        `GENOTYPE`: ~ 30 seconds. 

We have now used the `nextflow log` utility to easily extract a snapshot of the resource usage for all of our processes, including the CPU and memory requested (`cpus`, `memory`), how much of it was utilised (`%cpu`, `%mem`), and how long it ran (`duration`, `walltime`).

Next we will make this automated and reproducible by adding the `trace{}` directive to our configuration file. This will save all this information to a text file `runInfo/` to keep our launch directory neat.

## 2.2.3 Resource reporting with trace files 

[Trace reports](https://nextflow.io/docs/latest/reports.html#trace-file) can be customised to provide detailed records of each process executed within a pipeline. As we learnt earlier, this file is generated by adding the `-with-trace` flag to your execution command or including the `trace{}` directive in your custom configuration file. While the HTML reports we looked at earlier focus on visual summaries of the whole run, trace text files give you raw process-level data that you can slice, filter, and analyse as needed: perfect for working on the HPC! 


The Nextflow `trace` tool has a lot of functionality. Like `nextflow log`, not all fields are displayed by default. For a full list and description of trace fields available, users should refer to the Nextflow [trace documentation](https://nextflow.io/docs/latest/reports.html#trace-file). Note that the `log` and `trace` fields are not identical, but do have many shared fields. 

Next we will extract the same fields with trace as we did with `nextflow log`, but note the small difference in field names for percentage CPU usage (`%cpu` vs `pcpu`) and percentage memory usage (`%mem` vs `pmem`). As always in bioinformatics - read the docs! 


!!! example "Exercise: Enable trace reporting"

    1. Add the following to your `custom.config` file: 

    ```
    trace {
        enabled = true 
        overwrite = false 
        file = "./runInfo/trace-${params.timestamp}.txt"
        fields = 'name,status,exit,realtime,cpus,%cpu,memory,%mem,rss'
    }
    ```

    2. Save the config, then run your workflow again: 

    ```bash
    ./run.sh
    ```

    3. View your trace file:

    ```bash
    cat runInfo/trace-*.txt
    ```

    You should see something like this:

    === "Gadi (PBS)"

        | name                       | status    | exit | duration | realtime | cpus | %cpu  | memory | %mem | peak_rss |
        | -------------------------- | --------- | ---- | -------- | -------- | ---- | ----- | ------ | ---- | -------- |
        | ALIGN (1)                  | COMPLETED | 0    | 29.6s    | 1s       | 1    | 93.7% | 4 GB   | 0.0% | 95.8 MB  |
        | FASTQC (fastqc on NA12877) | COMPLETED | 0    | 34.5s    | 5s       | 1    | 76.2% | 4 GB   | 0.1% | 286.6 MB |
        | GENOTYPE (1)               | COMPLETED | 0    | 59.9s    | 45s      | 1    | 97.6% | 4 GB   | 0.5% | 950.3 MB |
        | JOINT_GENOTYPE (1)         | COMPLETED | 0    | 34.8s    | 16s      | 1    | 93.3% | 4 GB   | 0.3% | 508.8 MB |
        | STATS (1)                  | COMPLETED | 0    | 19.9s    | 0ms      | 1    | 73.4% | 4 GB   | 0.0% | 3.1 MB   |
        | MULTIQC                    | COMPLETED | 0    | 29.9s    | 4.7s     | 1    | 79.5% | 4 GB   | 0.0% | 97.2 MB  |

    === "Setonix (Slurm)"

        | name                       | status    | exit | duration | realtime | cpus | %cpu   | memory | %mem | peak_rss |
        | -------------------------- | --------- | ---- | -------- | -------- | ---- | ------ | ------ | ---- | -------- |
        | FASTQC (fastqc on NA12877) | COMPLETED | 0    | 13.8s    | 4s       | 1    | 135.0% | 2 GB   | 0.1% | 240.6 MB |
        | ALIGN (1)                  | COMPLETED | 0    | 13.8s    | 2s       | 1    | 100.1% | 2 GB   | 0.0% | 98.2 MB  |
        | GENOTYPE (1)               | COMPLETED | 0    | 39.9s    | 28s      | 1    | 164.8% | 2 GB   | 0.5% | 1.4 GB   |
        | JOINT_GENOTYPE (1)         | COMPLETED | 0    | 19.2s    | 8s       | 1    | 204.2% | 2 GB   | 0.2% | 466 MB   |
        | STATS (1)                  | COMPLETED | 0    | 14.9s    | 1s       | 1    | 45.2%  | 2 GB   | 0.0% | 2 MB     |
        | MULTIQC                    | COMPLETED | 0    | 19.9s    | 5.3s     | 1    | 62.4%  | 2 GB   | 0.0% | 78.6 MB  |

## 2.2.4 Trace file enhancement

When things go wrong with process execution, we often need to trawl through the work directory to view `.command.err` files or view the job with `qstat` or `sacct` to debug the error. With trace, we can add two very useful fields that help us quickly locate the work directory and the scheduler job ID for each process.

!!! example "Exercise: Add work directory and job ID to trace file"

    1. View the documentation on [trace fields](https://www.nextflow.io/docs/latest/reports.html#trace-fields), and locate the two fields that provide the:

        - `work/` directory path where the task was executed
        - job ID executed by a grid engine (i.e. the scheduler)

    2. Append these two field names to `trace.fields` in `custom.config`

    ??? question "Hint"

        ```
        trace {
            enabled = true
            overwrite = false
            file = "./runInfo/trace-${params.timestamp}.txt"
            fields = 'name,status,exit,duration,realtime,cpus,%cpu,memory,%mem,peak_rss,workdir,native_id'
        }
        ```
    3. Add the Nextflow `-resume` flag to your run script to avoid re-running completed processes:

    ??? question "Hint"

        ```
        nextflow run main.nf -profile slurm -c config/custom.config -resume
        ```
    

    4. Save the run script config and run the workflow:

    ```bash
    ./run.sh
    ```

    5. View the newly generated trace file inside the `runInfo/` folder

    ```bash
    cat runInfo/trace-<newest-timestamp>.txt
    ```

    6. Now view the previously generated trace file:

    ```bash
    cat runInfo/trace-<newest-timestamp>.txt
    ```

Note that despite running with `-resume`, the trace file is still generated afresh for the entire run, including full resource usage details as if the process was executed again. This is very handy for benchmarking and tracking resource usage over multiple runs, as the full workflow details are always included in the latest trace regardless of whether a task was cached or re-executed.

It is up to you how you want to configure your traces for your own pipelines and how much added information you require. The suggestions we have demonstrated in this lesson are a good starting point for most bioinformatics workflows.

## 2.2.5 Summary

In this section, we have enabled Nextflow's execution report and timeline features to generate HTML summary reports, explored how to extract custom information from the `nextflow log` command, and configured customised trace file reporting. 

These powerful profiling features can help you benchmark and debug your runs, monitor resource usage, and plan for upcoming compute needs. Applying these tools when you set up your own pipelines will help you build efficient and reliable workflows that make the best use of your HPC resources.

## 2.2.6 Code checkpoint

??? abstract "Show complete code at the end of Section 2.2"

    === "Gadi (PBS)"

        ```groovy title='config/custom.config'
        process {
            cpu = 1 // 'normalbw' queue = 128 GB / 28 CPU ~ 4.6 OR 9.1
            memory = 4.GB
        }

        // Name the reports according to when they were run
        params.timestamp = new java.util.Date().format('yyyy-MM-dd_HH-mm-ss')

        // Generate timeline-timestamp.html timeline report 
        timeline {
            enabled = true
            overwrite = false
            file = "./runInfo/timeline-${params.timestamp}.html"
        }

        // Generate report-timestamp.html execution report 
        report {
            enabled = true
            overwrite = false
            file = "./runInfo/report-${params.timestamp}.html"
        }

        trace {
            enabled = true 
            overwrite = false 
            file = "./runInfo/trace-${params.timestamp}.txt"
            fields = 'name,status,exit,realtime,cpus,%cpu,memory,%mem,rss,workdir,native_id'
        }
        ```

        
        ```groovy title="run.sh"
        #!/bin/bash

        module load nextflow/24.04.5
        module load singularity

        nextflow run main.nf -profile pbspro -c config/custom.config -resume
        ```

    === "Setonix (Slurm)"

        ```groovy title='config/custom.config'
        process {
            cpu = 1 // 'work' partition = 230 GB / 128 CPU ~ 1.8
            memory = 2.GB
        }

        // Name the reports according to when they were run
        params.timestamp = new java.util.Date().format('yyyy-MM-dd_HH-mm-ss')

        // Generate timeline-timestamp.html timeline report 
        timeline {
            enabled = true
            overwrite = false
            file = "./runInfo/timeline-${params.timestamp}.html"
        }

        // Generate report-timestamp.html execution report 
        report {
            enabled = true
            overwrite = false
            file = "./runInfo/report-${params.timestamp}.html"
        }

        trace {
            enabled = true 
            overwrite = false 
            file = "./runInfo/trace-${params.timestamp}.txt"
            fields = 'name,status,exit,realtime,cpus,%cpu,memory,%mem,rss,workdir,native_id'
        }
        ```

    
        ```groovy title="run.sh"
        #!/bin/bash

        module load nextflow/24.10.0
        module load singularity/4.1.0-slurm

        nextflow run main.nf -profile slurm -c config/custom.config -resume
        ```