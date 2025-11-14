# Pipeline monitoring and reporting

!!! info "Learning objectives"

    - Inspect workflow performance and resource usage utilising Nextflow's
    in-built monitoring features
    - Know which fields help determine efficiency and HPC resource usage
    - Compare the trade-offs between Nextflow's profiling features in comparison to unix tools such as `time` or `gprof`

In bioinformatics workflows, resource requirements are often not fixed. They can vary significantly depending on the size of input files, varying complexity of different genomic regions, the species you're working with, and sequencing format and depth. This means we can't assume CPU, memory, or time values will work the same for every sample. On HPC systems, where resources are shared and allocations may be charged, these differences matter. 

You cannot set and forget resource values for an entire workflow, you need to build in flexibility. This requires you to have visiblity over pipeline behaviour at the process level. Nextflow provides several monitoring and reporting tools that help you understand this behaviour including `nextflow log`, trace files, and reports. 

!!! note "Recall the trace file from Part 1"

    TODO add a reference to how we did this in part 1 and comment on how we will explore further here. Link out when we have page structure of Part 1 ready.

Now that our workflow is running without error on the scheduler, we will enable [Nextflow's reports](https://nextflow.io/docs/latest/reports.html). This allows us to view the resource usage of each process for our representative sample.


## 2.2.1 Using the execution report and timeline

Nextflow can produce an [execution report](https://nextflow.io/docs/latest/reports.html#execution-report) at the end of your workflow run that summarises all process execution metrics. Similarly, it can create a [execution timeline](https://nextflow.io/docs/latest/reports.html#execution-timeline). These reports summarise how your workflow ran, which processes were executed, and how long they took. They are very helpful during development, troubleshooting, and performance optimisation. These reports can be created when running the pipeline using the `-with-report` and `-with-timeline` flags or by adding the following to your configuration file:

```console
timeline {
    enabled = true
    overwrite = false
    file = "./runInfo/timeline.html"
    }

report {
    enabled = true
    overwrite = false
    file = "./runInfo/report.html"
        }
```

Both directives allow us to specify a custom file name and choose whether or not you overwrite the file for each run. Rather than overwriting the same report file every time we run the pipeline, we will add a small timestamp parameter that automatically labels each report with the exact date and time the workflow was launched. This makes it easier to track multiple runs, especially when you are iterating qucikly and comparing resource usage. 

!!! example "Exercise"

    Enable execution and timeline reports in our `custom.config` file. Copy the following into your `custom.config`:

    === "Gadi (PBSpro)"

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
  
    Run the workflow 
    ```bash
    ./run.sh
    ```
    Once the job completes, you should have a new folder called `runInfo`, where trace files are saved. Look at your first trace file:

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

We want these fields because they provide a great snapshot of the resource usage for all of our processes. These include the raw values for CPU and memory requested (`cpus`, `memory`), how much of it was utilised (`%cpu`, `%mem`), and how long it ran (`duration`, `walltime`).

Displaying all the usage information in a single file can avoid needing to benchmark each process manually, as we saw in Part 1.

The addition of timestamps and overwrite = false helps with benchmarking when you
need to compare settings before vs. after changing a configuration setting. 

We save this all into `runInfo/` to keep our launch directory neat.

## Customising the trace file

Currently, the trace file reports on the resources used per task. However,
when a tasks errors or produces an unexpected result, it is recommended to
view the work directory, or view the job run information with `qstat` or
`sacct`.

!!! example "Exercise"

    1. View the documentation on [trace fields](https://www.nextflow.io/docs/latest/reports.html#trace-fields)
    2. Locate the two field names that provides the:

        - directory path where the task was executed
        - the job ID when executed by a grid engine (scheduler)

    3. Append these to `trace.fields` in `nextflow.config`

    ??? question "Hint"

        ```groovy hl_lines="5"
        trace {
            enabled = true
            overwrite = false
            file = "./runInfo/trace-${params.timestamp}.txt"
            fields = 'name,status,exit,duration,realtime,cpus,%cpu,memory,%mem,peak_rss,workdir,native_id'
        }
        ```

    4. Save, and run the workflow:
    ```bash
    ./run.sh
    ```

    5. View the newly generate trace file under the `runInfo/` folder

    _Bonus: feel free to include several different fields and re-run your pipeline. However, ensure the fields between `name` through to `peak_rss` are included before proceeding to the next lessons._

These added fields help you track down the scheduler job and work directories
for debugging. It is up to you how you want to configure your traces for your own pipelines and how much added information you require.

## Summary

Nextflow's profiling features are powerful automations that scale nicely when
you have many processes in a pipeline. This saves the need for manual
benchmarking such as using `time`.

??? Code

    === "Gadi (PBS)"

        ```groovy title='config/custom.config' hl_lines="6-31"
        process {
            cpu = 1 // 'normalbw' queue = 128 GB / 28 CPU ~ 4.6 OR 9.1
            memory = 4.GB
        }

        params.timestamp = new java.util.Date().format('yyyy-MM-dd_HH-mm-ss')

        trace {
             enabled = true
             overwrite = false
             file = "./runInfo/trace-${params.timestamp}.txt"
             fields = 'name,status,exit,duration,realtime,cpus,%cpu,memory,%mem,peak_rss,workdir,native_id'
         }

        timeline {
            enabled = true
            overwrite = false
            file = "./runInfo/timeline-${params.timestamp}.html"
        }

        report {
            enabled = true
            overwrite = false
            file = "./runInfo/report-${params.timestamp}.html"
        }
    
        dag {
            enabled = true
            overwrite = false
            file = "./runInfo/dag-${params.timestamp}.html"
        }
        ```

    === "Setonix (Slurm)"

        ```groovy title='config/custom.config' hl_lines="6-31"
        process {
            cpu = 1 // 'work' partition = 230 GB / 128 CPU ~ 1.8
            memory = 2.GB
        }

        params.timestamp = new java.util.Date().format('yyyy-MM-dd_HH-mm-ss')

        trace {
             enabled = true
             overwrite = false
             file = "./runInfo/trace-${params.timestamp}.txt"
             fields = 'name,status,exit,duration,realtime,cpus,%cpu,memory,%mem,peak_rss,workdir,native_id'
         }

        timeline {
            enabled = true
            overwrite = false
            file = "./runInfo/timeline-${params.timestamp}.html"
        }

        report {
            enabled = true
            overwrite = false
            file = "./runInfo/report-${params.timestamp}.html"
        }
    
        dag {
            enabled = true
            overwrite = false
            file = "./runInfo/dag-${params.timestamp}.html"
        }
        ```