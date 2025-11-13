# Pipeline monitoring and reporting

!!! info "Learning objectives"

    - Inspect workflow performance and resource usage utilising Nextflow's
    in-built monitoring features
    - Know which fields help determine efficiency and HPC resource usage
    - Compare the trade-offs between Nextflow's profiling features in comparison to unix tools such as `time` or `gprof`

Once we get the workflow running without error on the scheduler, we need to enable Nextflow's reporting and monitoring functions. This allows us to view the resource requirements that each process uses, on our representative sample.

!!! example "Exercise"

    Enable all trace reporting available, with default/minimal settings.
    We use the configured `trace` from Part 1.

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
             fields = 'name,status,exit,duration,realtime,cpus,%cpu,memory,%mem,peak_rss'
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
             fields = 'name,status,exit,duration,realtime,cpus,%cpu,memory,%mem,peak_rss'
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

        TODO

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
