# Pipeline monitoring and reporting 

!!! info "Learning objectives"

    - Inspect workflow performance and resource usage utilising Nextflow's
    in-built monitoring features
    - Recognise the importance of benchmarking to request appropriate resources (e.g. CPUs, memory, time) when scheduling jobs on HPC
    - Know which fields help determine efficiency and HPC resource usage
    - Compare the trade-offs between Nextflow's profiling features in comparison to unix tools such as `time` or `gprof`

Once we get the workflow running without error on the scheduler, where can we optimise. 

!!! example "Exercise"

    Enable all trace reporting available, with default/minimal settings.
    We use the configured `trace` from Part 1.

    ```groovy title='nextflow.config'
     params.trace_timestamp = new java.util.Date().format('yyyy-MM-dd_HH-mm-ss')
     
     trace {
         enabled = true
         overwrite = false
         file = "./runInfo/trace-${params.trace_timestamp}.txt"
         fields = 'name,status,exit,duration,realtime,cpus,%cpu,memory,%mem,rss'
     }

    timeline { enabled = true }
    report { enabled = true }
    dag { enabled = true }
    ```

    Run the workflow using `run.sh`

Look at the report.html. Create a more informative trace file: https://www.nextflow.io/docs/latest/reports.html#trace-file.

Explain why we want these fields - tie in with benchmarking and HPC resource
allocation.

Addition of timestamp and overwrite = false - helps with benchmarking when you
need to compare settings before vs. after e.g. optimisation

## Customising the trace file

Currently, the trace file reports on the resources used per task. However,
when a tasks errors or produces an unexpected result, it is recommended to
view the work directory, or view the job run information with `qstat` or
`sacct`.

and select ones that are not yet included. Things such as the `workDir`bb

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
            file = "./runInfo/trace-${params.trace_timestamp}.txt"
            fields = 'name,status,exit,duration,realtime,cpus,%cpu,memory,%mem,rss,workdir,native_id'
        }
        ```

    4. Save, and run the workflow using `run.sh`
    5. View the newly generate trace file under the `runInfo/` folder

These added fields help you track down the scheduler job and work directories
for debugging.

## Summary

Nextflow's profiling features are powerful automations that scale nicely when
you have many processes in a pipeline. This saves the need for manual
benchmarking such as using `time`.
