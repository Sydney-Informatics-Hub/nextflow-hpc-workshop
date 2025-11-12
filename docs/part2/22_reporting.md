# Pipeline monitoring and reporting

!!! info "Learning objectives"

    - Inspect workflow performance and resource usage utilising Nextflow's
    in-built monitoring features
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
         fields = 'name,status,exit,duration,realtime,cpus,%cpu,memory,%mem,peak_rss'
     }

    // TODO: ADD path to runInfo
    timeline { enabled = true }
    report { enabled = true }
    dag { enabled = true }
    ```

    Run the workflow using
    ```bash
    ./run.sh
    ```
    Once the job completes, you should have a new folder called `runInfo`, where trace files are saved. Look at your first trace file:

    ```bash
    cat runInfo/trace-*.txt
    ```

    You should see something like this:

    | name                       | status     | exit | duration | realtime | cpus | %cpu   | memory | %mem | rss     |
    |----------------------------|------------|------|----------|----------|------|--------|--------|------|---------|
    | FASTQC (fastqc on NA12877) | COMPLETED  | 0    | 29.7s    | 4s       | 1    | 138.4% | 2 GB   | 0.1% | 251 MB  |
    | ALIGN (1)                  | COMPLETED  | 0    | 29.7s    | 1s       | 1    | 99.8%  | 2 GB   | 0.0% | 95 MB   |
    | GENOTYPE (1)               | COMPLETED  | 0    | 59.9s    | 30s      | 1    | 163.0% | 2 GB   | 0.5% | 1.3 GB  |
    | JOINT_GENOTYPE (1)         | COMPLETED  | 0    | 29.5s    | 7s       | 1    | 230.2% | 2 GB   | 0.1% | 400.9 MB|
    | STATS (1)                  | COMPLETED  | 0    | 29.8s    | 0ms      | 1    | 117.5% | 2 GB   | 0.0% | 2 MB    |
    | MULTIQC                    | COMPLETED  | 0    | 29.9s    | 4.3s     | 1    | 79.0%  | 2 GB   | 0.0% | 83.3 MB |

Explain why we want these fields - tie in with benchmarking and HPC resource
allocation.

Addition of timestamp and overwrite = false - helps with benchmarking when you
need to compare settings e.g. before vs. after configuration.

All saved into it's own folder for neatness

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
            file = "./runInfo/trace-${params.trace_timestamp}.txt"
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
for debugging.

## Summary

Nextflow's profiling features are powerful automations that scale nicely when
you have many processes in a pipeline. This saves the need for manual
benchmarking such as using `time`.
