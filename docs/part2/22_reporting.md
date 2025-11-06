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

!!! example "Exercise"

    TODO View trace fields, add some stuff that could be useful. e.g. jobid so you don't have to grep .nextflow.log

    TODO Run pipeline, view reports, particularly trace - what has changed?

    === "Gadi (PBS)"

TODO compare the trade-offs between Nextflow's profiling features in comparison to unix tools such as `time` or `gprof`
