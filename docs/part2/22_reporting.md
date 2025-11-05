# Pipeline monitoring and reporting 

!!! info "Learning objectives"

    - Inspect workflow performance and resource usage utilising Nextflow's
    in-built monitoring features
    - Recognise the importance of benchmarking to request appropriate resources (e.g. CPUs, memory, time) when scheduling jobs on HPC
    - Know which fields help determine efficiency and HPC resource usage
    - Compare the trade-offs between Nextflow's profiling features in comparison to unix tools such as `time` or `gprof`

Once we get the workflow running without error on the scheduler, where can we optimise. 

Look at the report.html. Create a more informative trace file: https://www.nextflow.io/docs/latest/reports.html#trace-file

!!! example "Exercise"

    TODO: Enable all trace reporting available, with default settings

!!! example "Exercise"
    
    Add the following to the end of your `nextflow.config`

    ```groovy title='nextflow.config'
     params.trace_timestamp = new java.util.Date().format('yyyy-MM-dd_HH-mm-ss')
     
     trace {
         enabled = true
         overwrite = false
         file = "./runInfo/trace-${params.trace_timestamp}.txt"
         fields = 'name,status,exit,duration,realtime,cpus,%cpu,memory,%mem,rss'
     }
    ```
   

!!! example "Exercise"

    TODO Run pipeline, view reports, particularly trace

!!! example "Exercise"

    TODO Customise trace. Run the pipeline again

