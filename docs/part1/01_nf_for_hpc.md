# Running Nextflow on HPC

!!! info "Learning objectives"

    - blah 
    - blah 

## What happens when a workflow runs on HPC? 

- Tasks map to jobs
- Executors: local vs Slurm vs PBS
- Staging work directories and `.command.sh`
- Parallelisation through job arrays / process parallelism
- **Key takeaway:** Nextflow does not run your workflow — the scheduler does!