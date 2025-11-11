# Scale to multiple samples

!!! info "Learning objectives"

    - Apply an optimised pipeline across multiple samples
    - Understand and compare sample-level and within-sample parallelism
    - Recall best practices for running multi-sample data with samplesheets
    
Now it has run successfully and efficiently on one, run on multiple samples in parallel. Explain sample-level paralellism. Refer back to smarter not harder lesson in HPC foundations. 

Compare to things like removing the need for [loop optimisations](https://pawsey.atlassian.net/wiki/spaces/US/pages/51925998/Loop+Optimisations)

This demonstrates another form of parallelisation, by sample.

!!! example "Exercise"

    Update your `run.sh` and change the `--samplesheet` param to use `samplesheets_full.csv`:

    === "Gadi (PBS)"

        ```groovy title="run.sh"
        #!/bin/bash

        module load nextflow/24.04.5
        module load singularity

        nextflow run main.nf -profile pbspro --pbspro_account vp91 -c conf/custom.config --samplesheet "samplesheet_full.csv"
        ```

    === "Setonix (Slurm)"

        ```groovy title="run.sh"
        #!/bin/bash

        module load nextflow/24.10.0
        module load singularity/4.1.0-slurm

        nextflow run main.nf -profile slurm --slurm_account courses01 -c conf/custom.config --samplesheet "samplesheet_full.csv"
        ```


    Save your script and re-run!
    
    ```bash
    ./run.sh
    ```

TODO: Show output