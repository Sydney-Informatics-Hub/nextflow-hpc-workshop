# 2.4 Layering Nextflow configurations

!!! info "Learning objectives"

    - Learn how Nextflow configurations can be layered on top of one another
    - Understand the configuration priorities
    - Learn how to fine-tune processes using the `withName` directive
    - Understand when you would want to fine-tune process requirements

We have now managed to configure and run an nf-core pipeline, but we have also seen that for our test dataset it isn't very well optimised: by default, the `sarek` pipeline requests many more CPUs and much more RAM than necessary for our purposes. This problem isn't just relevant to small test datasets; sometimes you might find when running a large dataset that the pipeline hasn't been optimised quite as well as you would like and is requesting fewer resources than you need. In these cases, what we would like to do is have a fine level of control over each process and the resources they request.

We have already seen how we can define a custom configuration file and layer it over the top of the default configuration using the `-c` option to Nextflow. We have also seen how the `sarek` pipeline defines the resources required by many of its processes within the `conf/base.config` file using `withLabel` and `withName` directives. In this final section of today's workshop, we will try to optimise the processes that we are running to more efficiently use the HPC resources by defining a new custom configuration file.
    
## 2.4.1 Configuration priorities

Before we get started, it's important to understand how Nextflow prioritises configuration files. Because we can provide configuration information at various levels and using multiple files, it is possible for some options in these places to overlap, and Nextflow needs to know which ones to give precedence to. Referring to the [Nextflow configuration documentation](https://www.nextflow.io/docs/latest/config.html), configuration files are prioritised in the following order, from lowest to highest priority:

1. The config file $HOME/.nextflow/config
2. The config file nextflow.config in the project directory (i.e. at the same level as the `main.nf` script)
3. The config file nextflow.config in the launch directory (i.e. the directory in which you run `nextflow run ...`)
4. Config files specified using the `-c <config-files>` option

Furthermore, when using the `-c` option, multiple configuration files can be provided, separated by commas, and are prioritised from lowest to highest in the order they are specified.

!!! example "Example: A very simple layered configuration"

    Consider the following (very basic) Nextflow file:

    ```groovy title="example.nf"
    params.value = "hello"

    workflow {
        println params.value
    }
    ```

    If we run this 'workflow', it will print the contents of the `value` parameter, i.e. "hello":

    ```bash
    nextflow run example.nf
    ```

    ```console title="Output"

    N E X T F L O W   ~  version 24.10.5

    Launching `example.nf` [tender_kay] DSL2 - revision: 573919f401

    hello
    ```

    Now suppose we create a `nextflow.config` file and set `value` to something different:

    ```groovy title="nextflow.config"
    params.value = "bye"
    ```

    Now, the workflow will print "bye"

    ```bash
    nextflow run example.nf
    ```

    ```console title="Output"

    N E X T F L O W   ~  version 24.10.5

    Launching `example.nf` [tender_kay] DSL2 - revision: 573919f401

    bye
    ```

    If we create another config file, define `params.value` in there, and layer it on top, that value will be used:

    ```groovy title="custom.config"
    params.value = "seeya"
    ```

    ```bash
    nextflow run example.nf -c custom.config
    ```

    ```console title="Output"

    N E X T F L O W   ~  version 24.10.5

    Launching `example.nf` [tender_kay] DSL2 - revision: 573919f401

    seeya
    ```

    And if we create a second custom config, define yet another value for `params.value`, and layer it on top as well, that value will be used:

    ```groovy title="another_custom.config"
    params.value = "ciao"
    ```

    ```bash
    nextflow run example.nf -c custom.config,another_custom.config
    ```

    ```console title="Output"

    N E X T F L O W   ~  version 24.10.5

    Launching `example.nf` [tender_kay] DSL2 - revision: 573919f401

    ciao
    ```

Process directives, such as CPU and memory requirements, can be configured in a number of ways, and these too are evaluated in a particular order by Nextflow. Briefly, they are prioritised in the following order, from lowest to highest priority:

1. Default process configuration settings in the configuration files (e.g. `process.cpu = 1`)
2. Process directives defined in the process definition
3. Process configuration settings within a matching `withLabel` selector
4. Process configuration settings within a matching `withName` selector

!!! example "Example: Configuring default resources and using process selectors"

    Consider the following `process {}` scope within a configuration file:

    ```groovy
    process {
        cpus = 4
        withLabel: hello { cpus = 8 }
        withName: bye { cpus = 16 }
    }
    ```

    This configuration will have the following consequences:

    - By default, all processes will be given 4 CPUs, unless their process definitions contain a `cpus` directive
    - Any process given the `label` "hello" will instead be given 8 CPUs
    - Any process named "bye" will be given 16 CPUs

## 2.4.2 Optimising `nf-core/sarek` for our data

!!! example "Exercise: Fine-tune `nf-core/sarek`"

    Start by creating a new blank file within the `config/` folder called `custom.config` and open it up in VSCode.

    We have eight distinct processes that we want to fine-tune:

    - `TABIX_BGZIPTABIX_INTERVAL_COMBINED`
    - `FASTQC`
    - `FASTP`
    - `BWAMEM1_MEM`
    - `MERGE_BAM`
    - `INDEX_MERGE_BAM`
    - `BAM_TO_CRAM_MAPPING`
    - `MULTIQC`

    From the trace file we received from the previous run of `sarek`, we saw that the processes were requesting between 1 and 24 CPUs, and up to 30 GB of memory for the `GATK4_CREATESEQUENCEDICTIONARY` process. For our example dataset, these values are overkill. Instead, we can get away with just 1-2 CPUs and 1-2GB of memory for each task. We'll also give each task just 2 minutes to complete, which is more than enough time.

    Let's translate this into the Nextflow configuration format:

    ```groovy title="config/custom.config"
    process {

        withName: 'TABIX_BGZIPTABIX_INTERVAL_COMBINED' {
            cpus = 1
            memory = 1.GB
            time = 2.min
        }

        withName: 'FASTQC' {
            cpus = 1
            memory = 1.GB
            time = 2.min
        }

        withName: 'FASTP' {
            cpus = 2
            memory = 2.GB
            time = 2.min
        }

        withName: 'BWAMEM1_MEM' {
            cpus = 1
            memory = 1.GB
            time = 2.min
        }

        withName: 'MERGE_BAM' {
            cpus = 1
            memory = 1.GB
            time = 2.min
        }

        withName: 'INDEX_MERGE_BAM' {
            cpus = 1
            memory = 1.GB
            time = 2.min
        }

        withName: 'BAM_TO_CRAM_MAPPING' {
            cpus = 1
            memory = 1.GB
            time = 2.min
        }

        withName: 'MULTIQC' {
            cpus = 1
            memory = 1.GB
            time = 2.min
        }

    }
    ```

    We are now giving most tasks 1 CPU and 1 GB of memory, except for `FASTP`, which we'll give 2 CPUs and 2 GB of memory. When running `sarek`, giving `FASTP` multiple CPUs causes it to also split the FASTQ up, thereby implementing a scatter-gather pattern. So, in this case, we are scattering our input FASTQ into two chunks and aligning them individually.

    Now that we have our custom configuration file created, we need to update our run script once again and add the new file to the `-c` option:

    === "Gadi (PBS)"

        ```bash title="run.sh" linenums="1" hl_lines="15"
        #!/bin/bash

        module load nextflow/24.04.5
        module load singularity

        nextflow run sarek/main.nf \
            --input ../data/fqs/samplesheet.single.csv \
            --fasta ../data/ref/Hg38.subsetchr20-22.fasta \
            --fasta_fai ../data/ref/Hg38.subsetchr20-22.fasta.fai \
            --dict ../data/ref/Hg38.subsetchr20-22.dict \
            --bwa ../data/ref \
            --step mapping \
            --skip_tools markduplicates,baserecalibrator,mosdepth,samtools \
            --outdir results \
            --no_intervals true \
            --igenomes_ignore true \
            -c config/gadi.config,config/custom.config \
            -resume
        ```

    === "Setonix (Slurm)"

        ```bash title="run.sh" linenums="1" hl_lines="15"
        #!/bin/bash

        module load nextflow/24.10.0
        module load singularity/4.1.0-slurm

        nextflow run sarek/main.nf \
            --input ../data/fqs/samplesheet.single.csv \
            --fasta ../data/ref/Hg38.subsetchr20-22.fasta \
            --fasta_fai ../data/ref/Hg38.subsetchr20-22.fasta.fai \
            --dict ../data/ref/Hg38.subsetchr20-22.dict \
            --bwa ../data/ref \
            --step mapping \
            --skip_tools markduplicates,baserecalibrator,mosdepth,samtools \
            --outdir results \
            --no_intervals true \
            --igenomes_ignore true \
            -c config/setonix.config,config/custom.config \
            -resume
        ```

    And now we're ready to re-run the pipeline!

    ```bash
    ./run.sh
    ```

    After a few minutes, the pipeline should finish. We can again inspect the trace file from the run to see how much CPU and memory was requested and used:

    ```bash
    # Your trace file will have a unique name based on the time it was run
    cut -f 1,6,8,10 runInfo/trace-2025-11-18_14-15-16.txt
    ```

    ```console title="Output"
    name    cpus    memory  peak_rss
    NFCORE_SAREK:SAREK:FASTP (test_sample1-all)     2       2 GB    2 MB
    NFCORE_SAREK:SAREK:FASTQC (test_sample1-all)    1       1 GB    218.8 MB
    NFCORE_SAREK:PREPARE_INTERVALS:TABIX_BGZIPTABIX_INTERVAL_COMBINED (no_intervals)        1       1 GB    2 MB
    NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:BWAMEM1_MEM (test_sample1)  1       1 GB    2 MB
    NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON:BWAMEM1_MEM (test_sample1)  1       1 GB    2 MB
    NFCORE_SAREK:SAREK:BAM_MERGE_INDEX_SAMTOOLS:MERGE_BAM (test_sample1)    1       1 GB    4 MB
    NFCORE_SAREK:SAREK:BAM_MERGE_INDEX_SAMTOOLS:INDEX_MERGE_BAM (test_sample1)      1       1 GB    2 MB
    NFCORE_SAREK:SAREK:BAM_TO_CRAM_MAPPING (test_sample1)   1       1 GB    18.3 MB
    NFCORE_SAREK:SAREK:MULTIQC      1       1 GB    686.6 MB
    ```

    We can see that `FASTQC` and `MULTIQC` are now using a much larger proportion of the memory assigned to them (in this case 200-700MB out of a total 1GB), so we are much more efficiently using our resources for these jobs. The other jobs are still only using a few MB of memory to run, but for such small jobs there isn't too much utility in optimising these any further; from a cost perspective, these jobs are already under-utilising memory per CPU, so there would be no benefit to reducing the request furhter; and we could also start running into failures due to fluctuations in the memory used between runs.

## 2.4.3 Scaling up to multiple samples

    Now that we have a fully-functioning run script and custom configuration, we can try scaling up to multiple samples.

    !!! example "Exercise: Run mapping on multiple samples"

        Update the `run.sh` script to use the full samplesheet with all three test samples:

        === "Gadi (PBS)"

            ```bash title="run.sh"
            #!/bin/bash

            module load nextflow/24.04.5
            module load singularity

            nextflow run sarek/main.nf \
                --input ../data/fqs/samplesheet.fq.csv \
                --fasta ../data/ref/Hg38.subsetchr20-22.fasta \
                --fasta_fai ../data/ref/Hg38.subsetchr20-22.fasta.fai \
                --dict ../data/ref/Hg38.subsetchr20-22.dict \
                --bwa ../data/ref \
                --step mapping \
                --skip_tools markduplicates,baserecalibrator,mosdepth,samtools \
                --outdir results \
                --no_intervals true \
                --igenomes_ignore true \
                -c config/gadi.config \
                -resume
            ```

        === "Setonix (Slurm)"

            ```bash title="run.sh"
            #!/bin/bash

            module load nextflow/24.10.0
            module load singularity/4.1.0-slurm

            nextflow run sarek/main.nf \
                --input ../data/fqs/samplesheet.fq.csv \
                --fasta ../data/ref/Hg38.subsetchr20-22.fasta \
                --fasta_fai ../data/ref/Hg38.subsetchr20-22.fasta.fai \
                --dict ../data/ref/Hg38.subsetchr20-22.dict \
                --bwa ../data/ref \
                --step mapping \
                --skip_tools markduplicates,baserecalibrator,mosdepth,samtools \
                --outdir results \
                --no_intervals true \
                --igenomes_ignore true \
                -c config/setonix.config \
                -resume
            ```

        Go ahead and re-run the script. The `-resume` flag will mean that the previously-run tasks for the first sample (`FASTQC`, `FASTP`, `BWAMEM1_MEM`, etc.) will not be re-run, but instead their outputs will be reused. Only the new samples will be run through these processes. `MULTIQC` will re-run at the end of the pipeline as it needs to summarise the results from all three samples.
