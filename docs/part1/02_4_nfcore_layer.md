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

    We have four distinct processes that we want to fine-tune:

    - `TABIX_BGZIPTABIX_INTERVAL_COMBINED`
    - `GATK4_CREATESEQUENCEDICTIONARY`
    - `GATK4_MARKDUPLICATES`
    - `MULTIQC`

    From the trace file we received from the previous run of `sarek`, we saw that the processes were requesting between 1 and 6 CPUs, and up to 36 GB of memory for the `GATK4_CREATESEQUENCEDICTIONARY` process. For our example dataset, these values are overkill. Instead, we can get away with just 1 CPU and 1GB of memory for each task. We'll also give each task just 2 minutes to complete, which is more than enough time.

    !!! note "Some tools are greedy!"

        You might have seen from our previous trace file that the `GATK4_MARKDUPLICATES` processes used several gigabytes of memory each (as reported in the `rss` column). So how come we can give them just 1GB now? GATK is a bit of a greedy tool and will often expand to use up lots of memory if is allowed to, so the values reported by the trace file aren't necessarily representative of how much memory the tool really needs. Optimising resources for tools therefore requires a bit of trial and error.

    Let's translate this into the Nextflow configuration format:

    ```groovy title="config/custom.config"
    process {

        withName: 'TABIX_BGZIPTABIX_INTERVAL_COMBINED' {
            cpus = 1
            memory = 1.GB
            time = 2.min
        }

        withName: 'GATK4_CREATESEQUENCEDICTIONARY' {
            cpus = 1
            memory = 1.GB
            time = 2.min
        }

        withName: 'GATK4_MARKDUPLICATES' {
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

    Now that we have our custom configuration file created, we need to update our run script one final time and add the new file to the `-c` option:

    === "Gadi"

        ```bash title="run.sh" linenums="1" hl_lines="15"
        #!/bin/bash

        module load nextflow/24.04.5
        module load singularity

        nextflow run sarek/main.nf \
            --input ../data/bams/samplesheet.csv \
            --fasta ../data/ref/Hg38.subsetchr20-22.fasta \
            --fasta_fai ../data/ref/Hg38.subsetchr20-22.fasta.fai \
            --step markduplicates \
            --skip_tools baserecalibrator,mosdepth,samtools \
            --outdir results \
            --no_intervals true \
            --igenomes_ignore true \
            -c config/hpc.config,config/custom.config \
            -resume
        ```

    === "Setonix"

        ```bash title="run.sh" linenums="1" hl_lines="15"
        #!/bin/bash

        module load nextflow/24.10.0
        module load singularity/4.1.0-slurm

        nextflow run sarek/main.nf \
            --input ../data/bams/samplesheet.csv \
            --fasta ../data/ref/Hg38.subsetchr20-22.fasta \
            --fasta_fai ../data/ref/Hg38.subsetchr20-22.fasta.fai \
            --step markduplicates \
            --skip_tools baserecalibrator,mosdepth,samtools \
            --outdir results \
            --no_intervals true \
            --igenomes_ignore true \
            -c config/hpc.config,config/custom.config \
            -resume
        ```

    And now we're ready to re-run the pipeline!

    ```bash
    ./run.sh
    ```

    After a few minutes, the pipeline should finish. We can again inspect the trace file from the run to see how much memory was requested and used:

    ```bash
    # Your trace file will have a unique name based on the time it was run
    cat runInfo/trace-2025-11-18_14-15-16.txt
    ```

    ```console title="Output"
    name	status	exit	duration	realtime	cpus	%cpu	memory	%mem	rss
    NFCORE_SAREK:PREPARE_INTERVALS:TABIX_BGZIPTABIX_INTERVAL_COMBINED (no_intervals)	COMPLETED	0	58s	0ms	1	40.9%	1 GB	0.0%	3.1 MB
    NFCORE_SAREK:PREPARE_GENOME:GATK4_CREATESEQUENCEDICTIONARY (Hg38.subsetchr20-22.fasta)	COMPLETED	0	1m 12s	10s	1	78.7%	1 GB	0.2%	326.5 MB
    NFCORE_SAREK:SAREK:BAM_MARKDUPLICATES:GATK4_MARKDUPLICATES (test_sample3)	COMPLETED	0	1m 9s	9s	1	93.6%	1 GB	0.5%	602.8 MB
    NFCORE_SAREK:SAREK:BAM_MARKDUPLICATES:GATK4_MARKDUPLICATES (test_sample1)	COMPLETED	0	1m 19s	12s	1	81.9%	1 GB	0.4%	597.3 MB
    NFCORE_SAREK:SAREK:BAM_MARKDUPLICATES:GATK4_MARKDUPLICATES (test_sample2)	COMPLETED	0	1m 11s	8s	1	94.5%	1 GB	0.5%	600.2 MB
    NFCORE_SAREK:SAREK:MULTIQC	COMPLETED	0	1m 20s	8.7s	1	67.4%	1 GB	0.2%	420.9 MB
    ```

    We can see that most of the processes are now using a larger proportion of the memory assigned to them (in this case 300-600MB out of a total 1GB), so we are much more efficiently using our resources. We could probably fine-tune this even further, but we'd get diminishing returns and risk some samples failing.