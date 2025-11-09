# 2.1 - nf-core

!!! info "Learning objectives"

    - Understand what nf-core is
    - Identify where to find and access nf-core pipelines
    - Recognise the basic structure and key components of an nf-core workflow

We have now seen how a Nextflow pipeline can be configured to run on an HPC. In tomorrow's section of the workshop, we will further explore optimising a custom Nextflow pipeline to efficiently utilise HPC resources. However, whenever you are considering building a workflow, it is always important to check whether a suitable tool already exists - after all, the goal of Nextflow is to build reproducible workflows, and we shouldn't re-invent the wheel if we don't have to! For the rest of this session, we will be looking at the **nf-core** project, which aims to address this very issue and provide a collection of bioinformatics Nextflow pipelines.

## 2.1.1 What is nf-core?

![nf-core logo](../assets/nf-core-logo.png)

nf-core is a community-driven effort to develop and curate open-source bioinformatics workflows built with Nextflow.

The project has a standardised set of best practices, guidelines, and templates for building modular, scalable, and portable bioinformatics workflows. Every workflow is well-documented, tested, and designed to work across multiple platforms, including cloud and HPC.

The key Features of nf-core workflows are:

- Documentation: nf-core workflows have extensive documentation covering installation, usage, and description of output files to ensure that you won’t be left in the dark.
- CI Testing: Every time a change is made to the workflow code, nf-core workflows use continuous-integration testing to ensure that nothing has broken.
- Stable Releases: nf-core workflows use GitHub releases to tag stable versions of the code and software, making workflow runs totally reproducible.
- Packaged software: Pipeline dependencies are automatically downloaded and handled using Docker, Singularity, Conda, or other software management tools. There is no need for any software installations.
- Portablility and reproducibility: nf-core workflows follow best practices to ensure maximum portability and reproducibility. The large community makes the workflows exceptionally well-tested and easy to execute.
- Cloud-ready: nf-core workflows are tested on AWS after every major release. You can even browse results live on the website and use outputs for your own benchmarking.

nf-core is published in Nature Biotechnology: [Nat Biotechnol 38, 276–278 (2020). Nature Biotechnology](https://www.nature.com/articles/s41587-020-0439-x)

## 2.1.2 Where to find nf-core pipelines

The nf-core website - [https://nf-co.re](https://nf-co.re) - hosts a list of all of the current nf-core pipelines, as well as their documentation, information for developers and users, and links to community forums, training sessions, and more.

The full list of available piplines are available at [https://nf-co.re/pipelines/](https://nf-co.re/pipelines/), and at the time of writing, this includes 84 stable workflows, 44 more that are under development, and 12 that have been archived.

Each workflow has a dedicated page that includes expansive documentation that is split into 6 sections:

- Introduction: An introduction and overview of the workflow
- Usage: Documentation and descriptions of how to execute the workflow
- Parameters: Documentation for all workflow parameters
- Output: Descriptions and examples of the expected output files
- Results: Example output files generated from the full test dataset on AWS
- Releases: Workflow version history

All of the workflow code is hosted on GitHub under the `nf-core` project. For example, today we will be working with the `sarek` pipeline, which is hosted on GitHub at [https://github.com/nf-core/sarek](https://github.com/nf-core/sarek). Helpfully, unless you are actively developing workflow code yourself, you usually won’t need to clone the workflow code from GitHub; instead, you can use Nextflow’s built-in functionality to pull a workflow:

```bash
nextflow pull nf-core/<pipeline>
```

Nextflow's `run` command will also automatically pull the workflow if it was not already available locally:

```bash
nextflow run nf-core/<pipeline>
```

By default, Nextflow will pull the default git branch of the pipeline unless a specific version is specified with the `-revision` or `-r` flag.

**Note** that in today's workshop, for full transparency and consistency, we **will** be pulling the pipeline directly from GitHub rather than using the shortcut method above.

## 2.1.5 Introducing today's pipeline: `nf-core/sarek`

As mentioned above, the rest of today's workshop will be focussing on the `nf-core/sarek` pipeline. This is a large workflow dedicated to performing variant calling on genome sequencing data. It is highly configurable and incoroprates a wide variety of tools and methods for detecting both germline and somatic variants.

![The structure of the nf-core/sarek pipeline](../assets/sarek_subway.png)

The pipeline consists of three main stages:

- Pre-processing: This stage takes in sequencing FASTQ files and performs quality control, read trimming, genome alignment, marking of duplicate reads, and recalibrating base quality scores.
- Variant calling: This is the main workhorse of the pipeline and is the part that actually detects germline and/or somatic variants in the sequencing data. It uses several tools to achieve this, including `deepvariant`, `bcftools`, and the `GATK` suite.
- Annotation: This stage takes the variants that were called and annotates them using public databases. Annotations include the predicted effects of the variants, such as potential frameshift or stop-gain mutations, and their clinical relevance.

As you can imagine, this is a very complex pipeline with lots of options and parameters, and could take a long time to run in its default mode. Luckily, it also includes a wide range of options for selecting which parts of the pipeline to run or skip. We will be making heavy use of these options today so that we can just run a small section of the workflow. We will be providing the pipeline with pre-aligned BAM files and just running the `markduplicates` part of the pre-processing stage.

![The markduplicates subsection of the nf-core/sarek pipeline](../assets/sarek_markduplicates.png)

!!! note "An explanation of the markduplicates workflow"

    For those that are unfamiliar with variant calling, be assured that no prior knowledge of variant calling is required for this workshop. Our focus today is on running and configuring an existing nf-core pipeline, not on what it is doing under the hood. Nevertheless, it is useful to have a large-picture overview of the workflow that we will be running and what each of the files are for.

    Our input to the pipeline will be a BAM file. This is a binary version of the SAM file format (which stands for Sequence Alignment Map). BAM/SAM files are created by an **alignment** process that takes in FASTQ files (containing the raw sequencing reads and their quality scores as determined by the sequencing instrument) and maps the reads within to a reference genome. The output BAM/SAM file contains information about each sequencing read, including what the sequence is, what the quality of each base is, where the the reads map to in the genome, and how well they map to that position (i.e. how much they vary from the reference genome). The `sarek` pipeline can either start with a pre-aligned BAM file like this, or it can alternatively start from scratch with FASTQ files and perform the alignment itself.

    In order to keep this workshop as simple and straight-forward as possible, we will only be running the `markduplicates` part of the `sarek` pipeline. This is an important quality control step that takes our aligned sequencing reads in the BAM format and marks which ones are duplicates of another read. Duplicates occur in sequencing libraries that were generated using PCR to amplify the DNA. This amplification process is incredibly useful at creating usable sequencing libraries from small amounts of DNA, but it also generates lots of duplicate reads that can interfere with our analyses. When we are trying to detect genomic variants in a sample, we want to have a measure of how confident we are in those variant calls. Duplicate reads can lead us to being over-confident in our variant calls, as they make it look like we have lots of evidence for a particular variant, when in fact they all originate from the same original DNA fragment. As such, we mark duplicate reads so that they can be skipped by the variant caller.

    Finally, the workflow we run today will generate a final MultiQC report that summarises the output of the `markduplicates` process and displays it in a handy HTML web page. If you were running other parts of the workflow, such as the variant calling stage, many of the outputs of those processes would also be captured in the same MultiQC report.

## 2.1.3 nf-core pipeline structure

At the [start of today](./01_0_intro.md#102-setup-the-project-space), we ran some setup scripts that downloaded the `nf-core/sarek` pipeline. You should see the `sarek/` folder in your current directory:

```bash
ls -1
```

```console title="Output" hl_lines="3"
config/
README.md
sarek/
scripts/
singularity/
```

Let's take a few minutes to explore the `sarek/` directory and understand the strucutre of a typical nf-core pipeline.

!!! example "Exercise: review the `nf-core/sarek` pipeline structure"

    First, let's list the directory structure of the `sarek` pipeline:

    ```bash
    cd sarek/
    ls -1
    ```

    ```console title="Output"
    assets/
    bin/
    CHANGELOG.md
    CITATIONS.md
    CODE_OF_CONDUCT.md
    conf/
    docs/
    LICENSE
    main.nf
    modules/
    modules.json
    nextflow_schema.json
    nextflow.config
    nf-test.config
    README.md
    singularity/
    subworkflows/
    tests/
    tower.yml
    workflows/
    ```

    First, note that the `singularity/` directory isn't part of the `sarek` pipeline, but is instead something we created at the start of the day in our setup script. This holds the singularity images that we will be using later on.

    Next, you'll see that there are quite a lot of files and folders. Some of these are simply housekeeping like the `LICENCE`, `CHANGELOG.md`, `CITATIONS.md`, and `CODE_OF_CONDUCT.md` files. The most important files and folders that we are interested in are described below:

    | File/folder | Description |
    | ----------- | ----------- |
    | main.nf | This is the entry point to workflow and the actual script that needs to be run to execute the pipeline |
    | workflows/ | This folder houses the main `sarek` workflow |
    | subworkflows/ | This folder houses smaller, self-contained, modular workflows that are called within the greater `sarek` pipeline |
    | modules/ | This houses Nextflow files that define individual processes. Breaking up your processes into modules like this is considered best practice. |
    | conf/ | This houses configuration files specific to various modules in the pipeline |
    | nextflow.config | This is the default configuration file and defines default parameters and profiles, and imports the module-specific configuration files defined in `conf/` |
    | bin/ | The `bin/` folder can house executable scripts that can be called directly from your Nextflow processes. In `sarek`, the only executable script in here is a python script to generate licence messages. |

    Diving a little deeper, you will see that the `modules/` and `subworkflows/` directories both contain two sub-folders: `local/` and `nf-core/`:

    ```bash
    ls modules subworkflows
    ```

    ```console title="Output"
    modules:
    local   nf-core

    subworkflows:
    local   nf-core
    ```

    Any modules and sub-workflows that originated from the [nf-core GitHub repository](https://github.com/nf-core/modules) will be placed in the `nf-core` sub-folder, and typically represent common tools and workflows that are useful in many pipelines. The `local` sub-folders, on the other hand, are for processes that were developed specificially for the current pipeline.

As you can see, nf-core pipelines can be quite complex. This is due to the attempt to heavily standardise and modularise these workflows. There are significant benefits to this, primarily in that it gives everyone a common template to work from and helps to break down the workflows into smaller, more manageable and maintainable chunks. However, it can also make it difficult to analyse and troubleshoot the code when developing them.

!!! example "Exercise: review the configuration files"

    Let's next have a look at the default configuration file: `nextflow.config`. Open it up in VSCode. You will see that it is quite long: 481 lines!

    ![sarek nextflow.config preview](./figs/02_01_nextflow_config.png)

    At line 9, you will see a `params` section with numerous parameters being defined, along with default settings for each one.

    Further down, at line 142, there is a `profiles` section. Profiles are a convenient way for setting groups of related parameters and options for running a pipeline at once, and are selected at the command line with the `nextflow run -profile <profile_name> ...` syntax. For example, lines 182-191 define the `singularity` profile, which enables singularity and disables other, mutually-exclusive backends like docker and conda.

    ![sarek profile definitions](./figs/02_01_sarek_profiles.png)

    Most of the rest of the `nextflow.config` file is dedicated to importing additional configuration files from the `conf/` folder; in this way, configurations are modularised like the pipeline itself, and help to break up the large set of parameters and options into more easily readable and maintainable chunks. We'll quickly take a look at one of these imported configuration files: on line 140, the `conf/base.config` file is imported, which defines a lot of the default resource requirements for the pipeline. In VSCode, open up `conf/base.config`:

    ![sarek base config](./figs/02_01_sarek_base_config.png)

    At the top, the file defines the default number of CPUs, amount of RAM, and wall time required by each process:

    ```groovy linenums="11"
    process {

        // TODO nf-core: Check the defaults for all processes
        cpus   = { 1      * task.attempt }
        memory = { 6.GB   * task.attempt }
        time   = { 4.h    * task.attempt }
    ```

    Note the use of the curly braces and the `* task.attempt` in each line. This is an example of a **dynamically calculated resource**. Nextflow keeps track of how many times a given task has been attempted; this can be used in case of a failure to increase the required resources on the next try. In this case, by default, jobs will be given 1 CPU, 6GB of memory, and 4h of walltime on the first try; on the second try, they will get 2 CPUs, 12GB of memory, and 8h of time. Line 20 sets `maxRetries = 1`, so this will only happen the one time, and failures won't cause an infinite loop of retries.

    Further down in the file, we see examples of overriding the default resources based on process names and labels. Nextflow processes can be specifically configured based on their name by using the `process.withName` scope. For example:

    ```groovy linenums="66"
    withName: 'FASTP'{
        cpus   = { 12   * task.attempt }
        memory = { 4.GB * task.attempt }
    }
    ```

    This chunk of code specifies that the `FASTQC` process will get 12 CPUs and 4GB of memory (multiplied by the attempt number). This overrides the defaults defined at the top of the file. Since `time` isn't specified here, the default run time of 4 hours is kept the same.

    The `process.withLabel` scope similarly lets us override parameters for any process that has a given label assigned to it. Labels are a Nextflow directive that can be applied to any process, and multiple labels can be applied to a single process. This is very useful for tasks that have very similar resource requirements: instead of individually specifying their resources with the `process.withName` scope, we can simply assign each process the same label, and then define the resources for that label. For example, any process that was labelled as `process_medium` will get 6 CPUs, 36GB memory, and 8h walltime:

    ```groovy linenums="42"
    withLabel:process_medium {
        cpus   = { 6     * task.attempt }
        memory = { 36.GB * task.attempt }
        time   = { 8.h   * task.attempt }
    }
    ```

Now we have seen what a typical nf-core pipeline looks like, in the next section we will write a script to try and run it!

!!! note "Move back to the parent directory"

    Before we move on, be sure to move back to the partent directory (`part1/`):

    ```bash
    cd ../
    ```
