# 1.5 Intro to nf-core

Now that we understand the fundamentals of running jobs on an HPC, we want to combine this with the reproducibility of Nextflow to allow us to consistently process data through a workflow at scale. To do this, we will introduce a community-driven Nextflow project called [nf-core](https://nf-co.re/).

## 1.5.1 What is nf-core?

![nf-core logo](/docs/assets/nf-core.png)

nf-core is a community-driven project to develop and curate open-source analysis workflows built with Nextflow.

The project has a standardised set of best practices, guidelines, and templates for building modular, scalable, and portable bioinformatics workflows. Every workflow is well-documented, tested, and designed to work across multiple platforms, including cloud and HPC.

The key Features of nf-core workflows are:

- Documentation: nf-core workflows have extensive documentation covering installation, usage, and description of output files to ensure that you won’t be left in the dark.
- CI Testing: Every time a change is made to the workflow code, nf-core workflows use continuous-integration testing to ensure that nothing has broken.
- Stable Releases: nf-core workflows use GitHub releases to tag stable versions of the code and software, making workflow runs totally reproducible.
- Packaged software: Pipeline dependencies are automatically downloaded and handled using Docker, Singularity, Conda, or other software management tools. There is no need for any software installations.
- Portablility and reproducibility: nf-core workflows follow best practices to ensure maximum portability and reproducibility. The large community makes the workflows exceptionally well-tested and easy to execute.
- Cloud-ready: nf-core workflows are tested on AWS after every major release. You can even browse results live on the website and use outputs for your own benchmarking.

nf-core is published in Nature Biotechnology: [Nat Biotechnol 38, 276–278 (2020). Nature Biotechnology](https://www.nature.com/articles/s41587-020-0439-x)

## 1.5.2 Available nf-core workflows

There are currently 139 workflows (October 2025) available as part of nf-core. These workflows are at various stages of development with 84 released, 43 under development, and 12 archived.

The [nf-core website](https://nf-co.re/pipelines/) has a full list of workflows, as well as their documentation, which can be explored.

Each workflow has a dedicated page that includes expansive documentation that is split into 6 sections:

- Introduction: An introduction and overview of the workflow
- Usage: Documentation and descriptions of how to execute the workflow
- Parameters: Documentation for all workflow parameters
- Output: Descriptions and examples of the expected output files
- Results: Example output files generated from the full test dataset on AWS
- Releases: Workflow version history

Unless you are actively developing workflow code, you don’t need to clone the workflow code from GitHub and can use Nextflow’s built-in functionality to pull and a workflow:

```bash
nextflow pull nf-core/<pipeline>
```

Nextflow run will also automatically pull the workflow if it was not already available locally:

```bash
nextflow run nf-core/<pipeline>
```

By default, Nextflow will pull the default git branch of the pipeline unless a specific version is specified with the `-revision` or `-r` flag.

## 1.5.3 nf-core tools

The nf-core project also has a command-line utility suite called **nf-core tools** that has been developed to aid using, developing, and testing nf-core workflows. We will make use of this tool today.

!!! note "Installing nf-core tools"

    During the workshop, the nf-core tools suite will be pre-installed on your system.

    To install nf-core tools yourself on another computer, you can do so through both the [Python Package Index (PyPI)](https://pypi.org/project/nf-core/) or [Bioconda](https://anaconda.org/bioconda/nf-core):

    ```bash title="Install nf-core tools via PyPI"
    pip install nf-core
    ```

    ```bash title="Instal nf-core tools via Bioconda"
    conda install -c bioconda nf-core
    ```
