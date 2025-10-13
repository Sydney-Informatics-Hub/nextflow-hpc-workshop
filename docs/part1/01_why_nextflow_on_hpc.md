# 1.1 Why Nextflow on HPC?

## 1.1.0 The scenario

You are a bioinformatician working in a busy genomics lab. The lab group has started an exciting new genomics project and will be analysing lots of RNA sequencing samples, with each sample having sizeable sequencing depth. This means that the dataset that will be generated will be quite large - hundreds of gigabytes of data - and processing everything on a laptop or local desktop machine will not be feasible. Instead, you need to be able to process this huge dataset in a timely manner, while also ensuring reproducibility and consistency of the analysis workflow between samples. For these reasons, you decide to use Nextflow on a HPC system, enabling parallel processing of each sample.

## 1.1.1 The current pipeline

For most of this session, we will be working with the [`nf-core/rnaseq` Nextflow pipeline](https://nf-co.re/rnaseq). This is a community-developed workflow that handles everything from initial quality control of the RNA sequencing data through to mapping to the genome and quantification of transcripts. As we will see shortly, this pipeline works fine on a local machine for small datasets, but with larger datasets it will take too long to run locally, and may even exceed memory and storage limits.

## 1.1.2 Running `nf-core/rnaseq` locally on a tiny dataset



## 1.1.3 Attempting to run a larger dataset

