#!/bin/bash

module load nextflow

nextflow run nf-core/sarek -r 3.6.0 \
    --input samplesheet.fq.csv \
    --step mapping \
    --skip_tools markduplicates,baserecalibrator,mosdepth,samtools \
    --outdir results \
    --no_intervals true \
    --bwa ../../../nextflow-on-hpc-materials/data/ref \
    --fasta ../../../nextflow-on-hpc-materials/data/ref/Hg38.subsetchr20-22.fasta \
    --fasta_fai ../../../nextflow-on-hpc-materials/data/ref/Hg38.subsetchr20-22.fasta.fai \
    --igenomes_ignore true \
    -c gadi.config,custom.config \
    -resume
