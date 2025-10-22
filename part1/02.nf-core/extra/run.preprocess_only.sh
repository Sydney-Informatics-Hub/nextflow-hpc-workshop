#!/bin/bash

module load nextflow
module load singularity

nextflow run -r 3.5.0 nf-core/sarek \
    --input samplesheet.fq.csv \
    --step mapping \
    --skip_tools markduplicates,baserecalibrator,mosdepth,samtools \
    --outdir results \
    --no_intervals true \
    --bwa ref/BWAIndex \
    --fasta ref/Homo_sapiens_assembly38.20-22.fasta \
    --fasta_fai ref/Homo_sapiens_assembly38.20-22.fasta.fai \
    --igenomes_ignore true \
    -c gadi.config,custom.config \
    -resume
