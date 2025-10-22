#!/bin/bash

module load nextflow/24.04.5
module load singularity

nextflow run -r 3.5.0 nf-core/sarek \
    --input samplesheet.bam.csv \
    --step markduplicates \
    --skip_tools baserecalibrator,mosdepth,samtools \
    --outdir results \
    --no_intervals true \
    --fasta ref/Homo_sapiens_assembly38.20-22.fasta \
    --fasta_fai ref/Homo_sapiens_assembly38.20-22.fasta.fai \
    --igenomes_ignore true \
    -c gadi.config,custom.config \
    -resume
