#!/bin/bash

module load nextflow
module load singularity

export NXF_SINGULARITY_CACHEDIR="$PWD/singularity"

set -euo pipefail

nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results \
    --fasta $PWD/data/ref/fasta/chr22.fa \
    --gtf $PWD/data/ref/genes/chr22.gtf \
    --skip_markduplicates \
    --save_trimmed true \
    --save_unaligned true \
    -c custom.config \
    -profile singularity

echo DONE