#!/bin/bash

SAMPLE_ID="tiny"
READS_1="data/${SAMPLE_ID}.R1.fq"
READS_2="data/${SAMPLE_ID}.R2.fq"

mkdir -p "results/fastqc_${SAMPLE_ID}_logs"
fastqc \
    --outdir "results/fastqc_${SAMPLE_ID}_logs" \
    --format fastq ${READS_1} ${READS_2}