#!/bin/bash
#PBS -q normal
#PBS -P er01
#PBS -l ncpus=1
#PBS -l mem=1GB
#PBS -l walltime=00:02:00
#PBS -l storage=scratch/er01
#PBS -l wd

set -euo pipefail

mkdir -p results/fastq_split_${SAMPLE_ID}
FQ1=data/${SAMPLE_ID}.R1.fq
FQ2=data/${SAMPLE_ID}.R2.fq

N_LINES=$(awk -v chunks="${CHUNKS}" 'NR % 4 == 1 { nreads += 1 } END { if (nreads % chunks == 0) { chunkreads = nreads / chunks } else { chunkreads = (nreads + chunks - (nreads % chunks)) / chunks }; printf("%0.f\n", chunkreads * 4) }' ${FQ1})

split -l ${N_LINES} -d --additional-suffix=.R1.fq ${FQ1} results/fastq_split_${SAMPLE_ID}/${SAMPLE_ID}.split_
split -l ${N_LINES} -d --additional-suffix=.R2.fq ${FQ2} results/fastq_split_${SAMPLE_ID}/${SAMPLE_ID}.split_

echo DONE
