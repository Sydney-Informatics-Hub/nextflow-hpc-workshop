#!/usr/bin/bash

for SAMPLE_ID in gut liver lung; do

	qsub \
		-N fastqc \
		-l wd \
		-l storage=gdata/za08+scratch/za08 \
		-l walltime=1:00:00 \
		-l mem=1GB \
		-l ncpus=1 \
		-v SAMPLE_ID="${SAMPLE_ID}" \
		fastqc.gadi.sh

done