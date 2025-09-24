#!/bin/bash

for SAMPLE_ID in gut liver lung; do

	sbatch \
		--account=$PAWSEY_PROJECT \
		--job-name=fastqc \
		--output='fastqc_%j.log' \
		--nodes=1 \
		--ntasks=1 \
		--cpus-per-task=2 \
		--mem=8GB \
		--partition=work \
		--time=00:10:00 \
		--mail-user='your@email.com' \
		--mail-type='END,FAIL'
		fastqc.setonix.sh ${SAMPLE_ID}

done