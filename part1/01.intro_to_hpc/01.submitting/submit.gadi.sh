#!/usr/bin/bash

qsub \
	-N fastqc \
	-l wd \
	-l storage=gdata/za08+scratch/za08 \
	-l walltime=1:00:00 \
	-l mem=1GB \
	-l ncpus=1 \
	fastqc.sh