#!/usr/bin/bash

SCRIPTDIR=$(cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd)

# Run the Index script
qsub -P za08 $SCRIPTDIR/00_index.sh

# Run the fastqc script
# This will only start running after the previous one is completed
current_jobs=`qselect -N index | cut -d'.' -f1 | xargs -I % sh -c '{ echo -n  % ; echo -n ":" ; }' | sed 's/:$//' `
#echo $current_jobs
qsub -P za08 -W depend=afterany:$current_jobs -v SAMPLE_ID=gut $SCRIPTDIR/01_fastqc.sh

# Run the quant script
# This will only start running after the previous one is completed
current_jobs=`qselect -N fastqc | cut -d'.' -f1 | xargs -I % sh -c '{ echo -n  % ; echo -n ":" ; }' | sed 's/:$//' `
#echo $current_jobs
qsub -P za08 -W depend=afterany:$current_jobs -v SAMPLE_ID=gut $SCRIPTDIR/02_quant.sh

# Run the multiqc script
# This will only start running after the previous one is completed
current_jobs=`qselect -N quant | cut -d'.' -f1 | xargs -I % sh -c '{ echo -n  % ; echo -n ":" ; }' | sed 's/:$//' `
#echo $current_jobs
qsub -P za08 -W depend=afterany:$current_jobs -v SAMPLE_ID=gut $SCRIPTDIR/03_multiqc.sh
