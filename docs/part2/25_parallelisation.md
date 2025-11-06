# Parallelisation

!!! info "Learning objectives"

    - Consider the limitations of parallelisation and scenarios where it should
    not be applied
    - Differentiate between multithreading and scatter-gather approaches to
    parallelisation
    - Implement parallelisation approaches (multithreading, multiprocessing)
    while preserving biological correctness


## Multithreading 

!!! example "Exercise"

    TODO Recall whether FASTQC() -t ${task.cpus} will make is more efficient or not

!!! example "Exercise"

    TODO In script processes GENOTYPE() and JOINT_GENOTYPE(), add --native-pair-hmm-threads ${task.cpus}

## Multi-processing with scatter-gather

!!! example "Exercise"

    TODO Add modules for SPLIT_FASTQ(), ALIGN_CHUNK(), MERGE_BAM()
    
!!! example "Exercise"

    TODO Update the channels in main.nf

## Bonus: Dynamic resourcing

!!! Exercise

    TODO error strategies with e.g. maxretries * task.memory
