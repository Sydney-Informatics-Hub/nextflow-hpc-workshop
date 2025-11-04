# Assigning process resources

!!! info "Learning objectives"

    - Identify how to size resource requests appropriately for each process
    in a workflow
    - Apply resource-aware design principles to improve job efficiency and
    reproducibility
    - Optimise processes for time, noting the trade-offs with cost (SU usage)

## Configuring process usage

https://sydney-informatics-hub.github.io/template-nf-guide/notebooks/modules.html

modules

withName

withLabel

!!! example "Exercises"

    TODO View trace, configure resources for each process in nextflow.config


## Configuring java heap sizes

!!! example "Exercises"

    TODO Update GENOTYPE and JOINT_GENOTYPE processes with -Xmx${tasks.memory}

Other things to consider - when writing custom R or Python scripts, writing
them efficiently. Utilising things like vectorisation, libraries such as numpy
etc.
