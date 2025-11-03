# 2.4 Layering Nextflow configurations

!!! info "Learning objectives"

    - Learn how Nextflow configurations can be layered on top of one another
    - Understand the configuration priorities
    - Learn how to fine-tune processes using the `withName` directive
    - Understand when you would want to fine-tune process requirements
    
## 2.4.1 Configuration priorities

TODO: Add explanation of configuration layers and the priorities of each layer.

- Config files > process definitions > process defaults
- Config files:
    - `-c` > `nextflow.config`
        - `-c`: last > first
    - `withName` > `withLabel` > defaults
        - `withLabel`: last > first

## 2.4.2 Optimising `nf-core/sarek` for our data

!!! example "Exercise: Fine-tune `nf-core/sarek`"

    TODO:

    - Add new blank file: `custom.config`
    - Add `process {}` scope
    - Add `withName` scopes for each process we run in the workflow, each with `cpus`, `memory`, and `time` directives specified
    - Run the workflow again, hopefully it runs quicker this time!