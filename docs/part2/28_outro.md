# Summary and takeaways

TODO “Best practices for running Nextflow on HPC” Points:

- read the user guide for the HPC you are running on  
- test on one small representative sample first

TODO Usage tips when doing things yourself:

- where to run your main job (not login node)  
- considering network access for pulling singularity containers
- etc.

Before optimising a task, consider:

- Does the tool support threads? → Use multithreading.
- Can the data be split meaningfully without compromising the biology? → Use scatter-gather.
- Will more CPUs or memory reduce queue time? → Check HPC architecture.
- Do I need to process widely variable data sizes? → Dynamic CPU/memory scaling later?