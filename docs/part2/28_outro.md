# Summary and takeaways

Optimising Nextflow pipelines on HPC involves more than just making it "run". It’s about making it run efficiently and being a responsible user on a shared system. Below are a few key takeaways when you apply to your own workflows.

## Get familiar with your HPC environment

Each HPC has unique hardware, limits, charging models, and ways to work. Before you run anything:

- Familiarse yourself with the user documentation
- Know where to run jobs (not on the login node!)
- Check for network access, especially if pulling containers or pipelines on-the-go
- Attend your HPCs next training event
- **Reach out to your HPC support team** - they are your best source of optimisation advice

## Start small, then scale

Before launching full-scale analyses:

- Run the workflow on a small, representative sample
- Ensure the pipeline works end-to-end with resources configured
- Use this opportunity to benchmark runtime and memory, cpu usage

## Use Nextflow's built-in profiling features

Use the trace and timeline outputs to guide configuration decisions. Review these before and after tuning.

- Identify long-running or ineficient processes
- Spot over-provisioned tasks (e.g. jobs using a small % of allocated memory)
- Adjust as needed, **remembering to configure to fit your HPC queues/partitions**

## Choose the right optimisation strategy

Does the tool support multithreading?
→ Increase CPUs, add ${task.cpus} to the `process`

Can data be split into chunks without breaking biology?
→ Use scatter-gather patterns

Would more resources help speed up processing time?
→ Match requests to node sizes (e.g. 2GB/core)

Is your input data heterogeneous?
→ Consider dynamic resource configuration 

## Iterate and adapt

Efficiency is context-specific. What works for one sample or one HPC may not work for another.

- Run benchmarks iteratively, especially if scaling to 10s - 1000s of samples
- Keep your pipeline reproducible and flexible across systems - **separate _what_ the workflow runs to _how_ and _where_ it should run**
- Document and comments any tuning decisions to support reproducibility

## Resources

### Developed by us

* [SIH Nextflow template](https://github.com/Sydney-Informatics-Hub/template-nf)
* [SIH Nextflow template guide](https://sydney-informatics-hub.github.io/template-nf-guide/)
* [SIH Nextflow for the life sciences 2025](https://sydney-informatics-hub.github.io/hello-nextflow-2025/)
* [SIH Customising nf-core workshop](https://sydney-informatics-hub.github.io/customising-nfcore-workshop/)
* [SIH Introduction to NCI Gadi - Optimisation](https://sydney-informatics-hub.github.io/training.gadi.intro/07-Optimisation/index.html)
* [Australian BioCommons Seqera Platform Service](https://www.biocommons.org.au/seqera-service)
* [NCI Gadi nf-core instutitonal config](https://nf-co.re/configs/nci_gadi/)
* [Pawsey Setonix nf-core instutitional config](https://nf-co.re/configs/pawsey_setonix/)

### Developed by others

* [Nextflow training](https://training.nextflow.io/)
* [Nextflow patterns](https://nextflow-io.github.io/patterns/index.html)
* [Nextflow blog](https://www.nextflow.io/blog.html)
* [Nextflow coding best practice recommendations](https://carpentries-incubator.github.io/Pipeline_Training_with_Nextflow/07-Nextflow_Best_Practice/index.html)
* [Seqera community forums](https://community.seqera.io/)