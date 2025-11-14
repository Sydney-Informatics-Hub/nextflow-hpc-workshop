
=== "Gadi (PBSpro)"

    ```bash
    cat conf/pbspro.config
    ```
    ```groovy title="conf/pbspro.config"
    params.pbspro_account = ""

    process {
      executor = 'pbspro'
      queue = 'normalbw'
      clusterOptions = "-P ${params.pbspro_account}"
      module = 'singularity'
    }

    singularity {
      enabled = true
      autoMounts = true
      cacheDir = "${projectDir}/singularity"
    }
    ```

=== "Setonix (Slurm)"

    ```bash
    cat conf/slurm.config
    ```
    ```groovy title="conf/slurm.config"
    params.slurm_account = ""

    process {
      executor = 'slurm'
      queue = 'work'
      clusterOptions = "--account=${params.slurm_account}"
      module = 'singularity/4.1.0-slurm'
    }

    singularity {
      enabled = true
      autoMounts = true
      cacheDir = "${projectDir}/singularity"
    }
    ```

This setup makes it easy for us to run the same workflow in different environments. By the end of Part 2, we’ll have extended these configs to better reflect the characteristics of each system, improving efficiency without touching the workflow logic itself.