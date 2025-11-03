# 1.4 HPC architecture for workflows

!!! info "Learning objectives"

    - Describe the roles of login nodes, compute nodes, and shared storage in HPC systems.
    - Distinguish between login nodes and compute nodes and know where workflows execute.
    - Explain how shared filesystems are accessed and managed safely in workflows.
    - Define job resource requirements (CPU, memory, walltime) and their impact on scheduling.

While HPCs can look intimidating, their architecture follows a simple structure that supports large-scale computation through shared resources. From a workflow perspective, this architecture means there are a few important realities to accept: work is not run interactively, resources must be requested rather than assumed and everything is governed by shared access.

![](figs/00_hpc_architecture.png)

## Login nodes

When a user connects to an HPC, they first land on a login node. This is a shared access point used to prepare work, not perform computations. From here, users submit jobs to the scheduler, monitor their progress and organise their project directories. The login node exists only to coordinate access to the system, and because it is shared by many people at once, it must not be overloaded with computational tasks.

## Compute nodes

The real work happens on the compute nodes. These are powerful machines with many CPU cores, large amounts of memory and fast access to storage. Workflows do not run directly on them; instead, the scheduler assigns workflow tasks to available compute nodes based on the resources requested. This separation between the login node and compute nodes allows users to interact with the system while computation is queued and executed elsewhere.

HPC systems divide compute resources into queues (sometimes called partitions in Slurm). Each queue represents a group of compute nodes with specific hardware characteristics, limits, and intended usage.

Queues help balance the system between different types of work, for example, short interactive tasks, long-running simulations, or large parallel jobs. Each queue enforces policies on maximum runtime (walltime), maximum cores or memory per job, job priority and eligible users or projects.

Because queues are mapped to distinct sets of compute nodes, requesting the right queue for your job is important for both performance and fairness.

Queue names and limits differ across infrastructure providers but some common examples include:

| Queue type        | Description                           | Typical limits                                         |
| ----------------- | ------------------------------------- | ------------------------------------------------------ |
| Normal/work       | Default queue for most batch jobs     | Moderate walltime (e.g. 24–48 h), general-purpose CPUs |
| Express/short     | Prioritised for rapid turnaround      | Small jobs, short walltime                             |
| Large/high memory | For jobs requiring many CPUs or nodes | Long walltime, higher resource requests                |
| GPU               | Access to GPU-enabled nodes           | GPU-specific workloads only                            |

For details about queue limits and scheduling policies on systems used today, see:

- [Setonix (Slurm) – Running Jobs on Setonix (Pawsey)](https://pawsey.atlassian.net/wiki/spaces/US/pages/51929058/Running+Jobs+on+Setonix)
- [Gadi (PBS Pro) – Queue Limits (NCI)](https://opus.nci.org.au/spaces/Help/pages/236881198/Queue+Limits)

!!! example "Exercise 1.4.1: Wheres my Jobs running ? - Login vs compute nodes"

    First, run the `hostname` command where you're working from, this is the **login node**, and its name varies depending on the infrastructure provider.

    === "Gadi (PBS)"

        ```bash
        hostname
        ```
        your will see the name of the node you ran your command on
        ```console
        gadi-login-05.gadi.nci.org.au
        ```
        Gadi uses descriptive hostnames that include the infrastructure name (gadi-login-XX).


    === "Setonix (Slurm)"

        ```bash
        hostname
        ```
        your will see the name of the node you ran your command on
        ```console
        setonix-01
        ```
        Setonix opts for simpler hostnames, often just setonix-XX.

    Next, submit a simple job to the scheduler to run the same command on a **compute node**.

    === "Gadi (PBS)"

        ```bash
        echo hostname | qsub
        ```
        PBS will return a job ID like:
        ```console
        123456.gadi-pbs
        ```

    === "Setonix (Slurm)"

        ```bash
        sbatch --wrap="hostname"
        ```
        ```console
        Submitted batch job 123456
        ```
    Once the job finishes, check the output file to confirm it ran on a compute node. Node names will differ from login nodes and vary by infrastructure.

    === "Gadi (PBS)"

        ```bash
        cat STDIN.o123456
        ```
        ```console
        gadi-cpu-clx-0534.gadi.nci.org.au

        ======================================================================================
                    Resource Usage on 2025-11-03 11:54:04:
        Job Id:             123456.gadi-pbs
        Project:            [project_id]
        Exit Status:        0
        Service Units:      0.00
        NCPUs Requested:    1                      NCPUs Used: 1
                                                CPU Time Used: 00:00:00
        Memory Requested:   500.0MB               Memory Used: 7.09MB
        Walltime requested: 00:01:00            Walltime Used: 00:00:01
        JobFS requested:    100.0MB                JobFS used: 0B
        ======================================================================================
        ```
        Note: On PBS systems like Gadi, a job file is automatically generated with detailed resource usage information and appended to the end of the job output file. This is not common on Slurm systems.


    === "Setonix (Slurm)"

        ```bash
        cat slurm-123456.out
        ```
        ```console
        nid002024
        ```
        Slurm does not generate a detailed job file by default—only a basic output file unless configured otherwise.


    This confirms that compute jobs are executed on separate **compute nodes**, not the **login nodes**.

## Shared storage

All nodes are connected to a shared parallel filesystem. This is a large, high-speed storage system where input data, reference files and workflow outputs are kept. Because it is shared across all users, it enables collaborative research and scalable workflows. However, it also introduces constraints around file organisation and performance, which is why workflows must be careful about how they read and write data here.

!!! example "Exercise 1.4.2: Hello from the other side - shared filesystems"

    Both login and compute nodes share the same file systems (e.g., `/scratch`).
    You can demonstrate this by writing to a file from the login node and then appending to it from a compute node.

    First, on the **login node**, create a simple file:

    ```bash
    echo "hello from " $(hostname) > ./shared_storage_system.txt
    cat shared_storage_system.txt
    ```
    ```console
    hello from gadi-login-05.gadi.nci.org.au
    ```
    This created a file in the current directory (./shared_storage_system.txt) and wrote the hostname of the login node into it.

    Next, submit a job that appends a new line from the compute node:

    === "Gadi (PBS)"

        ```bash
        echo "echo 'hello from' \$(hostname) >> $(pwd)/shared_storage_system.txt" | qsub
        ```

    === "Setonix (Slurm)"

        ```bash
        sbatch --wrap="echo 'hello from' \$(hostname) >> ./shared_storage_system.txt"
        ```

    This submitted a job to a compute node, which appended its own hostname to the same file.

    After the job completes, check the contents of the file again:

    === "Gadi (PBS)"

        ```bash
        cat shared_storage_system.txt
        ```
        ```console
        hello from gadi-login-05.gadi.nci.org.au
        hello from gadi-cpu-clx-0134.gadi.nci.org.au
        ```

    === "Setonix (Slurm)"

        ```bash
        cat shared_storage_system.txt
        ```
        ```console
        hello from setonix-01
        hello from nid002041
        ```

    The first entry in the file reflects the login node’s hostname, while the second entry corresponds to the compute node on which the job ran. Repeating the job may produce a different compute node name. This demonstrates that both login and compute nodes can access and modify the same file, confirming that they share a common storage system.

!!! warning "Overwriting Files in Shared Filesystems"

    Because the filesystem is shared, multiple jobs or users writing to the same file, especially if it has a common name, can accidentally overwrite each other’s data or cause race conditions.

    In the example above, we safely added a line to a file using a method that appends without deleting what's already there.

    However, in real workflows, most tools and scripts write to files in a way that replaces the entire file. On shared systems like /scratch, where many users and jobs access the same space, this can lead to conflicts or data loss.

    It’s good practice to:

    - Use unique filenames (e.g., include job ID or timestamp).
    - Avoid simultaneous writes unless explicitly managed.
    - Add file checks and error handling

## Job scheduler

At the centre of everything is the job scheduler. Rather than allowing users to run programs directly, HPCs rely on a scheduling system (e.g. Slurm or PBS Pro) to manage fair access to shared compute resources. When a job is submitted, it enters a queue where the scheduler decides when and where it will run. Jobs are matched to compute nodes based on requested resources like CPU, memory and runtime. Understanding how the scheduler behaves is essential for designing workflows that run efficiently.

Schedulers like PBS Pro and Slurm use queues to group jobs that share similar resource and policy constraints. When you submit a job, it’s placed in the appropriate queue, and the scheduler continuously evaluates all queued jobs to decide which can start next.

System load — high demand increases queue wait time
![](figs/00_HPC_scheduler_tetris.png)

!!! note "Understanding Job Scheduling with Tetris"

    Think of the scheduler like a giant game of **Tetris**, where every job you submit has its own unique “shape.”
    When you submit a batch job, you describe the resources it needs — this defines the shape of your job piece.
    The scheduler’s task is to fit all these different pieces together as efficiently as possible across the available compute nodes.

    Each job’s shape is determined by three key factors:

    - **CPU** – how many processor cores it needs
    - **Memory** – how much RAM it requires
    - **Walltime** – how long it is allowed to run

    Once submitted, your job enters a **queue**, much like a waiting line for compute resources. But unlike a simple first-come-first-served queue, the scheduler constantly reshuffles and fits jobs together — like sliding Tetris blocks — to maximise system usage.

    The order in which jobs run depends on several factors:

    - **Job priority** — determined by project, queue, and fair-share usage
    - **Requested resources** — smaller jobs can often “slot in” sooner
    - **Queue limits** — different queues prioritise short, long, or interactive jobs

    Getting the shape right matters. Overestimating your needs makes your job harder to fit, while underestimating can cause it to fail. Common outcomes include:

    - **Longer queue times** – large, awkwardly shaped jobs wait for space
    - **Wasted capacity** – unused cores or memory that could have run other jobs

    Just like in Tetris, the scheduler aims to **fill every gap** and keep the system running smoothly.
    Small, well-shaped jobs often fall neatly into open spaces, while larger ones wait for the perfect fit.
