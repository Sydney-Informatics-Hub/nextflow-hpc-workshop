# 1.2 Submitting jobs to an HPC

As mentioned in the previous section, working on an HPC is a bit different to working locally on a laptop or desktop. Instead of running jobs interactively via the command line, you are instead required to write a script that you then submit to the scheduler. The script is run on one of the compute nodes, and then after it completes, you can look through the outputs and the log files to ensure that everything has run successfully.

In this section, we will guide you through writing a simple script and submitting it to your HPC scheduler. In the process, we will explore how resources can be requested both on the command line and within the script itself, and how we can capture the standard ouptut and error streams to monitor the job and check that it has successfully run.

## 1.2.1 A simple FASTQC script

For all of the following examples, we will be using a simple script that runs the `FASTQC` software. This is a very common bioinformatics program that reads through one or more FASTQ files and constructs a quality report, looking at various metrics such as read length, quality scores, GC content, and more.

The basic FASTQC script that we will build on looks like this:

```bash title="fastqc.sh"
#!/usr/bin/bash

SAMPLE_ID="gut"
READS_1="data/ggal/${SAMPLE_ID}_1.fq"
READS_2="data/ggal/${SAMPLE_ID}_2.fq"

mkdir -p "results/fastqc_${SAMPLE_ID}_logs"
fastqc \
	--outdir "results/fastqc_${SAMPLE_ID}_logs" \
	--format fastq ${READS_1} ${READS_2}
```

!!! note

    Create a new blank file, copy and paste the script above, and save it as `fastqc.sh`

Let's break this down into it's components.

```bash
#!/usr/bin/bash
```

The first line starts with what's known as a "shebang": `#!`. This is a special comment that is used by the operating system to determine how to execute the script. In this case, we tell the OS to run the script using `bash` (which is located at `/usr/bin/bash` in the file system).

```bash
SAMPLE_ID="gut"
READS_1="data/ggal/${SAMPLE_ID}_1.fq"
READS_2="data/ggal/${SAMPLE_ID}_2.fq"
```

The next three lines define some bash variables: `SAMPLE_ID`, `READS_1` and `READS_2`. These are used to make the script a little more flexible and to prevent duplicating long paths that might introduce errors.

```bash
mkdir -p "results/fastqc_${SAMPLE_ID}_logs"
```

In the next line, we create an output directory to store our FASTQC reports in. The directory has a variable name based on the value stored in the `SAMPLE_ID` variable.

```bash
fastqc \
	--outdir "results/fastqc_${SAMPLE_ID}_logs" \
	--format fastq ${READS_1} ${READS_2}
```

Finally, we run the `fastqc` command, providing the output directory as well as the two input FASTQ files (stored in the variables `READS_1` and `READS_2`).

As it currently exists, this script isn't quite ready for use on an HPC. First, we need to make sure that the `fastqc` command is installed and available to use. If we were running this command on a local machine, we could simply download FASTQC from its website and install it, then run the script above. But it often isn't so straight forward on an HPC, as regular users don't have the ability to install software globally (for good reason!), and so we need other ways of accessing the software we need. There are two main approaches to this, which we will discuss more in depth in [a later section](./04_software_on_hpc.md), but for now, we will use what is called the `module` system to load a pre-installed version of FASTQC. We will add a single line after the shebang line and before defining our bash variables. This will differ slightly depending on the HPC system you are using, but in both cases the new line looks like `module load fastqc/<version>`. Select your HPC system in the code block below to see the specific line you will need to add to your script.

!!! note

    Throughout this workshop, we will use content tabs to deliver code that is specific to one system or another. Simply choose the system you are using ("Gadi" or "Setonix") and follow the code example in that box.

    === "Gadi"

        ```
        This is code specific to NCI/Gadi
        ```

    === "Setonix"

        ```
        This is code specific to Pawsey/Setonix
        ```

=== "Gadi"

    ```bash title="fastqc.sh" hl_lines="3"
    #!/usr/bin/bash

    module load fastqc/0.12.1

    SAMPLE_ID="gut"
    READS_1="data/ggal/${SAMPLE_ID}_1.fq"
    READS_2="data/ggal/${SAMPLE_ID}_2.fq"

    mkdir -p "results/fastqc_${SAMPLE_ID}_logs"
    fastqc \
        --outdir "results/fastqc_${SAMPLE_ID}_logs" \
        --format fastq ${READS_1} ${READS_2}
    ```

=== "Setonix"

    ```bash title="fastqc.sh" hl_lines="3"
    #!/usr/bin/bash

    module load fastqc/0.11.9--hdfd78af_1

    SAMPLE_ID="gut"
    READS_1="data/ggal/${SAMPLE_ID}_1.fq"
    READS_2="data/ggal/${SAMPLE_ID}_2.fq"

    mkdir -p "results/fastqc_${SAMPLE_ID}_logs"
    fastqc \
        --outdir "results/fastqc_${SAMPLE_ID}_logs" \
        --format fastq ${READS_1} ${READS_2}
    ```

## 1.2.2 Submitting the job

Now comes the exciting part! We will submit our new FASTQC job to the HPC scheduler and watch it run!

Gadi and Setonix use two different scheduling systems with different commands and syntax, but the general concept is the same for both. On Gadi, we have the PBSPro scheduler, and we submit jobs using the `qsub` command; on Setonix, we have the SLURM scheduler and submit jobs with `sbatch`. Both commands take several parameters, most of which have equivalents in both schedulers. The most important parameters you need to know about are for specifying:

- Your project - important for determining access permissions and billing purposes.
- The scheduler queue on which to run - important for determining resource allocations.
- The number of CPUs required.
- The amount of memory required.
- The amount of time (also called "walltime") required by the job.

=== "Gadi"

    ```bash
    qsub \
        -P $PROJECT \
        -N fastqc \
        -q normal \
        -l ncpus=1 \
        -l mem=1GB \
        -l walltime=00:10:00 \
        -l storage=scratch/$PROJECT \
        -l wd \
        fastqc.sh
    ```

=== "Setonix"

    ```bash
    sbatch \
        --account=$PAWSEY_PROJECT \
        --job-name=fastqc \
        --partition=work \
        --nodes=1 \
        --ntasks=1 \
        --cpus-per-task=1 \
        --mem=1GB \
        --time=00:10:00 \
        fastqc.sh
    ```

Once again, let's break these commands down into their components:

=== "Gadi"

    ```bash
    -P $PROJECT
    ```

    This specifies the project that should be used for running the job. This determines what project will get billed for the resources used, and will also determine how many resources in total can be requested.

    ```bash
    -N fastqc
    ```

    This gives the job a name that can be used to identify it among potentially multiple simultaneously running jobs.

    ```bash
    -q normal
    ```

    This specifies the **queue** in which the job should run. Different queues have different resource limits and as such are intended for different types of jobs. Some queues are for resource-intensive jobs, while others are for smaller, less-intensive jobs. The `normal` queue is a good default queue to use for standard jobs.

    ```bash
    -l ncpus=1
    ```

    This specifies how many CPUs the job requires; in our case, just 1 CPU is needed for FASTQC.

    ```bash
    -l mem=1GB
    ```

    This specifies the amount of memory or RAM we require. As this is a small job, we just require 1GB of memory.

    ```bash
    -l walltime=00:10:00
    ```

    This specifies the maximum amount of time the job can run for. If the job takes any longer than what is specified here, it will be killed, so we want to make sure that we give it enough time to finish.

    ```bash
    -l storage=scratch/$PROJECT
    ```

    Here we tell Gadi that we want to mount the scratch storage space for our project. By default, the network storage systems are not mounted and are unavailable to our jobs. We need to specify which file systems we want to mount with the `storage=<filesystem>/<project>` syntax.

    ```bash
    -l wd
    ```

    This tells Gadi to run our job in the current working directory that we are in. All commands that are executed within the script will be relative to our current directory.

    ```bash
    fastqc.sh
    ```

    Finally, we provide the path to the script to submit.

=== "Setonix"

    ```bash
    --account=$PAWSEY_PROJECT
    ```

    This specifies the project that should be used for running the job. This determines what project will get billed for the resources used, and will also determine how many resources in total can be requested.

    ```bash
    --job-name=fastqc
    ```

    This gives the job a name that can be used to identify it among potentially multiple simultaneously running jobs.

    ```bash
    --partition=work
    ```

    This specifies the **queue** in which the job should run (called "partitions" on Setonix). Different partitions have different resource limits and as such are intended for different types of jobs. Some partitions are for resource-intensive jobs, while others are for smaller, less-intensive jobs. The `work` partition is a good default choice for standard jobs.

    ```bash
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=1
    ```

    This specifies the number of nodes and CPUs required for the job. For a simple small job like FASTQC, we only require 1 compute node. We are also only running 1 task on this node - FASTQC. Finally, we only require 1 CPU for FASTQC. As such, we set all three of these parameters to `1`.

    ```bash
    -l mem=1GB
    ```

    This specifies the amount of memory or RAM we require. As this is a small job, we just require 1GB of memory.

    ```bash
    --time=00:10:00
    ```

    This specifies the maximum amount of time the job can run for. If the job takes any longer than what is specified here, it will be killed, so we want to make sure that we give it enough time to finish.

    ```bash
    fastqc.sh
    ```

    Finally, we provide the path to the script to submit.

    Additionally, note that the script will be run within our current working directory, and that all commands executed within the script will be relative to our current directory.

Go ahead and run the command for your system!

## 1.2.3 Monitoring a job's progress

When running jobs on an HPC, especially long-running jobs, it's useful to be able to monitor its progress and check whether it is queued, running, or complete. All HPC systems will have a way to do exactly this. On the command line, run the following command for your system to see your queue list:

=== "Gadi"

    ```bash
    qstat -u $USER
    ```

    You should see an output something like the following:

    ```console title="Output"
    gadi-pbs: 
                                                                    Req'd  Req'd   Elap
    Job ID               Username Queue    Jobname    SessID NDS TSK Memory Time  S Time
    -------------------- -------- -------- ---------- ------ --- --- ------ ----- - -----
    123456789.gadi-pbs   usr123   normal   fastqc        --    1   1  1024m 00:10 Q   -- 
    ```

=== "Setonix"

    ```bash
    squeue -u $USER
    ```

    You should see an output something like the following:

    ```console title="Output"
    JOBID        USER ACCOUNT                   NAME   EXEC_HOST ST     REASON START_TIME       END_TIME  TIME_LEFT NODES   PRIORITY       QOS
    12345678 username pawsey1234                fastqc       n/a PD       None N/A                   N/A      10:00     1      75423    normal
    ```

At first, the job will be in a **queued** stage, waiting to start. Once there is a suitable slot for the job to run, it will move into a **running** stage, and will continue until it finishes, or (hopefully not) fails. The above script will take only a minute or so to complete, at which point the job will move into a **completed** stage, during which the scheduler will wrap things up, including writing final outputs to the log files. After that, the job will disappear from the queue list.

## 1.2.4 Inspecting the output

Once the FASTQC job has completed, you should see a few new files. First, there will be a new `results/` folder containing the outputs of the `fastqc` command:

```bash
tree results
```

```console title="Output"
results
└── fastqc_gut_logs
    ├── gut_1_fastqc.html
    ├── gut_1_fastqc.zip
    ├── gut_2_fastqc.html
    └── gut_2_fastqc.zip
```

You will also see one or two new files that contain the standard output and standard error streams from the `fastqc` command. This differs slightly between systems:

=== "Gadi"

    On Gadi, the standard output and error streams are split across two files, typically named as `<job name>.[o|e]<job id>`. Note that the error file doesn't necessarily only contain error messages; it can also contain less urgent warning messages and other messages that might otherwise clutter up the main output file. The standard error file looks like this:

    ```bash
    cat fastqc.e123456789
    ```

    ```console title="Output"
    Started analysis of gut_1.fq
    Approx 30% complete for gut_1.fq
    Approx 65% complete for gut_1.fq
    Started analysis of gut_2.fq
    Approx 30% complete for gut_2.fq
    Approx 65% complete for gut_2.fq
    ```

    The standard output file looks like this:

    ```bash
    cat fastqc.o123456789
    ```

    ```console title="Output"
    Analysis complete for gut_1.fq
    Analysis complete for gut_2.fq

    ======================================================================================
                    Resource Usage on 2025-11-18 12:00:00:
    Job Id:             123456789.gadi-pbs
    Project:            ab01
    Exit Status:        0
    Service Units:      0.00
    NCPUs Requested:    1                      NCPUs Used: 1               
                                            CPU Time Used: 00:00:04        
    Memory Requested:   1.0GB                 Memory Used: 144.06MB        
    Walltime requested: 00:10:00            Walltime Used: 00:00:05        
    JobFS requested:    100.0MB                JobFS used: 0B              
    ======================================================================================
    ```

    Note that on Gadi, the standard output file also contains a helpful final summary of the job and its resource usage. This can be very useful for benchmarking purposes and determining the optimal resources to request. Note that the `Exit Status` line shows `0`, which indicates a successful run; if there were errors, a non-zero number will usually show here.

=== "Setonix"

    On Setonix, the standard output and error streams are combined into a single file, much like how you would see both streams print to your terminal if running interactively. The file will be typically named as `slurm-<job id>.out`:

    ```bash
    cat slurm-12345678.out
    ```

    ```console title="Output"
    Started analysis of gut_1.fq
    Approx 30% complete for gut_1.fq
    Approx 65% complete for gut_1.fq
    Analysis complete for gut_1.fq
    Started analysis of gut_2.fq
    Approx 30% complete for gut_2.fq
    Approx 65% complete for gut_2.fq
    Analysis complete for gut_2.fq
    ```

## 1.2.5 Simplifying job submission

You might think that it is a bit of a hassle to write the entire job submission command every time with all of those parameters - and you'd be right! Thankfully, many schedulers have a handy trick that lets you insert those parameters into the script itself, allowing you to significantly shorten and simplify the job submission command. This is done by adding special commented lines at the top of the script, like so:

=== "Gadi"

    ```bash title="fastqc.sh" hl_lines="2-9"
    #!/usr/bin/bash
    #PBS -P ab01
    #PBS -N fastqc
    #PBS -q normal
    #PBS -l ncpus=1
    #PBS -l mem=1GB
    #PBS -l walltime=00:10:00
    #PBS -l storage=scratch/ab01
    #PBS -l wd

    module load fastqc/0.12.1

    SAMPLE_ID="gut"
    READS_1="data/ggal/${SAMPLE_ID}_1.fq"
    READS_2="data/ggal/${SAMPLE_ID}_2.fq"

    mkdir -p "results/fastqc_${SAMPLE_ID}_logs"
    fastqc \
        --outdir "results/fastqc_${SAMPLE_ID}_logs" \
        --format fastq ${READS_1} ${READS_2}
    ```

    Note that when specifying your project with `#PBS -P ab01`, you need to explicity write out your project ID, and cannot use the environemnt variable `$PROJECT`. Similarly, you also need to supply the project ID explicity when specifying the storage resource: `#PBS -l storage=scratch/ab01`.

    Now, you can simply run the following command to submit your job:

    ```bash
    qsub fastqc.sh
    ```

=== "Setonix"

    ```bash title="fastqc.sh" hl_lines="2-9"
    #!/usr/bin/bash
    #SBATCH --account=pawsey1234
    #SBATCH --job-name=fastqc
    #SBATCH --partition=work
    #SBATCH --nodes=1
    #SBATCH --ntasks=1
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=1GB
    #SBATCH --time=00:10:00

    module load fastqc/0.11.9--hdfd78af_1

    SAMPLE_ID="gut"
    READS_1="data/ggal/${SAMPLE_ID}_1.fq"
    READS_2="data/ggal/${SAMPLE_ID}_2.fq"

    mkdir -p "results/fastqc_${SAMPLE_ID}_logs"
    fastqc \
        --outdir "results/fastqc_${SAMPLE_ID}_logs" \
        --format fastq ${READS_1} ${READS_2}
    ```

    Note that when specifying your account with `#SBATCH --acount=pawsey1234`, you need to explicity write out your project ID, and cannot use the environemnt variable `$PAWSEY_PROJECT`.

    Now, you can simply run the following command to submit your job:

    ```bash
    sbatch fastqc.sh
    ```
