# 2.2 Running nf-core pipelines

!!! info "Learning objectives"

    - Build a basic run script for a nf-core pipeline
    - Run an nf-core pipeline
    - Understand the need for HPC-specific configurations
    
## 2.2.1 Prepraring the environment

As mentioned in the [previous section](./02_1_nfcore_intro.md), we will be running a small section of the overall `sarek` workflow: the `markduplicates` stage. This takes in an alignment `.bam` file and a reference `.fasta` file as input, using a **samplesheet**: a `.csv` file that describes each sample, its name, and where to find its associated `.bam` file.

!!! example "Exercise: review the nf-core/sarek inputts"

    In the parent directory to your current working directory, there is a `data/` folder containing all the input files we need for this workshop. Within the `bams/` sub-directory is a pre-built CSV samplesheet (`samplesheet.csv`). Open it up in VSCode and review the contents:

    ```csv
    patient,sample,bam,bai
    NA12877,test_sample1,../data/bams/NA12877_chr20-22.bam,../data/bams/NA12877_chr20-22.bam.bai
    NA12878,test_sample2,../data/bams/NA12878_chr20-22.bam,../data/bams/NA12878_chr20-22.bam.bai
    NA12889,test_sample3,../data/bams/NA12889_chr20-22.bam,../data/bams/NA12889_chr20-22.bam.bai
    ```

    The first line is a **header line**: it describes the "column names" of the table that the CSV represents. We have four columns in this file:

    - `patient`: An identifier for the individual the sample came from. Note that more than one sample may originate from an individual, such as in the case of paired tumor and normal samples.
    - `sample`: A unique identifier for the sample. Each line must have a unique value for this column.
    - `bam`: The path to the BAM alignment file. The path can be relative to the working directory; we will run this from the `part1` directory, so the relative path is `../data/bams/*.bam`
    - `bai`: The path to the matched index file for the BAM file. Indexes are vital for being able to quickly search a large sequencing dataset for reads at a specific location in the genome.

    We can also see that we have three samples defined in our samplesheet. Samplesheets make it very easy to add new samples to our workflow run, simply by adding a new line and filling out the relevant columns.

    Let's now take a look at the BAM files themselves. These are binary files, so we can't view them directly, although we can use a tool called `samtools` to get a quick glance at their contents:

    === "Gadi"

        ```bash
        module load samtools
        samtools view ../data/bams/NA12878_chr20-22.bam | head -n 3
        ```

    === "Setonix"

        ```bash
        module load samtools/<TODO version>
        samtools view ../data/bams/NA12878_chr20-22.bam | head -n 3
        ```

    This will print out the first three lines of the BAM file for sample `NA12878_chr20-22`:

    ```console title="Output"
    HSQ1004:134:C0D8DACXX:1:1105:9389:51835 163     chr20   70768   6       9S81M1D11M      =       70907   231     TTTGGATGAGGCAATCAGACAAGAGAAAGAAATAAAGATCATTTAAATAGGAAGAGAAGAAGTTAAACTATCCCTGTTGGCAGATGACATATCCTATATCT      @C@FFFFFGHAHHGIJIIJJIJJJJJJJJIJJJJIGJJJIJJIIGGIIEIGHIFJIJJIGGIIIGGHIJJJGGHHHEEBEEFEDEDDEDCCDCDEDDDECC   MD:Z:9G9A8T0G3C0C4G0A7G5C14T6C4^G11        PG:Z:MarkDuplicates     RG:Z:C0D8DACXX.1        NM:i:13 AS:i:25 XS:i:36 MQ:i:22 MC:Z:5S92M4S
    HSQ1004:134:C0D8DACXX:1:1105:9389:51835 83      chr20   70907   22      5S92M4S =       70768   -231    CTGATAACTTCAGCAAAGTCTCAGGATACAAAATCAATGTGCAAAAATCATTAACATTCTTATACACCAACAGTCAAGCCAGGAGCCAAATCAGGAACACA      DDCEEDCEECDEEFFFDFFHGGFEHEIIJJJJIJJIJJIGGGIGJJIIIIIIGDIGIGJIGHF?GCIIHGJJJIJJJJJJJJJJJJJJHHGHHFFFFFCCC   MD:Z:32A3T0G7C2G5C21A15    PG:Z:MarkDuplicates     RG:Z:C0D8DACXX.1        NM:i:7  AS:i:57 XS:i:55 MQ:i:6  MC:Z:9S81M1D11M
    HSQ1004:134:C0D8DACXX:1:2101:9406:116923        163     chr20   81255   60      101M    =       81496   342     AATCATACTAGCTCTCATCAGATTGAAATGGCTGAAATGACAGACATAGAATTAATGATCTGGAAGGTAAGGAAGCTCAAGAACATTCAGGAGAAAGTTGA      ?B<D;:B?2AF?AGHIBACHH@4C3EHH@BGCDDDGHH@FBD>FB?GHHIE4?*9*??BBBFGE8BE)=@@@GA(.7=AE3AD@>);)6@>CDABD(---;   MD:Z:53C47 PG:Z:MarkDuplicates     RG:Z:C0D8DACXX.1        NM:i:1  AS:i:96 XS:i:29 MQ:i:60 MC:Z:101M
    ```

    You can see that these files have a tabular structure as well; briefly, the important columns to note are:

    - 1: The read name
    - 3: The chromosome it aligns to
    - 4: The position on the chromosome it aligns to
    - 5: The mapping quality, or how well the sequence aligns to the genome
    - 9: The length of the aligned sequence (with negative numbers indicating that the sequence aligns to the reverse DNA strand)
    - 10: The sequence itself
    - 11: The per-base quality of the sequence, encoded as ASCII characters

    Finaly, we'll take a look at the FASTA file. This is a plain text file that holds the full sequence of the reference genome that the sequencing data has been aligned to. You can open this file up in the VSCode editor: it's located at `data/ref/Hg38.subsetchr20-22.fasta`. You'll notice that the first line starts with a `>` character, followed by the name of the sequence that follows it; in this case, `chr20`. You'll also notice that the next 600 lines are just full of `N` characters; this is a "mask" character, and is used to essentially block out regions of the genome that are unreliably sequenced, as is the case with the chromosome ends. It isn't until line 602 that we see some actual A's, C's, G's and T's. In total, there are more than half a million lines in the file, representing the sequences of the human chromosomes 20 to 22. We have subset the genome to just these smallest chromosomes to speed up the example analyses in this workshop.

## 2.2.2 Write a simple run script

!!! example "Exercise: Create a run script for nf-core/sarek"

    In the current working directory, create a new blank file and name it `run.sh`. You can do this via VSCode's interface by right-clicking on the current directory (`part1`) in the explorer side bar, clicking on "New File...", writing "run.sh" and hitting the Enter key. You can also do this via the terminal:

    ```bash
    touch run.sh
    ```

    You should also make the file executable by running the following command:

    ```bash
    chmod +x run.sh
    ```

    Next, we'll start the run command by adding the following lines:

    === "Gadi"

        ```bash title="run.sh"
        #!/bin/bash

        module load nextflow/24.04.5
        module load singularity

        nextflow run sarek/main.nf
        ```

    === "Setonix"

        ```bash title="run.sh"
        #!/bin/bash

        module load nextflow/24.10.0
        module load singularity/4.1.0-slurm

        nextflow run sarek/main.nf
        ```

    So far, this is pretty straight-forward: we first tell the OS that we're running a bash script with the `#!/bin/bash` comment. Then, we load the nextflow and singularity modules; these are necessary for ensuring that the `nextflow` command is available to us, and that Nextflow is able to use Singularity for running containers. Finally, we run the `main.nf` Nextflow script within the `sarek` repository.

    As it exists, this won't run: the `sarek` pipeline requires you to give it a samplesheet as input, and we also need to configure a few other parameters to ensure that we only run our small `markduplicates` stage of the pipeline.

    Add a space and a backslash (` \`) to the last line, indicating that the command continues on the following line, then add the following lines to define the inputs:

    ```bash title="run.sh" linenums="6"
    nextflow run sarek/main.nf \
        --input ../data/bams/samplesheet.csv \
        --fasta ../data/ref/Hg38.subsetchr20-22.fasta \
        --fasta_fai ../data/ref/Hg38.subsetchr20-22.fasta.fai \
    ```

    We've provided here the samplesheet CSV file to the `--input` parameter, as well as the FASTA file and its paired index (FAI) file to the `--fasta` and `--fasta_fai` parameters.

    Next, we'll specify that we want to run the pipeline from the `markduplicates` step, and skip several downstream tools, including `baserecalibrator`, `mosdepth`, and `samtools`. This will have the effect of just running the `markduplicates` step, followed by a final run of MultiQC; in total we will run only four distinct processes, which will take only a few minutes for our dataset:

    ```bash title="run.sh" linenums="10"
        --step markduplicates \
        --skip_tools baserecalibrator,mosdepth,samtools \
    ```

    Wrapping up, we will tell the pipeline where to place our outputs, as well as configure a few additional parameters to make things run a bit smoother in our small example:

    ```bash title="run.sh" linenums="12"
        --outdir results \
        --no_intervals true \
        --igenomes_ignore true \
        -resume
    ```

    The `--no_intervals true` parameter tells `sarek` that we don't want to worry about splitting our data up into distinct genomic intervals. This is a very useful feature for large datasets, as it lets us parallelise large genomic sequencing data into chunks that can be processed simultaneously - and takes advantage of the parallel nature of HPCs! However, in our very small example here, it would actually cause things to run slower by adding more jobs and waiting for them to start and finish.

    The next line, `--igenomes_ignore true` stops the pipeline from downloading some additional files from public databases; again, this can be useful in a proper run, but for our purposes it is more of a nuisance.

    Finally, the `-resume` flag tells Nextflow to use the outputs of previously successful tasks where it is appropriate, rather than running them again every time we run the pipeline.
    
    At the end, the `run.sh` script should look like the following:

    === "Gadi"

    ```bash title="run.sh" linenums="1"
    #!/bin/bash

    module load nextflow/24.04.5
    module load singularity

    nextflow run sarek/main.nf \
        --input ../data/bams/samplesheet.csv \
        --fasta ../data/ref/Hg38.subsetchr20-22.fasta \
        --fasta_fai ../data/ref/Hg38.subsetchr20-22.fasta.fai \
        --step markduplicates \
        --skip_tools baserecalibrator,mosdepth,samtools \
        --outdir results \
        --no_intervals true \
        --igenomes_ignore true \
        -resume
    ```

    === "Setonix"

    ```bash title="run.sh" linenums="1"
    #!/bin/bash

    module load nextflow/24.10.0
    module load singularity/4.1.0-slurm

    nextflow run sarek/main.nf \
        --input ../data/bams/samplesheet.csv \
        --fasta ../data/ref/Hg38.subsetchr20-22.fasta \
        --fasta_fai ../data/ref/Hg38.subsetchr20-22.fasta.fai \
        --step markduplicates \
        --skip_tools baserecalibrator,mosdepth,samtools \
        --outdir results \
        --no_intervals true \
        --igenomes_ignore true \
        -resume
    ```

## 2.2.3 Running the pipeline

!!! example "Exercise: Run the workflow"

    We've created a full run script for our pipeline, so let's try it out and see what happens!

    ```bash
    ./run.sh
    ```

    ??? question "Result..."

        You should find that the pipeline quickly fails!

        ```console title="Output"
        TODO: error message
        ```
        
        But why did this happen? You should see in the error output a message that one or more tools couldn't be found:
        
        ```console title="Output"
        TODO: error message
        ```
        
        It's failing because we haven't yet told Nextflow to use Singularity, and each task that runs is trying to find the software installed on the system and is failing to do so.

        There's another problem with our run though: we haven't actually told Nextflow to use the HPC. Instead, it's trying to run everything on the login node: and that is not a good idea at all. Once again, login nodes have limited resources and should never be used for running anything computationally intensive. At best the jobs will fail due to not enough resources being available, and at worst we will slow down the entire login node for everyone else using the system!

        In the next section, we will build up a configuration file that tells Nextflow to both use Singularity and to talk to the HPC's scheduler to submit jobs to the compute nodes