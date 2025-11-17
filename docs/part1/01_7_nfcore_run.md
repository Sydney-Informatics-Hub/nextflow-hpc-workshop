# 1.7 Running nf-core on HPC

!!! info "Learning objectives"

    - Build a basic run script for a nf-core pipeline
    - Run an nf-core pipeline
    - Understand the need for HPC-specific configurations
    
## 1.7.1 An overview of the inputs

As mentioned in the [previous section](./02_1_nfcore_intro.md), we will be running a small section of the overall `sarek` workflow: the `mapping` stage. This stage requires a few distinct inputs:

- FASTQ files: These are plain text files that contain the sequences of the reads that came from the sequencing machine.
- A **samplesheet**: A `.csv` file that describes each sample, its name, and where to find its associated FASTQ files.
- A FASTA file: This holds the actual sequence of the chromosomes in the reference genome.
- A BWA index: This is a group of files that is used by the alignment tool `bwa mem` to efficiently find where in the genome a particular sequencing read best maps to.

We have two samplesheets prepared for today's workshop, both within the `../data/fqs` directory. The first is `samplesheet.fq.csv`. You can display it on the command line with:

```bash
cat ../data/fqs/samplesheet.fq.csv
```

!!! note "Viewing the samplesheet in VSCode"

    Unfortunately, because the `data/` folder is in the parent directory to where we are working, we can't see it in the VSCode explorer. However, if you wish to open the samplesheet in the VSCode editor, you can still do so via the terminal using the `code` command:

    ```bash
    code ../data/fqs/samplesheet.fq.csv
    ```

The samplesheet looks like this:

```csv title="samplesheet.fq.csv"
patient,sample,lane,fastq_1,fastq_2
NA12877,test_sample1,all,../data/fqs/NA12877_chr20-22.R1.fq.gz,../data/fqs/NA12877_chr20-22.R2.fq.gz
NA12878,test_sample2,all,../data/fqs/NA12878_chr20-22.R1.fq.gz,../data/fqs/NA12878_chr20-22.R2.fq.gz
NA12889,test_sample3,all,../data/fqs/NA12889_chr20-22.R1.fq.gz,../data/fqs/NA12889_chr20-22.R2.fq.gz
```

The first line is a **header line**: it describes the "column names" of the table that the CSV represents. We have four columns in this file:

- `patient`: An identifier for the individual the sample came from. Note that more than one sample may originate from an individual, such as in the case of paired tumor and normal samples.
- `sample`: A unique identifier for the sample. Each line must have a unique value for this column.
- `lane`: This identifies the lane of the flow cell that the sequencing read originated from.
- `fastq_1`: The path to the read 1 FASTQ file for that sample.
- `fastq_2`: The path to the read 2 FASTQ file for that sample.

We can also see that we have three samples defined in our samplesheet. Samplesheets make it very easy to add new samples to our workflow run, simply by adding a new line and filling out the relevant columns.

The second samplesheet we have prepared is `samplesheet.single.csv`:

```csv title="samplesheet.single.csv"
patient,sample,lane,fastq_1,fastq_2
NA12877,test_sample1,all,../data/fqs/NA12877_chr20-22.R1.fq.gz,../data/fqs/NA12877_chr20-22.R2.fq.gz
```

As you can see, this contains just a single sample. For the inital part of this lesson, we will work with this samplesheet in order to keep things simple and quick.

!!! note "Samplesheet columns for sarek"

    The `sarek` pipeline can use a lot of different tools and start and stop at various stages. This means there are different columns required in the samplesheet depending on what you are running. Check with the [usage documentation](https://nf-co.re/sarek/3.5.0/docs/usage/) to determine the correct columns for your run. In our case, we are using the [minimal configuration](https://nf-co.re/sarek/3.5.0/docs/usage/#start-with-mapping---step-mapping-default) when starting with the `mapping` stage.

??? example "Additional content: A deeper look at the input files"

    For those who are interested, here we look a bit closer at the structure of these input files.

    ### FASTQ files
    
    FASTQ files are plain text files that are structured in a particular way to capture the sequences of the reads in our sample as well as the quality of those sequences. Each read is represented by **four lines** in a FASTQ file:
    
    ```title="Two reads from a FASTQ file"
    @HSQ1004:134:C0D8DACXX:1:1105:9389:51835
    TGTGTTCCTGATTTGGCTCCTGGCTTGACTGTTGGTGTATAAGAATGTTAATGATTTTTGCACATTGATTTTGTATCCTGAGACTTTGCTGAAGTTATCAG
    +
    CCCFFFFFHHGHHJJJJJJJJJJJJJJIJJJGHIICG?FHGIJGIGIDGIIIIIIJJGIGGGIJJIJJIJJJJIIEHEFGGHFFDFFFEEDCEECDEECDD
    @HSQ1004:134:C0D8DACXX:1:2101:9406:116923
    AGGTTTTATTTTTAAATTTTTTTCTTTAGTTTTGTCTGACTGTGTTGATTCAAAGGACTGATCTTTGAGCTCTGAGATTCTTTCCTCAACTTGATCTATTC
    +
    @@?DBDDEBFFFFHIIIIDGHIIGEIIIG9DGGGIIGIIHGHFCBFC>GHGHGHCGEGGGIIEGIG;@=EHCEE@?:BEDADE@@@@;@ACACDCC>CD;@
    ```

    Each set of four lines is defined as follows:

    1. The read name. This often includes information about the flow cell it came from.
    2. The sequence of the read. This is the actual As, Cs, Gs and Ts of the read.
    3. A `+` symbol. This line isn't normally used for anything and serves mostly as a visual separator between the sequence and its quality scores.
    4. The quality scores for the read. The characters on this line represent quality scores for the corresponding position in the sequence. The quality scores are encoded as [ASCII characters](https://www.ascii-code.com/), starting at ASCII code 33 (`!`) for a value of `0` and increasing from there. The quality scores themselves are called [**Phred quality scores**](https://en.wikipedia.org/wiki/Phred_quality_score), which are logarithmically related to the probability of an error at that base (higher scores being of better quality). This encoding scheme is referred to as **Phred-33**.

    We are working with **paired-end** sequencing data, where every sequencing read has two **mates** originating from either end of the original DNA fragment. These are typically separated into two ordered FASTQ files, one containing the first mates (or "read 1") and the other containing the second mates (or "read 2").
    
    FASTQ files are also often compressed as **gzip** files. This helps to reduce storage space for large sequencing samples.

    ### FASTA files
    
    A FASTA file is another a plain text file type that holds the full sequence of the reference genome that the sequencing data will be aligned to. The FASTA format is much simpler than FASTQ. Each sequence starts with a line beginning with a `>` character. This line provides the name for the sequence, e.g. `chr20`. The actual sequence starts on the following line and can be split into any number of lines. The sequence ends when either a new one begins (with a new line starting with `>`) or the file ends.
    
    ```title="An example FASTA file"
    >seq1
    ATCGGTCA
    >seq2
    AATTGGCC
    GGCCTTAA
    >seq3
    AATGGATA
    CTTGATAC
    CCTG
    ```
    
    The FASTA file we are working with today (`Hg38.subsetchr20-22.fasta`) contains three sequences: human chromosomes 20-22. These are the smallest autosomal chromosomes, and were picked for this workshop to ensure speedy alignment. You'll notice that the first 600 lines of the `chr20` sequence are just full of `N` characters; this is a "mask" character, and is used to essentially block out regions of the genome that are unreliably sequenced, as is the case with the chromosome ends. It isn't until line 602 that we see some actual A's, C's, G's and T's. In total, there are more than half a million lines in the file.

    ```title="Hg38.subsetchr20-22.fasta"
    >chr20
    NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    ...
    ATTATCCATGAGGAGGGAGTGCAGACAAAGCAAAGAAGGAGGATGTTTGGAGAGGGGTAGTCTTTGAGTGGAGCCTTTAGGGATGAGAAGGGTGAATTGA
    GATATACCGGGAAAGTAGAAAAGATAAAACACGTATTAAAT
    >chr21
    NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    ```

    The FASTA file is typically associated with two additional helper files: the FASTA index (`.fai`) and dictionary (`.dict`) files. Both contain information about the sequences in the main FASTA file and are used by tools for quickly finding sequences within it, especially when the genome is very large.

    ### BWA Index

    The job of `bwa mem` is to take sequencing reads in the FASTQ format and try to align them to a reference genome. For all but the most trivially-sized genomes, this is a very computationally-intensive task. Therefore, `bwa mem` makes use of an index to more efficiently look up where in the genome specific sequences map to. BWA's index format consists of a directory containing multiple files:

    ```title="An example BWA index"
    Homo_sapiens_assembly38.20-22.fasta.amb  Homo_sapiens_assembly38.20-22.fasta.pac
    Homo_sapiens_assembly38.20-22.fasta.ann  Homo_sapiens_assembly38.20-22.fasta.sa
    Homo_sapiens_assembly38.20-22.fasta.bwt
    ```

    Note that different alignment tools will use their own specific index format; this format is specific to `bwa mem` and can't be used with other tools.

## 1.7.2 Write a simple run script

To run the mapping stage of `sarek`, we need to execute the `nextflow run` command and include several parameters. We will use a run script to compose this command bit-by-bit to keep things organised.

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

    === "Gadi (PBSpro)"

        ```bash title="run.sh"
        #!/bin/bash

        module load nextflow/24.04.5
        module load singularity

        nextflow run sarek/main.nf \
        ```

    === "Setonix (Slurm)"

        ```bash title="run.sh"
        #!/bin/bash

        module load nextflow/24.10.0
        module load singularity/4.1.0-slurm

        nextflow run sarek/main.nf \
        ```

    So far, this is pretty straight-forward: we first tell the OS that we're running a bash script with the `#!/bin/bash` comment. Then, we load the nextflow and singularity modules; these are necessary for ensuring that the `nextflow` command is available to us, and that Nextflow is able to use Singularity for running containers. Finally, we run the `main.nf` Nextflow script within the `sarek` repository.

    As it exists, this won't run: the `sarek` pipeline requires you to give it a samplesheet as input, and we also need to configure a few other parameters to ensure that we only run our small `mapping` stage of the pipeline. That's why we've ended the `nextflow run` line with a backslash (` \`), indicating that the command continues on the following line.
    
    Add the following lines to define the inputs to the `mapping` stage:

    ```bash title="run.sh" linenums="7"
        --input ../data/fqs/samplesheet.single.csv \
        --fasta ../data/ref/Hg38.subsetchr20-22.fasta \
        --fasta_fai ../data/ref/Hg38.subsetchr20-22.fasta.fai \
        --dict ../data/ref/Hg38.subsetchr20-22.dict \
        --bwa ../data/ref \
    ```

    We've provided here the samplesheet CSV file to the `--input` parameter, as well as the FASTA file and its paired index (FAI) and dictionary files to the `--fasta`, `--fasta_fai`, and `--dict` parameters. We also specified the path to the BWA index directory with `--bwa`.

    !!! note "Finding more information about the sarek parameters"

        The [sarek parameters documentation](https://nf-co.re/sarek/3.5.0/parameters/) contains information about all of the parameters supported by `sarek` and the valid arguments to them. The [usage page](https://nf-co.re/sarek/3.5.0/docs/usage/) in the docs also contains useful information about which parameters are required by the various stages of the pipeline.

    Next, we'll specify that we want to run the pipeline from the `mapping` step, and skip several downstream tools, including `markduplicates`, `baserecalibrator`, `mosdepth`, and `samtools`. This will have the effect of just running the `mapping` step, followed by a final run of MultiQC; in total we will run only eight distinct processes, which will take only a few minutes for our dataset:

    ```bash title="run.sh" linenums="12"
        --step mapping \
        --skip_tools markduplicates,baserecalibrator,mosdepth,samtools \
    ```

    Wrapping up, we will tell the pipeline where to place our outputs, as well as configure a few additional parameters to make things run a bit smoother in our small example:

    ```bash title="run.sh" linenums="14"
        --outdir results \
        --no_intervals true \
        --igenomes_ignore true
    ```

    The `--no_intervals true` parameter tells `sarek` that we don't want to worry about splitting our data up into distinct genomic intervals. Intervals are a very useful feature for large datasets, as it lets us parallelise large genomic sequencing data into chunks that can be processed simultaneously - and takes advantage of the parallel nature of HPCs! However, in our very small example here, it would actually cause things to run slower by adding more jobs and waiting for them to start and finish.

    The next line, `--igenomes_ignore true` stops the pipeline from downloading some additional files from public databases; again, this can be useful in a proper run, but for our purposes it is more of a nuisance.

    At the end, the `run.sh` script should look like the following:

    === "Gadi (PBSpro)"

        ```bash title="run.sh" linenums="1"
        #!/bin/bash

        module load nextflow/24.04.5
        module load singularity

        nextflow run sarek/main.nf \
            --input ../data/fqs/samplesheet.single.csv \
            --fasta ../data/ref/Hg38.subsetchr20-22.fasta \
            --fasta_fai ../data/ref/Hg38.subsetchr20-22.fasta.fai \
            --dict ../data/ref/Hg38.subsetchr20-22.dict \
            --bwa ../data/ref \
            --step mapping \
            --skip_tools markduplicates,baserecalibrator,mosdepth,samtools \
            --outdir results \
            --no_intervals true \
            --igenomes_ignore true
        ```

    === "Setonix (Slurm)"

        ```bash title="run.sh" linenums="1"
        #!/bin/bash

        module load nextflow/24.10.0
        module load singularity/4.1.0-slurm

        nextflow run sarek/main.nf \
            --input ../data/fqs/samplesheet.single.csv \
            --fasta ../data/ref/Hg38.subsetchr20-22.fasta \
            --fasta_fai ../data/ref/Hg38.subsetchr20-22.fasta.fai \
            --dict ../data/ref/Hg38.subsetchr20-22.dict \
            --bwa ../data/ref \
            --step mapping \
            --skip_tools markduplicates,baserecalibrator,mosdepth,samtools \
            --outdir results \
            --no_intervals true \
            --igenomes_ignore true
        ```

!!! question "How are you going?"

    If you're following along so far, let us know by reacting on zoom with a **":material-check:{ .check } Yes"**.
    
    If you're running into any issues, please react with a **":material-close:{ .close } No"** and we can help out before we move on to the next section.

## 1.7.3 Not quite ready!

We've created a full run script for our pipeline, but we won't run it just yet, because we haven't actually told Nextflow to use the HPC scheduler. If we were to run the workflow, it would attempt to run on the login node where we are currently working!

In addition, we haven't configured Nextflow to use Singularity yet, so Nextflow will assume all of the necessary software is installed and ready to use - which it isn't.

In the next section, we will build up a configuration file that tells Nextflow to both use Singularity and to talk to the HPC's scheduler to submit jobs to the compute nodes.