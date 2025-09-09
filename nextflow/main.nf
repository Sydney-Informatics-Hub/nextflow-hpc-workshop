// import processes to be run in the workflow
include { FETCHSRAFASTQS; DOWNLOADFASTQS } from './modules/download'
include { STARINDEX; ALIGNREADS; INDEXBAM; SPLITBAM; SPLITGTF; COUNT } from './modules/count'

// pipeline input parameters
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
params.gtf_file = "$projectDir/data/ggal/transcriptome.gtf"
params.star_index = ""
params.ref_name = params.star_index ? file(params.star_index).baseName : file(params.transcriptome_file).baseName
params.reads = "$projectDir/data/samplesheet.csv"
params.outdir = "results"

/*
 * define the `INDEX` process that creates a binary index
 * given the transcriptome file
 */
process INDEX {

    container "quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0"
    publishDir params.outdir, mode: 'copy'

    input:
    path transcriptome

    output:
    path 'salmon_index'

    script:
    """
    salmon index -t $transcriptome -i salmon_index
    """
}

process FASTQC {

    tag "fastqc on ${sample_id}"
    container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    mkdir -p "fastqc_${sample_id}_logs"
    fastqc --outdir "fastqc_${sample_id}_logs" --format fastq $reads_1 $reads_2 -t $task.cpus
    """
}

process QUANTIFICATION {

    tag "salmon on ${sample_id}"
    container "quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0"
    publishDir params.outdir, mode: 'copy'
    
    input:
    path salmon_index
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    path "$sample_id"

    script:
    def reads_arg = reads_2 ? "-1 ${reads_1} -2 ${reads_2}" : "-r ${reads_1}"
    """
    salmon quant --libType=U -i $salmon_index $reads_arg -o $sample_id
    """
}

process MULTIQC {

    container "quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0"
    publishDir params.outdir, mode: 'copy'

    input:
    path "*"

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc .
    """
}

// Define the workflow
workflow {

    // Get the STAR index
    if (params.star_index == "") {
        // Run the index step with the transcriptome parameter
        ref_data = Channel.of([ params.ref_name, file(params.transcriptome_file), file(params.gtf_file) ])
        STARINDEX(ref_data)
        ref_index = STARINDEX.out.star_index
    } else {
        // Get the user-supplied STAR index instead
        ref_index = Channel.fromPath(params.star_index)
            .first()  // Ensure we have a dataflow value (value channel) rather than a channel (queue channel)
    }

    // Define the fastqc input channel
    reads_in = Channel.fromPath(params.reads)
        .splitCsv(header: true)
        .map { row -> {
            def fastq_2_file = row.fastq_2 == '' ? [] : file(row.fastq_2)
            [row.sample, file(row.fastq_1), fastq_2_file] 
        }}

    // Run the fastqc step with the reads_in channel
    FASTQC(reads_in)

    // Run the quantification step with the index and reads_in channels
    ALIGNREADS(ref_index, reads_in)

    // Index the BAM file
    aligned_reads = ALIGNREADS.out.aligned_bam
    INDEXBAM(aligned_reads)

    // Split the BAM file into one BAM per chromosome
    indexed_reads = aligned_reads.join(INDEXBAM.out.indexed_bam)
    SPLITBAM(indexed_reads)

    // Split the GTF into per-chromosome files
    gtf_file = Channel.fromPath(params.gtf_file)
    SPLITGTF(gtf_file)

    // Prepare the input channel for COUNT
    split_gtfs = SPLITGTF.out.split_gtfs
        .flatten()
        .map { gtf -> {
            def chrom = gtf.baseName
            [ gtf, chrom ]
        }}
    
    split_reads = SPLITBAM.out.split_bams
        .transpose()
        .map { sample_id, bam_file -> {
            def chrom = bam_file.baseName.tokenize('.')[-1]
            [ sample_id, chrom, bam_file ]
        }}
        .combine(split_gtfs, by: 1)
        .map { chrom, sample_id, bam_file, gtf -> {
            [ sample_id, chrom, bam_file, gtf ]
        }}

    // Perform counting
    COUNT(split_reads)

    // Define the multiqc input channel
    count_summaries = COUNT.out.count_summary
        .map { it[1] }
    multiqc_in = FASTQC.out[0]
        .mix(count_summaries)
        .collect()

    /*
    * Generate the analysis report with the 
    * outputs from FASTQC and QUANTIFICATION
    */ 
    MULTIQC(multiqc_in)

}
