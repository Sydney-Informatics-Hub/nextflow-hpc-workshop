process SPLITBAM {
    container "quay.io/biocontainers/samtools:1.22.1--h96c455f_0"

    input:
    tuple val(sample_id), path(bam, stageAs: "input/*"), path(bai, stageAs: "input/*")

    output:
    tuple val(sample_id), path("${sample_id}.*.bam"), emit: split_bams

    script:
    """
    # Get list of chromosomes/contigs
    samtools idxstats $bam | cut -f1 | grep '^chr' > chromosomes.txt
    
    # Split BAM by chromosome
    while read chr; do
        samtools view -b $bam "\$chr" > ${sample_id}.\$chr.bam
    done < chromosomes.txt
    """
}

process SPLITGTF {
    container "quay.io/biocontainers/gawk:5.3.1"

    input:
    path annotation, stageAs: "input/*"

    output:
    path "*.gtf", emit: split_gtfs

    script:
    """
    # Split GTF by chromosome in a single pass
    awk '
    /^chr/ {
        chr = \$1
        print \$0 > chr ".gtf"
    }' $annotation
    """
}

process COUNT {
    container "quay.io/biocontainers/subread:2.1.1--h577a1d6_0"
    memory { 4.GB * Math.ceil(sortedBam.size() / 1024 ** 3) * task.attempt }
    time { 1.h * Math.ceil(sortedBam.size() / 1024 ** 3) * task.attempt }

    input:
    tuple val(sample_id), val(chr), path(sortedBam), path(annotation)

    output:
    tuple val(sample_id), path("${sample_id}.${chr}.gene_counts.txt"), emit: gene_counts
    tuple val(sample_id), path("${sample_id}.${chr}.gene_counts.txt.summary"), emit: count_summary

    script:
    """
    featureCounts -a $annotation -p --countReadPairs -T $task.cpus -o ${sample_id}.${chr}.gene_counts.txt $sortedBam
    """
}

process COMBINECHROMCOUNTS {
    container "quay.io/biocontainers/pandas:2.2.1"

    input:
    tuple val(sample_id), path(count_files, stageAs: 'counts/*'), path(summary_files, stageAs: 'summaries/*')

    output:
    tuple val(sample_id), path("${sample_id}.gene_counts.txt"), emit: gene_counts
    tuple val(sample_id), path("${sample_id}.gene_counts.txt.summary"), emit: count_summary

    script:
    """
    combine_featurecounts.py \\
        --sample ${sample_id} \\
        --counts counts/* \\
        --summaries summaries/* \\
        --output-counts ${sample_id}.gene_counts.txt.summary \\
        --output-sumary ${sample_id}.gene_counts.txt.summary
    """
}