process STARINDEX {
    container "quay.io/biocontainers/star:2.7.11b--h5ca1c30_7"

    input:
    tuple val(ref_name), path(reference), path(annotation)

    output:
    path "${ref_name}", emit: star_index

    script:
    """
    # Create STAR index directory
    mkdir -p ${ref_name}
    
    # Generate STAR genome index
    STAR --runMode genomeGenerate \\
         --genomeDir ${ref_name} \\
         --genomeFastaFiles $reference \\
         --runThreadN $task.cpus \\
         --sjdbGTFfile $annotation
    """
}

process ALIGNREADS {
    container "quay.io/biocontainers/star:2.7.11b--h5ca1c30_7"

    input:
    path star_index
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: aligned_bam

    script:
    """
    # Align reads with STAR
    STAR --genomeDir $star_index \\
         --outSAMtype BAM SortedByCoordinate \\
         --outFileNamePrefix aligned_ \\
         --runThreadN $task.cpus \\
         --readFilesIn $reads_1 $reads_2
    
    # Rename output file to match expected output
    mv aligned_Aligned.sortedByCoord.out.bam ${sample_id}.bam
    """
}

process INDEXBAM {
    container "quay.io/biocontainers/samtools:1.22.1--h96c455f_0"

    input:
    tuple val(sample_id), path(bam, stageAs: "input/*")

    output:
    tuple val(sample_id), path("${sample_id}.bam.bai"), emit: indexed_bam

    script:
    """
    samtools index -o ${sample_id}.bam.bai $bam
    """
}

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

    input:
    tuple val(sample_id), val(chr), path(sortedBam), path(annotation)

    output:
    tuple val(sample_id), path("${sample_id}.${chr}.gene_counts.txt"), emit: gene_counts
    tuple val(sample_id), path("${sample_id}.${chr}.gene_counts.txt.summary"), emit: count_summary

    script:
    """
    featureCounts -a $annotation -p --countReadPairs -o ${sample_id}.${chr}.gene_counts.txt $sortedBam
    """
}