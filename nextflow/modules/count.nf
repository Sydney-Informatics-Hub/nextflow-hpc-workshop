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
    tuple val(sample_id), path(fastqs)
    path star_index

    output:
    path "${sample_id}.bam"

    script:
    """
    # Align reads with STAR
    STAR --genomeDir $star_index \\
         --outSAMtype BAM SortedByCoordinate \\
         --outFileNamePrefix aligned_ \\
         --runThreadN $task.cpus \\
         --readFilesIn $fastqs
    
    # Rename output file to match expected output
    mv aligned_Aligned.sortedByCoord.out.bam ${sample_id}.bam
    """
}

process INDEXBAM {
    container "quay.io/biocontainers/samtools:1.22.1--h96c455f_0"

    input:
    tuple val(sample_id), path(bam, stageAs: "${sample_id}.bam")

    output:
    tuple val(sample_id), path("${sample_id}.bam.bai")

    script:
    """
    samtools index $bam
    """
}

process SPLITBAM {
    container "quay.io/biocontainers/samtools:1.22.1--h96c455f_0"

    input:
    tuple val(sample_id), path(bam, stageAs: "input/*"), path(bai, stageAs: "input/*")

    output:
    path "${sample_id}.*.bam"

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

process COUNT {
    container "quay.io/biocontainers/subread:2.1.1--h577a1d6_0"

    input:
    tuple val(sample_id), path(sortedBam)
    path annotation

    output:
    path "${sample_id}.gene_counts.txt"

    script:
    """
    featureCounts -a $annotation -o ${sample_id}.gene_counts.txt $sortedBam
    """
}

process CONCATCOUNTS {
    container "quay.io/biocontainers/gawk:5.3.1"

    input:
    tuple val(sample_id), path(countsFiles)

    output:
    path "${sample_id}.gene_counts.txt"

    script:
    """
    # Merge count files
    awk 'BEGIN { print "gene", "count" } $0 !~ /^(#|Geneid)/ { print $1, $7 }' ${countsFiles.join(" ")} > ${sample_id}.gene_counts.txt
    """
}