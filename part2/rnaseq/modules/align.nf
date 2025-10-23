process ALIGNREADS {
    container "quay.io/biocontainers/star:2.7.11b--h5ca1c30_7"
    memory { 4.GB * Math.ceil(reads_1.size() / 1024 ** 3) * task.attempt }
    time { 1.h * Math.ceil(reads_1.size() / 1024 ** 3) * task.attempt }

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
