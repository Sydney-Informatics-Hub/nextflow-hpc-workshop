// pipeline input parameters
params.samplesheet = "$projectDir/samplesheet.csv"
params.ref_prefix = "$projectDir/ref/Homo_sapiens_assembly38.20-22"
params.ref_fasta = "${params.ref_prefix}.fasta"
params.ref_fai = "${params.ref_prefix}.fasta.fai"
params.ref_dict = "${params.ref_prefix}.dict"
params.bwa_index = "$projectDir/ref/BWAIndex"
params.bwa_index_name = "Homo_sapiens_assembly38.20-22.fasta"
params.cohort_name = "cohort"
params.split_n = 4
params.outdir = "results"

process FASTQC {

    tag "fastqc on ${sample_id}"
    container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    path "fastqc_${sample_id}", emit: qc_out

    script:
    """
    mkdir -p "fastqc_${sample_id}"
    fastqc -t $task.cpus --outdir "fastqc_${sample_id}" --format fastq $reads_1 $reads_2
    """

}

process ALIGN {

    container "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:1bd8542a8a0b42e0981337910954371d0230828e-0"
    publishDir "${params.outdir}/alignment"

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)
    tuple val(ref_name), path(bwa_index)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: aligned_bam

    script:
    """
    bwa mem -t $task.cpus -R "@RG\\tID:${sample_id}\\tPL:ILLUMINA\\tPU:${sample_id}\\tSM:${sample_id}\\tLB:${sample_id}\\tCN:SEQ_CENTRE" ${bwa_index}/${ref_name} $reads_1 $reads_2 | samtools sort -O bam -o ${sample_id}.bam
    samtools index ${sample_id}.bam
    """

}

process GENOTYPE {

    container "quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_1"
    publishDir "${params.outdir}/genotyping"

    input:
    tuple val(sample_id), path(bam), path(bai)
    tuple path(ref_fasta), path(ref_fai), path(ref_dict)

    output:
    tuple val(sample_id), path("${sample_id}.g.vcf.gz"), path("${sample_id}.g.vcf.gz.tbi"), emit: gvcf

    script:
    """
    gatk --java-options "-Xmx4g" HaplotypeCaller -R $ref_fasta -I $bam -O ${sample_id}.g.vcf.gz -ERC GVCF
    """

}

process JOINT_GENOTYPE {

    container "quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_1"
    publishDir "${params.outdir}/genotyping"

    input:
    tuple val(cohort_id), path(gvcfs), path(gvcfs_idx)
    tuple path(ref_fasta), path(ref_fai), path(ref_dict)

    output:
    tuple val(cohort_id), path("${cohort_id}.vcf.gz"), path("${cohort_id}.vcf.gz.tbi"), emit: vcf

    script:
    variant_params = gvcfs.collect { f -> "--variant ${f}" }.join(" ")
    """
    gatk --java-options "-Xmx4g" CombineGVCFs -R $ref_fasta $variant_params -O cohort.g.vcf.gz
    gatk --java-options "-Xmx4g" GenotypeGVCFs -R $ref_fasta -V cohort.g.vcf.gz -O cohort.vcf.gz
    """

}

process STATS {

    container "quay.io/biocontainers/bcftools:1.22--h3a4d415_1"
    publishDir "${params.outdir}/genotyping"

    input:
    tuple val(cohort_id), path(cohort_vcf), path(cohort_vcf_idx)

    output:
    path "bcftools_stats.txt", emit: stats_out

    script:
    """
    bcftools stats $cohort_vcf > bcftools_stats.txt
    """

}

process MULTIQC {

    container "quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0"
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path "*"

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", emit: data

    script:
    """
    multiqc .
    """

}

// Add additional processes for implementing a scatter/gather pattern

process SPLIT_FASTQ {

    tag "split fastqs for ${sample_id}"
    container "quay.io/biocontainers/fastp:1.0.1--heae3180_0"

    input:
    tuple val(sample_id), path(reads_1), path(reads_2), val(n)

    output:
    tuple val(sample_id), path("*.${sample_id}.R1.fq"), path("*.${sample_id}.R2.fq"), emit: split_fq

    script:
    """
    fastp -Q -L -A -i $reads_1 -I $reads_2 -o ${sample_id}.R1.fq -O ${sample_id}.R2.fq -s $n
    """

}

process ALIGN_CHUNK {

    container "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:1bd8542a8a0b42e0981337910954371d0230828e-0"

    input:
    tuple val(sample_id), val(chunk), path(reads_1), path(reads_2)
    tuple val(ref_name), path(bwa_index)

    output:
    tuple val(sample_id), val(chunk), path("${sample_id}.${chunk}.bam"), path("${sample_id}.${chunk}.bam.bai"), emit: aligned_bam

    script:
    """
    bwa mem -t $task.cpus -R "@RG\\tID:${sample_id}\\tPL:ILLUMINA\\tPU:${sample_id}\\tSM:${sample_id}\\tLB:${sample_id}\\tCN:SEQ_CENTRE" ${bwa_index}/${ref_name} $reads_1 $reads_2 | samtools sort -O bam -o ${sample_id}.${chunk}.bam
    samtools index ${sample_id}.${chunk}.bam
    """

}

process MERGE_BAMS {

    container "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:1bd8542a8a0b42e0981337910954371d0230828e-0"
    publishDir "${params.outdir}/alignment"

    input:
    tuple val(sample_id), path(bams), path(bais)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: aligned_bam

    script:
    """
    samtools cat ${bams} | samtools sort -O bam -o ${sample_id}.bam
    samtools index ${sample_id}.bam
    """

}

// Define the workflow
workflow {

    // Define the fastqc input channel
    reads = Channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> {
            // def strandedness = row.strandedness ? row.strandedness : 'auto'
            [ row.sample, file(row.fastq_1), file(row.fastq_2) ] 
        }}

    bwa_index = Channel.fromPath(params.bwa_index)
        .map { idx -> [ params.bwa_index_name, idx ] }
        .first()
    ref = Channel.of( [ file(params.ref_fasta), file(params.ref_fai), file(params.ref_dict) ] )
        .first()

    // Run the fastqc step with the reads_in channel
    FASTQC(reads)

    // Perform a scatter/gather workflow by splitting the FASTQs into multiple chunks,
    // aligning each chunk separately, then merging the BAMs together at the end

    // Split FASTQs for each sample
    reads_to_split = reads
        .map { it + [ params.split_n ] }
    SPLIT_FASTQ(reads_to_split)

    // Extract the chunk ID from the split FASTQs and run ALIGN_CHUNK
    split_fqs_r1 = SPLIT_FASTQ.out.split_fq
        .map { sample_id, fqs_1, _fqs_2 -> [ sample_id, fqs_1 ] }
        .transpose()
        .map { sample_id, fq1 -> {
            def chunk_id = fq1.baseName.tokenize(".")[0]
            [ sample_id, chunk_id, fq1 ]
        } }
    split_fqs_r2 = SPLIT_FASTQ.out.split_fq
        .map { sample_id, _fqs_1, fqs_2 -> [ sample_id, fqs_2 ] }
        .transpose()
        .map { sample_id, fq2 -> {
            def chunk_id = fq2.baseName.tokenize(".")[0]
            [ sample_id, chunk_id, fq2 ]
        } }
    split_fqs = split_fqs_r1.join(split_fqs_r2, by: [0, 1])
    ALIGN_CHUNK(split_fqs, bwa_index)

    // Gather all BAM chunks for a sample and run MERGE_BAMS
    gathered_bams = ALIGN_CHUNK.out.aligned_bam
        .groupTuple()
        .map { sample_id, _chunks, bams, bais -> [ sample_id, bams, bais ] }
    MERGE_BAMS(gathered_bams)

    // Run genotyping with aligned bam and genome reference
    GENOTYPE(MERGE_BAMS.out.aligned_bam, ref)

    // Gather gvcfs and run joint genotyping
    all_gvcfs = GENOTYPE.out.gvcf
        .map { _sample_id, gvcf, gvcf_idx -> [ params.cohort_name, gvcf, gvcf_idx ] }
        .groupTuple()
    JOINT_GENOTYPE(all_gvcfs, ref)

    // Get VCF stats
    STATS(JOINT_GENOTYPE.out.vcf)

    // Collect summary data for MultiQC
    multiqc_in = FASTQC.out.qc_out
        .mix(STATS.out.stats_out)
        .collect()

    /*
    * Generate the analysis report with the 
    * outputs from fastqc and bcftools stats
    */ 
    MULTIQC(multiqc_in)

}
