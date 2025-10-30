#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters
params.reference = "ref/Homo_sapiens_assembly38.20-22.fasta"
params.bwa_index = "ref/BWAIndex/Homo_sapiens_assembly38.20-22.fasta"
params.outdir = "results"
params.chunks = 4
params.reads = "data/samplesheet.csv"

process SPLIT_FASTQ {
    tag "$sample_id"
    publishDir "${params.outdir}/fastq_split_${sample_id}", mode: 'copy'
    
    cpus 1
    memory '1 GB'
    time '2m'
    
    input:
    tuple val(sample_id), path(fq1), path(fq2)
    
    output:
    tuple val(sample_id), path("*.R1.fq"), path("*.R2.fq")
    
    script:
    """
    N_LINES=\$(awk -v chunks="${params.chunks}" \\
        'NR % 4 == 1 { nreads += 1 } \\
         END { \\
           if (nreads % chunks == 0) { \\
             chunkreads = nreads / chunks \\
           } else { \\
             chunkreads = (nreads + chunks - (nreads % chunks)) / chunks \\
           }; \\
           printf("%0.f\\n", chunkreads * 4) \\
         }' ${fq1})
    
    split -l \${N_LINES} -d --additional-suffix=.R1.fq ${fq1} ${sample_id}.split_
    split -l \${N_LINES} -d --additional-suffix=.R2.fq ${fq2} ${sample_id}.split_
    """
}

// Process 3: BWA alignment
process ALIGN {
    tag "${sample_id}_${split_id}"
    
    cpus 1
    memory '1 GB'
    time '5m'
    
    input:
    tuple val(sample_id), val(split_id), path(fq1), path(fq2), path(reference), path(bwa_index)
    
    output:
    tuple val(sample_id), path("${split_id}.bam")
    
    script:
    """
    bwa mem -t ${task.cpus} \\
        -R "@RG\\tID:${sample_id}\\tPL:ILLUMINA\\tPU:${sample_id}\\tSM:${sample_id}\\tLB:${sample_id}\\tCN:SEQ_CENTRE" \\
        ${reference} ${fq1} ${fq2} | \\
        samtools view -b -o ${split_id}.bam
    """
}

process MERGE_BAM {
    tag "$sample_id"
    publishDir "${params.outdir}/merge_${sample_id}", mode: 'copy'
    
    cpus 1
    memory '1 GB'
    time '5m'
    
    input:
    tuple val(sample_id), path(bams)
    
    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai")
    
    script:
    """
    samtools cat ${bams} | samtools sort -O bam -o ${sample_id}.bam
    samtools index ${sample_id}.bam
    """
}

process GENOTYPE {
    tag "$sample_id"
    publishDir "${params.outdir}/geno_${sample_id}", mode: 'copy'
    
    cpus 4
    memory '4 GB'
    time '10m'
    
    input:
    tuple val(sample_id), path(bam), path(bai), path(reference), path(ref_index), path(ref_dict)
    
    output:
    tuple val(sample_id), path("${sample_id}.g.vcf.gz"), path("${sample_id}.g.vcf.gz.tbi")
    
    script:
    """
    gatk --java-options "-Xmx4g" HaplotypeCaller \\
        -R ${reference} \\
        -I ${bam} \\
        -O ${sample_id}.g.vcf.gz \\
        -ERC GVCF
    """
}

process JOINT_GENOTYPE {
    publishDir "${params.outdir}/joint_geno", mode: 'copy'
    
    cpus 4
    memory '4 GB'
    time '10m'
    
    input:
    path(gvcfs)
    path(gvcf_indices)
    path(reference)
    path(ref_index)
    path(ref_dict)
    
    output:
    tuple path("cohort.vcf.gz"), path("cohort.vcf.gz.tbi")
    
    script:
    def gvcf_args = gvcfs.collect { "--variant $it" }.join(' ')
    """
    gatk --java-options "-Xmx4g" CombineGVCFs \\
        -R ${reference} \\
        ${gvcf_args} \\
        -O cohort.g.vcf.gz
    
    gatk --java-options "-Xmx4g" GenotypeGVCFs \\
        -R ${reference} \\
        -V cohort.g.vcf.gz \\
        -O cohort.vcf.gz
    """
}

process VCF_STATS {
    publishDir "${params.outdir}/vcf_stats", mode: 'copy'
    
    cpus 1
    memory '1 GB'
    time '5m'
    
    input:
    tuple path(vcf), path(vcf_index)
    
    output:
    path "bcftools_stats.txt"
    
    script:
    """
    bcftools stats ${vcf} > bcftools_stats.txt
    """
}

process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    
    cpus 1
    memory '1 GB'
    time '5m'
    
    input:
    path(fastqc_files)
    path(stats_file)
    
    output:
    path "multiqc_report.html"
    path "multiqc_data"
    
    script:
    """
    multiqc .
    """
}

workflow {

    // read in fastq reads from samplesheet
    reads_in = Channel.fromPath(params.reads)
        .splitCsv(header: true)
        .map {
            row -> [row.sample, file(row.fastq_1), file(row.fastq_2)]
        }
    
    // Reference files
    reference = Channel.fromPath(params.reference)
    ref_index = Channel.fromPath("${params.reference}.fai")
    ref_dict = Channel.fromPath(params.reference.replaceAll(/\.fasta$/, '.dict'))
    bwa_index = Channel.fromPath("${params.bwa_index}*").collect()
   
    //FASTQC(reads_in)
    
    // Split FASTQ files
    SPLIT_FASTQ(reads_in)
    
    // Prepare alignment inputs by splitting read pairs
    split_pairs = SPLIT_FASTQ.out
        .transpose()
        .map { sample_id, fq1s, fq2s -> tuple(sample_id, fq1s, fq2s) }
        .transpose()
        .map { sample_id, fq1, fq2 -> 
            def split_id = fq1.baseName.replaceAll(/\.R1$/, '')
            tuple(sample_id, split_id, fq1, fq2)
        }
        .combine(Channel.of(reference))
        .combine(bwa_index)
    
    // Align reads
    ALIGN(split_pairs)
    
    // Merge BAM files per sample
    merged_bams = ALIGN.out
        .groupTuple()
        .map { sample_id, bams -> tuple(sample_id, bams) }
    
    MERGE_BAM(merged_bams)
    
    // Call variants per sample
    geno_input = MERGE_BAM.out
        .combine(Channel.of(reference))
        .combine(Channel.of(ref_index))
        .combine(Channel.of(ref_dict))
    
    GENOTYPE(geno_input)
    
    // Collect all GVCFs
    all_gvcfs = GENOTYPE.out.map { sample_id, gvcf, idx -> gvcf }.collect()
    all_indices = GENOTYPE.out.map { sample_id, gvcf, idx -> idx }.collect()
    
    // Joint genotyping
    JOINT_GENOTYPE(
        all_gvcfs,
        all_indices,
        reference,
        ref_index,
        ref_dict
    )
    
    // Calculate stats
    VCF_STATS(JOINT_GENOTYPE.out)
    
    // Generate MultiQC report
    MULTIQC(
        FASTQC.out.collect(),
        VCF_STATS.out
    )
}
