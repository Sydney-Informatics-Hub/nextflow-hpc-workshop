process FETCHSRAFASTQS {
    container "ncbi/sra-tools:3.2.1"

    input:
    tuple val(sample_id), val(sra_id)

    output:
    tuple val(sample_id), path("${sra_id}*.fastq.gz")

    script:
    """
    fasterq-dump --split-files --threads $task.cpus $sra_id

    for fq in ${sra_id}*.fastq; do
        gzip -c \$fq > \$fq.gz
    done
    """
}


process DOWNLOADFASTQS {
    container "quay.io/biocontainers/wget:1.21.4"

    input:
    tuple val(sample_id), val(url)

    output:
    tuple val(sample_id), path("${sample_id}.fastq.gz")

    script:
    """
    wget -O $(basename $url) $url

    GZIPPED=\$(file $(basename $url) | grep gzip | wc -l)
    if [ \$GZIPPED -eq 0 ]; then
        gzip -c $(basename $url) > ${sample_id}.fastq.gz
    else
        mv $(basename $url) ${sample_id}.fastq.gz
    fi
    """
}