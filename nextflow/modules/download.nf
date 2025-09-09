process DOWNLOADGTF {
    container "quay.io/biocontainers/wget:1.21.4"

    input:
    val url

    output:
    path "annotation.gtf", emit: gtf_file

    script:
    gzipped = url.endsWith('.gz')
    dlpath = gzipped ? "annotation.gtf.gz" : "annotation.gtf"
    """
    wget -O $dlpath $url

    if $gzipped
    then
        gunzip -c $dlpath > annotation.gtf
    fi
    """
}
