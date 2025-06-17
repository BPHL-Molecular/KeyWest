process stats {
    publishDir "${params.output}/${prefix}/stats", mode: 'copy'
    
    input:
    tuple path(bam_file), val(prefix)
    
    output:
    tuple path("${prefix}_stats.txt"), val(prefix), emit: stats
    tuple path("${prefix}_flagstat.txt"), val(prefix), emit: flagstat
    tuple path("${prefix}_coverage.txt"), val(prefix), emit: coverage
    tuple path("${prefix}_depth.txt"), val(prefix), emit: depth
    
    script:
    """
    # Generate stats
    samtools stats ${bam_file} > ${prefix}_stats.txt
    samtools flagstat ${bam_file} > ${prefix}_flagstat.txt
    samtools coverage ${bam_file} > ${prefix}_coverage.txt
    samtools depth ${bam_file} > ${prefix}_depth.txt
    """
}