process stats {
    publishDir "${params.output}/${prefix}/stats", mode: 'copy'
    
    input:
    tuple path(bam_file), val(prefix)
    
    output:
    path "${prefix}_stats.txt"
    path "${prefix}_flagstat.txt"
    path "${prefix}_coverage.txt"
    path "${prefix}_depth.txt"
    
    script:
    """
    # Generate stats
    samtools stats ${bam_file} > ${prefix}_stats.txt
    samtools flagstat ${bam_file} > ${prefix}_flagstat.txt
    samtools coverage ${bam_file} > ${prefix}_coverage.txt
    samtools depth ${bam_file} > ${prefix}_depth.txt
    """
}