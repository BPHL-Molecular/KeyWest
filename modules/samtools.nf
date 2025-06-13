process samtools {
    publishDir "${params.output}/${prefix}/samtools", mode: 'copy'

    
    input:
    tuple path(sam_file), val(prefix)
    
    output:
    tuple path("${prefix}.sorted.bam"), val(prefix), emit: bam
    
    script:
    """
    # Convert SAM to BAM and sort
    samtools view -@ 1 -bS ${sam_file} | samtools sort -@ 1 -o ${prefix}.sorted.bam
    
    # Index BAM
    samtools index -@ 1 ${prefix}.sorted.bam
    """
}