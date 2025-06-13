process bwa {
    publishDir "${params.output}/${prefix}/bwa", mode: 'copy'
    
    input:
    tuple path(r1), path(r2), val(prefix)
    
    output:
    tuple path("${prefix}.sam"), val(prefix)
    
    script:
    def reference = params.reference
    """
    # Copy reference to work directory
    cp ${reference} reference.fna.gz
    
    # Index with BWA-MEM2
    bwa-mem2 index reference.fna.gz
    
    # Run BWA-MEM2 alignment
    bwa-mem2 mem reference.fna.gz ${r1} ${r2} > ${prefix}.sam
    """
}