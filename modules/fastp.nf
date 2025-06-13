process fastp {
    publishDir "${params.output}/${prefix}/fastp", mode: 'copy'

    
    input:
    tuple path(r1), path(r2), val(prefix)
    
    output:
    tuple path("${prefix}_1.clean.fastq.gz"), path("${prefix}_2.clean.fastq.gz"), val(prefix), emit: reads
    path("${prefix}_fastp.html"), emit: reports
    
    script:
    """
    fastp -i ${r1} -I ${r2} \\
          -o ${prefix}_1.clean.fastq.gz -O ${prefix}_2.clean.fastq.gz \\
          -h ${prefix}_fastp.html \\
          --detect_adapter_for_pe \\
          --correction --cut_right \\
          --thread ${task.cpus} \\
          --length_required 36
    """
}