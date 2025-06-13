process ntm_profiler {
    publishDir "${params.output}/${prefix}/ntm_profiler", mode: 'copy'
    
    input:
    tuple path(read1), path(read2), val(prefix)
    
    output:
    tuple path("${prefix}.results.csv"), val(prefix), emit: report
    
    script:
    """
    ntm-profiler profile \\
        --read1 ${read1} \\
        --read2 ${read2} \\
        --csv \\
        --prefix ${prefix} \\
        --db_dir ${params.NTM_profiler_db} \\
        -d .
    """
}