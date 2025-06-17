process ntm_profiler {
    publishDir "${params.output}/${prefix}/ntm_profiler", mode: 'copy'
    
    input:
    tuple path(read1), path(read2), val(prefix)
    
    output:
    tuple path("${prefix}.results.csv"), val(prefix), emit: report
    tuple path("${prefix}.sub_species.results.csv"), val(prefix), emit: sub_species
    
    script:
    """
    # Run standard NTM-Profiler analysis
    ntm-profiler profile \\
        --read1 ${read1} \\
        --read2 ${read2} \\
        --csv \\
        --prefix ${prefix} \\
        --db_dir ${params.NTM_profiler_db} \\
        -d .
    
    # Run sub-species analysis with external species database
    ntm-profiler profile \\
        --read1 ${read1} \\
        --read2 ${read2} \\
        --csv \\
        --external_species_db sub_ntmdb \\
        --species_only \\
        --prefix ${prefix}.sub_species \\
        -d .
    """
}