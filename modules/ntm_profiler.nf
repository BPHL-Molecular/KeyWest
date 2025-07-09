process ntm_profiler {
    errorStrategy 'ignore'  // Add this line - ignore failures and continue
    
    publishDir "${params.output}/${prefix}/ntm_profiler", mode: 'copy'
    
    input:
    tuple path(assembly), val(prefix)  // Fixed: Added 'path()' declaration
    
    output:
    tuple path("${prefix}.results.csv"), val(prefix), emit: report, optional: true
    tuple path("${prefix}.sub_species.results.csv"), val(prefix), emit: sub_species, optional: true
    tuple path("${prefix}.vcf.gz"), val(prefix), emit: vcf, optional: true
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    """
    # Run standard NTM-Profiler analysis
    ntm-profiler profile \\
        --fasta ${assembly} \\
        --csv \\
        --prefix ${prefix} \\
        --db_dir ${params.NTM_profiler_db} \\
        -d .
    
    # Check for 'multiple species' message
    if grep -q 'Multiple species found' "${prefix}.results.json"; then
        echo "[INFO] Multiple species detected for ${prefix}. Extracting top hit and rerunning with --species_filter"
    
        top_species=\$(awk -F',' 'NR==2 {print \$1}' "${prefix}.results.csv")
    
        ntm-profiler profile \\
            --fasta ${assembly} \\
            --csv \\
            --resistance_db Mycobacterium_abscessus \\
            --prefix ${prefix} \\
            --db_dir ${params.NTM_profiler_db} \\
            -d .
    fi
    
    # Run sub-species analysis with external species database using FASTA
    ntm-profiler profile \\
        --fasta ${assembly} \\
        --csv \\
        --external_species_db sub_ntmdb \\
        --species_only \\
        --prefix ${prefix}.sub_species \\
        -d .
    """
}