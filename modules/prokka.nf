process prokka {
    publishDir "${params.output}/${prefix}/prokka", mode: 'copy'
    container 'staphb/prokka:1.14.6'
    
    input:
    tuple path(assembly_file), val(prefix)
    
    output:
    tuple path("${prefix}.*"), val(prefix), emit: prokka_results
    
    script:
    """
    # Run Prokka for genome annotation (output directly to current directory)
    prokka --outdir . --prefix ${prefix} --force --centre X --compliant ${assembly_file}
    """
}