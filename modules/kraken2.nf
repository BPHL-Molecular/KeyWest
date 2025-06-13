process kraken2 {
    publishDir "${params.output}/${prefix}/kraken2", mode: 'copy'

  input:
    tuple path(r1), path(r2), val(prefix)
    
    output:
    tuple path("${prefix}_kraken2_report.txt"), val(prefix), emit: report
    tuple path("${prefix}_kraken2_output.txt"), val(prefix), emit: output
    
    script:
    def kraken2_db = params.kraken2_db ?: "/kraken2-db"  // Default or from params
    """
    kraken2 --db ${kraken2_db} \\
            --paired ${r1} ${r2} \\
            --output ${prefix}_kraken2_output.txt \\
            --report ${prefix}_kraken2_report.txt \\
            
    """
}