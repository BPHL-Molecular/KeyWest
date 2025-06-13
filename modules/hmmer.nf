// modules/trun_label.nf
process trun_label {
    publishDir "${params.output}/${prefix}/trun_label", mode: 'copy'
    
    input:
    tuple path(prokka_faa), val(prefix)
    
    output:
    tuple path("${prefix}_trun_label_results.txt"), val(prefix), emit: results
    tuple path("${prefix}_hmmer_output.txt"), val(prefix), emit: hmmer_output
    
    script:
    def trun_label_script = params.trun_label_script ?: "trun_label.py"
    """
    # Load required modules
    module load hmmer
    module load python
    
    # Run HMMER search
    hmmsearch \\
        --tblout ${prefix}_hmmer_output.txt \\
        --domtblout ${prefix}_hmmer_dom.txt \\
        --cpu ${task.cpus} \\
        ${params.hmm_database} \\
        ${prokka_faa} > ${prefix}_hmmer_search.out
    
    # Run Python script for trun-label analysis
    python ${trun_label_script} \\
        --hmmer_output ${prefix}_hmmer_output.txt \\
        --input_faa ${prokka_faa} \\
        --prefix ${prefix} \\
        --output ${prefix}_trun_label_results.txt
    """
}