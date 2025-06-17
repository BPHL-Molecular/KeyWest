process unicycler {
    publishDir "${params.output}/${prefix}/unicycler", mode: 'copy'
    
    input:
    tuple path(read1), path(read2), val(prefix)
    
    output:
    tuple path("${prefix}_unicycler.fasta"), val(prefix), emit: assembly
    path("${prefix}_unicycler_assembly.gfa"), emit: gfa
    path("${prefix}_unicycler.log"), emit: log
    
    script:
    """
    
    # Run Unicycler for short reads only
    unicycler \\
        -1 ${read1} \\
        -2 ${read2} \\
        -o unicycler_output \\
        --threads ${task.cpus} \\
        --verbosity 2 \\
        --keep 0
    
    # Copy and rename output files
    cp unicycler_output/assembly.fasta ${prefix}_unicycler.fasta
    cp unicycler_output/assembly.gfa ${prefix}_unicycler_assembly.gfa
    cp unicycler_output/unicycler.log ${prefix}_unicycler.log
    
    # Check if assembly was successful
    if [ ! -s ${prefix}_unicycler.fasta ]; then
        echo "ERROR: Unicycler assembly failed or produced empty output"
        exit 1
    fi
    
    # Print assembly stats
    echo "Assembly completed for ${prefix}"
    echo "Assembly file: ${prefix}_unicycler.fasta"
    grep -c ">" ${prefix}_unicycler.fasta | xargs echo "Number of contigs:"
    """
}