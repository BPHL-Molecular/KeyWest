process hmmer {
    publishDir "${params.output}/${prefix}/hmmer", mode: 'copy'
    
    input:
    tuple path(assembly), val(prefix)
    
    output:
    tuple path("${prefix}_hit_card_blast.txt"), val(prefix), emit: summary
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    """
    # Run nhmmer search
    nhmmer \\
        --tblout ${prefix}_hmmer_blast.tblout \\
        -A ${prefix}_hmmer_blast.stockholm \\
        ${params.erm41} \\
        $assembly
    # Check if hits were found
    if grep -q "^[^#]" ${prefix}_hmmer_blast.tblout; then
        echo "Hits found, extracting sequences and running BLAST..."
        
        # Extract hit sequences from Stockholm alignment
        esl-reformat \\
            -o ${prefix}_hit_sequence.fasta \\
            fasta \\
            ${prefix}_hmmer_blast.stockholm
        
        # Run BLASTX against CARD database
        if [ -f "${prefix}_hit_sequence.fasta" ] && [ -s "${prefix}_hit_sequence.fasta" ]; then
            blastx \\
                -query ${prefix}_hit_sequence.fasta \\
                -db ${params.card_database} \\
                -outfmt 6 \\
                -max_target_seqs 10 \\
                -evalue 1e-5 \\
                > ${prefix}_hit_card_blast.txt
        else
            echo "No sequences extracted from Stockholm file"
            touch ${prefix}_hit_card_blast.txt
        fi
    else
        echo "No hits found in nhmmer search"
        touch ${prefix}_hit_card_blast.txt
    fi
    """
    
    stub:
    """
    touch ${prefix}_hit_card_blast.txt
    """
}