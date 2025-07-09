process hmmer {
    module 'hmmer'
    module 'ncbi_blast'
    
    publishDir "${params.output}/${prefix}/hmmer", mode: 'copy'
    
    input:
    tuple path(assembly), val(prefix)
    
    output:
    tuple path("${prefix}_hit_card_blast.txt"), val(prefix), emit: summary
    tuple path("${prefix}_hmmer_blast.tblout"), val(prefix), emit: table, optional: true
    tuple path("${prefix}_hmmer_blast.stockholm"), val(prefix), emit: stockholm, optional: true
    tuple path("${prefix}_hit_sequence.fasta"), val(prefix), emit: sequences, optional: true
    tuple path("${prefix}_hmmer_analysis_log.txt"), val(prefix), emit: log, optional: true
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    """
    # Create analysis log file
    exec > >(tee ${prefix}_hmmer_analysis_log.txt) 2>&1
    
    echo "=========================================="
    echo "HMMER RESISTANCE GENE ANALYSIS"
    echo "Sample: ${prefix}"
    echo "Query: ${params.erm41}"
    echo "Assembly: $assembly"
    echo "Started: \$(date)"
    echo "=========================================="
    echo ""
    
    # Run nhmmer search
    echo "Running nhmmer search..."
    nhmmer \\
        --tblout ${prefix}_hmmer_blast.tblout \\
        -A ${prefix}_hmmer_blast.stockholm \\
        ${params.erm41} \\
        $assembly
    
    echo ""
    echo "=== HMMER TABLE HITS ==="
    if [ -f "${prefix}_hmmer_blast.tblout" ]; then
        echo "HMMER table output file: ${prefix}_hmmer_blast.tblout"
        echo "Number of lines: \$(wc -l < ${prefix}_hmmer_blast.tblout)"
        echo ""
        echo "Table content (excluding comments):"
        grep -v "^#" ${prefix}_hmmer_blast.tblout || echo "No hits found"
    else
        echo "No table output file generated"
    fi
    
    echo ""
    echo "=== STOCKHOLM FILE CONTENT ==="
    if [ -f "${prefix}_hmmer_blast.stockholm" ]; then
        echo "Stockholm file: ${prefix}_hmmer_blast.stockholm"
        echo "File size: \$(wc -l < ${prefix}_hmmer_blast.stockholm) lines"
        echo ""
        echo "Stockholm content:"
        cat ${prefix}_hmmer_blast.stockholm
    else
        echo "No Stockholm file generated"
    fi
    
    # Check if hits were found
    if grep -q "^[^#]" ${prefix}_hmmer_blast.tblout; then
        echo ""
        echo "=== SEQUENCE EXTRACTION ==="
        echo "Hits found in nhmmer search, checking Stockholm file..."
        
        # More flexible check for Stockholm sequence content
        if [ -s "${prefix}_hmmer_blast.stockholm" ] && grep -q -E "^[^#/=].*[ATCGN-]" ${prefix}_hmmer_blast.stockholm; then
            echo "Valid sequences found in Stockholm file, extracting..."
            
            # Extract hit sequences from Stockholm alignment
            if esl-reformat \\
                -o ${prefix}_hit_sequence.fasta \\
                fasta \\
                ${prefix}_hmmer_blast.stockholm; then
                
                echo "Sequences extracted successfully!"
                echo ""
                echo "=== EXTRACTED FASTA SEQUENCES ==="
                if [ -f "${prefix}_hit_sequence.fasta" ]; then
                    echo "FASTA file: ${prefix}_hit_sequence.fasta"
                    echo "Number of sequences: \$(grep -c '^>' ${prefix}_hit_sequence.fasta)"
                    echo ""
                    echo "FASTA content:"
                    cat ${prefix}_hit_sequence.fasta
                else
                    echo "No FASTA file created"
                fi
                
                # Verify the FASTA file has actual sequences
                if [ -f "${prefix}_hit_sequence.fasta" ] && [ -s "${prefix}_hit_sequence.fasta" ] && grep -q "^>" ${prefix}_hit_sequence.fasta; then
                    echo ""
                    echo "=== BLAST ANALYSIS ==="
                    echo "Running BLAST on extracted sequences..."
                    
                    blastx \\
                        -query ${prefix}_hit_sequence.fasta \\
                        -db ${params.card_database} \\
                        -outfmt 6 \\
                        -max_target_seqs 10 \\
                        -evalue 1e-5 \\
                        -num_threads ${task.cpus} \\
                        > ${prefix}_hit_card_blast.txt
                    
                    echo "BLAST completed!"
                    echo ""
                    echo "=== BLAST RESULTS ==="
                    if [ -s "${prefix}_hit_card_blast.txt" ]; then
                        echo "BLAST results file: ${prefix}_hit_card_blast.txt"
                        echo "Number of hits: \$(wc -l < ${prefix}_hit_card_blast.txt)"
                        echo ""
                        echo "BLAST output (top 20 hits):"
                        echo "Query_ID\tSubject_ID\tIdentity%\tAlignment_Length\tMismatches\tGaps\tQ_Start\tQ_End\tS_Start\tS_End\tE-value\tBit_Score"
                        head -20 ${prefix}_hit_card_blast.txt
                        echo ""
                        echo "Best hits summary:"
                        head -5 ${prefix}_hit_card_blast.txt | while IFS='\t' read -r qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore; do
                            echo "  - \$sseqid (Identity: \$pident%, E-value: \$evalue)"
                        done
                    else
                        echo "No BLAST hits found (empty results file)"
                    fi
                else
                    echo "No valid sequences in extracted FASTA file"
                    touch ${prefix}_hit_card_blast.txt
                fi
            else
                echo "Failed to extract sequences from Stockholm file (esl-reformat error)"
                touch ${prefix}_hit_card_blast.txt
            fi
        else
            echo "Stockholm file is empty or contains no valid sequences"
            echo "File size: \$(wc -l < ${prefix}_hmmer_blast.stockholm) lines"
            touch ${prefix}_hit_card_blast.txt
        fi
    else
        echo ""
        echo "=== NO HITS FOUND ==="
        echo "No hits found in nhmmer search"
        touch ${prefix}_hit_card_blast.txt
    fi
    
    # Ensure output file always exists for downstream processes
    if [ ! -f "${prefix}_hit_card_blast.txt" ]; then
        echo "Creating empty output file as final safety check"
        touch ${prefix}_hit_card_blast.txt
    fi
    
    # Create empty placeholder files if they don't exist
    if [ ! -f "${prefix}_hit_sequence.fasta" ]; then
        touch ${prefix}_hit_sequence.fasta
    fi
    
    # Final summary
    echo ""
    echo "=========================================="
    echo "ANALYSIS SUMMARY FOR ${prefix}"
    echo "=========================================="
    echo "HMMER table output: ${prefix}_hmmer_blast.tblout (\$(wc -l < ${prefix}_hmmer_blast.tblout 2>/dev/null || echo 0) lines)"
    echo "Stockholm alignment: ${prefix}_hmmer_blast.stockholm (\$(wc -l < ${prefix}_hmmer_blast.stockholm 2>/dev/null || echo 0) lines)"
    echo "Extracted sequences: ${prefix}_hit_sequence.fasta (\$(grep -c '^>' ${prefix}_hit_sequence.fasta 2>/dev/null || echo 0) sequences)"
    echo "BLAST results: ${prefix}_hit_card_blast.txt (\$(wc -l < ${prefix}_hit_card_blast.txt 2>/dev/null || echo 0) hits)"
    echo "Analysis log: ${prefix}_hmmer_analysis_log.txt"
    echo "Completed: \$(date)"
    echo "=========================================="
    """
}