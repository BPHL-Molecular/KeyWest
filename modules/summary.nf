process summary {
    publishDir "${params.output}/${prefix}/summary", mode: 'copy'
    
    input:
    tuple path(assembly), path(coverage_file), path(ntm_report), path(ntm_sub_species), path(kraken_report), path(hmmer_blast), val(prefix)
    
    output:
    tuple path("${prefix}_sample_summary.csv"), val(prefix), emit: summary
    
    script:
    """
    # Calculate assembly stats
    awk '/^>/ {
        if(seq) {
            lengths[++count] = length(seq)
            total += length(seq)
        }
        seq=""
        next
    } 
    {seq = seq \$0} 
    END {
        if(seq) {
            lengths[++count] = length(seq)
            total += length(seq)
        }
        
        # Sort lengths in descending order
        for(i=1; i<=count; i++) {
            for(j=i+1; j<=count; j++) {
                if(lengths[i] < lengths[j]) {
                    temp = lengths[i]
                    lengths[i] = lengths[j]
                    lengths[j] = temp
                }
            }
        }
        
        # Calculate N50
        target50 = total * 0.5
        sum = 0
        n50 = 0
        for(i=1; i<=count; i++) {
            sum += lengths[i]
            if(sum >= target50 && n50 == 0) {
                n50 = lengths[i]
                break
            }
        }
        
        # Store values
        printf "%d,%d,%d,%d,%d,%d\\n", count, total, lengths[1], lengths[count], int(total/count), n50 > "assembly_stats.tmp"
        
    }' ${assembly}
    
    # Parse coverage stats from samtools coverage file (skip header)
    tail -n +2 ${coverage_file} | awk -F'\\t' '
    BEGIN {
        total_length = 0
        total_covered_bases = 0
        weighted_depth = 0
        num_contigs = 0
    }
    {
        num_contigs++
        seq_length = \$3 - \$2 + 1
        covbases = \$5
        coverage = \$6
        meandepth = \$7
        
        total_length += seq_length
        total_covered_bases += covbases
        weighted_depth += (meandepth * seq_length)
    }
    END {
        overall_coverage = (total_covered_bases / total_length) * 100
        overall_depth = weighted_depth / total_length
        
        # Store coverage stats
        printf "%.4f,%.4f\\n", overall_coverage, overall_depth > "coverage_stats.tmp"
    }'
    
    # Parse NTM-Profiler sub-species report - use sub-species as main species result
    species=\$(awk '
        /^Species report/ { in_species = 1; next }
        /^Resistance report/ { in_species = 0 }
        in_species && NF > 0 && !/^-/ && !/^Species,/ {
            # First data line after header
            split(\$0, fields, ",")
            if (fields[1] != "" && fields[1] != "Species") {
                print fields[1]
                exit
            }
        }
    ' ${ntm_sub_species})

    accession=\$(awk '
        /^Species report/ { in_species = 1; next }
        /^Resistance report/ { in_species = 0 }
        in_species && NF > 0 && !/^-/ && !/^Species,/ {
            split(\$0, fields, ",")
            if (fields[2] != "" && fields[2] != "Accession") {
                print fields[2]
                exit
            }
        }
    ' ${ntm_sub_species})

    ani=\$(awk '
        /^Species report/ { in_species = 1; next }
        /^Resistance report/ { in_species = 0 }
        in_species && NF > 0 && !/^-/ && !/^Species,/ {
            split(\$0, fields, ",")
            if (fields[3] != "" && fields[3] != "Ani") {
                print fields[3]
                exit
            }
        }
    ' ${ntm_sub_species})

    # Parse resistance genes from main NTM report
    resistance_genes=\$(awk '
        /^Resistance genes report/ { in_genes = 1; next }
        /^Resistance variants report/ { in_genes = 0 }
        in_genes && NF > 0 && !/^-/ {
            if (genes == "") genes = \$0
            else genes = genes ";" \$0
        }
        END { 
            if (genes == "") print "None"
            else print genes
        }
    ' ${ntm_report})

    # Parse resistance variants from main NTM report
    resistance_variants=\$(awk '
        /^Resistance variants report/ { in_variants = 1; next }
        /^Other variants report/ { in_variants = 0 }
        in_variants && NF > 0 && !/^-/ {
            if (variants == "") variants = \$0
            else variants = variants ";" \$0
        }
        END { 
            if (variants == "") print "None"
            else print variants
        }
    ' ${ntm_report})
    
    # Parse Kraken2 report - find highest percentage species (excluding unclassified and root)
    kraken_species=\$(awk '{
        # Skip unclassified, root, and very low percentages
        if (\$1 > 0.1 && \$4 == "S" && \$1 > max_percent) {
            max_percent = \$1
            max_species = ""
            # Extract species name from the end of the line
            for(i=6; i<=NF; i++) {
                if(i==6) max_species = \$i
                else max_species = max_species " " \$i
            }
            # Clean up leading/trailing spaces
            gsub(/^[ \\t]+|[ \\t]+\$/, "", max_species)
        }
    } END {
        if (max_species == "") print "Unknown"
        else print max_species
    }' ${kraken_report})
    
    kraken_percentage=\$(awk '{
        if (\$1 > 0.1 && \$4 == "S" && \$1 > max_percent) {
            max_percent = \$1
        }
    } END {
        if (max_percent == "") print "0.00"
        else print max_percent
    }' ${kraken_report})
    
    # Parse HMMER BLAST output to check for Erm(41) and get identity percentage
    erm41_result=\$(awk '{
        # Check if the subject contains Erm(41) - case insensitive
        if (tolower(\$2) ~ /erm\\(41\\)/) {
            # Found Erm(41), store the identity percentage (column 3)
            identity = \$3
            found_erm41 = 1
            exit
        }
    } END {
        if (found_erm41 == 1) print "Yes," identity
        else print "No,NA"
    }' ${hmmer_blast})
    
    # Split the result
    erm41_functional=\$(echo \$erm41_result | cut -d',' -f1)
    erm41_identity=\$(echo \$erm41_result | cut -d',' -f2)
    
    # Handle empty values
    if [ -z "\$species" ]; then species="Unknown"; fi
    if [ -z "\$accession" ]; then accession="Unknown"; fi
    if [ -z "\$ani" ]; then ani="0"; fi
    if [ -z "\$kraken_species" ]; then kraken_species="Unknown"; fi
    if [ -z "\$kraken_percentage" ]; then kraken_percentage="0.00"; fi
    if [ -z "\$erm41_functional" ]; then erm41_functional="No"; fi
    if [ -z "\$erm41_identity" ]; then erm41_identity="NA"; fi
    
    # Read calculated values using proper bash syntax
    IFS=',' read num_contigs total_length largest_contig smallest_contig avg_length n50 < assembly_stats.tmp
    IFS=',' read overall_coverage overall_depth < coverage_stats.tmp
    
    # Extract reference genome name from params
    ref_name=\$(basename ${params.reference} | sed 's/\\.[^.]*\$//')
    
    # Create CSV output - using sub-species results in main species columns
    cat > ${prefix}_sample_summary.csv << EOF
Sample_ID,Reference_Genome,NTM_Profiler_Species,Accession,ANI,Kraken2_Species,Kraken2_Percentage,Total_Contigs,Total_Length_bp,Largest_Contig_bp,Smallest_Contig_bp,Average_Contig_Length_bp,N50_bp,Overall_Coverage_percent,Overall_Mean_Depth,Resistance_Genes,Resistance_Variants,erm41_functional,erm41_identity_percent
${prefix},\$ref_name,\$species,\$accession,\$ani,\$kraken_species,\$kraken_percentage,\$num_contigs,\$total_length,\$largest_contig,\$smallest_contig,\$avg_length,\$n50,\$overall_coverage,\$overall_depth,\$resistance_genes,\$resistance_variants,\$erm41_functional,\$erm41_identity
EOF
    
    # Display summary to console
    echo ""
    echo "=== SUMMARY FOR ${prefix} ==="
    echo "NTM-Profiler (Sub-species): \$species | ANI: \$ani% | Accession: \$accession"
    echo "Kraken2: \$kraken_species (\$kraken_percentage%)"
    echo "Contigs: \$num_contigs | Coverage: \${overall_coverage}% | Avg Depth: \${overall_depth}x"
    echo "Assembly: \$total_length bp | Largest contig: \$largest_contig bp"
    echo "Resistance genes: \$resistance_genes"
    echo "Resistance variants: \$resistance_variants"
    echo "Erm(41) functional: \$erm41_functional (\$erm41_identity% identity)"
    echo "CSV saved: ${prefix}_sample_summary.csv"
    
    # Clean up temporary files
    rm -f assembly_stats.tmp coverage_stats.tmp
    """
}