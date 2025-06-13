process pilon {
    publishDir "${params.output}/${prefix}/pilon", mode: 'copy'
    memory '32 GB'  // Increase memory allocation
    
    input:
    tuple path(bam_file), val(prefix)
    
    output:
    tuple path("${prefix}_pilon.fasta"), val(prefix)
    
    script:
    def reference = params.reference
    """
    # Load required modules
    module load pilon
    module load java
    
    # Copy and properly decompress reference file
    cp ${reference} reference_original
    
    # Force decompression
    if [[ ${reference} == *.gz ]]; then
        echo "Decompressing gzipped reference..."
        gunzip -c reference_original > reference.fna
    else
        echo "Reference not gzipped, copying..."
        cp reference_original reference.fna
    fi
    
    # Run Pilon with explicit Java memory settings
    java -Xmx28G -jar \$(which pilon | xargs dirname)/../share/pilon*/pilon*.jar \\
        --genome reference.fna \\
        --frags ${bam_file} \\
        --output ${prefix}_pilon \\
        --outdir . \\
        --changes \\
        --vcf
    """
}