process multiqc {
    publishDir "${params.output}/multiqc", mode: 'copy'

    
    input:
    path('*')
    
    output:
    path("*")  // All files created
    
    script:
    """
    multiqc . \\
        --title "NTM Pipeline Results" \\
        --filename multiqc_report.html \\
        --force
    
    # Move all data files to current directory (flatten structure)
    if [ -d multiqc_report_data ]; then
        mv multiqc_report_data/* .
        rmdir multiqc_report_data
    fi
    """
}