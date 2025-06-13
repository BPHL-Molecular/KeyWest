#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Include module definitions
include { fastp } from './modules/fastp'
include { ntm_profiler } from './modules/ntm_profiler'
include { kraken2 } from './modules/kraken2'
include { bwa } from './modules/bwa'
include { samtools } from './modules/samtools'
include { pilon } from './modules/pilon'
include { prokka } from './modules/prokka'
//include { trun_label } from './modules/trun_label'
include { stats } from './modules/stats'
include { multiqc } from './modules/multiqc'
//include { summary } from './modules/summary'

// Parameter validation
if (!params.input_dir) {
    error "Missing required parameter: --input_dir"
}
if (!params.reference) {
    error "Missing required parameter: --reference"
}
if (!params.output) {
    error "Missing required parameter: --output"
}
if (!file(params.input_dir).exists()) {
    error "Input directory does not exist: ${params.input_dir}"
}
if (!file(params.reference).exists()) {
    error "Reference file does not exist: ${params.reference}"
}

workflow {
    // Create paired FASTQ channel
    Channel
        .fromFilePairs("${params.input_dir}/*_{1,2}.fastq.gz", size: 2)
        .ifEmpty { error "No FASTQ pairs found in ${params.input_dir}" }
        .map { sample_id, reads -> 
            tuple(reads[0], reads[1], sample_id)
        }
        .set { paired_fastqs }
    
    // Pipeline steps
    // Step 1: Run fastp for QC and trimming
    fastp_out = fastp(paired_fastqs)
    trimmed_fastqs = fastp_out.reads
    fastp_reports = fastp_out.reports
    
    // Step 2: Run NTM-Profiler on trimmed reads
    ntm_profiler_out = ntm_profiler(trimmed_fastqs)
    ntm_reports = ntm_profiler_out.report
    
    // Step 3: Run Kraken2 for taxonomic classification (on trimmed reads)
    kraken2_out = kraken2(trimmed_fastqs)
    kraken2_reports = kraken2_out.report
    
    // Step 4: Run BWA on trimmed reads
    sam_files = bwa(trimmed_fastqs)
    
    // Step 5: Run Samtools (convert SAM to BAM)
    samtools_out = samtools(sam_files)
    bam_files = samtools_out.bam
    
    // Step 6: Run Pilon for reference-guided assembly using BAM files
    pilon_assembly = pilon(bam_files)
    
    // Step 7: Run Prokka on Pilon assembly
    prokka(pilon_assembly)
    
    // Step 8: Extract .faa files from Prokka output for trun-label
    //prokka_faa = prokka_out.prokka_results
        //.map { files, prefix -> 
            //def faa_file = files.find { it.toString().endsWith('.faa') }
            //return tuple(faa_file, prefix)
        //}
    
    // Step 9: Run trun-label on Prokka .faa files
    //trun_label_out = trun_label(prokka_faa)
    
    
    // Step 11: Generate Stats from BAM files
    stats_out = stats(bam_files)
    
// Step 12: Run MultiQC to aggregate all reports
    // Extract files from channels and collect them
    fastp_files = fastp_reports
    ntm_files = ntm_reports.map { file, prefix -> file }
    kraken2_files = kraken2_reports.map { file, prefix -> file }  
    stats_files = stats_out
    
    // Combine all report files for MultiQC
    multiqc_inputs = Channel.empty()
        .mix(fastp_files)
        .mix(ntm_files)
        .mix(kraken2_files)
        .mix(stats_files)
        .collect()
    
    multiqc(multiqc_inputs)
}

// Workflow completion message
workflow.onComplete {
    log.info """
    Pipeline completed at: ${workflow.complete}
    Execution status: ${workflow.success ? 'OK' : 'failed'}
    Execution duration: ${workflow.duration}
    Output directory: ${params.output}
    """
}