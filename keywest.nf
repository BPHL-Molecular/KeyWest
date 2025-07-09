#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Include module definitions
include { fastp } from './modules/fastp'
include { kraken2 } from './modules/kraken2'
include { bwa } from './modules/bwa'
include { samtools } from './modules/samtools'
include { unicycler } from './modules/unicycler'
include { ntm_profiler } from './modules/ntm_profiler'
include { hmmer } from './modules/hmmer'
include { stats } from './modules/stats'
include { summary } from './modules/summary'
include { multiqc } from './modules/multiqc'

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
if (!params.card_database) {
    error "Missing required parameter: --card_database"
}
if (!params.erm41) {
    error "Missing required parameter: --erm41"
}
if (!params.NTM_profiler_db) {
    error "Missing required parameter: --NTM_profiler_db"
}
if (!file(params.input_dir).exists()) {
    error "Input directory does not exist: ${params.input_dir}"
}
if (!file(params.reference).exists()) {
    error "Reference file does not exist: ${params.reference}"
}
if (!file(params.erm41).exists()) {
    error "ERM41 query file does not exist: ${params.erm41}"
}
if (!file(params.NTM_profiler_db).exists()) {
    error "NTM-Profiler database directory does not exist: ${params.NTM_profiler_db}"
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
    
    // Step 2: Run Kraken2 for taxonomic classification (on trimmed reads)
    kraken2_out = kraken2(trimmed_fastqs)
    kraken2_reports = kraken2_out.report
    
    // Step 3: Run BWA on trimmed reads
    sam_files = bwa(trimmed_fastqs)
    
    // Step 4: Run Samtools (convert SAM to BAM)
    samtools_out = samtools(sam_files)
    bam_files = samtools_out.bam
    
    // Step 5: Run unicycler for assembly
    unicycler_out = unicycler(trimmed_fastqs)
    unicycler_assembly = unicycler_out.assembly
    
    // Step 6: Run NTM-Profiler on assembly (MOVED AFTER UNICYCLER)
    ntm_profiler_out = ntm_profiler(unicycler_assembly)
    ntm_reports = ntm_profiler_out.report
    ntm_sub_species = ntm_profiler_out.sub_species
    ntm_vcfs = ntm_profiler_out.vcf
       
    // Step 7: Run HMMER analysis on unicycler assembly files
    hmmer_out = hmmer(unicycler_assembly)
    hmmer_summaries = hmmer_out.summary
    
    // Step 8: Generate Stats from BAM files
    stats_out = stats(bam_files)
    
    // Step 9: Generate Summary for each sample
    // Transform channels to [prefix, file] format for joining
    unicycler_ch = unicycler_assembly.map { assembly, prefix -> [prefix, assembly] }
    stats_ch = stats_out.coverage.map { coverage_file, prefix -> [prefix, coverage_file] }
    ntm_ch = ntm_reports.map { report, prefix -> [prefix, report] }
    ntm_sub_ch = ntm_sub_species.map { sub_report, prefix -> [prefix, sub_report] } 
    kraken2_ch = kraken2_reports.map { report, prefix -> [prefix, report] }
    hmmer_ch = hmmer_summaries.map { blast_file, prefix -> [prefix, blast_file] }
    
    // Join all channels by prefix - this ensures ALL processes finish before summary
    summary_input = unicycler_ch
        .join(stats_ch)
        .join(ntm_ch)
        .join(ntm_sub_ch)
        .join(kraken2_ch)
        .join(hmmer_ch)
        .map { prefix, assembly, coverage_file, ntm_report, ntm_sub_species, kraken_report, hmmer_blast ->
     	    tuple(assembly, coverage_file, ntm_report, ntm_sub_species, kraken_report, hmmer_blast, prefix)
	}    
    summary_out = summary(summary_input)
    summary_reports = summary_out.summary
    
    // Step 10: Run MultiQC to aggregate all reports
    fastp_files = fastp_reports
    kraken2_files = kraken2_reports.map { file, prefix -> file }  
    
    // Use the stats channel from your updated stats process
    stats_files = stats_out.stats.map { file, prefix -> file }
    
    // Combine all report files for MultiQC
    multiqc_inputs = Channel.empty()
        .mix(fastp_files)
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