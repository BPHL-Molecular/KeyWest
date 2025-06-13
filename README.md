# MABSC Identification and Resistance Prediction Pipeline

This modular [Nextflow](https://www.nextflow.io/) pipeline identifies *Mycobacterium abscessus complex* (MABSC), predicts drug resistance, performs identification, and assembles genomes from Illumina paired-end reads. It uses Docker/Singularity containers and generates summary results.

Features
- Contamination detection using `Kraken2`
- Adapter trimming using `fastp`
- Species & resistance prediction via `NTMprofiler`
- Assembly with `bwa`, `samtools`
- Resistance and virulence gene detection via `SnpEff`, `AMRfinder+`
- Outputs per-sample results and a final summary table

---
Raw FASTQ ? fastp ? NTM-Profiler (species ID/resistance)
                 ?
              Kraken2 (taxonomy)
                 ?
             BWA ? Samtools ? Pilon ? Prokka ? Trun-Label (resistance)
                                  ?
                            stats ? MultiQC
output/
+-- sample1/
¦   +-- fastp/           # QC + trimming
¦   +-- ntm_profiler/    # NTM species identification/Resistance prediction
¦   +-- kraken2/         # Taxonomic classification
¦   +-- bwa/            # Alignment
¦   +-- samtools/       # BAM files
¦   +-- pilon/          # Assembly
¦   +-- prokka/         # Annotation
¦   +-- trun_label/     # Resistance prediction
¦   +-- stats/          # Statistics
+-- multiqc/            # Report

Requirements

- Nextflow `>=22.10.0`
- Docker or Singularity
- Conda 
---

Input

Paired-end FASTQ files placed under `data/`:

```
data/
â”œâ”€â”€ sample1_1.fastq.gz
â”œâ”€â”€ sample1_2.fastq.gz
â”œâ”€â”€ sample2_1.fastq.gz
â”œâ”€â”€ sample2_2.fastq.gz
```

---

Configuration

Modify paths in `params.yaml` for database locations:


---

Run the Pipeline


---

Output

Each sample will have:
- QC reports
- Kraken2 contamination reports
- NTMprofiler results
- Assembly FASTA
- Gene annotations

The pipeline will generate a final summary file with key results across all samples.

---

Conda Environment 

```bash
conda env create -f env.yaml
```


---

## ðŸ“ž Contact

Developed by YH. For questions or contributions, please open an issue.
