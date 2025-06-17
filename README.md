# KeyWest - MABSC Identification and Resistance Prediction Pipeline

This Nextflow pipeline identifies *Mycobacterium abscessus complex* (MABSC), predicts drug resistance, performs species identification, and assembles genomes from Illumina paired-end reads. It uses Docker/Singularity containers and generates comprehensive summary results with Erm(41) resistance gene detection.

## Requirements
- Nextflow 23.04.0+
- Singularity 
- Conda 

---

## Input

Paired-end FASTQ files placed under `fastqs/`:
```
sample1_1.fastq.gz
sample1_2.fastq.gz
sample2_1.fastq.gz
sample2_2.fastq.gz
```

**For lab samples format as `prefix-date-xxx-xx_S1_L001_R1_001.fastq.gz`:**
```bash
# Use resource/rename.sh to convert lab sample names
bash resource/rename.sh
```

---

## Setup

### Check Parameters
Verify `params.yaml` has all required paths:
```yaml
input_dir: "/path/to/fastqs"
output: "/path/to/output"
reference: "/path/to/reference.fasta"
erm41: "/path/to/erm41.fasta"
card_database: "/path/to/card_db"
kraken2_db: "/path/to/kraken2_db"
NTM_profiler_db: "/path/to/NTM-profiler"
```

### Check FASTQ Files
Ensure files follow naming convention:
```
sample1_1.fastq.gz & sample1_2.fastq.gz
sample2_1.fastq.gz & sample2_2.fastq.gz
```

### Build Conda Environment
```bash
conda env create -f resource/keywest.yml
conda activate keywest
```

---

## Run the Pipeline

### Local Execution
```bash
nextflow run main.nf -params-file params.yaml
```

### HiPerGator Usage
```bash
sbatch keywest.sh
```

---

## Output Structure

```
output/
+-- {sample}/
¦   +-- fastp/              # Quality control reports
¦   +-- ntm_profiler/       # Species identification
¦   +-- kraken2/            # Taxonomic classification
¦   +-- bwa/                # Read alignment
¦   +-- samtools/           # BAM processing
¦   +-- unicycler/          # Genome assembly
¦   +-- hmmer/              # Erm(41) resistance gene analysis
¦   +-- stats/              # Coverage statistics
¦   +-- summary/            # Comprehensive CSV report
¦       +-- {sample}_sample_summary.csv
+-- multiqc_report.html     # Aggregate quality report
```


## Pipeline Workflow

1. **fastp** - Quality control and trimming
2. **NTM-Profiler** - Mycobacterial species identification
3. **Kraken2** - Taxonomic classification
4. **BWA** - Read mapping to reference
5. **Samtools** - BAM file processing
6. **Unicycler** - Genome assembly
7. **HMMER** - Erm(41) resistance gene detection + CARD BLAST
8. **Stats** - Coverage and mapping statistics
9. **Summary** - Comprehensive CSV report generation
10. **MultiQC** - Aggregate quality reporting

---

## Database Requirements

### Required Databases
- **Reference genome**: Mycobacterial reference (FASTA)
- **Erm(41) sequences**: Query genes for resistance detection
- **CARD database**: Comprehensive Antibiotic Resistance Database
- **Kraken2 database**: Taxonomic classification
- **NTM-Profiler database**: Mycobacterial species identification. (The NTM-Profiler database for BPHL use by HiperGator has been modified by author.)

### Database Setup 
(For Non-HiPerGator Users)
```bash
# CARD Database
wget https://card.mcmaster.ca/latest/data
tar -xf data
makeblastdb -in card_protein_homolog_model.fasta -dbtype prot -out card_db

# Kraken2 Database  
wget https://genome-idx.s3.amazonaws.com/kraken/minikraken2_v2_8GB_201904.tgz
tar -xzf minikraken2_v2_8GB_201904.tgz

# NTM-Profiler
git clone https://github.com/jodyphelan/NTM-Profiler.git
```

---

## Troubleshooting

(For Non-HiPerGator Users)

Please install the following modules manually:

**HMMER**: https://github.com/EddyRivasLab/hmmer/tree/master
```bash
wget http://eddylab.org/software/hmmer/hmmer.tar.gz
tar zxf hmmer.tar.gz
cd hmmer-3.4
./configure --prefix /your/install/path
make && make install
```

**NCBI BLAST+**: https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html
```bash
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-*-x64-linux.tar.gz
tar -xzf ncbi-blast-*-x64-linux.tar.gz
```

**BWA-MEM2**: https://github.com/bwa-mem2/bwa-mem2
```bash
git clone https://github.com/bwa-mem2/bwa-mem2
cd bwa-mem2
make
```
