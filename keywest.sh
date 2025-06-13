#!/usr/bin/bash
#SBATCH --account=bphl-umbrella
#SBATCH --qos=bphl-umbrella-b
#SBATCH --job-name=keywest
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=300gb
#SBATCH --time=48:00:00
#SBATCH --output=keywest.%j.out
#SBATCH --error=keywest.err
#SBATCH --mail-user=<EMAIL>
#SBATCH --mail-type=FAIL,END

module load singularity
module load nextflow
APPTAINER_CACHEDIR=./
export APPTAINER_CACHEDIR
module load bwa-mem2
module load pilon

nextflow run main.nf -params-file params.yaml