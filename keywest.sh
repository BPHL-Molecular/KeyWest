#!/usr/bin/bash
#SBATCH --account=bphl-umbrella
#SBATCH --qos=bphl-umbrella-b
#SBATCH --job-name=keywest
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=300gb
#SBATCH --time=48:00:00
#SBATCH --output=keywest.%j.out
#SBATCH --error=keywest.err
#SBATCH --mail-user=yi.huang@flhealth.gov
#SBATCH --mail-type=FAIL,END

module load singularity
module load nextflow/24.10.5
module load bwa-mem2
module load ncbi_blast
module load hmmer

APPTAINER_CACHEDIR=./
export APPTAINER_CACHEDIR
NXF_SINGULARITY_CACHEDIR=./
export NXF_SINGULARITY_CACHEDIR


nextflow run keywest.nf -params-file params.yaml

# Check if pipeline completed successfully
if [ $? -eq 0 ]; then
    echo "Pipeline completed successfully. Combining summary files..."
    
    # Find and combine summary files
    SUMMARY_FILES=($(find output -path "*/summary/*_sample_summary.csv" 2>/dev/null))
    
    if [ ${#SUMMARY_FILES[@]} -gt 0 ]; then
        echo "Found ${#SUMMARY_FILES[@]} summary files to combine"
        head -1 "${SUMMARY_FILES[0]}" > summary.csv
        find output -path "*/summary/*_sample_summary.csv" -exec tail -n +2 {} \; >> summary.csv
        echo "Combined summary saved to: summary.csv"
        echo "Total samples: $(($(wc -l < summary.csv) - 1))"
    else
        echo "Warning: No summary files found!"
    fi
    
    # Move files to output directory and rename with timestamp
    echo "Moving files and adding timestamp..."
    mv summary.csv ./output
    mv *.out ./output
    mv *.err ./output
    
    # Rename output directory with timestamp
    TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
    mv output output_${TIMESTAMP}
    echo "Final output directory: output_${TIMESTAMP}"
    
else
    echo "Pipeline failed. Check error logs."
    exit 1
fi

echo "Job completed at $(date)"
