for file in ./fastqs/*_R[12]_*.fastq.gz; do
    base=$(basename "$file")
    prefix=$(echo "$base" | sed -E 's/^([A-Za-z0-9]+)-.*_R[12]_.*\.fastq\.gz$/\1/')
    read_num=$(echo "$base" | sed -E 's/^.*_R([12])_.*\.fastq\.gz$/\1/')
    mv "$file" "./fastqs/${prefix}_${read_num}.fastq.gz"
done