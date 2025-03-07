# check integrity of fastq files
# katherine carbeck
# 5 march 2025

#!/bin/bash

#define output file
OUTPUT_FILE="gzip_check_results.txt"

# clear previous results
> "$OUTPUT_FILE"

# loop through all compressed fastq files
for file in *.fastq.gz *.fq.gz; do
    if [[ -f "$file" ]]; then 
        echo "Checking: $file" | tee -a "$OUTPUT_FILE"
        gunzip -t "$file" 2>>"$OUTPUT_FILE"
        if [[ $? -eq 0 ]]; then
            echo "$file: OK" | tee -a "$OUTPUT_FILE"
        else
            echo "$file: CORRUPTED or TRUNCATED" | tee -a "$OUTPUT_FILE"
        fi
        echo "-----------------" >> "$OUTPUT_FILE"
    fi
done

echo "check complete! results logged in $OUTPUT_FILE"

# LC_T3_171252-3 - corrupted R2 AT SOURCE (google drive)
# 171252-3_R2_001.fastq.gz
