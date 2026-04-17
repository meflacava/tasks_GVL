#!/bin/bash
#SBATCH --job-name=countReads
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH -p bmh

### Count the # lines in fastq files to get # sequenced reads for a specific list of samples

# --- Configuration ---
SAMPLE_LIST="" #text file containing list of samples to target (e.g., Sample1 to target file Sample1_R1.fastq.gz)
SEARCH_ROOT="" #path to fastq files
OUTPUT_FILE="readCounts_R1.txt"

# Clear previous output if it exists
> "$OUTPUT_FILE"

echo "Starting read count process at $(date)"
echo "Searching in: $SEARCH_ROOT"
echo "Reading samples from: $SAMPLE_LIST"

# Check if sample list exists
if [ ! -f "$SAMPLE_LIST" ]; then
    echo "ERROR: Sample list file '$SAMPLE_LIST' not found."
    exit 1
fi

# Loop through each sample ID in the text file
while IFS= read -r sample_id || [ -n "$sample_id" ]; do
    # Skip empty lines
    [ -z "$sample_id" ] && continue

    echo "Processing: $sample_id"

    # Construct the filename pattern we are looking for
    # Files look like: Ht01-62_1995_F08_R1.fastq.gz
    target_pattern="${sample_id}_R1.fastq.gz"

    # Find the file recursively
    # -type f: only files
    # -name: exact match for the constructed pattern
    # -print -quit: stop searching after the first match (speed optimization)
    file_path=$(find "$SEARCH_ROOT" -type f -name "$target_pattern" -print -quit 2>/dev/null)

    if [ -z "$file_path" ]; then
        echo "  WARNING: File not found for $sample_id ($target_pattern). Skipping."
        echo -e "${sample_id}\tNOT_FOUND" >> "$OUTPUT_FILE"
        continue
    fi

    # Count reads: zcat -> wc -l -> awk divide by 4
    read_count=$(zcat "$file_path" | wc -l | awk '{print int($1/4)}')

    # Append result to output file
    echo -e "${sample_id}\t${read_count}" >> "$OUTPUT_FILE"

done < "$SAMPLE_LIST"

echo "Process completed at $(date)."
echo "Results saved to: $OUTPUT_FILE"
