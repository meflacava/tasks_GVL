#!/bin/bash

### Determine the number of properly paired reads successfully aligned to the genome, and add it to a file containing the number of reads sequenced for each sample (run countReads_R1.sh to generate readCounts_R1.txt). Also calculate the proportion of aligned reads


module load samtools

INPUT="readCounts_R1.txt"
BAM_DIR="" #path to BAM files
OUT="propReadsAligned.txt"

echo "id n.seq n.align prop.align" > "$OUT"

while read sample total; do
    #double the number of R1 reads to get the total number of reads (R1+R2)
    total=$((total*2))
    
    # Find the BAM file starting with the sample ID
    # We look for files matching: ${sample}*.bam
    bam_file=$(find "$BAM_DIR" -maxdepth 1 -type f -name "${sample}*.bam" -print -quit 2>/dev/null)
    
    if [ -z "$bam_file" ]; then
        align="NA"
        prop="NA"
    else
        # Count aligned reads (-f 2 = properly paired reads; -F 1284 = filter out unmapped reads (4), secondary reads (256), and duplicate reads (1024))
        align=$(samtools view -c -f 2 -F 1284 "$bam_file" 2>/dev/null)
        
        if [ -z "$align" ] || [ "$align" -eq 0 ]; then
            align="NA"
            prop="NA"
        else
            prop=$(awk -v a="$align" -v t="$total" 'BEGIN {printf "%.6f", a/t}')
        fi
    fi
    
    echo -e "$sample\t$total\t$align\t$prop" >> "$OUT"
done < "$INPUT"
