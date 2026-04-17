# tasks_GVL

This repo contains various scripts and tools for lab and data analysis for the GVL


## Scripts:

seq_stats scripts: (run in order listed below)
  - countReads_R1.sh = count the number of sequenced reads for a specific list of samples
  - readStats.py = print the total number of sequenced reads across all samples, and the average, min and max number of reads per sample
  - addPropAligned.sh = calculate the proportion of sequenced reads that were properly paired and aligne to genome
  - propStats.py = print the average, min and max proportion of aligned reads per sample


LCT_HybridPanel_AlleleCounting.R = data manipulation and hybrid analysis using a species-specific allele counting method for the Lahontan cutthroat trout; input = Standard Biotools Genotyping Software SNP genotypes in the format of detailed table results, output = csv with hybrid proportions

LCT_HybridPanel_CreateStructureFile.R = reformat SNP genotypes from detailed table results output by Standard Biotools Genotyping Software to format needed to run STRUCTURE analysis

## Other tools:

Template_PlateMapConversion.xlsx = convert your long-form sample list to 8x12 plate format
