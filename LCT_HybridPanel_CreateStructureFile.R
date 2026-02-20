#Convert software output to structure format
#Created 2/29/2024 by MEFL, updated 7/23/2024 by MEFL

#### How to use this script ####

# Set input and output file names as variables in the next code block, then
#  run the entire script to create a STRUCTURE-formatted file

# After creating you structure file, open file in a text editor and delete all rows
#  except for the samples you want to run structure on (keep your positive controls
#  and any samples that did not come out as 100% one species in your allele counting
#  analysis)



#### SET 2 VARIABLES ####

## 1. Input file name
#  - Export called genotypes from software in format "Detailed table results" and list name of csv here
#  - List file name here with FULL PATH to its location on your computer (pull from Box or from local computer)
#example: file_in <- "Box/GVL/LCT/X9/LCT_HybPanel_Plate1_20240226/HybPlate1_20240226_DetailedTableResults.csv"
file_in <- "in.csv"


## 2. Desired STRUCTURE file name
#  - List desired output STRUCTURE file name ending with .stru with FULL PATH
#example: file_out <- "Desktop/HybridPlate1_20240226.stru"
file_out <- "out.csv"



#### Import data ####

## Import detailed table results csv from genotyping software
input <- read.csv(file_in, #import csv file designated in code block above
                  skip=15) #import starting from row 16 actual data

## Remove NTC rows and blank assays
input <- input[!grepl("NTC",input$Type,fixed=T),]
input <- input[!grepl("Blank",input$Assay,fixed=T),]

## Convert missing data to N:N
input$Converted[input$Converted=="No Call"] <- "N:N"
input$Converted[input$Converted=="Invalid"] <- "N:N"

## Split genotype into alleles
input$allele1 <- substr(input$Converted,1,1)
input$allele2 <- substr(input$Converted,3,3)


##### Locus data ####

## Subset input data to desired loci
#  - the 95 loci in the LCT hybrid panel have been reduced to a set of 58 loci (14 LCT, 32 RBT, 12 YCT)
#    that passed the validation test - only these loci should be analyzed
#  - manual list of loci below for easy use of this R script, but alternatively can import text file list of loci
#loci <- read.table("LCT_Final58loci_2024jul23.txt")[,1] #import text file list of loci
loci <- c("LCT_06", "LCT_11", "LCT_14", "LCT_25", "LCT_29", "LCT_35", "LCT_46", "LCT_60", "LCT_61", "LCT_70", "LCT_75", "LCT_84",
          "LCT_87", "LCT_92", "RBT_01", "RBT_07", "RBT_09", "RBT_17", "RBT_21", "RBT_22", "RBT_27", "RBT_28", "RBT_30", "RBT_31",
          "RBT_32", "RBT_33", "RBT_40", "RBT_41", "RBT_43", "RBT_44", "RBT_45", "RBT_54", "RBT_56", "RBT_62", "RBT_65", "RBT_71",
          "RBT_74", "RBT_76", "RBT_80", "RBT_81", "RBT_82", "RBT_88", "RBT_89", "RBT_93", "RBT_94", "RBT_95", "YCT_02", "YCT_10",
          "YCT_12", "YCT_42", "YCT_47", "YCT_49", "YCT_52", "YCT_55", "YCT_66", "YCT_77", "YCT_78","YCT_79")
input <- input[input$Assay %in% loci,]



#### Create structure file ####

## Convert long table to wide genotype matrix
wide_data <- reshape(input[,c("Assay","Name","allele1","allele2")],idvar="Name",timevar="Assay",direction="wide")


## Define a function to convert nucleotide calls to numerical codes
nucleotide_to_code <- function(genotype) {
  ifelse(genotype == "A", 1,
         ifelse(genotype == "T", 2,
                ifelse(genotype == "C", 3,
                       ifelse(genotype == "G", 4,
                              ifelse(genotype == "N", -9, NA)))))  # Encode missing data
}

## Apply the conversion function to each element of the genotype data
wide_data_numeric <- wide_data
wide_data_numeric[,-1] <- apply(wide_data[,-1], MARGIN = 2, FUN = function(col) nucleotide_to_code(col))
#View(wide_data_numeric)

## Add location data (required column to run structure, but it won't use the info so just = 1 for all)
wide_data_numeric <- cbind("Name" = wide_data_numeric[, 1], "loc" = 1, wide_data_numeric[, -1])

## Sort rows alphabetically, just to stay organized
wide_data_numeric <- wide_data_numeric[order(wide_data_numeric$Name),]

## Export structure-formatted file (excluding header)
write.table(wide_data_numeric,file_out,row.names=F,quote=F,col.names=F)

## Check number of loci and number of samples to input as parameters in structure run
#nrow(wide_data_numeric) #number of samples
#ncol(wide_data_numeric)/2 - 1 #number of loci
