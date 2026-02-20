## LCT hybrid panel - data manipulation and hybrid analysis using species-specific allele counting
#Create 2/29/2024 by MEFL, modified 7/23/2024 by MEFL

#### How to use this script ####

# Set input and output file names as variables in the next code block, then
#  run the entire script to calculate admixture proportions based on counting
#  species-specific alleles at loci for LCT, RBT, and YCT



#### SET 2 VARIABLES ####

## 1. Input file name
#  - Export called genotypes from software in format "Detailed table results" and list name of csv here
#  - List file name here with FULL PATH to its location on your computer (pull from Box or from local computer)
#example: file_in <- "Box/GVL/LCT/X9/LCT_HybPanel_Plate1_20240226/HybPlate1_20240226_DetailedTableResults.csv"
file_in <- "LCT_HybridPanel_Plate6_20240319_DetailedTableResults.csv"


## 2. Output file name
#  - List desired output file name with FULL PATH to GVL/LCT/Results folder on Box (or save to your computer, then load to Box)
#example: file_out <- "Box/GVL/LCT/Results/HumboldtBasin_2024/AlleleCountingResults_ByPlate/LCTHybridPanel_Plate1_20240226.csv"
file_out <- "NewResults.csv"



#### Import & re-format data ####

##### Genotype data ####

## Import detailed table results csv from genotyping software
input <- read.csv(file_in, #import csv file designated in code block above
                  skip=15) #import starting from row 16 where actual data is

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


## List of species-specific alleles at each locus
#  - manual allele info below for easy use of this R script, but alternatively can import csv allele info
#  - note that manual list only has alleles for 58 loci, whereas csv has alleles for all 95 loci
#allele <- read.csv("SpeciesAlleleList.csv") #import species-specific allele info
allele <- data.frame(Assay=loci,
                    LCT_allele=c("G", "A", "A", "T", "A", "T", "T", "T", "A", "T", "G", "T", "G", "T", "T", 
                                 "T", "A", "G", "C", "T", "A", "G", "C", "A", "C", "C", "A", "T", "C", "T", 
                                "G", "G", "A", "T", "G", "G", "C", "C", "G", "C", "A", "A", "T", "C", "G", 
                                "G", "G", "G", "C", "G", "C", "C", "T", "G", "T", "G", "C","A"),
                    RBT_allele=c("A", "G", "G", "C", "G", "G", "C", "C", "G", "C", "T", "C", "A", "G", "A", 
                                "G", "G", "A", "G", "G", "G", "A", "G", "G", "T", "T", "G", "G", "T", "G", 
                                "A", "A", "C", "G", "A", "A", "T", "T", "A", "T", "T", "G", "G", "T", "C", 
                                "T", "G", "G", "C", "G", "C", "C", "T", "G", "T", "G", "C","A"),
                    YCT_allele=c("A", "G", "G", "C", "G", "N", "C", "C", "G", "C", "T", "C", "A", "G", "T", 
                                "T", "A", "G", "C", "T", "A", "G", "C", "A", "C", "C", "A", "T", "C", "T", 
                                "G", "G", "A", "T", "G", "G", "C", "C", "G", "C", "A", "A", "T", "C", "G", 
                                "G", "A", "T", "T", "T", "T", "G", "A", "T", "C", "A", "A","G"))





#### Count species-specific alleles ####

## Add columns to input for each species
input$LCTcount <- NA
input$RBTcount <- NA
input$YCTcount <- NA

## Loop through input to count species-specific alleles
# 2 = both alleles match the species-specific allele
# 1 = one allele matches
# 0 = neither matches 
# NA = missing data

for (i in 1:nrow(input)){
  if (input$allele1[i]=="N"){} #skip rows with missing data
  else{
    if (grepl("LCT",input$Assay[i])){ #LCT loci
      input$LCTcount[i] <- (input$allele1[i] == allele$LCT_allele[allele$Assay==input$Assay[i]]) +
        (input$allele2[i] == allele$LCT_allele[allele$Assay==input$Assay[i]])
    }
    if (grepl("RBT",input$Assay[i])){ #RBT loci
      input$RBTcount[i] <- (input$allele1[i] == allele$RBT_allele[allele$Assay==input$Assay[i]]) +
        (input$allele2[i] == allele$RBT_allele[allele$Assay==input$Assay[i]])
    }
    if (grepl("YCT",input$Assay[i])){ #YCT loci
      input$YCTcount[i] <- (input$allele1[i] == allele$YCT_allele[allele$Assay==input$Assay[i]]) +
        (input$allele2[i] == allele$YCT_allele[allele$Assay==input$Assay[i]])
    }
  }
}


## Create table of samples to sum species-specific allele counts
results <- data.frame(sample=sort(unique(input$Name)),
                   LCTcount=NA,RBTcount=NA,YCTcount=NA, #sum number of species-specific alleles
                   LCTnonmiss=NA,RBTnonmiss=NA,YCTnonmiss=NA, #sum number of non-missing alleles genotyped
                   PropMissing=NA, #calc % missing data
                   LCTprop_AlleleCount=NA,RBTprop_AlleleCount=NA,YCTprop_AlleleCount=NA) #calc admixture prop for each species

for (i in 1:nrow(results)){
  #sum number of species-specific alleles
  results$LCTcount[i] <- sum(input$LCTcount[input$Name==results$sample[i]],na.rm=T)
  results$RBTcount[i] <- sum(input$RBTcount[input$Name==results$sample[i]],na.rm=T)
  results$YCTcount[i] <- sum(input$YCTcount[input$Name==results$sample[i]],na.rm=T)
  
  #sum number of non-missing alleles genotyped
  results$LCTnonmiss[i] <- sum(!is.na(input$LCTcount[input$Name==results$sample[i]])) * 2
  results$RBTnonmiss[i] <- sum(!is.na(input$RBTcount[input$Name==results$sample[i]])) * 2
  results$YCTnonmiss[i] <- sum(!is.na(input$YCTcount[input$Name==results$sample[i]])) * 2
}

for (i in 1:nrow(results)){
  #calc missing data %
  results$PropMissing[i] <- 1 - (rowSums(results[i,c("LCTnonmiss","RBTnonmiss","YCTnonmiss")]) / max(rowSums(results[,c("LCTnonmiss","RBTnonmiss","YCTnonmiss")])))
  
  #calc admixture prop for each species (# species-specific alleles / # non-missing alleles genotyped)
  # EXCLUDING samples with >50% missing data
  if (rowSums(results[i,c("LCTnonmiss","RBTnonmiss","YCTnonmiss")]) > max(rowSums(results[,c("LCTnonmiss","RBTnonmiss","YCTnonmiss")]))/2) {
    results$LCTprop_AlleleCount[i] <- results$LCTcount[i]/results$LCTnonmiss[i]
    results$RBTprop_AlleleCount[i] <- results$RBTcount[i]/results$RBTnonmiss[i]
    results$YCTprop_AlleleCount[i] <- results$YCTcount[i]/results$YCTnonmiss[i]
  }
}
#View(results)



#### Export results ####
write.csv(results[,c("sample","PropMissing","LCTprop_AlleleCount","RBTprop_AlleleCount","YCTprop_AlleleCount")],
          file_out, #name output based on variable name set at beginning of this script
          row.names=F,quote=F)
#SAVE copy of this results csv to Box folder: GVL/LCT/Results/AlleleCountingResults_ByPlate/
#COPY results into GVL/LCT/Results/Results_CompiledAcrossRuns.xlsx, excluding positive control samples
