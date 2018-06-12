#######################################
#
# cpval_column_check.R
# Brett Ford
# Created 20180612
#
# This script checks to see if columns in cpval files are equivalent and can
# therefore be concatenated into one
#
# Usage: Rscript cpval_column_check.R
#
######################################

#Load required libraries
library(readr)

#List files from original dryad files
#This requires that the dryad repository be downloaded
#Not present on the github repository, due to large memory requirements
dryad_path <- "/Users/brettford/Desktop/Northeastern/coding/forester_simulation_code/dryad"

#Get cpval files
Filenames.cpval <- list.files(path=paste0(dryad_path,"/SummaryFiles/"), pattern="Cpval")

#Get number of files; there should be 240
length(Filenames.cpval)

#Check to make sure the right files were listed
head(Filenames.cpval)

#Create empty dataframe to keep track of mismatches
column_check_results <- data.frame(output="Output")

#Suppress warnings from read.delim before starting. The function
#has trouble distinguishing between variable types in columns
suppressMessages(read.delim)

for(i in Filenames.cpval){
  
  #Read table was not working for me, so I switched to read.delim
  cpval_table <- read_delim(paste0(dryad_path, "/SummaryFiles/", i), " ", escape_double = FALSE, trim_ws = TRUE)
  
  #Create new, empty columns to keep track of mismatches between columns
  cpval_table$SNPIncluded_UseSNP <- "NA"
  cpval_table$SNP.LA2_LA.ALL <- "NA"
  
  for(j in 1:nrow(cpval_table)){
    if(cpval_table$SNPIncluded[j]==cpval_table$UseSNP[j]){
      cpval_table$SNPIncluded_UseSNP[j] <- "TRUE"
    }else{
      cpval_table$SNPIncluded_UseSNP[j] <- "FALSE"
    }
  }
  for(k in 1:nrow(cpval_table)){
    if(cpval_table$SNP.LA2[k]==cpval_table$LA.All[k]){
      cpval_table$SNP.LA2_LA.ALL[k] <- "TRUE"
    }else{
      cpval_table$SNP.LA2_LA.ALL[k] <- "FALSE"
    }
  }
  if(sum(grepl("FALSE", cpval_table$SNPIncluded_UseSNP))>0){
    mismatch_statement <- data.frame(output=paste("Mismatch identified between SNPIncluded and UseSNP in",i))
    column_check_results <- rbind(column_check_results, mismatch_statement)  
  }else{
    nomismatch_statement <- data.frame(output=paste("No mismatch identified between SNPIncluded and UseSNP in",i))
    column_check_results <- rbind(column_check_results, nomismatch_statement)
  }
  if(sum(grepl("FALSE", cpval_table$SNP.LA2_LA.ALL))>0){
    mismatch_statement <- data.frame(output=paste("Mismatch identified between SNP.LA2 and LA.ALL in",i))
    column_check_results <- rbind(column_check_results, mismatch_statement)  
  }else{
    nomismatch_statement <- data.frame(output=paste("No mismatch identified between SNP.LA2 and LA.ALL in",i))
    column_check_results <- rbind(column_check_results, nomismatch_statement)
  }
}

