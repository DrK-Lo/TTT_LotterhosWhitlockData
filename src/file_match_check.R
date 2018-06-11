#######################################
#
# file_match_check.R
# Brett Ford
# Created 20180611
#
# This script checks to make sure allele frequencies between variations 
# of the original LFMM files are equal. It only checks the first 10 columns
# between files
#
# Usage: Rscript file_match_check.R
#
######################################

#List files from original dryad files
#This requires that the dryad repository be downloaded
#Not present on the github repository, due to large memory requirements
dryad_path <- "/Users/brettford/Desktop/Northeastern/coding/forester_simulation_code/dryad"

#Get lfmm files
Filenames.lfmm <- list.files(path=paste0(dryad_path,"/SimFilesLFMM/"), pattern="NumPops=90.*6.lfmm|NumPops=90.*20.lfmm")

#Get number of files; there should be 72
length(Filenames.lfmm)
#Check to make sure the right files were listed
head(Filenames.lfmm)
logical_results <- data.frame(output="Output")
for(i in Filenames.lfmm){
  filename_df <- data.frame(output=i)
  logical_results <- rbind(logical_results, filename_df)
  print(logical_results)
  #Read in original file
  original_file <- read.table(paste0(dryad_path,"/SimFilesLFMM/",i), header = FALSE, sep=" ")
  #Change 2s to 1s which is how randomized files are formatted
  original_file[original_file==2] <- 1
  
  #Read in randomized file
  randomized_file <- read.table(paste0("../forester_simfiles/", i), header = FALSE, sep=" ")
  if(grepl("20.lfmm", i)){
    logical_table <- data.frame(matrix(NA, nrow = (nrow(randomized_file)/20), ncol = 10))
    for(j in 1:10){
    r1=1
      for(k in 1:(nrow(randomized_file)/20)){
      original_freq <- sum(original_file[r1:(20*k),j]==1)/20
      randomized_freq <- sum(randomized_file[r1:(20*k),j]==1)/20
        if(original_freq==randomized_freq){
        logical_table[k,j]<-"TRUE"
        }else{
            logical_table[k,j]<-"FALSE"
          }
      r1=r1+20
      }
    }
    if(sum(grepl("FALSE", logical_table))>0){
      mismatch_statement <- data.frame(output=paste("Mismatch identified in ",i))
      logical_results <- rbind(logical_results, mismatch_statement)
    }else{
      nomismatch_statement <- data.frame(output=paste("No mismatch identified in ",i))
      logical_results <- rbind(logical_results, nomismatch_statement)
    }
  }
  if(grepl("6.lfmm", i)){
    logical_table <- data.frame(matrix(NA, nrow = (nrow(randomized_file)/6), ncol = 10))
    for(j in 1:10){
    r1=1
      for(k in 1:(nrow(randomized_file)/6)){
        original_freq <- sum(original_file[r1:(6*k),j]==1)/6
        randomized_freq <- sum(randomized_file[r1:(6*k),j]==1)/6
        if(original_freq==randomized_freq){
          logical_table[k,j]<-"TRUE"
        }else{
          logical_table[k,j]<-"FALSE"
        }
        r1=r1+6
      }
    }
    if(sum(grepl("FALSE", logical_table))>0){
      mismatch_statement <- data.frame(output=paste("Mismatch identified in ",i))
      logical_results <- rbind(logical_results, mismatch_statement)  
    }else{
      nomismatch_statement <- data.frame(output=paste("No mismatch identified in ",i))
      logical_results <- rbind(logical_results, nomismatch_statement)
      }
  }
  
}

write.table(logical_results, file = paste0("../forester_results/","logical_results.txt"), sep = "\t", row.names=FALSE,col.names = FALSE)
