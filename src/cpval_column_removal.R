#######################################
#
# cpval_column_removal.R
# Brett Ford
# Created 20180612
#
# This script removes columns with extraneous or repetitive information and
# renames Bayescan results with new naming convention
#
# Usage: Rscript cpval_column_removal.R
#
######################################

#Load required libraries
library(readr)

#List cpval files with data appended from Forester et al. 2018
root_path <- "/Users/brettford/Desktop/Northeastern/coding/forester_simulation_code/forester_sim_code/TTT_LotterhosWhitlockData/"
forester_simfiles_path <- paste0(root_path, "forester_simfiles/")
forester_results_path <- paste0(root_path, "forester_results/")

#Get cpval files
Filenames.cpval <- list.files(path=forester_results_path, pattern="Cpval")

#Get number of files; there should be 72
length(Filenames.cpval)

#Check to make sure the right files were listed
head(Filenames.cpval)

#Create empty dataframe to keep track of mismatches
i <- "1R_P90_1351142954_453_1_NumPops=90_NumInd=20Bayenv2LFMMpca.Cpval"
for(i in Filenames.cpval){
  
  #Read table was not working for me, so I switched to read.delim
  cpval_table <- read.table(paste0(forester_results_path, i), header=TRUE, sep=" ")
  cpval_table <- cpval_table[,c("SNPnames", "LA.All", "Corr.All", "FST.All", 
                                 "s_high", "ENVI_ID", "demog", "Set", "p.end",
                                 "SNPIncluded", "filename", "He.LS", "He_samp",
                                 "log.bf", "rho", "xtx", "cols", "pch1", "rda_veganv2.4.3_loading",
                                 "dbrda_veganv2.4.3_loading")]
  #Replace column names with new naming convention
  names(cpval_table)[names(cpval_table)=="log.bf"] <- "bayescan_bayescan2.1_log.bf"
  names(cpval_table)[names(cpval_table)=="rho"] <- "bayescan_bayescan2.1_rho"
  names(cpval_table)[names(cpval_table)=="xtx"] <- "bayescan_bayescan2.1_xtx"
  
  #Check to make sure column names look right
  names(cpval_table)
  
  #Write table to new file (replacing old files)
  write.table(cpval_table, file=paste0(forester_results_path, i), sep= " ", col.names = TRUE, row.names = FALSE)
}