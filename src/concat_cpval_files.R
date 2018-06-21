######################################
#
# concat_cpval_files.R
# Brett Ford
# Created 20180314
#
# This script takes the original Lotterhos-Whitlock Cpval files
# and appends the data for the dbrda and rda analyses
#
# This assumes that ordination files have SNPnames
#
# Usage: Rscript concat_cpval_files.R
#
######################################

#Set directories- You should just have to change these if the repository was cloned
root_path <- root_path <- "/Users/brettford/Desktop/Northeastern/coding/forester_simulation_code/forester_sim_code/TTT_LotterhosWhitlockData/"
dryad_path <- "/Users/brettford/Desktop/Northeastern/coding/forester_simulation_code/dryad"

#Set subdirectories
forester_results_path <- paste0(root_path, "forester_results/")
ord_path <- paste0(forester_results_path, "ordination_results/")

#Load libraries
library(plyr)

#Unlist ordination files output from running Brenna's script
Filenames.ORD <- list.files(path=ord_path, pattern="ORD_locus_stats") 

for (i in 1:length(Filenames.ORD)) {               
  sim1 <- read.table(paste0(ord_path, Filenames.ORD[[i]]), header = T, sep = " ")
  sim1 <- sim1[,2:(ncol(sim1)-1)]

  #Column names should include: method_software&version_stat_logical(if applicable)
  colnames(sim1) <- c("SNPnames", "rda_veganv2.4-3_loading", "rda_veganv2.4-3_empP", "rda_veganv2.4-3_empQ", "rda_veganv2.4-3_empP_log", "rda_veganv2.4-3_SD2.5_log", 
                               "rda_veganv2.4-3_SD3.0_log", "dbrda_veganv2.4-3_loading", "dbrda_veganv2.4-3_empP", "dbrda_veganv2.4-3_empQ", "dbrda_veganv2.4-3_empP_log",
                               "dbrda_veganv2.4-3_SD2.5_log", "dbrda_veganv2.4-3_SD3.0_log")
  
  #Get basename of Cpval file from ordination files
  cpval_file <- gsub(Filenames.ORD[[i]], pattern="ORD_locus_stats_", replacement="")
  cpval_file <- gsub(cpval_file, pattern=".txt$", replacement="")
  
  #Read in Lotterhos Whitlock Cpval file
  Cpvaltable <- read.table(paste0(dryad_path, "/SummaryFiles/", cpval_file, "Bayenv2LFMMpca.Cpval"), header = T, sep = " ")

  # Create a lookup table from ordination file
  lookup <- unique(sim1)
  # Use join to append the lookup values to the Cpvaltable, based on SNPnames
  plyr1 <- join(Cpvaltable, lookup, by = "SNPnames")

  #write table to file (should write a total of 72 files)
  write.table(plyr1, file = paste0(forester_results_path, cpval_file, "_Bayenv2LFMMpca.Cpval"), sep = " ", row.names = F)
}


