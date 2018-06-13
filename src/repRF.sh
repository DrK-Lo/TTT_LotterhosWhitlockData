#!/bin/bash
set -e
set -u
set -o pipefail

######################################################################################
#
# repRF.sh
# Brett Ford
# Created 20180306
#
# This script runs multiple random forest R scripts across mutiple nodes on a cluster
#
# Usage: bash repRF.sh
#
#######################################################################################

mypath="/Users/brettford/Desktop/Northeastern/coding/forester_simulation_code/forester_sim_code/TTT_LotterhosWhitlockData/src"
resultspath="/Users/brettford/Desktop/Northeastern/coding/forester_simulation_code/forester_sim_code/TTT_LotterhosWhitlockData/forester_results"
cd $mypath

#Run across 60 nodes
start=1
finish=60

#Code to run across specific simulations
#declare -a redo=(1 3 5 7 9 11 13 23 47 55 65 67 71)

##############
#### run R script
#############
echo "Running R scripts"
for i in $(seq $start $finish)
#for i in "${redo[@]}"
do
	echo $i
	Rscript --vanilla BF_rf_sim_results.R ${i} > ${resultspath}"/R_out/"${i}"_RF_R.out" 2> ${resultspath}"/R_error/"${i}"_RF_R.error" & echo $!
	sleep 1m
done
