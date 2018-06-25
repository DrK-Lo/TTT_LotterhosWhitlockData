#!/bin/bash
set -e
set -u
set -o pipefail

######################################################################################
#
# repORD.sh
# Brett Ford
# Created 20180302
#
# This script runs multiple ordination R scripts across mutiple nodes on a cluster
#
# Usage: bash repORD.sh
#
#######################################################################################

# Specify directories
mypath="/home/br.ford/br.ford_remote/r_projects/forester_sim_code/TTT_LotterhosWhitlockData/src"
outpath="/home/br.ford/br.ford_remote/r_projects/forester_sim_code/TTT_LotterhosWhitlockData/forester_results"
cd $mypath

#Run across 36 processors
start=55
finish=72
echo $start $finish


##############
#### run R script
#############
echo "Running R scripts"
for i in $(seq $start $finish)
do
        echo $i
        Rscript --vanilla ordination_bylocus.R ${i} > ${outpath}"/R_out/"${i}"_ORD_R.out" 2> ${outpath}"/R_error/"${i}"_ORD_R.error" & echo $!
        sleep 1m
done
