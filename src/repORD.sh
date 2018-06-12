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

# Specify directory
mypath="/Users/brettford/Desktop/Northeastern/coding/forester_simulation_code/forester_sim_code/TTT_LotterhosWhitlockData/src"
cd $mypath

#Run across 60 nodes
start=1
finish=60
echo $start $finish


##############
#### run R script
#############
echo "Running R scripts"
for i in $(seq $start $finish)
do
        echo $i
        Rscript --vanilla ordination_bylocus.R ${i} > ${mypath}"/R_out/"${i}"_ORD_R.out" 2> ${mypath}"/R_error/"${i}"_ORD_R.error" & echo $!
        sleep 1m
done
