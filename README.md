# Simulations from Lotterhos and Whitlock 2015 and subsequent studies

The relative power of genome scans to detect local adaptation depends on sampling design and statistical method. Molecular Ecology, https://doi.org/10.1111/mec.13100

The simulations are from 90 randomly sampled populations. 10000 independent loci were simulated (9900 neutral and 100 under different strengths of selection).

These loci were analyzed by Lotterhos and Whitlock 2015 (LW2015), and subsequently by Forester et al. (2018; https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.14584). In this respository, we organized the results from both these studies, and in some cases re-ran some of the packages that had major updates since the LW2015 study. (NEED TO ADD DETAILS ABOUT WHAT WAS RE-RUN AND UPDATED)

Description of each folder or file and the data contained within:

* *"SchemeRandom1.txt"*: gives the x,y coordinates of each population in the data. These are plotted in *"SchemeRandom1_SampleLocs_90pops.pdf"*
  * "PopID" 
  * "X_Pops" 
  * "Y_Pops" 
  * "R90" 
  * "R60" 
  * "R30" 
  * "NSRangeTrans_30"

* "*simfiles*": the original simulations downloaded from the Dryad repository for LW2015- DELETED

* "*results*": the original results downloaded from the Dryad repository for LW2015- DELETED

* "*forester_simfiles*": This folder contains the LW2015 simulated data that was used in Forester et al. 2017. Because LW 2015 sampled allele frequencies, Forester et al. randomly sampled genotypes from those allele frequencies. So these are the datasets that should be used for subsequent analysis. This directory also includes the environment files for all simulations.

* "*forester_results*" This folder contains the reduced and appended cpval files. Columns from the original cpval were deleted due to repetitive information, and the results from Forester et al 2018 were appended to the remaining columns. Results from additional tests should be appended to these files. 
