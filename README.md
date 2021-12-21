# CMCMsim: Analysis of case-mother - control-mother simulated data

This repository contains the code to reproduce the simulations presented in Table 1 of the manuscript

**Methods and software to analyze gene-environment interactions under a case-mother control-mother design with partially missing child genotype**

by Alexandre Bureau, Yuang Tian, Patrick Levallois, Yves Giguere, Jinbo Chen, Hong Zhang

## Execution instructions

- Set the number of cores on your machine for parallel computing in the header of the scripts Data_Generation_dep.R and Data_Generation_indep.R
- Execute the scripts Data_Generation_dep.R and Data_Generation_indep.R, using for instance
 
        Rscript Data_Generation_dep.R
        Rscript Data_Generation_indep.R
        
  This will create the objects *CMCM.dep01.res* and *CMCM.indep01.res* in the R environment *.RData* file of the folder where you executed the scripts.
- Execute the script Simulation_processing.R using for instance

        Rscript Simulation_processing.R
        
  This will create the text file res.tex which was the basis for Table 1 of the manuscript.
