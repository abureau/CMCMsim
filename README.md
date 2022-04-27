# CMCMsim: Analysis of case-mother - control-mother simulated data

This repository contains the code to reproduce the simulation results presented in Tables 1 and 2 and Figures 1 and 2 of the manuscript

**Methods and software to analyze gene-environment interactions under a case-mother control-mother design with partially missing child genotype**

by Alexandre Bureau, Yuang Tian, Patrick Levallois, Yves Giguere, Jinbo Chen, Hong Zhang

## Execution instructions

- Set the number of cores on your machine for parallel computing in the header of the scripts. 
- Execute the scripts Data_Generation_dep.R, Data_Generation_depnull.R, Data_Generation_dep02.R, Data_Generation_dep02MNAR.R, Data_Generation_indep.R and Data_Generation_indepnull.R, using for instance
 
        Rscript Data_Generation_dep.R
         
  This will create objects containing the results in the R environment *.RData* file of the folder where you executed the scripts.
- Execute the scripts Simulation_processing_Haplin.R and Simulation_processing_Haplin_null.R using for instance

        Rscript Simulation_processing_Haplin.R
        
  This will create the text file res.tex which was the basis for Table 1 of the manuscript.

- Execute the script display_bias_power.R using for instance

        Rscript display_bias_power.R

  This will create the text file res_power_raw.tex which was the basis for Table 2 and the pdf files bias_dep1.pdf and emp_se_dep1.pdf which are Figures 1 and 2 of the manuscript.
  