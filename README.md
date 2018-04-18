# ventricular_myocyte_models

This folder contains basic code for different ventricular myocyte models. All equations, initial conditions, and parameters were taken from the original papers or online repositories. 

Please note the following:
1) If you are downloading only a single model, please also download the find_APD.m function. 
APD is calculated as difference between time of stimulus and time when voltage returns to -75 mV. 

2) The stimulus protocol can be changed by altering the variables in Step 2. 

3) The human models have the option of simulating different cell types (epi,M,endo). The default setting is endo and can be altered by changing the "celltype" variable in Step 1. 
