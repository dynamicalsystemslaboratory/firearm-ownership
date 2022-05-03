Matlab data tables and codes in this folder calibrate a spatial econometric model that estimates monthly state-level firearm ownership in the United States from data on background checks and suicides with firearms. We also include scripts that calibrate a non-spatial model and explore models with alternative weight schemes (reported in the Supporting Information). To run the analyses, place all files in the same directory and run "main.m". The expected output of the main script is a table called "firearm_ownership_monthly" that contains those estimations. Model calibration and inference take a couple of minutes to compelete.

# Data

- BC.mat: A table containing the monthly counts of background checks in each state between January 1999 and March 2021. (source: FBI NICS)
- D.mat: A symmetric matrix summarizing the distances between the centroids of states, measured in kilometers. (source: N. Abaid et al. "The effect of geography and citizen behavior on motorvehicle deaths in the United States," PLoS ONE 10 (2015) e0123339)
- S.mat: A table containing the monthly counts of suicides and suicides with firearms as well as the fraction of suicides with firearms, in each state between January 1999 and December 2019. (source: CDC Wonder)
- data.mat: A table summarizing yearly data on national- and state-level population size, firearm ownership, background checks, and suicides. (sources: US Census Bureau, Gallup, FBI NICS and CDC Wonder)
- A.mat: A table containing the area of each state. (source: US Census Bureau)
- GDP.mat: A table containing the gross domestic product of each state between 1997-2020. (source: Federal Reserve Economic Data)
- B.mat: A table containing the lengths of borders between each state (source: U.S. Geological Survey)

# Codes

- main.m: The main script that calls all other tables and functions to calibrate the model and infer monthly firearm ownership between January 2000 and October 2019.
- main_non_spatial.m: A script that calls the relevant tables and functions to calibrate a non-spatial model (W=0) and infer monthly firearm ownership between January 2000 and October 2019.
- test_Ws.m: A script that studies alternative models with difference spatial weight matrix. Each cell in the script 
- normw.m: A function that row-normalizes square matrices.
- calibrate_model.m: A function that calibrates the parameters of spatial models. This function was written by James P. LeSage and obtained from www.spatial-econometrics.com.
- invpd.m, lndetmc.m, f2_sar.m, hessian.m, f_sar.m: functions that are called by calibrate_model.m and support the estimation of parameters. These functions were written by James P. LeSage and obtained from www.spatial-econometrics.com.
