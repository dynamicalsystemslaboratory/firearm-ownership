Matlab data tables and codes in this folder calibrate a spatial econometric model that estimates monthly state-level firearm ownership in the United States from data on background checks and suicides with firearms. To run the analysis, place all files in the same directory and run "main.m".

# Data

- BC.mat: A table containing the monthly counts of background checks in each state between January 1999 and March 2021. (source: FBI NICS)
- D.mat: A symmetric matrix summarizing the distances between the centroids of states, measured in kilometers. (source: N. Abaid et al. "The effect of geography and citizen behavior on motorvehicle deaths in the United States," PLoS ONE 10 (2015) e0123339)
- S.mat: A table containing the monthly counts of suicides and suicides with firearms as well as the fraction of suicides with firearms, in each state between January 1999 and December 2019. (source: CDC Wonder)
- data.mat: A table summarizing yearly data on national- and state-level population size, firearm ownership, background checks, and suicides. (sources: US Census Bureau, Gallup, FBI NICS and CDC Wonder)

# Codes

- main.m: The main script that calls all other tables and functions to calibrate the model and infer monthly firearm ownership between January 2000 and October 2019.
- normw.m: A function that row-normalizes square matrices.
- Calibrate_model.m: A function that 
- 
