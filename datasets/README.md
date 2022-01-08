The data sets in this folder contain the data for replication of the study.


## csv files

- BC.csv: A table containing the raw counts of background checks in each state between January 1999 and March 2021. (source: FBI NICS)
- BCC.csv: A table containing the raw counts of background checks in each state between January 1999 and March 2021, divided by the population size in the state that month. (source: FBI NICS and US Census Bureau)
- BC_2000_2017_sa_dt.csv: A table containing the processed counts of background checks in each state between January 2000 and December 2017, after seasonal adjustment and detrending.
- BCC_2000_2017_sa_dt.csv: A table containing the processed counts of background checks per capita in each state between January 2000 and December 2017, after seasonal adjustment and detrending.
- FO.csv: A table comtaining the model's output for firearm ownership (fraction of firearm owners) in each state between January 2000 and October 2019.
- FO_2000_2017_sa_dt.csv: A table containing processed firearm ownership in each state between January 2000 and October 2017, after seasonal adjustment and detrending.
- MO_1999_2017.csv: A table summarizing the number of articles on regulation of firearms published by the New York Times and Washington Post between January 1999 and December 2017. (source: ProQuest)
- MS_1999_2017.csv: A table summarizing the mass shooting events that took place between January 1999 and December 2017. (source: Mother Jones)
- S.csv: A table containing the fraction of suicides committed with firearms in each state between November 1999 and October 2019. (source: CDC Wonder)
- S_2000_2017.csv: A table containing the processed fraction of suicides committed with firearms in each state between January 2000 and December 2017, after seasonal adjustment and detrending.


## mat files

- A.mat: A table containing the area of each state. (source: US Census Bureau)
- BC.mat: A table containing the monthly counts of background checks in each state between January 1999 and March 2021. (source: FBI NICS)
- BCC.mat: A table containing the monthly counts of background checks in each state between January 1999 and March 2021, divided by the population size in the state that month. (source: FBI NICS and US Census Bureau)
- D.mat: A symmetric matrix summarizing the distances between the centroids of states, measured in kilometers. (source: N. Abaid et al. "The effect of geography and citizen behavior on motorvehicle deaths in the United States," PLoS ONE 10 (2015) e0123339)
- data.mat: A table summarizing yearly data on national- and state-level population size, firearm ownership, background checks, and suicides. (sources: US Census Bureau, Gallup, FBI NICS and CDC Wonder)
- GDP.mat: A table containing the gross domestic product of each state between 1997-2020. (source: Federal Reserve Economic Data)
- S.mat: A table containing the monthly counts of suicides and suicides with firearms as well as the fraction of suicides with firearms, in each state between January 1999 and December 2019. (source: CDC Wonder)


## scripts

- plot_time_series.m: A script that plots the time series of each variable and reproduces figures 1, 2, and 3 in the manuscript.
