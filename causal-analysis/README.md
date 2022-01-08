Matlab data tables and codes in this folder perform four causal analyses between three variables at a time: background checks/mass shootings/media output on firearm regulations; background checks per capita/mass shootings/media output on firearm regulations; suicides with firearms/mass shootings/media output on firearm regulations; and firearm ownership/mass shootings/media output on firearm regulations. Codes to perform delay analyses where time series are shifted with respect to others are also included. 

To run the analyses, place all files in the same directory and run "main.m". Within each section of the main script, the expected output is a table called "TE_results" that summarizes the transfer entropy between each pair of variables, the 95th percentile of the associated surrogate distribution, and the p-value for significance. In addition, a table called "state_TE" is generated, summarizing the results on a state level. Each analysis takes a few hours to run, as surrogate distributions are created, containing 50,000 values each. In order to test the code briefly with fewer iterations, one can change the value of num_iterations in line 12 before.

# Data

- nature.mat: A table containing the time series for national-level background checks, mass shootings, media output on firearm laws and regulations, media output on unemployment, and media output on shootings excluding firearm laws and regulations (source: M. Porfiri et al. 2019, "Media coverage and firearm acquisition in the aftermath of a mass shooting," Nature Human Behaviour 3, 913–921)
- BC_2000_2017_sa_dt.mat: A table containing national- and state-level time series of background checks between January 2000 and December 2017, following seasonal adjustment and detrending (source: M. Porfiri et al. 2019, "Media coverage and firearm acquisition in the aftermath of a mass shooting," Nature Human Behaviour 3, 913–921).
- BCC_2000_2017_sa_dt.mat: A table containing national- and state-level time series of background checks per capita between January 2000 and December 2017, following seasonal adjustment and detrending.
- data.mat: A table summarizing yearly data on national- and state-level population size, firearm ownership, background checks, and suicides. (sources: US Census Bureau, Gallup, FBI NICS and CDC Wonder)
- FO_2000_2017_sa_dt.mat: A table containing processed firearm ownership in each state between January 2000 and October 2017, after seasonal adjustment and detrending.
- FO_no_W_2000_2017_sa_dt.mat: A table containing the firearm ownership prediction of a non-spatial model in each state between January 2000 and October 2017, after seasonal adjustment and detrending.
- S_2000_2017_sa_dt.mat: A table containing national- and state-level time series of the fraction of suicides committed with firearms between January 2000 and December 2017, following seasonal adjustment and detrending.


# Codes

- main.m: The main script that calls all other tables and functions to compute transfer entropy in the different triads.
- compute_TE.m: A function that computes the observed transfer entropy between pairs of variables and shuffles time series to create surrogate distributions.
- conditional_TE.m: A function that computes conditional transfer entropy while preserving structures in the data.
- delay_analysis.m: A script that calls compute_delay_TE.m to compute conditional transfer entropy with delays for every pair of variables under consideration in this study. The script also plots the figures in figures 5 and 6 in the Supporting Information.
- compute_delay_TE.m: A function that computes the observed transfer entropy between pairs of variables with a delay on the source variable and on the variable conditioned upon. It calls the conditional_TE.m function.
