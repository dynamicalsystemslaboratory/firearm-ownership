Matlab data tables and codes in this folder perform four causal analyses between three variables at a time: background checks/mass shootings/media output on firearm regulations; background checks per capita/mass shootings/media output on firearm regulations; suicides with firearms/mass shootings/media output on firearm regulations; and firearm ownership/mass shootings/media output on firearm regulations. To run the analyses, place all files in the same directory and run "main.m".

# Data

- nature.mat: A table containing the time series for national-level background checks, mass shootings, media output on firearm laws and regulations, media output on unemployment, and media output on shootings excluding firearm laws and regulations (source: M. Porfiri et al. 2019, "Media coverage and firearm acquisition in the aftermath of a mass shooting," Nature Human Behaviour 3, 913–921)
- restrictiveness.mat: A table containing firearm law restrictiveness scores for each state, based on the number of firearm laws implemented in that state, out of 133 possible laws (source: M. Porfiri et al. 2019, "Media coverage and firearm acquisition in the aftermath of a mass shooting," Nature Human Behaviour 3, 913–921).
- BC_2000_2017_sa_dt.mat: A table containing national- and state-level time series of background checks between January 2000 and December 2017, following seasonal adjustment and detrending (source: M. Porfiri et al. 2019, "Media coverage and firearm acquisition in the aftermath of a mass shooting," Nature Human Behaviour 3, 913–921).
- BCC_2000_2017_sa_dt.mat: A table containing national- and state-level time series of background checks per capita between January 2000 and December 2017, following seasonal adjustment and detrending.
- S_2000_2017_sa_dt.mat: A table containing national- and state-level time series of the fraction of suicides committed with firearms between January 2000 and December 2017, following seasonal adjustment and detrending.


# Codes

- main.m: The main script that calls all other tables and functions to compute transfer entropy in the different triads.
- compute_TE.m: A function that computes the observed transfer entropy between pairs of variables and shuffles time series to create surrogate distributions.
- conditional_TE.m: A function that computes conditional transfer entropy while preserving structures in the data.

