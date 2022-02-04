# SUITHER-calibration

This is the repository for the code used in the project of *Numerical Analysis for Partial Differential Equations* course held by Professor A. Quarteroni at Politecnico di Milano during academic year 2020/2021. The tutors of this work are Professor N. Parolini and Professor L. Dedè. \\
    We adapt a Bayesian AR Dirichlet-Multinomial model for COVID-19 data to the scenario of a differential model for modelling epidemic waves (SUITHER model). From one side, this goes in the direction to find a new way to calibrate it. From the other side, it represents in general a good tool for validating results and getting insights about the dynamics of the epidemics, not otherwise inferable through data.


## Code structure
    - Script.R: R-script containing all the steps to build the data set from raw data, estimating the missing time series, and to run the sampler
    - MCMC: folder containing the sampling function and a folder Initialization tables containing all the code needed to initialize table plus the file initial_tables.RData with the used initial tables
    -Posterior inference: folder contanining a script for performing posterior inference and a script for MCMC diagnostics plus the estimated values and the pics produced

    
    
