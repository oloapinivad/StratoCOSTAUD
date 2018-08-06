# StratoCostaud v0.1
## COherent STructures AUtomatic Detector for Stratocumulus LES simulations

July 2018

by P. Davini (CNR-ISAC, p.davini@isac.cnr.it)

----------------
## WHAT IS StratoCostaud

**StratoCostaud** is a set of R scripts aimed at detecting the coherent structures as updraft, downdraft and entrainment in Large Eddy Simulations (LES) data. To be used, you will need the two boundary layer and free tropospheric passive scalars as described in Davini et al. (2017). It has been created for UCLA-LES output so variable names and data structure must be arranged if you want to use a different source dataset. 

The method is based on the octants procedure which divides airmasses according to their physical properties and it is based on the two scalars and on vertical velocity. While for the boundary layer scalar and for the vertical velocity a simple standard anomalies from the slab-average is used, a more complex procedure is used to define a threshold for the free tropospheric scalar. This is done maximizing a cost function which describes the difference of the integrated heat flux carried by the entraiment and by the downdraft octant. For each octant statistics and various plots as in Davini et al. (2017) are provided.

----------------

## MAIN NOTES & REFERENCES

Be aware that this is a free scientific tool, then it may not be free of bugs. Please report any issue at p.davini@isac.cnr.it

In the unlikely hazardous case you decide to use these scripts please cite **StratoCOSTAUD** but more importantly please cite:

*Davini P., F. D’Andrea, S. Park, P. Gentine (2017). Coherent structures in large-eddy simulations of a non-precipitating stratocumulus-topped boundary layer Journal of the Atmospheric Sciences, 74(12), 4117–4137, DOI: 10.1175/JAS-D-17-0050.1*

----------------

## SOFTWARE REQUIREMENTS

- a. R version > 3.0
- b. R packages ncdf4, SDMTools, spam and fields
- c. Compiling environment (gcc) to compile R packages

R packages parallel, doParallel and foreach are optional if you want to enable parallel execution. 

-----------------

## EXECUTION TIMES

The CPU and memory cost for **StratoCOSTAUD** can be heavy considering the maximization of the cost function. 
On a 32-core Intel(R) Xeon(R) Gold @ 2.10GHz with 187 GB RAM, running with 4 cores in parallel mode with foreach package it takes about 120seconds to complete the maximizationon for a single time step in a 256x256x210 domain.  Using more than 4 cores seems pointless due to bad scaling.
About the twice the time for the execution of the whole `profile_computation.R` script.

-----------------

## HOW TO

* `config.R`: a prelimiary version of the config file is introduced. This should be edited to set the directories of input and outputs. However, in each file you will need to set the correct CODEDIR in order to find the scripts.
* `profile_computation.R`: this is the core script, which analyze a single timestep of the model extracting the required variable and saving in R format the profiles of the main variables according to each octant.
* `averaging_profiles.R`: this provides the final figures averaging all the outputs from the profile_computation.R file.
* `octant_figures.R`: this uses the threshold defined by the profile_computation.R to provide 2D horizontal and vertical sections as well as scatterplots. 

-----------------

## HISTORY

*v0.1 Jul 2018*
- First version on Github
- Cleanup of the code, ~30% faster bruteforce.fast() function   
- Zero order config file
