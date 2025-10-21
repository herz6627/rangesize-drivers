# rangesize-drivers

# ADD publication info

This repository contains all code associated with Herzog XXXX and associated file structure. All code is within the ./R/ directory. The output, figs, and data directories are empty place holder directories to allow users to easily see how the data was structured when writing the code. Files for the ./data/ directory can be found at the following links:

## APT011

## AWE012

## CGR023

## KFH011

## PBG011

## PPS011

## PVC021


The R scripts are numbered with assending numeric prefixes to indicate order that scripts should be run in. For example, 1_clean-PVC-data.R should be run before all other R scripts, followed by 2_range-size-estimation.R, etc..

Included R files are described below:

## 1_clean-PVC-data.R
Cleans and formats datasets.
## 2_range-size-estimation.R
Download and clean occurrence records, then estimate range size. Depends on functions in 2_func_mcp-sh.R
Specific iDigBio and GBIF records used to estimate range size can be found at XXXX
## 2_func_mcp-sh.R
Custom code largely based on the XXX R package. 
## 3_gh-range-cor.R
Estimate number of grasshopper consumers per plant species and correllate number of grasshopper consumers by plant range size.
## 4_analyses-sep-by-dat.R
Create linear models for colonization, change in cover, and extirpation using environmental drivers for each Konza plant species. Reliant on function in 4_functions.R
## 4_functions.R
Helper functions to estimate environmental effects on plant population dynamics.
## 5_additional-analyses.R
Using coefficients estimated in 4_ and range sizes estimated in 2_, find any linear trends in driver effect size and range size. Dependent on functions in 5_functions.R
## 5_functions.R
Helper functions to clean data, run linear models and make figures.
