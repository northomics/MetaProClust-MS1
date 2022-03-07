# MetaProClust-MS1
*MetaProClust-MS1: Gut microbiome metaproteomics clustering using rapid MS1 profiling*

## R notebooks

MPC-MS1 is presented with example data through a series of [R notebooks](https://github.com/northomics/MetaProClust-MS1/tree/main/Datasets1_2_3-drug_microbiome_interactions) that can be run through [RStudio](https://rstudio.com). The code can be run through an R notebook with a user's own data, or can be customized using modular scripts. 

The R notebooks include all relevent code with R scripts and calls for Python programs using the data from our manuscript. The code for all figures is also included. 

## Code

Command line scripts called by R notebooks can be found in [`bin/`](https://github.com/northomics/MetaProClust-MS1/tree/main/bin). This also includes two Python virtual environments required for [PRECISE ICA matrix decomposition](https://github.com/SBRG/precise-db) and k-medoid clustering.
