# RosR ChIP-chip processing

ChIP-chip processing pipeline for Tonner et al. "A regulatory hierarchy controls the dynamic transcriptional response to extreme oxidative stress in archaea."

* This repository contains the code necessary to recereate table S2 from the manuscript
* Raw data file can be downloaded from GEO: GSE58696

## Dependencies

* R
* Sweave - builtin into R environment
* Rmarkdown - http://rmarkdown.rstudio.com/


## Steps to recreate table:

* download raw data from GEO
* run preprocessing_0258_all_none_densityLoess_max100.Rnw with the command 'R CMD Sweave preprocessing_0258_all_none_densityLoess_max100.Rnw'
* run rosr-timcourse-analysis.Rmd with the R-markdown package
