# RosR ChIP-chip processing

ChIP-chip processing pipeline for Tonner et al. "A regulatory hierarchy controls the dynamic transcriptional response to extreme oxidative stress in archaea."

* This repository contains the code necessary to recereate table S2 from the manuscript
* Raw data file can be downloaded from GEO: GSE58696

## Dependencies

* R (version 2.15) - knitr, MeDiChI (available here: http://cran.r-project.org/src/contrib/Archive/MeDiChI/)


## Steps to recreate table:

* download raw data from GEO, create a folder called 'data/chip/' and unpack into that directory
* run the command 'R CMD Sweave preprocessing.Rnw'
* Download table S5 from Sharma et al. BMC Genomics 2012, and save the tab 'all\_data\_H2O2' as data/rosr\_h2o2\_exp.csv 
* from R, run 'knit2html(rosr-timcourse-analysis.Rmd)' with the knitr package
