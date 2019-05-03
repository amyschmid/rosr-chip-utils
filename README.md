# RosR ChIP-chip processing

ChIP-chip processing pipeline for Tonner et al. "A regulatory hierarchy controls the dynamic transcriptional response to extreme oxidative stress in archaea."

* This part of the repository contains the code necessary to recereate table S2 from the manuscript
* Raw data file can be downloaded from GEO: GSE58696

## Dependencies

* R (version 2.15.2):
* knitr
* limma - http://www.bioconductor.org/packages/release/bioc/html/limma.html
* maDB - http://www.bioconductor.org/packages//2.7/bioc/html/maDB.html
* MeDiChI - http://cran.r-project.org/src/contrib/Archive/MeDiChI/
* RColorBrewer
* KernSmooth
* yaml
* outliers
* LSD


## Steps to recreate table:

* download raw data from GEO, create a folder called 'data/chip/' and unpack into that directory
* run the command 'R CMD Sweave preprocessing.Rnw'
* Download table S5 from Sharma et al. BMC Genomics 2012, and save the tab 'all\_data\_H2O2' as data/rosr\_h2o2\_exp.csv 
* from R, run 'knit2html(rosr-timcourse-analysis.Rmd)' with the knitr package

# RosR ChIP-seq processing and comparison across multiple species of halophilic archaea.
ChIP-seq processing pipeline for Darnell and Hackley et al. "Functional divergence of the RosR transcription factor across species of halophilic archaea."

* This part of the repository contains the code necessary to recereate table X from the manuscript
* Raw data file can be downloaded from GEO: GSEXXXX

## Dependences

## Steps to recreate table X and figure Y:
