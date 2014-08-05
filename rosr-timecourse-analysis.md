Analysis of RosR Timcourse ChIP data
===================================


```r
rm(list = ls())
```


# Load Libraries and Data


```r
library(ggplot2)
source("lib/medichi.utils.R")
```

```
## MeDiChI utils. (c) David J Reiss, ISB.
## Please email dreiss@systemsbiology.org for questions or comments.
## Loading function 'get.strongest.hits'
## Loading function 'get.data'
## Loading function 'compare.chip.hits'
## Loading function 'get.genes.hit'
## Loading function 'chip.hits.to.sif'
## Loading function 'plot.significant.peaks'
## Loading function 'fraction.in.coding.rgns'
## Loading function 'medichi.clone.files'
## Loading function 'medichi.to.genome.browser' (thanks, Tie!)
```

```r
source("lib/chipchip.utils.R")
data("halo.lowres", package = "MeDiChI")
operons = read.delim("data/nrc1_operons.csv")
gene.coords$Gene_Name = gene.coords$canonical_Name

data_dir = ""
experiment = "rosr/"

targets = c("noH2O2.fits.Rdata", "wH2O2_10m.fits.Rdata", "wH2O2_20m.fits.Rdata", 
    "wH2O2_60m.fits.Rdata", "0258_22.fits.Rdata", "0258_11.fits.Rdata", "0258_13.fits.Rdata")

targets = sapply(targets, function(x) paste(data_dir, experiment, x, sep = ""))

datasets = c("noH2O2", "wH2O2_10m", "wH2O2_20m", "wH2O2_60m", "0258_22", "0258_11", 
    "0258_13")

conditions = c("noH2O2", "wH2O2_10m", "wH2O2_20m", "wH2O2_60m", "noH2O2", "noH2O2", 
    "noH2O2")
```



```r

hits = combineHits(targets, datasets, conditions, 0.2, 250)
```

```
## N.COEFFS = 389 
## N.COEFFS = 152 
## N.COEFFS = 205 
## N.COEFFS = 382 
## N.COEFFS = 673 
## N.COEFFS = 515 
## N.COEFFS = 232
```

```r
hits.combined = combineGeneHitsByDataset(hits, datasets)
hits.condition = combineGeneHitsByCondition(hits, conditions)
```



```r
# Convert NA pvalues to 10^(-.5) ~ .31 and NA intensities to the minimum
# value
pval.cols = grep("pval", colnames(hits.condition))
intens.cols = grep("intens", colnames(hits.condition))
intens.avg.cols = grep("intens.avg", colnames(hits.condition))
for (col in pval.cols) {
    hits.condition[is.na(hits.condition[, col]), col] = 10^(-0.5)
}
for (col in intens.cols) {
    hits.condition[is.na(hits.condition[, col]), col] = min(hits.condition[, 
        intens.cols], na.rm = T)
}

# use rows with at least one pvalue <.05
hits.condition = hits.condition[apply(hits.condition[, pval.cols], 1, function(x) {
    any(x < 0.05)
}), ]

# add operons to conditional hits
hits.condition = duplicateOperonRows(hits.condition, operons)

# save table
write.table(hits.condition, file = "rosr/table_s2.csv", sep = "\t", quote = F, 
    row.names = F)
```

