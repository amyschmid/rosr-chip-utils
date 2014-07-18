Analysis of RosR Timcourse ChIP data
===================================


```r
rm(list = ls())
```


# Load Libraries and Data


```r
library(ggplot2)
source("chipchip.utils.R")
```

```
## Loading required package: BiocGenerics
```

```
## Attaching package: 'BiocGenerics'
```

```
## The following object(s) are masked from 'package:stats':
## 
## xtabs
```

```
## The following object(s) are masked from 'package:base':
## 
## anyDuplicated, cbind, colnames, duplicated, eval, Filter, Find, get,
## intersect, lapply, Map, mapply, mget, order, paste, pmax, pmax.int, pmin,
## pmin.int, Position, rbind, Reduce, rep.int, rownames, sapply, setdiff,
## table, tapply, union, unique
```

```
## Loading required package: Biobase
```

```
## Welcome to Bioconductor
## 
## Vignettes contain introductory material; view with 'browseVignettes()'. To
## cite Bioconductor, see 'citation("Biobase")', and for packages
## 'citation("pkgname")'.
```

```
## Loading required package: lars
```

```
## Loaded lars 1.1
```

```
## Loading required package: quadprog
```

```
## Loading required package: corpcor
```

```
## Loading required package: Matrix
```

```
## Loading required package: lattice
```

```
## Attaching package: 'Matrix'
```

```
## The following object(s) are masked from 'package:stats':
## 
## toeplitz
```

```
## Loading MeDiChI, version 0.4.0 (Tue Oct 12 11:20:57 2010)
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
data("halo.lowres", package = "MeDiChI")
operons = read.delim("data/nrc1_operons.csv")
gene.coords$Gene_Name = gene.coords$canonical_Name

data_dir = "preprocessing/output/"
experiment = "0258_none_densityLoess_max100/"

targets = c("noH2O2.fits.Rdata", "wH2O2_10m.fits.Rdata", "wH2O2_20m.fits.Rdata", 
    "wH2O2_60m.fits.Rdata", "0258_22.fits.Rdata", "0258_11.fits.Rdata", "0258_12.fits.Rdata", 
    "0258_13.fits.Rdata")

targets = sapply(targets, function(x) paste(data_dir, experiment, x, sep = ""))

datasets = c("noH2O2", "wH2O2_10m", "wH2O2_20m", "wH2O2_60m", "0258_22", "0258_11", 
    "0258_12", "0258_13")

conditions = c("noH2O2", "wH2O2_10m", "wH2O2_20m", "wH2O2_60m", "noH2O2", "noH2O2", 
    "noH2O2", "noH2O2")
```



```r

hits = combineHits(targets, datasets, conditions, 0.2, 250)
```

```
## N.COEFFS = 404 
## N.COEFFS = 144 
## N.COEFFS = 213 
## N.COEFFS = 392 
## N.COEFFS = 630 
## N.COEFFS = 480 
## N.COEFFS = 349 
## N.COEFFS = 220
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
```


# Compare results with previous analysis


```r

# load table from Amy sent 7/10
old.hits = read.table("results/tables/ST1_07102013.csv", sep = "\t", header = T)
old.genes = unlist(strsplit(as.character(old.hits$genes..not.operon.expanded.), 
    " "))
old.genes = unique(old.genes)

all.genes = unique(c(old.genes, as.character(hits.condition$Gene)))
comparison = data.frame(gene = all.genes, group = ifelse(all.genes %in% old.genes & 
    all.genes %in% hits.condition$Gene, "both", ifelse(all.genes %in% old.genes, 
    "old", "new")))

# load ma2c hits from previous pipeline
load("results/data/ma2c_genes.RData")

comparison = cbind(comparison, in.ma2c = comparison$gene %in% ma2c_genes)

c <- ggplot(comparison, aes(group, fill = in.ma2c))
c + geom_bar()
```

![plot of chunk compare to previous hits](figure/compare to previous hits.png) 





There are 679 unique genes hit across all datasets (1644 total hits).


```r

h2o2_expression = cbind(read.table("data/rosr/expression/1a_h2o2.csv", header = T, 
    sep = "\t"), cluster = "1a")
h2o2_expression = rbind(h2o2_expression, cbind(read.table("data/rosr/expression/1b_h2o2.csv", 
    sep = "\t", header = T), cluster = "1b"))
h2o2_expression = rbind(h2o2_expression, cbind(read.table("data/rosr/expression/2a_h2o2.csv", 
    header = T, sep = "\t"), cluster = "2a"))
h2o2_expression = rbind(h2o2_expression, cbind(read.table("data/rosr/expression/2b_early_h2o2.csv", 
    header = T, sep = "\t"), cluster = "2b_early"))
h2o2_expression = rbind(h2o2_expression, cbind(read.table("data/rosr/expression/2b_late_h2o2.csv", 
    header = T, sep = "\t"), cluster = "2b_late"))
h2o2_expression = rbind(h2o2_expression, cbind(read.table("data/rosr/expression/3_h2o2.csv", 
    header = T, sep = "\t"), cluster = "3"))
h2o2_expression = rbind(h2o2_expression, cbind(read.table("data/rosr/expression/4a_h2o2.csv", 
    header = T, sep = "\t"), cluster = "4a"))
h2o2_expression = rbind(h2o2_expression, cbind(read.table("data/rosr/expression/4b_h2o2.csv", 
    header = T, sep = "\t"), cluster = "4b"))
h2o2_expression_dynamic = cbind(read.table("data/rosr/expression/dynamic_h2o2.csv", 
    header = T, sep = "\t"), Standard = "", cluster = "dynamic")
colnames(h2o2_expression_dynamic)[3] <- "genbank.gi"
colnames(h2o2_expression_dynamic)[7] <- "protein.fxn"
h2o2_expression = rbind(h2o2_expression, h2o2_expression_dynamic[, c(1, 2, 3, 
    4, 6, 5, 7, 9)])

deg = h2o2_expression

hits.combined = duplicateOperonRows(hits.combined, operons)
hits.allGenes = hits.combined

hits.combined = hits.combined[hits.combined$Gene %in% deg$ORF, ]
hits.combined = merge(hits.combined, deg, by.x = "Gene", by.y = "ORF")
hits.allGenes = merge(hits.allGenes, deg, by.x = "Gene", by.y = "ORF", all.x = T)
# intens = hits.combined[,grep('*.intens',names(hits.combined))]
# intens[is.na(intens)]=0
# 
# pval = hits.combined[,grep('*.p.val',names(hits.combined))]
# pval[is.na(pval)] = 1
```


There are 134 hits from differentially expressed genes.

# Visualize DEG ChIP timecourse


```r
expression = read.table("data/rosr/expression/allGenes_h2o2.csv", header = T)

expression_data = data.frame(gene = expression$gene)
for (tp in c("X.40", "X.20", "X0", "X10", "X20", "X40", "X60", "X80")) {
    # log ratio intensity
    ratio = rowMeans(expression[, grep(paste(tp, "_ura3*", sep = ""), colnames(expression))]) - 
        rowMeans(cbind(expression[, grep(paste(tp, "_0258*", sep = ""), colnames(expression))], 
            expression[, grep(paste(tp, "_0258*", sep = ""), colnames(expression))]))
    expression_data = cbind(expression_data, data.frame(ratio))
    colnames(expression_data)[ncol(expression_data)] = tp
}
```




```r
hits.condition = hits.condition[hits.condition$Gene %in% deg$ORF, ]

matplot(c(0, 10, 20, 60), t(hits.condition[, intens.avg.cols]), type = "l", 
    col = "black", xlab = "time", ylab = "ChIP intensity", lty = 1)
```

![plot of chunk visualize deg chip](figure/visualize deg chip.png) 



```r

k = 6
clusters = kmeans(hits.condition[, intens.avg.cols], k, iter.max = 50)

par(mfrow = c(k, 1))
for (i in 1:k) {
    matplot(c(0, 10, 20, 60), t(hits.condition[clusters$cluster == i, intens.avg.cols]), 
        type = "l", col = "grey", xlab = "time", ylab = "ChIP intensity", lty = 1)
    title(paste(clusters$size[i], "genes in cluster."))
    par(new = T)
    plot(c(0, 10, 20, 60), colMeans(hits.condition[clusters$cluster == i, intens.avg.cols]), 
        type = "l", lty = 1, lwd = 3, col = "black", xlab = "", ylab = "")
}
```

![plot of chunk cluster chip data](figure/cluster chip data.png) 



```r

# combine ge and chip data
combined.data = merge(expression_data, cbind(hits.condition$Gene, hits.condition[, 
    intens.avg.cols]), by.x = "gene", by.y = "hits.condition$Gene")
combined.data[, -1] = scale(combined.data[, -1])

k = 5
clusters = kmeans(combined.data[, -1], k, iter.max = 50)

par(mfrow = c(k, 1))
for (i in 1:k) {
    matplot(c(-40, -20, 0, 10, 20, 40, 60, 80), t(combined.data[clusters$cluster == 
        i, 2:9]), type = "l", col = "red", xlab = "time", ylab = "ChIP intensity", 
        lty = 1)
    par(new = T)
    matplot(c(0, 10, 20, 60), t(hits.condition[clusters$cluster == i, intens.avg.cols]), 
        type = "l", col = "grey", xlab = "time", ylab = "ChIP intensity", lty = 1, 
        xlim = c(-40, 80))
    # title(paste(clusters$size[i],'genes in cluster.')) par(new=T)
    # plot(c(0,10,20,60),colMeans(hits.condition[clusters$cluster==i,intens.avg.cols]),type='l',lty=1,lwd=3,col='black',xlab='',ylab='')
}
```

![plot of chunk cluster chip and deg](figure/cluster chip and deg.png) 

