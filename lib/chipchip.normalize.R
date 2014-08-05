library(LSD)
library(limma)
library(maDB)
library(KernSmooth)

modifiedLoess <-
  function(raw,cutoff=.2,kern_size=250)
  {
    raw$R = raw$R + 1e-9
    raw$G = raw$G + 1e-9
    M = log(raw$R/raw$G)
    A = log(raw$R*raw$G)/2
    weights = rep(0,length(M))

    #build a 2d kernel of the MA plot
    kern = bkde2D(cbind(A,M),c(var(A),var(M)),c(kern_size,kern_size))
    a_bw = kern$x1[2]-kern$x1[1]
    m_bw = kern$x2[2]-kern$x2[1]
    maxrow <- which(apply(kern$fhat==max(kern$fhat),1,sum)!=0)
    maxcol <- which(apply(kern$fhat==max(kern$fhat),2,sum)!=0) 
    s = kern$fhat[maxrow,maxcol]*a_bw*m_bw
    cur_max = max(kern$fhat)
    while(s<cutoff)
    {
      
      cur_max = max(kern$fhat[kern$fhat<cur_max])
      maxrow <- which(apply(kern$fhat==cur_max,1,sum)!=0)
      maxcol <- which(apply(kern$fhat==cur_max,2,sum)!=0)
      s = s + kern$fhat[maxrow,maxcol]*a_bw*m_bw
      
      weights[M>kern$x2[maxcol] & M<kern$x2[maxcol+1] & A>kern$x1[maxrow] & A<kern$x1[maxrow+1]] = 1
    }
    
    norm = normalizeWithinArrays(raw, raw$printer, method="loess",weights=weights,iterations=20)
    norm
}

gcNormalize<-function(data,isNorm=F)
{
  gc.bins = sapply(data$genes$Sequence,countGC)
  gc.bins = gc.bins+1 #offset 0 for indexing
  if(isNorm)
  {
    cy5 = .5*(data$M + 2*data$A)
    cy3 = cy5 - data$A
  }
  else
  {
    cy5 = log(data$R,2)[,1]
    cy3 = log(data$G,2)[,1]
  }
  
  cy5bins = list()
  cy3bins = list()
  for(bin in unique(gc.bins))
  {
    cy5bins[bin] = 0
    cy3bins[bin] = 0
  }
  for(i in 1:length(cy5))
  {
    cy5bins[[gc.bins[i]]] = c(cy5bins[[gc.bins[i]]],cy5[i])
    cy3bins[[gc.bins[i]]] = c(cy3bins[[gc.bins[i]]],cy3[i])
  }
  for(bin in unique(gc.bins))
  {
    cy5bins[[bin]] = cy5bins[[bin]][-1]
    cy3bins[[bin]] = cy3bins[[bin]][-1]
  }

  for(i in 1:length(cy5))
  { 
    if(length(cy5bins[[gc.bins[i]]])==1)
    {
      next; #no normalization for length one bin
    }
    
    tCy5 = cy5bins[[gc.bins[i]]]
    tCy3 = cy3bins[[gc.bins[i]]]
    gcCov = cov(tCy5,tCy3)
    cy5mean = mean(tCy5)
    cy5std = var(tCy5)^.5
    cy3mean = mean(tCy3)
    cy3std = var(tCy3)^.5
    
    cy5[i] = (cy5[i] - cy5mean)/(cy5std^2+cy3std^2-2*gcCov)^.5
    cy3[i] = (cy3[i] - cy3mean)/(cy5std^2+cy3std^2-2*gcCov)^.5
  }
  #global signal scaling to one
  gstd = var(cy5-cy3)^.5
  cy5 = cy5/gstd
  cy3 = cy3/gstd
  
  if(isNorm)
  {
    data$M = cy5-cy3
    data$A = (cy5+cy3)/2
  }
  else
  {
    data$R = 2^cy5
    data$G = 2^cy3
  }
  data$Select = select
  data$GCbins = gc.bins
  data$cy5bins = cy5bins
  data$cy3bins = cy3bins
  
  data
}

countGC<-function(seq)
{
  seq = gsub("A","",seq)
  seq = gsub("T","",seq)
  nchar(seq)
}

#ma plotting
#r,g are norm of the format where each row is a genome location, each column is a different microarray
heatmapMA<-
function(r,g,title="HeatMap MA Plot",applyLog = T)
{
	if(applyLog)
		heatmaplot(log2(r), log2(g),main=title)
	else
		heatmaplot(r, g,main=title)
	title(xlab="A",ylab="M (R/G)")
}

custom_heatmapMA<-
function (x, y, r, pch = 19, cexplot = 0.5, ncol = 30, grid = 100, 
    colpal = "heat", xlim = NULL, ylim = NULL, xlab = "M-values", 
    ylab = "A-values", main = "heatmaplot", cor = FALSE, alpha = NULL, 
    add.contour = FALSE, nlevels = 10, color.contour = "black", 
    add.line = TRUE, color.line = "darkgrey", lwd = 2, greyscale = FALSE, 
    ...) 
{
    if (!is.vector(x) | !is.vector(y)) 
        stop("First two argument must be vectors !")
    if (length(x) != length(y)) 
        stop("Data vectors must be of the same length !")
    m = r
    a = (x + y)/2
    heatscatter(a, m, pch = pch, cexplot = cexplot, ncol = ncol, 
        grid = grid, colpal = colpal, xlim = xlim, ylim = ylim, 
        main = main, cor = cor, alpha = alpha, add.contour = add.contour, 
        nlevels = nlevels, color.contour = color.contour, greyscale = greyscale, 
        ...)
    if (add.line) {
        abline(h = 0, col = color.line, lwd = lwd)
	}
}