library(LSD)
library(limma)
library(maDB)
library(KernSmooth)

normalize.results<-
  function(file, out, log=F, cols = c(1,2))
  {
    data = read.table(file,skip=2,header=F)
    data = data[,cols]
    if(log)
      data = 2^(data)
    
    png(out,height=1000,width=1000)
    layout(matrix(c(1,1,2,2,1,1,2,2,0,3,3,0,0,3,3,0), 4, 4, byrow=TRUE), respect=TRUE)
    par(cex=1.4)
    heatmapMA(data[,1],data[,2]); abline(0,0,col="green"); abline(0,9999999,col="green");
    heatscatter(log(data[,2]),log(data[,1]),xlab="Input",ylab="IP"); abline(0,1,col="green");
    h <- hist(log(data[,1]/data[,2]),breaks=1000,plot=F)
    plot(h$mid,h$counts,xlab="Normalized log-ratio",ylab="Counts",type='l')
    dev.off()
  }

modifiedLoess <-
  function(raw,cutoff=.2,kern_size=250)
  {
    
    M = log(raw$R/raw$G)
    A = (log(raw$R)+log(raw$G))/2
    weights = rep(0,length(M))

    #build a 2d kernel of the MA plot
    kern = bkde2D(cbind(A,M),c(var(A),var(M)),c(kern_size,kern_size))
    a_bw = kern$x1[2]-kern$x1[1]
    m_bw = kern$x2[2]-kern$x2[1]
    maxrow <- which(apply(kern$fhat==max(kern$fhat),1,sum)!=0)
    maxcol <- which(apply(kern$fhat==max(kern$fhat),2,sum)!=0) 
    s = kern$fhat[maxrow,maxcol]*a_bw*m_bw
    cur_max = max(kern$fhat)
    prog = cutoff/10
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
  select = data$genes$ControlType==0
  #gc.bins = sapply(data$genes$Sequence[select],countGC)
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
  #cy5 = cy5[select]
  #cy3 = cy3[select]

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
#     data$M[select] = cy5-cy3
#     data$A[select] = (cy5+cy3)/2
    data$M = cy5-cy3
    data$A = (cy5+cy3)/2
  }
  else
  {
#     data$R[select] = 2^cy5
#     data$G[select] = 2^cy3
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

# normalize.results<-
#   function(rfile, nfile, out)
#   {
#     raw = read.table(rfile,skip=2,header=F)
#     norm = read.table(nfile,skip=2,header=F)
#     #norm = normalize.loess(2^(norm[,c(4,5)]))
#     norm = 2^(norm[,c(4,5)])
#     
#     png(out,height=1000,width=1000)
#     par(mfrow=c(2,2))
#     heatmapMA(raw[,3],raw[,4]); abline(0,0); abline(0,99999999999);
#     heatscatter(log(raw[,4]),log(raw[,3]),xlab="Input",ylab="IP");
#     heatmapMA(norm[,1],norm[,2]); abline(0,0); abline(0,99999999999);
#     heatscatter(log(norm[,2]),log(norm[,1]),xlab="Input",ylab="IP");
#     dev.off()
#   }

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

normalizeCompareMA<-
function(cc,ccNorm,ccds,ccdsNorm,title="Normalize Comparison",method="",...)
{
	par(mfrow=c(2,2))
	heatmapMA(cc[,1],cc[,2],paste(title,"UnNormalized Sample 1"),...)
	heatmapMA(ccNorm[,1],ccNorm[,2],paste(title,method,"Normalized Sample 1"),...)
	heatmapMA(ccds[,1],ccds[,2],paste(title,"UnNormalized Sample 2 (Dye-Swap)"),...)
	heatmapMA(ccdsNorm[,1],ccdsNorm[,2],paste(title,method,"Normalized Sample 2 (Dye-Swap)"),...)
}

plotDensity<-
function(df,title="",group=F,condition=F)
{
	if(group)
	{
		qplot(ratios, ..density.., data = df, geom = "freqpoly", binwidth = .01, colour = normalized,cex=5) + labs(title = title)
	}
	else if(condition)
	{
		qplot(ratios, ..density.., data = df, geom = "freqpoly", binwidth = .01, colour = conditions) + labs(title = title)
	}
	else
	{
		qplot(ratios, ..density.., data = df, geom = "freqpoly", binwidth = .01, colour = experiments) + labs(title = title)
	}
}

ratioDensityPlot<-
function(ratios,title)
{
	h<-hist(ratios,breaks=1000,plot=F)
	plot(h$mid,h$counts,xlab="Ratio Density",ylab="Counts",main=title,type='l')
}

enrichment<-
function(ip,input,applyLog=F)
{
	if(applyLog)
		enrichment = log(ip)-log(input)
	else
		enrichment = ip-input
	pos = enrichment[enrichment>0]
	neg = abs(enrichment[enrichment<0])
	return (sum(pos)/sum(neg))
}


loessNormalize<-
function(procR,procG)
{
	return (normalize.loess(cbind(procR,procG)))
}

loessNormalizeComparison<-
function(rawFile,rawdsFile,outfile,raw=NULL,rawds=NULL,description="",rCol="rMedianSignal",gCol="gMedianSignal")
{
	cat("Loading raw files...\n")
	if(is.null(raw)){
		raw<-loadRawCc(rawFile)}
	if(is.null(rawds))
		rawds<-loadRawCc(rawdsFile)
	coords<-loadCoords("data/Chr_coords.txt")
	
	cat("Cleaning raw files...\n")
	proc<-cleanRawCc(coords,raw)
	procds<-cleanRawCc(coords,rawds,T)
	
	cat("Drawing graphs...\n")
	jpeg(outfile,height=800,width=800)
	normalizeCompare(cbind(proc[,rCol],proc[,gCol]),cbind(procds[,gCol],procds[,rCol]),loessNormalize,description,"Loess")
	dev.off()
}

ma2cNormalizeCompare<-
function(rawfile,normfile,rawdsfile,normdsfile,outfile="",description="")
{

	print("Raw")
	raw = read.table(rawfile,header=T)
	print("Norm")
	norm = read.table(normfile,header=T)
	rawds = read.table(rawdsfile,header=T)
	normds = read.table(normdsfile,header=T)
	
	cat (paste("MA2C comparison of",description,"\n"))
	
	enrichment = log(raw$IP) - log(raw$IN)
	pos = enrichment[enrichment>0]
	neg = abs(enrichment[enrichment<0])
	cat (paste("Raw enrichment:",sum(pos)/sum(neg),"\n"))
	
	enrichment = norm$IP - norm$IN
	pos = enrichment[enrichment>0]
	neg = abs(enrichment[enrichment<0])
	cat (paste("Normalized enrichment:",sum(pos)/sum(neg),"\n"))
	
	enrichment = log(rawds$IP) - log(rawds$IN)
	pos = enrichment[enrichment>0]
	neg = abs(enrichment[enrichment<0])
	cat (paste("Raw Dye Swap enrichment:",sum(pos)/sum(neg),"\n"))
	
	enrichment = normds$IP - normds$IN
	pos = enrichment[enrichment>0]
	neg = abs(enrichment[enrichment<0])
	cat (paste("Normalized Dye Swap enrichment:",sum(pos)/sum(neg),"\n"))
	
	if(! outfile == "")
	{
		jpeg(outfile,width=800,height=800)
		normalizeCompareMA(log2(subset(raw,select=c(IP,INPUT))), subset(norm,select=c(IP,INPUT)), log2(subset(rawds,select=c(IP,INPUT))), subset(normds,select=c(IP,INPUT)),description,"MA2C",applyLog = F)
		dev.off()
	}
}



ma2cMABatch<-
function(targets,datasets)
{
	for(val in targets)
	{
		raw = as.character(datasets[datasets$id==val & datasets$type=="raw" & datasets$chip=="reg","file"])
		norm = as.character(datasets[datasets$id==val & datasets$type=="norm" & datasets$chip=="reg","file"])
		rawds = as.character(datasets[datasets$id==val & datasets$type=="raw" & datasets$chip=="ds","file"])
		normds = as.character(datasets[datasets$id==val & datasets$type=="norm" & datasets$chip=="ds","file"])
		
		cat(raw,norm,rawds,normds,val,"\n")
		cat(paste(val,"_ma2c.jpeg",sep=""))
		ma2cNormalizeCompare(raw,norm,rawds,normds,paste(val,"_ma2c.jpeg",sep=""),description=as.character(val))
	}
}

normalizeCompare<-
function(cc,ccds,f,title="",functionDescription="")
{
	ccNorm=f(cc[,1],cc[,2])
	ccdsNorm=f(ccds[,1],ccds[,2])
	normalizeCompareMA(cc,ccNorm,ccds,ccdsNorm,title,functionDescription)
}

