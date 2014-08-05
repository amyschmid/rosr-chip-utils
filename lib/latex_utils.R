#Sweave and latex
latex_escape <- function(filename) {
	pretty = paste("\\\\verb|", filename,"|")
	return(pretty);
}

sanatize<-function(dirty){
	return(gsub('_', "\\\\\\\\_", dirty))
}

table <- function(filenames) {
	x ="Input Files ";
	for(i in 1:length(filenames)) {
		x=paste(x, " &", latex_escape(files[i]), "\\\\\\\\", sep="")
	}
	return(x);
}
# 
# getPairedMAFiguresBg <- function(tmpDir, config)  {
# 	drawCommand=""
# 	n=config$count
# 	if(config$quiet$backgroundCorrection){
# 		return(" ")
# 	}
# 	for(i in 0:(n-1)) {
# 		FILENAME=paste(tmpDir,"/maAfterBg", i, ".png", sep="" )
# 		drawCommand = paste(drawCommand, 
# 				"\\\\begin{center}\n",
# 				"\\\\includegraphics[width=6in]{", FILENAME, "}\n",
# 				"\\\\caption{MA plot of ", sanatize(config$targets$Name[i]),
# 				" after background correction}\n",
# 				"\\\\end{center}\\\\pagebreak\n",
# 				sep="")
# 	}
# 	return(drawCommand)
# }
# 
# getMAFiguresWAN <- function(tmpDir, config)  {
# 	if(config$vsn | config$quiet$withinArrayNormalization) {
# 		return(" ")
# 	}
# 	drawCommand=""
# 	for(i in 1:config$count) {
# 		FILENAME=paste(config$files$tmpDir,"/maAfterWAN", i-1, ".png", sep="" )
# 		drawCommand = paste(drawCommand, 
# 				"\\\\begin{center}\n",
# 				"\\\\includegraphics[width=6in]{", FILENAME, "}\n",
# 				"\\\\caption{", sanatize(config$targets$Name[i]), 
# 				" after within array normalization by ", 
# 				config$algorithms$withinArray,
# 				"}\n \\\\end{center}\\\\pagebreak\n",
# 				sep="")
# 	}
# 	return(drawCommand)
# }
# 
# getPairedMAFiguresBAN <- function(config)  {
# 	tmpDir = config$files$tmpDir
# 	n      = config$count
# 	method = config$algorithms$betweenArray
# 	
# 	drawCommand=""
# 	for(i in 0:(n-1)) {
# 		FILENAME=paste(tmpDir,"/maAfterBAN", i, ".png", sep="" )
# 		drawCommand = paste(drawCommand, 
# 				"\\\\begin{figure}[htb!]\n",
# 				"\\\\centering\n",
# 				"\\\\includegraphics[width=6in]{", FILENAME, "}\n",
# 				"\\\\caption{MA plot of ",sanatize(config$targets$Name[i]), 
# 				" after between array normalization by ", method, "}\n",
# 				"\\\\end{figure}\\\\pagebreak\n",
# 				sep="")
# 	}
# 	return(drawCommand)
# }

getMedichiFigures <- function(config) {
  drawCommand=""
  for(condition in config$conditions){
    drawCommand = paste(drawCommand,
                        "\n\n",
                        "\\\\subsection{", sanatize(condition),"}\n")
    FILENAME=paste(config$experiment,"/", condition, "_medichi_hits_chr.png",sep="")
    drawCommand = paste(drawCommand, 
                        "\\\\begin{figure}[H]\n",
                        "\\\\centering\n",
                        "\\\\includegraphics[width=6in]{", FILENAME, "}\n",
                        "MeDiChI peaks of ",sanatize(condition), 
                        " on the main chromosome\n",
                        "\\\\end{figure}\n",
                        sep="")
    
    FILENAME=paste(config$experiment,"/", condition, "_medichi_hits_pNRC100.png",sep="")
    drawCommand = paste(drawCommand, 
                        "\\\\begin{figure}[H]\n",
                        "\\\\centering\n",
                        "\\\\includegraphics[width=6in]{", FILENAME, "}\n",
                        "MeDiChI peaks of ",sanatize(condition), 
                        " on pNRC100\n",
                        "\\\\end{figure}\n",
                        sep="")
    
    FILENAME=paste(config$experiment,"/", condition, "_medichi_hits_pNRC200.png",sep="")
    drawCommand = paste(drawCommand, 
                        "\\\\begin{figure}[H]\n",
                        "\\\\centering\n",
                        "\\\\includegraphics[width=6in]{", FILENAME, "}\n",
                        "MeDiChI peaks of ",sanatize(condition), 
                        " on pNRC200\n",
                        "\\\\end{figure}\\\\pagebreak\n",
                        sep="")
  }
  return(drawCommand)
}

#For the document
getBackgroundText <- function(config) {
	return(config$bgText$main)
}

getMetadataText <- function(config){
  return(paste(
    "This is the analysis for the ",
    gsub("_",".",config$experiment),
    " experiments. ",sep=""));
}

# getBgMaPlotText <- function(config) {
# 	if(config$quiet$backgroundCorrection) {
# 		return(config$bgText$quiet)
# 	} else {
# 		return(config$bgText$verbose)
# 	}
# }

getExperimentResults <-function(config,Slides.raw,Slides.raw.bg,Slides.norm,Slides.norm.ban) {
  setRaw = newMadbSet(Slides.raw)
  setBg = newMadbSet(Slides.raw.bg)
  setWa = newMadbSet(Slides.norm)
  setBa = newMadbSet(Slides.norm.ban)
  
  outstr = ""
  for(i in 1:config$count)
  {
    FILENAME<-paste(config$experiment,"/",config$experiment, "_", i, ".png", sep="" )
    png(FILENAME,width=1500,height=1500)
    par(mfrow=c(2,2))
    drawMA(setRaw, r=2*i-1, g=2*i, 
           colramp=colorRampPalette(rev(brewer.pal(9,"Blues")[2:9])),main="Raw Data",cex=1.5)
    drawMA(setBg, r=2*i-1, g=2*i, 
           colramp=colorRampPalette(rev(brewer.pal(9,"Blues")[2:9])),main="Background Corrected",cex=1.5)
    drawMA(setWa, r=2*i-1, g=2*i, 
           colramp=colorRampPalette(rev(brewer.pal(9,"Blues")[2:9])),main="Within Array Normalization",cex=1.5)
    drawMA(setBa, r=2*i-1, g=2*i, 
           colramp=colorRampPalette(rev(brewer.pal(9,"Blues")[2:9])),main="Between Array Normalization",cex=1.5)
    dev.off()
    
    outstr= paste(outstr,
      "\\\\subsection{", sanatize(config$slides[[i]]$name)," Results}",
      "\\\\begin{center}\n",
      "\\\\includegraphics[scale=.6]{",  FILENAME,"}\\\\\\\\ \n",
      "Results of all processing steps; background correction ",
      "within array normalization, and between array normalization.",
      "\\\\end{center}\n\\\\pagebreak\n",
                  "\n\n",
                  sep="");
  }
  return(outstr)
}

getBoxplotCommandRaw <-function(config) {
  if(config$vsn) {return(" ")}
  return(paste(
    "\\\\begin{figure}[H]\n",
    "\\\\centering\n",
    "\\\\includegraphics{",	config$files$expDir, "/boxPlotRaw.png}\n",
    "\\\\caption{Boxplots of the signal intensities of each signal ",
    "channel of the microarrays. Raw data before any correction or normalization.}\n",
    "\\\\end{figure}\n\\\\pagebreak\n",sep=""));
}


getBoxplotCommandBg <-function(config) {
	if(config$vsn) {return(" ")}
	return(paste(
					"\\\\begin{figure}[H]\n",
					"\\\\centering\n",
					"\\\\includegraphics{",	config$files$expDir, "/boxPlotBg.png}\n",
					"\\\\caption{Boxplots of the signal intensities of each signal ",
					"channel of the microarrays. Raw data after background correction.}",
					"\\\\end{figure}\n\\\\pagebreak\n",sep=""));
}

getWithinArrayNormText <-function(config){
	if(config$vsn) {
		return(config$wanText$vsn)
	}
	return(gsub("xxx", config$algorithms$withinArray, config$wanText$other))
}

getBoxplotCommandWA <-function(config) {
	if(config$vsn) {return(" ")}
	return(paste(
          "\\\\begin{figure}[H]\n",
          "\\\\centering\n",
					"\\\\includegraphics{", config$files$expDir,
					"/boxPlotWA.png}\n \\\\caption{Boxplots of the signal intensities ",
					"of each signal channel of the microarrays after within array ",
					"normalization.}\\\\end{figure}\\\\pagebreak",sep=""));
}

getWanText <- function(config) {
	if(config$quiet$withinArrayNormalization) {
		return(config$wanText$quiet)
	} else {
		return(config$wanText$verbose)
	}
}

getMedichiText <- function(config) {
  return(paste("Run MeDiChI with paramters max.steps=",config$medichi$max.steps,
        ", fit.res=",config$medichi$fit.res,", n.boot=",config$medichi$n.boot,
        ", boot.sample.opt=",config$medichi$boot.sample.opt,".",sep=""))
}

#From carma web
if( !isGeneric("drawBoxplot") )
	setGeneric("drawBoxplot", function(x,new.plot=TRUE,xlim=NULL,at=0,exclude.flagged=FALSE,...)
				standardGeneric("drawBoxplot"))

setMethod("drawBoxplot","MadbSet",
		function(x,new.plot=TRUE,xlim=NULL,at=0,exclude.flagged=FALSE,...){
			if(exclude.flagged){
				W <- getWeights(x)
				for(i in 1:ncol(exprs(x))){
					exprs(x)[W[,i]==0,i] <- NA
				}
			}
			boxplots(x=exprs(x),new.plot=new.plot,xlim=xlim,at=at,...)
		}
)

boxplots <- function(x,new.plot=TRUE,xlim=NULL,at=0,col=NULL,log2.transform=FALSE,ylim=NULL,...){
	data <- x
#	if(class(x)=="exprSet" | class(x)=="EexprSet"){
#		data <- exprs(x)
#	}
	if(log2.transform){
		data <- log2(data)
	}
	if(is.null(ncol(data))){
		data <- as.matrix(data)
	}
	## remove those nasty infinite values
	for(i in 1:ncol(data)){
		data[is.infinite(data[,i]),i] <- NA
	}
	
	CN <- colnames(data)
	if(is.null(CN)){
		CN <- 1:ncol(data)
	}
	if(is.null(col)){
		col=rep(0,ncol(data))
	}
	if(length(col)!=ncol(data)){
		col=rep(col[1],ncol(data))
	}
	if(new.plot){
		if(is.null(xlim)){
			xlim=c(0.5,(ncol(data)+0.5))
		}
		par(omi=0.7*c(min(2,max(strwidth(CN,units="i"))),0,0,0),cex.axis=0.7)
		if(is.null(ylim)){
			ylim=c(min(data,na.rm=TRUE),max(data,na.rm=TRUE))
		}
		plot(1,1,pch=NA,ylab="expression",xaxt="n",ylim=ylim,xlim=xlim,xlab=NA)
	}
	axis(side=1,labels=CN,at=(1+at):(ncol(data)+at),las=3)
	for(i in 1:ncol(data)){
		boxplot(data[,i],at=(i+at),add=TRUE,col=col[i],range=0,...)
	}
}
