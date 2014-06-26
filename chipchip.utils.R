library(affy)
library(ggplot2)
require(MeDiChI)

##################
##Datset functions
##################

allDatasets<-
function()
{
	return (read.table("data/datasets.txt",header=T))
}

datasets.agilent<-
function()
{
	data = allDatasets()
	return (data[data$source=="agilent",])
}

data.bed<-
function(file)
{
	bed<-read.table(file)
	names(bed)<-c("chr","start","end","name","score")
	bed
}

save.bed<-
function(bed,ofile)
{
	write.table(subset(bed,select=c("chr","start","end","name","score")),ofile,row.names=F,col.names=F,sep="\t",quote=F)
}

loadMA2C<-
function(target)
{
	return(read.table(target,header=T,skip=1))
}

load_ma2c_peaks<-
function(condition,conditionLocation)
{
	conditionFile = paste(conditionLocation,condition,"_peaks.xls",sep="")
	ma2c = read.table(conditionFile,header=T)
	names(ma2c)[8]<-"ma2c_PVALUE"
	convert = 10^(ma2c[,8]/-10)
	ma2c[,8] = convert
	ma2c
}

load_medichi_peaks<-
function(datasets,medichi_folder,cutoff=.05)
{
	peaks = data.frame()
	for( data in datasets )
	{
		#load the MeDiChI data
		load(paste(medichi_folder,"fits_",data,".RData",sep="")) #load 'fits' from memory, peaks stored in 'fits'
	  #load(paste(medichi_folder,file,sep="")) #load 'fits' from memory, peaks stored in 'fits'
		
		#find overlap
		medichi_coeffs = get.strongest.hits(fits,cutoff) #get peaks with p value less than cutoff
		pos = medichi_coeffs[,"position"]
		chr = ifelse(rownames(medichi_coeffs)=="Chr",1,ifelse(rownames(medichi_coeffs)=="pNRC100","PNRC100","PNRC200"))
		pvals = medichi_coeffs[,"p.value <"]
		intens = medichi_coeffs[,"intensity"]
		names=paste("peak",1:nrow(medichi_coeffs),"dataset",as.character(data),sep="_")
		peaks<-rbind(peaks,data.frame(Position=pos,Chr=chr,Pvalue=pvals,Intensity=intens,Name=names,Dataset=data))
	}
	peaks
}



load.medichi.peaks<-
  function(config,cutoff=.05)
  {
    peaks = data.frame()
    conditions = unique(sapply(config$slides,function(x){x$condition}))
    for( condition in conditions )
    {
      #load the MeDiChI data
      #load(paste(medichi_folder,"fits_",data,".RData",sep="")) #load 'fits' from memory, peaks stored in 'fits'
      load(paste(config$files$outDir,"/",config$experiment,"/",condition,".fits.Rdata",sep="")) #load 'fits' from memory, peaks stored in 'fits'
      
      #find overlap
      medichi_coeffs = get.strongest.hits(fits,cutoff) #get peaks with p value less than cutoff
      pos = medichi_coeffs[,"position"]
      chr = ifelse(rownames(medichi_coeffs)=="Chr",1,ifelse(rownames(medichi_coeffs)=="pNRC100","PNRC100","PNRC200"))
      pvals = medichi_coeffs[,"p.value <"]
      intens = medichi_coeffs[,"intensity"]
      names=paste("peak",1:nrow(medichi_coeffs),as.character(condition),sep="_")
      peaks<-rbind(peaks,data.frame(Position=pos,Chr=chr,Pvalue=pvals,Intensity=intens,Name=names,Dataset=condition))
    }
    peaks
  }

###########################
## Fisher's pvalue combination function
fishers<-function(x)
{
	1-pchisq(-2*sum(log(x)),2*length(x))
}

###########################
## Pipeline Functions

pipeline.annotate.peaks<-
function(peaks)
{
	intersect <- annotate.features(updated.coords,
		peaks$intersect$start,
		peaks$intersect$end,
		peaks$intersect$chr,
		peaks$intersect[,c(-1,-2,-3)],
		names(peaks$intersect[,c(-1,-2,-3)]),
		upstream_window=1000,
		plasmid_map=c("1","PNRC100","PNRC200"))
	
	medichi <- annotate.features(updated.coords,
		peaks$medichi_only$start,
		peaks$medichi_only$end,
		peaks$medichi_only$chr,
		peaks$medichi_only[,c(-1,-2,-3)],
		names(peaks$medichi_only[,c(-1,-2,-3)]),
		upstream_window=1000,
		plasmid_map=c("1","PNRC100","PNRC200"))
	
	ma2c <- annotate.features(updated.coords,
		peaks$ma2c_only$Start,
		peaks$ma2c_only$End,
		peaks$ma2c_only$Chr,
		peaks$ma2c_only[,c(-1,-2,-3)],
		names(peaks$ma2c_only[,c(-1,-2,-3)]),
		upstream_window=1000,
		plasmid_map=c("1","PNRC100","PNRC200"))
	
	ret = list(intersect,medichi,ma2c)
	names(ret)<-c("intersect","medichi","ma2c")
	ret
}

pipeline.align.timecourse.peaks<-
function(noH2O2,wH202_10m,wH2O2_20m,wH2O2_60m)
{
	ret = data.frame()
	for(i in 1:nrow(noH2O2["intersect"]))
	{
		for(k in 1:3)
		{
			if(k==1)
				target = wH2O2_10m
			else if(k==2)
				target = wH2O2_20m
			else
				target = wH2O2_60m
				
			currentMin = 987654321
			hitRow = NULL
			for(j in 1:nrow(target["intersect"]))
			{
				dists = abs(noH2O2["intersect"]$Position-target["intersect"]$Position)
				m = min(dists)
				if(m < currentMin)
				{
					currentMin = m
					
					hit = which(dists<=m)
					hit = hit[1]
					hitRow = target["intersect"][hit,]
				}
			}
			for(j in 1:nrow(target["medichi_only"]))
			{
				dists = abs(noH2O2["intersect"]$Position-target["medichi_only"]$Position)
				m = min(dists)
				if(m < currentMin)
				{
					currentMin = m
					
					hit = which(dists<=m)
					hit = hit[1]
					hitRow = target["medichi_only"][hit,]
				}
			}
			for(j in 1:nrow(target["ma2c_only"]))
			{
				dists = abs(noH2O2["intersect"]$Position-target["ma2c_only"]$Position)
				m = min(dists)
				if(m < currentMin)
				{
					currentMin = m
					
					hit = which(dists<=m)
					hit = hit[1]
					hitRow = target["medichi_only"][hit,]
				}
			}
		}
	}
}

pipeline.select.peaks<-
function(medichi,ma2c,maxGap=120)
{
	medichi_chr = ifelse(medichi$chr=="1",1,ifelse(medichi$chr=="pNRC100",2,3))
	medichi_start = medichi$start
	medichi_end = medichi$end
	medichi_pval = medichi$pvalue
	medichi_dataset = as.character(medichi$datasets)
	ma2c_chr = ifelse(ma2c$Chr=="1",1,ifelse(ma2c$Chr=="PNRC100",2,3))
	ma2c_start = ma2c$Start
	ma2c_end = ma2c$End
	ma2c_fdr = ma2c$FDR
  ma2c_score = ma2c$MA2Cscore
	
	medichi_only = data.frame()
	ma2c_only = data.frame()
	intersect = data.frame()
	ma2c_ignore = c()
	for(i in 1:nrow(medichi))
	{
		if(medichi_pval[i]>.05)
			next;
			
		select = which(ma2c_fdr<.25 & medichi_chr[i] == ma2c_chr & ((medichi_start[i]>=ma2c_start-maxGap & medichi_start[i]<=ma2c_end+maxGap) | (medichi_end[i]>=ma2c_start-maxGap & medichi_end[i]<=ma2c_end+maxGap) ))
		if(length(select)>1)
			select = select[1]
		if(length(select)==0) 
		{
			uniqueDatasets = length(unique(strsplit(medichi_dataset[i],",")[[1]]))
			if(uniqueDatasets>2)
				medichi_only = rbind(medichi_only,medichi[i,]);
		}
		else
		{
			combo = cbind(medichi[i,],ma2c[select,])
			intersect = rbind(intersect,combo)
			ma2c_ignore = c(ma2c_ignore,select)
		}
	}
	for(i in 1:nrow(ma2c))
	{
		if(i %in% ma2c_ignore || ma2c_fdr[i]>0 || ma2c_score[i]<1.5)
			next;
		ma2c_only = rbind(ma2c_only,ma2c[i,])
	}
	ret = list(intersect,medichi_only,ma2c_only)
	names(ret) <- c("intersect","medichi_only","ma2c_only")
	ret
}

pipeline.combine.medichi.peaks<-
function(peaks,maxDist,maxLength,max_iter=-1)
{
	pos = peaks$Position #$medichi_pos
	pval = peaks$Pvalue #medichi_pval
	intensity = peaks$Intensity #medichi_intensity
	name = peaks$Name #medichi_name
	dataset = peaks$Dataset
	
	inf = 987654321
	#calculate the adjacency matrix for a set of peaks
	calcDist<-function(pos,name)
	{
		adj = data.frame()
		for(i in 1:length(pos))
		{
			dist<-abs(pos[i]-pos)
			adj <- rbind(adj,dist)
		}
		colnames(adj)<-name
		rownames(adj)<-name
		adj
	}
	ret = data.frame()
	for(chr in levels(peaks$Chr))
	{
		print(chr)
		#find peaks on this chromosome
		select = peaks$Chr == chr

		if(sum(select)==0)
			next
			
		#calculate distances
		adj<-calcDist(pos[select],name[select])
		
		#Combine peaks
		combined = matrix(rep(0,nrow(adj)*nrow(adj)),nrow=nrow(adj))
		for(i in 1:nrow(adj))
		{
			combined[i,i] = 1
			adj[i,i] = inf
		}
		currentMin = min(adj)
		iter = 1
		while(currentMin<maxDist)
		{
			print(sum(adj<=maxDist)/2)
			for(i in 1:nrow(adj))
			{
				for(j in 1:nrow(adj))
				{
					if(i==j || adj[i,j] > currentMin )
					#if(i==j || adj[i,j] > maxDist )
						next;
			
					peak1 = pos[select]*combined[i,]
					peak1 = peak1[peak1>0]
			
					peak2 = pos[select]*combined[j,]
					peak2 = peak2[peak2>0]
			
					combinedSize = max(peak1,peak2)-min(peak1,peak2)
					#print(peak1)
					#print(peak2)
					#print(combinedSize)
			
					#can we combine these two?
					if(combinedSize < maxLength)
					{
						combo = combined[i,]+combined[j,]
						
						for(i in which(combo>0))
							combined[i,]=combo
						
						#combined[i,] = combo
						#combined[j,] = combo
						#print(i)
						#print(j)
						#print(which(combo>0))
						#print(combined[i,])
						#print(combined[j,])
						
						#don't look at genes that are already combined
						adj[which(combo>0),which(combo>0)]=inf
					}
					else
					{
						#print("fail")
						#print(combinedSize)
						#cat("\n")
						adj[which(combined[i,]>0),j] = inf
						adj[which(combined[j,]>0),i] = inf
					}
				}
			}
			currentMin=min(adj)
			iter<-iter+1
			if(max_iter>0 & iter>max_iter)
				break
		}
		spos = pos[select]
		spval = pval[select]
		sintensity = intensity[select]
		sdata = dataset[select]
		combined = unique(combined)
		for(i in 1:nrow(combined))
		{
			start = min(spos[which(combined[i,]>0)])
			end = max(spos[which(combined[i,]>0)])
			if(start==end)
			{
				start = start-15
				end = end+15
			}
			cpval = fishers(spval[which(combined[i,]>0)])
			cint = mean(sintensity[which(combined[i,]>0)])
			cdata = paste(sdata[which(combined[i,]>0)],collapse=",")
			
			#cat(paste("Combining",paste(which(combined[i,]>0),collapse=","),"\n"))
			#cat(paste("Position:",paste(spos[which(combined[i,]>0)],collapse=" ")," (",cpos,")\n"))
			
			ret = rbind(ret,data.frame(
				chr,
				start,
				end,
				pvalue = cpval,
				intensity = cint,
				datasets = cdata
			))
		}
	}
	ret
}

expandOperons <- function(genes,operons)
{
  operon.names = sapply(as.character(operons$operon),function(x){strsplit(x," ")})
  inOperon <- function(x){sapply(operon.names,function(z) x %in% z)}
  
  sapply(genes,function(x)
    ifelse(any(inOperon(x)),
      operon.names[inOperon(x)],
      x)
    )
}

duplicateOperonRows <- function(table,operons,geneCol="Gene"){
  expand = expandOperons(table[,geneCol],operons)
  
  for(i in 1:length(expand)){
    op = expand[[i]]
    if(length(op)>1){ #in an operon
      for(o in op){
       if(! o %in% table[,geneCol]){
         newrow = table[i,]
         newrow[,geneCol]=o
         table = rbind(table,newrow)
       } 
      }      
    }
  }
  table
}

intersectChipDeg<-function(hits,deg,operons)
{
  
  gene.hits = as.character(hits$gene)
  gene.column = which(names(hits)=="gene")
  operon.names = sapply(as.character(operons$operon),function(x){strsplit(x," ")})
  orig.length = nrow(hits)
  
  #add operon hits
  for(i in 1:orig.length){
    
    search = sapply(operon.names,function(x){gene.hits[i] %in% x}) #check for gene i in all operons
    
    if(any(search)){ #this gene is in an operon
      
      for(op in operon.names[search]){
        
        for(g in op){
          
          annot = gene.coords[gene.coords$canonical_Name==g,] # get the gene's annotation
          
          #create a new row for this gene in hits
          new = hits[i,-gene.column]
          new = cbind(g,new)
          
          names(new) <- names(hits)
          
          hits = rbind(hits,new)
          
        }
        
      }
      
    }
    
  }
  #return hits that match ORFs of DEGs
  hits[sapply(hits$gene,function(x){x %in% deg$ORF}),]
}



# For a list of target files containing MeDiChI fits,
# combine the list of genes hit in each dataset.
# Downstream hits are removed, each row is labeled with the
# dataset it comes from
combineHits<-function(targets,dataset=c(),condition=c(),cutoff=.05,dist=500)
{
  out = data.frame()
  
  if(length(dataset)==0)
    dataset = 1:length(targets)
  
  if(length(condition)==0)
    condition = 1:length(targets)
  
  if("fits" %in% ls())
    rm("fits")
  
  for(i in 1:length(targets))
  {
    target = targets[i]
    load(target)
    if(!"fits" %in% ls()){
      cat(paste("Error: fits not found in",target))
      next;
    }
    
    hits = get.genes.hit(fits,p.cutoff=cutoff,dist.cut=dist)
    hits$Pk.intens = as.double(as.character(hits$Pk.intens))
    hits$Pk.p.val = as.double(as.character(hits$Pk.p.val))
    hits = hits[!hits$Where=="downstream",]
  
    
    if(nrow(hits)==0)
      next;
    
    out = rbind(out,cbind(hits,dataset[i],condition[i]))
      
    rm("fits")
  }
  
  n = names(hits)
  names(out) <- c(n,"dataset","condition")
  out
}

# DEPRICIATED: use combineHitsByDataset or combineHitsByCondition
#
# Create a table of genes hit at each timepoint.
# One row for each gene, with duplicate rows for multiple
# zero time point hits.
# combineTimepoints<-function(targets,dataset,cutoff=.05,dist=500){
#   
#   hits = combineHits(targets,dataset,cutoff=cutoff,dist=dist)
#   data=data.frame()
#   for(gene in unique(hits$Gene)){
#     rows = which(hits$Gene==gene)
#     sets = hits$dataset[rows]
#     ind  = which(dataset %in% sets)
#     
#     if("wH2O2_10m" %in% sets){
#       ind = which(hits$dataset=="wH2O2_10m" & hits$Gene==gene)
#       nonzero = data.frame(hits$Pk.coord[ind],hits$Pk.intens[ind],hits$Pk.p.val[ind])
#     }
#     else{
#       nonzero = data.frame(NA,NA,NA)
#     }
#     
#     if("wH2O2_20m" %in% sets){
#       ind = which(hits$dataset=="wH2O2_20m" & hits$Gene==gene)
#       nonzero = cbind(nonzero,hits$Pk.coord[ind],hits$Pk.intens[ind],hits$Pk.p.val[ind])
#     }
#     else{
#       nonzero = cbind(nonzero,NA,NA,NA)
#     }
#     
#     if("wH2O2_60m" %in% sets){
#       ind = which(hits$dataset=="wH2O2_60m" & hits$Gene==gene)
#       nonzero = cbind(nonzero,hits$Pk.coord[ind],hits$Pk.intens[ind],hits$Pk.p.val[ind])
#     }
#     else{
#       nonzero = cbind(nonzero,NA,NA,NA)
#     }
# 
#     names(nonzero) <- c("10m_h2o2_coord","10m_h2o2_intens","10m_h2o2_pval","20m_h2o2_coord","20m_h2o2_intens","20m_h2o2_pval","60m_h2o2_coord","60m_h2o2_intens","60m_h2o2_pval")
#     reached = F
#     for(s in unique(sets)){
#       if(!s=="wH2O2_10m" && !s=="wH2O2_20m" && !s=="wH2O2_60m"){
#         reached = T
#         ind = which(hits$dataset==s & hits$Gene==gene)
#         
#         nrow = cbind(gene,s,hits$Pk.coord[ind],hits$Pk.intens[ind],hits$Pk.p.val[ind],nonzero)
#         names(nrow) <- c("gene","0m_h2o2_dataset","0m_h2o2_coord","0m_h2o2_intens","0m_h2o2_pval",names(nonzero))
#         data = rbind(data,nrow)
#       }
#     }
#     if(!reached){
#       nrow = cbind(gene,NA,NA,NA,NA,nonzero)
#       names(nrow) <- c("gene","0m_h2o2_dataset","0m_h2o2_coord","0m_h2o2_intens","0m_h2o2_pval",names(nonzero))
#       data = rbind(data,nrow)
#     }
#   }
#   data
# }

  combineGeneHitsByDataset<-function(hits,datasets)
{
  out = hits[hits$dataset==datasets[1],]
  
  
  rename = names(hits)
  rename = rename[! rename == "GeneName"] #all columns other than Gene
  
  names(out)[2:ncol(out)] = paste(datasets[1],rename,sep=".") #set column names to match dataset
  
  dataset.col = which(rename=="dataset")
  out = out[,-(ncol(out)-length(rename)+dataset.col)] # remove dataset column
  
  for(i in 2:length(datasets)){
    out = merge(out,hits[hits$dataset==datasets[i],],by="GeneName",all=T)
    names(out)[(ncol(out)-length(rename)+1):ncol(out)] = paste(datasets[i],rename,sep=".") #set new column names to match dataset
    out = out[,-(ncol(out)-length(rename)+dataset.col)] # remove dataset column
  }
  out
}
 
# Given a list of conditions, combine gene hits within each
# condition using Fisher's pvalue combination and then 
# combine columns between different conditions.
combineGeneHitsByCondition<-function(hits,conditions)
{
  #build a seperate dataframe for each condition that combines all hits to each unique gene 
  cond.hits = list()
  for(cond in unique(conditions)){ # for each unique condition
    cond.hits[[cond]] = data.frame() #create a new dataframe for this condition
    for(gene in unique(hits$Gene[hits$condition==cond])){ #look at each unique gene in this condition
      select = hits$condition==cond & hits$Gene==gene # all rows in this condition with this gene
      
      #combine pvalues
      pvals = hits$Pk.p.val[select]
      combo.pval = fishers(pvals)
      dist = ifelse(hits$Where[select] == "upstream",hits$Distance[select],-hits$Distance[select]) # set distance relative to gene start
      
      #get the chromosome for this gene
      chr = hits$Chr[select][1]
      
      #find average, max, and min for peak coord, distance, and intensity
      avg.coord = mean(hits$Pk.coord[select])
      min.coord = min(hits$Pk.coord[select])
      max.coord = max(hits$Pk.coord[select])
      
      avg.dist = mean(dist)
      min.dist = min(dist)
      max.dist = max(dist)
      
      avg.intens = mean(hits$Pk.intens[select])
      min.intens = min(hits$Pk.intens[select])
      max.intens = max(hits$Pk.intens[select])
      
      #build a new row
      newrow = data.frame(Gene=gene,Pk.pval=combo.pval,Chr=chr,
                          Pk.coord.avg=avg.coord,Pk.coord.min=min.coord,Pk.coord.max=max.coord,
                          Pk.intens.avg=avg.intens,Pk.intens.min=min.intens,Pk.intens.max=max.intens,
                          Dist.avg=avg.dist,Distance.min=min.dist,Distance.max=max.dist)
      
      #add row to this condition
      cond.hits[[cond]] = rbind(cond.hits[[cond]],newrow)
    }
    
    #add 'dataset' label for this condition, used in next step
    cond.hits[[cond]] = cbind(cond.hits[[cond]],dataset=cond)
  }
  
  #build dataframe with all hits from all conditions
  cond.combined = data.frame()
  for(cond in cond.hits){
    cond.combined = rbind(cond.combined,cond)
  }
  
  #build table combining matching genes across conditions
  cond.combined = combineGeneHitsByDataset(cond.combined,names(cond.hits))
  
  cond.combined
}



# Link generic features (requiring only position and chromosome location) to the provided set of genes
# Assuming genes follow dataframe structure of Koiede et al. annotation (updated.coords_ORFs.RData)
# Additional data columns can be supplied (feature_misc & feature_misc_names) that will be appended to output dataframe
# TODO: add operon code
annotate.features<-
function(genes, feature_start, feature_end, feature_chr, feature_misc=NULL, feature_misc_names=NULL, operons=NULL, upstream_window=120, downstream_window=120,feature_size=60,plasmid_map=NULL)
{
	if(length(feature_start)==0)
		return(data.frame())
	#feature_start = feature_pos-feature_size/2
	#feature_end = feature_pos+feature_size/2
	if(!is.null(plasmid_map))
	{
		feature_chr = sapply(feature_chr,function(x){which(plasmid_map==x)})
	}
	
	# Organize gene vectors
	# Gene start is always the lower value (even in reverse strand genes)	
	gene_strand = ifelse(genes$Strand=="+",1,-1)
	gene_start = ifelse(genes$Start<genes$End,genes$Start,genes$End)
	gene_end = ifelse(genes$Start<genes$End,genes$End,genes$Start)
	gene_chr <- ifelse(genes$Chr=="Chr",1,ifelse(genes$Chr=="pNRC100",2,3))
	gene_id <- as.character( genes$vng )
	gene_name <- as.character( genes$Gene )
	is.for <- gene_strand == 1
	#names( gene_start ) <- names( gene_end ) <- names( is.for ) <- names( gene_chr ) <- names( gene_name ) <- gene_id

	add_hits<-
	function(orig,feature,hits,dists,where)
	{
		rbind( orig, cbind(data.frame( Gene=gene_id[ hits ], Distance=round( abs( dists[ hits ] ) ), Where=where, Chr=gene_chr[ hits ],
				     Start=feature_start[feature]),End=feature_end[feature], feature_misc[feature,]) )
	}

	out <- data.frame()
	for( i in 1:length( feature_start ) ) {
		#find upstream
		dists <- ifelse(is.for, gene_start - feature_end[i], feature_start[i] - gene_end) 
		hits <- which( dists>= 0 & abs(dists) <= upstream_window & gene_chr == feature_chr[i] )
		if(length(hits)>0)
		{
			out <- add_hits(out,i,hits,dists,"upstream")
		}
		
		#find downstream
		dists <- ifelse(is.for, feature_start[i] - gene_end, gene_start - feature_end[i]) 
		hits <- which( dists>= 0 & abs(dists) <= downstream_window & gene_chr == feature_chr[i] )
		if(length(hits)>0)
		{
			out <- add_hits(out,i,hits,dists,"downstream")
		}
		
		#find coding
		dists <- ifelse(is.for, gene_start - feature_end[i], feature_start[i] - gene_end) 
		hits = which(gene_chr == feature_chr[i] & ( (gene_start<=feature_end[i] & gene_end>=feature_end[i]) | (gene_start<=feature_start[i] & gene_end>=feature_start[i]) ) )
		if(length(hits)>0)
		{
			out <- add_hits(out,i,hits,dists,"coding")
		}
	}
	if ( nrow( out ) > 0 ) colnames( out ) <- c( "Gene", "Distance", "Where", "Chr",
		                              "Feature.start","Feature.end", feature_misc_names)

	out$Chr = ifelse(out$Chr==1,"1",ifelse(out$Chr==2,"pNRC100","pNRC200"))
	
	#TODO: Some of these genes should be kept
	out = out[out$Gene!="",] #throwout unnamed genes
		                              
	#return(unique(out))
	#handle operons
	old_out = out
	total = 0
	for(operon in as.character(operons$operon))
	{
		names = strsplit(operon," ")[[1]]
		if(!any(gene_id %in% names))
		{
			cat(paste("No genes found! (",names,")\n",sep=""))
			next;
		}
		
		is.for = gene_strand[gene_id %in% names]==1
		if(any(is.for)!=all(is.for))
		{
			cat(paste("Operon",operon,"not same direction!\n"))
		}
		#if(all(is.for))
		#{
		#	min_start = min(gene_start[gene_id %in% names])
		#	first = which(gene_start==min_start)
		#}
		#else
		#{
		#	max_end = max(gene_end[gene_id %in% names])
		#	first = which(gene_start==max_end)
		#}
		
		if(any(names %in% old_out$Gene))
		{
			for(gene in names)
			{
				others = which(names!=gene) 
				index = which(gene_id %in% names[others]) #global index of other genes in operon
				sub_gene_id = gene_id[index]
				sub_gene_start = gene_start[index]
				sub_gene_end = gene_end[index]
				sub.is.for = gene_strand[index]==1
				
				peaks = old_out[old_out$Gene==gene,]
				if(nrow(peaks)==0) next;
				
				sub_peaks_pos = peaks$Pk.coord
				sub_peaks_gene = peaks$Gene
				sub_peaks_where = peaks$Where
				sub_peaks_chr = peaks$Chr
				sub_peaks_intens = peaks$Pk.intens
				sub_peaks_pval = peaks$Pk.p.val
				sub_peaks_ma2c = peaks$MA2C_peak
				sub_peaks_dataset = peaks$Dataset
			
				for(i in 1:nrow(peaks))
				{
					dists = ifelse(sub.is.for, sub_gene_start - sub_peaks_pos[i], sub_peaks_pos[i] - gene_end)
					where = ifelse(sub.is.for,
							#forward strand case
							ifelse(sub_gene_start > sub_peaks_pos[i] + probe.size,
								"upstream",
								ifelse(sub_gene_end < sub_peaks_pos[i] - probe.size,
									"downstream",
									"coding"
									)
								),
							#reverse strand case
							ifelse(sub_gene_end < sub_peaks_pos[i] - probe.size,
								"upstream",
								ifelse(sub_gene_start > sub_peaks_pos[i] + probe.size,
									"downstream",
									"coding")
							
							)
						 )
					
					out <- rbind( out, data.frame( Gene=sub_gene_id, Distance=round( abs( dists) ), 
					     Where= where,
					     Chr=rep(sub_peaks_chr[i],length(index)),
					     Pk.coord=round( rep( sub_peaks_pos[i], length( index ) ) ),
					     Pk.intens=rep( sub_peaks_intens[i], length( index ) ),
					     Pk.p.val=rep( sub_peaks_pval[ i], length( index ) ),
					     MA2C_peak=rep(sub_peaks_ma2c[i],length(index)),
					     Dataset=rep(sub_peaks_dataset[i],length(index)) ) )
				}
				total = total + length(index)*nrow(peaks)
			}
		}
	}
	cat(paste(total,"features added from operons.\n"))
	out = unique(out)
	
	#combine pvals
	combined = data.frame()
	for(gene in unique(out$Gene))
	{
		hits = out[out$Gene == gene,]
		annotations = genes[genes$vng==gene,]
	}
	
	out
}

intersect.peaks.expression<-
function(hits,expression,min.hits)
{
	hit_count = ddply(hits,"Gene","nrow")
	hit_count = hit_count[order(hit_count$nrow,decreasing=T),]
	selected_genes = as.character(hit_count[sapply(strsplit(as.character(hit_count$Gene),'VNG'),length) < 3 & hit_count$nrow>=min.hits,]$Gene)
	select = expression$ORF %in% selected_genes
	overlap = expression[select,]
	#overlap<-cbind(overlap,hits=hit_count[hit_count$Gene %in% selected_genes & hit_count$nrow>=min.hits,"nrow"])
	overlap
}

peaks.neighbor.stats<-
function(peaks,window=120,probe.size=60,probe.delta=30)
{
	peak_start = peaks$medichi_pos-probe.size/2
	peak_end = peaks$medichi_pos+probe.size/2
	peak_mid = peaks$medichi_pos
	peak_chr = ifelse(peaks$chr=="Chr",1,ifelse(peaks$chr=="pNRC100",2,3))
	peak_score = peaks$medichi_intensity
	peak_pval = peaks$medichi_pval
	peak_ma2c = as.character(peaks$ma2c_name)
	peak_dataset = as.integer(peaks$dataset)
	hits = which(!is.na(peaks$ma2c_start)) #which land in ma2c regions?
	misses = which(is.na(peaks$ma2c_start)) #which don't?
	
	stats = data.frame()
	for(i in misses)
	{
		deselect = c(rep(TRUE,i-1),FALSE,rep(TRUE,nrow(peaks)-i)) #don't compare i to itself
		dist = apply(cbind(abs(peak_start[i]-peak_start),abs(peak_end[i]-peak_end)),1,min)
		hits = dist[deselect & dist<=window & peak_chr[i] == peak_chr & is.na(peaks$ma2c_start)]
		if(length(hits)==0)next;
		stats = rbind(stats,data.frame(avg=mean(hits),variance=var(hits),neighbors=length(hits)))
	}
	stats
}

#Run MeDiChI on target Agilent ChIP-chip files
pipeline.agilent.medichi<-
function(outdir="",targets=NULL)
{
	data("halo.hires",package="MeDiChI")

	if (outdir!="" && !file.exists(outdir))
		dir.create(outdir)
	
	datasets = datasets.agilent()
	ids = as.integer(levels(factor(datasets$id)))
	if(!is.null(targets))
	{
		ids = targets
	}
	
	for(id in ids)
	{
		select=datasets$id==id
		files=levels(datasets$file)[datasets$file[select]]
		dyeswap=datasets$dye_swap[select]
		cat(paste("Processing dataset",id,"\n"))
		cat("loading data\n")
		data=data.agilent.medichi(files,dyeswap)
		cat("finding peaks\n")
		fits = deconv.entire.genome(data,max.steps=100, fit.res=30, n.boot=10,boot.sample.opt="residual",kernel=kernel.halo.hires)
		outfile=paste(outdir,"/fits_",id,".RData",sep="")
		cat(paste("Saving to",outfile,"\n"))
		save( fits, file=outfile )
		outfile=paste(outdir,"/medichi_",id,".bed",sep="")
		cat(paste("Saving BED file to",outfile,"\n"))
		medichi.to.bed(fits,outfile)
		
		cat("\n")
	}
}

data.agilent.medichi<-
function(targets,ds=F)
{
	data=data.frame()
	files=data.frame(targets=targets,ds=ds,stringsAsFactors=F)
	for (i in 1:nrow(files))
	{
		target=files[i,1]
		d=files[i,2]
		raw = data.agilent.raw(target)
		if(!d)
			norm = loessNormalize(raw[,"rMedianSignal"],raw[,"gMedianSignal"])
		else
			norm = loessNormalize(raw[,"gMedianSignal"],raw[,"rMedianSignal"])
		newdata = cbind(raw[,c("chr_element","average_coord")],norm[,1]/norm[,2])
		data=rbind(data,newdata)
	}
	names(data)<-c("chr","coord","intensity")
	data
}

data.agilent.raw<-
function(target,ds=FALSE)
{
	loadCoords<-
	function(target)
	{
		return(read.delim (target, sep = "\t"))
	}

	loadRawCc<-
	function(target)
	{
		return(read.delim (target, sep = "\t", header = TRUE, skip = 9))
	}

	cleanRawCc<-
	function(coords,raw,ds=FALSE)
	{
		raw<-cbind(coords,raw)
		raw<-subset(raw,ControlType==0)
		raw<-subset (raw, select = c(chr_element, average_coord, gMedianSignal, rMedianSignal, gMeanSignal, rMeanSignal, gBGMedianSignal, rBGMedianSignal, gBGMeanSignal, rBGMeanSignal))
		rawChr<-subset(raw,chr_element == "Chr")
		raw100<-subset(raw,chr_element=="pNRC100")
		raw200<-subset(raw,chr_element=="pNRC200")
		full<-rbind(rawChr,raw100,raw200)
		return (full)
	}

	raw <- loadRawCc(target)
	coords<- loadCoords("data/Chr_coords.txt")
	return(cleanRawCc(coords,raw,ds))
}


ma2cDataFrame<-
function(rawDatasets,normDatasets)
{	
	rawFiles=rawDatasets$file
	normFiles=normDatasets$file
	combinedDatasets = rbind(rawDatasets,normDatasets)
	
	#build vector of observed log ratios
	ratios = vector()
	lengths = vector()
	for (target in rawFiles)
	{
		raw<-loadMA2C(target)
		ratios = c(ratios,log(raw$IP)-log(raw$IN))
		lengths = c(lengths,nrow(raw))
	}
	for (target in normFiles)
	{
		norm<-loadMA2C(target)
		ratios = c(ratios,norm$N)
		lengths = c(lengths,nrow(norm))
	}
	
	#build factor of experiment names
	factors = vector()
	targets = c(rawFiles,normFiles)
	for(ind in 1:length(targets))
	{
		factors = c(factors,rep(targets[ind],lengths[ind]))
	}
	experiments = factor(factors,targets)
	
	#build of factor marking experiments as normalized/un-normalized
	normalized=vector()
	for(ind in 1:length(rawFiles))
	{
		normalized = c(normalized,rep(F,lengths[ind]))
	}
	for(ind in 1:length(normFiles))
	{
		normalized = c(normalized,rep(T,lengths[ind+length(rawFiles)]))
	}
	normalized = factor(normalized)
	
	#build factor of conditions, dates
	conditions = combinedDatasets$condition
	newCond = vector()
	dates = vector()
	for(ind in 1:length(lengths))
	{
		cond = levels(conditions)[combinedDatasets[ind,"condition"]]
		newCond = c(newCond,rep(cond,lengths[ind]))
		dates = c(dates,rep(combinedDatasets[ind,"date"],lengths[ind]))
	}
	conditions = newCond
	dates = factor(dates)
	
	#build dataframe
	data<-data.frame(ratios,experiments,normalized,conditions,dates)
	return (data)
}

medichi.to.bed<-
function(fits.list,filename="medichi.bed",cutoff=.05)
{
	map <- c( "Chr", "pNRC200", "pNRC100" )
  	if ( class( fits.list ) == "chip.deconv.entire.genome" ) fits.list <- list( fits.list )
  	
  	all<-NULL
  	ind = 0
	for ( fits in fits.list ) {
		for ( where in map ) {
			loc = where
			if(loc=="Chr")
				loc = "1"
			else if(loc=="pNRC200")
				loc="PNRC200"
			else
				loc="PNRC100"
			
			select=coef(fits$fits.fin[[ where ]])[,3]<cutoff
			if(sum(select)==0) #empty
			{
				next;
			}
			else if(sum(select)==1) #simplified to vector
			{
				peak.isb <- fits$fits.fin[[ where ]][[ 1 ]]$coeffs[select,]
				peak.isb.start <- as.integer( peak.isb[1]-30)
				peak.isb.end <- as.integer( peak.isb[1]+30)
				names = paste("peak",ind+1:sum(select)-1,sep="")
				peak.isb <- cbind(loc, peak.isb.start, peak.isb.end, names, peak.isb[2] , coef(fits$fits.fin[[ where ]])[,3])
				colnames( peak.isb ) <- c("LOC","START", "END", "name", "SCORE", "PVALUE")
			}
			else
			{
				peak.isb <- fits$fits.fin[[ where ]][[ 1 ]]$coeffs[select,]
				peak.isb.start <- as.integer( peak.isb[,1]-30)
				peak.isb.end <- as.integer( peak.isb[,1]+30)
				names = paste("peak",ind+1:sum(select)-1,sep="")
				peak.isb <- cbind(loc, peak.isb.start, peak.isb.end, names, peak.isb[,2], coef(fits$fits.fin[[ where ]])[select,3])
				colnames( peak.isb ) <- c("LOC","START", "END", "name", "SCORE", "PVALUE")
			}
			ind=ind+sum(select)
			all<-rbind(all,peak.isb)
		}
	}
	if(!filename=="") #don't actually save the file
	{
		write.table ( subset(all,select=c(LOC,START,END,name,SCORE)), file= filename, col.names=F, row.names=F, quote=F, sep="\t")
	}
	data.frame(all)
}

