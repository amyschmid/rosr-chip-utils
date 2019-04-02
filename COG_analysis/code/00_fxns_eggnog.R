##eggNOGs custom functions

#useful shorthand
'%ni%' <- Negate('%in%')
colMax <- function(data) sapply(data, max, na.rm = TRUE)
colMin <- function(data) sapply(data, min, na.rm = TRUE)

#uses reutils packages to export fasta sequences to working directory
getseqs <- function(x, database = "protein", filename) {
  #x must be a vector of acc or gi numbers.
  require(reutils)
  path <- file.path(getwd(), filename)
  efetch(x, database, "fasta", "text", outfile = filename)
  print(paste("sequences can be found at ", path))
}

#checks column for duplicates and returns the number of duplicates if present.
any.dupes <- function(x,colname) {
  y <- duplicated(x[[colname]])
  if (nrow(x[y,]) == 0) {
    print(paste("there are 0 duplicates of", length(x[[colname]]), "entries in", colname))
  } else {
    z <- nrow(x[y,])
    print(paste("there are", z, "duplicates in", colname, "out of", length(x[[colname]]), "entries"))
  }
}

#checks that all querys are accounted for in the keyfile
check.gene.list <- function(x,y) {
  #x is gene identifies in cluster files/genelist
  #y is the corresponding column in the key file
  xiny <- length(na.omit(x[x %in% y]))
  yinx <-  length(na.omit(y[y %in% x]))
  if (xiny == yinx) {
    print("all protins are included in the Keyfile")
  }
}

#splits eggNOGs column into individual columns... if multiple entries in NOG database per protein, returns cat version
getNOGs <- function(x, colname, sep = " ", factors = FALSE) {
  lis <- list()
  for (i in 1:length(x[[colname]])) {
    vec <- c(str_extract(x[[colname]][i], "[a-zA-Z0-9]*@halNOG"),
             str_extract(x[[colname]][i], "[a-zA-Z0-9]*@eurNOG"),
             capture.output(cat(unlist(str_extract_all(x[[colname]][i], "[a-zA-Z0-9]*@arNOG")), sep = sep)),
             capture.output(cat(unlist(str_extract_all(x[[colname]][i], "[a-zA-Z0-9]*@NOG")), sep = sep)))
    lis[[i]] <- vec
  }
  df <- as.data.frame(do.call(rbind, lis), stringsAsFactors = FALSE)
  colnames(df) <- c("halNOG", "eurNOG", "arNOG", "NOG")
  df <- cbind("query" = x[["query"]], df, "COG_category" = x[["COG Cat."]], "HMM_description" = x[["eggNOG HMM Desc."]])
  return(df)
}

#function randomly samples without replacement the number NOGs in the peak lists for each species. 
#function then calculates/returns overlap. sample.unique and sample.shared functions represent altered distributions to sample from
sample.NOGs <- function(x, c, samples, nog = "NOGs", list.NOGs = FALSE, seed = TRUE) {
  #x is a vector of df names corresponding to GETNOG outputs of each species (whole genome)
  #c is a vector of the number unique NOGs in each peak.list
  #assumes x and c are in same order
  #samples = the number of times to sample
  
  if (seed == TRUE) {
    set.seed(1234) 
    #this ensures reproducibility when running in multple chunks / notebooks. 
    #If I move inside the for loop, each sampling instances is the same (since it effectively resets the seed every iteration)
    #moved to argument in case I want to change the seed outside of the function. 
  }
  
  #create vectors of unique NOGs
  one <- na.omit(eval(as.name(x[1]))[[nog]])
  two <- na.omit(eval(as.name(x[2]))[[nog]])
  three <- na.omit(eval(as.name(x[3]))[[nog]])
  
  sampled.nogs <- list() #to store overlap numbers
  nogs.intersection <- list() 
  
  for (i in 1:samples){
    
    assign(paste("r.", names(c)[1], sep = ""), sample(one, c[1])) #sample vector of NOGs w/o replacement
    assign(paste("r.", names(c)[2], sep = ""), sample(two, c[2]))
    assign(paste("r.", names(c)[3], sep = ""), sample(three, c[3]))
    #did it the above way for consistency. Now there are three vectors that are samples of the three genomes, called r.hbt, r.hvol, and r.hmed. 
    #Not generlaizable to other genomes currently; tmp vectors names are variable but currently hard coded in.
    
    #calculate overlaps between the three samples and append to list
    hv.hm <- (na.omit(intersect(r.hvol, r.hmed))) #these contain actual NOGs
    hb.hm <- (na.omit(intersect(r.hbt, r.hmed))) 
    hb.hv <- (na.omit(intersect(r.hbt, r.hvol)))
    hb.hv.hm <- intersect(hv.hm, r.hbt) #create all shared
    
    #check intersection code
    #if (length(hb.hv.hm) > 0) {
    #  print(length(hb.hv.hm) == length(intersect(hb.hv, r.hmed)))
    #  print(length(hb.hv.hm) == length(intersect(hb.hm, r.hvol)))}
    
    #redefine pairwise to omit those that are in all shared
    hv.hm <- hv.hm[hv.hm %ni% hb.hv.hm]
    hb.hv <- hb.hv[hb.hv %ni% hb.hv.hm]
    hb.hm <- hb.hm[hb.hm %ni% hb.hv.hm]
    
    #named vector passed names to df as colnames
    overlap <- c("hbt&hvol&hmed"= length(hb.hv.hm),
                 "hbt&hmed" = length(hb.hm),
                 "hbt&hvol" = length(hb.hv),
                 "hvol&hmed" = length(hv.hm), 
                 "hbt" = c[[1]]-length(hb.hv.hm)-length(hb.hm)-length(hb.hv),
                 "hvol"= c[[3]]-length(hb.hv.hm)-length(hv.hm)-length(hb.hv),
                 "hmed"= c[[2]]-length(hb.hv.hm)-length(hb.hm)-length(hv.hm))
    sampled.nogs[[i]] <- overlap
    if (list.NOGs == TRUE) {
      nogs.intersection[[i]] <- hb.hv.hm
    }
  }
  df <- as.data.frame(do.call(rbind, sampled.nogs))
  if (list.NOGs == TRUE) {
    df$intersection.NOGs <- nogs.intersection
  }
  return(df)
}

#function randomly samples unique NOGs from genomes without replacement and calculates/returns the overlap in DF format.
sample.unique.NOGs <- function(x, c, samples, nog = "NOGs", list.NOGs = FALSE, seed = TRUE) {
  #x is a vector of df names corresponding to GETNOG outputs of each species (whole genome)
  #c is a vector of the number unique NOGs in each peak.list
  #assumes x and c are in same order
  #samples = the number of times to sample
  
  if (seed == TRUE) {
    set.seed(1234) 
    #this ensures reproducibility when running in multple chunks / notebooks. 
    #If I move inside the for loop, each sampling instances is the same (since it effectively resets the seed every iteration)
    #moved to argument in case I want to change the seed outside of the function. 
  }
  
  #create vectors of unique NOGs
  one <- na.omit(unique(eval(as.name(x[1]))[[nog]]))
  two <- na.omit(unique(eval(as.name(x[2]))[[nog]]))
  three <- na.omit(unique(eval(as.name(x[3]))[[nog]]))
  
  sampled.nogs <- list() #to store overlap numbers
  nogs.intersection <- list() 
  
  for (i in 1:samples){
    
    assign(paste("r.", names(c)[1], sep = ""), sample(one, c[1])) #sample vector of NOGs w/o replacement
    assign(paste("r.", names(c)[2], sep = ""), sample(two, c[2]))
    assign(paste("r.", names(c)[3], sep = ""), sample(three, c[3]))
    #did it the above way for consistency. Now there are three vectors that are samples of the three genomes, called r.hbt, r.hvol, and r.hmed. 
    #Not generlaizable to other genomes currently; tmp vectors names are variable but currently hard coded in.
    
    #calculate overlaps between the three samples and append to list
    hv.hm <- (na.omit(intersect(r.hvol, r.hmed))) #these contain actual NOGs
    hb.hm <- (na.omit(intersect(r.hbt, r.hmed))) 
    hb.hv <- (na.omit(intersect(r.hbt, r.hvol)))
    hb.hv.hm <- intersect(hv.hm, r.hbt) #create all shared
    
    #check intersection code
    #if (length(hb.hv.hm) > 0) {
    #  print(length(hb.hv.hm) == length(intersect(hb.hv, r.hmed)))
    #  print(length(hb.hv.hm) == length(intersect(hb.hm, r.hvol)))}
    
    #redefine pairwise to omit those that are in all shared
    hv.hm <- hv.hm[hv.hm %ni% hb.hv.hm]
    hb.hv <- hb.hv[hb.hv %ni% hb.hv.hm]
    hb.hm <- hb.hm[hb.hm %ni% hb.hv.hm]
    
    #named vector passed names to df as colnames
    overlap <- c("hbt&hvol&hmed"= length(hb.hv.hm),
                 "hbt&hmed" = length(hb.hm),
                 "hbt&hvol" = length(hb.hv),
                 "hvol&hmed" = length(hv.hm), 
                 "hbt" = c[[1]]-length(hb.hv.hm)-length(hb.hm)-length(hb.hv),
                 "hvol"= c[[3]]-length(hb.hv.hm)-length(hv.hm)-length(hb.hv),
                 "hmed"= c[[2]]-length(hb.hv.hm)-length(hb.hm)-length(hv.hm))
    sampled.nogs[[i]] <- overlap
    if (list.NOGs == TRUE) {
      nogs.intersection[[i]] <- hb.hv.hm
    }
  }
  df <- as.data.frame(do.call(rbind, sampled.nogs))
  if (list.NOGs == TRUE) {
    df$intersection.NOGs <- nogs.intersection
  }
  return(df)
}

#function randomly samples from only NOGs that are conserved in all three species without replacement and calculates/returns the overlap in DF format.
sample.shared.NOGs <- function(x, c, samples, nog = "NOGs", list.NOGs = FALSE, seed = TRUE) {
  #x is a vector of df names corresponding to GETNOG outputs of each species (whole genome)
  #c is a vector of the number unique NOGs in each peak.list
  #assumes x and c are in same order
  #samples = the number of times to sample
  
  if (seed == TRUE) {
    set.seed(1234) 
    #this ensures reproducibility when running in multple chunks / notebooks. 
    #If I move inside the for loop, each sampling instances is the same (since it effectively resets the seed every iteration)
    #moved to argument in case I want to change the seed outside of the function. 
  }
  
  #create vectors of unique NOGs
  one <- na.omit(eval(as.name(x[1]))[[nog]])
  two <- na.omit(eval(as.name(x[2]))[[nog]])
  three <- na.omit(eval(as.name(x[3]))[[nog]])
  
  #master vector of NOGs that are in all genomes
  shared.nogs <- intersect(intersect(one, two), three)
  sampled.nogs <- list() #to store overlap numbers
  nogs.intersection <- list() 
  
  for (i in 1:samples){
    
    assign(paste("r.", names(c)[1], sep = ""), sample(shared.nogs, c[1])) #sample vector of NOGs w/o replacement
    assign(paste("r.", names(c)[2], sep = ""), sample(shared.nogs, c[2]))
    assign(paste("r.", names(c)[3], sep = ""), sample(shared.nogs, c[3]))
    #did it the above way for consistency. Now there are three vectors that are samples of the three genomes, called r.hbt, r.hvol, and r.hmed. 
    #Not generlaizable to other genomes currently; tmp vectors names are variable but currently hard coded in.
    
    #calculate overlaps between the three samples and append to list
    hv.hm <- (na.omit(intersect(r.hvol, r.hmed))) #these contain actual NOGs
    hb.hm <- (na.omit(intersect(r.hbt, r.hmed))) 
    hb.hv <- (na.omit(intersect(r.hbt, r.hvol)))
    hb.hv.hm <- intersect(hv.hm, r.hbt) #create all shared
    
    #check intersection code
    #if (length(hb.hv.hm) > 0) {
    #  print(length(hb.hv.hm) == length(intersect(hb.hv, r.hmed)))
    #  print(length(hb.hv.hm) == length(intersect(hb.hm, r.hvol)))}
    
    #redefine pairwise to omit those that are in all shared
    hv.hm <- hv.hm[hv.hm %ni% hb.hv.hm]
    hb.hv <- hb.hv[hb.hv %ni% hb.hv.hm]
    hb.hm <- hb.hm[hb.hm %ni% hb.hv.hm]
    
    #named vector passed names to df as colnames
    overlap <- c("hbt&hvol&hmed"= length(hb.hv.hm),
                 "hbt&hmed" = length(hb.hm),
                 "hbt&hvol" = length(hb.hv),
                 "hvol&hmed" = length(hv.hm), 
                 "hbt" = c[[1]]-length(hb.hv.hm)-length(hb.hm)-length(hb.hv),
                 "hvol"= c[[3]]-length(hb.hv.hm)-length(hv.hm)-length(hb.hv),
                 "hmed"= c[[2]]-length(hb.hv.hm)-length(hb.hm)-length(hv.hm))
    sampled.nogs[[i]] <- overlap
    if (list.NOGs == TRUE) {
      nogs.intersection[[i]] <- hb.hv.hm
    }
  }
  df <- as.data.frame(do.call(rbind, sampled.nogs))
  if (list.NOGs == TRUE) {
    df$intersection.NOGs <- nogs.intersection
  }
  return(df)
}

#TESTS FOR FUNCTIONAL ENRICHMENT
#performs hypergeometric test on provided subset of genes/proteins relative to the genome
#code adapted from Keely Dulmage
nogtest <- function(namelist,nogfile,pvalue, cutoff = 5) {
  #namelist is a vector of protein on gene names you wnat to test for enrichment
  #nogfile is the genome-wide GETNOG output
  #p-value is significance threshold desired
  #cutoff preset prevents functional categories with less than the designated number of genes/proteins being displayed 
  
  nogs <- nogfile[nogfile[["query"]] %in% namelist,]
  clust <-  table(nogs[["COG_category"]])
  resm <- matrix(0,length(clust),3) #create 0 matrix
  res <- data.frame(resm)  #make 0 df
  rownames(res) <- names(clust)
  colnames(res) <- c("probability", "expected","count")
  all <- table(nogfile[["COG_category"]][nogfile[["COG_category"]] %in% nogs[["COG_category"]]])
  tot <- sum(table(nogfile[["COG_category"]]))
  print(tot); print(all); print(clust)
  for (i in 1:length(clust)){   #calc expected frequencies and pval by hypergeo and append to DF
    
    res[i,1] <- signif(phyper(clust[i], all[i], tot-all[i], nrow(nogs),lower.tail=F), digits = 4)
    res[i,2] <- signif(all[i]*(nrow(nogs)/tot), digits = 4)
    res[i,3] <- clust[i]
  }
  fin <- subset(res, probability <= pvalue & count >= cutoff)
  fin$COG_category <- rownames(fin)
  fin <- fin[, c("COG_category", "probability", "expected", "count")]
  return(fin)
}

#Use the following function to look at the genes in your cluster associated with a particular COG.
nogset= function(namelist,nogfile, cog.category) {
  subset(nogfile, is.element(nogfile$query, namelist) & is.element(nogfile$COG_category, cog.category)==TRUE)
}

#Splits multidomain proteins with multiple functional categories and returns a table of category sums by domain. 
  #MissingasS = TRUE converts any "NA" cog.cat into "S"
all.domains <- function(x, missingasS = FALSE) {
  #x is the a freq table/DF of cog cats
  require("plyr")
  df.list <- list()
  for (i in 1:nrow(x)) {
    a <- unlist(strsplit(as.character(x$COG_category[i]), ", "))
    b <- rep(x$n[i], length(a))
    df <- data.frame("COG_category" = a, "n" = b)
    df.list[[i]] <- df
  }
  bind_rows(df.list) -> df
  if (missingasS == TRUE) {df$COG_category[is.na(df$COG_category)] <- "S" }
  df <- ddply(df, "COG_category", numcolwise(sum))
  return(df)
}

#converts any "NA" cog.cat into "S"
missingasS <- function(x, colname) {
  x[[colname]][is.na(x[[colname]])] <- "S"
  sums <- ddply(x, colname, numcolwise(sum))
  return(sums)
}

#coverts frequency table to percentages
freq.as.percent <- function(x, colname) {
  x[[colname]] <- x[["n"]]/sum(x[["n"]])*100
  x[["n"]] <- NULL
  return(x)
}

sig.list <-  function(x, coglist, clusterlen) {
  require("dplyr")
  df.list <- list()
  for (i in 1:length(coglist)) {
    df <- nogset(x[["query"]], x, coglist[i])
    if (i <= clusterlen) {df$cluster <- rep(1, nrow(df))} 
    else {df$cluster <- rep(2, nrow(df))}
    df.list[[i]] <- df
  }
  fin <- bind_rows(df.list)
  return(fin)
}

#Not used, jaccard function adapatation 
jaccard.RH <- function(a, b, similarity = T) {
  overlap <- intersect(a,b)
  all <- intersection(a,b)
  
  if (similarity == T) {
    if (length(overlap) == 0) {return(0)} else {return(length(overlap) / length(all))}
  } else {
    if (length(overlap) == 0) {return(1)} else {return(1 -(length(overlap) / length(all))) }
  }
} 