# bsdAnalysisFunc.R
# Function definition source file for bsdAnalysis.R
# By Lee Pang - Institute for Systems Biology

# History:
# 20090818	Lee Pang	Created
# 20090819	Lee Pang	Added:
#							getLayout()
#							rmsd()
#							bsdMergeReplicates()
# 20090827	Lee Pang	Fixes:
#							bsdProcess()
#								now has option to omit converting time
#							getMuMax(), getMuMaxTime()
#								were silently grabbing bsd$at
#								getMuMax() doesn't require it - removed
#								getMuMaxTime() has it specified as argument
#						Added:
#							getMuPeaks()
#								supercedes getMuMax() and getMuMaxTime()

## BGN FUNCTION : SIMP
## Simple function that implement's Simpson's rule
## retrieved from: https://stat.ethz.ch/pipermail/r-help/2008-July/166422.html
simp <- function (y, a = NULL, b = NULL, x = NULL, n = 200)
{
     if (is.null(a) | is.null(b)) {
         if (is.null(x))
             stop("No x values provided to integrate over.\n")
     }
     else {
         x <- c(a, b)
     }
     fff <- 1
     if (length(x) == 2) {
         if (x[1] == x[2])
             return(0)
         if (x[2] < x[1]) {
             fff <- -1
             x <- rev(x)
         }
         x <- seq(x[1], x[2], length = n)
         if (is.function(y))
             y <- y(x)
         else {
             cat("y must be a function when x is\n")
             cat("of length equal to 2.\n")
             stop("Bailing out.\n")
         }
         equisp <- TRUE
     }
     else {
         if (is.function(y))
             y <- y(x)
         else if (length(y) != length(x))
             stop("Mismatch in lengths of x and y.\n")
         s <- order(x)
         x <- x[s]
         ddd <- diff(x)
         if (any(ddd == 0))
             stop("Gridpoints must be distinct.\n")
         equisp <- isTRUE(all.equal(diff(ddd), rep(0, length(ddd) - 1)))
         y <- y[s]
     }
     n <- length(x) - 1
     if (equisp) {
         old.op <- options(warn = -1)
         on.exit(options(old.op))
         M <- matrix(y, nrow = n + 2, ncol = 4)[1:(n - 2), ]
         h <- x[2] - x[1]
         fc <- h * c(-1, 13, 13, -1)/24
         aa <- apply(t(M) * fc, 2, sum)
         a1 <- h * sum(y[1:3] * c(5, 8, -1))/12
         an <- h * sum(y[(n - 1):(n + 1)] * c(-1, 8, 5))/12
         return(fff * sum(c(a1, aa, an)))
     }
     m <- n%/%2
     i <- 1:(m + 1)
     a <- x[2 * i] - x[2 * i - 1]
     i <- 1:m
     b <- x[2 * i + 1] - x[2 * i]
     o <- (a[i] * b + 2 * a[i] * a[i] - b * b)/(6 * a[i])
     p <- (a[i] + b)^3/(6 * a[i] * b)
     q <- (a[i] * b + 2 * b * b - a[i] * a[i])/(6 * b)
     k <- numeric(n + 1)
     k[1] <- o[1]
     i <- 1:(m - 1)
     k[2 * i] <- p[i]
     k[2 * i + 1] <- q[i] + o[-1]
     if (n > 2 * m) {
         aa <- a[m + 1]
         bb <- b[m]
         den <- 6 * bb * (bb + aa)
         k[2 * m] <- p[m] - (aa^3)/den
         k[2 * m + 1] <- q[m] + (aa^3 + 4 * bb * aa^2 + 3 * aa *
             bb^2)/den
         k[2 * m + 2] <- (2 * bb * aa^2 + 3 * aa * bb^2)/den
     }
     else {
         k[2 * m] <- p[m]
         k[2 * m + 1] <- q[m]
     }
     fff * sum(k * y)
}
## END FUNCTION

## BGN FUNCTION PEAKDET
peakdet <- function(v, delta, x = NULL) {
#PEAKDET Detect peaks in a vector
#        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
#        maxima and minima ("peaks") in the vector V.
#        MAXTAB and MINTAB consists of two columns. Column 1
#        contains indices in V, and column 2 the found values.
#      
#        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
#        in MAXTAB and MINTAB are replaced with the corresponding
#        X-values.
#
#        A point is considered a maximum peak if it has the maximal
#        value, and was preceded (to the left) by a value lower by
#        DELTA.

# Eli Billauer, 3.4.05 (Explicitly not copyrighted).
# This function is released to the public domain; Any use is allowed.

# 2009-08-15, Lee Pang, function ported to R

	maxtab = NULL
	mintab = NULL

	v = c(v) # Just in case this wasn't a proper vector

	if (is.null(x)) {
	  x = 1:length(v)
	} else { 
	  x = c(x)
	  if (length(v)!= length(x)) {
		stop('Input vectors v and x must have same length')
	  }
	}
	  
	if (length(c(delta))>1) stop('Input argument DELTA must be a scalar');
	if (delta <= 0) stop('Input argument DELTA must be positive')

	mn = Inf
	mx = -Inf
	mnpos = NaN
	mxpos = NaN

	lookformax = TRUE

	for (i in 1:length(v)) {
	  this = v[i]
	  if (this > mx) {
		mx = this
		mxpos = x[i]
		}
	  if (this < mn) {
		mn = this
		mnpos = x[i]
		}
	  
	  if (lookformax) {
		if (this < mx-delta) {
		  maxtab = rbind(maxtab, c(mxpos, mx))
		  mn = this
		  mnpos = x[i]
		  lookformax = FALSE
		}
	  } else {
		if (this > mn+delta) {
		  mintab = rbind(mintab, c(mnpos, mn))
		  mx = this
		  mxpos = x[i]
		  lookformax = TRUE
		}
	  }
	}
	
	if (!is.null(maxtab)) colnames(maxtab) = c('position', 'value')
	if (!is.null(mintab)) colnames(mintab) = c('position', 'value')
	
	return(list(maxima=maxtab, minima=mintab))
}
## END FUNCTION


## BGN FUNCTION bsdProcess
bsdProcess <- function(dataFile, convert.time=TRUE, do.smoothing=TRUE, do.normalize=TRUE, path.corr=2, smooth.win=0.1) {
	# strip out the time strings and convert to something mathematically useable
	# import the raw data file as a data.frame.
	x = read.csv(dataFile, row.names=NULL, as.is=TRUE)
	
	if (convert.time) {
		# Time strings are formatted as h:mm:ss
		# split on ":" and convert numeric values to hours accordingly
		at = unlist(
			lapply(
				lapply(
					strsplit(x$Time, ":"), 
					as.numeric), 
				function(x){
					sum(x*c(1,1/60,1/3600))
				}
			)
		)
	} else {
		# assume that time is already in a usable numeric format
		at = x$Time
	}
	
	# convert the rest of the OD600 data to a matrix so it will be easier to work
	# with mathematically
	x = as.matrix(x[-1])

	# smooth the data using lowess smoothing with a window 10% the length of data.
	if (do.smoothing) {
		xs = apply(x, 2, function(z){return(lowess(at, z, f = smooth.win)$y)})
	} else {
		xs = x
	}
	
	# pathlength correction: approximately 2 to return
	# OD600 values comparable to standard 1cm lengths
	xs = xs*path.corr
	
	# normalize the data by the first point.
	if (do.normalize) {
		xn = apply(xs, 2, function(x){return(x / x[1])})
	} else {
		xn = xs
	}
	
	# calculate the first derivative of the natural log transform
	# to extract instantaneous rates.
	# peaks in this profile are the maximum characteristic rates for
	# each growth phase.
	dX = diff(log(xn))
	dT = matrix(rep(diff(at), dim(dX)[2]), nrow=dim(dX)[1], ncol=dim(dX)[2])
	mu = dX / dT
	
	return(list(at=at, x=x, xs=xs, xn=xn, mu=mu))
}
## END FUNCTION

## BGN FUNCTION getMuMax
# retrieves the maximum growth rate
# delta is the threshold for peak detection - e.g. for delta = 0.1 if a data
# point is 10% of the range greater than it's neighbors it will be flagged as
# a local maxima.  See peakdet() for details.
getMuMax <- function(mu, delta=0.1) {
	return(
		apply(
			mu,
			2,
			function(x){
				sort(
					peakdet(x, diff(range(x))*delta)$maxima[,"value"],
					decreasing=TRUE
				)[1]
			}
		)
	)
}
## END FUNCTION

## BGN FUNCTION getMuPeaks
# Retrieves peaks and their coordinates (time or indices) from instantaneous
# growth rate data.  If a time vector is not specified in 'at', indices are
# returned.  Rate data specified in matrix form is processed on a column
# basis.  Threshold for peak detection 'delta' defaults to 10% of the data
# range in each column.
getMuPeaks <- function(mu, at=NULL, delta=0.1) {
	return(
		apply(mu, 2,
			# on each column c
			function(c){
				# detect peaks
				pd.mx = peakdet(c, diff(range(c))*delta, x=at[-1])$maxima
				# order by decreasing magnitude
				pd.mx = pd.mx[order(pd.mx[,'value'], decreasing=TRUE),]
	}))
}
## END FUNCTION

## BGN FUNCTION getMuMaxTime
# retrieves the time or vector index that the maximum growth rate occured
# delta is the threshold for peak detection - e.g. for delta = 0.1 if a data
# point is 10% of the range greater than it's neighbors it will be flagged as
# a local maxima.  See peakdet() for details.
getMuMaxTime <- function(at, mu, delta=0.1) {
	return(
		apply(
			mu,
			2,
			function(x){
				pd = peakdet(x, diff(range(x))*delta, x=at[-1]) 
				return(pd$maxima[
					sort(pd$maxima[,"value"], decreasing=TRUE,index.return=TRUE)$ix[1],
					"position"
					])
			}
		)
	)
}
## END FUNCTION

## BGN FUNCTION getLayout
# retrieves the layout of the bioscreen data either as a matrix or as a list
# returns unique names and replicate well ids
getLayout <- function(layoutFile, type='matrix') {
	if (type == 'matrix') {
		lyt = as.matrix(read.delim(layoutFile, row.names=1, as.is=T))
		# convert to list format
		
		# prepare the column names - remove the 'X' that was assigned via
		# import as a data.frame
		colnames(lyt) = sub('X', '', colnames(lyt))
		
		# create a matrix of identical size with colIDs
		CIDs = t(array(rep(as.numeric(colnames(lyt)), dim(lyt)[1]), dim(lyt)[2:1]))
		
		# create a matrix of identical size with rowIDs
		RIDs = array(rep(as.numeric(rownames(lyt)), dim(lyt)[2]), dim(lyt))
		
		# merge via addition to create wellIDs
		WIDs = CIDs + RIDs
		
		# convert the wellIDs to character strings
		WIDs = apply(WIDs, 2, as.character)
		
		# stitch the column of wellIDs with the column of culture names list
		# to form the list layout format
		lyt = as.matrix(cbind(c(WIDs), c(lyt)))
		colnames(lyt) = c('Well', 'Name')
		
	} else if (type == 'list') {
		lyt = as.matrix(read.delim(layoutFile, as.is=T))
	} else {
		stop('unrecognized layout type')
	}
	
	lyt.names = unique(lyt[,'Name'])
	lyt.reps = vector('list', length(lyt.names))
	names(lyt.reps) = lyt.names
	for (i in 1:length(lyt.names)) {
		lyt.reps[[i]] = lyt[which(lyt[,'Name'] == unique(lyt[,'Name'])[i]), 'Well']
	}
	
	return(list(Names=lyt.names, Reps=lyt.reps))
}
## END FUNCTION

## BGN FUNCTION rmsd
# calculates the normalized root mean squared deviation between two vectors
rmsd <- function(a, b, normalized = TRUE) {
	value = sqrt(mean((a - b)^2))
	
	if (normalized) {
		value = value / diff(range(a, b))
	}
	
	return(value)
}
## END FUNCTION

## BGN FUNCTION bsdMergeReplicates
bsdMergeReplicates <- function(x, lyt) {
	xbar = NULL
	rep.rmsd = vector('list', length(lyt$Names))
	for (i in 1:length(lyt$Names)) {
		# calculate the mean profile of all replicates
		colIDs = lyt$Reps[[lyt$Names[i]]]
		colIDs = unlist(lapply(colIDs, function(x){paste(c('Well', x), collapse=".")}))
		
		if (min(dim(x)) > 1) {
			xbar = cbind(xbar, apply(x[,colIDs], 1, mean))
			# calculate rmsd values for individual replicate profile relative to the mean profile
			rep.rmsd[[i]] = apply(x[,colIDs], 2, function(a) rmsd(a, xbar[, i]))
			
		} else {
			# avoid the error generated by matices with only one row - ugh!
			xbar = cbind(xbar, mean(x[,colIDs]))
			rep.rmsd[[i]] = rmsd(x[,colIDs], xbar[,i])
		}
	}
	
	# assign appropriate column/object names
	colnames(xbar) = lyt$Names
	names(rep.rmsd) = lyt$Names
	
	return(list(xbar=xbar, rmsd=rep.rmsd))
}

## END FUNCTION