##try( source( "chip.deconv.R" ) )
require( MeDiChI )

cat( "MeDiChI utils. (c) David J Reiss, ISB.\nPlease email dreiss@systemsbiology.org for questions or comments.\n" )

cat( "Loading function 'get.strongest.hits'\n" )
## Note p.cutoff can be a p-value, or an intensity (if a non-integer) or a peak count (if an integer)
get.strongest.hits <- function( fits, p.cutoff=0.05 ) {
  if ( ! is.matrix( fits ) ) coeffs <- coef( fits )
  else coeffs <- fits ## Can pass a 2- or 3- column matrix (of coeffs) instead of a fit object.

  has.p <- ncol( coeffs ) == 3

  if ( p.cutoff < 1 ) { ## It's a p-value
    if ( ! has.p ) coeffs <- coeffs[ order( coeffs[ ,2 ], decreasing=T ), ,drop=F ][ 1:( round( nrow( coeffs ) ) *
                                                                             ( 1 - p.cutoff ) ), ,drop=F ]
    else coeffs <- coeffs[ coeffs[ ,3 ] <= p.cutoff, ,drop=F ]
  } else if ( round( p.cutoff ) == p.cutoff ) { ## It's a peak count
    p.cutoff <- min( p.cutoff, nrow( coeffs ) )
    coeffs <- coeffs[ order( coeffs[ ,2 ], decreasing=T ), ,drop=F ][ 1:p.cutoff, ,drop=F ]
  } else if ( is.numeric( p.cutoff ) ) { ## It's an intensity cutoff
    coeffs <- coeffs[ coeffs[ ,2 ] <= p.cutoff, ,drop=F ]
  }
  if ( nrow( coeffs ) <= 0 ) stop( "No hits at that cutoff!" )
  if ( has.p ) coeffs <- coeffs[ order( coeffs[ ,3 ] ), ,drop=F ]
  else coeffs <- coeffs[ order( coeffs[ ,2 ], decreasing=T ), ,drop=F ]
  cat( "N.COEFFS =", nrow( coeffs ), "\n" )
  invisible( coeffs )
}

cat( "Loading function 'get.data'\n" )
get.data <- function( object, ... ) {
  if ( ! is.null( object$data ) ) return( object$data )
  data <- NULL
  for ( chr in names( object$fits.fin ) ) {
    tmp <- object$fits.fin[[ chr ]][[ 1 ]]$data 
    rownames( tmp ) <- rep( chr, nrow( tmp ) )
    data <- rbind( data, tmp )
  }
  data
}

## For a given fit object (using its kernel), get the fit for a given set of coefficients
## across a given coordinates (x.in - default is the object's data coords)
get.fit.for.coeffs <- function( obj, coeffs, x.in=NA, hi.res=1, verbose=F, ... ) {
  if ( is.na( x.in[ 1 ] ) ) x.in <- get.data( obj )[ ,1 ]
  kernel <- obj$kernel; if ( is.null( obj$kernel ) ) kernel <- obj$fits.fin[[ 1 ]][[ 1 ]]$kernel
  coe <- coeffs; coe[ ,1 ] <- coe[ ,1 ] - min( x.in ) + 1
  xx <- make.predictor.matrix( x.in - min( x.in ) + 1, kernel, fit.res=hi.res,
                              good.posns.hires=unique( round( coe[ ,1 ] ) ), 
                              sparse=F, verbose=verbose ) ##fit.bg=fit.bg,
  colnames( xx ) <- as.character( unique( round( coeffs[ ,1 ] ) ) )
  rownames( xx ) <- as.character( x.in )
  pred <- ( xx %*% coeffs[ ,2 ] )[ ,1 ]
  pred
}

cat( "Loading function 'compare.chip.hits'\n" )
compare.chip.hits <- function( fits.1, fits.2, p.cutoffs=c(0.05,0.05), n.rnd=100, within=c(50,100,200),
                              plot.dist.cut=2000, plot.intens=F, reverse=F, main="" ) {

  if ( length( p.cutoffs ) == 1 ) p.cutoffs <- rep( p.cutoffs, 2 )
  coeffs.1 <- get.strongest.hits( fits.1, p.cutoffs[ 1 ] )
  coeffs.2 <- get.strongest.hits( fits.2, p.cutoffs[ 2 ] )
  
  if ( reverse ) {
    tmp <- coeffs.1
    coeffs.1 <- coeffs.2
    coeffs.2 <- tmp
    rm( tmp )
  }

  cat( "COEFFS: FIT1 =", nrow( coeffs.1 ), "; FIT2 =", nrow( coeffs.2 ), "\n" )
  
  dists <- outer( coeffs.1[ ,1 ], coeffs.2[ ,1 ], "-" )
  same.chr <- outer( rownames( coeffs.1 ), rownames( coeffs.2 ), "==" )
  rownames( dists ) <- rownames( same.chr ) <- rownames( coeffs.1 )
  colnames( dists ) <- colnames( same.chr ) <- rownames( coeffs.2 )

  dists.tmp <- round( dists ); dists.tmp[ ! same.chr ] <- 1e15
  n.in <- sapply( within, function( i ) sum( apply( abs( dists.tmp ), 2, min ) <= i, na.rm=T ) )
  cat( "PAIRED:", n.in, "\n" )

  paired <- apply( abs( dists.tmp ), 2, which.min )
  paired.dists <- apply( dists.tmp, 2, function( i ) sign( i[ which.min( abs( i ) ) ] ) * min( abs( i ) ) )
  paired.good <- abs( paired.dists ) <= max( within )

  out.tab <- cbind( coeffs.1[ paired[ paired.good ], ], coeffs.2[ paired.good, ], paired.dists[ paired.good ] )
  ##out.tab <- cbind( coeffs.1[ pairs[ ,1 ], ], coeffs.2[ pairs[ ,2 ], ], distance=dists[ pairs ] )
  ##rownames( out.tab ) <- rownames( coeffs.1 )[ pairs[ ,1 ] ]
  colnames( out.tab ) <- c( "coord.1", "intens.1", "p.val.1", "coord.2", "intens.2", "p.val.2", "distance" )
  if ( reverse ) colnames( out.tab ) <- c( "coord.2", "intens.2", "p.val.2", "coord.1", "intens.1", "p.val.1", "distance" )
  out.tab <- out.tab[ order( abs( out.tab[ ,"distance" ] ) ), ]

  if ( ! is.na( plot.dist.cut ) && any( abs( dists ) <= plot.dist.cut ) ) {
    par( mfrow=c( 2, 1 ) )
    x.hist <- hist( dists[ abs( dists ) <= plot.dist.cut ], breaks=50, plot=F, xlab="Distance", ylab="Count" )
    x.hist <- hist( dists[ abs( dists ) <= plot.dist.cut ], breaks=50, freq=F, main=main,
                   ylim=range( x.hist$density ), xlab="Distance", ylab="Count" )
  }
  
  if ( n.rnd > 0 ) {
    n.better <- within * 0
    all.dists.rnd <- numeric()
    for ( iter in 1:n.rnd ) {
      dists.rnd <- numeric()
      for ( chr in unique( rownames( coeffs.1 ) ) ) {
        coe.2 <- coeffs.2[ rownames( coeffs.2 ) == chr, 1 ]
        if ( length( coe.2 ) <= 0 ) next
        poss.coords.1 <- fits.1$fits.fin[[ chr ]][[ 1 ]]$data[ ,1 ]
        tmp <- round( runif( length( coe.2 ), min=min( poss.coords.1 ), max=max( poss.coords.1 ) ) )
        dists.rnd <- c( dists.rnd, as.vector( outer( coe.2, tmp, "-" ) ) )
      }

      all.dists.rnd <- c( all.dists.rnd, dists.rnd[ abs( dists.rnd ) <= plot.dist.cut + 100 ] )

      n.in.rnd <- sapply( within, function( i ) sum( abs( dists.rnd ) <= i ) )
      n.better <- n.better + as.integer( n.in.rnd >= n.in )
      if ( iter %% 100 == 0 || iter == n.rnd ) cat( "RANDOM TEST:", iter, "|", within, "|", n.in, "|",
                                 n.in.rnd, "|", n.better, "\n" )
    }

    if ( ! is.na( plot.dist.cut ) && any( abs( dists ) <= plot.dist.cut ) &&
        any( abs( all.dists.rnd ) <= plot.dist.cut ) )
      hist( all.dists.rnd[ abs( all.dists.rnd ) <= plot.dist.cut ], breaks=50, main="Random", freq=F,
           ylim=range( x.hist$density ), xlab="Distance", ylab="Count" )

    ##return( invisible( list( dists, dists.rnd ) ) )
  }
  if ( plot.intens ) {
    plot( out.tab[ ,2 ], out.tab[ ,5 ], xlab='Intensity 1', ylab='Intensity 2', pch=19, cex=0.7 )
    cat( "Intensity correlation:", cor( out.tab[ ,2 ], out.tab[ ,5 ], method="spearman" ), "\n" )
    cat( "              p-value:", cor.test( out.tab[ ,2 ], out.tab[ ,5 ] )$p.value, "\n" )
  }
  invisible( list( out.tab=out.tab, coeffs.1=coeffs.1, coeffs.2=coeffs.2, dists=dists ) )
}

cat( "Loading function 'get.genes.hit'\n" )
get.genes.hit <- function( fits, coeffs=NULL, p.cutoff=0.05, dist.cut=2000 ) {
  
  if ( is.null( coeffs ) ) coeffs <- get.strongest.hits( fits, p.cutoff )
  else coeffs <- get.strongest.hits( coeffs, p.cutoff )

  coords <- get.gene.coords.halo()

  is.for <- coords$Orientation == "For"
  start <- as.integer( as.vector( coords$Start ) )
  end <- as.integer( as.vector( coords$Stop ) )
  wheres <- as.character( coords$where )
  genes <- as.character( coords$canonical_Name )
  name <- as.character( coords$Gene_Name )
  chr <- as.character( coords$where )
  names( start ) <- names( end ) <- names( is.for ) <- names( wheres ) <- names( name ) <- names( chr ) <- genes
  
  out <- data.frame()
  for( i in 1:nrow( coeffs ) ) {
    dists <- coeffs[ i, 1 ] - start
    hits <- which( abs( dists ) <= dist.cut & wheres == rownames( coeffs )[ i ] )
    if ( length( hits ) <= 0 ) next

    where.hit <- 
        ifelse( dists[ hits ] <= 0 & is.for[ hits ], "upstream",
            ifelse( dists[ hits ] >= 0 & ! is.for[ hits ], "upstream",
                ifelse( dists[ hits ] > 0 & is.for[ hits ] & dists[ hits ] <= end[ hits ], "coding",
                    ifelse( dists[ hits ] < 0 & ! is.for[ hits ] & dists[ hits ] >= end[ hits ], "coding",
                        ifelse( dists[ hits ] > end[ hits ] & is.for[ hits ], "downstream",
                            ifelse( dists[ hits ] < end[ hits ] & ! is.for[ hits ], "downstream", "" ) ) ) ) ) )

    out <- rbind( out, data.frame( name[ hits ], round( abs( dists[ hits ] ) ), where.hit, chr[ hits ],
                             round( rep( coeffs[ i, 1 ], length( hits ) ) ),
                             sprintf( "%.3f", rep( coeffs[ i, 2 ], length( hits ) ) ),
                             sprintf( "%.5f", rep( coeffs[ i, 3 ], length( hits ) ) ) ) )
  }
  if ( nrow( out ) > 0 ) colnames( out ) <- c( "Gene", "Distance", "Where", "Chr",
                                              "Pk.coord", "Pk.intens", "Pk.p.val" )
  out
}

cat( "Loading function 'chip.hits.to.sif'\n" )
chip.hits.to.sif <- function( fits.list, p.cutoff=0.05, dist.cut=2000, org=NULL, dir.out=NULL,
                             append=T ) {
  if ( class( fits.list ) == "chip.deconv.entire.genome" ) fits.list <- list( TF=fits.list )

  sifs <- list()
  for ( tf.name in names( fits.list ) ) {
    cat( "IP =", tf.name, "\n" )
    fits <- fits.list[[ tf.name ]]
    coeffs <- get.strongest.hits( fits, p.cutoff )
    coords <- get.gene.coords.halo()

    is.for <- coords$Orientation == "For"
    start <- as.integer( as.vector( coords$Start ) )
    end <- as.integer( as.vector( coords$Stop ) )
    wheres <- as.character( coords$where )
    genes <- as.character( coords$canonical_Name )
    names( start ) <- names( end ) <- names( is.for ) <- names( wheres ) <- genes

    sif <- edge.atts.i <- edge.atts.pos <- edge.atts.d <- edge.atts.p <- edge.atts.where <- character()
    for( i in 1:nrow( coeffs ) ) {
      dists <- coeffs[ i, 1 ] - start
      hits <- which( abs( dists ) <= dist.cut & wheres == rownames( coeffs )[ i ] )
      if ( length( hits ) <= 0 ) next
      sif <- c( sif, paste( tf.name, "pd", genes[ hits ] ) )
      edge.atts.i <- c( edge.atts.i, paste( tf.name, "(pd)", genes[ hits ], "=", coeffs[ i, 2 ] ) )
      edge.atts.pos <- c( edge.atts.pos, paste( tf.name, "(pd)", genes[ hits ], "=", round( coeffs[ i, 1 ] ) ) )
      edge.atts.d <- c( edge.atts.d, paste( tf.name, "(pd)", genes[ hits ], "=", abs( dists[ hits ] ) ) )
      where.hit <- 
        ifelse( dists[ hits ] <= 0 & is.for[ hits ], "upstream",
            ifelse( dists[ hits ] >= 0 & ! is.for[ hits ], "upstream",
                ifelse( dists[ hits ] > 0 & is.for[ hits ] & dists[ hits ] <= end[ hits ], "coding",
                    ifelse( dists[ hits ] < 0 & ! is.for[ hits ] & dists[ hits ] >= end[ hits ], "coding",
                        ifelse( dists[ hits ] > end[ hits ] & is.for[ hits ], "downstream",
                            ifelse( dists[ hits ] < end[ hits ] & ! is.for[ hits ], "downstream", "" ) ) ) ) ) )
      edge.atts.where <- c( edge.atts.where, paste( tf.name, "(pd)", genes[ hits ], "=", where.hit ) )
      if ( ncol( coeffs ) == 3 )
        edge.atts.p <- c( edge.atts.p, paste( tf.name, "(pd)", genes[ hits ], "=", coeffs[ i, 3 ] ) )
    }

    out <- list( sif=sif, edge.atts.i=edge.atts.i, edge.atts.d=edge.atts.d, edge.atts.pos=edge.atts.pos,
                edge.atts.where=edge.atts.where )
    if ( length( edge.atts.p ) > 1 ) out$edge.atts.p <- edge.atts.p
    sifs[[ tf.name ]] <- out
  }
  
  if ( ! is.null( dir.out ) ) {
    tf.names <- names( sif )
    sif <- unlist( lapply( sifs, "[[", "sif" ) )
    i.eda <- c( "intensity", unlist( lapply( sifs, "[[", "edge.atts.i" ) ) )
    d.eda <- c( "distance", unlist( lapply( sifs, "[[", "edge.atts.d" ) ) )
    p.eda <- c( "p.value", unlist( lapply( sifs, "[[", "edge.atts.p" ) ) )
    pos.eda <- c( "coordinate", unlist( lapply( sifs, "[[", "edge.atts.pos" ) ) )
    where.eda <- c( "where.hit", unlist( lapply( sifs, "[[", "edge.atts.where" ) ) )
    IPed.noa <- paste( names( sifs ), "= TRUE" )

    load( "data/halo.coords.RData" )
    cnames.noa <- c( "commonName", paste( tf.names, "=", tf.names ),
                    paste( halo.coords$canonical_Name, "=", halo.coords$Gene_Name ) )

    if ( ! file.exists( dir.out ) ) dir.create( dir.out )
    cat( sif, file=paste( dir.out, "/chip_chip_net.sif", sep="" ), append=append, sep="\n", collapse="\n" )
    cat( i.eda, file=paste( dir.out, "/intensity.eda", sep="" ), append=append, sep="\n", collapse="\n" )
    cat( d.eda, file=paste( dir.out, "/distance.eda", sep="" ), append=append, sep="\n", collapse="\n" )
    cat( p.eda, file=paste( dir.out, "/p.value.eda", sep="" ), append=append, sep="\n", collapse="\n" )
    cat( pos.eda, file=paste( dir.out, "/coordinate.eda", sep="" ), append=append, sep="\n", collapse="\n" )
    cat( where.eda, file=paste( dir.out, "/where.hit.eda", sep="" ), append=append, sep="\n", collapse="\n" )
    cat( "was.IPed", IPed.noa, file=paste( dir.out, "/was.IPed.noa", sep="" ), append=append, sep="\n", collapse="\n" )
    cat( cnames.noa, file=paste( dir.out, "/commonNames.noa", sep="" ), append=append, sep="\n", collapse="\n" )
  }    
  
  invisible( sifs )
}

## Re-fit kernel to one or more peaks centered at "centers" using non-linear least squares
re.fit.one.peak <- function( data, centers, intensities, kernel, where=NA, window=1500, 
                            n.boot=1, verbose=F, plot=F ) {
  get.kernel <- function( par ) {
    kern <- kernel
    kern[ ,1 ] <- kern[ ,1 ] + par[ "center" ]
    kern[ ,2 ] <- kern[ ,2 ] * par[ "height" ] + par[ "offset" ]
    kern
  }

  get.fit.y <- function( par, in.x ) {
    y <- yy <- NA
    for ( i in 0:( length( par ) / 3 - 1 ) ) {
      kern <- get.kernel( par[ i*3 + 1:3 ] )
      appr <- approx( kern[ ,1 ], kern[, 2 ], in.x )
      if ( is.na( y ) ) y <- appr$y else y <- y + appr$y
      y <- appr$y
      y[ is.na( y ) ] <- 0
      if ( is.na( yy ) ) yy <- y else yy <- yy + y
    }
    yy
  }
  
  iter <- -1
  func <- function( par, ... ) { ## par is center, height, zero-offset (repeated for multiple peaks)
    yy <- get.fit.y( par, dat[ ,1 ] )
    rss <- sum( ( dat[ ,2 ] - yy )^2 )
    if ( verbose && ( iter <<- iter + 1 ) %% 100 == 0 ) cat( iter + 1, par, rss, "\n" )
    rss
  }

  orig.dat <- data[ data[ ,1 ] >= min( centers ) - window & data[ ,1 ] <= max( centers ) + window, ]
  if ( ! is.na( where ) && ! is.null( rownames( orig.dat ) ) ) orig.dat <- orig.dat[ rownames( orig.dat ) == where, ]
  xx <- min( orig.dat[ ,1 ] ):max( orig.dat[ ,1 ] )
  n.par <- length( centers )
  start.par <- as.vector( rbind( center=centers, height=intensities, offset=rep( 0, n.par ) ) )
  names( start.par ) <- rep( c( "center", "height", "offset" ), n.par )
  
  out.pars <- list()
  for ( boot in 1:n.boot ) {
    dat <- orig.dat
    
    if ( boot > 1 ) {
      ft <- get.fit.y( out.pars[[ 1 ]]$par, dat[ ,1 ] )
      resids <- dat[ ,2 ] - ft
      resids <- resids * sample( c( -( sqrt( 5 ) - 1 ) / 2, ( sqrt( 5 ) + 1 ) / 2 ), length( resids ),
                            prob=c( ( sqrt( 5 ) - 1 ) / ( 2 * sqrt( 5 ) ), ( sqrt( 5 ) + 1 ) / ( 2 * sqrt( 5 ) ) ),
                                replace=T )
      dat[ ,2 ] <- dat[ ,2 ] + resids
    }
  
    if ( plot ) { par( mfrow=c( 2, 1 ) );
                  plot( dat ); lines( xx, get.fit.y( start.par, xx ), col="red" );
                  for ( i in 0:(n.par-1) ) {
                    lines( get.kernel( start.par[ i*3 + 1:3 ] ), col="green" );
                    lines( rep( start.par[ i*3 + 1 ], 2 ), c( 0, start.par[ i*3 + 2 ] ), col="blue" ) } }
  
    x <- optim( start.par, func, method="L-BFGS-B", 
               lower=rep( c( min( dat[ ,1 ] ), 0, min( dat[ ,2 ] ) ), n.par ),
               upper=rep( c( max( dat[ ,1 ] ), max( dat[ ,2 ] ) * 1.5, max( dat[ ,2 ] ) * 0.1 ), n.par ),
               control=list( parscale=rep( c( 1000, 1, 0.1 ), n.par ) ) )
               
    if ( plot ) { plot( dat ); lines( xx, get.fit.y( x$par, xx ), col="red" );
                  for ( i in 0:(n.par-1) ) {
                    lines( get.kernel( x$par[ i*3 + 1:3 ] ), col="green" )
                    lines( rep( x$par[ i*3 + 1 ], 2 ), c( 0, x$par[ i*3 + 2 ] ), col="blue" ) } }
    cat( boot, "\n" ); print( x$par )
    out.pars[[ boot ]] <- x
  }

  if ( plot && n.boot > 1 ) {
    dat <- orig.dat
    pars <- sapply( out.pars, "[[", "par" )
    plot( dat )
    for ( x in out.pars ) {
      lines( xx, get.fit.y( x$par, xx ), col="red" );
      for ( i in 0:(n.par-1) ) lines( get.kernel( x$par[ i*3 + 1:3 ] ), col="green" )
    }
    centers <- as.vector( pars[ rownames( pars ) == "center", ] )
    dens <- density( centers, bw=5, weights=rep( 10, length( centers ) ) )
    lines( dens$x, dens$y, col="blue" )
  }
  out.pars
}

cat( "Loading function 'plot.significant.peaks'\n" )
plot.significant.peaks <- function( fits, p.cutoff=0.05, window=20000, ... ) {
  coeffs <- get.strongest.hits( fits, p.cutoff )
  par( mfrow=c( 2, 1 ) )
  for ( i in 1:nrow( coeffs ) ) {
    cat( coeffs[ i, ], "\n" )
    plot( fits$fits.fin[[ rownames( coeffs )[ i ] ]], center=coeffs[ i, 1 ], window=window,
         main=paste( "P-value =", coeffs[ i, 3 ] ), plot.genes=T, ... )
  }
}

cat( "Loading function 'fraction.in.coding.rgns'\n" )
fraction.in.coding.rgns <- function( fits, p.cutoff=0.05, slop=0, chr.in=NA, ... ) {
  coeffs <- get.strongest.hits( fits, p.cutoff )
  get.fraction.in.coding.rgns( coeffs, slop=slop, chr=chr.in, ... )
}

get.fraction.in.coding.rgns <- function( coeffs, max.npeaks=NA, slop=0, chr=NA ) {
  ##slop <- 20 ## Allow +/- 20 bp "slop" in gene start/stop sites -- see if this does anything.
  coords <- get.gene.coords.halo()
  if ( is.na( max.npeaks ) || max.npeaks > nrow( coeffs ) ) max.npeaks <- nrow( coeffs )
  if ( is.na( chr ) ) chr <- unique( as.character( coords$where ) )
  out <- NULL
  out.chr <- character()
  hits.in <- hits.out <- bg.in <- bg.out <- total.len <- 0
  for ( w in chr ) {
    coe <- coeffs[ rownames( coeffs ) == w, ,drop=F ]
    if ( nrow( coe ) <= 0 ) next
    if ( nrow( coe ) > max.npeaks ) coe <- coe[ order( coe[ ,2 ], decreasing=TRUE ), ][ 1:max.npeaks, ]
    max.np <- nrow( coe )
    cc <- coords[ as.character( coords$where ) == w, ]
    coo <- rep( FALSE, max( c( cc$Start, cc$Stop ) ) )
    for ( i in 1:nrow( cc ) ) coo[ cc$Start[ i ]:cc$Stop[ i ] ] <- TRUE
    if ( slop > 0 ) {
      for ( i in 1:slop ) {
        c.up <- c( FALSE, coo[ 1:( length( coo ) - 1 ) ] )
        c.down <- c( coo[ 2:length( coo ) ], FALSE )
        coo <- coo & c.up & c.down
      }
    }
    ##if ( verbose ) cat( w, sum( coo ) / length( coo ), sum( ! coo ) / length( coo ), "\t" )
    coe[ ,1 ] <- round( coe[ ,1 ] )
    hits.in <- hits.in + sum( coo[ coe[ ,1 ] ], na.rm=T )
    hits.out <- hits.out + sum( ! coo[ coe[ ,1 ] ], na.rm=T )
    bg.in <- bg.in + sum( coo, na.rm=T )
    bg.out <- bg.out + sum( ! coo, na.rm=T )
    total.len <- total.len + length( coo )
##     if ( verbose ) {
##       cat( hits.in, hits.out, "\t" )
##       cat( pbinom( sum( coo[ coe[ ,1 ] ] ), nrow( coe ), sum( coo ) / length( coo ) ), "" )
##       cat( pbinom( sum( ! coo[ coe[ ,1 ] ] ), nrow( coe ), sum( ! coo ) / length( coo ), lower=F ), "\n" )
##     }
##     out.chr <- c( out.chr, w )
##     out <- rbind( out, c( max.np, slop, sum( coo ) / length( coo ), sum( ! coo ) / length( coo ),
##                          hits.in, hits.out,
##                          pbinom( sum( coo[ coe[ ,1 ] ] ), nrow( coe ), sum( coo ) / length( coo ) ),
##                          pbinom( sum( ! coo[ coe[ ,1 ] ] ), nrow( coe ), sum( ! coo ) / length( coo ), lower=F ) ) )
  }

  out <- c( nrow( coeffs ), slop, bg.in / total.len, bg.out / total.len,
           hits.in / nrow( coeffs ), hits.out / nrow( coeffs ),
           log10( pbinom( hits.in, nrow( coeffs ), bg.in / total.len ) ),
           log10( pbinom( hits.out, nrow( coeffs ), bg.out / total.len, lower=F ) ) )
  
  names( out ) <- c( "N. peaks", "Slop", "Expected Out", "Expected In", "Observed Out", "Observed In",
                       "log10-P-value Out", "log10-P-value In" )
  t( t( out ) )
}

get.coding.rgns <- function( slop=0, chr=NA, dirs="BOTH" ) {
  coords <- get.gene.coords.halo()
  if ( is.na( chr ) ) chr <- c( "Chr", "pNRC100", "pNRC200" )
  out <- list()
  for ( w in chr ) {
    cc <- coords[ as.character( coords$where ) == w, ]
    if ( dirs == "FORWARD" ) cc <- cc[ as.character( cc$Orientation ) == "For", ]
    else if ( dirs == "REVERSE" ) cc <- cc[ as.character( cc$Orientation ) == "Rev", ]
    coo <- rep( FALSE, max( c( cc$Start, cc$Stop ) ) )
    for ( i in 1:nrow( cc ) ) coo[ cc$Start[ i ]:cc$Stop[ i ] ] <- TRUE
    if ( slop != 0 ) {
      for ( i in 1:slop ) {
        c.up <- c( FALSE, coo[ 1:( length( coo ) - 1 ) ] )
        c.down <- c( coo[ 2:length( coo ) ], FALSE )
        if ( slop > 0 ) coo <- coo & c.up & c.down
        else if ( slop < 0 ) coo <- coo | c.up | c.down
      }
    }
    out[[ w ]] <- coo
  }
  out
}

#### HALO-SPECIFIC FUNCTIONS HERE!!!

load.data.halo <- function( data, where=NA, wind=NA, verbose=F,
                           data.file=c( isb="data/tfbD_cmyc_Nim_IP_vs_tfbD_cmyc_Nim_wce.clone.gz",
                             nimb="data/tfbd_chip_OID308/GFF/63616_ratio.gff.gz" ),
                           merge.reps=FALSE, cor.cutoff=0.4, ... ) {
  in.window <- function( x, w ) { x >= w[ 1 ] & x <= w[ 2 ] }
  chr.map <- c( Chr="NC_002607", pNRC100="NC_001869", pNRC200="NC_002608" )
  rev.chr.map <- c( NC_002607="Chr", NC_001869="pNRC100", NC_002608="pNRC200" )
  if ( length( data.file ) == 2 ) data.file <- data.file[ data ]
  if ( data == "nimb" ) {
    ## Load tfbD data set
    ##rdata.file <- paste( paste( strsplit( data.file, "/" )[[ 1 ]][ 1:2 ], collapse="/" ), "RData", sep="." )
    rdata.file <- gsub( ".gff.gz", ".RData", data.file )
    if ( file.exists( rdata.file ) ) {
      if ( verbose ) cat( "LOADING NIMB DATA FILE", rdata.file, "\n" ) 
      load( rdata.file ) ## loads 'nimb.chip'
    } else {
      if ( verbose ) cat( "LOADING NIMB DATA FILE", data.file, "\n" ) 
      nimb.chip <- read.delim( gzfile( data.file ), head=F )
      nimb.chip <- nimb.chip[ ,c( 1, 3:6, 9 ) ]
      colnames( nimb.chip ) <- c( "SEQ_ID", "RATIO_NAME", "POSN_START", "POSN_END", "LOG2_RATIO", "PROBE_ID" )
      if ( verbose ) cat( "SAVING NIMB DATA FILE AS", rdata.file, "\n" ) 
      save( nimb.chip, file=rdata.file )
    }

    nimb.in <- nimb.chip
    if ( ! is.na( wind[ 1 ] ) ) nimb.in <- nimb.chip[ in.window( nimb.chip[ ,"POSN_START" ], wind ) |
                                                in.window( nimb.chip[ ,"POSN_END" ], wind ), ]
    nimb.chr <- rev.chr.map[ as.character( nimb.in[ , "SEQ_ID" ] ) ]
    if ( ! is.na( where ) ) nimb.in <- nimb.in[ nimb.chr == where, ]
    nimb.x <- apply( nimb.in[ , c( "POSN_START", "POSN_END" ) ], 1, mean )
    ## Sort the x coords and y values convert y values to intensities (from log2 intensities)
    ## (fwd and rev strand are "replicates" for chip-chip)
    x <- nimb.x
    y <- nimb.in[ , "LOG2_RATIO" ]
    x.inds <- order( x )
    x <- x[ x.inds ]
    y <- 2^y[ x.inds ] ##sapply( seq( 1, length( y ), by=2 ), function( i ) mean( tmp[ i:(i+1) ] ) )
    chr <- rev.chr.map[ as.character( nimb.in[ , "SEQ_ID" ] ) ][ x.inds ]
    rm( nimb.in, nimb.x, nimb.chr, x.inds )
  }

  else if ( data == "isb" ) {
    isb.chip <- NULL
    for ( i in 1:length( data.file ) ) {
      if ( verbose ) cat( "LOADING ISB DATA FILE", data.file[ i ], "\n" )
      if ( length( grep( "gz", data.file[ i ] ) ) >= 1 ) tmp <- read.delim( gzfile( data.file[ i ] ) )
      else tmp <- read.delim( data.file[ i ] )
      colnames( tmp )[ 1 ] <- "GENE_NAME"
      tmp <- tmp[ ! is.na( tmp$chromosome_start_position ) & ! is.na( tmp$chromosome_end_position ), ]
      if ( merge.reps ) { ## Use muRATIO combined ratios from replicates via vera/sam
        if ( verbose ) cat( "Using VERA/SAM MERGED REPS!\n" )
        tmp <- tmp[ , c( "GENE_NAME", "lambda", "muRATIO", "chromosome_start_position",
                        "chromosome_end_position" ) ]
        tmp[ ,"muRATIO" ] <- 10^tmp[ ,"muRATIO" ]
        isb.chip <- rbind( isb.chip, tmp )
      } else { ## Use all replicate ratios individually!!!
        pseudo <- 1
        x.cols <- grep( "^X\\d", colnames( tmp ), perl=T, val=T )
        y.cols <- grep( "^Y\\d", colnames( tmp ), perl=T, val=T )
        f.cols <- grep( "^F\\d", colnames( tmp ), perl=T, val=T )
        if ( verbose ) cat( "Using", length( x.cols ), "INDIVIDUAL UN-MERGED REPS!\n" )
        tmp2 <- data.frame()
        for ( j in 1:length( x.cols ) ) {
          ## Flags: I means included, O means outlier (although we dont know what that means)
          okay <- as.character( tmp[ ,f.cols[ j ] ] ) == "I" ##| as.character( tmp[ ,f.cols[ j ] ] ) == "IO"
          ttmp <- tmp[ okay, ]
          muRATIO <- ( as.numeric( as.vector( ttmp[ ,x.cols[ j ] ] ) ) + pseudo ) /
                     ( as.numeric( as.vector( ttmp[ ,y.cols[ j ] ] ) ) + pseudo )
          tmp.cor <- cor( muRATIO, 10^as.vector( ttmp$muRATIO ), method="spearman" )
          if ( tmp.cor < cor.cutoff )
            muRATIO <- ( as.numeric( as.vector( ttmp[ ,y.cols[ j ] ] ) ) + pseudo ) /
                       ( as.numeric( as.vector( ttmp[ ,x.cols[ j ] ] ) ) + pseudo )
          ##cat(j,range(muRATIO),cor( muRATIO, 10^as.numeric( as.vector( ttmp$muRATIO ) ) ),"\n")
          tmp.cor <- cor( muRATIO, 10^as.numeric( as.vector( ttmp$muRATIO ) ), method="spearman" )
          if ( tmp.cor < cor.cutoff ) {
            cat( "SKIPPING replicate", j, abs( tmp.cor ), "\n" )
            ##plot( muRATIO, 10^as.numeric( as.vector( ttmp$muRATIO ) ) )
            next
          }
          if ( verbose ) cat( "Using replicate", j, abs( tmp.cor ), "\n" )
          tmp2 <- rbind( tmp2, data.frame( GENE_NAME=as.character( ttmp$GENE_NAME ), muRATIO=muRATIO,
                                                    chromosome_start_position=ttmp$chromosome_start_position,
                                                    chromosome_end_position=ttmp$chromosome_end_position ) )
        }
        isb.chip <- rbind( isb.chip, tmp2 )
      }
    }

    isb.in <- isb.chip
    if ( ! is.na( wind[ 1 ] ) ) isb.in <- isb.chip[ in.window( isb.chip[ ,"chromosome_start_position" ], wind ) |
                                              in.window( isb.chip[ ,"chromosome_end_position" ], wind ), ]
    isb.chr <- sapply( strsplit( as.character( isb.in[ ,"GENE_NAME" ] ), "_", fixed=T ), "[[", 1 )
    if ( ! is.na( where ) ) isb.in <- isb.in[ isb.chr == where, ]
    ##isb.x <- apply( isb.in[ , c( "chromosome_start_position", "chromosome_end_position" ) ], 1, mean )
    isb.x <- apply( cbind( as.integer( as.character( isb.in$chromosome_start_position ) ),
                          as.integer( as.character( isb.in$chromosome_start_position ) ) ), 1, mean )
    isb.srt <- order( isb.x )
    x <- isb.x[ isb.srt ]
    ##y <- 10^isb.in[ isb.srt, "muRATIO" ]
    y <- isb.in[ isb.srt, "muRATIO" ]
    chr <- sapply( strsplit( as.character( isb.in[ isb.srt, "GENE_NAME" ] ), "_", fixed=T ), "[[", 1 )
    rm( isb.in, isb.chr, isb.x, isb.srt )
  }
  out <- cbind( x, y )
  rownames( out ) <- chr
  out <- out[ ! is.na( out[ ,1 ] ), ]
  out
}

cat( "Loading function 'medichi.clone.files'\n" )
medichi.clone.files <- function( files, fit.res=30, n.boot=10, verbose=F, ... ) {  
  data <- NULL
  for ( ff in files ) {
    cat( "LOADING:", ff, "\n" )
    data.tmp <- try( load.data.halo( "isb", data.file=ff, verbose=T, cor.cutoff=0.4, ... ) )
    if ( class( data.tmp ) == "try-error" ) { cat( "UH OH (1) - cant load data!\n" ); next }
    if ( nrow( data.tmp ) <= 5 ) { cat( "UH OH (2) - cant load data!\n" ); next }
    if ( unique( rownames( data.tmp ) ) > 3 )
      data.tmp <- data.tmp[ rownames( data.tmp ) %in% c( "Chr", "pNRC100", "pNRC200" ), ]
    if ( nrow( data.tmp ) <= 5 ) { cat( "UH OH (2) - cant load data!\n" ); next }
    data <- rbind( data, data.tmp )
  }
  if ( is.null( data ) || nrow( data ) <= 5 ) { cat( "UH OH (3) - cant load data!\n" ); next }

  ##cm.fit <- chip.deconv( data, where="Chr", fit.res=10, center=650000, wind=20000, max.steps=100, n.boot=10,
  ##                      kernel=kernel.halo.lowres, verbose=T ); plot( cm.fit, boot="scaled.prob", plot.genes=T )

  try( data( "halo.lowres", package="MeDiChI" ) ) ## Load the deconv. kernel and gene.coords (halo)
  fits <- deconv.entire.genome( data=data, max.steps=100, fit.res=fit.res, n.boot=n.boot,
                               boot.sample="residual", kernel=kernel.halo.lowres, verbose=verbose, ... )
  if ( class( fits ) == "try-error" ) { cat( "UH OH - cant load data!\n" ); print( fits ); next }
  ##cat( "Saving fits to", rdata.file, "\n" )
  ##save( data, fits, file=rdata.file, compress=T )
  fits
}

cat( "Loading function 'medichi.to.genome.browser' (thanks, Tie!)\n" )
medichi.to.genome.browser <- function( fits.list, file.prefix="medichi.output", path.to.template="./genomeBrowser" ) {
  options( scipen=999 )
  map <- c( "Chr", "pNRC200", "pNRC100" )
  if ( class( fits.list ) == "chip.deconv.entire.genome" ) fits.list <- list( fits.list )
  if ( ! file.exists( file.prefix ) ) dir.create( file.prefix )

  Lu <-  paste( file.prefix, "data", "peak" )
  for ( fits in fits.list ) {
    for ( where in map ) {
      if ( ! file.exists( paste( file.prefix, where, sep="/" ) ) ) dir.create( paste( file.prefix, where, sep="/" ) )
      dat.isb <- fits$fits.fin[[ where ]][[ 1 ]]$data
      dat.isb.start <- as.integer( dat.isb[,1]-250 )
      dat.isb.end <- as.integer( dat.isb[,1]+250 )
      dat.isb <- cbind( dat.isb.start, dat.isb.end , dat.isb[,2])
      colnames( dat.isb ) <- c("START", "END", file.prefix)
      print( dim(dat.isb) )
      u = dim(dat.isb)[1]
      write.table ( dat.isb, file= paste( file.prefix, "/", where, "/", basename( file.prefix ), "_isb_",
                               where, ".tsv", sep=""), row.names=F, quote=F, sep="\t")
      
      peak.isb <- fits$fits.fin[[ where ]][[ 1 ]]$coeffs
      peak.isb.start <- as.integer( peak.isb[,1])
      peak.isb.end <- as.integer( peak.isb[,1])
      peak.isb <- cbind( peak.isb.start, peak.isb[,2] )
      colnames( peak.isb ) <- c("POSITION", file.prefix)
      print( dim(peak.isb) )  
      v = dim(peak.isb)[1]
      L <- c( paste(file.prefix, "_", "isb_", where, ".tsv", sep=""), u,v )
      Lu = rbind (Lu, L)
      write.table ( peak.isb, file= paste( file.prefix, "/", where, "/", "peak_", basename( file.prefix ),
                                "_isb_", where, ".tsv", sep=""), row.names=F, quote=F, sep="\t")
    }
  }
  ##write.table(Lu, file= paste(file.prefix, "/length_", file.prefix, ".txt", sep=""), col.names=F, row.names=F,
  ##            quote=F, sep="\t")
  tpath <- path.to.template
  file.copy( paste( tpath, "Chr/chromosome_coordinates.txt", sep="/" ),
            paste( file.prefix, "Chr/chromosome_coordinates.txt", sep="/" ) )
  file.copy( paste( tpath, "pNRC200/pnrc200_coordinates.txt", sep="/" ),
            paste( file.prefix, "pNRC200/pnrc200_coordinates.txt", sep="/" ) )
  file.copy( paste( tpath, "pNRC100/pnrc100_coordinates.txt", sep="/" ),
            paste( file.prefix, "pNRC100/pnrc100_coordinates.txt", sep="/" ) )
  file.copy( paste( tpath, "Chr/rna.tsv", sep="/" ), paste( file.prefix, "Chr/rna.tsv", sep="/" ) )
  file.copy( paste( tpath, "Chr/halo1.png", sep="/" ), paste( file.prefix, "Chr/halo1.png", sep="/" ) )
  file.copy( paste( tpath, "pNRC200/halo2.png", sep="/" ), paste( file.prefix, "pNRC200/halo2.png", sep="/" ) )
  file.copy( paste( tpath, "pNRC100/halo3.png", sep="/" ), paste( file.prefix, "pNRC100/halo3.png", sep="/" ) )
  dataset.file <- readLines( paste( tpath, "/halo.dataset.template", sep="" ) )
  dataset.file <- gsub( "LENGTH_CHR", Lu[ 2, 2 ], dataset.file )
  dataset.file <- gsub( "LENGTH_PEAKS_CHR", Lu[ 2, 3 ], dataset.file )
  dataset.file <- gsub( "LENGTH_PNRC200", Lu[ 3, 2 ], dataset.file )
  dataset.file <- gsub( "LENGTH_PEAKS_PNRC200", Lu[ 3, 3 ], dataset.file )
  dataset.file <- gsub( "LENGTH_PNRC100", Lu[ 4, 2 ], dataset.file )
  dataset.file <- gsub( "LENGTH_PEAKS_PNRC100", Lu[ 4, 3 ], dataset.file )
  dataset.file <- gsub( "FILENAME", basename( file.prefix ), dataset.file )
  cat( dataset.file, file=paste( file.prefix, "/", basename( file.prefix ), ".halo.dataset", sep="" ),
      sep="\n", collapse="\n" )
  cat( "\nSaved files to directory", getwd(), "/", file.prefix, "\n", sep="" )
  cat( "\nTo open in Genome Browser, just use 'File -> Load local data set' to open the file:\n",
      getwd(), "/", file.prefix, "/", basename( file.prefix ), ".halo.dataset\n\n", sep="" )
  invisible( Lu )
}

get.gene.coords.halo <- function() {
  if ( ! exists( "gene.coords" ) ) try( data( "halo.lowres", package="MeDiChI" ) )
  if ( exists( "gene.coords" ) ) return( gene.coords )
  try( load( "data/halo.coords.RData" ) ) ## loads "halo.coords"
  halo.coords
}

## ## Remove multiple adjacent coeffs right next to each other. By "coeffs", here, they are actually
## ## output of get.chip.hits()  (see zzz-halo.R)
## filter.coeffs <- function( coeffs, res=40 ) {
##   out <- NULL
##   if ( is.null( rownames( coeffs ) ) ) rownames( coeffs ) <- rep( "XXX", nrow( coeffs ) )
##   for ( w in unique( rownames( coeffs ) ) ) {
##     ##cat( w, "\n" )
##     coe <- coeffs[ rownames( coeffs ) == w, ,drop=F ]
##     coe <- coe[ order( coe[ ,1 ] ), ,drop=F ]
##     neighs <- c( ( coe[ ,1 ] + res ) %in% coe[ ,1 ], FALSE )
##     tmp <- NULL; i <- 1
##     while ( i <= nrow( coe ) ) {
##       if ( neighs[ i ] ) {
##         j <- i
##         while( all( neighs[ j ] ) ) j <- c( j, j[ length( j ) ] + 1 )
##         keep <- j[ which.max( coe[ j, 2 ] ) ]
##         tmp <- rbind( tmp, coe[ keep, ] )
##         i <- max( j ) + 1
##       } else {
##         tmp <- rbind( tmp, coe[ i, ] )
##         i <- i + 1
##       }
##     }
##     rownames( tmp ) <- rep( w, nrow( tmp ) )
##     out <- rbind( out, tmp )
##   }
##   out
## }    

## ## Input should be, e.g. fits$fits.fin from deconv.entire.genome().
## get.chip.hits <- function( obj, p.cutoff=0.9, boot.type="prob", no.boot=F, hi.res=NA, smooth=T ) {
##   coeffs <- NULL
##   in.hi.res <- hi.res
##   if ( ! is.null( obj$coeffs ) ) obj <- list( obj=obj )
##   for ( where in names( obj ) ) {
##     is.boot <- ! is.null( obj[[ where ]][[ 1 ]]$args ) && obj[[ where ]][[ 1 ]]$args$n.boot > 1
##     if ( ! is.boot || no.boot ) { ## Get all hits that are above p.cutoff fraction of the hits
##       cc <- obj[[ where ]]$coeffs
##       if ( ! is.na( p.cutoff ) && require( lattice ) )
##         cc <- cc[ cc[ ,2 ] >= lattice:::fast.quantile( cc[ ,2 ], p.cutoff ), ,drop=F ]
##     } else { ## get all locs that have a boostrap prob above p.cutoff
##       if ( is.na( hi.res ) ) hi.res <- obj[[ where ]][[ 1 ]]$args$fit.res
##       cc <- get.chip.boot.probs( obj[[ where ]], boot.results=boot.type, hi.res=hi.res, smooth=smooth )[ ,c(1,3) ]
##       if ( smooth && ! is.na( in.hi.res ) && hi.res > 1 ) {
##         tmp <- tapply( cc[ ,2 ], round( cc[ ,1 ] / hi.res ) * hi.res, mean, na.rm=T )
##         cc <- cbind( as.numeric( names( tmp ) ), tmp )
##       } else {
##         cc[ ,1 ] <- round( cc[ ,1 ] )
##       }
##       if ( ! is.na( p.cutoff ) ) cc <- cc[ cc[ ,2 ] >= p.cutoff, ,drop=F ]
##     }
##     if ( nrow( cc ) <= 0 ) next
##     rownames( cc ) <- rep( where, nrow( cc ) )
##     coeffs <- rbind( coeffs, cc )
##   }
##   unique( coeffs )
## }

## Compare two chip hit sets (should be output of filter.coeffs(get.chip.hits())... see zzz-halo.R)
## Moved to medichi.utils.R
## compare.two.peaksets <- function( coeffs1, coeffs2, pos.cutoff=500, max.npeaks=c(nrow(coeffs1),nrow(coeffs2)) ) {
##   max.npeaks[ 1 ] <- min( nrow( coeffs1 ), max.npeaks[ 1 ] )
##   max.npeaks[ 2 ] <- min( nrow( coeffs2 ), max.npeaks[ 2 ] )
##   coeffs1 <- coeffs1[ order( coeffs1[ ,2 ], decreasing=T ), ][ 1:max.npeaks[ 1 ], ]
##   coeffs2 <- coeffs2[ order( coeffs2[ ,2 ], decreasing=T ), ][ 1:max.npeaks[ 2 ], ]
##   out <- NULL
##   for ( w in unique( rownames( coeffs1 ) ) ) {
##     cc1 <- coeffs1[ rownames( coeffs1 ) == w, ,drop=F ]
##     cc2 <- coeffs2[ rownames( coeffs2 ) == w, ,drop=F ]
##     for ( i in 1:nrow( cc1 ) ) {
##       f1 <- cc1[ i, ]
##       tmp <- which( abs( cc2[ ,1 ] - f1[ 1 ] ) == min( abs( cc2[ ,1 ] - f1[ 1 ] ), na.rm=T ) )
##       if ( length( tmp ) <= 0 ) next
##       f2 <- cc2[ tmp, ]
##       if ( abs( f1[ 1 ] - f2[ 1 ] ) > pos.cutoff ) next
##       cc2[ tmp, ] <- NA
##       out <- rbind( out, c( f1, f2 ) )
##       rownames( out )[ nrow( out ) ] <- w
##     }
##   }
##   tmp <- out[ abs( out[ ,1 ] - out[ ,3 ] ) < pos.cutoff, ,drop=F ]
##   med <- median( abs( tmp[ ,1 ] - tmp[ ,3 ] ), na.rm=T )
##   mad <- mad( tmp[ ,1 ] - tmp[ ,3 ], na.rm=T )
##   c( med=med, mad=mad, nr1=nrow( coeffs1 ), nr2=nrow( coeffs2 ) )
## }

if ( FALSE ) {
  source("medichi.utils.R")
  load("data/chip_data/trh3.fits.RData")
  fits.3=fits
  load("data/chip_data/trh4.fits.RData")
  fits.4=fits
  load("data/chip_data/trh2.fits.RData")
  fits.2=fits
  source("medichi.utils.R")
  pdf("trh.distances.pdf")
  compare.chip.hits(fits.2,fits.4,p.cut=0.5,n.rnd=100,plot.dist.cut=1000,main="Trh2, Trh4")
  compare.chip.hits(fits.3,fits.4,p.cut=0.5,n.rnd=100,plot.dist.cut=1000,main="Trh3, Trh4")
  compare.chip.hits(fits.3,fits.2,p.cut=0.5,n.rnd=100,plot.dist.cut=1000,main="Trh2, Trh3")
  graphics.off()
}

if ( FALSE ) { ## For Deep's  VNG1179C data
  for ( i in c( list.files( "data/deep/VNG1179C_20060615/", patt="clone.gz", full=T ),
               list.files( "data/deep/VNG1179C_20060711/", patt="clone.gz", full=T ) ) ) {
    if ( file.exists( gsub( "clone.gz", "fits.RData", i ) ) ) next
    fits <- medichi.clone.files( i )
    save( fits, file=gsub( "clone.gz", "fits.RData", i ) )
  }

  pdf()
  for ( i in c( list.files( "data/deep/VNG1179C_20060615/", patt="clone.gz", full=T ),
               list.files( "data/deep/VNG1179C_20060711/", patt="clone.gz", full=T ) ) ) {
    load( gsub( "clone.gz", "fits.RData", i ) )
    cat( i, "\n" )
    print( fraction.in.coding.rgns( fits, p.cutoff=0.05 ) )
    medichi.to.genome.browser( fits, file=gsub( "clone.gz", "gb", i ) )
    plot( fits$fits.fin$Chr, center=528500, window=20000, plot.genes=T, main=basename( i ) )
  }
  graphics.off()
}

